#ifndef MATH_WLSM_H
#define MATH_WLSM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <array>
#include <assert.h>

#include "defs.h"
#include "math/curve.h"
#include "gsl/minimizer.h"

namespace prg
{
    /*!
     * \brief Описание полинома k-го порядка p0 + p1*x + p2*x^2 + ... + pk*x^k
     * \note полином k-го порядка имеет k+1 параметров, причем параметры полинома не храняться
     *  а заимствуются снаружи путем передачи адреса начала массива параметров p[]
     */
    struct Polinome
    {
        Unsigned power; ///< порядок полинома

        Polinome(Unsigned power=0) : power(power) {}

        /*!
         * \brief подсчитать значение полинома в точке, находящейся на расстоянии x от точки разложения
         * \param p параметры полинома
         * \param x расстояние от точки разложения полинома
         */
        double operator()(const double p[], double x)
        {
            double sum = p[power];
            for (unsigned i=power; i>0; i--) sum = sum * x + p[i-1];
            return sum;
        }

        /*!
         * \brief подсчитать значение полинома[его производных k-го порядка] в точке,
         *  находящейся на расстоянии x от точки разложения
         * \param p параметры полинома
         * \param x расстояние от точки разложения полинома
         */
        double operator()(unsigned k, const double p[], double x)
        {
            if (k > power) return 0;

            double coef = prg::factorial(power) / prg::factorial(power-k);
            double sum = p[power] * coef;
            for (unsigned i=power-k; i>0; i--)
            {
                coef *= double(i) / (i + k);
                sum = sum * x + p[i+k-1] * coef;
            }
            return sum;
        }

        /*!
         * \brief расчет одного члена функции потерь (дисперсия в заданной точке) при подборе полинома
         * \note используется для WLSM и имеет название совместимое с GSL
         *
         * \param p коэффициенты полинома
         * \param x расстояние от точки данных xi до точки разложения
         * \param y значение в точке данных yi
         * \param w вес точки данных {xi,yi}
         * \return вклад точки в функцию потерь метода WLSM
         */
        double f(const double p[], double x, double y, double w)
        {
            double ff = this->operator()(p, x) - y;
            return sqr(ff) * w;
        }

        /*!
         * \brief расчет одного члена функции потерь и его градиента (дисперсия в заданной точке) при подборе полинома
         * \note название функции обеспечивает совместимость с GSL
         *
         * \param p коэффициенты полинома
         * \param g[out] вклад в градиент функции потерь от точки {xi,yi}
         * \param x расстояние от точки данных xi до точки разложения
         * \param y значение в точке данных yi
         * \param w вес точки данных {xi,yi}
         * \return вклад точки в функцию потерь метода WLSM
         */
        double fdf(const double p[], double g[], double x, double y, double w)
        {
            double ff = this->operator()(p, x) - y;
            double fw = ff * w;
            for (unsigned i=0; i<=power; i++) { g[i] += fw; fw *= x; }
            return sqr(ff) * w;
        }
    };

    struct WLSM_item
    {
        Unsigned power = 2; ///< мощность полинома
        Unsigned iterations = 40; ///< число итераций подгонки полинома
        flagged_value<double> hx = {0.01, true}; ///< шаг сглаживания
        flagged_value<double> width = {0.05, true}; ///< окно сглаживания
        interval<double> xrange={inf<double>(),-inf<double>(0)}; ///< интервал сглаживания
        flagged_value<double> left_indent = {0.1, true}; ///< отступ слева при отсутствии задания интервала
        flagged_value<double> right_indent = {0.1, true}; ///< отступ справа при отсутствии задания интервала
        bool print_flag=false; ///< флаг печати процесса подгонки
    };

    inline std::string to_string(const WLSM_item &item)
    {
        std::string s;
        s += "power\t : " + std::to_string(item.power) + "\n";
        s += "iterations\t : " + std::to_string(item.iterations) + "\n";
        s += "hx\t : " + to_string(item.hx) + "\n";
        s += "width\t : " + to_string(item.width) + "\n";
        s += "xrange\t : " + to_string(item.xrange) + "\n";
        s += "left_indent\t : " + to_string(item.left_indent) + "\n";
        s += "right_indent\t : " + to_string(item.right_indent) + "\n";
        s += "print_flag\t : " + to_string(item.print_flag) + "\n";
        return s;
    }

    /*!
     * \brief Промежуточный объект, который является связкой с минимизаторами GSL и обеспечивает
     *  реализацию Weighted Least Square Method для любого полинома.
     */
    struct WLSM : public Polinome
    {
        using Polinome::power;

        const double wcoef = 4; // все |x-xx| >= wcoef * w будут игнорированы
            // при этом погрешность при 3w равна ~1%, при 4w - 0.03%

        const double *X = nullptr, *Y = nullptr; // массивы исходных данных
        Unsigned n = 0; // размерность массивов исходных данных
        double xx; // точка разложения полинома
        double width; // ширина окна (коэффициент падения веса точки)

    protected:

        friend prg::Minimizer<WLSM>; // обеспечим доступ минимизатора к функциям f(), fdf()

        /*!
         * \brief выполнить минимизацию коэффициентов полинома в заданной точке
         * \note сглаженное значение в точке равно коэффициенту P[0]
         * \note дисперсия сглаженного значения в произвольной точке не оценивается, только в точке данных
         * \note дисперсия сглаженного значения в i-той точке данных равна sqr(P[0] - Y[i])
         *
         * \param iteration число итераций
         * \param print_flag печатать или нет путь итерирования
         * \return минимально возможное значение функции потерь в заданной точке
         */
        flagged_value<double> run(double xx, Unsigned iterations, double *p, bool print_flag=false)
        {
            this->xx = xx; // точка разложения
            prg::Minimizer<WLSM> minimizer(*this);
            minimizer.get_params().print_flag = print_flag;
            return minimizer.run(iterations, &p[0], &p[power+1]);
        }

        /*!
         * \brief расчет значения функции потерь для заданного полинома
         * \note функция используется GSL минимизаторами и имеет совместимое с ними имя с целью оптимизации
         *  коэффициентов полинома для его лучшей подгонки к внешним данным (массивы X,Y)
         * \param p коэффициенты полинома
         * \return значение функции потерь
         */
        double f(const double p[])
        {
            double maxwidth = wcoef * width;
            Unsigned ix1 = std::distance(&X[0], std::lower_bound(&X[0], &X[n], xx-maxwidth));
            Unsigned ix2 = std::distance(&X[0], std::upper_bound(&X[0], &X[n], xx+maxwidth));

            double sum = 0;
            double wsum = 0;
            for (Unsigned i=ix1; i<ix2; i++)
            {
                double w = exp(-0.5 * sqr((X[i] - xx) / width));
                sum += Polinome::f(p, X[i]-xx, Y[i], w);
                    // разложения делаем около точки xx, потому использум x-xx
                wsum += w;
            }
            return sum <= 0 ? 0 : 0.5 * sum / wsum;
        }

        /*!
         * \brief расчет значения функции потерь и ее градиентов для заданного полинома
         * \note функция используется GSL минимизаторами и имеет совместимое с ними имя с целью оптимизации
         *  коэффициентов полинома для его лучшей подгонки к внешним данным (массивы X,Y)
         * \param p коэффициенты полинома
         * \return значение функции потерь
         */
        double fdf(double g[], const double p[])
        {
            for (unsigned i=0; i<=power; i++) g[i] = 0;

            double maxwidth = wcoef * width;
            Unsigned ix1 = std::distance(&X[0], std::lower_bound(&X[0], &X[n], xx-maxwidth));
            Unsigned ix2 = std::distance(&X[0], std::upper_bound(&X[0], &X[n], xx+maxwidth));

            double sum = 0;
            double wsum = 0;
            for (Unsigned i=ix1; i<ix2; i++)
            {
                double w = exp(-0.5 * sqr((X[i] - xx) / width));
                sum += Polinome::fdf(p, g, X[i]-xx, Y[i], w);
                wsum += w;
            }
            return sum <= 0 ? 0 : 0.5 * sum / wsum;
        }

    public:

        /// получение диапазона данных, в котором можно сглаживать внешнюю функцию
        double xmin() const { return X ? X[0] : inf<double>(); }
        double xmax() const { return X ? X[n-1] : -inf<double>(); }

        /*!
         * \brief сгладить стохастическую функцию f
         * \param k порядок производной
         * \param item параметры сглаживания
         * \param ix колонка, где содержится независимая координата
         * \param iy колонка, где содержится зависимая координата
         */
        curve<2> smooth(Unsigned k, const curve<2> &f, const WLSM_item &item=WLSM_item());
        curve<3> smooth_with_errors(Unsigned k, const curve<2> &f, const WLSM_item &item=WLSM_item());

        shared<curve<2>> smooth(Unsigned k, shared<curve<2>> f, const WLSM_item &item=WLSM_item())
            { return shared<curve<2>>(new curve<2>(smooth(k, *f, item))); }

       shared<curve<3>> smooth_with_errors(Unsigned k, shared<curve<2>> f, const WLSM_item &item=WLSM_item())
            { return shared<curve<3>>(new curve<3>(smooth_with_errors(k, *f, item))); }

        template <Unsigned n> curve<2> smooth(Unsigned k, const curve<n> &f, WLSM_item item, Unsigned ix=0, Unsigned iy=1)
        {
            curve<2> g; for (auto e : f) g.push_back({e[ix], e[iy]});
            return smooth(k, g, item);
        }

        /*!
         * \brief сглаживание объекта любого класса при условии наличия в нем функций [resize(), push_back()]
         * \param c контейнер (обычно std::vector)
         * \param op унарная операция, возвращающая std::pair<double, double>{x,y}
         */
        template <typename Container, typename Unary>
        void smooth(Unsigned k, WLSM_item item, Container &c, Unary op);

    };

    inline curve<2> WLSM::smooth(Unsigned k, const curve<2> &f, const WLSM_item &item)
    {
        if (f.size() <= item.power) return f;

        curve<2> out;
        Polinome::power = item.power;
        std::vector<double> P(item.power + 1, 0); // коэффициенты полинома

        //---------------------------------------------------------------------------------------------
        //      задание функций для оптимизации регрессии методом наименьших квадратов
        //---------------------------------------------------------------------------------------------
        std::vector<double> X, Y; for (auto e : f) { X.push_back(e[0]); Y.push_back(e[1]); }

        n = X.size();
        this->X = X.data();
        this->Y = Y.data();

        double xa = 0, xb = 0;
        if (item.xrange.invalid()) { xa = xmin();  xb = xmax(); }
        else { xa = item.xrange.left(); xb = item.xrange.right(); }

        double len = xb - xa; // стартовая длина интервала
        double da = item.left_indent.value();  if (item.left_indent.flag()) da *= len;
        double db = item.right_indent.value();  if (item.right_indent.flag()) db *= len;
        xa -= da;
        xb += db;
        len = xb - xa; // текущая длина интервала с отступами

        width = item.width.flag() == false ? item.width.value() : item.width.value() * len;
            // выбор между абсолютным значением (0) и относительным (1)

        double hx = item.hx.flag() == false ? item.hx.value() :  item.hx.value() * len;

        for (double x=xa; x<=xb; x+=hx)
        {
            auto r = run(x, item.iterations, &P[0], item.print_flag); // оптимизация полинома
            double y = (*this)(k, &P[0], 0); // значение сглаженной функции или ее производной
            if (r.flag()) out.push_back({x, y});
        }
        return out;
    }

    inline curve<3> WLSM::smooth_with_errors(Unsigned k, const curve<2> &f, const WLSM_item &item)
    {
        curve<3> out;
        if (f.size() <= item.power)
        {
            for (auto e : f) out.push_back({e[0], e[1], 0.});
            return out;
        }

        Polinome::power = item.power;
        std::vector<double> P(item.power + 1, 0); // коэффициенты полинома

        //---------------------------------------------------------------------------------------------
        //      задание функций для оптимизации регрессии методом наименьших квадратов
        //---------------------------------------------------------------------------------------------
        std::vector<double> X, Y; for (auto e : f) { X.push_back(e[0]); Y.push_back(e[1]); }

        n = X.size();
        this->X = X.data();
        this->Y = Y.data();

        double xa = 0, xb = 0;
        if (item.xrange.invalid()) { xa = xmin();  xb = xmax(); }
        else { xa = item.xrange.left(); xb = item.xrange.right(); }

        double len = xb - xa; // стартовая длина интервала
        double da = item.left_indent.value();  if (item.left_indent.flag()) da *= len;
        double db = item.right_indent.value();  if (item.right_indent.flag()) db *= len;
        xa -= da;
        xb += db;
        len = xb - xa; // текущая длина интервала с отступами

        width = item.width.flag() == false ? item.width.value() : item.width.value() * len;
            // выбор между абсолютным значением (0) и относительным (1)

        double hx = item.hx.flag() == false ? item.hx.value() :  item.hx.value() * len;

        for (double x=xa; x<=xb; x+=hx)
        {
            auto r = run(x, item.iterations, &P[0], item.print_flag); // оптимизация полинома
            double sigma = std::sqrt(2 * r.value());
                //(!) вычисляется суммарная взвешанная дисперсия в окне

            double y = (*this)(k, &P[0], 0); // значение сглаженной функции или ее производной
            if (r.flag()) out.push_back({x, y, sigma});
        }
        return out;
    }

    template <typename Container, typename Unary>
    void WLSM::smooth(Unsigned k, WLSM_item item, Container &container, Unary op)
    {
        if (container.size() <= item.power) return;

        Polinome::power = item.power;
        std::vector<double> P(power + 1, 0); // коэффициенты полинома

        //---------------------------------------------------------------------------------------------
        //                  считывание данных и параметров сглаживания
        //---------------------------------------------------------------------------------------------
        std::vector<double> X, Y;
        for (auto first=container.begin(); first!=container.end(); ++first)
        {
            auto value = *first; // доступ к элементу неизвестного типа
            paired<double> v = op(value);  // обеспечение извлечения {x,y}
            X.push_back(v.first);
            Y.push_back(v.second);
        }
        container.resize(0); // очистка контейнера

        n = X.size();
        this->X = X.data();
        this->Y = Y.data();

        double xa = 0, xb = 0;
        if (item.xrange.invalid()) { xa = xmin();  xb = xmax(); }
        else { xa = item.xrange.left(); xb = item.xrange.right(); }

        double len = xb - xa; // стартовая длина интервала
        double da = item.left_indent.value();  if (item.left_indent.flag()) da *= len;
        double db = item.right_indent.value();  if (item.right_indent.flag()) db *= len;
        xa -= da;
        xb += db;
        len = xb - xa; // текущая длина интервала с отступами

        width = item.width.flag() == false ? item.width.value() : item.width.value() * len;
            // выбор между абсолютным значением (0) и относительным (1)

        double hx = item.hx.flag() == false ? item.hx.value() :  item.hx.value() * len;
        for (double x=xa; x<=xb; x+=hx)
        {
            auto r = run(x, item.iterations, &P[0]); // оптимизация полинома
            double y = (*this)(k, &P[0], 0); // значение сглаженной функции или ее производной
            if (r.flag()) container.push_back({x, y}); // построение сглаженной функции
        }
    }

} // namespace prg
#endif // MATH_WLSM_H
