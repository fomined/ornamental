#ifndef MATH_WLSM_H
#define MATH_WLSM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <array>
#include <assert.h>
#include <defs.h>

#include "prg.defs.h"
#include "prg.curve.h"
#include "prg.math.h"

#include "gsl/prg.minimizer.h"
#include "gsl/prg.rootfinder.h"
#include "gsl/prg.interpolation.h"

#include "defs.h"

/*!
 * \brief Описание полинома k-го порядка p0 + p1*x + p2*x^2 + ... + pk*x^k
 * \note полином k-го порядка имеет k+1 параметров, причем параметры полинома не храняться
 *  а заимствуются снаружи путем передачи адреса начала массива параметров p[]
 */
struct Polinome
{
    unsigned power; ///< порядок полинома

    Polinome(unsigned power) : power(power) {}

    /// подсчитать значение полинома в точке, находящейся на расстоянии x от точки разложения
    double operator()(const double p[], double x)
    {
        double sum = p[power];
        for (unsigned i=power; i>0; i--) sum = sum * x + p[i-1];
        return sum;
    }

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

    /// подсчитать значение полинома[его производных k-го порядка] в точке, находящейся
    /// на расстоянии x от точки разложения
    double derivation(unsigned k, const double p[], double x)
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
     * \brief функция, использующаяся для WLSM подбора оптимального полинома, обеспечивающего
     *  минимум дисперсии для случайных данных
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
     * \brief функция, использующаяся для WLSM подбора оптимального полинома, обеспечивающего
     *  минимум дисперсии для случайных данных
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


/*!
 * Структура, обеспечивающая реализацию Weighted Least Square Method для любого полинома.
 */
struct WLSM : protected Polinome
{
    using Polinome::power;
    using Polinome::derivation;

    const double wcoef = 4; // все |x-xx| >= wcoef * w будут игнорированы
        // при этом погрешность при 3w равна ~1%, при 4w - 0.03%

    const std::vector<point<2>> *Y; // массив исходных данных
    double xx; // позиция, относительно которой отсчитываются веса точек
    double width; // ширина окна (коэффициент падения веса точки)

    WLSM(unsigned power, const std::vector<point<2>> *data, double x, double w)
        : Polinome(power), Y(data), xx(x), width(w) {}

    double operator()(const double p[], double x) { return Polinome::operator()(p, x); }
    double operator()(unsigned k, const double p[], double x) { return Polinome::operator()(k, p, x); }

    double f(const double p[])
    {
        unsigned n = unsigned(Y->size());

        double sum = 0;
        double wsum = 0;

        Unsigned ix = Unsigned(std::distance(Y->begin(), std::lower_bound(Y->begin(), Y->end(), xx,
            [](const point<2> &a, double x) { return a[0] < x; }
        )));
        if (ix == n) ix = n - 1; // иначе нельзя сойти с индекса вниз

        double maxwidth = wcoef * width;

        Unsigned ix1 = ix;
        while (ix1 > 0 && std::abs((*Y)[ix1][0] - xx) < maxwidth) ix1--;

        Unsigned ix2 = ix;
        while (ix2 < n && std::abs((*Y)[ix2][0] - xx) < maxwidth) ix2++;

        for (Unsigned i=ix1; i<ix2; i++)
        {
            double x = (*Y)[i][0];
            double y = (*Y)[i][1];
            double w = exp(-0.5 * sqr((x - xx) / width));
            sum += Polinome::f(p, x-xx, y, w);
                // разложения делаем около точки xx, потому использум x-xx
            wsum += w;
        }
        //return sum <= 0 ? 0 : 0.5 * sum / (wsum - power - 1);
        return sum <= 0 ? 0 : 0.5 * sum / wsum;
    }

    double fdf(double g[], const double p[])
    {
        unsigned n = unsigned(Y->size());
        for (unsigned i=0; i<=power; i++) g[i] = 0;

        double sum = 0;
        double wsum = 0;

        Unsigned ix = Unsigned(std::distance(Y->begin(), std::lower_bound(Y->begin(), Y->end(), xx,
            [](const point<2> &a, double x) { return a[0] < x; }
        )));
        if (ix == n) ix = n - 1; // иначе нельзя сойти с индекса вниз

        double maxwidth = wcoef * width;

        Unsigned ix1 = ix;
        while (ix1 > 0 && std::abs((*Y)[ix1][0] - xx) < maxwidth) ix1--;

        Unsigned ix2 = ix;
        while (ix2 < n && std::abs((*Y)[ix2][0] - xx) < maxwidth) ix2++;

        for (Unsigned i=ix1; i<ix2; i++)
        {
            double x = (*Y)[i][0];
            double y = (*Y)[i][1];
            double w = exp(-0.5 * sqr((x - xx) / width));
            sum += Polinome::fdf(p, g, x-xx, y, w);
                // разложения делаем около точки xx, потому использум x-xx
            wsum += w;
        }
        //return sum <= 0 ? 0 : 0.5 * sum / (wsum - power - 1);
        return sum <= 0 ? 0 : 0.5 * sum / wsum;
    }

    /*!
     * \brief выполнить минимизацию коэффициентов полинома и получить дисперсию
     * \param iteration число итераций
     * \param print_flag печатать или нет путь итерирования
     * \param first,last массив парамеров
     * \return средняя дисперсия при заданной ширине окна в заданной точке
     */
    template <typename Iterator>
    flagged_value<double> run(unsigned iteration, Iterator first, Iterator last,
        bool print_flag=false)
    {
        prg::Minimizer<WLSM> minimizer(*this);
        minimizer.get_params().print_flag = print_flag;
        return minimizer.run(iteration, first, last);
    }

    /*!
     * \brief выполнить минимизацию коэффициентов полинома и получить дисперсию данных
     * \param iteration число итераций
     * \param print_flag печатать или нет путь итерирования
     * \param first,last массив парамеров
     * \return средняя дисперсия данных при заданной ширине окна в заданной точке
     */
    template <typename Iterator>
    flagged_value<double> rmsd(unsigned iteration, Iterator first, Iterator last,
        int print_flag=-1)
    {
        prg::Minimizer<WLSM> minimizer(*this);
        minimizer.get_params().print_flag = print_flag;

        flagged_value<double> r = minimizer.run(iteration, first, last);
        return flagged_value<double>(std::sqrt(2 * r.first), r.second);
    }

    flagged_value<double> run(unsigned iteration, bool print_flag=false)
    {
        std::vector<double> p(power+1, 0);
        return run(iteration, &p[0], &p[power+1], print_flag);
    }
};

/*!
 * \brief простое (без связи разложений в соседних точках грида) сглаживание стохастической кривой
 * \param k порядок производной сглаженной кривой (k=0 - сама сглаженная функция)
 * \param power ранг полинома, интерполирующего данные
 * \param fdata стохастические данные
 * \param width ширина окна учета данных
 * \param hx шаг сетки из опорных точек, в которых используется полином (желательно hg < width)
 * \return
 */
template <unsigned n>
curve_t wslm_simple_smooth(unsigned k, unsigned power, const curve<n> &fdata, unsigned ix, unsigned iy,
    double width, double hx)
{
    if (fdata.size() <= 1) return curve_t();

    curve<2> data;
    for (auto e : fdata) data.push_back(point<2>{e[ix], e[iy]});

//    assert(data.size() > 1);
    assert(prg::is_sorted_by_column<2>(data, 0));

    //---------------------------------------------------------------------------------------------
    //      задание функций для оптимизации регрессии методом наименьших квадратов
    //---------------------------------------------------------------------------------------------
    WLSM lsm(power, &data, 0, width);

    double xxa = data.front()[ix];
    double xxb = data.back()[ix];
    double xxbb = xxb + 0.3 * hx;

    curve_t tmp;
    for (double x=xxa; x<xxbb; x+=hx)
    {
        if (x > xxb) x = xxb;
        lsm.xx = x;

        std::vector<double> p(power+1, 0);
        //-----------------------------------------------------------------------------------------
        //                              оптимизация параметров
        //-----------------------------------------------------------------------------------------
        prg::Minimizer<WLSM> minimizer(lsm);
        minimizer.get_params().print_flag = -1;

        flagged_value<double> r = minimizer.rmsd(40, &p[0], &p[power+1]);
        if (!r.second) continue;

        tmp.push_back(point<CDIM>{x, lsm(k, &p[0], 0), r.first});
    }
    return tmp;
}

#endif // MATH_WLSM_H
