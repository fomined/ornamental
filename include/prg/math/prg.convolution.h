#ifndef PRG_CONVOLUTION
#define PRG_CONVOLUTION

#include <array>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <random>
#include <memory>

#include "prg.defs.h"
#include "prg.maxmin.h"
#include "prg.point.h"

namespace prg
{

    typedef std::vector<point<2>> function; // функция y(x)
    typedef std::shared_ptr<function> pfunction; // указатель на функцию y(x)

    /// печать функции в виде [x, y]
    inline void print(const function &v)
    {
        for (auto e : v) std::printf("%lf %lf\n", e[0], e[1]);
        std::cout.flush();
    }
    inline void print(const pfunction &v) { print(*v); }

    /*!
     * \brief вероятность эквивалентности a==0 [1 - да, 0 - нет]
     * \param a сравниваемое значение
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return [0,1] уровень эквивалентности
     */
    inline double zero_grade(double a, double sigma=confidence_interval_sigma)
    {
        assert(sigma > 0);
        double v = (std::abs(a) <= sigma * confidence_interval_multiplier) ?
            exp(-0.5*sqr(a/sigma)) : 0. /* строгий 0 */;
            // возврат строгого нуля (0.) позволяет делат сравнение результата ret != 0
        return v;
    }

    /*!
     * \brief вероятность эквивалентности чисел
     * \param a,b сравниваемые значения
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return [0,1] уровень эквивалентности
     */
    inline double equal_grade(double a, double b, double sigma=confidence_interval_sigma)
    { return zero_grade(a - b, sigma); }

    /*!
     * \brief вероятность циклической эквивалентности чисел
     * \param L длина цикла
     * \param a,b сравниваемые значения
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return [0,1] уровень эквивалентности
     */
    inline double beltway_equal_grade(double a, double b, double L,
        double sigma=confidence_interval_sigma)
    {
        a = mod(a, L); b = mod(b, L); if (a > b) std::swap(a, b);
        double rmin = std::min(b - a, a + L - b);
        return zero_grade(rmin, sigma);
    }

    //--------------------------------------------------------------------------------------------
    //                              интегральные функции
    //--------------------------------------------------------------------------------------------
    /*!
     * \brief находит область данных, которые близки к заданному значению и находит их взвешенную
     *  сумму и сумму весов
     * \note возврат взвешенной суммы и суммы весов позволяет при необходимости иметь как сумму,
     *  так и вычислить взвешенное среднее
     * \param t целевое знчение
     * \param w функция, в которой делается поиск
     * \param L область цикличности
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return сумма значений функции во всех близких точках (с учетом цикличности) + сумма весов
     */
    inline std::pair<double, double> beltway_equal_weigted_sum(double t, pfunction w, double L,
        double sigma=confidence_interval_sigma)
    {
        if (w->empty()) return {0, 0};

        typedef point<2> point;
        t = mod(t, L); // приведение к циклической области

        // основная область попавшая в интервал [0, L]
        auto lowerit = std::lower_bound(w->begin(), w->end(), t,
            [sigma](const point &p, double t){ return less(p[0], t, sigma); });
        auto upperit = std::upper_bound(w->begin(), w->end(), t,
            [sigma](double t, const point &p){ return less(t, p[0], sigma); });

        // найдем сумму значений функции и сумму весов значений в основной области
        double sum = 0, sumw = 0;
        for (auto it=lowerit; it!=upperit; ++it)
        {
            double x = (*it)[0], y = (*it)[1];
            double w = beltway_equal_grade(t, x, L, sigma);
            sum += y * w;
            sumw += w;
        }


        // найдем сумму значений функции и сумму весов значений за правым краем основной области
        if (upperit == w->end())
        {
            upperit = std::upper_bound(w->begin(), w->end(), t - L,
                [sigma](double t, const point &p){ return less(t, p[0], sigma); }
            );

            for (auto it=w->begin(); it!=upperit; ++it)
            {
                double x = (*it)[0], y = (*it)[1];
                double w = beltway_equal_grade(t, x, L, sigma);
                sum += y * w;
                sumw += w;
            }
        }

        // найдем сумму значений функции и сумму весов значений перед левым краем основной области
        if (lowerit == w->begin())
        {
            lowerit = std::lower_bound(w->begin(), w->end(), t + L,
                [sigma](const point &p, double t){ return less(p[0], t, sigma); }
            );
            for (auto it=lowerit; it!=w->end(); ++it)
            {
                double x = (*it)[0], y = (*it)[1];
                double w = beltway_equal_grade(t, x, L, sigma);
                sum += y * w;
                sumw += w;
            }
        }
        return {sum, sumw};
    }

    /*!
     * \brief находит область данных, которые близки к заданному значению
     * \param t целевое знчение
     * \param w функция, где делается поиск
     * \param L область цикличности
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return выделенная область (с учетом цикличности)
     */
    inline pfunction beltway_equal(double t, pfunction w, double L,
        double sigma=confidence_interval_sigma)
    {
        pfunction f(new function);
        if (w->empty()) return f;

        typedef point<2> point;
        t = mod(t, L); // приведение к циклической области

        auto lowerit = std::lower_bound(w->begin(), w->end(), t,
            [sigma](const point &p, double t){ return less(p[0], t, sigma); }
        );

        auto upperit = std::upper_bound(w->begin(), w->end(), t,
            [sigma](double t, const point &p){ return less(t, p[0], sigma);
            }
        );
        f->assign(lowerit, upperit);
            // область попавшая в интервал [0, L]

        if (upperit == w->end())
        {
            upperit = std::upper_bound(w->begin(), w->end(), t - L,
                [sigma](double t, const point &p){ return less(t, p[0], sigma); }
            );
            f->insert(f->end(), w->begin(), upperit);
                // добавляем область в интервале [0, upperit(t - L)]
        }

        if (lowerit == w->begin())
        {
            lowerit = std::lower_bound(w->begin(), w->end(), t + L,
                [sigma](const point &p, double t){ return less(p[0], t, sigma); }
            );
            f->insert(f->end(), lowerit, w->end());
                // добавляем область в интервале [lowerit(t+L), L]
        }
        return f;
    }

    /*!
     * \brief находит циклическую корреляцию двух функций
     * \param w1,w2 заданные функции
     * \param L область цикличности
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return корреляция заданных двух функций
     */
    inline pfunction beltway_correlation(pfunction w1, pfunction w2, double L,
        double sigma=confidence_interval_sigma)
    {
        if (w1->empty() || w2->empty()) return pfunction();

        typedef point<2> point;
        //----------------------------------------------------------------------------------------
        //          найдем смещения, при которых совпадает хотя бы одна линия
        //----------------------------------------------------------------------------------------
        std::set<double> ww;
        for (auto e : *w1)
        {
            double m1 = e[0];
            for (auto e2 : *w2)
            {
                double m2 = e2[0];
                double dm = mod(m2 - m1, L);
                auto pit = std::equal_range(ww.begin(), ww.end(), dm,
                    [sigma](double a, double b) { return less(a, b, sigma); });
                if (pit.first == pit.second) { ww.insert(dm);
                }
                else {
                }
            }
        }

        //----------------------------------------------------------------------------------------
        //          для найденных смещений подсчитаем мощность совпадения линий
        //----------------------------------------------------------------------------------------

        pfunction tmp(new function);
        // первый цикл обязательно по разностям, так как всплески есть именно в них
        for (auto &dm : ww)
        {
            double w = 0;
            for (auto e1 : *w1)
            {
                double m1 = dm + e1[0];
                auto weightsum = beltway_equal_weigted_sum(m1, w2, L, sigma);
                w += e1[1] * weightsum.first;
            }
            tmp->push_back(point{dm, w});
        }
        return tmp;
    }

    /*!
     * \brief находит циклическую конволюцию двух функций
     * \param w1,w2 заданные функции
     * \param L область цикличности
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return корреляция заданных двух функций
     */
    inline pfunction beltway_convolution(pfunction w1, pfunction w2, double L,
        double sigma=confidence_interval_sigma)
    {
        if (w1->empty() || w2->empty()) return pfunction();

        typedef point<2> point;
        //----------------------------------------------------------------------------------------
        //          найдем смещения, при которых совпадает хотя бы одна линия
        //----------------------------------------------------------------------------------------
        std::set<double> ww;
        for (auto e : *w1)
        {
            double m1 = e[0];
            for (auto e2 : *w2)
            {
                double m2 = e2[0];
                double dm = mod(m2 + m1, L);

                auto pit = std::equal_range(ww.begin(), ww.end(), dm,
                    [sigma](double a, double b) { return less(a, b, sigma); });
                if (pit.first == pit.second) ww.insert(dm);
            }
        }

        //----------------------------------------------------------------------------------------
        //          для найденных смещений подсчитаем мощность совпадения линий
        //----------------------------------------------------------------------------------------

        pfunction tmp(new function);
        // первый цикл обязательно по разностям, так как всплески есть именно в них
        for (auto &dm : ww)
        {
            double w = 0;
            for (auto e1 : *w1)
            {
                double m1 = dm - e1[0];
                auto weightsum = beltway_equal_weigted_sum(m1, w2, L, sigma);
                w += e1[1] * weightsum.first;
            }
            tmp->push_back(point{dm, w});
        }
        return tmp;
    }

    /*!
     * \brief находит циклическую автокорреляцию функции
     * \param w заданная функция
     * \param L область цикличности
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return автокорреляция функции
     */
    inline pfunction beltway_autocorrelation(pfunction w, double L,
        double sigma=confidence_interval_sigma)
    { return beltway_correlation(w, w, L, sigma); }

    /*!
     * \brief находит циклическую конволюцию двух функций
     * \param w заданная функция
     * \param L область цикличности
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return автоконволюция функции
     */
    inline pfunction beltway_autoconvolution(pfunction w, double L,
        double sigma=confidence_interval_sigma)
    { return beltway_convolution(w, w, L, sigma); }


    /*!
     * \brief находит циклическую автокорреляцию функции (только в точках самой функции)
     * \param f заданная функция
     * \param L область цикличности
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return автокорреляция заданной функции (только в точках самой функции)
     */
    inline pfunction beltway_restricted_autocorrelation(pfunction f, double L,
        double sigma=confidence_interval_sigma)
    {
        if (f->empty()) return pfunction();
        typedef point<2> point;

        pfunction tmp(new function);
        for (auto e : *f)
        {
            double x = e[0];
            double w = 0;
            for (auto ee : *f)
            {
                double xx = ee[0];
                auto weightsum = beltway_equal_weigted_sum(x + xx, f, L, sigma);
                w += sqrt(ee[1] * weightsum.first);
            }
            tmp->push_back(point{x, w});
        }
        return tmp;
    }

    /*!
     * \brief находит циклическую автоконволюцию функции (только в точках самой функции)
     * \param f заданная функция
     * \param L область цикличности
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return автоконволюция заданной функции (только в точках самой функции)
     */
    inline pfunction beltway_restricted_autoconvolution(pfunction f, double L,
        double sigma=confidence_interval_sigma)
    {
        if (f->empty()) return pfunction();
        typedef point<2> point;

        pfunction tmp(new function);
        for (auto e : *f)
        {
            double x = e[0];
            double w = 0;
            for (auto ee : *f)
            {
                double xx = ee[0];
                auto weightsum = beltway_equal_weigted_sum(x - xx, f, L, sigma);
                w += sqrt(ee[1] * weightsum.first);
            }
            tmp->push_back(point{x, w});
        }
        return tmp;
    }

    /*!
     * \brief находит пересечение функций (только в точках самой функции)
     * \note полагается, что области определения функций совпадают
     * \param f1,f2 заданные функции
     * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
     * \return автоконволюция заданной функции (только в точках самой функции)
     */
    inline pfunction restricted_intersection(pfunction f1, pfunction f2,
        double=confidence_interval_sigma)
    {
        if (f1->empty() || f2->empty()) return pfunction();
        typedef point<2> point;

        pfunction tmp(new function);
        for (Unsigned i1=0,i2=0; i1<f1->size(); i1++, i2++)
        {
            double x1 = (*f1)[i1][0];
            double w = std::min((*f1)[i1][1], (*f2)[i2][1]);
            tmp->push_back(point{x1, w});
        }
        return tmp;
    }

    /*!
     * \brief находит число ненулевых элементов в функции
     * \param f заданная функция
     */
    inline Unsigned unzero_count(pfunction f)
    {
        if (f->empty()) return 0;
        Unsigned cnt = 0;
        for (auto e : *f) if (e[1] != 0.) cnt++;
        return cnt;
    }

}

#endif // PRG_CONVOLUTION
