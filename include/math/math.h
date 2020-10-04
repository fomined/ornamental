#ifndef PRG_MATH
#define PRG_MATH

/**
 * @file хранение стандартных сокращений и макросов
 */

#include <cmath>
#ifdef BOOST_USED
    #include <boost/math/special_functions/factorials.hpp>
    #include <boost/math/special_functions/binomial.hpp>
#endif
#include "defs.h"

namespace prg
{

    /// неопределенное значение величины
    constexpr double undef_value = infinity<double>();

    /// функция f(x) = sign(x), возвращает: -1 в случае |x| != x и 1 если |x| == x
    template <typename T> inline constexpr int sign(T a) { return a < 0 ? -1 : 1; }

    /// значение для которого |x| == x
    template <typename T> inline constexpr bool positive(T a) { return a < 0 ? false : true; }

    /// значение для которого |x| != x
    template <typename T> inline constexpr bool negative(T a) { return a < 0 ? true : false; }


#ifdef BOOST_USED
    /*!
     * \brief получить факториал $n!$
     * \param x аргумент функции [тип unsigned согласно boost::]
     * \note возвращаемое значение всегда double, иначе переполнение уже при малых n
     */
    inline double factorial(Unsigned x) { return boost::math::factorial<double>(x); }

    /*!
     * \brief получить биномиальный коэффициент $C^m_n$
     * \note возвращаемое значение всегда double, иначе переполнение уже при малых n
     */
    inline double binomial(Unsigned n, Unsigned m)
        { return n >= m ? boost::math::binomial_coefficient<double>(n, m) : 0; }
#else
    inline double factorial(Unsigned x) { return x <= 1 ? 1 : double(x) * factorial(x-1); }
#endif

    /*!
     * \brief значение функции Пуассона $e^-m m^k / k!$
     * \param mean среднее значение
     * \param k точка, в которой считается функция
     */
    inline double poisson(double mean, Unsigned k)
        { return std::exp(-mean) * std::pow(mean, k) / factorial(k); }

    /// функция сигмодиды f(x) = 1/(1+exp(-(x-x0)/s)
    inline double sigmoid(double x, double x0, double sigma) { return 1./(1 + exp(-(x-x0)/sigma)); }

    /// сигмодида и ее производные по параметрам {f, df/dx, df/dx0, df/dsigma}
    inline std::array<double, 4> grad_sigmoid(double x, double x0, double sigma)
    {
        double _1sigma = 1. / sigma;
        double m = (x - x0) * _1sigma;
        double ex = exp(-m);
        double f = 1. / (1 + ex);
        double ff = f * (f - 1) * _1sigma;
        std::array<double, 4> a = {{ f, -ff, ff, ff * m }};
        return a;
    }

    /// функция обратной сигмодиды f(x) = 1/(1+exp((x-x0)/s)
    inline double rsigmoid(double x, double x0, double sigma) { return 1./(1 + exp((x-x0)/sigma)); }

    /// сигмодида и ее производные по параметрам {f, df/dx, df/dx0, df/dsigma}
    inline std::array<double, 4> grad_rsigmoid(double x, double x0, double sigma)
    {
        double _1sigma = 1. / sigma;
        double m = (x - x0) * _1sigma;
        double ex = exp(m);
        double f = 1. / (1 + ex);
        double ff = f * (f - 1) * _1sigma;
        std::array<double, 4> a = {{ f, ff, -ff, -ff * m }};
        return a;
    }

    /// функция двойной сигмодиды f(x) = s(x) + rs(x) - 1;
    inline double sigmoid2s(double x, double x0, double sigma0, double x1, double sigma1)
    { return sigmoid(x, x0, sigma0) + rsigmoid(x, x1, sigma1) - 1; }

    /// 2-ная сигмодида и ее производные по параметрам {f, df/dx, df/dx1, df/ds1, df/dx2, df/ds2}
    inline std::array<double, 6> grad_sigmoid2s(double x, double x1, double sigma1, double x2, double sigma2)
    {
        std::array<double, 4> a = grad_sigmoid (x, x1, sigma1);
        std::array<double, 4> b = grad_rsigmoid(x, x2, sigma2);
        std::array<double, 6> aa = {{ a[0] + b[0] - 1, a[1] + b[1], a[2], a[3], b[2], b[3] }};
        return aa;
    }

    /// функция двойной сигмодиды f(x) = s(x) * rs(x);
    inline double sigmoid2p(double x, double x0, double sigma0, double x1, double sigma1)
    { return sigmoid(x, x0, sigma0) * rsigmoid(x, x1, sigma1); }

    /// 2-ная сигмодида и ее производные по параметрам {f, df/dx, df/dx1, df/ds1, df/dx2, df/ds2}
    inline std::array<double, 6> grad_sigmoid2p(double x, double x1, double sigma1, double x2, double sigma2)
    {
        std::array<double, 4> a = grad_sigmoid (x, x1, sigma1);
        std::array<double, 4> b = grad_rsigmoid(x, x2, sigma2);
        std::array<double, 6> aa =
        {{
            a[0] * b[0], a[1] * b[0] + a[0] * b[1],
            a[2] * b[0], a[3] * b[0],
            a[0] * b[2], a[0] * b[3]
        }};
        return aa;
    }

    /// функция двойной сигмодиды (0 - как сумма, 1 - как произведение
    inline double sigmoid2(double x, double x0, double sigma0, double x1, double sigma1, int flag=0)
    { return flag == 0 ? sigmoid2s(x, x0, sigma0, x1, sigma1) : sigmoid2p(x, x0, sigma0, x1, sigma1); }

    /// 2-ная сигмодида и ее производные по параметрам {f, df/dx, df/dx1, df/ds1, df/dx2, df/ds2}
    inline std::array<double, 6> grad_sigmoid2(double x, double x1, double sigma1, double x2, double sigma2, int flag=0)
    {
        return flag == 0
            ? grad_sigmoid2s(x, x1, sigma1, x2, sigma2)
            : grad_sigmoid2p(x, x1, sigma1, x2, sigma2);
    }

}
#endif // PRG_MATH
