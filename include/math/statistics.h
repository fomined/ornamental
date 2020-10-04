#ifndef MATH_STATISTICS_H
#define MATH_STATISTICS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <array>

#include "defs.h"
#include "math/math.h"
#include "math/curve.h"
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/statistics/bivariate_statistics.hpp>
#include <boost/math/statistics/t_test.hpp>
#include <boost/math/statistics/anderson_darling.hpp>
#include <boost/math/statistics/ljung_box.hpp>

namespace prg
{
    namespace bs = boost::math::statistics;

    inline std::string pvalue_as_string(double x)
    {
        if (x < 0.000001)
        {
            char buf[128];
            std::sprintf(buf, "%.2e", x);
            return buf;
        }
        return std::to_string(x);
    }

    inline double normalCDF(double x, double mu, double sigma)
        { return 0.5 * (1 + std::erf((x - mu) * M_SQRT1_2 / sigma)); }

    inline double normalCDF(double x) { return 0.5 * (1 + std::erf((x) * M_SQRT1_2)); }

    inline double average(const std::vector<double> &v) { return bs::mean(v); }
    inline double median(const std::vector<double> &v) { std::vector<double> t(v); return bs::median(t); }
        // vector is modified by bs::median()

    /// returns the population variance ( 1/(n-1) instead 1/n)
    inline double variance(const std::vector<double> &v) { return bs::variance(v); }
    inline double sample_variance(const std::vector<double> &v) { return bs::sample_variance(v); }
    inline double standard_derivation(const std::vector<double> &v) { return std::sqrt(variance(v)); }

    inline double skewness(const std::vector<double> &v) { return bs::skewness(v); }
    inline double kurtosis(const std::vector<double> &v) { return bs::kurtosis(v); }

    inline double covariance(const std::vector<double> &u, const std::vector<double> &v)
        { return bs::covariance(u, v); }

    /// returns Pearson correlation coefficient
    inline double correlation_coefficient(const std::vector<double> &u, const std::vector<double> &v)
        { return bs::correlation_coefficient(u, v); }

    /// returns [test statistic, p-value]
    inline auto one_sample_t_test(double sample_mean, double sample_variance, double num_samples, double assumed_mean)
        { return bs::one_sample_t_test(sample_mean, sample_variance, num_samples, assumed_mean); }

    const Unsigned anderson_darling_boost = 0; // from BOOST library
    const Unsigned anderson_darling_sym = 1; // from PAST tool tutorial (!) DON'T WORK
    const Unsigned anderson_darling_asym = 2; // PAST tool tutorial (alternative) = boost

    template <Unsigned key>
    double anderson_darling_normality_statistic(const std::vector<double> &, double, double);

    template <>
    inline double anderson_darling_normality_statistic<anderson_darling_boost>(const std::vector<double> &v,
        double mean, double stddev)
    {
        assert(!v.empty());
        std::vector<double> t(v);
        std::sort(t.begin(), t.end());

        return bs::anderson_darling_normality_statistic(t, mean, stddev);
    }

    template <>
    inline double anderson_darling_normality_statistic<anderson_darling_sym>(const std::vector<double> &v,
        double mean, double stddev)
    {
        assert(!v.empty());
        std::vector<double> t(v);
        std::sort(t.begin(), t.end());

        Unsigned n = t.size();
        double A2 = 0;
        for (Unsigned i=0; i<n; i++)
        {
            double Yi = (t[i] - mean ) / stddev;
            double Yni = (t[n-i] - mean ) / stddev;
            double fi = normalCDF(Yi);
            double fni = normalCDF(Yni);
            A2 -= (2*i + 1) * (std::log(fi) + std::log(1 - fni));
        }
        return A2/n - n;
    }

    template <>
    inline double anderson_darling_normality_statistic<anderson_darling_asym>(const std::vector<double> &v,
        double mean, double stddev)
    {
        assert(!v.empty());
        std::vector<double> t(v);
        std::sort(t.begin(), t.end());

        Unsigned n = t.size();
        double A2 = 0;
        for (Unsigned i=0; i<n; i++)
        {
            double Yi = (t[i] - mean ) / stddev;
            double fi = normalCDF(Yi);
            A2 -= (2*i + 1) * std::log(fi) + (2 * (n-i) - 1) * std::log(1 - fi);
        }
        return A2/n - n;
    }

    /// value of anderson_darling needs to be adjusted for small sample sizes
    inline auto adjusted_anderson_darling_statistic(double A2, Unsigned n)
        { return A2 * (1 + 0.75 / n + 2.25 / sqr(n)); }

    /// pvalue from adjusted_anderson_darling_normality_statistic
    inline double stephens_pvalue(double aad)
    {
        double pvalue = 0;
        if (aad >= 0.6 ) pvalue = std::exp(1.2937 - 5.709 * aad + 0.0186 * sqr(aad));
        else if (aad >= 0.34) pvalue = std::exp(0.9177 - 4.279 * aad - 1.38 * sqr(aad));
        else if (aad >= 0.2 ) pvalue = 1 - std::exp(-8.318 + 42.796 * aad - 59.938 * sqr(aad));
        else pvalue = 1 - std::exp(-13.436 + 101.14 *aad - 223.73 *sqr(aad));
        return pvalue;
    }

    template <Unsigned n>
    inline auto anderson_darling_normality_test(const std::vector<double> &v,
        double mean, double stddev)
    {
        double A2 = anderson_darling_normality_statistic<n>(v, mean, stddev);
        A2 = adjusted_anderson_darling_statistic(A2, v.size());
        double pvalue = stephens_pvalue(A2);
        return std::make_pair(A2, pvalue);
    }

    /*!
     * \brief ljung_box test
     * \note The Ljung-Box test is used to test if residuals from a fitted model have unwanted
     *      autocorrelation. If autocorrelation exists in the residuals, then presumably a model
     *      with more parameters can be fitted to the original data and explain more of
     *      the structure it contains. H0 = "residuals -> normal distrubution"
     * \param v
     * \param lags
     * \param fit_dof
     * \return [test statistic, p-value]
     */
    inline auto ljung_box(const std::vector<double> &v, int64_t lags = -1, int64_t fit_dof = 0)
        { return bs::ljung_box(v, lags, fit_dof); }

    template <Unsigned N1, Unsigned N2>
    inline double determination_coef(shared<curve<N1>> fit, shared<curve<N2>> data)
    {
        std::sort(fit->begin(), fit->end(), [](auto a, auto b) { return a[0] < b[0]; });
        std::sort(data->begin(), data->end(), [](auto a, auto b) { return a[0] < b[0]; });

        assert(fit->front()[0] <= data->front()[0]);
        assert(fit->back()[0] >= data->back()[0]);

        double ya = 0; Unsigned n = 0;
        for (auto e : *data) { ya += e[1]; n++; } ya /=n;
            // среднее значение функции на интервале

        double dispE = 0; // дисперсия полная
        for (auto e : *data) { dispE += sqr(e[1] - ya); }
        if (dispE == 0) return 1;

        double dispT = 0; // дисперсия необъясненная
        for (auto e : *data)
        {
            double x = e[0];
            double y = e[1];
            auto it = std::lower_bound(fit->begin(), fit->end(), x,
                [](auto a, double x){ return a[0] < x; }
            );
            if (it == fit->end()) continue;

            double y1 = (*it)[1];
            dispT += sqr(y - y1);

        }
        return 1 - dispT/ dispE;
    }

    enum test_enum {ONE_SIDED_TEST, TWO_SIDED_TEST};

    /// получение одно/двух стороннего p-значения
    inline double pvalue(double x, double mean, double sigma, test_enum e=TWO_SIDED_TEST)
    {
        double p = 1. - std::erf((x - mean)/std::sqrt(2*sqr(sigma)));
        return e == ONE_SIDED_TEST ? 0.5 * p : p;
    }

    /**
     * @brief t_statistic используется для проверки значимости линейного коэффициента регрессии r(X,Y)
     * \note H0 коэффициент корреляции незначим, и не отличается от нуля.
     * @param r sample correlation coefficient
     * @param nfreedom число степеней свободы
     * @return t_statistic
     */
    inline double t_statistic(double r, Signed nfreedom)
    {
        assert(nfreedom > 2);
        return r * (nfreedom - 2.) / std::sqrt(1. - r*r);
    }

    inline auto correlation_test(double rxy, Signed nfreedom)
    {
        double t = t_statistic(rxy, nfreedom);
        double pv = pvalue(t, 0, 1, TWO_SIDED_TEST);
        return std::make_pair(t, pv);
    }

    /// моменты распределения случайной величины ([0] среднее, [1] дисперсия, ..)
    struct moment_t
    {
        enum { n = 2 };
        std::array<double, n> a;

        moment_t(double x=0, double s=0) { a[0] = x; a[1] = s; }
        const double &operator[](size_t i) const { assert(i<n); return a[i]; }
        double &operator[](size_t i) { assert(i<n); return a[i]; }
    };

    inline std::ostream &operator<<(std::ostream &os, const moment_t &m)
        { os << "{average, sigma} = {" << m[0] << "," << m[1] << "}\n"; return os; }

    /*!
     * \brief получение статистических центральных моментов (среднего и дисперсии) для массива
     */
    inline moment_t statistics(const std::vector<double> &v)
    {
        size_t n = v.size();

        double av = 0;
        for (size_t i=0; i<n; i++) av += v[i];
        av /= n;

        double s = 0;
        for (size_t i=0; i<n; i++) s += sqr(v[i] - av);
        s = sqrt(s / n);

        return moment_t(av, s);
    }

    /*!
     * \brief получение статистических центральных моментов (среднего и дисперсии) для массива
     * \param w веса элементов
     */
    inline moment_t statistics(const std::vector<double> &v, const std::vector<double> &w)
    {
        size_t cnt = v.size();

        double av = 0, sumw = 0;
        for (size_t i=0; i<cnt; i++) { av += w[i] * v[i]; sumw += w[i]; }
        assert(sumw > 0);
        av /= sumw;

        double s = 0;
        for (size_t i=0; i<cnt; i++) { s += w[i] * sqr(v[i] - av); }
        s = sqrt(s / sumw);

        return moment_t(av, s);
    }

    class Histogram : protected std::vector<Unsigned>
    {
        typedef std::vector<Unsigned> base_t;
        double a, b, dx; // границы и ширина корзины

        flagged_value<Unsigned> index(double v) const
        {
            if (a <= v && v < b)
            {
                Unsigned k = Unsigned((v - a) / dx);
                assert(k < size());
                return {k, true};
            }
            return {0, false};
        }

        double value(Unsigned k) const { return k * dx; }

    public :

        Histogram(double a, double b, double dx)
        : base_t(Unsigned((b - a) / dx) + 1, 0), a(a), b(b), dx(dx) {}

        Histogram(double a, double b, Unsigned n)
        : base_t(n, 0), a(a), b(b), dx((b - a)/n) {}

        /// вставка только тех значений, которые находятся в диапазоне [a, b)
        bool insert(double v) { auto f = index(v); if (f.second) at(f.first)++; return f.second; }

        void print()
        {
            std::printf("--------- histogram start ---------\n");
            for (Unsigned i=0; i<size(); i++)
            {
                std::printf("%lf %ld\n", value(i), at(i));
            }
            std::printf("--------- histogram end ---------\n");
            FLUSH
        }
    };
}
#endif // MATH_STATISTICS_H








//#ifdef BOOST_1_74
//#else

//    /// получение среднего для массива
//    inline double average(const std::vector<double> &v)
//    {
//        size_t n = v.size();
//        double sum = 0; for (size_t i=0; i<n; i++) sum += v[i];
//        return sum / n;
//    }

//    /// получение медианного значения для массива
//    inline double median(const std::vector<double> &v)
//    {
//        std::vector<double> vv(v);
//        std::sort(vv.begin(), vv.end());

//        Unsigned n = v.size();
//        if (n % 2 == 1) return vv[n / 2 + 1];

//        return vv[n/2 - 1] + vv[n/2];
//    }

//    /// получение дисперсии для массива
//    /// !делать по умалчиванию mean=0 нельзя, так как провоцирует ошибки забывания
//    inline double variance(const std::vector<double> &v, double mean)
//    {
//        size_t n = v.size();
//        double sum = 0; for (size_t i=0; i<n; i++) sum += sqr(v[i] - mean);
//        return sum / n;
//    }

//    enum test_enum {ONE_SIDED_TEST, TWO_SIDED_TEST};

//    /// получение одно/двух стороннего p-значения
//    inline double pvalue(double x, double mean, double sigma, test_enum e=TWO_SIDED_TEST)
//    {
//        double p = 1. - std::erf((x - mean)/std::sqrt(2*sqr(sigma)));
//        return e == ONE_SIDED_TEST ? 0.5 * p : p;
//    }

//    /// моменты распределения случайной величины ([0] среднее, [1] дисперсия, ..)
//    struct moment_t
//    {
//        enum { n = 2 };
//        std::array<double, n> a;

//        moment_t(double x=0, double s=0) { a[0] = x; a[1] = s; }
//        const double &operator[](size_t i) const { assert(i<n); return a[i]; }
//        double &operator[](size_t i) { assert(i<n); return a[i]; }
//    };

//    inline std::ostream &operator<<(std::ostream &os, const moment_t &m)
//        { os << "{average, sigma} = {" << m[0] << "," << m[1] << "}\n"; return os; }

//    /*!
//     * \brief получение статистических центральных моментов (среднего и дисперсии) для массива
//     */
//    inline moment_t statistics(const std::vector<double> &v)
//    {
//        size_t n = v.size();

//        double av = 0;
//        for (size_t i=0; i<n; i++) av += v[i];
//        av /= n;

//        double s = 0;
//        for (size_t i=0; i<n; i++) s += sqr(v[i] - av);
//        s = sqrt(s / n);

//        return moment_t(av, s);
//    }

//    /*!
//     * \brief получение статистических центральных моментов (среднего и дисперсии) для массива
//     * \param w веса элементов
//     */
//    inline moment_t statistics(const std::vector<double> &v, const std::vector<double> &w)
//    {
//        size_t cnt = v.size();

//        double av = 0, sumw = 0;
//        for (size_t i=0; i<cnt; i++) { av += w[i] * v[i]; sumw += w[i]; }
//        assert(sumw > 0);
//        av /= sumw;

//        double s = 0;
//        for (size_t i=0; i<cnt; i++) { s += w[i] * sqr(v[i] - av); }
//        s = sqrt(s / sumw);

//        return moment_t(av, s);
//    }

//    class Histogram : protected std::vector<Unsigned>
//    {
//        typedef std::vector<Unsigned> base_t;
//        double a, b, dx; // границы и ширина корзины

//        flagged_value<Unsigned> index(double v) const
//        {
//            if (a <= v && v < b)
//            {
//                Unsigned k = Unsigned((v - a) / dx);
//                assert(k < size());
//                return {k, true};
//            }
//            return {0, false};
//        }

//        double value(Unsigned k) const { return k * dx; }

//    public :

//        Histogram(double a, double b, double dx)
//        : base_t(Unsigned((b - a) / dx) + 1, 0), a(a), b(b), dx(dx) {}

//        Histogram(double a, double b, Unsigned n)
//        : base_t(n, 0), a(a), b(b), dx((b - a)/n) {}

//        /// вставка только тех значений, которые находятся в диапазоне [a, b)
//        bool insert(double v) { auto f = index(v); if (f.second) at(f.first)++; return f.second; }

//        void print()
//        {
//            std::printf("--------- histogram start ---------\n");
//            for (Unsigned i=0; i<size(); i++)
//            {
//                std::printf("%lf %ld\n", value(i), at(i));
//            }
//            std::printf("--------- histogram end ---------\n");
//            FLUSH
//        }
//    };
//#endif
