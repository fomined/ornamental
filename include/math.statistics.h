#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <array>

#include "defs.h"
#include "math.defs.h"

/// получение среднего для массива
inline double average(const std::vector<double> &v)
{
    size_t n = v.size();
    double sum = 0; for (size_t i=0; i<n; i++) sum += v[i];
    return sum / n;
}

/// получение дисперсии для массива
/// !делать по умалчиванию mean=0 нельзя, так как провоцирует ошибки забывания
inline double variance(const std::vector<double> &v, double mean)
{
    size_t n = v.size();
    double sum = 0; for (size_t i=0; i<n; i++) sum += sqr(v[i] - mean);
    return sum / n;
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

#endif // STATISTICS_H
