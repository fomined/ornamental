#ifndef PRG_MATH
#define PRG_MATH

/**
 * @file хранение стандартных сокращений и макросов
 */

#include <cmath>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include "prg.defs.h"

namespace prg
{

    /*!
     * \brief получить факториал $n!$
     * \param x аргумент функции [тип unsigned согласно boost::]
     * \note возвращаемое значение всегда double, иначе переполнение уже при малых n
     */
    inline double factorial(unsigned x) { return boost::math::factorial<double>(x); }

    /*!
     * \brief получить биномиальный коэффициент $C^m_n$
     * \note возвращаемое значение всегда double, иначе переполнение уже при малых n
     */
    inline double binomial(unsigned n, unsigned m)
        { return boost::math::binomial_coefficient<double>(n, m); }


    /*!
     * \brief значение функции Пуассона $e^-m m^k / k!$
     * \param mean среднее значение
     * \param k точка, в которой считается функция
     */
    inline double poisson(double mean, size_t k)
        { return std::exp(-mean) * std::pow(mean, k) / factorial(unsigned(k)); }

}
#endif // PRG_MATH
