/**
* \file _random.h
* \brief Файл содержит генераторы случайных чисел и распределений.
*/
#ifndef RANDOM_H
#define RANDOM_H

#include <type_traits>
#include <ctime>
#include <random>
#include "prg.iterator.h"

/// если RANDGEN не определен, то задается по умалчиванию std::mt19937
#ifndef RANDGEN
    #define RANDGEN std::mt19937
#endif

/*!
 * \brief возвращает генератор случайных чисел (он создается при при первом обращении)
 * \note в _DEBUG режиме инициализация генератора 0, иначе std::time(0)
 * \return ссылка на ГСЧ
 */
inline static RANDGEN &randgen()
{
#ifdef _DEBUG
    static RANDGEN instance(0);
#else
    static RANDGEN instance(std::time(0));
#endif
    return instance;
}

namespace prg
{

	/**
    * \brief Выдача случайного значения в диапазоне [0, max].
    * \param max максимальное значение (включено в диапазон)
	*/
    template <class T, typename std::enable_if<std::is_integral<T>::value>::type * = nullptr>
    inline T ran0(T max)
	{
        std::uniform_int_distribution<T> dist(0, max);
        return dist(randgen());
	}

    /**
    * \brief Выдача случайного значения в диапазоне [min, max].
    * \param min минимальное значение (включено в диапазон)
    * \param max максимальное значение (включено в диапазон)
    */
    template <class T, typename std::enable_if<std::is_integral<T>::value>::type * = nullptr>
    inline T ran0(T min, T max)
    {
        std::uniform_int_distribution<T> dist(min, max);
        return dist(randgen());
    }

    /**
    * \brief Выдача случайного значения в заданном диапазоне.
    * \param max максимальное значение (включено в диапазон)
    * \param min минимальное значение (включено в диапазон)
    */
    template <class T, typename std::enable_if<std::is_floating_point<T>::value>::type * = nullptr>
    inline T ran0(T max)
    {
        std::uniform_real_distribution<T> dist(0., max);
        return dist(randgen());
    }

    /**
    * \brief Выдача случайного значения в заданном диапазоне.
    * \param max максимальное значение (включено в диапазон)
    * \param min минимальное значение (включено в диапазон)
    */
    template <class T, typename std::enable_if<std::is_floating_point<T>::value>::type * = nullptr>
    inline T ran0(T min, T max)
    {
        std::uniform_real_distribution<T> dist(min, max);
        return dist(randgen());
    }


    /*!
     * \brief выбрать число в диапазоне [0, n-1] с вероятностью пропорциональной весу числа
     * \param n число возможных результатов
     * \param w веса результатов
     */
    struct piecewise_constant_distribution : public std::piecewise_constant_distribution<>
    {
        typedef std::piecewise_constant_distribution<> base_t;

        template <typename Container>
        piecewise_constant_distribution(const Container &w)
        : base_t(prg::real_iterator(0.1), prg::real_iterator(w.size()), w.begin()) {}

        piecewise_constant_distribution(size_t n, double w[])
        : base_t(prg::real_iterator(0.1), prg::real_iterator(n), w) {}
    };
}
#endif
