#ifndef RANDOM_H
#define RANDOM_H

/**
* \file _random.h
* \brief Файл содержит генераторы случайных чисел и распределений.
*/

#include <type_traits>
#include <ctime>
#include <random>
#include "defs.h"

/*!
 * \brief возвращает генератор случайных чисел (он создается при при первом обращении)
 * \note инициализация генератора всегда std::time(0), если необходимо генерировать
 *  одинаковые последовательности для отладки, то необходимо использовать вызов ranseed(int)
 * \return ссылка на ГСЧ
 */
inline static std::mt19937_64 &randgen()
{
    static std::random_device randev;
    static std::mt19937_64 tmp(randev());
    return tmp;
}

/*!
 * \brief сгенерировать случайое число в качестве seed для генератора
 * \note используется, чтобы знать, какое seed использовать 
 *  для воспроизведения того же поведения при отладке
 */
inline Unsigned get_ranseed()
{
    static std::uniform_int_distribution<Unsigned> dist(0, std::numeric_limits<Unsigned>::max());
    auto tmp = dist(randgen());
    std::printf("genered seed = %ld\n", tmp); FLUSH;
    return tmp;
}

#define rangen randgen()

/// переустановка генератора случайных чисел в DEBUG режим
#define set_ranseed(v) randgen().seed(v)

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
        return dist(rangen);
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
        return dist(rangen);
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
        return dist(rangen);
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
        return dist(rangen);
    }
}
#endif
