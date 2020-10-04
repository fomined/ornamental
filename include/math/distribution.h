#ifndef PRG_DISTRIBUTION_H
#define PRG_DISTRIBUTION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tuple>
#include <vector>
#include <random>
#include <memory>
#include <algorithm>

#include "defs.h"

namespace  prg
{

    class numeric_distribution
    {
        typedef std::piecewise_linear_distribution<> distribution_type;
        typedef std::shared_ptr<distribution_type> pointer_type;
        typedef distribution_type::result_type result_type;

        pointer_type pd; // указатель на распределение

    public:

        /// пустое распределение
        numeric_distribution() {}

        /// стандартная инициализация распределения
        template <typename Iterator>
        numeric_distribution(Iterator xa, Iterator xb, Iterator w)
        : pd(new distribution_type(xa, xb, w)) {}

        /// стандартная инициализация распределения
        template <typename Iterator>
        void init(Iterator first, Iterator last, Iterator w)
        {
            pd = pointer_type(new distribution_type(first, last, w));
        }

        /*!
         * \brief создание генератора по массиву данных
         * \param first,last итераторы по данным
         * \param op операция, возвращаюшая по данным пару {x, w}
         */
        template <typename Iterator, typename UnaryOperation>
        void init(Iterator first, Iterator last, UnaryOperation op)
        {
            std::vector<double> x, w;
            for (; first!=last; ++first)
            {
                auto v = op(*first);
                x.push_back(v.first);
                w.push_back(v.second);
            }
            pd = pointer_type(new distribution_type(x.begin(), x.end(), w.begin()));
        }

        /*!
         * \brief создание генератора по массиву данных
         * \param first,last итераторы по данным
         * \param op операция, возвращаюшая по данным пару {x, w}
         */
        template <typename Iterator, typename UnaryOperation>
        numeric_distribution(Iterator first, Iterator last, UnaryOperation op)
        { init(first, last, op); }

        /// генерация одного числа по заданному распределению
        template <class Generator>
        result_type operator()(Generator& g) { return (*pd)(g); }

        /// генерация массива чисел по заданному распределению
        template <class Generator>
        std::vector<result_type> operator()(Generator& g, Unsigned n)
        {
            std::vector<result_type> tmp;
            for (Unsigned i=0; i<n; i++) tmp.push_back((*pd)(g));
            return tmp;
        }

        /// генерация массива чисел по заданному распределению с упорядочивание согласно образцу f
        template <class Generator>
        std::vector<result_type> operator()(Generator& g, const std::vector<result_type> &f)
        {
            Unsigned n = f.size();
            //------------------------------------------------------------------------------------
            //          генерация заданного числа значений по имеющемуся распределению
            //------------------------------------------------------------------------------------
            std::vector<result_type> tmp;
            for (Unsigned i=0; i<n; i++) tmp.push_back((*pd)(g));

            //------------------------------------------------------------------------------------
            //                  упорядочение значений согласно образцу
            //------------------------------------------------------------------------------------
            typedef std::tuple<Unsigned, double, double> tuple_type;
            std::vector<tuple_type> ff;

            for (Unsigned i=0; i<n; i++) ff.push_back(tuple_type{i, f[i], 0});
            std::sort(ff.begin(), ff.end(),
                [](const tuple_type &a, const tuple_type &b){ return std::get<1>(a) < std::get<1>(b); }
            ); // упорядочили образец

            std::sort(tmp.begin(), tmp.end()); // упорядочили сгенерированную выборку

            for (Unsigned i=0; i<n; i++) std::get<2>(ff[i]) = tmp[i];
                // сопоставили упорядоченные образец и выборку

            std::sort(ff.begin(), ff.end(),
                [](const tuple_type &a, const tuple_type &b){ return std::get<0>(a) < std::get<0>(b); }
            ); // вернулись к упорядочению образца

            for (Unsigned i=0; i<n; i++) tmp[i] = std::get<2>(ff[i]);
                // перенесли данные, упорядоченные по образцу

            return tmp;
        }
    };


} // namespace prg
#endif // PRG_DISTRIBUTION_H
