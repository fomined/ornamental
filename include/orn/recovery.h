#ifndef ORN_RECOVERY_H
#define ORN_RECOVERY_H

#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <random>

#include "math/statistics.h"
#include "math/math.h"
#include "math/random.h"
#include "orn/orn.defs.h"
#include "orn/time.h"

namespace prg
{

    inline std::pair<double,double> intrainterval_localization_by_directMLM(
        std::vector<rough_value<Time>*> &lb,
        std::vector<rough_value<Time>*> &rb)
    {
        // проходим по интервалу возможных значений с некоторым шагом (прежпочтительней 1)
        // и для каждого значения считаем дисперсию (отклонение от от заданного значения),
        // позиция с минимальной дисперсией означает точку максимальной вероятности

        const double min_phase_duration = 1; // минимально возможная длина фазы

        Unsigned n = lb.size();
        const double dt = 0.1;

        //------------------------------------------------------------------------------------------
        //                      получение границ диапазона
        //------------------------------------------------------------------------------------------

        double minv = inf<double>(), maxv = 0;
        for (Unsigned i=0; i<n; i++)
        {
            auto range = rb[i]->range() - lb[i]->range();
            double min = std::max(range.left().to_double(), min_phase_duration);
            double max = range.right().to_double();

            if (minv > min) minv = min;
            if (maxv < max) maxv = max;
        }

        //------------------------------------------------------------------------------------------
        //              сканирование диапазона с поиском минимума дисперсии
        //------------------------------------------------------------------------------------------
        double minaver = inf<double>(), minsigma = inf<double>();
        for (double v=minv; v<maxv; v+=dt)
        {
            std::vector<double> tmp;
            for (Unsigned i=0; i<n; i++)
            {
                auto range = rb[i]->range() - lb[i]->range();
                double min = std::max(range.left().to_double(), min_phase_duration);
                double max = range.right().to_double();

                if (v < min) tmp.push_back(min);
                else if (v > max) tmp.push_back(max);
                else tmp.push_back(v);
            }
            auto m = statistics(tmp);
            if (minsigma > m[1]) { minaver = m[0]; minsigma = m[1]; }
        }

        //------------------------------------------------------------------------------------------
        //                      установка оптимальных значений
        //------------------------------------------------------------------------------------------
        for (Unsigned i=0; i<n; i++)
        {
            auto range = rb[i]->range() - lb[i]->range();
            double min = std::max(range.left().to_double(), min_phase_duration);
            double max = range.right().to_double();

            double dt = minaver; // устанавливаемое значение длины фазы для конкретного цветка
            if (dt < min) dt = min;
            else if (dt > max) dt = max;

            if (dt < range.minimal_value())
            {
                lb[i]->value() = lb[i]->range().maximal_value();
                rb[i]->value() = rb[i]->range().minimal_value();
            }
            else if (dt > range.maximal_value())
            {
                lb[i]->value() = lb[i]->range().minimal_value();
                rb[i]->value() = rb[i]->range().maximal_value();
            }
            else
            {
                Time s = lb[i]->right(); // предположительный старт
                Time f = s + dt; // предположительный финиш
                if (rb[i]->include(f))
                {
                    // установим на конец стартового интервала неопределенности, это самая вероятная
                    // дата начала фазы, так как в этот день происходит наблюдение
                    lb[i]->value() = s;
                    rb[i]->value() = f;
                }
                else
                {
                    // установим на конец финишного интервала неопределенности, это также вероятная
                    // дата начала фазы (но следующей), других вариантов нет, если не учитывать, что
                    // фаза начинается и кончается внутри интервалов неопределенности
                    rb[i]->value() = rb[i]->right();
                    lb[i]->value() = rb[i]->value() - dt;
                }
            }
        }

        return {minaver, minsigma};
    }

    /*!
     * \brief установить значения rough values, обеспечивающих минимальную дисперсию
     * \param b1,b2,b3 [out] лучшие значения rough values
     * \param l1 заданный интервал b2 - b1
     * \param l2 заданный интервал b3 - b2
     * \return дисперсия
     */
    inline double best_localization_for_directMLM(
        rough_value<Time>* b1, rough_value<Time>* b2, rough_value<Time>* b3,
        double l1, double l2)
    {
        Time min1 = b1->minimal_value(), max1 = b1->maximal_value();
        Time min2 = b2->minimal_value(), max2 = b2->maximal_value();
        Time min3 = b3->minimal_value(), max3 = b3->maximal_value();

        double min_sigma = inf<double>();
        bool flag = false;

        for (Time t2=max2; t2>=min2; t2=t2-1)
        {
            Time maxx1 = std::min(max1, t2-1); // max for 1st phase
            Time minn3 = std::max(min3, t2+1); // min for 3st phase
            for (Time t1=maxx1; t1>=min1; t1=t1-1)
            {
                double s1 = sqr(std::abs(l1 - (t2-t1)));
                for (Time t3=max3; t3>=minn3; t3=t3-1)
                {
                    double sigma = s1 + sqr(std::abs(l2 - (t3-t2)));
                    if (min_sigma > sigma)
                    {
                        min_sigma = sigma;
                        b1->value() = t1;
                        b2->value() = t2;
                        b3->value() = t3;
                        flag = true;
                    }
                }
            }
        }
        if (!flag)
        {
            b2->value() = b3->value() - 1;
            if (b2->value() < min2) b2->value() = min2;

            b1->value() = b2->value() - 1;
            if (b1->value() < min1) b1->value() = min1;

            Time t3 = b3->value();
            Time t2 = b2->value();
            Time t1 = b1->value();

            double s1 = sqr(std::abs(l1 - (t2-t1)));
            min_sigma = s1 + sqr(std::abs(l2 - (t3-t2)));
        }
        return min_sigma;
    }

    inline std::pair<double,double> intrainterval_localization_by_directMLM(
        std::vector<rough_value<Time>*> &b1,
        std::vector<rough_value<Time>*> &b2,
        std::vector<rough_value<Time>*> &b3)
    {
        // проходим по интервалу возможных значений с некоторым шагом (прежпочтительней 1)
        // и для каждого значения считаем дисперсию (отклонение от от заданного значения),
        // позиция с минимальной дисперсией означает точку максимальной вероятности

        const double min_phase_duration = 1; // минимально возможная длина фазы

        Unsigned n = b1.size();
        const double dt = 0.5;

        //------------------------------------------------------------------------------------------
        //                      получение границ диапазона
        //------------------------------------------------------------------------------------------

        double minv1 = inf<double>(), maxv1 = 0;
        double minv2 = inf<double>(), maxv2 = 0;
        for (Unsigned i=0; i<n; i++)
        {
            auto range1 = b2[i]->range() - b1[i]->range();
            double min1 = std::max(range1.left().to_double(), min_phase_duration);
            double max1 = range1.right().to_double();
            if (minv1 > min1) minv1 = min1;
            if (maxv1 < max1) maxv1 = max1;

            auto range2 = b3[i]->range() - b2[i]->range();
            double min2 = std::max(range2.left().to_double(), min_phase_duration);
            double max2 = range2.right().to_double();
            if (minv2 > min2) minv2 = min2;
            if (maxv2 < max2) maxv2 = max2;
        }

        //------------------------------------------------------------------------------------------
        //              сканирование диапазона с поиском минимума дисперсии
        //------------------------------------------------------------------------------------------

        double min_sigma = inf<double>(), sum_sigma = 0;
        double o1 = inf<double>(), o2 = inf<double>();

        for (double v1=minv1; v1<maxv1; v1+=dt)
        {
            for (double v2=minv2; v2<maxv2; v2+=dt)
            {
                sum_sigma = 0;
                for (Unsigned i=0; i<n; i++)
                {
                    sum_sigma += best_localization_for_directMLM(b1[i], b2[i], b3[i], v1, v2);
                }
                if (min_sigma > sum_sigma)
                {
                    min_sigma = sum_sigma;
                    o1 = v1;
                    o2 = v2;
                }
            }
        }

        //------------------------------------------------------------------------------------------
        //                      установка оптимальных значений
        //------------------------------------------------------------------------------------------
        for (Unsigned i=0; i<n; i++)
        {
            best_localization_for_directMLM(b1[i], b2[i], b3[i], o1, o2);
        }

        return {o1, o2};
    }

    /*!
     * \brief установить случайные значения rough values
     * \param b1,b2,b3 [out] значения rough values
     */
    inline void random_localization(rough_value<Time>* b1, rough_value<Time>* b2, rough_value<Time>* b3)
    {
        Time min1 = b1->minimal_value(), max1 = b1->maximal_value();
        Time min2 = b2->minimal_value(), max2 = b2->maximal_value();
        Time min3 = b3->minimal_value(), max3 = b3->maximal_value();

        std::uniform_int_distribution<> dist1(min1.to_signed(), max1.to_signed());
        std::uniform_int_distribution<> dist2(min2.to_signed(), max2.to_signed());
        std::uniform_int_distribution<> dist3(min3.to_signed(), max3.to_signed());

        Unsigned iterations = 100;
        while (iterations-- > 0)
        {
            b1->value() = dist1(randgen());
            b2->value() = dist2(randgen());
            b3->value() = dist3(randgen());

            if (b1->value() < b2->value()
                && b2->value() < b3->value()) break;
        }
    }

    inline void intrainterval_localization_by_random(
        std::vector<rough_value<Time>*> &b1,
        std::vector<rough_value<Time>*> &b2,
        std::vector<rough_value<Time>*> &b3)
    {
        Unsigned n = b1.size();
        for (Unsigned i=0; i<n; i++)
            random_localization(b1[i], b2[i], b3[i]);
    }

}
#endif // ORN_RECOVERY_H
