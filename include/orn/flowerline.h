#ifndef FLOWERLINE_H
#define FLOWERLINE_H

#include <string>
#include <vector>
#include <string>
#include <set>
#include <tuple>
#include <memory>

#include "container/vector.h"
#include "math/math.h"
#include "math/statistics.h"

#include "os/os.h"

#include "orn/orn.defs.h"
#include "orn/time.h"
#include "orn/phase.h"
#include "orn/flower.h"

/*!
 * \file типа хранения записи flowerline и основные функторы его обработки
 */
namespace prg
{

    /// основной тип данных хранения начала и конца фазы для отдельного цветка
    class flowerline
    {
        rough_value<Time> start_interval_[PHASE_CNT]; // начало каждой фазы развития (b, f, o, s, d)

    public:

        FlowerIndex fsp;

        /// тип записи
        bool dynamic_flag() const { return fsp.flag(); }
        bool empty(Phase e) const { return start_interval_[to_unsigned(e)].range().empty(); }

        /// ВНИМАНИЕ для MLM необходим возврат указателя (модификация на месте, вместо перекопирования)
        /// доступ по значению enum, чтобы избежать явного его преобразования к индексу
        /// возврат указателя, чтобы обеспечить прямой доступ на модификацию значения
        /// (возврат ссылки провоцирует ошибки, так как для константного объекта происходит
        /// создание временного объекта)
        const rough_value<Time> *start_interval(Phase e) const { return &start_interval_[to_unsigned(e)]; }
        rough_value<Time> *start_interval(Phase e) { return &start_interval_[to_unsigned(e)]; }

        /*!
         * \brief расчет длины фазы (по "центральным" от грубых значений)
         * \return false, если какая либо фаза не установлена
         */
        flagged_value<double> duration_of_phase(Phase e) const
        {
            if (dynamic_flag())
            {
                auto t1 = start_interval(e)->value();
                auto t2 = start_interval(e+1)->value();
                return {t2 - t1, t1.valid() && t2.valid()};
            }
            return {0, false};
        }

        /// доступ к "центральному значению"
        Time value(Phase e) const { assert(dynamic_flag()); return start_interval(e)->value(); }
        Time &value(Phase e) { assert(dynamic_flag()); return start_interval(e)->value(); }

        /// получение фазы на заданную дату
        Phase phase(Time t)
        {
            assert(dynamic_flag());
            if (t < value(BUD)) return UNDEF;
            if (t < value(FLOWER)) return BUD;
            if (t < value(OVARY)) return FLOWER;
            if (t < value(SEED)) return OVARY;
            if (t < value(DEATH)) return SEED;
            return DEATH;
        }

        /// тестовая печать объекта
        std::string print() const
        {
            std::string msg;

            std::string fsps = fsp.to_string();
            msg += to_string("%8s : ", fsps.c_str());

            if (!dynamic_flag()) return msg;

            for (auto e : { BUD, FLOWER, OVARY, SEED, DEATH })
            {
                if (!empty(e))
                {
                    auto v = start_interval_[to_unsigned(e)];
                    Time x = v.value();
                    Unsigned x1 = v.range().left().to_month_day();
                    Unsigned x2 = v.range().right().to_month_day();

                    std::string s1 = std::to_string(x1);
                    if (s1.size() == 1) s1 = " " + s1;
                    std::string s2 = std::to_string(x2);
                    if (s2.size() == 1) s2 = " " + s2;

                    msg += x.to_string() + '(' + s1 + ".." + s2 + "]  ";
                }
            }
            return msg;
        }
    };

    inline bool operator< (const flowerline &a, const flowerline &b) { return a.fsp < b.fsp; }
    inline std::ostream &operator<< (std::ostream &os, const flowerline &f) { os << f.print(); return os; }


    //----------------------------------------------------------------------------------------------
    //                          обработка входных данных
    //----------------------------------------------------------------------------------------------
    /// число цветков на разных осях (и суммарно)
    class flower_count_functor
    {
        Unsigned axis;
    public :
        enum { n = 2 }; // число возвращаемых значений
        typedef flagged_value<curve<n>::value_type> return_type; // возвращаемый тип данных

        /*!
         * \brief инициализация функтора
         * \param axis ось, для которой считается число цветков (счет с 1)
         * \note при axis=0 считается общее число цветков на всех осях
         */
        flower_count_functor(Unsigned axis) : axis(axis) {}

        /// операция расчета
        return_type operator() (const flowerline &e)
        {
            return_type r({}, false);
            if (e.dynamic_flag()) return r;

            auto n = e.fsp.axes();

            if (axis == 0) r.first[0] = double(n[0] + n[1] + n[2] + n[3]);
            else r.first[0] = double(n[axis-1]);
            r.first[1] = 1;
            r.second = true;
            if (r.first[0] > 300)
            {
std::cout << r.first[0] << "\n";
std::cout.flush();
            }
            return r;
        }
    };

    class axes2_count_functor
    {
    public :
        enum { n = 2 }; // число возвращаемых значений
        typedef flagged_value<curve<n>::value_type> return_type; // возвращаемый тип данных

        axes2_count_functor() {}

        /*!
         * \brief выдача событий (наличия оси заданного порядка у побега)
         * \return {sprout, 1}, {sprout, 1}, {sprout, 1} ..
         *  Сумма событий дает число осей заданного порядка у побега
         */
        return_type operator() (const flowerline &e)
        {
            return_type empty({}, false);
            if (!e.dynamic_flag()) return empty;

            FlowerIndex f = e.fsp;
            if (f.axis() != 2 || f.position() != 1) return empty;

            double sprout = f.sprout();
            return {{ sprout, 1.}, true};
        }
    };

    /// отношение числа цветков между заданными осями
    class flower_rate_functor
    {
        Unsigned xaxis; // независимая ось
        Unsigned yaxis; // зависимая ось
    public :
        enum { n = 2 }; // число возвращаемых значений
        typedef flagged_value<curve<n>::value_type> return_type; // возвращаемый тип данных

        /*!
         * \brief инициализация функтора
         * \n note обрабатывает только структурные записи
         * \param xaxis независимая ось, для которой считается число цветков (счет с 1)
         * \param yaxis зависимая ось, для которой считается число цветков (счет с 1)
         * \note при yaxis=0 считается общее число цветков на всех осях, больших xaxis
         */
        flower_rate_functor(Unsigned xaxis, Unsigned yaxis) : xaxis(xaxis), yaxis(yaxis) {}

        /// операция расчета
        return_type operator() (const flowerline &e)
        {
            return_type r({}, false);
            if (e.dynamic_flag()) return r;

            auto n = e.fsp.axes();
            r.first[0] = n[xaxis-1];

            if (yaxis == 0) for (Unsigned i=xaxis; i<4; i++) r.first[1] += n[i];
            else r.first[1] = n[yaxis-1];
            r.second = true;

            return r;
        }
    };

    class axis_density_functor
    {
        Unsigned axis; // ось цветения
    public :
        enum { n = 2 }; // число возвращаемых значений
        typedef flagged_value<curve<n>::value_type> return_type; // возвращаемый тип данных

        /*!
         * \brief инициализация функтора
         * \n note обрабатывает только записи динамики
         * \param axis ось цветения (счет с 1)
         */
        axis_density_functor(Unsigned axis) : axis(axis)  {}

        /// операция расчета
        return_type operator() (const flowerline &e)
        {
            return_type r({}, false);
            if (!e.dynamic_flag()) return r;

            FlowerIndex f = e.fsp;
            if (f.axis() != axis) return r;

            double x = f[1]; // позиция по высоте
            r.first = {x, 1};
            r.second = true;

            return r;
        }
    };

    /// длина любой фазы развития
    class phase_duration_functor
    {
        Phase phase; // идентификатор фазы
    public :
        enum { n = 2 }; // число возвращаемых значений
        typedef flagged_value<curve<n>::value_type> return_type; // возвращаемый тип данных

        phase_duration_functor(Phase phase) : phase(phase) {}
        return_type operator() (const flowerline &e)
        {
            return_type r({}, false);
            if (!e.dynamic_flag()) return r;

            auto dt = e.duration_of_phase(phase);
            r.first[0] = dt.value();
            r.first[1] = 1;
            r.second = dt.flag();
            return r;
        }
    };

    const Unsigned PHASE_START = 0;
    const Unsigned PHASE_MID = 1;
    const Unsigned PHASE_END = 2;

    /// вычленение начала/конца фазы
    class axis_phase_functor
    {
        Unsigned axis; // ось цветения
        Phase phase; // идентификатор фазы
        Unsigned flag; // начало[0], середина[1], конец[2] фазы
    public :
        enum { n = 2 }; // число возвращаемых значений
        typedef flagged_value<curve<n>::value_type> return_type; // возвращаемый тип данных

        /*!
         * \brief инициализация функтора
         * \n note обрабатывает только записи динамики
         * \param axis ось цветения (счет с 1)
         * \param phase фаза развития
         * \param flag (PHASE_START, PHASE_MID, PHASE_END)
         */
        axis_phase_functor(Unsigned axis, Phase phase, Unsigned flag=PHASE_START)
        : axis(axis), phase(phase), flag(flag) {}

        /*!
         * \brief операция извлечения данных
         * \return {h, time}
         */
        return_type operator() (const flowerline &e)
        {
            return_type r({}, false);
            if (!e.dynamic_flag()) return r;

            FlowerIndex f = e.fsp;
            if (f.axis() != axis) return r;

            double x = f[1]; // позиция по высоте
            double y = 0;

            if (flag == PHASE_START)
            {
                Time t = e.value(phase);
                y = t.to_double();
                r.second = t.valid();
            }
            if (flag == PHASE_END)
            {
                Time t = e.value(phase + 1);
                    // конец фазы = начало следующей
                y = t.to_double();
                r.second = t.valid();
            }
            if (flag == PHASE_MID)
            {
                Time t1 = e.value(phase);
                Time t2 = e.value(phase+1) - 1;
                y = 0.5 * (t2.to_double() + t1.to_double());
                r.second = t1.valid() && t2.valid();
            }
            r.first = {x, y};
            return r;
        }
    };

    /// вычленение и начала и конца фазы относительно заданной даты
    class extraction_phase_of_axis
    {
        Unsigned axis; // ось цветения
        Phase phase; // идентификатор фазы
        Time start; // starting point
    public :
        enum { n = 3 }; // число возвращаемых значений
        typedef flagged_value<curve<n>::value_type> return_type; // возвращаемый тип данных

        /*!
         * \brief инициализация функтора
         * \note обрабатывает только записи динамики
         * \param axis ось цветения (счет с 1)
         * \param phase фаза развития
         */
        extraction_phase_of_axis(Unsigned axis, Phase phase, Time start=0.)
        : axis(axis), phase(phase), start(start) {}

        /*!
         * \brief операция извлечения данных
         * \return {h, time_start, time_finish}
         */
        return_type operator() (const flowerline &e)
        {
            return_type r({}, false);
            if (!e.dynamic_flag()) return r;

            FlowerIndex f = e.fsp;
            if (f.axis() != axis) return r;

            double x = f[1]; // позиция по высоте
            Time t1 = e.value(phase);
            Time t2 = e.value(phase + 1) - 1;
                // конец фазы на день ранее, чем начало следующей

            r.second = t1.valid() && t2.valid();
                // после -start могут стать недействительными

            t1 = t1 - start;
            t2 = t2 - start;

            r.first = {x, t1.to_double(), t2.to_double()};

            return r;
        }
    };

    //----------------------------------------------------------------------------------------------
    //                          обработка выходных данных
    //----------------------------------------------------------------------------------------------
//    const WLSM_item default_WLSM_options =
//    {
//        2, // мощность полинома
//        40, // число итераций подгонки полинома
//        {0.01, true}, // шаг сглаживания
//        {0.05, true}, // окно сглаживания
//        {inf<double>(),-inf<double>(0)}, // интервал сглаживания
//        {0.1, true}, // отступ слева при отсутствии задания интервала
//        {0.1, true}, // отступ справа при отсутствии задания интервала
//        false // флаг печати процесса подгонки
//    };

//    /// вызов объекта, который выполняет сглаживание
//    inline shared<curve<2>> smooth(Unsigned k, shared<curve<2>> p, const WLSM_item &item=WLSM_item())
//        { return WLSM().smooth(k, p, item); }
//    inline shared<curve<3>> smooth_with_errors(Unsigned k, shared<curve<2>> p, const WLSM_item &item=WLSM_item())
//        { return WLSM().smooth_with_errors(k, p, item); }


    template <Unsigned, typename T>
    inline shared<curve<2>> transform(shared<curve<2>> p, T param);

    template <Unsigned, typename T>
    inline shared<curve<2>> transform(shared<curve<1>> p, T param);

    constexpr Unsigned MIN_BOUND    = 1;
    constexpr Unsigned MAX_BOUND    = 2;
    constexpr Unsigned UNITE        = 3; // (, double) точность объединения
    constexpr Unsigned ACCUMULATE   = 4; // (, double) не используется
    constexpr Unsigned FILL_GAPS    = 5; // (, std::array<double, 3>) шаг, левый и правый отступы (ОТНОСИТЕЛЬНЫЕ!)
    constexpr Unsigned TO_HEIGHTS   = 6;
    constexpr Unsigned TO_START_DAY = 7;
    constexpr Unsigned EXCLUDE_X    = 8;
    constexpr Unsigned SORT         = 9; // (, int) номер колонки
    constexpr Unsigned ABOVE        = 10;
    constexpr Unsigned TO_XRANGE    = 11;
    constexpr Unsigned TO_YRANGE    = 12;
    constexpr Unsigned SMOOTH       = 13;
    constexpr Unsigned SMOOTH0      = 14;
    constexpr Unsigned SMOOTH1      = 15;
    constexpr Unsigned TO_HISTOGRAM = 16;

    template <>
    inline shared<curve<2>> transform<MIN_BOUND>(shared<curve<2>> f, double eps)
    {
        double len = f->back()[0] - f->front()[0];

        eps *= len;
        std::sort(f->begin(), f->end(),
            [eps](auto a, auto b)
            {
                if (equal(a[0], b[0], eps)) return a[1] < b[1];
                return a[0] < b[0];
            }
        );
        auto last = std::unique(f->begin(), f->end(),
            [eps](auto a, auto b) { return equal(a[0], b[0], eps); }
        );
        f->erase(last, f->end());
        return f;
    }

    template <>
    inline shared<curve<2>> transform<MAX_BOUND>(shared<curve<2>> f, double eps)
    {
        double len = f->back()[0] - f->front()[0];

        eps *= len;
        std::sort(f->begin(), f->end(),
            [eps](auto a, auto b)
            {
                if (equal(a[0], b[0], eps)) return a[1] > b[1];
                return a[0] < b[0];
            }
        );
        auto last = std::unique(f->begin(), f->end(),
            [eps](auto a, auto b) { return equal(a[0], b[0], eps); }
        );
        f->erase(last, f->end());
        return f;
    }

    template <>
    inline shared<curve<2>> transform<UNITE>(shared<curve<2>> f, double eps)
    {
        double len = f->back()[0] - f->front()[0];
        eps *= len;

        Unsigned iprev = 0; double vprev = (*f)[0][0];
        for (Unsigned i=1; i<f->size(); i++)
        {
            if (equal((*f)[i][0], vprev, eps))
            {
                (*f)[iprev][1] += (*f)[i][1];
                (*f)[i][1] = 0;
            }
            else
            {
                vprev = (*f)[i][0];
                iprev = i;
            }
        }
        auto last = std::stable_partition(f->begin(), f->end(), [](auto a){ return a[1] > 0; });
        f->erase(last, f->end());
        return f;
    }

    template <>
    inline shared<curve<2>> transform<ACCUMULATE>(shared<curve<2>> f, double)
    {
        for (Unsigned i=1; i<f->size(); i++) (*f)[i][1] += (*f)[i-1][1];
        double norm_coef = 1. / f->back()[1];
        for (auto &e : *f) e[1] *= norm_coef;
        return f;
    }

    template <>
    inline shared<curve<2>> transform<FILL_GAPS>(shared<curve<2>> f, std::array<double, 3> param)
    {
        double len = f->back()[0] - f->front()[0];
        assert(len > 0);
        param[0] *= len; // hx
        param[1] *= len; // left indent
        param[2] *= len; // right indent

        double xa = f->front()[0] - param[1];
        double xb = f->back()[0] + param[2];

        shared<curve<2>> tmp(new curve<2>);
        Unsigned i = 0, n = f->size(); double prev = 0;
        for (double x=xa; x<=xb; x+=param[0])
        {
            while (i < n && (*f)[i][0] <= x) { prev = (*f)[i][1]; i++; }
            tmp->push_back({x, prev});
        }
        if (i != n) tmp->push_back({xb, f->back()[1]});

        //f = std::move(tmp);
        f->clear();
        for (auto e : *tmp) f->push_back(e);

        return f;
    }

    template <>
    inline shared<curve<2>> transform<TO_HEIGHTS>(shared<curve<2>> f, Unsigned n)
    {
        for (auto &e : *f) e[0] = relative_height(e[0], double(n));
        std::sort(f->begin(), f->end(), [](auto a, auto b) { return a[0] < b[0]; });
            // упорядочим по высоте для сглаживания
        return f;
    }

    template <>
    inline shared<curve<2>> transform<TO_START_DAY>(shared<curve<2>> f, Time day)
    {
        double t = day.to_double();
        for (auto &e : *f) e[1] = e[1] - t;
        return f;
    }

    template <>
    inline shared<curve<2>> transform<EXCLUDE_X>(shared<curve<2>> f, int)
    {
        for (auto &e : *f) { e[0] = e[1]; e[1] = 1; }
        return f;
    }

    template <>
    inline shared<curve<2>> transform<SORT>(shared<curve<2>> f, int i)
    {
        std::sort(f->begin(), f->end(), [i](auto a, auto b){ return a[i] < b[i]; });
        return f;
    }

    template <>
    inline shared<curve<2>> transform<ABOVE>(shared<curve<2>> f, double h)
    {
        for (auto &e : *f) { if (e[1] < h) e[1] = h; }
        return f;
    }

    template <>
    inline shared<curve<2>> transform<TO_XRANGE>(shared<curve<2>> f, Interval<double, 2> range)
    {
        double xa = range.minimal_value();
        double xb = range.maximal_value();

        auto last = std::stable_partition(f->begin(), f->end(), [xa, xb](auto a){ return xa <= a[0] && a[0] <=xb; });
        f->erase(last, f->end());

        return f;
    }

    template <>
    inline shared<curve<2>> transform<TO_YRANGE>(shared<curve<2>> f, Interval<double, 2> range)
    {
        double xa = range.minimal_value();
        double xb = range.maximal_value();

        auto last = std::stable_partition(f->begin(), f->end(), [xa, xb](auto a){ return xa <= a[1] && a[1] <=xb; });
        f->erase(last, f->end());

        return f;
    }

    template <>
    inline shared<curve<2>> transform<SMOOTH>(shared<curve<2>> f, int k)
    {
        auto g = smooth(k, f); // обязательное присвоение, т.к. smooth не меняет аргумент
        f->swap(*g);
        return f;
    }

    template <>
    inline shared<curve<2>> transform<TO_HISTOGRAM>(shared<curve<1>> f, int)
    {
        shared<curve<2>> tmp(new curve<2>);
        curve<2>::value_type v;
        for (auto e : *f) tmp->push_back({e[0], 1.});
        transform<SORT>(tmp, 0);
        transform<UNITE>(tmp, 0.01);
        transform<ACCUMULATE>(tmp, 0.);
        transform<FILL_GAPS>(tmp, std::array<double, 3>{0.025, 0.1, 0.1});
        transform<SMOOTH>(tmp, 1);
        transform<ABOVE>(tmp, 0.);
        return tmp;
    }

//    //---------------------------------------------------------------------------------------------
//    //                              micro-functions
//    //---------------------------------------------------------------------------------------------
//    /*!
//     * \brief extract lower|upper bound of ordered curve
//     * \note область определения делится на сегменты длиной hx и в каждом сегменте выбирается
//     *  верхний/нижний представитель
//     * \param f функция
//     * \param lower [true] lower bound, [false] upper bound
//     * \param hx длина сегмента
//     * \param repr part of selected points
//     */
//    template <Unsigned n> inline void extract_bound(shared<curve<n>> f, bool lower,
//        double hx=0.1, double repr=0.3)
//    {
//        range_assert(hx, 0.01, 1.);
//        range_assert(repr, 0., 1.);

//        double len = f->back()[0] - f->front()[0];
//        hx *= len;

//        auto first = f->begin();
//        auto last = first;

//        Unsigned k1 = 0, k2 = 0;
//        double xx = (*first)[0] + hx;
//        while (first != f->end())
//        {
//            while (last != f->end() && (*last)[0] <= xx) ++last;
//            k2 = std::distance(f->begin(), last);

//            if (lower) std::sort(&(*f)[k1], &(*f)[k2], [](auto a, auto b) { return a[1] < b[1]; });
//            else std::sort(&(*f)[k1], &(*f)[k2], [](auto a, auto b) { return a[1] > b[1]; });

//            Unsigned i = k1 + Unsigned(std::ceil(repr * (k2 - k1)));
//            for (; i<k2; ++i) (*f)[i][0] = inf<double>();

//            xx = xx + hx;
//            k1 = k2;
//            first = last;
//        }
//        std::sort(f->begin(), f->end(), [](auto a, auto b) { return a[0] < b[0]; });
//        f->erase(std::find_if(f->begin(), f->end(), [](auto a) { return a[0] == inf<double>(); }), f->end());
//    }

//    /*!
//     * \brief объединение близких точек (x[0] - cреднее, x[i != 0] - суммы)
//     * \return
//     */
//    template <Unsigned n> inline void unite(shared<curve<n>> f, double hx)
//    {
//        range_assert(hx, 0.01, 1.);
//        double len = f->back()[0] - f->front()[0];
//        hx *= len;

//        double invalid_value = f->back()[0] + len;
//        auto first = f->begin();
//        auto last = first;

//        Unsigned k1 = 0, k2 = 0;
//        double xx = (*first)[0] + hx;
//        while (first != f->end())
//        {
//            while (last != f->end() && (*last)[0] <= xx) ++last;
//            k2 = std::distance(f->begin(), last);

//            if (k1 != k2)
//            {
//                typename curve<n>::value_type sum = 0;
//                for (Unsigned i=k1; i<k2; ++i) { sum += (*f)[i]; (*f)[i] = invalid_value; }
//                (*f)[k1] = sum;  // set sums
//                (*f)[k1][0] *= 1./(k2 - k1); // set average
//            }
//            xx = xx + hx;
//            k1 = k2;
//            first = last;
//        }

//        auto it = std::stable_partition(f->begin(), f->end(), [invalid_value](auto a) { return a[0] < invalid_value; });
//        f->erase(it, f->end());
//    }

//    /*!
//     * \brief make sequence of partial sum
//     * \param f
//     * \param norm normalization to 1
//     * \param ndx
//     */
//    template <Unsigned n> inline void partial_sum(shared<curve<n>> f, bool norm, Unsigned ndx=1)
//    {
//        for (Unsigned i=1; i<f->size(); i++) (*f)[i][ndx] += (*f)[i-1][ndx];
//        if (norm)
//        {
//            double norm_coef = 1. / f->back()[ndx];
//            for (auto &e : *f) e[ndx] *= norm_coef;
//        }
//    }

//    /*!
//     * \brief fill gaps by previous values
//     * \param f curve
//     * \param hx step
//     */
//    template <Unsigned n> inline void fill_gaps(shared<curve<n>> f, double hx,
//        double left_indent, double right_indent)
//    {
//        double len = f->back()[0] - f->front()[0];
//        assert(len > 0);

//        hx *= len;
//        left_indent *= len;
//        right_indent *= len;
//        double xa = f->front()[0] - left_indent;
//        double xb = f->back()[0] + right_indent;

//        shared<curve<2>> tmp(new curve<2>);
//        Unsigned i = 0, nn = f->size(); double prev = 0;
//        for (double x=xa; x<=xb; x+=hx)
//        {
//            while (i < nn && (*f)[i][0] <= x) { prev = (*f)[i][1]; i++; }
//            tmp->push_back({x, prev});
//        }
//        if (i != nn) tmp->push_back({xb, f->back()[1]});
//        tmp->swap(*f);
//    }


//    template <Unsigned n> void to_height(shared<curve<n>> f, double h)
//    {
//        for (auto &e : *f) e[0] = relative_height(e[0], h);
//        std::sort(f->begin(), f->end(), [](auto a, auto b) { return a[0] < b[0]; });
//            // упорядочим, так как порядок изменился
//    }


//    template <Unsigned n> inline void to_start_day(shared<curve<n>> f, double day)
//        { for (auto &e : *f) e[1] = e[1] - day; }

//    template <Unsigned n> inline void to_relative(shared<curve<n>> f, shared<curve<n>> g)
//    {
//        for (auto &e : *f)
//        {
//            double h = e[0];
//            auto it = std::lower_bound(g->begin(), g->end(), h, [](auto a, auto h) { return a[0] < h; });
//            e[1] = e[1] - (*it)[1];
//        }
//    }

//    template <Unsigned n> inline void sort(shared<curve<n>> f, Unsigned i)
//        { std::sort(f->begin(), f->end(), [i](auto a, auto b){ return a[i] < b[i]; }); }

//    template <Unsigned n> inline void cut_below(shared<curve<n>> f)
//        { for (auto &e : *f) { if (e[1] < 0) e[1] = 0; } }

//    template <Unsigned n> inline void cut_left(shared<curve<n>> f)
//        { for (auto &e : *f) { if (e[0] < 0) e[1] = 0; } }

//    template <Unsigned n> inline void to_range(shared<curve<n>> f, Interval<double, 2> range, Unsigned i)
//    {
//        double a = range.minimal_value();
//        double b = range.maximal_value();

//        // сдвинем близкие к краям точки на края интервала, чтобы не потерять краевые точки
//        auto it = std::lower_bound(f->begin(), f->end(), a, [](auto a, auto x){ return a[0] < x; });
//        if (it != f->begin()) { --it; (*it)[0] = a; }

//        it = std::lower_bound(f->begin(), f->end(), b, [](auto a, auto x){ return a[0] < x; });
//        if (it != f->end()) { (*it)[0] = b; }

//        auto last = std::stable_partition(f->begin(), f->end(), [a, b, i](auto v){ return a <= v[i] && v[i] <= b; });
//        f->erase(last, f->end());
//    }

//    template <Unsigned n> inline void to_xrange(shared<curve<n>> f, Interval<double, 2> range)
//        { to_range(f, range, 0); }

//    template <Unsigned n> inline void to_yrange(shared<curve<n>> f, Interval<double, 2> range)
//        { to_range(f, range, 1); }

//    template <Unsigned n> inline void exclude_x(shared<curve<n>> f)
//    {
//        for (auto &e : *f) { e[0] = e[1]; e[1] = 1; }
//        return f;
//    }

//    template <Unsigned n> inline void to_normalized(shared<curve<n>> f, double norm_value)
//    {
//        double sum = 0;
//        for (auto e : *f) sum += e[1];
//        double norm_coef = norm_value / sum;
//        for (auto &e : *f) e[1] *= norm_coef;
//    }

//    template <Unsigned n> inline void round(shared<curve<n>> f)
//        { for (auto &e : *f) e[1] = std::round(e[1]); }

//    template <Unsigned n> inline void floor(shared<curve<n>> f)
//        { for (auto &e : *f) e[1] = std::floor(e[1]); }

//    template <Unsigned n> inline void ceil(shared<curve<n>> f)
//        { for (auto &e : *f) e[1] = std::ceil(e[1]); }


//    /// создать множество событий по разностям элементов пр одном значении x
//    template <Unsigned n> inline shared<curve<2>> to_beetween(shared<curve<n>> f, double eps)
//    {
//        double xprev = f->front()[0];
//        double yprev = f->front()[1];

//        shared<curve<n>> tmp(new curve<n>);
//        for (Unsigned i=1; i<f->size(); ++i)
//        {
//            double x = (*f)[i][0];
//            double y = (*f)[i][1];
//            if (equal(xprev, x, eps))
//            {
//                tmp->push_back({std::abs(y - yprev), 1});
//            }
//            xprev = x;
//            yprev = y;
//        }
//        return tmp;
//    }

//    template <Unsigned n>
//    inline shared<curve<2>> extract(shared<curve<n>> f, Unsigned ix, Unsigned iy)
//    {
//        shared<curve<2>> tmp(new curve<2>);
//        for (auto e : *f) tmp->push_back({e[ix], e[iy]});
//        return tmp;
//    }

//    template <Unsigned n>
//    inline shared<curve<2>> extract_mean(shared<curve<n>> f) { return extract(f, 0, 1); }

//    template <Unsigned n>
//    inline shared<curve<2>> extract_errors_range(shared<curve<n>> f)
//    {
//        shared<curve<2>> r(new curve<2>);

//        for (Unsigned i=0; i<f->size(); i++)
//        {
//            double x = (*f)[i][0];
//            double y = (*f)[i][1];
//            double s = (*f)[i][2];
//            r->push_back({x, y + s});
//        }
//        for (Unsigned i=f->size()-1; (Signed)i>=0; i--)
//        {
//            double x = (*f)[i][0];
//            double y = (*f)[i][1];
//            double s = (*f)[i][2];
//            r->push_back({x, y - s});
//        }
//        return r;
//    }

//    template <Unsigned n>
//    inline flagged_value<typename curve<n>::value_type> find_exact(shared<curve<n>> f, double x, double eps=0.001)
//    {
//        if (f->empty()) return {typename curve<n>::value_type(), false};
//        if (less(x, f->front()[0], eps)) return {f->front()[0], false};
//        if (greater(x, f->back()[0],eps)) return {f->back()[0], false};

//        auto it = std::lower_bound(f->begin(), f->end(), x,
//            [](auto a, double x) { return a[0] < x; }
//        );
//        assert (it != f->end());
//        if (equal((*it)[0], x, eps)) return { *it, true};
//        if (f->size() > 1 && equal((*(it-1))[0], x, eps)) return { *(it-1), true};

//        return {typename curve<n>::value_type(), false};
//    }

//    template <Unsigned n>
//    inline flagged_value<typename curve<n>::value_type> find_nearest(shared<curve<n>> f, double x)
//    {
//        if (f->empty()) return {typename curve<n>::value_type(), false};
//        if (x < f->front()[0]) return {f->front()[0], false};
//        if (x > f->back()[0] ) return {f->back()[0], false};

//        auto it = std::lower_bound(f->begin(), f->end(), x,
//            [](auto a, double x) { return a[0] < x; }
//        );
//        if (it == f->end()) --it;
//        return {*it, true};
//    }

//    /*!
//     * \brief make_cumulative_density_function
//     * \param f event function of the form {value1, 1}, {value2, 1}, ..
//     * \param hx
//     * \return
//     */
//    inline void to_cumulative_density_function(shared<curve<2>> f, double hx,
//        double left_indent, double right_indent)
//    {
//        sort(f, 0);
//        unite(f, hx);
//        partial_sum(f, true);
//        fill_gaps(f, hx, left_indent, right_indent);
//    }

//    inline shared<curve<2>> make_cumulative_density_function(shared<curve<1>> f, double hx,
//        double left_indent, double right_indent)
//    {
//        shared<curve<2>> tmp(new curve<2>);
//        for (auto e : *f) tmp->push_back({e[0], 1.});
//        to_cumulative_density_function(tmp, hx, left_indent, right_indent);
//        return tmp;
//    }
//    inline shared<curve<2>> make_cumulative_density_function(shared<curve<2>> f, double hx,
//        double left_indent, double right_indent)
//    {
//        shared<curve<2>> tmp(new curve<2>(*f));
//        to_cumulative_density_function(tmp, hx, left_indent, right_indent);
//        return tmp;
//    }

//    template <Unsigned n> inline void smooth(shared<curve<n>> f, int k)
//        { auto g = smooth(k, f); g->swap(*f); }

//    inline shared<curve<2>> make_histogram(shared<curve<2>> f, double hx=0.05)
//    {
//        auto tmp = make_cumulative_density_function(f, hx, 0.3, 0.3);
//        return smooth(1, tmp);
//    }

        /*!
         * \brief options to transform curve
         * \note all values are relative
         */
        struct transform_options
        {
            // совокупность преобразований для извлечения части данных
            flagged_value<double> min_bound = {0.01, false}; // извлечь нижнюю границу {eps, ..}
            flagged_value<double> max_bound = {0.01, false}; // извлечь верхнюю границу {eps, ..}

            // совокупность преобразований для получения распределений
            flagged_value<double> unite = {0.1, false}; // суммировать одинаковые значения
            flagged_value<double> accumulate = {0., false}; // создать массив частичных сумм и нормировать на значение
            flagged_value<double> left_indent = {0., false}; // вставить отступ слева
            flagged_value<double> right_indent = {0., false}; // вставить отступ справа
            flagged_value<double> fill_gaps = {0.01, false}; // заполнить пропуски предыдущими значениями (построение распределений)

            // совокупность преобразований для приведения координат представления графиков
            flagged_value<double> relative_height = {0, false}; // привести к относительной высоте
            flagged_value<Time> day_align_to = {0., false }; // привести к дате начала фазы
        };

        #define DEFAULT1_OPTIONS {0.01, false}, {0.01, false}
        #define DEFAULT2_OPTIONS {0.1, false}, {0., false}, {0., false}, {0., false}, {0.01, false}
        #define DEFAULT3_OPTIONS {0, false}, {0., false }

        inline transform_options default_distribution_transform_options()
        {
            transform_options option =
            {
                DEFAULT1_OPTIONS,
                {0.01, true}, {1., true}, {0.3, true}, {0.3, true}, {0.05, true},
                DEFAULT3_OPTIONS
            };
            return option;
        }

        inline transform_options default_axis_phase_transform_options(bool start, double height, Time day)
        {
            transform_options option =
            {
                {0.01, start}, {0.01, !start},
                DEFAULT2_OPTIONS,
                {height, true}, {day.to_double(), true}
            };
            return option;
        }

        inline transform_options default_axis_density_transform_options(double height)
        {
            transform_options option =
            {
                DEFAULT1_OPTIONS,
                {0.01, true}, {0., false}, {0., false}, {0., false}, {0.01, false},
                {height, true}, {0., false}
            };
            return option;
        }

        #undef DEFAULT1_OPTIONS
        #undef DEFAULT2_OPTIONS
        #undef DEFAULT3_OPTIONS


        inline shared<curve<2>> transform(shared<curve<2>> p, const transform_options &q)
        {
            shared<curve<2>> f(new curve<2>(*p));
            if (p->empty()) return f;

            double len = p->back()[0] - p->front()[0];
                // диапазон для определения абсолютных значений

            if (q.min_bound.flag())
            {
                double eps = q.min_bound.value() * len;
                std::sort(f->begin(), f->end(),
                    [eps](auto a, auto b)
                    {
                        if (equal(a[0], b[0], eps)) return a[1] < b[1];
                        return a[0] < b[0];
                    }
                );
                auto last = std::unique(f->begin(), f->end(),
                    [eps](auto a, auto b) { return equal(a[0], b[0], eps); }
                );
                f->erase(last, f->end());
            }

            if (q.max_bound.flag())
            {
                double eps = q.max_bound.value() * len;
                std::sort(f->begin(), f->end(),
                    [eps](auto a, auto b)
                    {
                        if (equal(a[0], b[0], eps)) return a[1] > b[1];
                        return a[0] < b[0];
                    }
                );
                auto last = std::unique(f->begin(), f->end(),
                    [eps](auto a, auto b) { return equal(a[0], b[0], eps); }
                );
                f->erase(last, f->end());
            }

            if (q.unite.flag())
            {
                double eps = q.unite.value() * len;

                Unsigned iprev = 0; double vprev = (*f)[0][0];
                for (Unsigned i=1; i<f->size(); i++)
                {
                    if (equal((*f)[i][0], vprev, eps))
                    {
                        (*f)[iprev][1] += (*f)[i][1];
                        (*f)[i][1] = 0;
                    }
                    else
                    {
                        vprev = (*f)[i][0];
                        iprev = i;
                    }
                }
                auto last = std::stable_partition(f->begin(), f->end(), [](auto a){ return a[1] > 0; });
                f->erase(last, f->end());
                std::cout << "1unit:" << *f << "\n"; std::cout.flush();
            }

            if (q.accumulate.flag())
            {
                for (Unsigned i=1; i<f->size(); i++) (*f)[i][1] += (*f)[i-1][1];

                if (q.accumulate.value() > 0) // нормировка под 1 на границе
                {
                    double norm_coef = q.accumulate.value() / f->back()[1];
                    for (auto &e : *f) e[1] *= norm_coef;
                }
                std::cout << "1acc:" << *f << "\n"; std::cout.flush();
            }

            if (q.fill_gaps.flag() && len > 0)
            {
                double hx = q.fill_gaps.value() * len;

                double da = q.left_indent.value();  if (q.left_indent.flag()) da *= len;
                double db = q.right_indent.value();  if (q.right_indent.flag()) db *= len;
                double xa = f->front()[0] - da;
                double xb = f->back()[0] + db;

                shared<curve<2>> tmp(new curve<2>);
                Unsigned i = 0, n = f->size(); double prev = 0;
                for (double x=xa; x<=xb; x+=hx)
                {
                    while (i < n && (*f)[i][0] <= x) { prev = (*f)[i][1]; i++; }
                    tmp->push_back({x, prev});
                }
                if (i != n) tmp->push_back({xb, f->back()[1]});
                f = tmp;
                std::cout << "1gap:" << *f << "\n"; std::cout.flush();
            }

            if (q.relative_height.flag()) // нормировка по высоте
            {
                double n = q.relative_height.value();
                for (auto &e : *f) e[0] = relative_height(e[0], n);
                std::sort(f->begin(), f->end(), [](auto a, auto b) { return a[0] < b[0]; });
                    // упорядочим по высоте для сглаживания
            }

            if (q.day_align_to.flag()) // смещение дней
            {
                double t = q.day_align_to.value().to_double();
                for (auto &e : *f) e[1] = e[1] - t;
            }

            return f;
        }


}
#endif // FLOWERLINE
