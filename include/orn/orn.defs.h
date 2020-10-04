#ifndef ORN_DEFS_H
#define ORN_DEFS_H

#include <unistd.h> // unix
#include <cstddef>
#include <numeric>
#include <array>
#include <type_traits>

#include "defs.h"
#include "math/random.h"
#include "math/WLSM.h"

namespace prg
{
    extern const char *species[];

    extern shared<class Flowerdome> rflowerdome;
    extern shared<class Flowerdome> wflowerdome;
    extern shared<class Sproutmodel> sproutmodel;

    extern std::string HOME_DIR; // домашнаяя дирктория пользователя
    extern std::string TMP_DIR; // директория временных файлов
    extern std::string WORK_DIR; // директория загрузки данных
    extern std::string OUTPUT_DIR; // выходная директория
    extern std::string LOCALE; // идентификатор локали

    extern std::string PDF_VIEW; // команда запуска просмотра pdf
    extern std::string SVG_TO_PDF; // команда конвертации svg файла в pdf

    extern Unsigned active_sprout;
    extern std::string active_model;
    extern Unsigned active_plot;
    const double default_width = 0.1; // относительная ширина сглаживания графиков

    extern std::map<std::string, shared<Flowerdome>> flowerdomes;
    extern std::map<std::string, shared<Sproutmodel>> sproutmodels;
    extern std::set<std::string> iflowerdomes;

    extern std::string read_flowerdome;
    extern std::string read_model;
    extern std::string write_flowerdome;
    extern std::string write_model;


    /// максимальное число побегов в одном файле
    const unsigned NUM_SPROUT = 3;
    const unsigned NUM_SPECIES = 9;
    const unsigned MAX_FLOWERING_TIME = 90;

    extern const char *species_name[];
    extern const std::string PROJECT_DIR;
    extern const std::string DATA_DIR;
    extern const std::string FLOWER_DIR;
    extern const std::string RAWDATA_DIR;
    extern const std::string MAP_FILE;

    /*!
     * \file файл стандартных определений используемых типов проекта и функции отладки
     */

    const size_t BUF_LEN = 2048; // длина записи в файле
    const int ITERATIONS = 150; // стандартное число итераций оптимизации
    const double THICKNESS = 0.1; // доля объектов, принадлежащих границе

    /// получить имя файла для заданного вида
    inline std::string make_ifilename(const std::string &dirname, unsigned species)
    { return dirname + species_name[species]; }

    /// получить имя файла для заданного вида
    inline std::string make_ifilename(const std::string &dirname, const std::string &species)
    { return dirname + species; }

    /// получить имя файла для заданного вида
    inline std::string make_ofilename(const std::string &dirname, unsigned species,
        const std::vector<unsigned> &sprouts)
    {
        std::string postfix("{" + std::to_string(sprouts[0]));
        for (unsigned i=1; i<sprouts.size(); i++)
            postfix += "," + std::to_string(sprouts[i]);
        postfix += "}";
        return dirname + std::string(species_name[species]) + postfix + ".pdf";
    }

    /// получить имя файла для заданного вида
    inline std::string make_ofilename(const std::string &dirname, const std::string &species,
        const std::vector<unsigned> &sprouts)
    {
        std::string postfix("{" + std::to_string(sprouts[0]));
        for (unsigned i=1; i<sprouts.size(); i++)
            postfix += "," + std::to_string(sprouts[i]);
        postfix += "}";
        return dirname + species + postfix + ".pdf";
    }

    /// литеральная строка со всеми пустыми символами и кавычками, удаляемая trim
    const char * const remove_symbols = " \"\t\n\r\f\v";

    /// выход из программы при отсутствии необходимого файла
    inline void no_file_exit [[noreturn]] (const std::string &filename)
    {
        std::printf("[ERROR] no file %s\n", filename.c_str());
        std::cout.flush();
        throw std::exception();
    }

    /// перевод доли [0,1] в проценты
    inline double to_percent(double v) { return v * 100; }

    /// выход из программы при отсутствии необходимого файла
    #define OPEN_FILE(file, filename) \
        std::ifstream file(filename); \
        if (!file) no_file_exit(filename); \

    template <typename T>
    void DEBUG_PRINT(T t) { std::cout << t << "\n"; std::cout.flush(); }

    template <typename T, typename U>
    void DEBUG_PRINT(T t, U u) { std::cout << t << " " << u << "\n"; std::cout.flush(); }

    template <typename T>
    void DEBUG_PRINT(const std::vector<T> &t)
    { for (auto e:t) std::cout << e << "\n"; std::cout.flush(); }

    template <typename T> void _DEBUG_PRINT(const T &) {}

    template <Unsigned n>
    inline shared<curve<n>> operator+(shared<curve<n>> a, shared<curve<n>> b)
    {
        shared<curve<n>> out(new curve<n>(*a));
        for (auto e : *b) out->push_back(e);
        std::sort(out->begin(), out->end(), [](auto a, auto b){ return a[0] < b[0]; });
        return out;
    }

    template <unsigned n>
    inline curve<n> lower_bound(const curve<n> &f, double portion=0.5, unsigned ix=0, unsigned iy=1)
    {
        if (f.empty()) return curve<n>();
        curve<n> tmp;

        std::map<double, std::vector<double>> m;
        for (auto e : f)
        {
            double h = e[ix];
            double v = e[iy];
            if (m.count(h) == 0) m[h] = std::vector<double>();
            m[h].push_back(v);
        }

        for (auto &e : m)
        {
            auto h = e.first;
            auto &vec = e.second;
            std::sort(vec.begin(), vec.end());

            double ts = vec.front();
            double te = vec.back();
            double tcut = ts + portion * (te - ts);
            auto it = std::upper_bound(vec.begin(), vec.end(), tcut,
                [](double tcut, double a){ return tcut < a; }
            );
            vec.erase(it, vec.end());
            for (auto e : vec) tmp.push_back({h, e});
        }
        assert(prg::is_sorted_by_column<n>(tmp, 0));

        return tmp;
    }

    //---------------------------------------------------------------------------------------------
    //                              micro-functions
    //---------------------------------------------------------------------------------------------
    /*!
     * \brief relative height of the flower on sprout
     * \param f flower 1st-order index
     * \param n number of flowers on the 1st axis
     * \return
     */
    inline double relative_height(double f, double fmax) { return (fmax - f) / (fmax - 1); }

    /// получить номер цветка 1-го порядка в зависимости от заданной высоты
    inline Unsigned number_from_relative_height(double h, double fmax)
        { return Unsigned(std::round(fmax - h * (fmax - 1))); }


    const WLSM_item default_WLSM_options =
    {
        2, // мощность полинома
        40, // число итераций подгонки полинома
        {0.01, true}, // шаг сглаживания
        {0.05, true}, // окно сглаживания
        {inf<double>(),-inf<double>(0)}, // интервал сглаживания
        {0.1, true}, // отступ слева при отсутствии задания интервала
        {0.1, true}, // отступ справа при отсутствии задания интервала
        false // флаг печати процесса подгонки
    };

    /// вызов объекта, который выполняет сглаживание
    inline shared<curve<2>> smooth(Unsigned k, shared<curve<2>> p, const WLSM_item &item=WLSM_item())
        { return WLSM().smooth(k, p, item); }

    inline shared<curve<3>> smooth_with_errors(Unsigned k, shared<curve<2>> p, const WLSM_item &item=WLSM_item())
        { return WLSM().smooth_with_errors(k, p, item); }

    /*!
     * \brief extract lower|upper bound of ordered curve
     * \note область определения делится на сегменты длиной hx и в каждом сегменте выбирается
     *  верхний/нижний представитель
     * \param f функция
     * \param lower [true] lower bound, [false] upper bound
     * \param hx длина сегмента
     * \param repr part of selected points
     */
    template <Unsigned n> inline void extract_bound(shared<curve<n>> f, bool lower,
        double hx=0.1, double repr=0.3)
    {
        range_assert(hx, 0.01, 1.);
        range_assert(repr, 0., 1.);

        double len = f->back()[0] - f->front()[0];
        hx *= len;

        auto first = f->begin();
        auto last = first;

        Unsigned k1 = 0, k2 = 0;
        double xx = (*first)[0] + hx;
        while (first != f->end())
        {
            while (last != f->end() && (*last)[0] <= xx) ++last;
            k2 = std::distance(f->begin(), last);

            if (lower) std::sort(&(*f)[k1], &(*f)[k2], [](auto a, auto b) { return a[1] < b[1]; });
            else std::sort(&(*f)[k1], &(*f)[k2], [](auto a, auto b) { return a[1] > b[1]; });

            Unsigned i = k1 + Unsigned(std::ceil(repr * (k2 - k1)));
            for (; i<k2; ++i) (*f)[i][0] = inf<double>();

            xx = xx + hx;
            k1 = k2;
            first = last;
        }
        std::sort(f->begin(), f->end(), [](auto a, auto b) { return a[0] < b[0]; });
        f->erase(std::find_if(f->begin(), f->end(), [](auto a) { return a[0] == inf<double>(); }), f->end());
    }

    /*!
     * \brief объединение близких точек (x[0] - cреднее, x[i != 0] - суммы)
     * \return
     */
    template <Unsigned n> inline void unite(shared<curve<n>> f, double hx)
    {
        range_assert(hx, 0.01, 1.);
        double len = f->back()[0] - f->front()[0];
        hx *= len;

        double invalid_value = f->back()[0] + ((len > 0) ? len : 1);
        auto first = f->begin();
        auto last = first;

        Unsigned k1 = 0, k2 = 0;
        double xx = (*first)[0] + hx;
        while (first != f->end())
        {
            while (last != f->end() && (*last)[0] <= xx) ++last;
            k2 = std::distance(f->begin(), last);

            if (k1 != k2)
            {
                typename curve<n>::value_type sum = 0;
                for (Unsigned i=k1; i<k2; ++i) { sum += (*f)[i]; (*f)[i] = invalid_value; }
                (*f)[k1] = sum;  // set sums
                (*f)[k1][0] *= 1./(k2 - k1); // set average
            }
            xx = xx + hx;
            k1 = k2;
            first = last;
        }

        auto it = std::stable_partition(f->begin(), f->end(), [invalid_value](auto a) { return a[0] < invalid_value; });
        f->erase(it, f->end());
    }

    /*!
     * \brief make sequence of partial sum
     * \param f
     * \param norm normalization to 1
     * \param ndx
     */
    template <Unsigned n> inline void partial_sum(shared<curve<n>> f, bool norm, Unsigned ndx=1)
    {
        for (Unsigned i=1; i<f->size(); i++) (*f)[i][ndx] += (*f)[i-1][ndx];
        if (norm)
        {
            double norm_coef = 1. / f->back()[ndx];
            for (auto &e : *f) e[ndx] *= norm_coef;
        }
    }

    /*!
     * \brief fill gaps by previous values
     * \param f curve
     * \param hx step
     */
    template <Unsigned n> inline void fill_gaps(shared<curve<n>> f, double hx,
        double left_indent, double right_indent)
    {
        double len = f->back()[0] - f->front()[0];
        len  = len > 0 ? len : 1;

        hx *= len;
        left_indent *= len;
        right_indent *= len;
        double xa = f->front()[0] - left_indent;
        double xb = f->back()[0] + right_indent;

        shared<curve<2>> tmp(new curve<2>);
        Unsigned i = 0, nn = f->size(); double prev = 0;
        for (double x=xa; x<=xb; x+=hx)
        {
            while (i < nn && (*f)[i][0] <= x) { prev = (*f)[i][1]; i++; }
            tmp->push_back({x, prev});
        }
        if (i != nn) tmp->push_back({xb, f->back()[1]});
        tmp->swap(*f);
    }


    template <Unsigned n> void to_height(shared<curve<n>> f, double h)
    {
        for (auto &e : *f) e[0] = relative_height(e[0], h);
        std::sort(f->begin(), f->end(), [](auto a, auto b) { return a[0] < b[0]; });
            // упорядочим, так как порядок изменился
    }


    template <Unsigned n> inline void to_start_day(shared<curve<n>> f, double day)
        { for (auto &e : *f) e[1] = e[1] - day; }

    template <Unsigned n> inline void to_relative(shared<curve<n>> f, shared<curve<n>> g)
    {
        for (auto &e : *f)
        {
            double h = e[0];
            auto it = std::lower_bound(g->begin(), g->end(), h, [](auto a, auto h) { return a[0] < h; });
            e[1] = e[1] - (*it)[1];
        }
    }

    template <Unsigned n> inline void sort(shared<curve<n>> f, Unsigned i)
    {
        if (f->empty()) return;
        std::sort(f->begin(), f->end(), [i](auto a, auto b){ return a[i] < b[i]; });
    }

    template <Unsigned n> inline void cut_below(shared<curve<n>> f)
    {
        if (f->empty()) return;
        for (auto &e : *f) { if (e[1] < 0) e[1] = 0; }
    }

    template <Unsigned n> inline void cut_below(shared<curve<n>> f, double x1, double x2)
    {
        if (f->empty()) return;
        for (auto &e : *f) { if (e[1] < 0 && e[0] >= x1 && e[0] <= x2) e[1] = 0; }
    }

    template <Unsigned n> inline void cut_up(shared<curve<n>> f, double up=1)
    {
        if (f->empty()) return;
        for (auto &e : *f) { if (e[1] > up) e[1] = up; }
    }

    template <Unsigned n> inline void cut_left(shared<curve<n>> f)
    {
        if (f->empty()) return;
        for (auto &e : *f) { if (e[0] < 0) e[1] = 0; }
    }

    template <Unsigned n> inline shared<curve<n>> cut_zerobounds(shared<curve<n>> f)
    {
        shared<curve<n>> tmp(new curve<n>);
        if (f->empty()) return tmp;

        Unsigned nn = f->size();
        Unsigned n1 = 0; while (n1 < nn && (*f)[n1][1] <= 0) n1++;
        Unsigned n2 = nn; while (n2 > 0 && (*f)[n2-1][1] <= 0) n2--;
//        if (n1 > 0) n1--;
//        if (n2 < nn) n2++;

        for (Unsigned i=n1; i<n2; i++) tmp->push_back((*f)[i]);
        return tmp;
    }

    template <Unsigned n> inline double first_nonzero(shared<curve<n>> f)
        { double x = -inf<double>(); for (auto e : *f) if (e[1] > 0) { x = e[0]; break; }  return x; }

    template <Unsigned n> inline void to_range(shared<curve<n>> f, Interval<double, 2> range, Unsigned i)
    {
        double a = range.minimal_value();
        double b = range.maximal_value();

        // сдвинем близкие к краям точки на края интервала, чтобы не потерять краевые точки
        auto it = std::lower_bound(f->begin(), f->end(), a, [](auto a, auto x){ return a[0] < x; });
        if (it != f->begin()) { --it; (*it)[0] = a; }

        it = std::lower_bound(f->begin(), f->end(), b, [](auto a, auto x){ return a[0] < x; });
        if (it != f->end()) { (*it)[0] = b; }

        auto last = std::stable_partition(f->begin(), f->end(), [a, b, i](auto v){ return a <= v[i] && v[i] <= b; });
        f->erase(last, f->end());
    }

    template <Unsigned n> inline void to_xrange(shared<curve<n>> f, Interval<double, 2> range)
        { to_range(f, range, 0); }

    template <Unsigned n> inline void to_yrange(shared<curve<n>> f, Interval<double, 2> range)
        { to_range(f, range, 1); }

    template <Unsigned n> inline void exclude_x(shared<curve<n>> f)
    {
        for (auto &e : *f) { e[0] = e[1]; e[1] = 1; }
        return f;
    }

    template <Unsigned n> inline void to_normalized(shared<curve<n>> f, double norm_value)
    {
        double sum = 0;
        for (auto e : *f) sum += e[1];
        double norm_coef = norm_value / sum;
        for (auto &e : *f) e[1] *= norm_coef;
    }

    template <Unsigned n> inline void round(shared<curve<n>> f)
        { for (auto &e : *f) e[1] = std::round(e[1]); }

    template <Unsigned n> inline void floor(shared<curve<n>> f)
        { for (auto &e : *f) e[1] = std::floor(e[1]); }

    template <Unsigned n> inline void ceil(shared<curve<n>> f)
        { for (auto &e : *f) e[1] = std::ceil(e[1]); }


    /// создать множество событий по разностям элементов пр одном значении x
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

    inline shared<curve<2>> to_beetween(shared<curve<2>> f, double eps)
    {
        double xprev = f->front()[0];
        double yprev = f->front()[1];

        shared<curve<2>> tmp(new curve<2>);
        for (Unsigned i=1; i<f->size(); ++i)
        {
            double x = (*f)[i][0];
            double y = (*f)[i][1];
            if (equal(xprev, x, eps))
            {
                tmp->push_back({std::abs(y - yprev), 1});
            }
            xprev = x;
            yprev = y;
        }
        return tmp;
    }

    template <Unsigned n>
    inline shared<curve<2>> extract(shared<curve<n>> f, Unsigned ix, Unsigned iy)
    {
        shared<curve<2>> tmp(new curve<2>);
        for (auto e : *f) tmp->push_back({e[ix], e[iy]});
        return tmp;
    }

    template <Unsigned n>
    inline shared<curve<2>> extract_mean(shared<curve<n>> f) { return extract(f, 0, 1); }

    template <Unsigned n>
    inline shared<curve<2>> extract_errors_range(shared<curve<n>> f)
    {
        shared<curve<2>> r(new curve<2>);

        for (Unsigned i=0; i<f->size(); i++)
        {
            double x = (*f)[i][0];
            double y = (*f)[i][1];
            double s = (*f)[i][2];
            r->push_back({x, y + s});
        }
        for (Unsigned i=f->size()-1; (Signed)i>=0; i--)
        {
            double x = (*f)[i][0];
            double y = (*f)[i][1];
            double s = (*f)[i][2];
            r->push_back({x, y - s});
        }
        return r;
    }

    template <Unsigned n>
    inline flagged_value<typename curve<n>::value_type> find_exact(shared<curve<n>> f, double x, double eps=0.001)
    {
        if (f->empty()) return {typename curve<n>::value_type(), false};
        if (less(x, f->front()[0], eps)) return {f->front()[0], false};
        if (greater(x, f->back()[0],eps)) return {f->back()[0], false};

        auto it = std::lower_bound(f->begin(), f->end(), x,
            [](auto a, double x) { return a[0] < x; }
        );
        assert (it != f->end());
        if (equal((*it)[0], x, eps)) return { *it, true};
        if (f->size() > 1 && equal((*(it-1))[0], x, eps)) return { *(it-1), true};

        return {typename curve<n>::value_type(), false};
    }

    template <Unsigned n>
    inline flagged_value<typename curve<n>::value_type> find_nearest(shared<curve<n>> f, double x)
    {
        if (f->empty()) return {typename curve<n>::value_type(), false};
        if (x < f->front()[0]) return {f->front()[0], false};
        if (x > f->back()[0] ) return {f->back()[0], false};

        auto it = std::lower_bound(f->begin(), f->end(), x,
            [](auto a, double x) { return a[0] < x; }
        );
        if (it == f->end()) --it;
        return {*it, true};
    }

    /*!
     * \brief получить минимальное значение y по заданной оси i
     * \param i номер оси y={1, 2}
     */
    template <Unsigned n>
    inline flagged_value<double> find_minimal(shared<curve<n>> f, Unsigned i)
    {
        if (f->empty()) return {0., false};
        double v = f->front()[i];

        for (auto e : *f) if (v > e[i]) v = e[i];
        return {v, true};
    }

    /*!
     * \brief получить максимальное значение y по заданной оси i
     * \param i номер оси y={1, 2}
     */
    template <Unsigned n>
    inline flagged_value<double> find_maximal(shared<curve<n>> f, Unsigned i)
    {
        if (f->empty()) return {0., false};
        double v = f->front()[i];

        for (auto e : *f) if (v < e[i]) v = e[i];
        return {v, true};
    }

    /*!
     * \brief make_cumulative_density_function
     * \param f event function of the form {value1, 1}, {value2, 1}, ..
     * \param hx
     * \return
     */
    inline void to_cumulative_density_function(shared<curve<2>> f, double hx,
        double left_indent, double right_indent)
    {
        sort(f, 0);
        unite(f, hx);
        partial_sum(f, true);
        fill_gaps(f, hx, left_indent, right_indent);
    }

    inline shared<curve<2>> make_cumulative_density_function(shared<curve<1>> f, double hx,
        double left_indent, double right_indent)
    {
        shared<curve<2>> tmp(new curve<2>);
        for (auto e : *f) tmp->push_back({e[0], 1.});
        to_cumulative_density_function(tmp, hx, left_indent, right_indent);
        return tmp;
    }
    inline shared<curve<2>> make_cumulative_density_function(shared<curve<2>> f, double hx,
        double left_indent, double right_indent)
    {
        shared<curve<2>> tmp(new curve<2>(*f));
        to_cumulative_density_function(tmp, hx, left_indent, right_indent);
        return tmp;
    }

    template <Unsigned n> inline void smooth(shared<curve<n>> f, int k)
        { auto g = smooth(k, f); g->swap(*f); }

    inline shared<curve<2>> make_histogram(shared<curve<2>> f, double hx=0.05)
    {
        auto tmp = make_cumulative_density_function(f, hx, 0.3, 0.3);
        return smooth(1, tmp);
    }

//    template <Unsigned N1, Unsigned N2>
//    inline shared<curve<N1>> minimize_dispersion(shared<curve<N1>> p1, shared<curve<N2>> f,
//        Unsigned steps=100)
//    {
//        Unsigned n = p1->size();
//        shared<curve<N1>> out(new curve<N1>(*p1));

//        for (Unsigned i=0; i<steps; i++)
//        {
//            std::uniform_int_distribution<> dist(0, n-2);
//            Unsigned i1 = dist(randgen());
//            Unsigned i2 = i1 + 1;

//            auto &e1 = (*p1)[i1];
//            auto &e2 = (*p1)[i2];
//            double x1 = e1[0];
//            double y1 = e1[1];
//            double x2 = e2[0];
//            double y2 = e2[1];

//            auto q1 = find_nearest(f, x1).value();
//            auto q2 = find_nearest(f, x2).value();
//            double disp_prev = sqr(q1[1] - e1[1]) + sqr(q2[1] - e2[1]);
//            double disp_new  = sqr(q1[1] - y2) + sqr(q2[1] - y1);
//            if (disp_new > disp_prev) continue;

//            e1[1] = y2;
//            e2[1] = y1;
//        }
//        return out;
//    }

    template <Unsigned N1, Unsigned N2>
    inline shared<curve<N1>> minimize_dispersion(shared<curve<N1>> p1, shared<curve<N2>> f,
        Unsigned steps=10)
    {
        Unsigned n = p1->size();
        shared<curve<N1>> out(new curve<N1>(*p1));

        for (Unsigned step=0; step<steps; step++)
        {
            for (Unsigned i1=0; i1<n-2; i1++)
            {
//            std::uniform_int_distribution<> dist(0, n-2);
//            Unsigned i1 = dist(randgen());
                Unsigned i2 = i1 + 1;

                auto &e1 = (*p1)[i1];
                auto &e2 = (*p1)[i2];
                double x1 = e1[0];
                double y1 = e1[1];
                double x2 = e2[0];
                double y2 = e2[1];

                auto q1 = find_nearest(f, x1).value();
                auto q2 = find_nearest(f, x2).value();
                double disp_prev = sqr(q1[1] - e1[1]) + sqr(q2[1] - e2[1]);
                double disp_new  = sqr(q1[1] - y2) + sqr(q2[1] - y1);
                if (disp_new > disp_prev) continue;

                e1[1] = y2;
                e2[1] = y1;
            }
        }
        return out;
    }

}
#endif // ORN_DEFS_H
