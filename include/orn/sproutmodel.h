#ifndef SPROUT_MODEL_H
#define SPROUT_MODEL_H

#include <string>
#include <vector>
#include <string>
#include <set>
#include <tuple>
#include <memory>

#include "container/vector.h"
#include "math/math.h"
#include "math/random.h"
#include "math/distribution.h"
#include "math/statistics.h"
#include "math/curve.h"
#include "os/os.h"

#include "orn/orn.defs.h"
#include "orn/time.h"
#include "orn/phase.h"
#include "orn/flower.h"
#include "orn/flowerline.h"
#include "orn/sprout.h"
#include "orn/flowerdome.h"
#include "orn/gnuplot.h"

/// определяется при тестировании правильности генерации распределений и регрессий
// #define TEST_PLOT
const Unsigned DISTRIBUTION_TEST_VALUES = 10;

namespace prg
{

    struct dependence_item
    {
        typedef void (Sproutmodel::*PF1)(dependence_item &, const Flowerdome &, const std::set<Unsigned> &);
        typedef void (Sproutmodel::*PF2)(dependence_item &, double);
        typedef std::string (Sproutmodel::*PF3)(dependence_item &, Unsigned);

        std::string title =""; // name of dependence
        PF1 pf1 = nullptr; // data preparation function
        PF2 pf2 = nullptr; // data smoothing function
        PF3 pf3 = nullptr; // test statistics function

        shared<numeric_distribution> dist = nullptr;
        std::string info = "";
        bool correlation_coef = false;
        bool determination_coef = false;
        bool normality_test = false;
        Unsigned normality_test_n = 100;

        shared<curve<2>> data1 = nullptr; // experimental data for 1st-order
        shared<curve<2>> data2 = nullptr; // experimental data for 2nd-order
        shared<curve<2>> data = nullptr; // sum of experimental data 1st + 2nd
        shared<curve<2>> ddata = nullptr; // diff of experimental data 1st - 2nd

        shared<curve<3>> smooth1 = nullptr;  // smoothed data with errors
        shared<curve<3>> smooth2 = nullptr;  // smoothed data with errors
        shared<curve<3>> smooth = nullptr;  // smoothed data with errors
        shared<curve<3>> dsmooth = nullptr;  // smoothed data with errors

        gnuplotParams gnuparam; // data plot parameters
    };

    /// построить имя переменной для распределения объектов(of) по параметру(by)
    #define DISTRIBUTION(of, by) distribution_var_##of##_##by

    /// построить имя перемнной для функции для регрессии параметра(of) от параметра(on)
    #define REGRESSION(of, on) regression_var_##of##_##on

    /// построить имя перемнной для функции для регрессии параметра(of) от параметра(on)
    #define CURVE(of, on) curve_var_##of##_##on

    /*!
     * \brief получить ширину сглаживания, обеспечивающую отсутствие артефактов
     * \note исходная ширина окна увеличивается до тех пор, пока гарантировано не охватит
     *  такое число точек, которое обеспечит гладкость
     * \param power мощность полинома
     * \param data набор данных
     * \param relative_width исходная ширина окна
     * \return относительная ширина окна
     */
    inline double smooth_width(Unsigned power, shared<curve<2>> data, double relative_width)
    {
        if (data->size() <= power + 1) return 1;

        double len = std::abs(data->back()[0] - data->front()[0]);
        double w = relative_width;
        for (Unsigned i=power; i<data->size(); i++)
        {
            double x1 = (*data)[i-power][0];
            double x2 = (*data)[i][0];
            double ww = std::abs(x2 - x1) / len;
            if (w < ww) w = ww;
        }
        return 0.5 * w;
//        return relative_width;
    }


    /*!
     * \brief сгладить распределение с заданной шириной окна
     * \param data случайные данные
     * \param width ширина окна
     * \return сглаженные данные
     */
    inline shared<curve<2>> smooth_distribution(shared<curve<2>> data, double width=default_width)
    {
        WLSM_item options(default_WLSM_options);
        width = smooth_width(options.power, data, width);
        options.width = {width, true};

        auto f = smooth(1, data, options);
        cut_below(f);
        cut_left(f);
        return f;
    }
    inline shared<curve<3>> distribution_smooth_with_errors(shared<curve<2>> data, double width=default_width)
    {
        WLSM_item options(default_WLSM_options);
        width = smooth_width(options.power, data, width);
        options.width = {width, true};

        auto f = WLSM().smooth_with_errors(1, data, options);
        cut_below(f);
        cut_left(f);
        return f;
    }

    /*!
     * \brief сгладить регрессию
     * \note результат сглаживания сводится к диапазону опредленеия данных
     * \param data случайные данные
     * \param width ширина сглаживания
     * \return сглаженная кривая со средней ошибкой сглаживания
     */
    inline shared<curve<3>> smooth_regression(shared<curve<2>> data, double width=default_width, bool correct_width=true)
    {
        if (data->empty()) return shared<curve<3>>(new curve<3>);

        WLSM_item options(default_WLSM_options);
        if (correct_width) width = smooth_width(options.power, data, width);
        options.width = {width, true};
        auto f = WLSM().smooth_with_errors(0, data, options);

        Interval<double,2> xrange(data->front()[0], data->back()[0]);
        to_xrange(data, xrange);
        to_xrange(f, xrange);
        return f;
    }

    inline std::string statistics_info(const std::vector<double> &tmp)
    {
        if (tmp.empty()) return "";

        Unsigned n = tmp.size();
        double mean = average(tmp);
        double med = median(tmp);
        double var = sample_variance(tmp);
        double stddev = std::sqrt(var);
        double skw = skewness(tmp);
        double krt = kurtosis(tmp);

        std::string msg = "\n";
        msg += "\tn = " + std::to_string(n) + "\n";
        msg += "\tmean = " + std::to_string(mean) + "\n";
        msg += "\tmedian = " + std::to_string(med) + "\n";
        msg += "\tvariance = " + std::to_string(var) + "\n";
        msg += "\tstddev = " + std::to_string(stddev) + "\n";
        msg += "\tskewness = " + std::to_string(skw) + "\n";
        msg += "\tkurtosis = " + std::to_string(krt) + "\n";
        msg += "\n";

        return msg;
    }

    inline std::string normality_test(shared<numeric_distribution> dist,
        Unsigned n=DISTRIBUTION_TEST_VALUES)
    {
    #ifdef NORMAL_STATISTICS_TEST
        std::vector<double> tmp;
        std::normal_distribution<> d{-20, 10};
        for (Unsigned i=0; i<n; i++) tmp.push_back(d(randgen()));
    #else
        std::vector<double> tmp = dist->operator()(randgen(), n);
    #endif

        std::string msg = "\n*** distribution normality test ***\n";
        msg += statistics_info(tmp);

        double mean = average(tmp);
        double var = sample_variance(tmp);
        double stddev = std::sqrt(var);

        auto [aad, pvalue] = anderson_darling_normality_test<0>(tmp, mean, stddev);

        msg += "\tA^2* (anderson darling test) = " + std::to_string(aad) + "\n";
        msg += "\tpvalue = ";
        if (pvalue > 0.00001) msg += std::to_string(pvalue);
        else
        {
            char buf[128];
            std::sprintf(buf, "%e", pvalue);
            msg += buf;
        }
        msg += "\n";
        msg += "\tnormality(0.05): \t";
        msg += pvalue > 0.05 ? "yes" : "no";
        return msg;
    }

    inline std::string regression_test(dependence_item &item, int64_t lags = -1, int64_t fit_dof = 1)
    {
        shared<curve<2>> data = item.data;
        shared<curve<3>> smooth = item.smooth;

        std::string msg = "\n";
        msg += "******************************************\n";
        msg += "***             statistics             ***\n";
        msg += "******************************************\n";
        if (item.correlation_coef)
        {
            std::vector<double> X,Y;
            for (auto e : *data) { X.push_back(e[0]); Y.push_back(e[1]); }

            double r = correlation_coefficient(X, Y);
            Unsigned n = data->size();
            auto [t, pvalue] = correlation_test(r, n);

            msg += "\tn = " + std::to_string(n) + "\n\n";
            msg += "\tr(X,Y) = " + std::to_string(r) + "\n";
            msg += "\tt (two-tailed t test) = " + std::to_string(t) + "\n";
            msg += "\tpvalue = " + pvalue_as_string(pvalue) + "\n";
            msg += "\tcorrelation exist(0.05): ";
            msg += pvalue > 0.05 ? "no" : "yes";
            msg += "\n\n";
        }

        if (item.determination_coef)
        {
            double R2 = determination_coef(smooth, data);
            msg += "\tR2 = " + std::to_string(R2) + "\n";
        }

        if (item.normality_test)
        {
            msg += "\n";
            msg += "******************************************\n";
            msg += "*** regression residues normality test ***\n";
            msg += "******************************************\n";


//#define LJUNG_TEST
#ifdef LJUNG_TEST
            std::vector<double> tmps;
            for (auto e : *smooth) tmps.push_back(e[1]);
            double aver = average(tmps);
            for (auto &e : *smooth) e[1] = aver;
#endif
            std::vector<double> residues;
            for (auto e : *data)
            {
                double x = e[0];
                double y = e[1];

                auto it = std::lower_bound(smooth->begin(), smooth->end(), x,
                    [](auto v, auto x){ return v[0] < x; });
                assert(it != smooth->end());

                double y1 = (*it)[1];
                residues.push_back(y1 - y);
            }
            auto [Q, pvalue] = ljung_box(residues, lags, fit_dof);
            // auto [Q, pvalue] = ljung_box(residues, 2, fit_dof);

            double mean = average(residues);
            double stddev = standard_derivation(residues);
            auto [A2, pvalue2] = anderson_darling_normality_test<0>(residues, mean, stddev);

            msg += statistics_info(residues);

            msg += "\tQ (ljung_box test) = " + std::to_string(Q) + "\n";
            msg += "\tpvalue = " + pvalue_as_string(pvalue) + "\n";
            msg += "\tno correlations(0.05): \t";
            msg += pvalue > 0.05 ? "yes" : "no";
            msg += "\n\n";

            msg += "\tA2 (anderson_darling test) = " + std::to_string(A2) + "\n";
            msg += "\tpvalue = " + pvalue_as_string(pvalue2) + "\n";
            msg += "\tnormality of rests(0.05): \t";
            msg += pvalue2 > 0.05 ? "yes" : "no";
            msg += "\n";
        }
        return msg;
    }

    inline std::string regression_test(shared<curve<2>> data, shared<curve<3>> smooth,
        dependence_item &item, int64_t lags = -1, int64_t fit_dof = 1)
    {
        std::string msg = "\n";
        if (item.correlation_coef)
        {
            std::vector<double> X,Y;
            for (auto e : *data) { X.push_back(e[0]); Y.push_back(e[1]); }

            double r = correlation_coefficient(X, Y);
            Unsigned n = data->size();
            auto [t, pvalue] = correlation_test(r, n);

            msg += "\tn = " + std::to_string(n) + "\n\n";
            msg += "\tr(X,Y) = " + std::to_string(r) + "\n";
            msg += "\tt (two-tailed t test) = " + std::to_string(t) + "\n";
            msg += "\tpvalue = " + pvalue_as_string(pvalue) + "\n";
            msg += "\tcorrelation exist(0.05): ";
            msg += pvalue > 0.05 ? "no" : "yes";
            msg += "\n\n";
        }

        if (item.determination_coef)
        {
            double R2 = determination_coef(smooth, data);
            msg += "\tR2 = " + std::to_string(R2) + "\n";
        }

        if (item.normality_test)
        {
            msg += "\n";
            msg += "--- regression residues normality test ---\n";

            std::vector<double> residues;
            for (auto e : *data)
            {
                double x = e[0];
                double y = e[1];

                auto it = std::lower_bound(smooth->begin(), smooth->end(), x,
                    [](auto v, auto x){ return v[0] < x; });
                assert(it != smooth->end());

                double y1 = (*it)[1];
                residues.push_back(y1 - y);
            }
            auto [Q, pvalue] = ljung_box(residues, lags, fit_dof);

            double mean = average(residues);
            double stddev = standard_derivation(residues);
            auto [A2, pvalue2] = anderson_darling_normality_test<0>(residues, mean, stddev);

            msg += statistics_info(residues);

            msg += "\tQ (ljung_box test) = " + std::to_string(Q) + "\n";
            msg += "\tpvalue = " + pvalue_as_string(pvalue) + "\n";
            msg += "\tno correlations(0.05): \t";
            msg += pvalue > 0.05 ? "yes" : "no";
            msg += "\n\n";

            msg += "\tA2 (anderson_darling test) = " + std::to_string(A2) + "\n";
            msg += "\tpvalue = " + pvalue_as_string(pvalue2) + "\n";
            msg += "\tnormality of rests(0.05): \t";
            msg += pvalue2 > 0.05 ? "yes" : "no";
            msg += "\n";
        }
        return msg;
    }
    /*!
     * \brief модель побега
     * \note модель состоит из совокупности распределений и регрессионных соотношений,
     * позволяющих по ним генерить образцы модельных побегов
     */
    class Sproutmodel
    {
    public:
        std::string info = ""; // some info
        Unsigned nflower1 = 0;
        Unsigned nflower2 = 0;

    protected:

        //------------------------------------------------------------------------------------------
        //              описание структуры побега (вплоть до осей 2-го порядка)
        //------------------------------------------------------------------------------------------
        shared<numeric_distribution> DISTRIBUTION(shoots, flowers1) // распределение побегов по числу цветков 1-го порядка
            = shared<numeric_distribution>(new numeric_distribution);
        shared<numeric_distribution> DISTRIBUTION(flowers, budding_time) // распределение цветков по времени бутонизации
            = shared<numeric_distribution>(new numeric_distribution);
        shared<numeric_distribution> DISTRIBUTION(flowers, flowering_time) // распределение цветков по времени цветения
            = shared<numeric_distribution>(new numeric_distribution);
        shared<numeric_distribution> DISTRIBUTION(flowers, ovaring_time) // распределение цветков по времени завязывания
            = shared<numeric_distribution>(new numeric_distribution);
        shared<numeric_distribution> DISTRIBUTION(flowers, between_time) // распределение пар цветков по времени между последовательными распусканиями на одной оси
            = shared<numeric_distribution>(new numeric_distribution);

        shared<curve<3>> REGRESSION(flowers2, flowers1); // зависимость числа цветков 2-го порядков от числа цветков 1-го порядка
        shared<curve<3>> REGRESSION(density2, height); // плотность цветков 2-го порядка по высоте

        //------------------------------------------------------------------------------------------
        //                      описание динамики цветения побега
        //------------------------------------------------------------------------------------------
        shared<curve<3>> REGRESSION(start1, height); // кривая старта цветения оси 1-го порядка (относительно начала цветения побега)
        shared<curve<3>> REGRESSION(start2, height); // кривая старта цветения оси 2-го порядка (относительно начала цветения терминалов)
        shared<curve<3>> REGRESSION(finish2, height); // кривая конца цветения оси 2-го порядка (относительно начала цветения терминалов)

        shared<curve<3>> CURVE(flowering1, time); // кривая цветения 1-го порядка
        shared<curve<3>> CURVE(flowering2, time); // кривая цветения 2-го порядка
        shared<curve<3>> CURVE(flowering, time); // кривая цветения 1-го + 2-го порядков

        typedef std::pair<Unsigned, gnuplotParams> view_type; // тип объекта, описывающего ту или иную зависимость и опции ее отрисовки
        std::map<Unsigned, gnuplotParams> dependencies; // список всех сформированных зависимостей
            // имя зависимости лежит в title of gnuplotParams

        Interval<double, 2> get_interval(shared<curve<2>> f) const { return {f->front()[0], f->back()[0]}; }
        Interval<double, 2> get_interval(shared<curve<3>> f) const { return {f->front()[0], f->back()[0]}; }

    public:

        typedef void (Sproutmodel::*PF)(const Flowerdome &, const std::set<Unsigned> &, double);
        std::map<Unsigned, std::pair<gnuplotParams, PF>> Xdependencies; // список всех сформированных зависимостей
            // имя зависимости лежит в title of gnuplotParams

        std::string get_info() const { return info; }

        //-----------------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------------
        std::vector<dependence_item> dependencies_items;

    #define ADD_SPECTRUM_LINES \
        shared<curve<2>> spectrum(new curve<2>); \
        for (auto e : *data) spectrum->push_back(e); \
        sort(spectrum, 0); \
        unite(spectrum, 0.001); \
        item.ddata = spectrum;

    #define PLOT_DISTRIBUTION \
        gnuparam.insert_data(f2, "u 1:2 w l lw 5 lc rgb 'blue' ");

    #define PLOT_SPECTRUM_LINES \
        flagged_value<double> max1 = find_maximal(item.smooth, 1); \
        flagged_value<double> max2 = find_maximal(item.ddata, 1); \
        assert(max1.flag() && max2.flag()); \
        double norm = max1 / max2; \
        norm *= 0.8; \
        for (auto &e : *item.ddata) e[1] *= norm; \
        gnuparam.insert_data(item.ddata, "u 1:(0):(0):2 w vectors nohead lw 5 lc rgb 'blue' ");

        //-----------------------------------------------------------------------------------------
        void shoots_by_flowers1(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            double hx = 0.1;
            double left_indent = 0.3;
            double right_indent = 0.3;

            auto data = flowerdome.foreach(sprouts, flower_count_functor(1));
            if (data->empty()) return;
            item.normality_test_n = data->size();

            ADD_SPECTRUM_LINES
            to_cumulative_density_function(data, hx, left_indent, right_indent);
            item.data = data;
        }

        void shoots_by_flowers1(dependence_item &item, double width)
        {
            item.smooth = distribution_smooth_with_errors(item.data, width);
            auto f = extract(item.smooth, 0, 1);
            auto f2 = cut_zerobounds(f);

            gnuplotParams &gnuparam = item.gnuparam;
            PLOT_DISTRIBUTION
            PLOT_SPECTRUM_LINES

            DISTRIBUTION(shoots, flowers1)->init(f2->begin(), f2->end(),
                [](auto v) { return std::make_pair(v[0], v[1]); }

            ); // инициализация распределения
        }

        std::string shoots_by_flowers1_test(dependence_item &, Unsigned n)
        {
            return normality_test(DISTRIBUTION(shoots, flowers1), n);
        }

        //-----------------------------------------------------------------------------------------
        auto flowers_by_phase_duration_time(Phase phase, dependence_item &item,
            const Flowerdome &flowerdome, const std::set<Unsigned> &sprouts)
        {
            double hx = 0.1;
            double left_indent = 0.3;
            double right_indent = 0.3;

            shared<curve<2>> data = flowerdome.foreach(sprouts, phase_duration_functor(phase));
            Unsigned data_n = data->size();

            if (!data->empty())
            {
                ADD_SPECTRUM_LINES
                to_cumulative_density_function(data, hx, left_indent, right_indent);
            }
            return std::make_pair(data_n, data);
        }

        void flowers_by_budding_time(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            auto [n, data] = flowers_by_phase_duration_time(BUD, item, flowerdome, sprouts);
            item.normality_test_n = n;
            item.data = data;
        }
        void flowers_by_budding_time(dependence_item &item, double width)
        {
            item.smooth = distribution_smooth_with_errors(item.data, width);
            auto f = extract(item.smooth, 0, 1);
            auto f2 = cut_zerobounds(f);

            gnuplotParams &gnuparam = item.gnuparam;
            PLOT_DISTRIBUTION
            PLOT_SPECTRUM_LINES

            DISTRIBUTION(flowers, budding_time)->init(f2->begin(), f2->end(),
                [](auto v) { return std::make_pair(v[0], v[1]); }
            );
        }

        std::string flowers_by_budding_time_test(dependence_item &, Unsigned n)
        {
            return normality_test(DISTRIBUTION(flowers, budding_time), n);
        }

        //-----------------------------------------------------------------------------------------
        void flowers_by_flowering_time(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            auto [n, data] = flowers_by_phase_duration_time(FLOWER, item, flowerdome, sprouts);
            item.normality_test_n = n;
            item.data = data;
        }

        void flowers_by_flowering_time(dependence_item &item, double width)
        {
            item.smooth = distribution_smooth_with_errors(item.data, width);
            auto f = extract(item.smooth, 0, 1);
            auto f2 = cut_zerobounds(f);

            gnuplotParams &gnuparam = item.gnuparam;
            PLOT_DISTRIBUTION
            PLOT_SPECTRUM_LINES

            DISTRIBUTION(flowers, flowering_time)->init(f2->begin(), f2->end(),
                [](auto v) { return std::make_pair(v[0], v[1]); }
            );
        }

        std::string flowers_by_flowering_time_test(dependence_item &, Unsigned n)
        {
            return normality_test(DISTRIBUTION(flowers, flowering_time), n);
        }

        //-----------------------------------------------------------------------------------------
        void flowers_by_ovaring_time(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            auto [n, data] = flowers_by_phase_duration_time(OVARY, item, flowerdome, sprouts);
            item.normality_test_n = n;
            item.data = data;
        }

        void flowers_by_ovaring_time(dependence_item &item, double width)
        {
            for (auto &e : *item.data) e[1] *= 5;
            item.smooth = distribution_smooth_with_errors(item.data, width);
            for (auto &e : *item.data) e[1] /= 5;
            for (auto &e : *item.smooth) e[1] /= 5;

            auto f = extract(item.smooth, 0, 1);
            auto f2 = cut_zerobounds(f);

            gnuplotParams &gnuparam = item.gnuparam;
            PLOT_DISTRIBUTION
            PLOT_SPECTRUM_LINES

            DISTRIBUTION(flowers, ovaring_time)->init(f2->begin(), f2->end(),
                [](auto v) { return std::make_pair(v[0], v[1]); }
            );
        }

        std::string flowers_by_ovaring_time_test(dependence_item &, Unsigned n)
        {
            return normality_test(DISTRIBUTION(flowers, ovaring_time), n);
        }

        //-----------------------------------------------------------------------------------------
        void flowers_by_between_time(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            shared<curve<2>> data(new curve<2>);
            for (auto sprout : sprouts)
            {
                auto tmp = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_START));
                if (tmp->empty()) continue;

                double search_eps = 0.01;
                std::sort(tmp->begin(), tmp->end(), [search_eps](auto a, auto b)
                    {
                        if (equal(a[0], b[0], search_eps)) return a[1] < b[1];
                        return a[0]< b[0];
                    }
                ); // сортировка для извлечения разниц

                data = data + to_beetween(tmp, search_eps);
            }
            item.data = data;
            item.normality_test_n = item.data->size();

            double hx = 0.1;
            double left_indent = 0.3;
            double right_indent = 0.3;

            ADD_SPECTRUM_LINES
            to_cumulative_density_function(item.data, hx, left_indent, right_indent);
        }

        void flowers_by_between_time(dependence_item &item, double width)
        {
            item.smooth = distribution_smooth_with_errors(item.data, width);
            auto f = extract(item.smooth, 0, 1);
            auto f2 = cut_zerobounds(f);

            gnuplotParams &gnuparam = item.gnuparam;
            PLOT_DISTRIBUTION
            PLOT_SPECTRUM_LINES
            DISTRIBUTION(flowers, between_time)->init(f2->begin(), f2->end(),
                [](auto v) { return std::make_pair(v[0], v[1]); }
            );
        }

        std::string flowers_by_between_time_test(dependence_item &, Unsigned n)
        {
            return normality_test(DISTRIBUTION(flowers, between_time), n);
        }

        //-----------------------------------------------------------------------------------------
        void flowers2_by_flowers1(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            item.data =  flowerdome.foreach(sprouts, flower_rate_functor(1, 2));
        }

        void flowers2_by_flowers1(dependence_item &item, double width)
        {
            item.smooth = smooth_regression(item.data, width);

            gnuplotParams &gnuparam = item.gnuparam;
            gnuparam.key = "top left";
            gnuparam.insert_data(extract_errors_range(item.smooth),
                "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder title 'δf' ");
            gnuparam.insert_data(extract_mean(item.smooth), "u 1:2 w l lw 5 lc 'blue' title 'f' ");
            gnuparam.insert_data(item.data, "u 1:2 w p pt 7 lc rgb 'blue' title '' ");

            REGRESSION(flowers2, flowers1) = item.smooth;
        }

        std::string flowers2_by_flowers1_test(dependence_item &item, Unsigned)
        {
            return regression_test(item);
        }

        //-----------------------------------------------------------------------------------------

        void density12_by_height(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            double hx = 0.01;

            shared<curve<2>> data(new curve<2>);  // 2nd order
            for (auto sprout : sprouts)
            {
                const Sprout &s = flowerdome.at(sprout);
                Unsigned nf1 = s.number_of_axis1_flowers();

                auto tmp = flowerdome.foreach({sprout}, axis_density_functor(2));
                if (tmp->empty()) continue;

                unite(tmp, hx);
                to_height(tmp, nf1);

                for (auto e : *tmp) data->push_back(e);
            }
            item.data = data;

            shared<curve<2>> data1(new curve<2>); // 1st order
            for (auto sprout : sprouts)
            {
                const Sprout &s = flowerdome.at(sprout);
                Unsigned nf1 = s.number_of_axis1_flowers();

                auto tmp = flowerdome.foreach({sprout}, axis_density_functor(1));
                if (tmp->empty()) continue;

                unite(tmp, hx);
                to_height(tmp, nf1);

                for (auto e : *tmp) data1->push_back(e);
            }
            item.data1 = data1;
        }

        void density12_by_height(dependence_item &item, double width)
        {
            item.smooth = smooth_regression(item.data, width);
            item.smooth1 = smooth_regression(item.data1, width);
            auto max = find_maximal(item.data, 1).value();

            gnuplotParams &gnuparam = item.gnuparam;
            gnuparam.title = "flowering density";
            gnuparam.key = "top right";
            gnuparam.set_xrange(0, max + 1);
            gnuparam.set_yrange(0, 1);

            gnuparam.insert_data(extract_errors_range(item.smooth1),
                "u 2:1 w filledcu closed fc 'red' fill transparent solid 0.2 noborder title 'δf_{ρ1}(h)' ");
            gnuparam.insert_data(extract_mean(item.smooth1), "u 2:1 w l lw 5 lc 'red' title 'f_{ρ1}(h)' ");
            gnuparam.insert_data(item.data1, "u 2:1 w p pt 7 lc rgb 'red' title '' ");

            gnuparam.insert_data(extract_errors_range(item.smooth),
                "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder title 'δf_{ρ2}(h)' ");
            gnuparam.insert_data(extract_mean(item.smooth), "u 2:1 w l lw 5 lc 'blue' title 'f_{ρ2}(h)' ");
            gnuparam.insert_data(item.data, "u 2:1 w p pt 7 lc rgb 'blue' title '' ");

            REGRESSION(density2, height) = item.smooth;
        }

        std::string density12_by_height_test(dependence_item &item, Unsigned)
        {
            return std::string("\n2nd order axis :\n") + regression_test(item);
        }
        //-----------------------------------------------------------------------------------------

        void start1_by_height(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            shared<curve<2>> data(new curve<2>);
            for (auto sprout : sprouts)
            {
                const Sprout &s = flowerdome.at(sprout);
                Unsigned nf1 = s.number_of_axis1_flowers();

                auto tmp = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
                if (tmp->empty()) continue;

                to_height(tmp, nf1);

                for (auto e : *tmp) data->push_back(e);
            }

            item.data = data;
        }

        void start1_by_height(dependence_item &item, double width)
        {
            item.smooth = smooth_regression(item.data, width);

            gnuplotParams &gnuparam = item.gnuparam;

            gnuparam.insert_data(extract_errors_range(item.smooth),
                "u 2:1 w filledcu closed fc 'red' fill transparent solid 0.2 noborder ");
            gnuparam.insert_data(extract_mean(item.smooth), "u 2:1 w l lw 5 lc 'red' ");
            gnuparam.insert_data(item.data, "u 2:1 w p pt 7 lc rgb 'red' ");

            REGRESSION(start1, height) = item.smooth;
        }

        std::string start1_by_height_test(dependence_item &item, Unsigned)
        {
            return regression_test(item);
        }

        void finish1_by_height(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            shared<curve<2>> data(new curve<2>);
            for (auto sprout : sprouts)
            {
                const Sprout &s = flowerdome.at(sprout);
                Unsigned nf1 = s.number_of_axis1_flowers();

                auto tmp = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_END));
                if (tmp->empty()) continue;

                to_height(tmp, nf1);

                for (auto e : *tmp) data->push_back(e);
            }
            item.data = data;
        }

        //-----------------------------------------------------------------------------------------
        void start2_by_height(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            shared<curve<2>> data(new curve<2>);
            for (auto sprout : sprouts)
            {
                const Sprout &s = flowerdome.at(sprout);
                Unsigned nf1 = s.number_of_axis1_flowers();

                auto tmp1 = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
                if (tmp1->empty()) continue;

                to_height(tmp1, nf1); // для получения относительного старта

                auto tmp2 = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_START));
                if (tmp2->empty()) continue;

                double search_eps = 0.01;
                double represent_part = 0.3;
                extract_bound(tmp2, true, search_eps, represent_part);

                to_height(tmp2, nf1);
                to_relative(tmp2, tmp1);

                for (auto e : *tmp2) data->push_back(e);
            }
            item.data = data;
        }

        void start2_by_height(dependence_item &item, double width)
        {
            item.smooth = smooth_regression(item.data, width);

            gnuplotParams &gnuparam = item.gnuparam;

            gnuparam.insert_data(extract_errors_range(item.smooth),
                "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ");
            gnuparam.insert_data(extract_mean(item.smooth), "u 2:1 w l lw 5 lc 'blue' ");
            gnuparam.insert_data(item.data, "u 2:1 w p pt 7 lc rgb 'blue' ");

            REGRESSION(start2, height) = item.smooth;
        }

        std::string start2_by_height_test(dependence_item &item, Unsigned)
        {
            return regression_test(item);
        }

        void start12_by_height(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            {
                shared<curve<2>> data(new curve<2>);
                for (auto sprout : sprouts)
                {
                    const Sprout &s = flowerdome.at(sprout);
                    Unsigned nf1 = s.number_of_axis1_flowers();

                    auto tmp = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
                    if (tmp->empty()) continue;

                    to_height(tmp, nf1);

                    for (auto e : *tmp) data->push_back(e);

                    Time ts = s.day_of_first_start(FLOWER);
                    for (auto &e : *data) e[1] = e[1] - ts;
                }

                item.data1 = data;
            }
            {
                shared<curve<2>> data(new curve<2>);
                for (auto sprout : sprouts)
                {
                    const Sprout &s = flowerdome.at(sprout);
                    Unsigned nf1 = s.number_of_axis1_flowers();

                    auto tmp1 = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
                    if (tmp1->empty()) continue;

                    to_height(tmp1, nf1); // для получения относительного старта

                    auto tmp2 = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_START));
                    if (tmp2->empty()) continue;

                    double search_eps = 0.01;
                    double represent_part = 0.3;
                    extract_bound(tmp2, true, search_eps, represent_part);

                    to_height(tmp2, nf1);
                    to_relative(tmp2, tmp1);

                    for (auto e : *tmp2) data->push_back(e);
                }
                item.data2 = data;
            }
        }

        void start12_by_height(dependence_item &item, double width)
        {
            {
                item.smooth1 = smooth_regression(item.data1, width);

                gnuplotParams &gnuparam = item.gnuparam;
                gnuparam.key = "top right";
                gnuparam.title = "flowering patterns";

                gnuparam.insert_data(extract_errors_range(item.smooth1),
                    "u 2:1 w filledcu closed fc 'red' fill transparent solid 0.2 noborder title 'δf_{t1}(h)' ");
                gnuparam.insert_data(extract_mean(item.smooth1), "u 2:1 w l lw 5 lc 'red' title 'f_{t1}(h)' ");
                gnuparam.insert_data(item.data1, "u 2:1 w p pt 7 lc rgb 'red' title '' ");

                REGRESSION(start1, height) = item.smooth1;
            }

            {
                item.smooth2 = smooth_regression(item.data2, width);
                gnuplotParams &gnuparam = item.gnuparam;

                gnuparam.insert_data(extract_errors_range(item.smooth2),
                    "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder title 'δf_{t2}(h)' ");
                gnuparam.insert_data(extract_mean(item.smooth2), "u 2:1 w l lw 5 lc 'blue' title 'f_{t2}(h)' ");
                gnuparam.insert_data(item.data2, "u 2:1 w p pt 7 lc rgb 'blue' title '' ");

                REGRESSION(start2, height) = item.smooth2;
            }
        }

        std::string start12_by_height_test(dependence_item &item, Unsigned)
        {
            std::string msg1 = "\n1-st order axis: \n" + regression_test(item.data1, item.smooth1, item);
            std::string msg2 = "2-nd order axis: \n" + regression_test(item.data2, item.smooth2, item);
            return msg1 + '\n' + msg2;
        }

        //-----------------------------------------------------------------------------------------
        void finish2_by_height(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            shared<curve<2>> data(new curve<2>);
            for (auto sprout : sprouts)
            {
                const Sprout &s = flowerdome.at(sprout);
                Unsigned nf1 = s.number_of_axis1_flowers();

                auto tmp1 = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
                if (tmp1->empty()) continue;

                to_height(tmp1, nf1); // для получения относительног старта

                auto tmp2 = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_END));
                if (tmp2->empty()) continue;

                double search_eps = 0.01;
                double represent_part = 0.3;
                extract_bound(tmp2, false, search_eps, represent_part);

                to_height(tmp2, nf1);
                to_relative(tmp2, tmp1);

                for (auto e : *tmp2) data->push_back(e);
            }
            item.data = data;
        }

        void finish2_by_height(dependence_item &item, double width)
        {
            item.smooth = smooth_regression(item.data, width);

            gnuplotParams &gnuparam = item.gnuparam;

            gnuparam.insert_data(extract_errors_range(item.smooth),
                "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ");
            gnuparam.insert_data(extract_mean(item.smooth), "u 2:1 w l lw 5 lc 'blue' ");
            gnuparam.insert_data(item.data, "u 2:1 w p pt 7 lc rgb 'blue' ");

            REGRESSION(finish2, height) = item.smooth;
        }

        std::string finish2_by_height_test(dependence_item &item, Unsigned)
        {
            return regression_test(item);
        }

        void abs_start2_by_height(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            shared<curve<2>> data(new curve<2>);
            for (auto sprout : sprouts)
            {
                const Sprout &s = flowerdome.at(sprout);
                Unsigned nf1 = s.number_of_axis1_flowers();

                auto tmp2 = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_START));
                if (tmp2->empty()) continue;

                double search_eps = 0.01;
                double represent_part = 0.3;
                extract_bound(tmp2, true, search_eps, represent_part);

                to_height(tmp2, nf1);
                for (auto e : *tmp2) data->push_back(e);
            }
            item.data = data;
        }

        void abs_finish2_by_height(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            shared<curve<2>> data(new curve<2>);
            for (auto sprout : sprouts)
            {
                const Sprout &s = flowerdome.at(sprout);
                Unsigned nf1 = s.number_of_axis1_flowers();

                auto tmp2 = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_END));
                if (tmp2->empty()) continue;

                double search_eps = 0.01;
                double represent_part = 0.3;
                extract_bound(tmp2, true, search_eps, represent_part);

                to_height(tmp2, nf1);
                for (auto e : *tmp2) data->push_back(e);
            }
            item.data = data;
        }

        //-----------------------------------------------------------------------------------------
        shared<curve<2>> development_by_time(extraction_phase_of_axis fn, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            shared<curve<2>> data(new curve<2>);

            Time tmin = 365.; Time tmax = 0.; bool flag = false;
            for (auto sprout : sprouts)
            {
                auto data = flowerdome.foreach({sprout}, fn);
                    // data = { {height, start, finish}, {..} ... }

                if (data->empty()) continue;

                auto min = find_minimal(data, 1); // получить дату самого раннего старта
                if (tmin > min.value()) tmin = min.value();

                auto max = find_maximal(data, 2); // получить дату самого позднего финиша
                if (tmax < max.value()) tmax = max.value();
                flag = true;
            }
            if (!flag) return data;

            Unsigned left_indent = 5;
            Unsigned right_indent = 5;

            Signed min_ndx = tmin.to_signed() - left_indent;
            Signed max_ndx = tmax.to_signed() + right_indent;

            Unsigned count = max_ndx - min_ndx + 1;
            prg::vector<double> sums(min_ndx, count, 0.);

            for (auto sprout : sprouts)
            {
                auto data = flowerdome.foreach({sprout}, fn);
                if (data->empty()) continue;

                for (auto e : *data)
                {
                    Signed t1 = e[1];
                    Signed t2 = e[2];
                    for (auto t=t1; t<=t2; t++) sums[t] += 1;
                }
            }

            for (Signed t=sums.minimal_index(); t<=sums.maximal_index(); t++)
                data->push_back({double(t), sums[t]});

            return data;
        }

        void flowering1_by_time(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            Time tstart = 365.; // starting date
            for (auto sprout : sprouts)
            {
                auto data = flowerdome.foreach({sprout}, extraction_phase_of_axis(1, FLOWER, 0.));
                if (data->empty()) continue;

                auto min = find_minimal(data, 1); // получить дату самого раннего старта
                if (tstart > min.value()) tstart = min.value();
            }

            item.data = development_by_time(extraction_phase_of_axis(1, FLOWER, tstart), flowerdome, sprouts);
            cut_below(item.data);
        }

        void flowering1_by_time(dependence_item &item, double width)
        {
            item.smooth = smooth_regression(item.data, width);

            auto mean = extract_mean(item.smooth);
            cut_below(mean);
            auto range = extract_errors_range(item.smooth);
            cut_below(range);

            gnuplotParams &gnuparam = item.gnuparam;

            gnuparam.insert_data(range,
                "u 1:2 w filledcu closed fc 'red' fill transparent solid 0.2 noborder ");
            gnuparam.insert_data(mean, "u 1:2 w l lw 5 lc 'red' ");
            gnuparam.insert_data(item.data, "u 1:2 w p pt 7 lc rgb 'red' ");

            CURVE(flowering1, time) = item.smooth;
        }

        std::string flowering1_by_time_test(dependence_item &item, Unsigned)
        {
            return regression_test(item);
        }

        //-----------------------------------------------------------------------------------------
        void flowering2_by_time(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            Time tstart = 365.; // date of beginning of flowering for the sprout
            for (auto sprout : sprouts)
            {
                auto data = flowerdome.foreach({sprout}, extraction_phase_of_axis(1, FLOWER, 0.));
                if (data->empty()) continue;

                auto min = find_minimal(data, 1); // получить дату самого раннего старта
                if (tstart > min.value()) tstart = min.value();
            }

            item.data = development_by_time(extraction_phase_of_axis(2, FLOWER, tstart), flowerdome, sprouts);
            cut_below(item.data);
        }

        void flowering2_by_time(dependence_item &item, double width)
        {
            item.smooth = smooth_regression(item.data, width);

            auto mean = extract_mean(item.smooth);
            cut_below(mean);
            auto range = extract_errors_range(item.smooth);
            cut_below(range);

            gnuplotParams &gnuparam = item.gnuparam;

            gnuparam.insert_data(range,
                "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ");
            gnuparam.insert_data(mean, "u 1:2 w l lw 5 lc 'blue' ");
            gnuparam.insert_data(item.data, "u 1:2 w p pt 7 lc rgb 'blue' ");

            CURVE(flowering2, time) = item.smooth;
        }

        std::string flowering2_by_time_test(dependence_item &item, Unsigned)
        {
            return regression_test(item);
        }

        void flowering_by_time(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            flowering1_by_time(item, flowerdome, sprouts);
            item.data1 = item.data;

            flowering2_by_time(item, flowerdome, sprouts);
            item.data2 = item.data;

            item.data = summarize<2>(item.data1, item.data2);
            item.ddata = subtract<2>(item.data1, item.data2);
        }

        void flowering_by_time(dependence_item &item, double width)
        {
            item.smooth1 = smooth_regression(item.data1, width);
            item.smooth2 = smooth_regression(item.data2, width);

//            item.smooth  = summarize<3>(item.smooth1, item.smooth2);
//            double eps = (item.smooth->back()[0] - item.smooth->front()[0]) / item.smooth->size();
//            item.dsmooth  = subtract<3>(item.smooth1, item.smooth2, eps);
            // produces artefacts

            item.smooth  = smooth_regression(item.data , width);
            item.dsmooth = smooth_regression(item.ddata, width);

            auto mean1 = extract_mean(item.smooth1);
            cut_below(mean1);

            double x = first_nonzero(mean1);
            auto mean2 = extract_mean(item.smooth2);
            cut_below(mean2);
            auto mean = extract_mean(item.smooth);
            cut_below(mean);
            auto dmean = extract_mean(item.dsmooth);
            cut_below(dmean, -100, x);

            auto range1 = extract_errors_range(item.smooth1);
            cut_below(range1);
            auto range2 = extract_errors_range(item.smooth2);
            cut_below(range2);

            gnuplotParams &gnuparam = item.gnuparam;
            gnuparam.key = "top right";
            gnuparam.title = "flowering curves";

            gnuparam.insert_data(range1, "u 1:2 w filledcu closed fc 'red' fill transparent solid 0.2 noborder title 'δf_1(t)' ");
            gnuparam.insert_data(mean1, "u 1:2 w l lw 5 lc 'red' title 'f_1(t)' ");
            gnuparam.insert_data(range2, "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder title 'δf_2(t)' ");
            gnuparam.insert_data(mean2, "u 1:2 w l lw 5 lc 'blue' title 'f_2(t)' ");

            gnuparam.insert_data(item.data1, "u 1:2 w p pt 7 lc rgb 'red' title '' ");
            gnuparam.insert_data(item.data2, "u 1:2 w p pt 7 lc rgb 'blue' title '' ");
            //gnuparam.insert_data(item.data, "u 1:2 w p pt 7 lc rgb 'gray' title '' ");
            gnuparam.insert_data(mean, "u 1:2 w l lw 5 lc 'black' title 'f(t)' ");
            gnuparam.insert_data(dmean, "u 1:2 w l lt 0 lw 5 lc 'black' title 'Δf(t)' ");

            CURVE(flowering, time) = item.smooth;
        }

        std::string flowering_by_time_test(dependence_item &item, Unsigned)
        {
            return regression_test(item);
        }

        void flowering_ranges(dependence_item &item, const Flowerdome &flowerdome,
            const std::set<Unsigned> &sprouts)
        {
            Time ts = flowerdome.day_of_fist_start(*sprouts.begin(), FLOWER);

            dependence_item start1_item;
            start1_by_height(start1_item, flowerdome, sprouts);
            for (auto &e : *start1_item.data) e[1] = e[1] - ts;
            item.data1 = start1_item.data; // используем для записи start1

            dependence_item finish1_item;
            finish1_by_height(finish1_item, flowerdome, sprouts);
            for (auto &e : *finish1_item.data) e[1] = e[1] - ts;
            item.data = finish1_item.data; // используем для записи finish1
            assert(item.data1->size() == item.data->size());

            dependence_item start2_item;
            abs_start2_by_height(start2_item, flowerdome, sprouts);
            for (auto &e : *start2_item.data) e[1] = e[1] - ts;
            item.data2 = start2_item.data; // используем для записи start2

            dependence_item finish2_item;
            abs_finish2_by_height(finish2_item, flowerdome, sprouts);
            for (auto &e : *finish2_item.data) e[1] = e[1] - ts;
            item.ddata = finish2_item.data; // используем для записи finish2
            assert(item.data2->size() == item.ddata->size());
        }

        void flowering_ranges(dependence_item &item, double width)
        {
            item.smooth1 = smooth_regression(item.data1, width); // start1
            item.smooth  = smooth_regression(item.data , width); // end1
            item.smooth2 = smooth_regression(item.data2, width); // start2
            item.dsmooth = smooth_regression(item.ddata, width); // end2

            auto mean1 = extract_mean(item.smooth1); // start1
            auto mean = extract_mean(item.smooth); // end1
            auto mean2 = extract_mean(item.smooth2); // start2
            auto dmean = extract_mean(item.dsmooth); // end2

            auto range1 = mean1;
            for (auto it = mean->rbegin(); it != mean->rend(); ++it) range1->push_back(*it);

            auto range2 = mean2;
            for (auto it = dmean->rbegin(); it != dmean->rend(); ++it) range2->push_back(*it);

            gnuplotParams &gnuparam = item.gnuparam;
            gnuparam.title = "flowering ranges";
            gnuparam.insert_data(range1, "u 2:1 w filledcu closed fc 'red' fill transparent solid 0.2 noborder title 'D_1(τ,h)' ");
            gnuparam.insert_data(range2, "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder title 'D_2(τ,h)' ");

            shared<curve<4>> data1(new curve<4>);
            assert(item.data1->size() == item.data->size());
            Unsigned n = item.data1->size();

            for (Unsigned i=0; i<n; i++)
            {
                double h1 = (*item.data1)[i][0];
                double t1 = (*item.data1)[i][1];
                double h2 = (*item.data)[i][0];
                double t2 = (*item.data)[i][1];
                assert(equal(h1, h2, 0.001));
                data1->push_back({h1, t1, 0, t2 - t1});
            }
            gnuparam.insert_data(data1, "u 2:1:4:3 w vectors lw 4 lc rgb 'red' title 'τ_1(h)' ");

            shared<curve<4>> data2(new curve<4>);
            assert(item.data2->size() == item.ddata->size());
            n = item.data2->size();

            for (Unsigned i=0; i<n; i++)
            {
                double h1 = (*item.data2)[i][0];
                double t1 = (*item.data2)[i][1];
                double h2 = (*item.ddata)[i][0];
                double t2 = (*item.ddata)[i][1];
                assert(equal(h1, h2, 0.001));
                data2->push_back({h1, t1, 0, t2 - t1});
            }
            gnuparam.insert_data(data2, "u 2:1:4:3 w vectors lw 4 lc rgb 'blue' title 'τ_2(h)' ");
            item.info = "no statistics information";
        }

        std::string flowering_ranges_test(dependence_item &, Unsigned)
        {
            return "\nno statistics information";
        }

        //-----------------------------------------------------------------------------------------

        /// получить случайное значение в точке x около заданной кривой
        double generate_at(shared<curve<3>> f, double x) const
        {
            typedef std::uniform_real_distribution<> distribution_type;
                // нормальное распределение дает большие выбросы за пределы области

            distribution_type dist;
            auto p = find_nearest(f, x);

            auto flag = p.flag();
            if (flag == false) return 0;

            auto v = p.value();

            distribution_type::param_type param(v[1] - v[2], v[1] + v[2]);
            dist.param(param);

            return dist(rangen);
        }

    public:

        /// генерить массив событий, содержащий число цветков 1-го порядка
        double generate_flower1() const
        {
            double v = (*DISTRIBUTION(shoots, flowers1))(randgen());
            return v < 1 ? 1 : v;
        }

        /// генерить массив событий, содержащий число цветков 1-го порядка
        double generate_between() const
        {
            double len = std::round((*DISTRIBUTION(flowers, between_time))(randgen()));
            return len < 0 ? 0 : len;
        }

        /// генерить массив событий, содержащий число цветков 1-го порядка
        double generate_duration(Phase phase) const
        {
            double len = 0;
            if (phase == BUD) len = (*DISTRIBUTION(flowers, budding_time))(randgen());
            if (phase == FLOWER) len = (*DISTRIBUTION(flowers, flowering_time))(randgen());
            if (phase == OVARY) len = (*DISTRIBUTION(flowers, ovaring_time))(randgen());
            len = std::round(len);
            return len < 1 ? 1 : len;
        }

        /// сгенерить массив числа цветков 1-го порядка (для тестирования)
        shared<curve<2>> generate_flower1(Unsigned n) const
        {
            shared<curve<2>> out(new curve<2>);
            for (Unsigned i=0; i<n; i++)
            {
                double v = generate_flower1();
                out->push_back({std::ceil(v), 1});
            }
            return out;
        }

        /// генерить массив, содержащий число цветков 2-го порядка для соответствующего массива числа цветков 1-го порядка
        double generate_flower2(double n1) const
        {
            double v = generate_at(REGRESSION(flowers2, flowers1), n1);
            return v < 0 ? 0 : v;
        }

        /// генерить массив числа цветков 2-го порядка (для тестирования)
        shared<curve<2>> generate_flower2(shared<curve<2>> n1) const
        {
            shared<curve<2>> out(new curve<2>);
            for (auto e : *n1)
            {
                double nx = e[0];
                double v = generate_flower2(nx);
                out->push_back({nx, v});
            }
            return out;
        }

        /// сгенерить плотность цветков на заданной высоте
        double generate_density2(double h) const { return generate_at(REGRESSION(density2, height), h); }

        /// сгенерить сортированный список высот для заданного числа цветков (1|2)-го порядка
        shared<curve<1>> generate_heights(double nf1, Unsigned axis)
        {
            shared<curve<1>> out(new curve<1>);
            for (Unsigned i=nf1; i>=1; i--) out->push_back(relative_height(i, nf1));

            if (axis == 2)
            {
                auto hrange = get_interval(REGRESSION(density2, height));
                to_range(out, hrange, 0); // отрисовка только != 1 значений
             }
            return out;
        }

        /// сгенерить плотность цветков по списку высот и полному числу цветков 2-го порядка
        shared<curve<2>> generate_density2(shared<curve<1>> heights, double nf2) const
        {
            shared<curve<2>> out(new curve<2>);
            double sum = 0;
            for (auto e : *heights)
            {
                double h = e[0];
                double v = generate_density2(h);
                out->push_back({h, v});
                sum += v;
            }
            if (sum != 0)
            {
                to_normalized(out, nf2);
                round(out);
                auto hrange = get_interval(REGRESSION(density2, height));
                to_range(out, hrange, 0); // отрисовка только != 1 значений
            }

        #ifdef TEST_PLOT
            auto f3 = smooth_regression(out, 0.2);
            gnuplotParams param;
            param.xlabel = "number of the 2nd-order flowers";
            param.ylabel = "relative high";
            param.insert_data(extract_errors_range(f3), "u 2:1 w filledcu closed");
            param.insert_data(extract_mean(f3), "u 2:1 w l");
            param.insert_data(out, "u 2:1 w p");
            plot(param);
        #endif

            return out;
        }

        /// сгенерить дату старта цветения цветков 1-го порядка на заданной высоте
        double generate_start1(double h) const { return generate_at(REGRESSION(start1, height), h); }

        /// сгенерить кривую старта цветения цветков 1-го порядка
        shared<curve<2>> generate_start1(shared<curve<1>> heights) const
        {
            shared<curve<2>> tmp1(new curve<2>);
            for (auto e : *heights)
            {
                double h = e[0];
                double v = generate_start1(h);
                tmp1->push_back({h, v});
            }
            round(tmp1);

            Unsigned steps = 500;
            auto out = minimize_dispersion(tmp1, REGRESSION(start1, height), steps);

        #ifdef TEST_PLOT
            auto f3 = smooth_regression(out, 0.2);
            gnuplotParams param;
            param.xlabel = "start of the 1st-order axis flowering";
            param.ylabel = "relative high";
            param.insert_data(extract_errors_range(f3), "u 2:1 w filledcu closed");
            param.insert_data(extract_mean(f3), "u 2:1 w l");
            param.insert_data(out, "u 2:1 w p");
            plot(param);
        #endif

            return out;
        }

        /// сгенерить дату старта цветения цветков 1-го порядка на заданной высоте
        double generate_start2(double h) const { return generate_at(REGRESSION(start2, height), h); }

        /// сгенерить кривую старта цветения цветков 2-го порядка
        shared<curve<2>> generate_start2(shared<curve<1>> heights) const
        {
            shared<curve<2>> tmp(new curve<2>);
            for (auto e : *heights)
            {
                double h = e[0];
                double v = generate_start2(h);
                tmp->push_back({h, v});
            }
            round(tmp);

            Unsigned steps = 100;
            auto out = minimize_dispersion(tmp, REGRESSION(start2, height), steps);

        #ifdef TEST_PLOT
            auto f3 = smooth_regression(out, 0.2);
            gnuplotParams param;
            param.xlabel = "relative start of the 2nd-order axis flowering";
            param.ylabel = "relative high";
            param.insert_data(extract_errors_range(f3), "u 2:1 w filledcu closed");
            param.insert_data(extract_mean(f3), "u 2:1 w l");
            param.insert_data(out, "u 2:1 w p");
            plot(param);
        #endif

            return out;
        }

    public:

        Sproutmodel() {}

        /*!
         * \brief build the model
         * \param flowerdome база данных по которой строятся общие распределения по побегам
         * \param sprout номер побега, по которому строятся модель конкретного побега
         */
        void build(const Flowerdome &flowerdome, Unsigned sprout, const std::string &comment="");

        /*!
         * \brief print the list of all the distributions and regressions
         * \note адаптировано под вывод в ncurses terminal
         */
        std::string view_list() const;

        /// число зависимостей
        Unsigned view_list_size() const { return dependencies_items.size(); }

        /*!
         * \brief view the distributions and regressions
         * \param идентификатор зависимости
         */
        std::string view(Unsigned id);
        std::string view(const std::vector<Unsigned> &ids);

        /*!
         * \brief генерация случайного образца модельного побега
         * \param restrictions ограничения по числу цветков различных порядков
         *  [0] по полному числу цветков
         *  [1] по числу цветков на 1й оси
         *  [2] по числу цветков на каждой 2й оси ...
         */
        shared<Sprout> generate_sprout(const std::vector<Interval<Unsigned, 1>> &restrictions);
        shared<Sprout> generate_sprout(Unsigned n1=inf<Unsigned>(), Unsigned n2=inf<Unsigned>(),
            const std::string &flowering_start_date="1-Jul", Unsigned flowering_days=10);
    };

    inline std::string Sproutmodel::view_list() const
    {
        if (dependencies_items.empty())
            return "No data yet. Build the model [model sprout_id]";

        std::string msg;
        for (Unsigned ndx=1; ndx<dependencies_items.size(); ndx++)
        {
            Unsigned id = ndx;
            std::string title = dependencies_items[ndx].title;
            msg += '[' + std::to_string(id) + "] " + title + "\n";
        }
        return msg;
    }

    inline std::string Sproutmodel::view(Unsigned id)
    {
        if (dependencies_items.empty())
            return "No data yet. Build the model [model sprout_id]";

        auto item = dependencies_items.at(id);

        gnuplotParams v = item.gnuparam;
        plot(v);
        return item.info;
    }

    inline std::string Sproutmodel::view(const std::vector<Unsigned> &ids)
    {
        bool flag = true;
        for (auto id : ids) if (dependencies.count(id) == 0) flag = false;
        if (!flag) return "No data yet. Build the model [model sprout_id]";

        gnuplotParams params = dependencies.at(ids.front());
        for (Unsigned i=1; i<ids.size(); i++)
        {
            gnuplotParams v = dependencies.at(ids[i]);
            if (v.xlabel != params.xlabel) return "fail";
            if (v.ylabel != params.ylabel) return "fail";

            for (Unsigned j=0; j<v.size(); j++)
            {
                auto data = v.get_data(j);
                auto options = v.get_options(j);
                params.insert_data(data, options);
            }
        }
        plot(params);
        return "OK";
    }

    inline void Sproutmodel::build(const Flowerdome &flowerdome, Unsigned sprout,
        const std::string &comment)
    {
        info = comment; // save external information about model

        nflower1 = flowerdome.sprout(sprout).number_of_axis1_flowers();
        nflower2 = flowerdome.sprout(sprout).number_of_axis2_flowers();

        auto sprouts = flowerdome.get_sprouts();

        dependence_item item;
        gnuplotParams &gnuparam = item.gnuparam;
        gnuparam.output = TMP_DIR + "/terminal.pdf";

        dependencies_items.push_back(item);  // empty for [0]

        //------------------------------------------------------------------------------------------
        //                      dependencies for the all sprouts given
        //------------------------------------------------------------------------------------------
        item.title = "distribution of shoots by the number of 1st-order flowers";
        item.correlation_coef = false;
        item.determination_coef = false;
        item.normality_test = true;
        item.pf1 = &Sproutmodel::shoots_by_flowers1;
        item.pf2 = &Sproutmodel::shoots_by_flowers1;
        item.pf3 = &Sproutmodel::shoots_by_flowers1_test;
        gnuparam.clear();
        gnuparam.title = item.title;
        gnuparam.xlabel = "number of the 1st-order flowers";
        gnuparam.ylabel = "probability density";
        item.info = "";
        (this->*(item.pf1))(item, flowerdome, sprouts);
        (this->*(item.pf2))(item, default_width);
        dependencies_items.push_back(item);

        item.title = "distribution of flowers by duration of budding";
        item.correlation_coef = false;
        item.determination_coef = false;
        item.normality_test = true;
        item.pf1 = &Sproutmodel::flowers_by_budding_time;
        item.pf2 = &Sproutmodel::flowers_by_budding_time;
        item.pf3 = &Sproutmodel::flowers_by_budding_time_test;
        gnuparam.clear();
        gnuparam.title = item.title;
        gnuparam.xlabel = "budding duration";
        gnuparam.ylabel = "probability density";
        item.info = "";
        (this->*(item.pf1))(item, flowerdome, sprouts);
        (this->*(item.pf2))(item, default_width);
        dependencies_items.push_back(item);

        item.title = "distribution of flowers by duration of flowering";
        item.correlation_coef = false;
        item.determination_coef = false;
        item.normality_test = true;
        item.pf1 = &Sproutmodel::flowers_by_flowering_time;
        item.pf2 = &Sproutmodel::flowers_by_flowering_time;
        item.pf3 = &Sproutmodel::flowers_by_flowering_time_test;
        gnuparam.clear();
        gnuparam.title = item.title;
        gnuparam.xlabel = "flowering duration";
        gnuparam.ylabel = "probability density";
        item.info = "";
        (this->*(item.pf1))(item, flowerdome, sprouts);
        (this->*(item.pf2))(item, default_width);
        dependencies_items.push_back(item);

        item.title = "distribution of flowers by duration of ovaring";
        item.correlation_coef = false;
        item.determination_coef = false;
        item.normality_test = true;
        item.pf1 = &Sproutmodel::flowers_by_ovaring_time;
        item.pf2 = &Sproutmodel::flowers_by_ovaring_time;
        item.pf3 = &Sproutmodel::flowers_by_ovaring_time_test;
        gnuparam.clear();
        gnuparam.title = item.title;
        gnuparam.xlabel = "ovaring duration";
        gnuparam.ylabel = "probability density";
        item.info = "";
        (this->*(item.pf1))(item, flowerdome, sprouts);
        (this->*(item.pf2))(item, default_width);
        dependencies_items.push_back(item);

        item.title = "distribution of intervals";
        item.correlation_coef = false;
        item.determination_coef = false;
        item.normality_test = true;
        item.pf1 = &Sproutmodel::flowers_by_between_time;
        item.pf2 = &Sproutmodel::flowers_by_between_time;
        item.pf3 = &Sproutmodel::flowers_by_between_time_test;
        gnuparam.clear();
        gnuparam.title = item.title;
        gnuparam.xlabel = "duration between 2nd-order flower openings";
        gnuparam.ylabel = "probability density";
        item.info = "";
        (this->*(item.pf1))(item, flowerdome, sprouts);
        (this->*(item.pf2))(item, default_width);
        dependencies_items.push_back(item);

        item.title = "regression between the number of the 1st-order flowers and the 2nd-order flowers";
        item.correlation_coef = true;
        item.determination_coef = true;
        item.normality_test = true;
        item.pf1 = &Sproutmodel::flowers2_by_flowers1;
        item.pf2 = &Sproutmodel::flowers2_by_flowers1;
        item.pf3 = &Sproutmodel::flowers2_by_flowers1_test;
        gnuparam.clear();
        gnuparam.title = item.title;
        gnuparam.xlabel = "number of the 1st-order flowers";
        gnuparam.ylabel = "number of the 2nd-order flowers";
        item.info = "";
        (this->*(item.pf1))(item, flowerdome, sprouts);
        (this->*(item.pf2))(item, default_width);
        dependencies_items.push_back(item);

        //------------------------------------------------------------------------------------------
        //              dependencies for the model of sprout given but not for all sprouts
        //------------------------------------------------------------------------------------------
//        item.title = "the number of flowers on the axes of the 2nd order, depending on the relative height of the axis";
//        item.correlation_coef = false;
//        item.determination_coef = true;
//        item.normality_test = true;
//        item.pf1 = &Sproutmodel::density2_by_height;
//        item.pf2 = &Sproutmodel::density2_by_height;
//        item.pf3 = &Sproutmodel::density2_by_height_test;
//        gnuparam.clear();
//        gnuparam.title = item.title;
//        gnuparam.xlabel = "number of the 2nd-order flowers";
//        gnuparam.ylabel = "h, relative height";
//        item.info = "";
//        (this->*(item.pf1))(item, flowerdome, {sprout});
//        (this->*(item.pf2))(item, default_width);
//        dependencies_items.push_back(item);

        item.title = "the number of flowers on axes, depending on height";
        item.correlation_coef = false;
        item.determination_coef = true;
        item.normality_test = true;
        item.pf1 = &Sproutmodel::density12_by_height;
        item.pf2 = &Sproutmodel::density12_by_height;
        item.pf3 = &Sproutmodel::density12_by_height_test;
        gnuparam.clear();
        gnuparam.title = item.title;
        gnuparam.xlabel = "number of flowers";
        gnuparam.ylabel = "h, relative height";
        item.info = "";
        (this->*(item.pf1))(item, flowerdome, {sprout});
        (this->*(item.pf2))(item, default_width);
        dependencies_items.push_back(item);

//        item.title = "pattern of flowering of the 1st-order axis";
//        item.correlation_coef = false;
//        item.determination_coef = true;
//        item.normality_test = true;
//        item.pf1 = &Sproutmodel::start1_by_height;
//        item.pf2 = &Sproutmodel::start1_by_height;
//        item.pf3 = &Sproutmodel::start1_by_height_test;
//        gnuparam.clear();
//        gnuparam.title = item.title;
////        gnuparam.xlabel = "start of the 1st-order axis flowering";
//        gnuparam.xlabel = "T, days from the beginning of flowering";
//        gnuparam.ylabel = "h, relative height";
//        item.info = "";
//        (this->*(item.pf1))(item, flowerdome, {sprout});
//        (this->*(item.pf2))(item, default_width);
//        dependencies_items.push_back(item);

//        item.title = "delays in the beginning of flowering of axes of the 2nd order";
//        item.correlation_coef = false;
//        item.determination_coef = true;
//        item.normality_test = true;
//        item.pf1 = &Sproutmodel::start2_by_height;
//        item.pf2 = &Sproutmodel::start2_by_height;
//        item.pf3 = &Sproutmodel::start2_by_height_test;
//        gnuparam.clear();
//        gnuparam.title = item.title;
//        gnuparam.xlabel = "relative start of the 2nd-order axis flowering";
//        gnuparam.ylabel = "h, relative height";
//        item.info = "";
//        (this->*(item.pf1))(item, flowerdome, {sprout});
//        (this->*(item.pf2))(item, default_width);
//        dependencies_items.push_back(item);

        item.title = "patterns of flowering of the 1st-order and 2-nd order axes";
        item.correlation_coef = false;
        item.determination_coef = true;
        item.normality_test = true;
        item.pf1 = &Sproutmodel::start12_by_height;
        item.pf2 = &Sproutmodel::start12_by_height;
        item.pf3 = &Sproutmodel::start12_by_height_test;
        gnuparam.clear();
        gnuparam.title = item.title;
        gnuparam.xlabel = "T, days from the beginning of flowering";
        gnuparam.ylabel = "h, relative height";
        item.info = "";
        (this->*(item.pf1))(item, flowerdome, {sprout});
        (this->*(item.pf2))(item, default_width);
        dependencies_items.push_back(item);



//        item.title = "delays in the finishing of flowering of axes of the 2nd order";
//        item.correlation_coef = false;
//        item.determination_coef = true;
//        item.normality_test = true;
//        item.pf1 = &Sproutmodel::finish2_by_height;
//        item.pf2 = &Sproutmodel::finish2_by_height;
//        item.pf3 = &Sproutmodel::finish2_by_height_test;
//        gnuparam.clear();
//        gnuparam.title = item.title;
//        gnuparam.xlabel = "relative finish of the 2nd-order axis flowering";
//        gnuparam.ylabel = "h, relative height";
//        item.info = "";
//        (this->*(item.pf1))(item, flowerdome, {sprout});
//        (this->*(item.pf2))(item, default_width);
//        dependencies_items.push_back(item);

        item.title = "flowering curves";
        item.correlation_coef = false;
        item.determination_coef = true;
        item.normality_test = true;
        item.pf1 = &Sproutmodel::flowering_by_time;
        item.pf2 = &Sproutmodel::flowering_by_time;
        item.pf3 = &Sproutmodel::flowering_by_time_test;
        gnuparam.clear();
        gnuparam.title = item.title;
        gnuparam.xlabel = "T, days from the beginning of flowering";
        gnuparam.ylabel = "n, open flowers";
        item.info = "";
        (this->*(item.pf1))(item, flowerdome, {sprout});
        (this->*(item.pf2))(item, default_width);
        dependencies_items.push_back(item);

        item.title = "flowering ranges";
        item.correlation_coef = false;
        item.determination_coef = false;
        item.normality_test = false;
        item.pf1 = &Sproutmodel::flowering_ranges;
        item.pf2 = &Sproutmodel::flowering_ranges;
        item.pf3 = &Sproutmodel::flowering_ranges_test;
        gnuparam.clear();
        gnuparam.title = item.title;
        gnuparam.xlabel = "T, days from the beginning of flowering";
        gnuparam.ylabel = "h, relative height";
        item.info = "";
        (this->*(item.pf1))(item, flowerdome, {sprout});
        (this->*(item.pf2))(item, default_width);
        dependencies_items.push_back(item);
    }

    inline shared<Sprout> Sproutmodel::generate_sprout(Unsigned flower1_number,
        Unsigned flower2_number, const std::string &flowering_start_date, Unsigned flowering_days)
    {
        Time sprout_flowering_start(flowering_start_date.c_str());
        std::uniform_int_distribution<> dayshift_dist(0, flowering_days);
        sprout_flowering_start = sprout_flowering_start + dayshift_dist(randgen());
            // случайная дата начала цветения побега

        shared<Sprout> sprout(new Sprout);

        //-----------------------------------------------------------------------------------------
        //                    генерация структуры побега
        //-----------------------------------------------------------------------------------------
        if (flower1_number == inf<Unsigned>()) flower1_number = generate_flower1();
        else if (flower1_number == 0) flower1_number = nflower1;

        if (flower2_number == inf<Unsigned>()) flower2_number = generate_flower2(flower1_number);
        else if (flower2_number == 0) flower2_number = nflower2;

        auto heights1 = generate_heights(flower1_number, 1); // создадим список высот 1-го порядка
        auto heights2 = generate_heights(flower1_number, 2); // создадим список высот 2-го порядка

        shared<curve<2>> density2 = generate_density2(heights1, flower2_number);

        // пересчитаем реальное число цветков 2-го порядка после генерации
        Unsigned sum = 0;
        for (auto d : *density2) sum += Unsigned(d[1]);
        //std::printf("number of 2nd order flowers reset : %ld -> %ld\n", flower2_number, sum);
        flower2_number = sum;

        flowerline line;

        line.fsp = FlowerIndex(1, false, std::array<unsigned, 4>{(unsigned)flower1_number, (unsigned)flower2_number, 0, 0});
        sprout->push_back(line); // структурная запись

        for (Unsigned k1=1; k1<=flower1_number; ++k1)
        {
            line.fsp = FlowerIndex(1, true, std::array<unsigned, 4>{(unsigned)k1, 0, 0, 0});
            sprout->push_back(line);
        }

        if (flower2_number != 0)
        {
            for (auto d : *density2)
            {
                Unsigned k1 = number_from_relative_height(d[0], flower1_number);
                Unsigned n = d[1];
                for (Unsigned k2=1; k2<=n; ++k2)
                {
                    line.fsp = FlowerIndex(1, true, std::array<unsigned, 4>{(unsigned)k1, (unsigned)k2, 0, 0});
                    sprout->push_back(line);
                }
            }
        }
        std::sort(sprout->begin(), sprout->end(), [](auto a, auto b){ return a.fsp < b.fsp; });

        //-----------------------------------------------------------------------------------------
        //                    генерация динамики цветения побега
        //-----------------------------------------------------------------------------------------
        shared<curve<2>> start1 = generate_start1(heights1);

        shared<curve<2>> start2;
        if (flower2_number != 0) start2 = generate_start2(heights2);

        for (auto e : *start1)
        {
            Unsigned k1 = number_from_relative_height(e[0], flower1_number);
            auto line = sprout->operator[](FlowerIndex(1, true, std::array<unsigned, 4>{(unsigned)k1, 0, 0, 0}));

            Time flower_start = e[1] + sprout_flowering_start.to_double();
            rough_value<Time> *range = line->start_interval(FLOWER);
            range->range() = { flower_start-1, flower_start};
            range->value() = flower_start;

            double dt = generate_duration(BUD);
            Time bud_start = flower_start - dt;
            range = line->start_interval(BUD);
            range->range() = { bud_start-1, bud_start };
            range->value() = bud_start ;

            dt = generate_duration(FLOWER);
            Time ovary_start = flower_start + dt;
            range = line->start_interval(OVARY);
            range->range() = { ovary_start-1, ovary_start };
            range->value() = ovary_start;

            dt = generate_duration(OVARY);
            Time seed_start = ovary_start + dt;
            range = line->start_interval(SEED);
            range->range() = { seed_start-1, seed_start };
            range->value() = seed_start;

            Time death = seed_start + 90;
            range = line->start_interval(DEATH);
            range->range() = { death-1, death };
            range->value() = death;
        }

        if (flower2_number != 0)
        {
            for (auto e : *start2)
            {
                double h = e[0];
                Unsigned k1 = number_from_relative_height(h, flower1_number);
                auto line1 = sprout->operator[](FlowerIndex(1, true, std::array<unsigned, 4>{(unsigned)k1, 0, 0, 0}));
                Time t1 = line1->value(FLOWER); // старт терминала боковой оси

                auto r = find_exact(density2, h, 0.001);
                assert(r.flag() == true);
                Unsigned n2 = Unsigned(r.value()[1]);
                    // число цветков на оси 2-го порядка

                Time flower_start = t1 + e[1]; // старт первого цветка боковой оси, отличного от терминала
                for (Unsigned k2=1; k2<=n2; ++k2)
                {
                    auto line = sprout->operator[](FlowerIndex(1, true, std::array<unsigned, 4>{(unsigned)k1, (unsigned)k2, 0, 0}));

                    rough_value<Time> *range = line->start_interval(FLOWER);
                    range->range() = { flower_start-1, flower_start };
                    range->value() = flower_start;

                    double dt = generate_duration(BUD);
                    Time bud_start = flower_start - dt;
                    range = line->start_interval(BUD);
                    range->range() = { bud_start-1, bud_start };
                    range->value() = bud_start;

                    dt = generate_duration(FLOWER);
                    Time ovary_start = flower_start + dt;
                    range = line->start_interval(OVARY);
                    range->range() = { ovary_start-1, ovary_start };
                    range->value() = ovary_start;

                    dt = generate_duration(OVARY);
                    Time seed_start = ovary_start + dt;
                    range = line->start_interval(SEED);
                    range->range() = { seed_start-1, seed_start };
                    range->value() = seed_start;

                    Time death = seed_start + 90;
                    range = line->start_interval(DEATH);
                    range->range() = { death-1, death };
                    range->value() = death;

                    flower_start += generate_between();
                }
            }
        }
        return sprout;
    }

#undef DISTRIBUTION
#undef REGRESSION

}
#endif // SPROUT_MODEL_H









//    /// построить имя функции для распределения объектов(of) по параметру(by)
//    #define MAKE_DISTRIBUTION(of, by) distribution_function_##of##_##by

//    /// построить имя переменной для распределения объектов(of) по параметру(by)
//    #define DISTRIBUTION(of, by) distribution_var_##of##_##by

//    /// построить имя функции для регрессии параметра(of) от параметра(on)
//    #define MAKE_REGRESSION(of, on) regression_function_##of##_##on

//    /// построить имя перемнной для функции для регрессии параметра(of) от параметра(on)
//    #define REGRESSION(of, on) regression_var_##of##_##on

//    /// построить имя функции для регрессии параметра(of) от параметра(on)
//    #define MAKE_CURVE(of, on) curve_function_##of##_##on

//    /// построить имя перемнной для функции для регрессии параметра(of) от параметра(on)
//    #define CURVE(of, on) curve_var_##of##_##on

//    const double default_width = 0.1;

//    /*!
//     * \brief сгладить распределение с заданной шириной окна
//     * \param data случайные данные
//     * \param width ширина окна
//     * \return сглаженные данные
//     */
//    inline shared<curve<2>> smooth_distribution(shared<curve<2>> data, double width=default_width)
//    {
//        WLSM_item options(default_WLSM_options);
//        options.width = {width, true};

//        auto f = smooth(1, data, options);
//        cut_below(f);
//        cut_left(f);
//        return f;
//    }
//    inline shared<curve<3>> distribution_smooth_with_errors(shared<curve<2>> data, double width=default_width)
//    {
//        WLSM_item options(default_WLSM_options);
//        options.width = {width, true};

//        auto f = WLSM().smooth_with_errors(1, data, options);
//        cut_below(f);
//        cut_left(f);
//        return f;
//    }

//    /*!
//     * \brief сгладить регрессию
//     * \note результат сглаживания сводится к диапазону опредленеия данных
//     * \param data случайные данные
//     * \param width ширина сглаживания
//     * \return сглаженная кривая со средней ошибкой сглаживания
//     */
//    inline shared<curve<3>> smooth_regression(shared<curve<2>> data, double width=0.2)
//    {
//        WLSM_item options(default_WLSM_options);
//        options.width = {width, true};
//        auto f = WLSM().smooth_with_errors(0, data, options);

//        Interval<double,2> xrange(data->front()[0], data->back()[0]);
//        to_xrange(data, xrange);
//        to_xrange(f, xrange);
//        return f;
//    }

//    /*!
//     * \brief модель побега
//     * \note модель состоит из совокупности распределений и регрессионных соотношений,
//     * позволяющих по ним генерить образцы модельных побегов
//     */
//    class Sproutmodel
//    {
//        //------------------------------------------------------------------------------------------
//        //              описание структуры побега (вплоть до осей 2-го порядка)
//        //------------------------------------------------------------------------------------------
//        mutable numeric_distribution DISTRIBUTION(shoots, flowers1); // распределение побегов по числу цветков 1-го порядка
//        mutable numeric_distribution DISTRIBUTION(flowers, budding_time); // распределение цветков по времени бутонизации
//        mutable numeric_distribution DISTRIBUTION(flowers, flowering_time); // распределение цветков по времени цветения
//        mutable numeric_distribution DISTRIBUTION(flowers, ovaring_time); // распределение цветков по времени завязывания
//        mutable numeric_distribution DISTRIBUTION(flowers, between_time); // распределение пар цветков по времени между последовательными распусканиями на одной оси

//        shared<curve<3>> REGRESSION(flowers2, flowers1); // зависимость числа цветков 2-го порядков от числа цветков 1-го порядка
//        shared<curve<3>> REGRESSION(density2, height); // плотность цветков 2-го порядка по высоте

//        //------------------------------------------------------------------------------------------
//        //                      описание динамики цветения побега
//        //------------------------------------------------------------------------------------------
//        shared<curve<3>> REGRESSION(start1, height); // кривая старта цветения оси 1-го порядка (относительно начала цветения побега)
//        shared<curve<3>> REGRESSION(start2, height); // кривая старта цветения оси 2-го порядка (относительно начала цветения терминалов)
//        shared<curve<3>> REGRESSION(finish2, height); // кривая конца цветения оси 2-го порядка (относительно начала цветения терминалов)

//        shared<curve<3>> CURVE(flowering1, time); // кривая цветения 1-го порядка
//        shared<curve<3>> CURVE(flowering2, time); // кривая цветения 2-го порядка
//        shared<curve<3>> CURVE(flowering, time); // кривая цветения 1-uj + 2-го порядков

//        typedef std::pair<Unsigned, gnuplotParams> view_type; // тип объекта, описывающего ту или иную зависимость и опции ее отрисовки
//        std::map<Unsigned, gnuplotParams> dependencies; // список всех сформированных зависимостей
//            // имя зависимости лежит в title of gnuplotParams

//        Interval<double, 2> get_interval(shared<curve<2>> f) const { return {f->front()[0], f->back()[0]}; }
//        Interval<double, 2> get_interval(shared<curve<3>> f) const { return {f->front()[0], f->back()[0]}; }

//    public:

//        typedef void (Sproutmodel::*PF)(const Flowerdome &, const std::set<Unsigned> &, double);
//        std::map<Unsigned, std::pair<gnuplotParams, PF>> Xdependencies; // список всех сформированных зависимостей
//            // имя зависимости лежит в title of gnuplotParams


//        void MAKE_DISTRIBUTION(shoots, flowers1)(const Flowerdome &flowerdome, Unsigned)
//        {
//            auto sprouts = flowerdome.get_sprouts();
//            double hx = 0.1;
//            double left_indent = 0.3;
//            double right_indent = 0.3;

//            shared<curve<2>> data = flowerdome.foreach(sprouts, flower_count_functor(1));
//            to_cumulative_density_function(data, hx, left_indent, right_indent);

//            auto f2 = smooth_distribution(data);
//            auto f3 = distribution_smooth_with_errors(data);

//            gnuplotParams gnuparam;
//            gnuparam.title = "distribution of shoots by the number of flowers of the 1st order";
//            gnuparam.xlabel = "number of the 1st-order flowers";
//            gnuparam.ylabel = "probability density";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(f2);
//            dependencies[1] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            DISTRIBUTION(shoots, flowers1).init(f2->begin(), f2->end(),
//                [](auto v) { return std::make_pair(v[0], v[1]); }
//            ); // инициализация распределения
//        }

//        void MAKE_DISTRIBUTION(flowers, budding_time)(const Flowerdome &flowerdome, Unsigned)
//        {
//            auto sprouts = flowerdome.get_sprouts();
//            double hx = 0.1;
//            double left_indent = 0.3;
//            double right_indent = 0.3;

//            shared<curve<2>> data = flowerdome.foreach(sprouts, phase_duration_functor(BUD));
//            to_cumulative_density_function(data, hx, left_indent, right_indent);

//            auto f2 = smooth_distribution(data);

//            gnuplotParams gnuparam;
//            gnuparam.title = "distribution of flowers by duration of budding";
//            gnuparam.xlabel = "budding duration";
//            gnuparam.ylabel = "probability density";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(f2);
//            dependencies[2] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            DISTRIBUTION(flowers, budding_time).init(f2->begin(), f2->end(),
//                [](auto v) { return std::make_pair(v[0], v[1]); }
//            ); // инициализация распределения
//        }

//        void MAKE_DISTRIBUTION(flowers, flowering_time)(const Flowerdome &flowerdome, Unsigned)
//        {
//            auto sprouts = flowerdome.get_sprouts();
//            double hx = 0.1;
//            double left_indent = 0.3;
//            double right_indent = 0.3;

//            shared<curve<2>> data = flowerdome.foreach(sprouts, phase_duration_functor(FLOWER));
//            to_cumulative_density_function(data, hx, left_indent, right_indent);

//            auto f2 = smooth_distribution(data);

//            gnuplotParams gnuparam;
//            gnuparam.title = "distribution of flowers by duration of flowering";
//            gnuparam.xlabel = "flowering duration";
//            gnuparam.ylabel = "probability density";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(f2);
//            dependencies[3] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            DISTRIBUTION(flowers, flowering_time).init(f2->begin(), f2->end(),
//                [](auto v) { return std::make_pair(v[0], v[1]); }
//            ); // инициализация распределения
//        }

//        void MAKE_DISTRIBUTION(flowers, ovaring_time)(const Flowerdome &flowerdome, Unsigned)
//        {
//            auto sprouts = flowerdome.get_sprouts();
//            double hx = 0.1;
//            double left_indent = 0.3;
//            double right_indent = 0.3;

//            shared<curve<2>> data = flowerdome.foreach(sprouts, phase_duration_functor(OVARY));
//            to_cumulative_density_function(data, hx, left_indent, right_indent);

//            auto f2 = smooth_distribution(data);

//            gnuplotParams gnuparam;
//            gnuparam.title = "distribution of flowers by duration of ovaring";
//            gnuparam.xlabel = "ovaring duration";
//            gnuparam.ylabel = "probability density";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(f2);
//            dependencies[4] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            DISTRIBUTION(flowers, ovaring_time).init(f2->begin(), f2->end(),
//                [](auto v) { return std::make_pair(v[0], v[1]); }
//            ); // инициализация распределения
//        }

//        void MAKE_DISTRIBUTION(flowers, between_time)(const Flowerdome &flowerdome, Unsigned)
//        {
//            auto sprouts = flowerdome.get_sprouts();
//            double hx = 0.1;
//            double left_indent = 0.3;
//            double right_indent = 0.3;

//            shared<curve<2>> data(new curve<2>);
//            for (auto sprout : sprouts)
//            {
//                auto tmp = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_START));
//                if (tmp->empty()) continue;

//                double search_eps = 0.01;
//                std::sort(tmp->begin(), tmp->end(), [search_eps](auto a, auto b)
//                    {
//                        if (equal(a[0], b[0], search_eps)) return a[1] < b[1];
//                        return a[0]< b[0];
//                    }
//                ); // сортировка для извлечения разниц

//                data += to_beetween(tmp, search_eps);
//            }
//            to_cumulative_density_function(data, hx, left_indent, right_indent);

//            auto f2 = smooth_distribution(data);

//            gnuplotParams gnuparam;
//            gnuparam.title = "distribution of flowers on axes of the 2nd order by time intervals between their openings";
//            gnuparam.xlabel = "duration between flower openings";
//            gnuparam.ylabel = "probability density";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(f2);
//            dependencies[5] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            DISTRIBUTION(flowers, between_time).init(f2->begin(), f2->end(),
//                [](auto v) { return std::make_pair(v[0], v[1]); }
//            ); // инициализация распределения
//        }

//        /// dependence of the number of flowers of the 2nd order on the number of flowers of the 2nd order
//        void MAKE_REGRESSION(flowers2, flowers1)(const Flowerdome &flowerdome, Unsigned)
//        {
//            auto sprouts = flowerdome.get_sprouts();
//            auto data = flowerdome.foreach(sprouts, flower_rate_functor(1, 2));

//            auto f3 = smooth_regression(data, 0.2);

//            gnuplotParams gnuparam;
//            gnuparam.title = "regression between the number of 1st-order flowers and the number of 2nd-order flowers";
//            gnuparam.xlabel = "number of the 1st-order flowers";
//            gnuparam.ylabel = "number of the 2nd-order flowers";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(f3), "u 1:2 w filledcu closed");
//            gnuparam.insert_data(extract_mean(f3), "u 1:2 w l");
//            gnuparam.insert_data(data, "u 1:2 w p");
//            dependencies[6] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            REGRESSION(flowers2, flowers1) = f3;
//        }

//        void MAKE_REGRESSION(density2, height)(const Flowerdome &flowerdome, Unsigned sprout)
//        {
//            double hx = 0.01;
//            const Sprout &s = flowerdome.at(sprout);
//            Unsigned nf1 = s.number_of_axis1_flowers();

//            auto data = flowerdome.foreach({sprout}, axis_density_functor(2));
//            unite(data, hx);
//            to_height(data, nf1);

//            auto f3 = smooth_regression(data, 0.2);

//            gnuplotParams gnuparam;
//            gnuparam.title = "the number of flowers on the axes of the 2nd order, depending on the relative height of the axis";
//            gnuparam.xlabel = "number of the 2nd-order flowers";
//            gnuparam.ylabel = "relative high";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(f3), "u 2:1 w filledcu closed");
//            gnuparam.insert_data(extract_mean(f3), "u 2:1 w l");
//            gnuparam.insert_data(data, "u 2:1 w p");
//            dependencies[7] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            REGRESSION(density2, height) = f3;
//        }

//        void MAKE_REGRESSION(start1, height)(const Flowerdome &flowerdome, Unsigned sprout)
//        {
//            const Sprout &s = flowerdome.at(sprout);
//            Unsigned nf1 = s.number_of_axis1_flowers();

//            auto data = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
//            to_height(data, nf1);

//            auto f3 = smooth_regression(data, 0.2);

//            gnuplotParams gnuparam;
//            gnuparam.title = "the order of blooming flowers on the axis of the 1st order";
//            gnuparam.xlabel = "start of the 1st-order axis flowering";
//            gnuparam.ylabel = "relative high";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(f3), "u 2:1 w filledcu closed");
//            gnuparam.insert_data(extract_mean(f3), "u 2:1 w l");
//            gnuparam.insert_data(data, "u 2:1 w p");
//            dependencies[8] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            REGRESSION(start1, height) = f3;
//        }

//        void MAKE_REGRESSION(start2, height)(const Flowerdome &flowerdome, Unsigned sprout)
//        {
//            const Sprout &s = flowerdome.at(sprout);
//            Unsigned nf1 = s.number_of_axis1_flowers();

//            auto data1 = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
//            to_height(data1, nf1); // для получения относительног старта

//            auto data2 = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_START));

//            double search_eps = 0.01;
//            double represent_part = 0.3;
//            extract_bound(data2, true, search_eps, represent_part);

//            to_height(data2, nf1);
//            to_relative(data2, data1);

//            auto f3 = smooth_regression(data2, 0.2);

//            gnuplotParams gnuparam;
//            gnuparam.title = "The delays in the beginning of flowering of axes of the 2nd order";
//            gnuparam.xlabel = "relative start of the 2nd-order axis flowering";
//            gnuparam.ylabel = "relative high";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(f3), "u 2:1 w filledcu closed");
//            gnuparam.insert_data(extract_mean(f3), "u 2:1 w l");
//            gnuparam.insert_data(data2, "u 2:1 w p");
//            dependencies[9] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            REGRESSION(start2, height) = f3;
//        }

//        void MAKE_REGRESSION(finish2, height)(const Flowerdome &flowerdome, Unsigned sprout)
//        {
//            const Sprout &s = flowerdome.at(sprout);
//            Unsigned nf1 = s.number_of_axis1_flowers();

//            auto data1 = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
//            to_height(data1, nf1); // для получения относительног старта

//            auto data2 = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_END));

//            double search_eps = 0.01;
//            double represent_part = 0.3;
//            extract_bound(data2, false, search_eps, represent_part);

//            to_height(data2, nf1);
//            to_relative(data2, data1);

//            auto f3 = smooth_regression(data2, 0.2);

//            gnuplotParams gnuparam;
//            gnuparam.title = "The delays in the finishing of flowering of axes of the 2nd order";
//            gnuparam.xlabel = "relative finish of the 2nd-order axis flowering";
//            gnuparam.ylabel = "relative high";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(f3), "u 2:1 w filledcu closed");
//            gnuparam.insert_data(extract_mean(f3), "u 2:1 w l");
//            gnuparam.insert_data(data2, "u 2:1 w p");
//            dependencies[10] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            REGRESSION(finish2, height) = f3;
//        }

//        void MAKE_CURVE(flowering1, time)(const Flowerdome &flowerdome, Unsigned sprout, double width=2.)
//        {
//            auto data = flowerdome.foreach({sprout}, extraction_phase_of_axis(1, FLOWER));
//            auto min = find_minimal(data, 1); // получить дату самого раннего старта
//            auto max = find_maximal(data, 2); // получить дату самого позднего финиша

////            typename prg::vector<double>::range_type range(Signed(min), Signed(max+1));
////            prg::vector<double> sums(range, 0.);
//            prg::vector<double> sums(Signed(min), Signed(max - min + 1), 0.);

//            for (auto e : *data)
//            {
//                Signed t1 = e[1];
//                Signed t2 = e[2];
//                for (auto t=t1; t<=t2; t++) sums[t] += 1;
//            }

//            shared<curve<2>> tmp(new curve<2>);
//            for (Signed t=sums.lower_index(); t<sums.upper_index(); t++)
//                tmp->push_back({double(t - min), sums[t]}); // приведем к относительной дате

//            Unsigned n = sums.size();
//            double hx = 1. / n;
//            double left_indent = 5. / n;
//            double right_indent = 5. / n;
//            width /= n;
//                // относительные значения

//            fill_gaps(tmp, hx, left_indent, right_indent);

//            auto f3 = smooth_regression(tmp, width);

//            gnuplotParams gnuparam;
//            gnuparam.title = "flowering curves [1st]";
//            gnuparam.xlabel = "T, days from the beginning of flowering";
//            gnuparam.ylabel = "%, open flowers";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(f3), "u 1:2 w filledcu closed fc 'red' fill transparent solid 0.2 noborder ");
//            gnuparam.insert_data(extract_mean(f3), "u 1:2 w l lc 'red' ");
//            gnuparam.insert_data(tmp, "u 1:2 w p pt 7 ");

//            dependencies[11] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            CURVE(flowering1, time) = f3;
//        }

//        void MAKE_CURVE(flowering2, time)(const Flowerdome &flowerdome, Unsigned sprout)
//        {
//            auto data = flowerdome.foreach({sprout}, extraction_phase_of_axis(2, FLOWER));
//            auto min = find_minimal(data, 1); // получить дату самого раннего старта
//            auto max = find_maximal(data, 2); // получить дату самого позднего финиша

////            typename prg::vector<double>::range_type range(Signed(min), Signed(max+1));
////            prg::vector<double> sums(range, 0.);
//            prg::vector<double> sums(Signed(min), Signed(max - min + 1), 0.);

//            for (auto e : *data)
//            {
//                Signed t1 = e[1];
//                Signed t2 = e[2];
//                for (auto t=t1; t<=t2; t++) sums[t] += 1;
//            }

//            shared<curve<2>> tmp(new curve<2>);
//            for (Signed t=sums.lower_index(); t<sums.upper_index(); t++)
//                tmp->push_back({double(t - min), sums[t]}); // приведем к относительной дате

//            Unsigned n = sums.size();
//            double hx = 1. / n;
//            double left_indent = 5. / n;
//            double right_indent = 5. / n;
//            double width = 2. / n;
//                // относительные значения

//            fill_gaps(tmp, hx, left_indent, right_indent);

//            auto f3 = smooth_regression(tmp, width);

//            gnuplotParams gnuparam;
//            gnuparam.title = "flowering curves [2nd]";
//            gnuparam.xlabel = "T, days from the beginning of flowering";
//            gnuparam.ylabel = "%, open flowers";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

////            "u 2:1 w filledcu closed fc 'red' fill transparent solid 0.2 noborder title '{/Symbol d}f_{t_1}(h)'",
////            "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder title '{/Symbol d}f_{t_2}(h)'",
////            "u 2:1 w l lc 'red' title 'f_{t_1}(h)' ",
////            "u 2:1 w l lc 'blue' title 'f_{t_2}(h)' ",
////            "u 2:1:4 w xerr pt 7 fc 'gray' title '' ",
////            "u 2:1:4 w xerr pt 7 fc 'red' title '' ",
////            "u 2:1:4 w xerr pt 7 fc 'blue' title '' "

//            gnuparam.insert_data(extract_errors_range(f3), "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ");
//            gnuparam.insert_data(extract_mean(f3), "u 1:2 w l lc 'blue' ");
//            gnuparam.insert_data(tmp, "u 1:2 w p pt 7 ");
//            dependencies[12] = gnuparam;

//        #ifdef TEST_PLOT
//            plot(gnuparam);
//        #endif

//            CURVE(flowering2, time) = f3;
//        }

//        //-----------------------------------------------------------------------------------------
//        //-----------------------------------------------------------------------------------------
//        //-----------------------------------------------------------------------------------------
//        struct dependence_item
//        {
//            typedef void (Sproutmodel::*PF1)(dependence_item &, const Flowerdome &, const std::set<Unsigned> &);
//            typedef void (Sproutmodel::*PF2)(dependence_item &, double);

//            std::string title =""; // name of dependence
//            PF1 pf1 = nullptr; // data preparation function
//            PF2 pf2 = nullptr; // data smoothing function

//            shared<curve<2>> data1 = nullptr; // experimental data for 1st-order
//            shared<curve<2>> data2 = nullptr; // experimental data for 2nd-order
//            shared<curve<2>> data = nullptr; // sum of experimental data 1st + 2nd
//            shared<curve<2>> ddata = nullptr; // diff of experimental data 1st - 2nd

//            shared<curve<3>> smooth1 = nullptr;  // smoothed data with errors
//            shared<curve<3>> smooth2 = nullptr;  // smoothed data with errors
//            shared<curve<3>> smooth = nullptr;  // smoothed data with errors
//            shared<curve<3>> dsmooth = nullptr;  // smoothed data with errors

//            gnuplotParams gnuparam; // data plot parameters
//        };
//        std::vector<dependence_item> dependencies_items;

//        //-----------------------------------------------------------------------------------------
//        void shoots_by_flowers1(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            double hx = 0.1;
//            double left_indent = 0.3;
//            double right_indent = 0.3;

//            auto data = flowerdome.foreach(sprouts, flower_count_functor(1));
//            if (data->empty()) return;

//            to_cumulative_density_function(data, hx, left_indent, right_indent);
//            item.data = data;
//        }

//        void shoots_by_flowers1(dependence_item &item, double width)
//        {
//            item.smooth = distribution_smooth_with_errors(item.data, width);
//            auto f2 = extract(item.smooth, 0, 1);
//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "distribution of shoots by the number of flowers of the 1st order";
//            gnuparam.xlabel = "number of the 1st-order flowers";
//            gnuparam.ylabel = "probability density";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";
//            gnuparam.insert_data(f2);

//            DISTRIBUTION(shoots, flowers1).init(f2->begin(), f2->end(),
//                [](auto v) { return std::make_pair(v[0], v[1]); }
//            ); // инициализация распределения
//        }

//        //-----------------------------------------------------------------------------------------
//        shared<curve<2>> flowers_by_phase_duration_time(Phase phase,
//            const Flowerdome &flowerdome, const std::set<Unsigned> &sprouts)
//        {
//            double hx = 0.1;
//            double left_indent = 0.3;
//            double right_indent = 0.3;

//            shared<curve<2>> data = flowerdome.foreach(sprouts, phase_duration_functor(phase));

//            if (!data->empty())
//                to_cumulative_density_function(data, hx, left_indent, right_indent);

//            return data;
//        }

//        void flowers_by_budding_time(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            item.data = flowers_by_phase_duration_time(BUD, flowerdome, sprouts);
//        }
//        void flowers_by_budding_time(dependence_item &item, double width)
//        {
//            item.smooth = distribution_smooth_with_errors(item.data, width);
//            auto f2 = extract(item.smooth, 0, 1);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "distribution of flowers by duration of budding";
//            gnuparam.xlabel = "budding duration";
//            gnuparam.ylabel = "probability density";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";
//            gnuparam.insert_data(f2);

//            DISTRIBUTION(flowers, budding_time).init(f2->begin(), f2->end(),
//                [](auto v) { return std::make_pair(v[0], v[1]); }
//            ); // инициализация распределения
//        }

//        //-----------------------------------------------------------------------------------------
//        void flowers_by_flowering_time(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            item.data = flowers_by_phase_duration_time(FLOWER, flowerdome, sprouts);
//        }

//        void flowers_by_flowering_time(dependence_item &item, double width)
//        {
//            item.smooth = distribution_smooth_with_errors(item.data, width);
//            auto f2 = extract(item.smooth, 0, 1);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "distribution of flowers by duration of flowering";
//            gnuparam.xlabel = "flowering duration";
//            gnuparam.ylabel = "probability density";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";
//            gnuparam.insert_data(f2);

//            DISTRIBUTION(flowers, flowering_time).init(f2->begin(), f2->end(),
//                [](auto v) { return std::make_pair(v[0], v[1]); }
//            ); // инициализация распределения
//        }

//        //-----------------------------------------------------------------------------------------
//        void flowers_by_ovaring_time(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            item.data = flowers_by_phase_duration_time(OVARY, flowerdome, sprouts);
//        }

//        void flowers_by_ovaring_time(dependence_item &item, double width)
//        {
//            item.smooth = distribution_smooth_with_errors(item.data, width);
//            auto f2 = extract(item.smooth, 0, 1);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "distribution of flowers by duration of ovaring";
//            gnuparam.xlabel = "ovaring duration";
//            gnuparam.ylabel = "probability density";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";
//            gnuparam.insert_data(f2);

//            DISTRIBUTION(flowers, ovaring_time).init(f2->begin(), f2->end(),
//                [](auto v) { return std::make_pair(v[0], v[1]); }
//            ); // инициализация распределения
//        }

//        //-----------------------------------------------------------------------------------------
//        void flowers_by_between_time(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            shared<curve<2>> data(new curve<2>);
//            for (auto sprout : sprouts)
//            {
//                auto tmp = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_START));
//                if (tmp->empty()) continue;

//                double search_eps = 0.01;
//                std::sort(tmp->begin(), tmp->end(), [search_eps](auto a, auto b)
//                    {
//                        if (equal(a[0], b[0], search_eps)) return a[1] < b[1];
//                        return a[0]< b[0];
//                    }
//                ); // сортировка для извлечения разниц

//                data += to_beetween(tmp, search_eps);
//            }
//            item.data = data;

//            double hx = 0.1;
//            double left_indent = 0.3;
//            double right_indent = 0.3;

//            to_cumulative_density_function(item.data, hx, left_indent, right_indent);
//        }

//        void flowers_by_between_time(dependence_item &item, double width)
//        {
//            item.smooth = distribution_smooth_with_errors(item.data, width);
//            auto f2 = extract(item.smooth, 0, 1);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "distribution of flowers on axes of the 2nd order by time intervals between their openings";
//            gnuparam.xlabel = "duration between flower openings";
//            gnuparam.ylabel = "probability density";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";
//            gnuparam.insert_data(f2);

//            DISTRIBUTION(flowers, between_time).init(f2->begin(), f2->end(),
//                [](auto v) { return std::make_pair(v[0], v[1]); }
//            ); // инициализация распределения
//        }

//        //-----------------------------------------------------------------------------------------
//        void flowers2_by_flowers1(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            item.data =  flowerdome.foreach(sprouts, flower_rate_functor(1, 2));
//        }

//        void flowers2_by_flowers1(dependence_item &item, double width)
//        {
//            item.smooth = smooth_regression(item.data, width);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "regression between the number of 1st-order flowers and the number of 2nd-order flowers";
//            gnuparam.xlabel = "number of the 1st-order flowers";
//            gnuparam.ylabel = "number of the 2nd-order flowers";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(item.smooth), "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.2 ");
//            gnuparam.insert_data(extract_mean(item.smooth), "u 1:2 w l lc 'blue' ");
//            gnuparam.insert_data(item.data, "u 1:2 w p pt 7 lc rgb 'blue' ");

//            REGRESSION(flowers2, flowers1) = item.smooth;
//        }
//        //-----------------------------------------------------------------------------------------

//        void density2_by_height(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            double hx = 0.01;

//            shared<curve<2>> data(new curve<2>);
//            for (auto sprout : sprouts)
//            {
//                const Sprout &s = flowerdome.at(sprout);
//                Unsigned nf1 = s.number_of_axis1_flowers();

//                auto tmp = flowerdome.foreach({sprout}, axis_density_functor(2));
//                if (tmp->empty()) continue;

//                unite(tmp, hx);
//                to_height(tmp, nf1);

//                for (auto e : *tmp) data->push_back(e);
//            }
//            item.data = data;
//        }

//        void density2_by_height(dependence_item &item, double width)
//        {
//            item.smooth = smooth_regression(item.data, width);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "the number of flowers on the axes of the 2nd order, depending on the relative height of the axis";
//            gnuparam.xlabel = "number of the 2nd-order flowers";
//            gnuparam.ylabel = "relative high";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(item.smooth), "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ");
//            gnuparam.insert_data(extract_mean(item.smooth), "u 2:1 w l lc 'blue' ");
//            gnuparam.insert_data(item.data, "u 2:1 w p pt 7 lc rgb 'blue' ");

//            REGRESSION(density2, height) = item.smooth;
//        }
//        //-----------------------------------------------------------------------------------------

//        void start1_by_height(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            shared<curve<2>> data(new curve<2>);
//            for (auto sprout : sprouts)
//            {
//                const Sprout &s = flowerdome.at(sprout);
//                Unsigned nf1 = s.number_of_axis1_flowers();

//                auto tmp = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
//                if (tmp->empty()) continue;

//                to_height(tmp, nf1);

//                for (auto e : *tmp) data->push_back(e);
//            }
//            item.data = data;
//        }

//        void start1_by_height(dependence_item &item, double width)
//        {
//            item.smooth = smooth_regression(item.data, width);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "the order of blooming flowers on the axis of the 1st order";
//            gnuparam.xlabel = "start of the 1st-order axis flowering";
//            gnuparam.ylabel = "relative high";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(item.smooth), "u 2:1 w filledcu closed fc 'red' fill transparent solid 0.2 noborder ");
//            gnuparam.insert_data(extract_mean(item.smooth), "u 2:1 w l lc 'red' ");
//            gnuparam.insert_data(item.data, "u 2:1 w p pt 7 lc rgb 'red' ");

//            REGRESSION(start1, height) = item.smooth;
//        }

//        //-----------------------------------------------------------------------------------------
//        void start2_by_height(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            shared<curve<2>> data(new curve<2>);
//            for (auto sprout : sprouts)
//            {
//                const Sprout &s = flowerdome.at(sprout);
//                Unsigned nf1 = s.number_of_axis1_flowers();

//                auto tmp1 = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
//                if (tmp1->empty()) continue;

//                to_height(tmp1, nf1); // для получения относительного старта

//                auto tmp2 = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_START));
//                if (tmp2->empty()) continue;

//                double search_eps = 0.01;
//                double represent_part = 0.3;
//                extract_bound(tmp2, true, search_eps, represent_part);

//                to_height(tmp2, nf1);
//                to_relative(tmp2, tmp1);

//                for (auto e : *tmp2) data->push_back(e);
//            }
//            item.data = data;
//        }

//        void start2_by_height(dependence_item &item, double width)
//        {
//            item.smooth = smooth_regression(item.data, width);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "The delays in the beginning of flowering of axes of the 2nd order";
//            gnuparam.xlabel = "relative start of the 2nd-order axis flowering";
//            gnuparam.ylabel = "relative high";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(item.smooth), "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ");
//            gnuparam.insert_data(extract_mean(item.smooth), "u 2:1 w l lc 'blue' ");
//            gnuparam.insert_data(item.data, "u 2:1 w p pt 7 lc rgb 'blue' ");

//            REGRESSION(start2, height) = item.smooth;
//        }

//        //-----------------------------------------------------------------------------------------
//        void finish2_by_height(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            shared<curve<2>> data(new curve<2>);
//            for (auto sprout : sprouts)
//            {
//                const Sprout &s = flowerdome.at(sprout);
//                Unsigned nf1 = s.number_of_axis1_flowers();

//                auto tmp1 = flowerdome.foreach({sprout}, axis_phase_functor(1, FLOWER, PHASE_START));
//                if (tmp1->empty()) continue;

//                to_height(tmp1, nf1); // для получения относительног старта

//                auto tmp2 = flowerdome.foreach({sprout}, axis_phase_functor(2, FLOWER, PHASE_END));
//                if (tmp2->empty()) continue;

//                double search_eps = 0.01;
//                double represent_part = 0.3;
//                extract_bound(tmp2, false, search_eps, represent_part);

//                to_height(tmp2, nf1);
//                to_relative(tmp2, tmp1);

//                for (auto e : *tmp2) data->push_back(e);
//            }
//            item.data = data;
//        }

//        void finish2_by_height(dependence_item &item, double width)
//        {
//            item.smooth = smooth_regression(item.data, width);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "The delays in the finishing of flowering of axes of the 2nd order";
//            gnuparam.xlabel = "relative finish of the 2nd-order axis flowering";
//            gnuparam.ylabel = "relative high";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(extract_errors_range(item.smooth), "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ");
//            gnuparam.insert_data(extract_mean(item.smooth), "u 2:1 w l lc 'blue' ");
//            gnuparam.insert_data(item.data, "u 2:1 w p pt 7 lc rgb 'blue' ");

//            REGRESSION(finish2, height) = item.smooth;
//        }

//        //-----------------------------------------------------------------------------------------
//        shared<curve<2>> development_by_time(extraction_phase_of_axis fn, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            Time tmin = 365.; Time tmax = 0.;
//            for (auto sprout : sprouts)
//            {
//                auto data = flowerdome.foreach({sprout}, fn);
//                    // data = { {height, start, finish}, {..} ... }

//                if (data->empty()) continue;

//                auto min = find_minimal(data, 1); // получить дату самого раннего старта
//                if (tmin > min.value()) tmin = min.value();

//                auto max = find_maximal(data, 2); // получить дату самого позднего финиша
//                if (tmax < max.value()) tmax = max.value();
//            }

//            Unsigned left_indent = 5;
//            Unsigned right_indent = 5;

//            Signed min_ndx = tmin.to_signed() - left_indent;
//            Signed max_ndx = tmax.to_signed() + right_indent;

//            Unsigned count = max_ndx - min_ndx + 1;
//            prg::vector<double> sums(min_ndx, count, 0.);

//            for (auto sprout : sprouts)
//            {
//                auto data = flowerdome.foreach({sprout}, fn);
//                if (data->empty()) continue;

//                for (auto e : *data)
//                {
//                    Signed t1 = e[1];
//                    Signed t2 = e[2];
//                    for (auto t=t1; t<=t2; t++) sums[t] += 1;
//                }
//            }

//            shared<curve<2>> data(new curve<2>);
//            for (Signed t=sums.minimal_index(); t<=sums.maximal_index(); t++)
//                data->push_back({double(t), sums[t]});

//            return data;
//        }

//        void flowering1_by_time(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            Time tstart = 365.; // starting date
//            for (auto sprout : sprouts)
//            {
//                auto data = flowerdome.foreach({sprout}, extraction_phase_of_axis(1, FLOWER, 0.));
//                if (data->empty()) continue;

//                auto min = find_minimal(data, 1); // получить дату самого раннего старта
//                if (tstart > min.value()) tstart = min.value();
//            }

//            item.data = development_by_time(extraction_phase_of_axis(1, FLOWER, tstart), flowerdome, sprouts);
//        }

//        void flowering1_by_time(dependence_item &item, double width)
//        {
//            item.smooth = smooth_regression(item.data, width);

//            auto mean = extract_mean(item.smooth);
//            cut_below(mean);
//            auto range = extract_errors_range(item.smooth);
//            cut_below(range);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "flowering curves [1st]";
//            gnuparam.xlabel = "T, days from the beginning of flowering";
//            gnuparam.ylabel = "%, open flowers";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            //            "u 2:1 w filledcu closed fc 'red' fill transparent solid 0.2 noborder title '{/Symbol d}f_{t_1}(h)'",
//            //            "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder title '{/Symbol d}f_{t_2}(h)'",
//            //            "u 2:1 w l lc 'red' title 'f_{t_1}(h)' ",
//            //            "u 2:1 w l lc 'blue' title 'f_{t_2}(h)' ",
//            //            "u 2:1:4 w xerr pt 7 fc 'gray' title '' ",
//            //            "u 2:1:4 w xerr pt 7 fc 'red' title '' ",
//            //            "u 2:1:4 w xerr pt 7 fc 'blue' title '' "

//            gnuparam.insert_data(range, "u 1:2 w filledcu closed fc 'red' fill transparent solid 0.2 noborder ");
//            gnuparam.insert_data(mean, "u 1:2 w l lc 'red' ");
//            gnuparam.insert_data(item.data, "u 1:2 w p pt 7 lc rgb 'red' ");

//            CURVE(flowering1, time) = item.smooth;
//        }

//        //-----------------------------------------------------------------------------------------
//        void flowering2_by_time(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            Time tstart = 365.; // date of beginning of flowering for the sprout
//            for (auto sprout : sprouts)
//            {
//                auto data = flowerdome.foreach({sprout}, extraction_phase_of_axis(1, FLOWER, 0.));
//                if (data->empty()) continue;

//                auto min = find_minimal(data, 1); // получить дату самого раннего старта
//                if (tstart > min.value()) tstart = min.value();
//            }

//            item.data = development_by_time(extraction_phase_of_axis(2, FLOWER, tstart), flowerdome, sprouts);
//        }

//        void flowering2_by_time(dependence_item &item, double width)
//        {
//            item.smooth = smooth_regression(item.data, width);

//            auto mean = extract_mean(item.smooth);
//            cut_below(mean);
//            auto range = extract_errors_range(item.smooth);
//            cut_below(range);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "flowering curves [2nd]";
//            gnuparam.xlabel = "T, days from the beginning of flowering";
//            gnuparam.ylabel = "%, open flowers";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(range, "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ");
//            gnuparam.insert_data(mean, "u 1:2 w l lc 'blue' ");
//            gnuparam.insert_data(item.data, "u 1:2 w p pt 7 lc rgb 'blue' ");

//            CURVE(flowering2, time) = item.smooth;
//        }

//        void flowering_by_time(dependence_item &item, const Flowerdome &flowerdome,
//            const std::set<Unsigned> &sprouts)
//        {
//            flowering1_by_time(item, flowerdome, sprouts);
//            item.data1 = item.data;

//            flowering2_by_time(item, flowerdome, sprouts);
//            item.data2 = item.data;

//            item.data = summarize<2>(item.data1, item.data2);
//            item.ddata = subtract<2>(item.data1, item.data2);
//        }

//        void flowering_by_time(dependence_item &item, double width)
//        {
//            item.smooth1 = smooth_regression(item.data1, width);
//            item.smooth2 = smooth_regression(item.data2, width);
//            item.smooth  = smooth_regression(item.data , width);
//            item.dsmooth = smooth_regression(item.ddata, width);

//            auto mean1 = extract_mean(item.smooth1);
//            cut_below(mean1);
//            auto mean2 = extract_mean(item.smooth2);
//            cut_below(mean2);
//            auto mean = extract_mean(item.smooth);
//            cut_below(mean);
//            auto dmean = extract_mean(item.dsmooth);

//            auto range1 = extract_errors_range(item.smooth1);
//            cut_below(range1);
//            auto range2 = extract_errors_range(item.smooth2);
//            cut_below(range2);

//            gnuplotParams &gnuparam = item.gnuparam;
//            gnuparam.title = "flowering curves [1st + 2nd]";
//            gnuparam.xlabel = "T, days from the beginning of flowering";
//            gnuparam.ylabel = "%, open flowers";
//            gnuparam.output = TMP_DIR + "/terminal.pdf";

//            gnuparam.insert_data(range1, "u 1:2 w filledcu closed fc 'red' fill transparent solid 0.2 noborder ");
//            gnuparam.insert_data(mean1, "u 1:2 w l lc 'red' ");
//            gnuparam.insert_data(range2, "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ");
//            gnuparam.insert_data(mean2, "u 1:2 w l lc 'blue' ");

//            gnuparam.insert_data(item.data1, "u 1:2 w p pt 7 lc rgb 'red' ");
//            gnuparam.insert_data(item.data2, "u 1:2 w p pt 7 lc rgb 'blue' ");
//            gnuparam.insert_data(item.data, "u 1:2 w p pt 7 lc rgb 'gray' ");
//            gnuparam.insert_data(mean, "u 1:2 w l lw 2 lc 'black' ");
//            gnuparam.insert_data(dmean, "u 1:2 w l lt 0 lw 2 lc 'black' ");

//            CURVE(flowering, time) = item.smooth;
//        }

//        //-----------------------------------------------------------------------------------------

//        /// получить случайное значение в точке x около заданной кривой
//        double generate_at(shared<curve<3>> f, double x) const
//        {
//            typedef std::uniform_real_distribution<> distribution_type;
//                // нормальное распределение дает большие выбросы за пределы области

//            distribution_type dist;
//            auto p = find_nearest(f, x);

//            auto flag = p.flag();
//            if (flag == false) return 0;

//            auto v = p.value();

//            distribution_type::param_type param(v[1] - v[2], v[1] + v[2]);
//            dist.param(param);

//            return dist(rangen);
//        }

//    public:

//        /// генерить массив событий, содержащий число цветков 1-го порядка
//        double generate_flower1() const
//        {
//            double v = DISTRIBUTION(shoots, flowers1)(randgen());
//            return v < 1 ? 1 : v;
//        }

//        /// генерить массив событий, содержащий число цветков 1-го порядка
//        double generate_between() const
//        {
//            double len = std::round(DISTRIBUTION(flowers, between_time)(randgen()));
//            return len < 0 ? 0 : len;
//        }

//        /// генерить массив событий, содержащий число цветков 1-го порядка
//        double generate_duration(Phase phase) const
//        {
//            double len = 0;
//            if (phase == BUD) len = DISTRIBUTION(flowers, budding_time)(randgen());
//            if (phase == FLOWER) len = DISTRIBUTION(flowers, flowering_time)(randgen());
//            if (phase == OVARY) len = DISTRIBUTION(flowers, ovaring_time)(randgen());
//            len = std::round(len);
//            return len < 1 ? 1 : len;
//        }

//        /// сгенерить массив числа цветков 1-го порядка (для тестирования)
//        shared<curve<2>> generate_flower1(Unsigned n) const
//        {
//            shared<curve<2>> out(new curve<2>);
//            for (Unsigned i=0; i<n; i++)
//            {
//                double v = generate_flower1();
//                out->push_back({std::ceil(v), 1});
//            }
//            return out;
//        }

//        /// генерить массив, содержащий число цветков 2-го порядка для соответствующего массива числа цветков 1-го порядка
//        double generate_flower2(double n1) const
//        {
//            double v = generate_at(REGRESSION(flowers2, flowers1), n1);
//            return v < 0 ? 0 : v;
//        }

//        /// генерить массив числа цветков 2-го порядка (для тестирования)
//        shared<curve<2>> generate_flower2(shared<curve<2>> n1) const
//        {
//            shared<curve<2>> out(new curve<2>);
//            for (auto e : *n1)
//            {
//                double nx = e[0];
//                double v = generate_flower2(nx);
//                out->push_back({nx, v});
//            }
//            return out;
//        }

//        /// сгенерить плотность цветков на заданной высоте
//        double generate_density2(double h) const { return generate_at(REGRESSION(density2, height), h); }

//        /// сгенерить сортированный список высот для заданного числа цветков (1|2)-го порядка
//        shared<curve<1>> generate_heights(double nf1, Unsigned axis)
//        {
//            shared<curve<1>> out(new curve<1>);
//            for (Unsigned i=nf1; i>=1; i--) out->push_back(relative_height(i, nf1));

//            if (axis == 2)
//            {
//                auto hrange = get_interval(REGRESSION(density2, height));
//                to_range(out, hrange, 0); // отрисовка только != 1 значений
//             }
//            return out;
//        }

//        /// сгенерить плотность цветков по списку высот и полному числу цветков 2-го порядка
//        shared<curve<2>> generate_density2(shared<curve<1>> heights, double nf2) const
//        {
//            shared<curve<2>> out(new curve<2>);
//            double sum = 0;
//            for (auto e : *heights)
//            {
//                double h = e[0];
//                double v = generate_density2(h);
//                out->push_back({h, v});
//                sum += v;
//            }
//            to_normalized(out, nf2);
//            round(out);

//            auto hrange = get_interval(REGRESSION(density2, height));
//            to_range(out, hrange, 0); // отрисовка только != 1 значений

//        #ifdef TEST_PLOT
//            auto f3 = smooth_regression(out, 0.2);
//            gnuplotParams param;
//            param.xlabel = "number of the 2nd-order flowers";
//            param.ylabel = "relative high";
//            param.insert_data(extract_errors_range(f3), "u 2:1 w filledcu closed");
//            param.insert_data(extract_mean(f3), "u 2:1 w l");
//            param.insert_data(out, "u 2:1 w p");
//            plot(param);
//        #endif

//            return out;
//        }

//        /// сгенерить дату старта цветения цветков 1-го порядка на заданной высоте
//        double generate_start1(double h) const { return generate_at(REGRESSION(start1, height), h); }

//        /// сгенерить кривую старта цветения цветков 1-го порядка
//        shared<curve<2>> generate_start1(shared<curve<1>> heights) const
//        {
//            shared<curve<2>> out(new curve<2>);
//            for (auto e : *heights)
//            {
//                double h = e[0];
//                double v = generate_start1(h);
//                out->push_back({h, v});
//            }
//            round(out);

//        #ifdef TEST_PLOT
//            auto f3 = smooth_regression(out, 0.2);
//            gnuplotParams param;
//            param.xlabel = "start of the 1st-order axis flowering";
//            param.ylabel = "relative high";
//            param.insert_data(extract_errors_range(f3), "u 2:1 w filledcu closed");
//            param.insert_data(extract_mean(f3), "u 2:1 w l");
//            param.insert_data(out, "u 2:1 w p");
//            plot(param);
//        #endif

//            return out;
//        }

//        /// сгенерить дату старта цветения цветков 1-го порядка на заданной высоте
//        double generate_start2(double h) const { return generate_at(REGRESSION(start2, height), h); }

//        /// сгенерить кривую старта цветения цветков 2-го порядка
//        shared<curve<2>> generate_start2(shared<curve<1>> heights) const
//        {
//            shared<curve<2>> out(new curve<2>);
//            for (auto e : *heights)
//            {
//                double h = e[0];
//                double v = generate_start2(h);
//                out->push_back({h, v});
//            }
//            round(out);

//        #ifdef TEST_PLOT
//            auto f3 = smooth_regression(out, 0.2);
//            gnuplotParams param;
//            param.xlabel = "relative start of the 2nd-order axis flowering";
//            param.ylabel = "relative high";
//            param.insert_data(extract_errors_range(f3), "u 2:1 w filledcu closed");
//            param.insert_data(extract_mean(f3), "u 2:1 w l");
//            param.insert_data(out, "u 2:1 w p");
//            plot(param);
//        #endif

//            return out;
//        }

//    public:

//        Sproutmodel() {}

//        /*!
//         * \brief build the model
//         * \param flowerdome база данных по которой строятся общие распределения по побегам
//         * \param sprout номер побега, по которому строятся модель конкретного побега
//         */
//        void build(const Flowerdome &flowerdome, Unsigned sprout);

//        /*!
//         * \brief print the list of all the distributions and regressions
//         * \note адаптировано под вывод в ncurses terminal
//         */
//        std::string view_list() const;


//        /*!
//         * \brief view the distributions and regressions
//         * \param идентификатор зависимости
//         */
//        std::string view(Unsigned id);
//        std::string view(const std::vector<Unsigned> &ids);

//        /*!
//         * \brief генерация случайного образца модельного побега
//         * \param restrictions ограничения по числу цветков различных порядков
//         *  [0] по полному числу цветков
//         *  [1] по числу цветков на 1й оси
//         *  [2] по числу цветков на каждой 2й оси ...
//         */
//        shared<Sprout> generate_sprout(const std::vector<Interval<Unsigned, 1>> &restrictions);
//        shared<Sprout> generate_sprout();

//    };

//    inline std::string Sproutmodel::view_list() const
//    {
//        if (dependencies_items.empty())
//            return "No data yet. Build the model [model sprout_id]";

//        std::string msg;
//        for (Unsigned ndx=1; ndx<dependencies_items.size(); ndx++)
//        {
//            Unsigned id = ndx;
//            std::string title = dependencies_items[ndx].title;
//            msg += '[' + std::to_string(id) + "] " + title + "\n";
//        }
//        return msg;
//    }

//    inline std::string Sproutmodel::view(Unsigned id)
//    {
//        if (dependencies_items.empty())
//            return "No data yet. Build the model [model sprout_id]";

//        gnuplotParams v = dependencies_items.at(id).gnuparam;
//        plot(v);
//        return "OK";
//    }

//    inline std::string Sproutmodel::view(const std::vector<Unsigned> &ids)
//    {
//        bool flag = true;
//        for (auto id : ids) if (dependencies.count(id) == 0) flag = false;
//        if (!flag) return "No data yet. Build the model [model sprout_id]";

//        gnuplotParams params = dependencies.at(ids.front());
//        for (Unsigned i=1; i<ids.size(); i++)
//        {
//            gnuplotParams v = dependencies.at(ids[i]);
//            if (v.xlabel != params.xlabel) return "fail";
//            if (v.ylabel != params.ylabel) return "fail";

//            for (Unsigned j=0; j<v.size(); j++)
//            {
//                auto data = v.get_data(j);
//                auto options = v.get_options(j);
//                params.insert_data(data, options);
//            }
//        }
//        plot(params);
//        return "OK";
//    }

//#undef DISTRIBUTION
//#undef REGRESSION


//        item.title = "flowering curves [1st]";
//        item.pf1 = &Sproutmodel::flowering1_by_time;
//        item.pf2 = &Sproutmodel::flowering1_by_time;
//        gnuparam.clear();
//        gnuparam.title = item.title;
//        gnuparam.xlabel = "T, days from the beginning of flowering";
//        gnuparam.ylabel = "n, open flowers";
//        (this->*(item.pf1))(item, flowerdome, {sprout});
//        (this->*(item.pf2))(item, default_width);
//        dependencies_items.push_back(item);

//        item.title = "flowering curves [2nd]";
//        item.pf1 = &Sproutmodel::flowering2_by_time;
//        item.pf2 = &Sproutmodel::flowering2_by_time;
//        gnuparam.clear();
//        gnuparam.title = item.title;
//        gnuparam.xlabel = "T, days from the beginning of flowering";
//        gnuparam.ylabel = "%, open flowers";
//        (this->*(item.pf1))(item, flowerdome, {sprout});
//        (this->*(item.pf2))(item, default_width);
//        dependencies_items.push_back(item);


        //        void density2_by_height(dependence_item &item, const Flowerdome &flowerdome,
        //            const std::set<Unsigned> &sprouts)
        //        {
        //            double hx = 0.01;

        //            shared<curve<2>> data(new curve<2>);
        //            for (auto sprout : sprouts)
        //            {
        //                const Sprout &s = flowerdome.at(sprout);
        //                Unsigned nf1 = s.number_of_axis1_flowers();

        //                auto tmp = flowerdome.foreach({sprout}, axis_density_functor(2));
        //                if (tmp->empty()) continue;

        //                unite(tmp, hx);
        //                to_height(tmp, nf1);

        //                for (auto e : *tmp) data->push_back(e);
        //            }
        //            item.data = data;
        //        }

        //        void density2_by_height(dependence_item &item, double width)
        //        {
        //            item.smooth = smooth_regression(item.data, width);

        //            gnuplotParams &gnuparam = item.gnuparam;

        //            gnuparam.insert_data(extract_errors_range(item.smooth),
        //                "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ");
        //            gnuparam.insert_data(extract_mean(item.smooth), "u 2:1 w l lw 5 lc 'blue' ");
        //            gnuparam.insert_data(item.data, "u 2:1 w p pt 7 lc rgb 'blue' ");

        //            REGRESSION(density2, height) = item.smooth;
        //        }

        //        std::string density2_by_height_test(dependence_item &item, Unsigned)
        //        {
        //            return regression_test(item);
        //        }

