#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include "defs.h"
#include "math.WLSM.h"
#include "math.gnuplot.h"
#include "flower.ensemble.h"
#include "plot.flowering.h"

/// при установке рисуются все кривые (1я ось, 2я ось, сумма и разность)
#define DRAW_ALL_CURVES

/// вывести только заданное число точек, чтобы не загромождать график
#define DRAW_ONLY 30

curve_t prepare_PDF(const curve_t &v1)
{
    //---------------------------------------------------------------------------------------------
    //                          рассчитаем функцию распределения
    //---------------------------------------------------------------------------------------------
    const unsigned nmid = 18; // число промежуточных точек
    double indent = std::max(0.25 * (v1.back()[0] - v1.front()[0]), double(v1.front()[0]));
//    double indentL = double(v1.front()[0]);
//    double indentR = 0.5 * (v1.back()[0] - v1.front()[0]);
    double indentL = indent;
    double indentR = indent;
    
    curve_t data;
    double prevx = v1.front()[0] - indentL, cnt = 0, dx = 0;
    for (unsigned i=0; i<v1.size(); i++)
    {
        dx = (v1[i][0] - prevx) / nmid;
        for (unsigned k=0; k<nmid; k++)
        {
            data.push_back({prevx, cnt});
            prevx += dx;
        }
        cnt += v1[i][1];
        prevx = v1[i][0];
    }
    dx = indentR / nmid;
    for (unsigned k=0; k<nmid; k++)
    {
        data.push_back({prevx, cnt});
        prevx += dx;
    }

    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    unsigned power = 2; // степень апроксимирующего полинома
    double hx = 0.125; // расстояния между точками при разложении по полиномам
    double wx1 = 0.125 * (data.back()[0] - data.front()[0]); // ширина окна усреднения

    double norm = 0;
    for (auto e : data) if (norm < e[1]) norm = e[1];
    norm = 1./norm;
    for (auto &e : data) { e[1] *= norm; e[2] *= norm; }

    //---------------------------------------------------------------------------------------------
    //                      расчет гистограммы цветения
    //---------------------------------------------------------------------------------------------
    auto f = wslm_simple_smooth<CDIM>(1, power, data, 0, 1, wx1, hx); // PDF

    //---------------------------------------------------------------------------------------------
    //                  нормировка на полное число цветков на площадь
    //---------------------------------------------------------------------------------------------
    norm = 0; for (auto e : f) norm += e[1];
    norm = 1. / (hx * norm);
    for (auto &e : f) { e[1] *= norm; e[2] *= norm; }
        // нормировали так, что полный интеграл равен 1

//    auto f = wslm_simple_smooth<CDIM>(0, power, data, 0, 1, wx1, hx); // CDF
//    norm = 1. / f.back()[1];
//    for (auto &e : f) { e[1] *= norm; e[2] *= norm; }

    return f;
}

void plot_countPDF(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2, const curve_t &vs)
{
    if (v1.empty() || v2.empty() || vs.empty()) return;

    //---------------------------------------------------------------------------------------------
    //                          рассчитаем функцию распределения
    //---------------------------------------------------------------------------------------------
    auto f1 = prepare_PDF(v1);
    auto f2 = prepare_PDF(v2);
    auto fs = prepare_PDF(vs);

    //---------------------------------------------------------------------------------------------
    //                  нормировка на полное число цветков на площадь
    //---------------------------------------------------------------------------------------------
    double max = 0; // максимум графиков
    for (auto e : f1) if (max < e[1]) max = e[1];
    for (auto e : f2) if (max < e[1]) max = e[1];
    for (auto e : fs) if (max < e[1]) max = e[1];

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::vector<curve_t> patterns { v1, v2, f1, f2, fs };
    std::vector<gnuplot_packet> options
    {
        "u 1:(0):(0):2 w vectors filled nohead lw 4 lt 0 lc 'red' title '' axis x1y2",
        "u 1:(0):(0):2 w vectors filled nohead lw 4 lt 0 lc 'blue' title '' axis x1y2",
        "u 1:2 w l lt 7 lw 3 fc 'red' title '{/Symbol Y}_1(n)' ",
        "u 1:2 w l lt 7 lw 3 fc 'blue' title '{/Symbol Y}_2(n)' ",
        "u 1:2 w l lt 0 lw 8 fc 'black' title '{/Symbol Y}(n)' "
    };

    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }

    //---------------------------------------------------------------------------------
    //                          старт реальной отрисовки
    //---------------------------------------------------------------------------------
    //pipe << "set key off";

    std::string title = species + ", PDFs ";
    pipe << gnuplot_packet("set title '" + title + "'");

    pipe << "set xlabel 'n, number of flowers'";
    pipe << "set xrange[-0.5:*]";

    gnuplot_packet pack("set yrange[0:");
    pack << 1.05 * max; // чуть выше для красоты
    pack << "]";
    pipe << pack;

    pipe << "set ytics nomirror";
    pipe << "set y2tics";
    pipe << "set y2range[0:*]";

    pipe << "set ylabel 'probability density'\n";
    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

void plot_countRF(const std::string &output, const std::string &species, const curve_t &data1)
{
    if (data1.empty()) return;

    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    unsigned power = 1; // степень апроксимирующего полинома
    double wx1 = 0.25 * (data1.back()[0] - data1.front()[0]); // ширина окна усреднения
    double hx = 0.5; // расстояния между точками при разложении по полиномам

    //---------------------------------------------------------------------------------------------
    //                      расчет гистограммы цветения
    //---------------------------------------------------------------------------------------------
    auto f1d = wslm_simple_smooth<CDIM>(0, power, data1, 0, 1, wx1, hx);
    auto closed1 = prg::closed_bound<CDIM>(f1d); // закрытые области, описывающие границы
    for (auto &e : closed1) if (e[1] < 0) e[1] = 0;

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::vector<curve_t> patterns
    {
        closed1,
        data1, f1d
    };

    std::vector<gnuplot_packet> options
    {
        "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder title '{/Symbol d}f_{n_2}(n_1)'",
        "u 1:2 w p pt 7 fc 'blue' title '\\{n_1,n_2\\}'",
        "u 1:2 w l lt 0 lw 8 fc 'blue' title 'f_{n_2}(n_1)' ",
    };

    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }

    //---------------------------------------------------------------------------------
    //                          старт реальной отрисовки
    //---------------------------------------------------------------------------------
    //pipe << "set key off";

    std::string title = species + ", RF";
    pipe << gnuplot_packet("set title '" + title);

    pipe << "set xlabel 'n_1, number of 1^{st} order flowers'";
    pipe << "set xrange[*:*]";

    pipe << "set yrange[*:*]";
    pipe << "set ylabel 'n_2, number of 2^{nd} order flowers'";
    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

void plot_axeRF(const std::string &output, const std::string &species,
    const curve_t &data1, const curve_t &data2, const std::string &xrange)
{
    if (data1.empty() || data2.empty()) return;

    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    unsigned power = 2; // степень апроксимирующего полинома
    double wy1 = .15; // ширина окна усреднения
    double wy2 = .4; // ширина окна усреднения
    double hy = 0.0125; // расстояния между точками при разложении по полиномам
    double thickness = THICKNESS; // доля объектов принадлежащих границе

    //---------------------------------------------------------------------------------
    //                      расчет данных для 1й оси
    //---------------------------------------------------------------------------------
    auto curve1 = wslm_simple_smooth<CDIM>(0, power, data1, 0, 1, wy1, hy);
        // кривая расцветания
    auto closed1 = prg::closed_bound<CDIM>(curve1);
        // закрытые области, описывающие границы и их погрешности

    //---------------------------------------------------------------------------------
    //                      расчет данных для 2й оси
    //---------------------------------------------------------------------------------
    auto ldata2 = lower_bound<CDIM>(data2, thickness);//fom

    auto curve2 = wslm_simple_smooth<CDIM>(0, power, ldata2, 0, 1, wy2, hy);
        // кривые расцветания на левой и правой границах, соответственно
    auto closed2 = prg::closed_bound<CDIM>(curve2);
        // закрытые области, описывающие границы и их погрешности

    //---------------------------------------------------------------------------------
    //              "облегчим" графики от части сырых данных
    //---------------------------------------------------------------------------------

    std::uniform_int_distribution<> d;
    typedef std::uniform_int_distribution<>::param_type param_type;
    curve_t data1_, data2_, ldata2_;

    if (DRAW_ONLY < data1.size())
    {
        d.param(param_type{0, int(data1.size())-1});
        for (Unsigned i=0; i<DRAW_ONLY; i++) data1_.push_back(data1[Unsigned(d(rangen))]);
        std::sort(data1_.begin(), data1_.end(),
            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
        );
    }
    else for (auto e : data1) data1_.push_back(e);

    if (DRAW_ONLY < data2.size())
    {
        d.param(param_type{0, int(data2.size())-1});
        for (Unsigned i=0; i<DRAW_ONLY; i++) data2_.push_back(data2[Unsigned(d(rangen))]);
        std::sort(data2_.begin(), data2_.end(),
            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
        );
    }
    else for (auto e : data2) data2_.push_back(e);

    if (DRAW_ONLY < ldata2.size())
    {
        d.param(param_type{0, int(ldata2.size())-1});
        for (Unsigned i=0; i<DRAW_ONLY; i++) ldata2_.push_back(ldata2[Unsigned(d(rangen))]);
        std::sort(ldata2_.begin(), ldata2_.end(),
            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
        );
    }
    else for (auto e : ldata2) ldata2_.push_back(e);

    //---------------------------------------------------------------------------------
    //                  обращаем данные по высоте
    //---------------------------------------------------------------------------------
    for (auto &e : closed1) e[0] = 1 - e[0];
    for (auto &e : closed2) e[0] = 1 - e[0];
    for (auto &e : curve1) e[0] = 1 - e[0];
    for (auto &e : curve2) e[0] = 1 - e[0];
    for (auto &e : data1_) e[0] = 1 - e[0];
    for (auto &e : data2_) e[0] = 1 - e[0];
    for (auto &e : ldata2_) e[0] = 1 - e[0];

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::vector<curve_t> patterns
    {
        closed1, closed2, // сперва закрытые области, чтобы избежать сокрытия данных
        curve1, curve2,
        data2_, data1_, ldata2_ // серые незначимы, потому допустимо перекрывать
    };

    std::vector<gnuplot_packet> options
    {
        "u 2:1 w filledcu closed fc 'red' fill transparent solid 0.2 noborder title '{/Symbol d}f_{t_1}(h)'",
        "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder title '{/Symbol d}f_{t_2}(h)'",
        "u 2:1 w l lc 'red' title 'f_{t_1}(h)' ",
        "u 2:1 w l lc 'blue' title 'f_{t_2}(h)' ",
        "u 2:1:4 w xerr pt 7 fc 'gray' title '' ",
        "u 2:1:4 w xerr pt 7 fc 'red' title '' ",
        "u 2:1:4 w xerr pt 7 fc 'blue' title '' "
    };
    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }

    pipe << "set ylabel 'h, rel.height of flower'\n";
//    pipe << "set key off";
    pipe << "set yrange[*:*]";

    std::string title = species;
    pipe << gnuplot_packet("set title '" + title);

    pipe << gnuplot_packet("set xrange" + xrange);

    pipe << "set xlabel 'T, days from the beginning of flowering'\n";
    //pipe << "set xtics ('0' 0, '10' 10, '20' 20, '30' 30, '40' 40, '50' 50, '60' 60)";
    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

void plot_densityRF(const std::string &output, const std::string &species,
    const curve_t &data1, const curve_t &data2)
{
    if (data1.empty() || data2.empty()) return;

    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    unsigned power = 2; // степень апроксимирующего полинома
    double hy = 0.0125; // расстояния между точками при разложении по полиномам
    double wyd = .15; // ширина окна усреднения для плотности

    //---------------------------------------------------------------------------------------------
    //                          расчет графика плотности
    //---------------------------------------------------------------------------------------------
    curve_t density1_data, density2_data;
    double max = 0;
    {
        typedef std::pair<double, unsigned> key_t;

        std::map<key_t, unsigned> sum;
        for (auto e : data1)
        {
            key_t key{e[0], e[5]};
            if (sum.count(key) == 0) sum[key] = 0;
            sum[key]++;
        }
        for (auto e : sum) density1_data.push_back({e.first.first, double(e.second)});

        sum.clear();
        for (auto e : data2)
        {
            key_t key{e[0], e[5]};
            if (sum.count(key) == 0) sum[key] = 0;
            sum[key]++;
        }
        for (auto e : sum) density2_data.push_back({e.first.first, double(e.second)});

        for (auto e : sum) if (max < e.second) max = e.second;
    }

    auto density1_curve = wslm_simple_smooth<CDIM>(0, power, density1_data, 0, 1, wyd, hy);
    auto density2_curve = wslm_simple_smooth<CDIM>(0, power, density2_data, 0, 1, wyd, hy);

    double thin = 0.05;
    for (auto &e : density1_curve) if (e[2] < thin) e[2] = thin;
    for (auto &e : density2_curve) if (e[2] < thin) e[2] = thin;
        // добавим толщину к линии
    
    auto density1_closed = prg::closed_bound<CDIM>(density1_curve);

    auto density2_closed = prg::closed_bound<CDIM>(density2_curve);
        // данные, кривая и закрытая область для описания плотности цветения 2й оси

//    density1_curve.insert(density1_curve.begin(), density1_curve.front());
//    density1_curve.front()[1] = 0;
//    density1_curve.insert(density1_curve.end(), density1_curve.back());
//    density1_curve.back()[1] = 0;
//        // вставка линий к 0 для красоты отображения
//    density2_curve.insert(density2_curve.begin(), density2_curve.front());
//    density2_curve.front()[1] = 0;
//    density2_curve.insert(density2_curve.end(), density2_curve.back());
//    density2_curve.back()[1] = 0;
//        // вставка линий к 0 для красоты отображения

    //---------------------------------------------------------------------------------
    //              "облегчим" графики от части сырых данных
    //---------------------------------------------------------------------------------
    std::uniform_int_distribution<> d;
    typedef std::uniform_int_distribution<>::param_type param_type;
    curve_t data1_, data2_, density2_data_;

    if (DRAW_ONLY < data1.size())
    {
        d.param(param_type{0, int(data1.size())-1});
        for (Unsigned i=0; i<DRAW_ONLY; i++) data1_.push_back(data1[Unsigned(d(rangen))]);
        std::sort(data1_.begin(), data1_.end(),
            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
        );
    }
    else for (auto e : data1) data1_.push_back(e);

    if (DRAW_ONLY < density2_data.size())
    {
        d.param(param_type{0, int(density2_data.size())-1});
        for (Unsigned i=0; i<DRAW_ONLY; i++) density2_data_.push_back(density2_data[Unsigned(d(rangen))]);
        std::sort(density2_data_.begin(), density2_data_.end(),
            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
        );
    }
    else for (auto e : density2_data) density2_data_.push_back(e);

    //---------------------------------------------------------------------------------
    //                  обращаем данные по высоте
    //---------------------------------------------------------------------------------
    for (auto &e : data1_) e[0] = 1 - e[0];
    for (auto &e : data2_) e[0] = 1 - e[0];
    for (auto &e : density1_curve) e[0] = 1 - e[0];
    for (auto &e : density2_curve) e[0] = 1 - e[0];
    for (auto &e : density2_data_) e[0] = 1 - e[0];
    for (auto &e : density1_closed) e[0] = 1 - e[0];
    for (auto &e : density2_closed) e[0] = 1 - e[0];

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::vector<curve_t> patterns
    {
        density1_closed, density2_closed,
        density1_curve, density2_curve,
        density2_data_
    };

    std::vector<gnuplot_packet> options
    {
        "u 2:1 w filledcu fc 'red' fill transparent solid 0.2 noborder title '{/Symbol d}f_{{/Symbol r}_1}(h)' ",
        "u 2:1 w filledcu fc 'blue' fill transparent solid 0.2 noborder title '{/Symbol d}f_{{/Symbol r}_2}(h)' ",
        "u 2:1 w l lt 0 lw 8 lc 'red' title 'f_{{/Symbol r}_1}(h)' ",
        "u 2:1 w l lt 0 lw 8 lc 'blue' title 'f_{{/Symbol r}_2}(h)' ",
        "u 2:1:3 w xerr pt 7 fc 'blue' title ''",
    };

    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }


    std::string title = species;
    pipe << gnuplot_packet("set title '" + title + "'");

    pipe << "set ylabel 'h, rel.height of flower'\n";
    //pipe << "set key off";
    pipe << "set yrange[-0.05:1.05]";

    //pipe << "set xrange[0.0:*]";
    gnuplot_packet pack("set xrange[0:");
    pack << max + 1 << "]";
    pipe << pack;

    pipe << "set xlabel '{/Symbol s}, flowering density'";
    pipe << "set xtics ('0' 0, '1' 1, '2' 2, '3' 3, '4' 4, '5' 5, '6' 6, '7' 7, '8' 8)";

    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

void plot_spanRF(const std::string &output, const std::string &species,
    const curve_t &data1, const curve_t &data2)
{
    if (data1.empty() || data2.empty()) return;

    double max = 0;
    for (auto e : data1) if (max < e[2]) max = e[2];
    for (auto e : data2) if (max < e[2]) max = e[2];

    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    unsigned power = 2; // степень апроксимирующего полинома
    double wy = .2; // ширина окна усреднения
    double hy = 0.0125; // расстояния между точками при разложении по полиномам
    //---------------------------------------------------------------------------------------------

    //---------------------------------------------------------------------------------------------
    //                      расчет графиков длительности цветения
    //---------------------------------------------------------------------------------------------
    std::set<double> heights; // набор уникальных высот
    for (auto e : data2) heights.insert(e[0]);

    // на 2й оси используем накопление по всем цветкам одной высоты
    curve_t flowerspan2;
    for (auto h : heights)
    {
        auto first = std::lower_bound(data2.begin(), data2.end(), h,
            [](const point<CDIM> &v, double h){ return v[0] < h; });
        auto last = std::upper_bound(data2.begin(), data2.end(), h,
            [](double h, const point<CDIM> &v){ return h < v[0]; });
        assert(first != last);

//        std::vector<double> tmp;
//        std::vector<double> eps; // погрешность длительности цветения
//        for (; first!=last; ++first)
//        {
//            tmp.push_back((*first)[2]);
//            eps.push_back((*first)[4]);
//        }

//        double mean = average(tmp);
//        double var = variance(tmp, mean);
//        double meane = average(eps);
//        double vare = variance(eps, meane);

//        double sigma2 = var + vare;
//        flowerspan2.push_back({h, 0, mean, 0, sqrt(sigma2)});
//            // порядок как у data1

        std::vector<double> tmp;
        double eps2 = 0; // погрешность длительности цветения
        for (; first!=last; ++first)
        {
            tmp.push_back((*first)[2]);
            eps2 += sqr((*first)[4]);
        }

        double mean = average(tmp);
        double var = variance(tmp, mean);

        double sigma2 = var + eps2;
        flowerspan2.push_back({h, 0, mean, 0, sqrt(sigma2)});
            // порядок как у data1
    }

    auto fcurve1 = wslm_simple_smooth<CDIM>(0, power, data1, 0, 2, wy, hy);
    auto fcurve2 = wslm_simple_smooth<CDIM>(0, power, flowerspan2, 0, 2, wy, hy);
        // кривые расцветания на левой и правой границах, соответственно

    double thin = 0.05;
    for (auto &e : fcurve1) if (e[2] < thin) e[2] = thin;
    for (auto &e : fcurve2) if (e[2] < thin) e[2] = thin;
        // добавим толщину к линии

    auto closed1 = prg::closed_bound<CDIM>(fcurve1);
    auto closed2 = prg::closed_bound<CDIM>(fcurve2);
        // закрытые области, описывающие границы и их погрешности

//    fcurve1.insert(fcurve1.begin(), fcurve1.front());
//    fcurve1.front()[1] = 0;
//    fcurve1.insert(fcurve1.end(), fcurve1.back());
//    fcurve1.back()[1] = 0;
//        // вставка линий к 0 для красоты отображения
//    fcurve2.insert(fcurve2.begin(), fcurve2.front());
//    fcurve2.front()[1] = 0;
//    fcurve2.insert(fcurve2.end(), fcurve2.back());
//    fcurve2.back()[1] = 0;
//        // вставка линий к 0 для красоты отображения


//    for (auto &e : fcurve1) e[1] += e[2];
//    for (auto &e : fcurve2) e[1] += e[2];
        // кривые расцветания на левой и правой границах вместе с погрешностью аппроксимации

    //---------------------------------------------------------------------------------
    //              "облегчим" графики от части сырых данных
    //---------------------------------------------------------------------------------
    std::uniform_int_distribution<> d;
    typedef std::uniform_int_distribution<>::param_type param_type;
    curve_t data1_, data2_, density2_data_;

    if (DRAW_ONLY < data1.size())
    {
        d.param(param_type{0, int(data1.size())-1});
        for (Unsigned i=0; i<DRAW_ONLY; i++) data1_.push_back(data1[Unsigned(d(rangen))]);
        std::sort(data1_.begin(), data1_.end(),
            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
        );
    }
    else for (auto e : data1) data1_.push_back(e);

    if (DRAW_ONLY < flowerspan2.size())
    {
        d.param(param_type{0, int(flowerspan2.size())-1});
        for (Unsigned i=0; i<DRAW_ONLY; i++) data2_.push_back(flowerspan2[Unsigned(d(rangen))]);
        std::sort(data2_.begin(), data2_.end(),
            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
        );
    }
    else for (auto e : flowerspan2) data2_.push_back(e);

    //---------------------------------------------------------------------------------
    //                  обращаем данные по высоте
    //---------------------------------------------------------------------------------
    for (auto &e : fcurve1) e[0] = 1 - e[0];
    for (auto &e : fcurve2) e[0] = 1 - e[0];
    for (auto &e : closed1) e[0] = 1 - e[0];
    for (auto &e : closed2) e[0] = 1 - e[0];
    for (auto &e : data1_) e[0] = 1 - e[0];
    for (auto &e : data2_) e[0] = 1 - e[0];

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::vector<curve_t> patterns
    {
        closed1, closed2, 
        fcurve1, fcurve2, 
        data1_, data2_
    };

    std::vector<gnuplot_packet> options
    {
        "u 2:1 w filledcu y1 fc 'red' fill transparent solid 0.2 noborder title 'f_{{/Symbol t}_1}(h)'",
        "u 2:1 w filledcu y1 fc 'blue' fill transparent solid 0.2 noborder title 'f_{{/Symbol t}_2}(h)'",
        "u 2:1 w l lw 8 lt 0 lc 'red' title '{/Symbol d}f_{{/Symbol t}_1}(h)'",
        "u 2:1 w l lw 8 lt 0 lc 'blue' title '{/Symbol d}f_{{/Symbol t}_2}(h)'",
        "u 3:1:5 w xerr pt 7 fc 'red' title '' ",
        "u 3:1:5 w xerr pt 7 fc 'blue' title ''"
    };

    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }

    pipe << "set ylabel 'h, rel.height of flower'\n";
    //pipe << "set key off";
    pipe << "set yrange[-0.05:1.05]";

    std::string title = species;
    pipe << gnuplot_packet("set title '" + title + "'");

    //pipe << "set xrange[0:*]";
    gnuplot_packet pack("set xrange[0:");
    pack << max + 2 << "]";
    pipe << pack;

    pipe << "set xlabel '{/Symbol t}, flowering time '";
    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

void plot_spanRF(const std::string &output, const std::string &species,
    const std::vector<curve_t> &data1, const std::vector<curve_t> &data2)
{
    if (data1.empty() || data2.empty()) return;
    Unsigned count = data1.size();

    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    unsigned power = 2; // степень апроксимирующего полинома
    double wy = .2; // ширина окна усреднения
    double hy = 0.0125; // расстояния между точками при разложении по полиномам

    double max = 0;
    for (Unsigned i=0; i<count; i++)
    {
        for (auto e : data1[i]) if (max < e[2]) max = e[2];
        for (auto e : data2[i]) if (max < e[2]) max = e[2];
    }

    //---------------------------------------------------------------------------------------------
    //                      расчет графиков длительности цветения
    //---------------------------------------------------------------------------------------------
    std::vector<curve_t> FCURVE1, FCURVE2, FCLOSED1, FCLOSED2;
    for (Unsigned i=0; i<count; i++)
    {
        std::set<double> heights; // набор уникальных высот
        for (auto e : data2[i]) heights.insert(e[0]);
    
        // на 2й оси используем накопление по всем цветкам одной высоты
        curve_t flowerspan2;
        for (auto h : heights)
        {
            auto first = std::lower_bound(data2[i].begin(), data2[i].end(), h,
                [](const point<CDIM> &v, double h){ return v[0] < h; });
            auto last = std::upper_bound(data2[i].begin(), data2[i].end(), h,
                [](double h, const point<CDIM> &v){ return h < v[0]; });
            assert(first != last);
    
            std::vector<double> tmp;
            for (; first!=last; ++first) tmp.push_back((*first)[2]);
    
            flowerspan2.push_back({h, 0, average(tmp), 0, 0});
                // порядок как у data1
        }
    
        auto fcurve1 = wslm_simple_smooth<CDIM>(0, power, data1[i], 0, 2, wy, hy);
        
        unsigned power2 = power; if (flowerspan2.size() <= 3) power2 = 1;
        auto fcurve2 = wslm_simple_smooth<CDIM>(0, power2, flowerspan2, 0, 2, wy, hy);
            // кривые расцветания на левой и правой границах, соответственно
    
        //---------------------------------------------------------------------------------
        //                  обращаем данные по высоте
        //---------------------------------------------------------------------------------
        for (auto &e : fcurve1) e[0] = 1 - e[0];
        for (auto &e : fcurve2) e[0] = 1 - e[0];

        auto closed1 = prg::closed_bound<CDIM>(fcurve1);
        auto closed2 = prg::closed_bound<CDIM>(fcurve2);
        
        FCURVE1.push_back(fcurve1);
        FCURVE2.push_back(fcurve2);
        FCLOSED1.push_back(closed1);
        FCLOSED2.push_back(closed2);
    }

    Unsigned sz1 = FCURVE1[0].size();
    curve_t AVERAGE1(sz1, {0, 0});

    for (Unsigned k=0; k<sz1; k++)
    {
        AVERAGE1[k][0] = FCURVE1[0][k][0];
        for (Unsigned i=0; i<count; i++) AVERAGE1[k][1] += FCURVE1[i][k][1];
        AVERAGE1[k][1] /= count;
    }
        
    Unsigned sz2 = FCURVE2[0].size();
    curve_t AVERAGE2(sz2, {0, 0});

    for (Unsigned k=0; k<sz2; k++)
    {
        AVERAGE2[k][0] = FCURVE2[0][k][0];
        for (Unsigned i=0; i<count; i++) AVERAGE2[k][1] += FCURVE2[i][k][1];
        AVERAGE2[k][1] /= count;
    }

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::vector<curve_t> patterns;
//    for (Unsigned i=0; i<count; i++) patterns.push_back(FCURVE1[i]);
//    for (Unsigned i=0; i<count; i++) patterns.push_back(FCURVE2[i]);
    for (Unsigned i=0; i<count; i++) patterns.push_back(FCLOSED1[i]);
    for (Unsigned i=0; i<count; i++) patterns.push_back(FCLOSED2[i]);
    patterns.push_back(AVERAGE1);
    patterns.push_back(AVERAGE2);

    std::vector<gnuplot_packet> options;
//    for (Unsigned i=0; i<count; i++) options.push_back("u 2:1 w l lc 'red' lw 1 title ''");
//    for (Unsigned i=0; i<count; i++) options.push_back("u 2:1 w l lc 'blue' lw 1 title ''");

    const char *opt1 = "u 2:1 w filledcu y1 fc 'red' fill transparent solid 0.01 noborder title ''";
    const char *opt2 = "u 2:1 w filledcu y1 fc 'blue' fill transparent solid 0.01 noborder title ''";
    for (Unsigned i=0; i<count; i++) options.push_back(opt1);
    for (Unsigned i=0; i<count; i++) options.push_back(opt2);
    options.push_back("u 2:1 w l lt 0 lc 'black' lw 8 title '{/Symbol d}f_{{/Symbol t}_1}(h)'");
    options.push_back("u 2:1 w l lt 0 lc 'black' lw 8 title '{/Symbol d}f_{{/Symbol t}_2}(h)'");

    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }

    pipe << "set ylabel 'h, rel.height of flower'\n";
    pipe << "set yrange[-0.05:1.05]";

    std::string title = species;
    pipe << gnuplot_packet("set title '" + title + "'");

    gnuplot_packet pack("set xrange[0:");
    pack << max + 2 << "]";
    pipe << pack;

    pipe << "set xlabel '{/Symbol t}, flowering time '";
    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

void plot_spanPDF(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2)
{
    if (v1.empty() || v2.empty()) return;

    curve_t v1_, v2_, fs;
    std::map<double, unsigned> map;
    for (auto e : v1)
    {
        unsigned dt = unsigned(e[2]);
        if (map.count(dt) == 0) map[dt] = 0;
        ++map[dt];
    }
    for (auto e : map)
    {
        v1_.push_back({e.first, double(e.second)});
        fs.push_back({e.first, double(e.second)});
    }

    map.clear();
    for (auto e : v2)
    {
        unsigned dt = unsigned(e[2]);
        if (map.count(dt) == 0) map[dt] = 0;
        ++map[dt];
    }
    for (auto e : map)
    {
        v2_.push_back({e.first, double(e.second)});
        fs.push_back({e.first, double(e.second)});
    }
    std::sort(fs.begin(), fs.end(),
        [](const point_t &a, const point_t &b){ return a[0] < b[0]; }
    );

    auto f1 = prepare_PDF(v1_);
    auto f2 = prepare_PDF(v2_);
    fs = prepare_PDF(fs);

    //---------------------------------------------------------------------------------------------
    //                  нормировка на полное число цветков на площадь
    //---------------------------------------------------------------------------------------------
    double max = 0; // максимум графиков
    for (auto e : f1) if (max < e[1]) max = e[1];
    for (auto e : f2) if (max < e[1]) max = e[1];
    for (auto e : fs) if (max < e[1]) max = e[1];

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::vector<curve_t> patterns
    {
        v1_, v2_, f1, f2, fs
    };

    std::vector<gnuplot_packet> options
    {
        "u 1:(0):(0):2 w vectors filled nohead lw 4 lt 0 lc 'red' title '' axis x1y2",
        "u 1:(0):(0):2 w vectors filled nohead lw 4 lt 0 lc 'blue' title '' axis x1y2",
        "u 1:2 w l lt 7 lw 3 fc 'red' title '{/Symbol Y}_1({/Symbol t})'",
        "u 1:2 w l lt 7 lw 3 fc 'blue' title '{/Symbol Y}_2({/Symbol t})'",
        "u 1:2 w l lt 0 lw 8 fc 'black' title '{/Symbol Y}({/Symbol t})'"
    };

    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }

    //---------------------------------------------------------------------------------
    //                          старт реальной отрисовки
    //---------------------------------------------------------------------------------
    //pipe << "set key off";

    std::string title = species + ", PDFs ";
    pipe << gnuplot_packet("set title '" + title + "'");

    pipe << "set xlabel '{/Symbol t}, flowering time '";
    pipe << "set xrange[0:*]";

    gnuplot_packet pack("set yrange[0:");
    pack << 1.05 * max; // чуть выше для красоты
    pack << "]";
    pipe << pack;

    pipe << "set ytics nomirror";
    pipe << "set y2tics";

    pipe << "set ylabel 'probability density'\n";
    pipe << "set y2range[0:*]";
    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

void plot_spanPDF(const std::string &output, const std::string &species,
    const std::vector<curve_t> &v1, const std::vector<curve_t> &v2)
{
    if (v1.empty() || v2.empty()) return;
    Unsigned count = v1.size();
    double max = 0; // максимум графиков
    
    std::vector<curve_t> F1, F2;
    for (Unsigned i=0; i<count; i++)
    {
        curve_t v1_, v2_, fs;
        std::map<double, unsigned> map;
        for (auto e : v1[i])
        {
            unsigned dt = unsigned(e[2]);
            if (map.count(dt) == 0) map[dt] = 0;
            ++map[dt];
        }
        for (auto e : map)
        {
            v1_.push_back({e.first, double(e.second)});
            fs.push_back({e.first, double(e.second)});
        }
    
        map.clear();
        for (auto e : v2[i])
        {
            unsigned dt = unsigned(e[2]);
            if (map.count(dt) == 0) map[dt] = 0;
            ++map[dt];
        }
        for (auto e : map)
        {
            v2_.push_back({e.first, double(e.second)});
            fs.push_back({e.first, double(e.second)});
        }
        std::sort(fs.begin(), fs.end(),
            [](const point_t &a, const point_t &b){ return a[0] < b[0]; }
        );
    
        auto f1 = prepare_PDF(v1_);
        auto f2 = prepare_PDF(v2_);
        fs = prepare_PDF(fs);
    
        //---------------------------------------------------------------------------------------------
        //                  нормировка на полное число цветков на площадь
        //---------------------------------------------------------------------------------------------
        for (auto e : f1) if (max < e[1]) max = e[1];
        for (auto e : f2) if (max < e[1]) max = e[1];
        for (auto e : fs) if (max < e[1]) max = e[1];
        
        F1.push_back(f1);
        F2.push_back(f2);
    }

    double xmax = 0;
    for (Unsigned i=0; i<count; i++) 
    {
        if (xmax < F1[i].back()[0]) xmax = F1[i].back()[0];
        if (xmax < F2[i].back()[0]) xmax = F2[i].back()[0];
    }
    
    Unsigned sz = 100;
    double dx = xmax / sz;
    
    curve_t AVERAGE1(sz, {0, 0}), AVERAGE2(sz, {0, 0});
    for (Unsigned k=0; k<sz; k++)
    {
        double x = dx * k;
        AVERAGE1[k][0] = x;
        AVERAGE2[k][0] = x;

        for (Unsigned i=0; i<count; i++) 
        {
            auto it = std::lower_bound(F1[i].begin(), F1[i].end(), x, 
                [](const point_t &a, double x) { return a[0] < x; }
            );
            if (it != F1[i].end()) AVERAGE1[k][1] += (*it)[1];

            it = std::lower_bound(F2[i].begin(), F2[i].end(), x, 
                [](const point_t &a, double x) { return a[0] < x; }
            );
            if (it != F2[i].end()) AVERAGE2[k][1] += (*it)[1];
        }
        AVERAGE1[k][1] /= count;
        AVERAGE2[k][1] /= count;
    }

    unsigned power = 3;
    double wx = 0.5;
    double hx = 0.05;
    AVERAGE1 = wslm_simple_smooth<CDIM>(0, power, AVERAGE1, 0, 1, wx, hx);
    AVERAGE2 = wslm_simple_smooth<CDIM>(0, power, AVERAGE2, 0, 1, wx, hx);
    
    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::vector<curve_t> patterns;
//    for (Unsigned i=0; i<count; i++) patterns.push_back(F1[i]);
//    for (Unsigned i=0; i<count; i++) patterns.push_back(F2[i]);

    for (Unsigned i=0; i<count; i++) { patterns.push_back(F1[i]); patterns.push_back(F2[i]); }
    patterns.push_back(AVERAGE1);
    patterns.push_back(AVERAGE2);
    
    std::string opt1("u 1:2 w l lt 7 lw 1 fc rgb "); 
    opt1 += gnu::Pink1;  
    opt1 += " title ''";      
               
    std::string opt2("u 1:2 w l lt 7 lw 1 fc rgb "); 
    opt2 += gnu::LightBlue1;  
    opt2 += " title ''";      
    
    std::vector<gnuplot_packet> options;
//    for (Unsigned i=0; i<count; i++) options.push_back(opt1);
//    for (Unsigned i=0; i<count; i++) options.push_back(opt2);

    for (Unsigned i=0; i<count; i++) { options.push_back(opt1); options.push_back(opt2); }
    options.push_back("u 1:2 w l lt 0 lw 8 fc 'red' title '{/Symbol Y}_1({/Symbol t})'");
    options.push_back("u 1:2 w l lt 0 lw 8 fc 'blue' title '{/Symbol Y}_2({/Symbol t})'");
    
//    {
//        "u 1:2 w l lt 7 lw 3 fc 'red' title '{/Symbol Y}_1({/Symbol t})'",
//        "u 1:2 w l lt 7 lw 3 fc 'blue' title '{/Symbol Y}_2({/Symbol t})'",
//    };

    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }

    //---------------------------------------------------------------------------------
    //                          старт реальной отрисовки
    //---------------------------------------------------------------------------------
    std::string title = species + ", PDFs ";
    pipe << gnuplot_packet("set title '" + title + "'");

    pipe << "set xlabel '{/Symbol t}, flowering time '";
    pipe << "set xrange[0:*]";

    gnuplot_packet pack("set yrange[0:");
    pack << 1.05 * max; // чуть выше для красоты
    pack << "]";
    pipe << pack;

    pipe << "set ytics nomirror";
    pipe << "set y2tics";

    pipe << "set ylabel 'probability density'\n";
    pipe << "set y2range[0:*]";
    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

void plot_phase_range(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2, const std::string &xrange)
{
    if (v1.empty() || v2.empty()) return;

    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    unsigned power = 3; // степень апроксимирующего полинома
    double wy = 0.15; // ширина окна усреднения границ для кривых
    double wy2 = 0.6; // ширина окна усреднения границ для 2й оси
    double hy = 0.01; // расстояния между точками при разложении по полиномам

    //---------------------------------------------------------------------------------------------
    //                      пространственные области цветения
    //---------------------------------------------------------------------------------------------
    curve_t data1, data2;
    for (auto e : v1) data1.push_back(e);
    for (auto e : v2) data2.push_back(e);

    double shifty = .5 / data1.size(); // смещение по вертикали при отрисовке векторов для 2й оси
    for (auto &e : data2) e[0] -= shifty;

    // построим правые границы цветения (старт + продолжительность)
    curve_t data1R, data2R;
    for (auto e : data1) data1R.push_back({e[0], e[1] + e[2], sqrt(sqr(e[3]) + sqr(e[4]))});
    for (auto e : data2) data2R.push_back({e[0], e[1] + e[2], sqrt(sqr(e[3]) + sqr(e[4]))});

    auto L1 = wslm_simple_smooth<CDIM>(0, power, data1,  0, 1, wy, hy);
    auto R1 = wslm_simple_smooth<CDIM>(0, power, data1R, 0, 1, wy, hy);
        // сглаженные правые и левые границы диапазона цветения 1й оси

    auto xL1 = prg::subtract<CDIM>(prg::select<CDIM>(L1, 0, 1), prg::select<CDIM>(L1, 0, 2));
    auto xR1 = prg::summarize<CDIM>(prg::select<CDIM>(R1, 0, 1), prg::select<CDIM>(R1, 0, 2));
        // уширенные на стандартное отклонение правые и левые границы диапазона цветения 1й оси
    auto range_closed1 = prg::closed_bound<CDIM>(xL1, xR1);
        // полная область цветения 1й оси

    auto BL2 = prg::lower_bound<CDIM>(data2.begin(), data2.end(), wy2, 0, 1);
    auto BR2 = prg::upper_bound<CDIM>(data2R.begin(), data2R.end(), wy2, 0, 1);
        // правые и левые границы диапазона цветения 2й оси

    auto L2 = wslm_simple_smooth<CDIM>(0, power, BL2, 0, 1, wy, hy);
    auto R2 = wslm_simple_smooth<CDIM>(0, power, BR2, 0, 1, wy, hy);

        // сглаженные правые и левые границы диапазона цветения 2й оси
    auto xL2 = prg::subtract<CDIM>(prg::select<CDIM>(L2, 0, 1), prg::select<CDIM>(L2, 0, 2));
    auto xR2 = prg::summarize<CDIM>(prg::select<CDIM>(R2, 0, 1), prg::select<CDIM>(R2, 0, 2));
        // уширенные на стандартное отклонение правые и левые границы диапазона цветения 2й оси
    auto range_closed2 = prg::closed_bound<CDIM>(xL2, xR2);

    //---------------------------------------------------------------------------------
    //              "облегчим" графики от части сырых данных
    //---------------------------------------------------------------------------------
    std::uniform_int_distribution<> d;
    typedef std::uniform_int_distribution<>::param_type param_type;
    curve_t data1_, data2_;

    if (DRAW_ONLY < data1.size())
    {
        d.param(param_type{0, int(data1.size())-1});
        for (Unsigned i=0; i<DRAW_ONLY; i++) data1_.push_back(data1[Unsigned(d(rangen))]);
        std::sort(data1_.begin(), data1_.end(),
            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
        );
    }
    else for (auto e : data1) data1_.push_back(e);

    if (DRAW_ONLY < data2.size())
    {
        d.param(param_type{0, int(data2.size())-1});
        for (Unsigned i=0; i<DRAW_ONLY; i++) data2_.push_back(data2[Unsigned(d(rangen))]);
        std::sort(data2_.begin(), data2_.end(),
            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
        );
    }
    else for (auto e : data2) data2_.push_back(e);

    //---------------------------------------------------------------------------------
    //                  обращаем данные по высоте
    //---------------------------------------------------------------------------------
    for (auto &e : range_closed1) e[0] = 1 - e[0];
    for (auto &e : range_closed2) e[0] = 1 - e[0];
    for (auto &e : data1_) e[0] = 1 - e[0];
    for (auto &e : data2_) e[0] = 1 - e[0];

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::vector<curve_t> patterns
    {
        range_closed1, range_closed2, // сглаженные области цветения от начала до конца
        data1_, data2_ // точки реальных данных старта цветения
    }; // сперва закрытые области, чтобы избежать сокрытия данных

    std::vector<gnuplot_packet> options
    {
        "u 2:1 w filledcu closed fc 'red' fill transparent solid 0.2 noborder title '\\{t,{/Symbol t}\\}_1(h)' ",
        "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder title '\\{t,{/Symbol t}\\}_2(h)'",
        "u 2:1:3:(0) w vectors filled lw 5 fc 'red' title '{/Symbol t}_1(h)' ",
        "u 2:1:3:(0) w vectors filled lw 5 fc 'blue' title '{/Symbol t}_2(h)' "
    };
    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }
    //pipe << "set key off";

    std::string title = species + ", flowering ranges";
    //---------------------------------------------------------------------------------
    //                          старт реальной отрисовки
    //---------------------------------------------------------------------------------
    pipe << gnuplot_packet("set title '" + title + "'");

    pipe << gnuplot_packet("set xrange" + xrange);
    pipe << "set xlabel 'T, days from the beginning of flowering'";

    pipe << "set yrange[-0.02:1.02]";
    pipe << "set ylabel 'h, rel.height'";
    pipe << "set ytics";

    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, 4);

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

void plot_phase_range(const std::string &output, const std::string &species,
    const std::vector<curve_t> &v1, const std::vector<curve_t> &v2, const std::string &xrange)
{
    if (v1.empty() || v2.empty()) return;

    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    unsigned power = 3; // степень апроксимирующего полинома
    double wy = 0.15; // ширина окна усреднения границ для кривых
    double wy2 = 0.6; // ширина окна усреднения границ для 2й оси
    double hy = 0.01; // расстояния между точками при разложении по полиномам

    //---------------------------------------------------------------------------------------------
    //                      пространственные области цветения
    //---------------------------------------------------------------------------------------------
    std::vector<curve_t> data1, data2;
    for (auto e : v1) data1.push_back(e);
    for (auto e : v2) data2.push_back(e);

    Unsigned count = v1.size();
    std::vector<curve_t> CLOSED1, CLOSED2;

    auto L = wslm_simple_smooth<CDIM>(0, power, data1[0],  0, 1, wy, hy);
    Unsigned sz = L.size();
    curve_t AL1(L), AR1(L), AL2(L), AR2(L);
    for (Unsigned i=0; i<sz; i++)
    {
        AL1[i][1] = AR1[i][1] = AL2[i][1] = AR2[i][1] = 0;
        AL1[i][2] = AR1[i][2] = AL2[i][2] = AR2[i][2] = 0;
    }

    for (Unsigned i=0; i<count; i++)
    {
        // построим правые границы цветения (старт + продолжительность)
        curve_t data1R, data2R;
        for (auto e : data1[i]) data1R.push_back({e[0], e[1] + e[2], sqrt(sqr(e[3]) + sqr(e[4]))});
        for (auto e : data2[i]) data2R.push_back({e[0], e[1] + e[2], sqrt(sqr(e[3]) + sqr(e[4]))});

        auto L1 = wslm_simple_smooth<CDIM>(0, power, data1[i],  0, 1, wy, hy);
        auto R1 = wslm_simple_smooth<CDIM>(0, power, data1R, 0, 1, wy, hy);
            // сглаженные правые и левые границы диапазона цветения 1й оси

        if (!L1.empty() && !R1.empty())
        {
            for (Unsigned j=0; j<L1.size(); ++j)
            {
                double y = L1[j][0];
                Unsigned k = Unsigned(std::round(y / hy));
                if (k < sz)
                {
                    AL1[k][1] += L1[j][1]; AL1[k][2] += 1;
                    AR1[k][1] += R1[j][1]; AR1[k][2] += 1;
                }
            }
        }


        auto range_closed1 = prg::closed_bound<CDIM>(L1, R1);
            // полная область цветения 1й оси

        auto BL2 = prg::lower_bound<CDIM>(data2[i].begin(), data2[i].end(), wy2, 0, 1);
        auto BR2 = prg::upper_bound<CDIM>(data2R.begin(), data2R.end(), wy2, 0, 1);
            // правые и левые границы диапазона цветения 2й оси

        auto L2 = wslm_simple_smooth<CDIM>(0, power, BL2, 0, 1, wy2, hy);
        auto R2 = wslm_simple_smooth<CDIM>(0, power, BR2, 0, 1, wy2, hy);
            // сглаженные правые и левые границы диапазона цветения 2й оси

        if (!L2.empty() && !R2.empty())
        {
            for (Unsigned j=0; j<L2.size(); ++j)
            {
                double y = L2[j][0];
                Unsigned k = Unsigned(std::round(y / hy));
                if (k < sz)
                {
                    AL2[k][1] += L2[j][1]; AL2[k][2] += 1;
                    AR2[k][1] += R2[j][1]; AR2[k][2] += 1;
                }
            }
        }

        auto range_closed2 = prg::closed_bound<CDIM>(L2, R2);

        //---------------------------------------------------------------------------------
        //                  обращаем данные по высоте
        //---------------------------------------------------------------------------------
        for (auto &e : range_closed1) e[0] = 1 - e[0];
        for (auto &e : range_closed2) e[0] = 1 - e[0];

        CLOSED1.push_back(range_closed1);
        CLOSED2.push_back(range_closed2);
    }

    for (auto &e : AL1) { e[0] = 1 - e[0]; if (e[2] > 1) e[1] /= e[2]; }
    for (auto &e : AR1) { e[0] = 1 - e[0]; if (e[2] > 1) e[1] /= e[2]; }
    for (auto &e : AL2) { e[0] = 1 - e[0]; if (e[2] > 1) e[1] /= e[2]; }
    for (auto &e : AR2) { e[0] = 1 - e[0]; if (e[2] > 1) e[1] /= e[2]; }

    auto it = std::stable_partition(AL2.begin(), AL2.end(), [](const point_t &a){ return a[2] >= 2; });
    AL2.erase(it, AL2.end());
    it = std::stable_partition(AR2.begin(), AR2.end(), [](const point_t &a){ return a[2] >= 2; });
    AR2.erase(it, AR2.end());

    std::sort(AL2.begin(), AL2.end(), [](const point_t &a, const point_t &b){ return a[0] < b[0]; });
    std::sort(AR2.begin(), AR2.end(), [](const point_t &a, const point_t &b){ return a[0] < b[0]; });
    AL2 = wslm_simple_smooth<CDIM>(0, power, AL2, 0, 1, wy2, hy);
    AR2 = wslm_simple_smooth<CDIM>(0, power, AR2, 0, 1, wy2, hy);

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::swap(CLOSED1.front(), CLOSED1.back());
    std::swap(CLOSED2.front(), CLOSED2.back());

    std::vector<curve_t> patterns;
    for (Unsigned i=0; i<count; i++) patterns.push_back(CLOSED1[i]);

//    patterns.push_back(AL1);
//    patterns.push_back(AR1);
//    curve_t LINE1U, LINE1D;
//    LINE1U.push_back(AL1.front());
//    LINE1U.push_back(AR1.front());
//    LINE1D.push_back(AL1.back());
//    LINE1D.push_back(AR1.back());
//    patterns.push_back(LINE1U);
//    patterns.push_back(LINE1D);

    for (Unsigned i=0; i<count; i++) patterns.push_back(CLOSED2[i]);
//    patterns.push_back(AL2);
//    patterns.push_back(AR2);
//    curve_t LINE2U, LINE2D;
//    LINE2U.push_back(AL2.front());
//    LINE2U.push_back(AR2.front());
//    LINE2D.push_back(AL2.back());
//    LINE2D.push_back(AR2.back());
//    patterns.push_back(LINE2U);
//    patterns.push_back(LINE2D);

    std::string op1 = "u 2:1 w filledcu closed fc 'red' fill transparent solid 0.005 noborder title ''";
    std::string op2 = "u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.005 noborder title ''";

    std::vector<gnuplot_packet> options;
    //options.push_back("u 2:1 w filledcu closed fc 'red' fill transparent solid 0.02 noborder title '\\{t,{/Symbol t}\\}_{m_1}(h)' ");

    for (Unsigned i=0; i<count; i++) options.push_back(op1);
//    options.push_back("u 2:1 w l lc 'red' lw 8 lt 0 title '' ");
//    options.push_back("u 2:1 w l lc 'red' lw 8 lt 0 title '' ");
//    options.push_back("u 2:1 w l lc 'red' lw 8 lt 0 title '' ");
//    options.push_back("u 2:1 w l lc 'red' lw 8 lt 0 title '' ");

    for (Unsigned i=0; i<count; i++) options.push_back(op2);
//    options.push_back("u 2:1 w l lc 'blue' lw 8 lt 0 title '' ");
//    options.push_back("u 2:1 w l lc 'blue' lw 8 lt 0 title '' ");
//    options.push_back("u 2:1 w l lc 'blue' lw 8 lt 0 title '' ");
//    options.push_back("u 2:1 w l lc 'blue' lw 8 lt 0 title '' ");

    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }
    //pipe << "set key off";

    std::string title = species + ", flowering ranges";
    //---------------------------------------------------------------------------------
    //                          старт реальной отрисовки
    //---------------------------------------------------------------------------------
    pipe << gnuplot_packet("set title '" + title + "'");

    pipe << gnuplot_packet("set xrange" + xrange);
    pipe << "set xlabel 'T, days from the beginning of flowering'";

    pipe << "set yrange[-0.02:1.02]";
    pipe << "set ylabel 'h, rel.height'";
    pipe << "set ytics";

    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

//void plot_phase_range(const std::string &output, const std::string &species,
//    const std::vector<curve_t> &v1, const std::vector<curve_t> &v2, const std::string &xrange)
//{
//    if (v1.empty() || v2.empty()) return;

//    //---------------------------------------------------------------------------------------------
//    //                          параметры интерполяции данных
//    //---------------------------------------------------------------------------------------------
//    unsigned power = 3; // степень апроксимирующего полинома
//    double wy = 0.15; // ширина окна усреднения границ для кривых
//    double wy2 = 0.6; // ширина окна усреднения границ для 2й оси
//    double hy = 0.01; // расстояния между точками при разложении по полиномам

//    //---------------------------------------------------------------------------------------------
//    //                      пространственные области цветения
//    //---------------------------------------------------------------------------------------------
//    std::vector<curve_t> data1, data2;
//    for (auto e : v1) data1.push_back(e);
//    for (auto e : v2) data2.push_back(e);

//    Unsigned count = v1.size();
//    std::vector<curve_t> CLOSED1, CLOSED2;

//    for (Unsigned i=0; i<count; i++)
//    {
//        // построим правые границы цветения (старт + продолжительность)
//        curve_t data1R, data2R;
//        for (auto e : data1[i]) data1R.push_back({e[0], e[1] + e[2], sqrt(sqr(e[3]) + sqr(e[4]))});
//        for (auto e : data2[i]) data2R.push_back({e[0], e[1] + e[2], sqrt(sqr(e[3]) + sqr(e[4]))});

//        auto L1 = wslm_simple_smooth<CDIM>(0, power, data1[i],  0, 1, wy, hy);
//        auto R1 = wslm_simple_smooth<CDIM>(0, power, data1R, 0, 1, wy, hy);
//            // сглаженные правые и левые границы диапазона цветения 1й оси

//        auto xL1 = prg::subtract<CDIM>(prg::select<CDIM>(L1, 0, 1), prg::select<CDIM>(L1, 0, 2));
//        auto xR1 = prg::summarize<CDIM>(prg::select<CDIM>(R1, 0, 1), prg::select<CDIM>(R1, 0, 2));
//            // уширенные на стандартное отклонение правые и левые границы диапазона цветения 1й оси
//        auto range_closed1 = prg::closed_bound<CDIM>(xL1, xR1);
//            // полная область цветения 1й оси

//        auto BL2 = prg::lower_bound<CDIM>(data2[i].begin(), data2[i].end(), wy2, 0, 1);
//        auto BR2 = prg::upper_bound<CDIM>(data2R.begin(), data2R.end(), wy2, 0, 1);
//            // правые и левые границы диапазона цветения 2й оси

//        auto L2 = wslm_simple_smooth<CDIM>(0, power, BL2, 0, 1, wy2, hy);
//        auto R2 = wslm_simple_smooth<CDIM>(0, power, BR2, 0, 1, wy2, hy);

//            // сглаженные правые и левые границы диапазона цветения 2й оси
//        auto xL2 = prg::subtract<CDIM>(prg::select<CDIM>(L2, 0, 1), prg::select<CDIM>(L2, 0, 2));
//        auto xR2 = prg::summarize<CDIM>(prg::select<CDIM>(R2, 0, 1), prg::select<CDIM>(R2, 0, 2));
//            // уширенные на стандартное отклонение правые и левые границы диапазона цветения 2й оси
//        auto range_closed2 = prg::closed_bound<CDIM>(xL2, xR2);

//        //---------------------------------------------------------------------------------
//        //                  обращаем данные по высоте
//        //---------------------------------------------------------------------------------
//        for (auto &e : range_closed1) e[0] = 1 - e[0];
//        for (auto &e : range_closed2) e[0] = 1 - e[0];

//        CLOSED1.push_back(range_closed1);
//        CLOSED2.push_back(range_closed2);
//    }


//    //---------------------------------------------------------------------------------
//    //              рисуем диаграмму с помощью вызова gnuplot
//    //---------------------------------------------------------------------------------
//    std::swap(CLOSED1.front(), CLOSED1.back());
//    std::swap(CLOSED2.front(), CLOSED2.back());

//    std::vector<curve_t> patterns;
//    for (Unsigned i=0; i<count; i++) patterns.push_back(CLOSED1[i]);
//    curve_t LINE1;
//    LINE1.push_back(CLOSED1.back().front());
//    LINE1.push_back(CLOSED1.back().back());
//    patterns.push_back(LINE1);

//    for (Unsigned i=0; i<count; i++) patterns.push_back(CLOSED2[i]);
//    curve_t LINE2;
//    LINE2.push_back(CLOSED2.back().front());
//    LINE2.push_back(CLOSED2.back().back());
//    patterns.push_back(LINE2);

//    std::vector<gnuplot_packet> options;
//    options.push_back("u 2:1 w filledcu closed fc 'red' fill transparent solid 0.02 noborder title '\\{t,{/Symbol t}\\}_{m_1}(h)' ");
//    for (Unsigned i=1; i<count-1; i++) options.push_back("u 2:1 w filledcu closed fc 'red' fill transparent solid 0.02 noborder title ''");
//    options.push_back("u 2:1 w l lc 'red' lw 8 lt 0 title '\\{t,{/Symbol t}\\}_{e_1}(h)' ");
//    options.push_back("u 2:1 w l lc 'red' lw 8 lt 0 title '' ");

////    std::vector<gnuplot_packet> options;
////    options.push_back("u 2:1 w l lc 'red' title '\\{t,{/Symbol t}\\}_{m_1}(h)' ");
////    for (Unsigned i=1; i<count-1; i++) options.push_back("u 2:1 w l lc 'red' title ''");
////    options.push_back("u 2:1 w l lc 'red' lw 8 lt 0 title '\\{t,{/Symbol t}\\}_{e_1}(h)' ");
////    options.push_back("u 2:1 w l lc 'red' lw 8 lt 0 title '' ");

//    options.push_back("u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.02 noborder title '\\{t,{/Symbol t}\\}_{m_2}(h)'");
//    for (Unsigned i=1; i<count-1; i++) options.push_back("u 2:1 w filledcu closed fc 'blue' fill transparent solid 0.02 noborder title ''");
//    options.push_back("u 2:1 w l lc 'blue' lw 8 lt 0 title '\\{t,{/Symbol t}\\}_{e_2}(h)' ");
//    options.push_back("u 2:1 w l lc 'blue' lw 8 lt 0 title '' ");

//    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

//    gnuplot pipe; // инициализация gnuplot
//    if (!output.empty())
//    {
//        pipe << std::string("set output '" + output + "'");
//        pipe << "set term pdf enhanced";
//    }
//    //pipe << "set key off";

//    std::string title = species + ", flowering ranges";
//    //---------------------------------------------------------------------------------
//    //                          старт реальной отрисовки
//    //---------------------------------------------------------------------------------
//    pipe << gnuplot_packet("set title '" + title + "'");

//    pipe << gnuplot_packet("set xrange" + xrange);
//    pipe << "set xlabel 'T, days from the beginning of flowering'";

//    pipe << "set yrange[-0.02:1.02]";
//    pipe << "set ylabel 'h, rel.height'";
//    pipe << "set ytics";

//    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

//    if (!output.empty())
//    {
//        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
//        std::cout.flush();
//    }
//}

void plot_phase_curve(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2, const std::string &xrange)
{
    if (v1.empty() || v2.empty()) return;

    const int indent = 5; // необходимый отступ слева и справа
    //------------------------------------------------------------------------------------
    //                      выравнивание и определение границ
    //------------------------------------------------------------------------------------
    curve_t data1_(v1), data2_(v2);

    double ts = 1000, te = -1000;
    for (auto e : data1_)
    {
        if (ts > e[1]) ts = e[1];
        if (te < e[1] + e[2]) te = e[1] + e[2];
    }
    for (auto e : data2_)
    {
        if (ts > e[1]) ts = e[1];
        if (te < e[1] + e[2]) te = e[1] + e[2];
    }
    int tstart = int(std::floor(ts));
    int tfinish = int(std::ceil(te));

    //------------------------------------------------------------------------------------
    //                      накопление кривых по ансамблю
    //------------------------------------------------------------------------------------
    unsigned ndays = unsigned(tfinish - tstart + 1);
    std::vector<unsigned> sum1(ndays, 0);
    std::vector<unsigned> sum2(ndays, 0);

    for (auto e : data1_)
    {
        int ta = int(std::round(e[1]));
        int tb = int(std::round(e[1] + e[2]));
        for (int i=ta; i<=tb; i++)
        {
            unsigned ndx = unsigned(i - tstart);
            assert(ndx < ndays);
            sum1[ndx]++;
        }
    }
    for (auto e : data2_)
    {
        int ta = int(std::round(e[1]));
        int tb = int(std::round(e[1] + e[2]));
        for (int i=ta; i<=tb; i++)
        {
            unsigned ndx = unsigned(i - tstart);
            assert(ndx < ndays);
            sum2[ndx]++;
        }
    }

    //------------------------------------------------------------------------------------
    //              выбор данных по интервалу и нормировка
    //------------------------------------------------------------------------------------
    Unsigned num_flowers = data1_.size() + data2_.size();
    curve_t data1, data2, datasum, datadif;

    // добавление indent, чтобы уйти от резких разрывов до и после графика
    for (int t=tstart-indent; t<tstart; t++)
    {
        data1.push_back({double(t-tstart), 0., 0.});
        data2.push_back({double(t-tstart), 0., 0.});
        datasum.push_back({double(t-tstart), 0., 0.});
        datadif.push_back({double(t-tstart), 0., 0.});
    }
    for (int t=tstart; t<=tfinish; t++)
    {
        data1.push_back({double(t-tstart), sum1[unsigned(t-tstart)] / double(num_flowers), 0});
        data2.push_back({double(t-tstart), sum2[unsigned(t-tstart)] / double(num_flowers), 0});
        datasum.push_back({double(t-tstart), (sum1[unsigned(t-tstart)] + sum2[unsigned(t-tstart)]) / double(num_flowers), 0});
        datadif.push_back({double(t-tstart), int(sum1[unsigned(t-tstart)] - sum2[unsigned(t-tstart)]) / double(num_flowers), 0});
    }
    for (int t=tfinish+1; t<tfinish+indent; t++)
    {
        data1.push_back({double(t-tstart), 0., 0.});
        data2.push_back({double(t-tstart), 0., 0.});
        datasum.push_back({double(t-tstart), 0., 0.});
        datadif.push_back({double(t-tstart), 0., 0.});
    }

    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    unsigned power = 3; // степень апроксимирующего полинома
    //double wx = 0.75; // ширина окна усреднения
    double wx = 1.25; // ширина окна усреднения
    double hx = 0.125; // расстояния между точками при разложении по полиномам

    auto data_curve1(data1), data_curve2(data2), data_curvesum(datasum), data_curvedif(datadif);
    for (auto &e : data_curve1) e[1] *= 100;
    for (auto &e : data_curve2) e[1] *= 100;
    for (auto &e : data_curvesum) e[1] *= 100;
    for (auto &e : data_curvedif) e[1] *= 100;
        // перевод в %

    auto curve1 = wslm_simple_smooth<CDIM>(0, power, data_curve1, 0, 1, wx, hx);
    auto curve2 = wslm_simple_smooth<CDIM>(0, power, data_curve2, 0, 1, wx, hx);
    auto curvesum = wslm_simple_smooth<CDIM>(0, power, data_curvesum, 0, 1, wx, hx);
    auto curvedif = wslm_simple_smooth<CDIM>(0, power, data_curvedif, 0, 1, wx, hx);
        // кривые расцветания на левой и правой границах, соответственно

    auto curve_closed1 = prg::closed_bound<CDIM>(curve1);
    auto curve_closed2 = prg::closed_bound<CDIM>(curve2);
    auto curve_closedsum = prg::closed_bound<CDIM>(curvesum);
    auto curve_closeddif = prg::closed_bound<CDIM>(curvedif);
        // закрытые области, описывающие границы и их погрешности

    curve_t zero;
    zero.push_back({ts, 0., 0.});
    zero.push_back({te, 0., 0.});

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------
    std::vector<curve_t> patterns
    {
    #ifdef DRAW_ALL_CURVES
        zero,
        //curve_closedsum, curve_closeddif,
        curve_closed1, curve_closed2,
        curve1, curve2, curvesum, curvedif
    #else
        curve_closed1, curve_closed2
    #endif
    };

    std::vector<gnuplot_packet> options
    {
    #ifdef DRAW_ALL_CURVES
        "u 1:2 w lines lc 'black' lt 0 title ''",
//        "u 1:2 w filledcu closed fc 'gray' fill transparent solid 0.9 noborder title '{/Symbol d}f(t)' ",
//        "u 1:2 w filledcu closed fc 'gray' fill transparent solid 0.9 noborder title '' ",
        "u 1:2 w filledcu closed fc 'red' fill transparent solid 0.4 noborder title '{/Symbol d}f_1(t)' ",
        "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.4 noborder title '{/Symbol d}f_2(t)' ",
        "u 1:2 w lines lc 'red' lw 3 title 'f_1(t)'",
        "u 1:2 w lines lc 'blue' lw 3 title 'f_2(t)'",
        "u 1:2 w lines lc 'black' lw 5 title 'f(t)'",
        "u 1:2 w lines lt 0 lc 'black' lw 5 title '{/Symbol D}f(t)'"
    #else
        "u 1:2 w filledcu closed fc 'red' fill transparent solid 0.2 noborder ",
        "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder "
    #endif
    };
    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }

    //pipe << "set key off";

    std::string title = species + ", flowering curves";
    //---------------------------------------------------------------------------------
    //                          старт реальной отрисовки
    //---------------------------------------------------------------------------------
    pipe << gnuplot_packet("set title '" + title);

    pipe << gnuplot_packet("set xrange" + xrange);
    pipe << "set xlabel 'T, days from the beginning of flowering'";

    pipe << "set yrange[*:*]";
    pipe << "set ylabel '%, open flowers'\n";

    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

void plot_phase_curve(const std::string &output, const std::string &species,
    const std::vector<curve_t> &v1, const std::vector<curve_t> &v2, const std::string &xrange)
{
    if (v1.empty() || v2.empty()) return;

    const int indent = 5; // необходимый отступ слева и справа
    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    unsigned power = 3; // степень апроксимирующего полинома
    double wx = 1.25; // ширина окна усреднения
    double hx = 0.125; // расстояния между точками при разложении по полиномам

    //------------------------------------------------------------------------------------
    //                      выравнивание и определение границ
    //------------------------------------------------------------------------------------

    Unsigned count = v1.size(); // число спектров
    std::vector<curve_t> data1_(v1), data2_(v2);

    double ts = 1000, te = -1000;
    for (Unsigned i=0; i<count; i++)
    {
        for (auto e : data1_[i])
        {
            if (ts > e[1]) ts = e[1];
            if (te < e[1] + e[2]) te = e[1] + e[2];
        }
        for (auto e : data2_[i])
        {
            if (ts > e[1]) ts = e[1];
            if (te < e[1] + e[2]) te = e[1] + e[2];
        }
    }
    int tstart = int(std::floor(ts));
    int tfinish = int(std::ceil(te));

    //------------------------------------------------------------------------------------
    //                      накопление кривых по ансамблю
    //------------------------------------------------------------------------------------
    unsigned ndays = unsigned(tfinish - tstart + 1);

    std::vector<curve_t> SUM;
    Unsigned ss = 0;

    for (Unsigned i=0; i<count; i++)
    {

        std::vector<unsigned> sum1(ndays, 0);
        std::vector<unsigned> sum2(ndays, 0);

        for (auto e : data1_[i])
        {
            int ta = int(std::round(e[1]));
            int tb = int(std::round(e[1] + e[2]));
            for (int i=ta; i<=tb; i++)
            {
                unsigned ndx = unsigned(i - tstart);
                assert(ndx < ndays);
                sum1[ndx]++;
                ss++;
            }
        }
        for (auto e : data2_[i])
        {
            int ta = int(std::round(e[1]));
            int tb = int(std::round(e[1] + e[2]));
            for (int i=ta; i<=tb; i++)
            {
                unsigned ndx = unsigned(i - tstart);
                assert(ndx < ndays);
                sum2[ndx]++;
                ss++;
            }
        }

        //------------------------------------------------------------------------------------
        //              выбор данных по интервалу и нормировка
        //------------------------------------------------------------------------------------
        Unsigned num_flowers = data1_[i].size() + data2_[i].size();
        curve_t data1, data2, datasum, datadif;

        // добавление indent, чтобы уйти от резких разрывов до и после графика
        for (int t=tstart-indent; t<tstart; t++)
        {
            data1.push_back({double(t-tstart), 0., 0.});
            data2.push_back({double(t-tstart), 0., 0.});
            datasum.push_back({double(t-tstart), 0., 0.});
            datadif.push_back({double(t-tstart), 0., 0.});
        }
        for (int t=tstart; t<=tfinish; t++)
        {
            data1.push_back({double(t-tstart), sum1[unsigned(t-tstart)] / double(num_flowers), 0});
            data2.push_back({double(t-tstart), sum2[unsigned(t-tstart)] / double(num_flowers), 0});
            datasum.push_back({double(t-tstart), (sum1[unsigned(t-tstart)] + sum2[unsigned(t-tstart)]) / double(num_flowers), 0});
            datadif.push_back({double(t-tstart), int(sum1[unsigned(t-tstart)] - sum2[unsigned(t-tstart)]) / double(num_flowers), 0});
        }
        for (int t=tfinish+1; t<tfinish+indent; t++)
        {
            data1.push_back({double(t-tstart), 0., 0.});
            data2.push_back({double(t-tstart), 0., 0.});
            datasum.push_back({double(t-tstart), 0., 0.});
            datadif.push_back({double(t-tstart), 0., 0.});
        }

        auto data_curve1(data1), data_curve2(data2), data_curvesum(datasum), data_curvedif(datadif);
        for (auto &e : data_curve1) e[1] *= 100;
        for (auto &e : data_curve2) e[1] *= 100;
        for (auto &e : data_curvesum) e[1] *= 100;
        for (auto &e : data_curvedif) e[1] *= 100;
            // перевод в %

        auto curve1 = wslm_simple_smooth<CDIM>(0, power, data_curve1, 0, 1, wx, hx);
        auto curve2 = wslm_simple_smooth<CDIM>(0, power, data_curve2, 0, 1, wx, hx);
        auto curvesum = wslm_simple_smooth<CDIM>(0, power, data_curvesum, 0, 1, wx, hx);
        auto curvedif = wslm_simple_smooth<CDIM>(0, power, data_curvedif, 0, 1, wx, hx);
            // кривые расцветания на левой и правой границах, соответственно

        auto curve_closed1 = prg::closed_bound<CDIM>(curve1);
        auto curve_closed2 = prg::closed_bound<CDIM>(curve2);
        auto curve_closedsum = prg::closed_bound<CDIM>(curvesum);
        auto curve_closeddif = prg::closed_bound<CDIM>(curvedif);
            // закрытые области, описывающие границы и их погрешности

        SUM.push_back(curvesum);
    }

    Unsigned sz = SUM[0].size();
    curve_t AVERAGE(sz, {0, 0});

    for (Unsigned k=0; k<sz; k++)
    {
        AVERAGE[k][0] = SUM[0][k][0];
        for (Unsigned i=0; i<count; i++) AVERAGE[k][1] += SUM[i][k][1];
        AVERAGE[k][1] /= count;
    }

    curve_t zero;
    zero.push_back({ts, 0., 0.});
    zero.push_back({te, 0., 0.});

    //---------------------------------------------------------------------------------
    //              рисуем диаграмму с помощью вызова gnuplot
    //---------------------------------------------------------------------------------

    std::vector<curve_t> patterns;
    for (auto e : SUM) patterns.push_back(e);
    patterns.push_back(AVERAGE);
        // ставим учетный на последнее место, чтобы не затирать при отрисовке

    std::vector<gnuplot_packet> options;
    for (Unsigned i=0; i<count-1; i++) options.push_back("u 1:2 w l lc 'gray' lw 1 title ''");
    options.push_back("u 1:2 w l lc 'gray' lw 1 title 'f_m(t)'");
    options.push_back("u 1:2 w l lt 0 lc 'black' lw 8 title '<f>(t)'");

    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

    gnuplot pipe; // инициализация gnuplot
    if (!output.empty())
    {
        pipe << std::string("set output '" + output + "'");
        pipe << "set term pdf enhanced";
    }

    std::string title = species + ", flowering curves";
    //---------------------------------------------------------------------------------
    //                          старт реальной отрисовки
    //---------------------------------------------------------------------------------
    pipe << gnuplot_packet("set title '" + title);

    pipe << gnuplot_packet("set xrange" + xrange);
    pipe << "set xlabel 'T, days from the beginning of flowering'";

    pipe << "set yrange[*:*]";
    pipe << "set ylabel '%, open flowers'\n";

    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

    if (!output.empty())
    {
        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
        std::cout.flush();
    }
}

////////////////////////////////////////////////////////////////////////////////
//void plot_spanPDF(const std::string &output, const std::string &species,
//    const curve_t &v1, const curve_t &v2)
//{
//    //---------------------------------------------------------------------------------------------
//    //                          параметры интерполяции данных
//    //---------------------------------------------------------------------------------------------
//    unsigned power = 3; // степень апроксимирующего полинома
//    double wx = 1.5; // ширина окна усреднения
//    double hx = 0.25; // расстояния между точками при разложении по полиномам

//    //---------------------------------------------------------------------------------------------
//    //                      расчет гистограммы цветения
//    //---------------------------------------------------------------------------------------------
//    unsigned num_flowers = unsigned(v1.size() + v2.size());
//    auto flowerspan1 = construct_flowerspan(v1, num_flowers);
//    auto flowerspan2 = construct_flowerspan(v2, num_flowers);

//    //---------------------------------------------------------------------------------------------
//    //                  сглаживание гистограммы цветения
//    //---------------------------------------------------------------------------------------------
//    auto fcurve1 = wslm_simple_smooth<CDIM>(0, power, flowerspan1, 0, 1, wx, hx);
//    auto fcurve2 = wslm_simple_smooth<CDIM>(0, power, flowerspan2, 0, 1, wx, hx);

//    for (auto &e : fcurve1) if (e[1] < 0) e[1] = 0;
//    for (auto &e : fcurve2) if (e[1] < 0) e[1] = 0;

//    auto closed1 = prg::closed_bound<CDIM>(fcurve1);
//    auto closed2 = prg::closed_bound<CDIM>(fcurve2);

//    //---------------------------------------------------------------------------------
//    //              рисуем диаграмму с помощью вызова gnuplot
//    //---------------------------------------------------------------------------------
//    std::vector<curve_t> patterns
//    {
//        closed1, closed2, flowerspan1, flowerspan2, flowerspan1, flowerspan2
//    };

//    std::vector<gnuplot_packet> options
//    {
//        "u 1:2 w filledcu closed fc 'red' fill transparent solid 0.2 noborder ",
//        "u 1:2 w filledcu closed fc 'blue' fill transparent solid 0.2 noborder ",
//        "u 1:2 w l lt 0 fc 'red' ",
//        "u 1:2 w l lt 0 fc 'blue' ",
//        "u 1:2 w p lt 7 fc 'red' ",
//        "u 1:2 w p lt 7 fc 'blue' "
//    };

//    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

//    gnuplot pipe; // инициализация gnuplot
//    if (!output.empty())
//    {
//        pipe << std::string("set output '" + output + "'");
//        pipe << "set term pdf enhanced";
//    }

//    //---------------------------------------------------------------------------------
//    //                          старт реальной отрисовки
//    //---------------------------------------------------------------------------------
//    pipe << "set key off";

//    std::string title = species + ", flowering lifespan histogram";
//    pipe << gnuplot_packet("set title '" + title);

//    pipe << "set xlabel '{/Symbol D}T, flowering time '";
//    pipe << "set xrange[-0.5:*]";

//    pipe << "set yrange[0:*]";
//    pipe << "set ylabel '%, flowers'\n";
//    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

//    if (!output.empty())
//    {
//        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
//        std::cout.flush();
//    }
//}

//void plot_countPDF(const std::string &output, const std::string &species,
//    const curve_t &v1, const curve_t &v2, const curve_t &vs)
//{
//    //---------------------------------------------------------------------------------------------
//    //                          рассчитаем функцию распределения
//    //---------------------------------------------------------------------------------------------
//    const unsigned nmid = 8; // число промежуточных точек
//    double indent1 = std::max(0.5 * (v1.back()[0] - v1.front()[0]), double(v1.front()[0]));
//    double indent2 = std::max(0.5 * (v2.back()[0] - v2.front()[0]), double(v2.front()[0]));
//    double indents = std::max(0.5 * (vs.back()[0] - vs.front()[0]), double(vs.front()[0]));

//    curve_t data1, data2, datas;
//    double prevx = v1.front()[0] - indent1, cnt = 0, dx = 0;
//    for (unsigned i=0; i<v1.size(); i++)
//    {
//        dx = (v1[i][0] - prevx) / nmid;
//        for (unsigned k=0; k<nmid; k++)
//        {
//            data1.push_back({prevx, cnt});
//            prevx += dx;
//        }
//        cnt++;
//        prevx = v1[i][0];
//    }
//    dx = indent1 / nmid;
//    for (unsigned k=0; k<nmid; k++)
//    {
//        data1.push_back({prevx, cnt});
//        prevx += dx;
//    }
//    double count1 = cnt;

//    prevx = v2.front()[0] - indent2; cnt = 0; dx = 0;
//    for (unsigned i=0; i<v2.size(); i++)
//    {
//        dx = (v2[i][0] - prevx) / nmid;
//        for (unsigned k=0; k<nmid; k++)
//        {
//            data2.push_back({prevx, cnt});
//            prevx += dx;
//        }
//        cnt++;
//        prevx = v2[i][0];
//    }

//    dx = indent2 / nmid;
//    for (unsigned k=0; k<nmid; k++)
//    {
//        data2.push_back({prevx, cnt});
//        prevx += dx;
//    }
//    double count2 = cnt;

//    prevx = vs.front()[0] - indents; cnt = 0; dx = 0;
//    for (unsigned i=0; i<vs.size(); i++)
//    {
//        dx = (vs[i][0] - prevx) / nmid;
//        for (unsigned k=0; k<nmid; k++)
//        {
//            datas.push_back({prevx, cnt});
//            prevx += dx;
//        }
//        cnt++;
//        prevx = vs[i][0];
//    }

//    dx = indents / nmid;
//    for (unsigned k=0; k<nmid; k++)
//    {
//        datas.push_back({prevx, cnt});
//        prevx += dx;
//    }
//    //---------------------------------------------------------------------------------------------
//    //                          параметры интерполяции данных
//    //---------------------------------------------------------------------------------------------
//    unsigned power = 2; // степень апроксимирующего полинома
//    double hx = 0.125; // расстояния между точками при разложении по полиномам
//    double wx1 = 0.125 * (data1.back()[0] - data1.front()[0]); // ширина окна усреднения
//    double wx2 = 0.125 * (data2.back()[0] - data2.front()[0]); // ширина окна усреднения
//    double wxs = 0.125 * (datas.back()[0] - datas.front()[0]); // ширина окна усреднения

//    double norm = 0;
//    for (auto e : data1) if (norm < e[1]) norm = e[1];
//    norm = 1./norm;
//    for (auto &e : data1) { e[1] *= norm; e[2] *= norm; }

//    norm = 0;
//    for (auto e : data2) if (norm < e[1]) norm = e[1];
//    norm = 1./norm;
//    for (auto &e : data2) { e[1] *= norm; e[2] *= norm; }

//    norm = 0;
//    for (auto e : datas) if (norm < e[1]) norm = e[1];
//    norm = 1./norm;
//    for (auto &e : datas) { e[1] *= norm; e[2] *= norm; }
//    //---------------------------------------------------------------------------------------------
//    //                      расчет гистограммы цветения
//    //---------------------------------------------------------------------------------------------
//    auto f1d = wslm_simple_smooth<CDIM>(1, power, data1, 0, 1, wx1, hx);
//    auto f2d = wslm_simple_smooth<CDIM>(1, power, data2, 0, 1, wx2, hx);
//    auto f2s = wslm_simple_smooth<CDIM>(1, power, datas, 0, 1, wxs, hx);

//    //---------------------------------------------------------------------------------------------
//    //                  нормировка на полное число цветков на площадь
//    //---------------------------------------------------------------------------------------------
//    norm = 0; for (auto e : f1d) norm += e[1]; norm = count1 / norm;
//    for (auto &e : f1d) { e[1] *= norm; e[2] *= norm; }
//    norm = 0; for (auto e : f2d) norm += e[1]; norm = count2 / norm;
//    for (auto &e : f2d) { e[1] *= norm; e[2] *= norm; }
//    norm = 0; for (auto e : f2s) norm += e[1]; norm = (count1 + count2) / norm;
//    for (auto &e : f2s) { e[1] *= norm; e[2] *= norm; }
//        // нормировали так, что полный интеграл равен числу цветков

//    norm = 0;
//    for (auto e : f1d) if (norm < e[1]) norm = e[1];
//    for (auto e : f2d) if (norm < e[1]) norm = e[1];
//    for (auto e : f2s) if (norm < e[1]) norm = e[1];
//    norm = 1./norm;
//    for (auto &e : f1d) { e[1] *= norm; e[2] *= norm; }
//    for (auto &e : f2d) { e[1] *= norm; e[2] *= norm; }
//    for (auto &e : f2s) { e[1] *= norm; e[2] *= norm; }
//        // нормировали так, что максимум одного из двух графиков равен 1

//    auto v1_ = v1, v2_ = v2;
//    norm = 0;
//    for (auto e : v1) if (norm < e[1]) norm = e[1];
//    for (auto e : v2) if (norm < e[1]) norm = e[1];
//    norm = .5 / norm;
//    for (auto &e : v1_) e[1] *= norm;
//    for (auto &e : v2_) e[1] *= norm;
//        // нормировали так, что максимум одного из двух графиков равен 0.75

//    //---------------------------------------------------------------------------------
//    //              рисуем диаграмму с помощью вызова gnuplot
//    //---------------------------------------------------------------------------------
//    std::vector<curve_t> patterns
//    {
//        v1_, v2_, f1d, f2d, f2s
//    };

//    std::vector<gnuplot_packet> options
//    {
//        "u 1:(0):(0):2 w vectors filled nohead lw 4 lt 0 lc 'red' title ''",
//        "u 1:(0):(0):2 w vectors filled nohead lw 4 lt 0 lc 'blue' title ''",
//        "u 1:2 w l lt 7 fc 'red' title '{/Symbol Y}_1(n)'",
//        "u 1:2 w l lt 7 fc 'blue' title '{/Symbol Y}_2(n)'",
//        "u 1:2 w l lt 0 lw 8 fc 'black' title '{/Symbol Y}(n)'"
//    };

//    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

//    gnuplot pipe; // инициализация gnuplot
//    if (!output.empty())
//    {
//        pipe << std::string("set output '" + output + "'");
//        pipe << "set term pdf enhanced";
//    }

//    //---------------------------------------------------------------------------------
//    //                          старт реальной отрисовки
//    //---------------------------------------------------------------------------------
//    //pipe << "set key off";

//    std::string title = species + ", PDFs ";
//    pipe << gnuplot_packet("set title '" + title + "'");

//    pipe << "set xlabel 'n, number of flowers'";
//    pipe << "set xrange[0:*]";

//    pipe << "set yrange[0:1.05]";
//    pipe << "set ylabel '%, flowers'\n";
//    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, patterns.size());

//    if (!output.empty())
//    {
//        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
//        std::cout.flush();
//    }
//}

//void plot_densityRF(const std::string &output, const std::string &species,
//    const curve_t &data1, const curve_t &data2)
//{
//    if (data1.empty() || data2.empty()) return;

//    //---------------------------------------------------------------------------------------------
//    //                          параметры интерполяции данных
//    //---------------------------------------------------------------------------------------------
//    unsigned power = 2; // степень апроксимирующего полинома
//    double hy = 0.0125; // расстояния между точками при разложении по полиномам
//    double wyd = .15; // ширина окна усреднения для плотности

//    //---------------------------------------------------------------------------------------------
//    //                          расчет графика плотности
//    //---------------------------------------------------------------------------------------------
//    curve_t density1_data, density2_data;
//    double max = 0;
//    {
//        typedef std::pair<double, unsigned> key_t;

//        std::map<key_t, unsigned> sum;
//        for (auto e : data1)
//        {
//            key_t key{e[0], e[5]};
//            if (sum.count(key) == 0) sum[key] = 0;
//            sum[key]++;
//        }
//		for (auto e : sum) density1_data.push_back({e.first.first, double(e.second)});

//        sum.clear();
//        for (auto e : data2)
//        {
//            key_t key{e[0], e[5]};
//            if (sum.count(key) == 0) sum[key] = 0;
//            sum[key]++;
//        }
//		for (auto e : sum) density2_data.push_back({e.first.first, double(e.second)});

//        for (auto e : sum) if (max < e.second) max = e.second;
//    }

//    auto density1_curve = wslm_simple_smooth<CDIM>(0, power, density1_data, 0, 1, wyd, hy);
//    auto density1_closed = prg::closed_bound<CDIM>(density1_curve);

//    auto density2_curve = wslm_simple_smooth<CDIM>(0, power, density2_data, 0, 1, wyd, hy);
//    auto density2_closed = prg::closed_bound<CDIM>(density2_curve);
//        // данные, кривая и закрытая область для описания плотности цветения 2й оси

//    //---------------------------------------------------------------------------------
//    //              "облегчим" графики от части сырых данных
//    //---------------------------------------------------------------------------------
//    std::uniform_int_distribution<> d;
//    typedef std::uniform_int_distribution<>::param_type param_type;
//    curve_t data1_, data2_, density2_data_;

//    if (DRAW_ONLY < data1.size())
//    {
//        d.param(param_type{0, int(data1.size())-1});
//        for (Unsigned i=0; i<DRAW_ONLY; i++) data1_.push_back(data1[Unsigned(d(rangen))]);
//        std::sort(data1_.begin(), data1_.end(),
//            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
//        );
//    }
//    else for (auto e : data1) data1_.push_back(e);

//    if (DRAW_ONLY < density2_data.size())
//    {
//        d.param(param_type{0, int(density2_data.size())-1});
//        for (Unsigned i=0; i<DRAW_ONLY; i++) density2_data_.push_back(density2_data[Unsigned(d(rangen))]);
//        std::sort(density2_data_.begin(), density2_data_.end(),
//            [](const point<CDIM> &a, const point<CDIM> &b){ return a[0] < b[0]; }
//        );
//    }
//    else for (auto e : density2_data) density2_data_.push_back(e);

//    //---------------------------------------------------------------------------------
//    //                  обращаем данные по высоте
//    //---------------------------------------------------------------------------------
//    for (auto &e : data1_) e[0] = 1 - e[0];
//    for (auto &e : data2_) e[0] = 1 - e[0];
//    for (auto &e : density1_curve) e[0] = 1 - e[0];
//    for (auto &e : density2_curve) e[0] = 1 - e[0];
//    for (auto &e : density2_data_) e[0] = 1 - e[0];

//    //---------------------------------------------------------------------------------
//    //              рисуем диаграмму с помощью вызова gnuplot
//    //---------------------------------------------------------------------------------
//    std::vector<curve_t> patterns
//    {
//        density1_curve, density2_curve,
//        density2_data_
//    };

//    std::vector<gnuplot_packet> options
//    {
//        "u 2:1 w filledcu y1 fc 'red' fill transparent solid 0.2 noborder title 'f_{{/Symbol r}_1}(h)' ",
//        "u 2:1 w filledcu y1 fc 'blue' fill transparent solid 0.2 noborder title 'f_{{/Symbol r}_2}(h)' ",
//        "u 2:1:3 w xerr pt 7 fc 'blue' title ''",
//    };

//    gnuplot_prepare<CDIM>(patterns); // подготовка временных файлов данных для вызова gnuplot

//    gnuplot pipe; // инициализация gnuplot
//    if (!output.empty())
//    {
//        pipe << std::string("set output '" + output + "'");
//        pipe << "set term pdf enhanced";
//    }


//    std::string title = species;
//    pipe << gnuplot_packet("set title '" + title + "'");

//    pipe << "set ylabel 'h, rel.height of flower'\n";
//    //pipe << "set key off";
//    pipe << "set yrange[-0.05:1.05]";

//    //pipe << "set xrange[0.0:*]";
//    gnuplot_packet pack("set xrange[0:");
//    pack << max + 1 << "]";
//    pipe << pack;

//    pipe << "set xlabel '{/Symbol s}, flowering density'";
//    pipe << "set xtics ('0' 0, '1' 1, '2' 2, '3' 3, '4' 4, '5' 5, '6' 6, '7' 7, '8' 8)";

//    gnuplot_insert_plot<CDIM>(pipe, patterns, options, 0, 3);

//    if (!output.empty())
//    {
//        printf("gnuplot_diagram '%s' is drawn\n", output.c_str());
//        std::cout.flush();
//    }
//}

