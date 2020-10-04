#ifndef ORN_GNUPLOT_H
#define ORN_GNUPLOT_H

#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include "orn/orn.defs.h"
#include "math/math.h"
#include "gnu/gnuplot.h"

/// при установке рисуются все кривые (1я ось, 2я ось, сумма и разность)
#define DRAW_ALL_CURVES

/// вывести только заданное число точек, чтобы не загромождать график
#define DRAW_ONLY 30

namespace prg
{
    struct gnuplotParams
    {
    private:
        std::vector<shared<curve<4>>> line_data;
        std::vector<gnuplot_packet> line_options;

        std::string xrange = "[*:*]"; ///< область оси
        std::string yrange = "[*:*]"; ///< область оси

    protected:

        std::string range(double a, double b)
            { return "[" + std::to_string(a) + ":" + std::to_string(b) + "]"; }

    public:

        std::string output = ""; ///< файл вывода (или ""[terminal])
        std::string title = ""; ///< заголовок графика
        std::string key = "off"; ///< нет вывода легенды (по умалчиванию)
        std::string xlabel = ""; ///< метка оси
        std::string ylabel = ""; ///< метка оси
        std::string ytics = ""; ///< метки оси

        void clear()
        {
            line_data.clear();
            line_options.clear();
            set_xrange();
            set_yrange();
        }

        void set_xrange(double a, double b) { xrange = range(a, b); }
        void set_yrange(double a, double b) { yrange = range(a, b); }
        void set_xrange() { xrange = "[*:*]"; }
        void set_yrange() { yrange = "[*:*]"; }

        /// вставить график (данные и опции рисования)
        void insert_data(shared<curve<2>> p,
            const gnuplot_packet &option="u 1:2 w l lt 7 lw 1 fc 'black' title '' ")
        {
            shared<curve<4>> p4(new curve<4>);
            for (auto e : *p) p4->push_back({e[0], e[1], 0, 0});
            line_data.push_back(p4);
            line_options.push_back(option);
        }

        void insert_data(shared<curve<4>> p,
            const gnuplot_packet &option="u 1:2 w l lt 7 lw 1 fc 'black' title '' ")
            { line_data.push_back(p); line_options.push_back(option); }

        /// число графиков
        Unsigned size() const { return line_data.size(); }
        bool empty() const { return line_data.empty(); }

        /// получить данные графика
        shared<curve<4>> get_data(Unsigned i) const { return line_data[i]; }
        gnuplot_packet get_options(Unsigned i) const { return line_options[i]; }

        const std::vector<shared<curve<4>>> &get_data() const { return line_data; }
        const std::vector<gnuplot_packet> &get_options() const { return line_options; }

        std::string get_xrange() const { return xrange; }
        std::string get_yrange() const { return yrange; }
    };

    inline void plot(const gnuplotParams &param)
    {
        // if (param.empty()) return;
        // закоментировано, чтобы рисовать пустой график

        //---------------------------------------------------------------------------------
        //              рисуем диаграмму с помощью вызова gnuplot
        //---------------------------------------------------------------------------------
        // подготовка временных файлов данных для вызова gnuplot
        for (Unsigned i=0; i<param.size(); i++)
            gnuplot_prepare(i, param.get_data(i));

        gnuplot pipe; // инициализация gnuplot
        if (!param.output.empty())
        {
            pipe << std::string("set output '" + param.output + "'");
            pipe << "set term pdf enhanced";
        }

        //---------------------------------------------------------------------------------
        //                          старт реальной отрисовки
        //---------------------------------------------------------------------------------
        pipe << gnuplot_packet("set key " + param.key);
        pipe << gnuplot_packet("set title '" + param.title + "'");

        pipe << gnuplot_packet("set xrange" + param.get_xrange());
        pipe << gnuplot_packet("set xlabel '" + param.xlabel + "'");

        pipe << gnuplot_packet("set yrange" + param.get_yrange());
        pipe << gnuplot_packet("set ylabel '" + param.ylabel + "'");

        gnuplot_insert_plot(pipe, param.get_data(), param.get_options());

        if (param.empty())
        {
            pipe << std::string("plot 0 w l lc 'black' title '' ");
        }
    }

} // namespace
#endif
