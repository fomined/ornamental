#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <random>

#include "prg.os.h"
#include "prg.random.h"

#include "test.defs.h"
#include "math.defs.h"
#include "math.WLSM.h"
#include "math.statistics.h"
#include "flower.generator.h"

void Flowermodel::insert_flower(const Flowerstamp &fs, const point<3> &v)
{
    flowers.push_back(fs);
    phases.push_back(v);

    unsigned axis = fs.axis_order();
    if (axis == 1) num_flowers1++;
    if (axis == 2) num_flowers2++;
}

curve_t Flowermodel::flowering_select(unsigned axis, unsigned char flags) const
{
    //------------------------------------------------------------------------------------
    //                  получение стартов для выравнивания
    //------------------------------------------------------------------------------------
    const double inf = 1000.;

    double sstart = inf; // старт побега для выравнивания
    std::vector<double> hstarts(num_flowers1+1, inf); // старты на разных высотах

    for (Unsigned i=0; i<num_flowers(); i++)
    {
        Flowerstamp fs = flowers[i];
        auto time = phases[i];
        if (fs.axis_order() != 1) continue; // старты задаются только по 1й оси

        double t = time[1]; // время старта цветка 1й оси
        unsigned j = fs.get(1); // высота цветка 1й оси

        if (sstart > t) sstart = t;
        if (hstarts[j] > t) hstarts[j] = t;
    }

    //------------------------------------------------------------------------------------
    //                  получение данных для выбранных цветков
    //------------------------------------------------------------------------------------
    curve_t tmp;

    for (Unsigned i=0; i<num_flowers(); i++)
    {
        Flowerstamp fs = flowers[i];
        auto time = phases[i];
        if (fs.axis_order() != axis) continue;

        double t = time[1]; // время старта цветка 1й оси
        double te = time[2]; // время отцветания цветка
        unsigned j = fs.get(1); // высота цветка

        if (flags & ALIGN_SPROUT) // выравнивание по побегу
        { t -= sstart; te -= sstart; }

        if ((flags & ALIGN_1AXIS) && axis == 2) // выравнивание по 1й оси
        { t -= hstarts[j]; te -= hstarts[j]; }

        double height = j;
        if (flags & NORM) // исправление высот
        {
            double norm = 1. / (num_flowers1 - 1);
            height = (j - 1) * norm;
        }

        tmp.push_back({height, t, te-t, 0., 0.});
    }

    if (flags & DAY_ALIGN_CODE) // выровнять границы диапазона (колонки 1, 2)
        tmp = prg::to_int<CDIM>(tmp, {1, 2});

    return tmp;
}

pFlowermodel Flowergenerator::generate_sprout(unsigned sprout, int num1, int num2)
{
    assert(!params.empty());

    std::normal_distribution<> flowernum_distribution;
    std::normal_distribution<> budspan_distribution;
    std::normal_distribution<> start_distribution;
    std::normal_distribution<> span_distribution;
    std::normal_distribution<> density_distribution;
    std::normal_distribution<> between_distribution;

    pFlowermodel pmodel(new Flowermodel());

    //--------------------------------------------------------------------------------------------
    //                          генерация числа цветков
    //--------------------------------------------------------------------------------------------
    unsigned num_flowers1 = unsigned(num1);
    while (int(num_flowers1) < 0)
    {
        flowernum_distribution.param(params.flowernum_params(1));
        num_flowers1 = unsigned(std::round(flowernum_distribution(rangen)));
    }

    unsigned num_flowers2 = unsigned(num2);
    while (int(num_flowers2) < 0)
    {
        flowernum_distribution.param(params.flowernum_params(2));
        num_flowers2 = unsigned(std::round(num_flowers1 * flowernum_distribution(rangen)));
            // так как в params.flowernum_params(2) хранится доля от числа num1
    }
    budspan_distribution.param(params.budspan_params());

    //--------------------------------------------------------------------------------------------
    //                          генерация и нормировка плотности
    //--------------------------------------------------------------------------------------------
    std::vector<double> density(num_flowers1 + 1, 0.);

    double sumd = 0;
    for (unsigned iflower=1; iflower<=num_flowers1; iflower++)
    {
        double h = double(iflower - 1) / (num_flowers1 - 1);

        auto param = params.density_params(h);
        if (!param) continue; // при малой высоте нет кривой

        density_distribution.param(param);

        double d = density_distribution(rangen);
        while (d < 0) d = density_distribution(rangen);
            // коррекция, чтобы не уйти в нефизическую область

        density[iflower] = d;
        sumd += d;
    }
    double norm = num_flowers2 / sumd;
    for (auto &e : density) e *= norm;
    for (auto &e : density) e = std::round(e);

    //--------------------------------------------------------------------------------------------
    //                          генерация времени старта побега
    //--------------------------------------------------------------------------------------------
    std::vector<double> tmp;
    for (unsigned i=0; i<50; i++)
    {
        double tmin = 1000.;
        for (unsigned iflower=1; iflower<=num_flowers1; iflower++)
        {
            double h = double(iflower - 1) / (num_flowers1 - 1);
            start_distribution.param(params.axis_params(h, 1));
            double ts = start_distribution(rangen);
            if (tmin > ts) tmin = ts;
        }
        tmp.push_back(tmin);
    }
    double tstart = average(tmp);

    //--------------------------------------------------------------------------------------------
    //                          генерация моделей осей
    //--------------------------------------------------------------------------------------------
    for (unsigned iflower=1; iflower<=num_flowers1; iflower++)
    {
        double h = double(iflower - 1) / (num_flowers1 - 1);

        //---------------------------------- генерация 1й оси ----------------------------------
        start_distribution.param(params.axis_params(h, 1));
        span_distribution.param(params.flowerspan_params(h, 1));
            // установили генераторы на параметры модельных кривых

        double ts = start_distribution(rangen);
        if (ts < tstart) ts = tstart;

        double dt = span_distribution(rangen);
        while (dt < 2) dt = span_distribution(rangen);
            // старт и длительность цветения на высоте h = average + среднекв.ошибка
            // + коррекция, чтобы не уйти в нефизическую область

        auto tb = budspan_distribution(rangen);
        pmodel->insert_flower(Flowerstamp{sprout, iflower, 0}, {ts - tb, ts, ts + dt});

        //---------------------------------- генерация 2й оси ----------------------------------
        if (density[iflower] == 0.) continue;

        //---------------------- левая граница  ----------------------
        start_distribution.param(params.axis_params(h, 2));
        span_distribution.param(params.flowerspan_params(h, 2));
        between_distribution.param(params.between_params(h, 2));
            // установили генераторы на параметры модельных кривых

        double ts2 = ts + start_distribution(rangen);
        dt = span_distribution(rangen); while (dt < 1) dt = span_distribution(rangen);

        for (unsigned i=1; i<=unsigned(density[iflower]); i++)
        {
            auto tb = budspan_distribution(rangen);
            pmodel->insert_flower(Flowerstamp{sprout, iflower, i}, {ts2 - tb, ts2, ts2 + dt});
            ts2 = ts2 + between_distribution(rangen);
            dt = span_distribution(rangen); while (dt < 1) dt = span_distribution(rangen);
        }
    }
    return pmodel;
}
