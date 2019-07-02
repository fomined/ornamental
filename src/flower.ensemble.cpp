#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "prg.os.h"
#include "prg.random.h"

#include "test.defs.h"
#include "math.defs.h"
#include "math.WLSM.h"
#include "math.statistics.h"
#include "data.recovery.h"
#include "flower.dome.h"
#include "flower.ensemble.h"

//#define USE_XMLM

#ifdef USE_XMLM   
std::pair<double, std::shared_ptr<species_ensemble>>
Flowerensemble::generate_ensemble(Unsigned num_iterations, bool) const
{
    //-----------------------------------------------------------------------------------------
    //                          рассчитаем суженные интервалы
    //-----------------------------------------------------------------------------------------
    Unsigned ne = bounds.size();
    
    auto return_values = intrainterval_localization_by_XMLM(ne, 3,
        &intervals[0], num_iterations);
    auto restored_values = return_values.second;

    //-----------------------------------------------------------------------------------------
    //                  вернем данные к формату возвращаемого значения
    //-----------------------------------------------------------------------------------------
    std::pair<double, std::shared_ptr<species_ensemble>>
        tmp(0, std::shared_ptr<species_ensemble>(new species_ensemble()));
    tmp.first = return_values.first;

    for (Unsigned i=0; i<ne; i++)
    {
        double t1 = restored_values[3*i    ];
        double t2 = restored_values[3*i + 1];
        double t3 = restored_values[3*i + 2];

        if (t1 >= t2) tmp.first = 10000; // игнорируем
        if (t2 >= t3) tmp.first = 10000; // игнорируем
        if (t1 >= t2 || t2 >= t3) std::printf("t1=%2.0lf t2=%2.0lf t3=%2.0lf\n", t1, t2, t3);
        
        tmp.second->push_back({t1, t2, t3});
    }

    return tmp;
}
#else
std::pair<double, std::shared_ptr<species_ensemble>>
Flowerensemble::generate_ensemble(Unsigned num_iterations, bool print_flag) const
{
    // typedef intraranges_uniform_generator generator_type;
    typedef intraranges_dispersion_generator generator_type;

    generator_type generator;
    //-----------------------------------------------------------------------------------------
    //                          рассчитаем суженные интервалы
    //-----------------------------------------------------------------------------------------
    Unsigned ne = bounds.size();
    
    auto restored_values = intrainterval_localization_by_MLM<generator_type>(ne, 3,
        &intervals[0], num_iterations, generator);

    //-----------------------------------------------------------------------------------------
    //                  вернем данные к формату возвращаемого значения
    //-----------------------------------------------------------------------------------------
    std::pair<double, std::shared_ptr<species_ensemble>>
        tmp(0, std::shared_ptr<species_ensemble>(new species_ensemble()));
    std::vector<double> spans;

    species_ensemble ensemble;
    for (Unsigned i=0; i<ne; i++)
    {
        double t1 = restored_values[3*i    ];
        double t2 = restored_values[3*i + 1];
        double t3 = restored_values[3*i + 2];

        if (t2 > t3) std::swap(t2, t3); //TODO заплатка
        spans.push_back(t3 - t2);
        tmp.second->push_back({t1, t2, t3});
    }
    //-----------------------------------------------------------------------------------------
    //                  сосчитаем логарифм вероятности ансамбля
    //-----------------------------------------------------------------------------------------
    double mu = average(spans);
    double sigma2 = variance(spans, mu);
    double N = spans.size();
    double p = 1. / std::pow(sigma2, 0.5 * N);
    tmp.first = sigma2;

    if (print_flag)
    {
        printf("dispersion = %e possibility = %e\n", sigma2, p); FLUSH;
    }
    return tmp;
}
#endif

void Flowerensemble::generate(Unsigned runs, Unsigned num_iterations, bool print_flag)
{
    for (unsigned i=0; i<runs; i++)
    {
        auto v = generate_ensemble(num_iterations, print_flag);
        concentrator.insert(v, [](const value_type &a, const value_type &b)
            { return std::abs(a.first - b.first); });
        //std::printf("%d run done.\n", i); FLUSH;
    }

    if (print_flag)
    {
        for (Unsigned i=0; i<bounds.size(); i++)
        {
            print_flowering_bounds(bounds, i);
            for (auto e : concentrator) print_flowering_phase(e.second, i);
        }
    }
}

interval<Timestamp> Flowerensemble::flowering_dates() const
{
    Timestamp ts(1000), te(0);

    auto first = concentrator.begin(), last = concentrator.end();
    for (; first!=last; ++first)
    {
        value_type record = first->second; // {variance, ptr->micro_ensemble}
        std::shared_ptr<exemple_type> example = record.second; // множество фаз цветения

        Unsigned num_flowers = example->size();
        for (Unsigned i=0; i<num_flowers; i++)
        {
            Flowerstamp fs = bounds[i].fsp;
            unsigned axis = fs.axis_order();

            auto r = (*example)[i];
            if (ts > r[1] && axis==1) ts = r[1];
            if (te < r[2]) te = r[2];
        }
    }
    return interval<Timestamp>(ts, te);
}

std::pair<double, double> Flowerensemble::flowering_start(unsigned sprout) const
{
    //------------------------------------------------------------------------------------
    //              получение номеров записей, удовлетворяющих фильтру
    //------------------------------------------------------------------------------------
    std::vector<unsigned> indexes; // индексы цветков, удовлетворяющих фильтру
    for (unsigned index=0; index<bounds.size(); index++)
    {
        Flowerstamp fs = bounds[index].fsp;
        if (fs.sprout() != sprout) continue;
        if (fs.axis_order() != 1) continue;
        indexes.push_back(index);
    }

    //------------------------------------------------------------------------------------
    //                  получение данных для выбранных цветков
    //------------------------------------------------------------------------------------
    curve_t tmp; double tstart = 1000., tsigma2 = 0.;
    for (auto i : indexes)
    {
        std::vector<double> ts;
        auto first = concentrator.begin(), last = concentrator.end();
        for (unsigned iensemble=0; first!=last; ++first, ++iensemble)
        {
            value_type record = first->second; // {variance, ptr->micro_ensemble}
            std::shared_ptr<exemple_type> example = record.second; // множество фаз цветения

            auto r = (*example)[i];
            ts.push_back(r[1]);
        }
        double mean = average(ts);
        double sigma2 = variance(ts, mean);
        if (tstart > mean) { tstart = mean; tsigma2 = sigma2; }
    }
    return {tstart, sqrt(tsigma2)};
}

std::pair<double, double> Flowerensemble::flowering_start(unsigned sprout, unsigned height) const
{
    std::vector<unsigned> indexes; // индексы цветков, удовлетворяющих фильтру

    int neighbours = 0;
    while (indexes.empty() && neighbours < 5)
    {
        for (unsigned index=0; index<bounds.size(); index++)
        {
            Flowerstamp fs = bounds[index].fsp;

            if (fs.sprout() != sprout) continue;
            if (fs.axis_order() != 1) continue;

            if (std::abs(int(fs.get(1)) - int(height)) > neighbours) continue;
                // усредняем по 2*neighbours ближайшим соседям

            indexes.push_back(index);
        }
        neighbours++;
    }
    assert(!indexes.empty());

    //------------------------------------------------------------------------------------
    //                  получение данных для выбранных цветков
    //------------------------------------------------------------------------------------
    curve_t tmp; double tstart = 1000., tsigma2 = 0.;
    for (auto i : indexes)
    {
        std::vector<double> ts;
        auto first = concentrator.begin(), last = concentrator.end();
        for (unsigned iensemble=0; first!=last; ++first, ++iensemble)
        {
            value_type record = first->second; // {variance, ptr->micro_ensemble}
            std::shared_ptr<exemple_type> example = record.second; // множество фаз цветения

            auto r = (*example)[i];
            ts.push_back(r[1]);
        }

        double mean= average(ts);
        double sigma2 = variance(ts, mean);
        if (tstart > mean) { tstart = mean; tsigma2 = sigma2; }
    }
    return {tstart, sqrt(tsigma2)};
}

unsigned Flowerensemble::flowering_height(unsigned sprout) const
{
    unsigned height = 0;
    for (auto e : bounds)
    {
        Flowerstamp fs = e.fsp;
        if (fs.sprout() != sprout) continue;

        unsigned pos = fs.get(1);
        if (height < pos) height = pos;
    }
    return height;
}

curve_t Flowerensemble::count_select(const std::vector<unsigned> &sprouts, unsigned axis) const
{
    std::map<unsigned, unsigned> num;

    for (auto e : sprouts)
    {
        auto n = dome.number_of_flowers(e, axis);
        if (num.count(n) == 0) num[n] = 0;
        num[n]++;
    }
    curve_t data;
    for (auto e : num) data.push_back({double(e.first), double(e.second)});
    return data;
}

curve_t Flowerensemble::phase_select(Unsigned sample, phase_enum_t phase, unsigned sprout, unsigned axis,
    unsigned char flags, unsigned pos) const
{
    assert(sample < concentrator.size());
    assert(sprout != 0); // побег должен быть определен

    unsigned phase_pos = 0; // BUD
    if (phase == FLOWER) phase_pos = 1;

    //------------------------------------------------------------------------------------
    //              получение номеров записей, удовлетворяющих фильтру
    //------------------------------------------------------------------------------------
    std::vector<unsigned> indexes; // индексы цветков, удовлетворяющих фильтру
    for (unsigned index=0; index<bounds.size(); index++)
    {
        Flowerstamp fs = bounds[index].fsp;

        if (fs.sprout() != sprout) continue;
        if (fs.axis_order() != axis) continue;
        if (pos != 0 && fs.axis_position() != pos) continue;

        indexes.push_back(index);
    }
    if (indexes.empty()) return curve_t();

    //------------------------------------------------------------------------------------
    //                  получение данных для выбранных цветков
    //------------------------------------------------------------------------------------
    curve_t tmp;
    for (auto i : indexes)
    {
        Flowerstamp fs = bounds[i].fsp;
        
        auto first = concentrator.begin();
        for (unsigned iensemble=0; iensemble!=sample; ++iensemble) ++first;
            // сдвиг на заданную позицию
        
        value_type record = first->second; // {variance, ptr->micro_ensemble}
        std::shared_ptr<exemple_type> example = record.second; // множество фаз цветения

        double height = double(fs.get(1));
            // высота как номер (временно, далее корректируем)

        auto r = (*example)[i];
        double ts = r[phase_pos];
        double dt = r[phase_pos+1] - r[phase_pos];

        tmp.push_back({height, ts, dt, 0., 0., double(fs.sprout())});
    }

    if (flags & ALIGN_SPROUT) // выравнивание по побегу (независимо от оси)
    {
        std::pair<double, double> m = flowering_start(sprout);
        for (auto &e : tmp) e[1] -= m.first;
    }

    if ((flags & ALIGN_1AXIS) && axis == 2) // выравнивание по 1й оси только для 2й оси
    {
        // подготовим заранее данные по стартам на разной высоте, чтобы уменьшить затраты
        unsigned max_height = flowering_height(sprout);

        std::vector<double> starts(max_height+1, 0); // старты на разных высотах
        for (unsigned ih=1; ih<=max_height; ih++)
        {
            std::pair<double, double> m = flowering_start(sprout, ih);
            starts[ih] = m.first;
        }
        for (auto &e : tmp)
        {
            unsigned ih = unsigned(std::round(e[0]));
            double start = starts[ih];
            e[1] -= start;
        }
    }

    if (flags & NORM) // исправление высот ПОСЛЕ выравнивания, так как исходные высоты используются
    {
        unsigned height = flowering_height(sprout);
        double norm = 1. / (height - 1);
        for (auto &e : tmp) e[0] = (e[0] - 1) * norm;
    }

    if (flags & DAY_ALIGN_CODE) // выровнять границы диапазона (колонки 1, 2)
        tmp = prg::to_int<CDIM>(tmp, {1, 2});

    return tmp;
}

curve_t Flowerensemble::phase_select(Unsigned sample, phase_enum_t phase, const std::vector<unsigned> &sprouts,
    unsigned axis, unsigned char flags, unsigned pos) const
{
    curve_t data; // точечные данные по цветению всех цветков
    for (auto sprout : sprouts)
    {
        auto tmp = phase_select(sample, phase, sprout, axis, flags, pos);
            // выбрать на побеге (sprout) цветки 1-го и 2-го порядка (все=0), выровнять (ALIGN)
            // по побегу и нормировать (NORM) на число цветков в оси 1-го порядка
        if (tmp.empty()) continue;

        data = prg::merge<CDIM>(data, tmp, 0);
            // объединить данные с сортировкой по колонке 0 (height)
    }
    return data;
}

curve_t Flowerensemble::phase_select(phase_enum_t phase, unsigned sprout, unsigned axis,
    unsigned char flags, unsigned pos) const
{
    assert(sprout != 0); // побег должен быть определен

    unsigned phase_pos = 0; // BUD
    if (phase == FLOWER) phase_pos = 1;

    //------------------------------------------------------------------------------------
    //              получение номеров записей, удовлетворяющих фильтру
    //------------------------------------------------------------------------------------
    std::vector<unsigned> indexes; // индексы цветков, удовлетворяющих фильтру
    for (unsigned index=0; index<bounds.size(); index++)
    {
        Flowerstamp fs = bounds[index].fsp;

        if (fs.sprout() != sprout) continue;
        if (fs.axis_order() != axis) continue;
        if (pos != 0 && fs.axis_position() != pos) continue;

        indexes.push_back(index);
    }
    if (indexes.empty()) return curve_t();

    //------------------------------------------------------------------------------------
    //                  получение данных для выбранных цветков
    //------------------------------------------------------------------------------------
    curve_t tmp;
    for (auto i : indexes)
    {
        Flowerstamp fs = bounds[i].fsp;
        
        std::vector<double> ts;
        std::vector<double> dt;

        auto first = concentrator.begin(), last = concentrator.end();
        for (unsigned iensemble=0; first!=last; ++first, ++iensemble)
        {
            value_type record = first->second; // {variance, ptr->micro_ensemble}
            std::shared_ptr<exemple_type> example = record.second; // множество фаз цветения

            auto r = (*example)[i];
            ts.push_back(r[phase_pos]);
            dt.push_back(r[phase_pos+1] - r[phase_pos]);
        }

        double height = double(fs.get(1));
            // высота как номер (временно, далее корректируем)

        double meants = average(ts);
        double sigmats = variance(ts, meants);
        double meandt = average(dt);
        double sigmadt = variance(dt, meandt);

        tmp.push_back({height, meants, meandt, sqrt(sigmats), sqrt(sigmadt), double(fs.sprout())});
    }

    if (flags & ALIGN_SPROUT) // выравнивание по побегу (независимо от оси)
    {
        std::pair<double, double> m = flowering_start(sprout);
        for (auto &e : tmp) e[1] -= m.first;
    }

    if ((flags & ALIGN_1AXIS) && axis == 2) // выравнивание по 1й оси только для 2й оси
    {
        // подготовим заранее данные по стартам на разной высоте, чтобы уменьшить затраты
        unsigned max_height = flowering_height(sprout);

        std::vector<double> starts(max_height+1, 0); // старты на разных высотах
        for (unsigned ih=1; ih<=max_height; ih++)
        {
            std::pair<double, double> m = flowering_start(sprout, ih);
            starts[ih] = m.first;
        }
        for (auto &e : tmp)
        {
            unsigned ih = unsigned(std::round(e[0]));
            double start = starts[ih];
            e[1] -= start;
        }
    }

    if (flags & NORM) // исправление высот ПОСЛЕ выравнивания, так как исходные высоты используются
    {
        unsigned height = flowering_height(sprout);
        double norm = 1. / (height - 1);
        for (auto &e : tmp) e[0] = (e[0] - 1) * norm;
    }

    if (flags & DAY_ALIGN_CODE) // выровнять границы диапазона (колонки 1, 2)
        tmp = prg::to_int<CDIM>(tmp, {1, 2});

    return tmp;
}

curve_t Flowerensemble::phase_select(phase_enum_t phase, const std::vector<unsigned> &sprouts,
    unsigned axis, unsigned char flags, unsigned pos) const
{
    curve_t data; // точечные данные по цветению всех цветков
    for (auto sprout : sprouts)
    {
        auto tmp = phase_select(phase, sprout, axis, flags, pos);
            // выбрать на побеге (sprout) цветки 1-го и 2-го порядка (все=0), выровнять (ALIGN)
            // по побегу и нормировать (NORM) на число цветков в оси 1-го порядка
        if (tmp.empty()) continue;

        data = prg::merge<CDIM>(data, tmp, 0);
            // объединить данные с сортировкой по колонке 0 (height)
    }

    return data;
}

curve_t construct_flowerspan(const curve_t &data, unsigned num_flowers)
{
    //---------------------------------------------------------------------------------------------
    //                          параметры интерполяции данных
    //---------------------------------------------------------------------------------------------
    const double indent = 3; // добавочный отступ при отрисовке

    //---------------------------------------------------------------------------------------------
    //                      расчет гистограммы цветения
    //---------------------------------------------------------------------------------------------
    curve_t flowerspan;
    std::map<double, unsigned> map; double tmin = 100., tmax = -100.;
    for (auto e : data)
    {
        unsigned dt = unsigned(std::round(e[2]));
        if (map.count(dt) == 0) { map[dt] = 0; if (tmin > dt) tmin = dt; if (tmax < dt) tmax = dt; }
        ++map[dt];
    }
    for (double i=tmin-indent; i<tmin; i+=1) flowerspan.push_back({i, 0});
    for (auto e : map) flowerspan.push_back({e.first, 100*double(e.second)/num_flowers });
    for (double i=tmax+1; i<tmax+indent; i+=1) flowerspan.push_back({i, 0});

    return flowerspan;
}

hidden_parameters_type Flowerensemble::get_hidden_parameters(const std::vector<unsigned> &sprouts) const
{
    hidden_parameters_type params;

    unsigned power = 3; // степень аппроксимирующего полинома
    double wy = .15; // ширина окна усреднения для 1й оси
    double wyb = .4; // ширина окна усреднения для границ 2й оси
    double wyd = .05; // ширина окна усреднения для плотности

    double hy = 0.0125; // расстояния между точками при разложении по полиномам
    double thickness = THICKNESS; // доля объектов принадлежащих границе

    std::vector<double> num1, num2; // число цветков 1го и 2го порядков
    curve_t data1; // точечные данные по расцветанию всех цветков 1-го порядка
    curve_t data2; // точечные данные по расцветанию всех цветков 2-го порядка
    curve_t data4; // точечные данные по цветанию всех цветков 2-го порядка

//    curve_t numc; numc.push_back({0, 0});
    for (auto sprout : sprouts)
    {
        auto tmp1 = phase_select(FLOWER, sprout, 1, NORM | ALIGN_SPROUT, 0);
        auto tmp2 = phase_select(FLOWER, sprout, 2, NORM | ALIGN_1AXIS, 0);
        auto tmp4 = phase_select(FLOWER, sprout, 2, NORM | ALIGN_1AXIS, 0);
            // относительно 1й оси, чтобы улучшить сглаживание

        num1.push_back(tmp1.size());
        num2.push_back(double(tmp2.size()) / tmp1.size());
            // доля, относительно числа цветков на 1й оси

        // объединить данные с сортировкой по колонке 0
        data1 = prg::merge<CDIM>(data1, tmp1, 0);
        data2 = prg::merge<CDIM>(data2, tmp2, 0);
        data4 = prg::merge<CDIM>(data4, tmp4, 0);
    }

    //----------------------- параметры числа цветков на осях -----------------------
    params.num_flowers[1][0] = average(num1);
    params.num_flowers[1][1] = sqrt(variance(num1, params.num_flowers[1][0]));
    params.num_flowers[2][0] = average(num2);
    params.num_flowers[2][1] = sqrt(variance(num2, params.num_flowers[2][0]));

    //----------------------- параметры 1й оси -----------------------
    params.axis[1] = wslm_simple_smooth<CDIM>(0, power, data1, 0, 1, wy, hy);

    //----------------------- параметры 2й оси -----------------------
    auto ldata2 = lower_bound<CDIM>(data4, thickness);
    params.axis[2] = wslm_simple_smooth<CDIM>(0, power, ldata2, 0, 1, wyb, hy);

    //----------------------- параметры плотности -----------------------
    curve_t density_data;
    {
        typedef std::pair<double, unsigned> key_t;

        std::map<key_t, unsigned> sum;
        for (auto e : data4)
        {
            key_t key{e[0], e[5]}; // с учетом того, к какому побегу e[5] принадлежит точка
            if (sum.count(key) == 0) sum[key] = 0;
            sum[key]++;
        }
		for (auto e : sum) density_data.push_back({e.first.first, double(e.second)});
    }

    params.density = wslm_simple_smooth<CDIM>(0, power, density_data, 0, 1, wyd, hy);

    //----------------------- длительность цветения -----------------------
    params.flowerspan[1] = wslm_simple_smooth<CDIM>(0, power, data1, 0, 2, wy, hy);

    curve_t flowerspan;
    std::set<double> heights; // набор уникальных высот
    for (auto e : data4) heights.insert(e[0]);

    for (auto h : heights)
    {
        auto first = std::lower_bound(data4.begin(), data4.end(), h,
            [](const point<CDIM> &v, double h){ return v[0] < h; });
        auto last = std::upper_bound(data4.begin(), data4.end(), h,
            [](double h, const point<CDIM> &v){ return h < v[0]; });
        assert(first != last);
            // на 2й оси используем накопление по всем цветкам одной высоты

        std::vector<double> tmp;
        for (; first!=last; ++first) tmp.push_back((*first)[2]);

        flowerspan.push_back({h, average(tmp), 0.});
    }
    params.flowerspan[2] = wslm_simple_smooth<CDIM>(0, power, flowerspan, 0, 1, wy, hy);

    //------------------- промежутки между цветением 2й оси --------------------

    std::sort(data4.begin(), data4.end(),
        [](const point<CDIM> &a, const point<CDIM> &b)
        {
            if (a[0] < b[0]) return true;
            if (a[0] > b[0]) return false;
            return a[1] < b[1];
        }
    );

    curve_t between;
    for (auto h : heights)
    {
        auto first = std::lower_bound(data4.begin(), data4.end(), h,
            [](const point<CDIM> &v, double h){ return v[0] < h; });
        auto last = std::upper_bound(data4.begin(), data4.end(), h,
            [](double h, const point<CDIM> &v){ return h < v[0]; });

        std::vector<double> tmp; double prev = (*first)[1];
        for (++first; first!=last; ++first)
        {
            tmp.push_back((*first)[1] - prev);
            prev = (*first)[1];
        }
        double dt = tmp.size() < 1 ? 0 : average(tmp);
        between.push_back({h, dt, 0.});
    }
    params.between = wslm_simple_smooth<CDIM>(0, power, between, 0, 1, wy, hy);

    params.ready_flag = true;
    return params;
}

void print_flowering_bounds(const std::vector<flowerline3_t> &bounds, Unsigned i)
{
    Flowerstamp flower = bounds[i].fsp;

    auto r1 = bounds[i].intervals[1];
    int tb1 = r1[0].to_int();
    int tb2 = r1[1].to_int();

    auto r2 = bounds[i].intervals[2];
    int tc1 = r2[0].to_int();
    int tc2 = r2[1].to_int();

    int indent = 5;
    std::string ranges(unsigned(tc2 + indent), ' ');
    for (int i=tb1; i<=tc2; i++)
    {
        char ch = ' ';
        if (i <= tb2)
        {
          if (i < tc1) ch = '-';
          else ch = '.';
        }
        else ch = '+';
        ranges[unsigned(i)] = ch;
    }
    printf("%8s: %s\n", flower.to_string().c_str(), ranges.c_str()); FLUSH;
}

void print_flowering_phase(const std::pair<double, std::shared_ptr<species_ensemble>> &ens,
    Unsigned i)
{
    int indent = 5;

    std::shared_ptr<species_ensemble> p = ens.second;
    auto r = (*p)[i];

    std::string current(unsigned(r[2] + indent), ' ');
    int t1 = int(std::round(r[1]));
    int t2 = int(std::round(r[2]));
    for (int i=t1; i<=t2; i++) current[unsigned(i)] = '*';
    printf("%8s: %s\n", "", current.c_str());
}

