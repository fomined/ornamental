#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <random>

#include "data.recovery.h"
#include "math.statistics.h"
#include "math.WLSM.h"
#include "gsl/prg.interpolation.h"

/*!
 * \brief генератор случайных дат внутри заданных диапазонов с контролем пересечений интервалов
 */
struct intraranges_date_generator
{
    typedef std::uniform_int_distribution<int> distribution_type;
    typedef distribution_type::param_type param_type;
    distribution_type distribution; // генератор распределения

    std::array<int, 3> operator()(const interval<double> *range)
    {
        int num_iteration = 300; // для избегания зацикливания

        // генерим значение в интервале смены BUD-FLOWER и от него в обе стороны, используя
        // нормальные распределения длин фаз, пытаемся попасть в интервалы *-BUD и FLOWER-OVARY
        std::array<int, 3> t;

        interval<double> r0 = range[0]; // интервал *-BUD
        interval<double> r1 = range[1]; // интервал BUD-FLOWER
        interval<double> r2 = range[2]; // интервал FLOWER-OVARY

        while (num_iteration > 0)
        {
            num_iteration--;

        #define PARAMS(range) param_type{int(std::round(range[0])), int(std::round(range[1]))}
        
            distribution.param(PARAMS(r1));
            t[1] = distribution(rangen);

            distribution.param(PARAMS(r0));
            t[0] = distribution(rangen);
            if (t[0] >= t[1]) continue;

            distribution.param(PARAMS(r2));
            t[2] = distribution(rangen);
            if (t[2] <= t[1]) continue;
            
        #undef PARAMS        
       
            break;
        };
        assert(num_iteration > 0);
        
        return t;
    }
};

/// k=0 (BUD), k=1(FLOWER)
curve_t make_histogram(Unsigned ne, const std::vector<int> &x, Unsigned k)
{
    assert(k == 0 || k == 1);
    
    std::map<int, int> hist;
    for (Unsigned i=0; i<ne; i++)
    {
        auto m = x[3*i + k+1] - x[3*i + k];
        if (hist.count(m) == 0) hist[m] = 0;
        hist[m]++;
    }

    curve_t data;    
    for (auto e : hist) data.push_back(point_t{double(e.first), double(e.second)});
    return data;
}

/// рассчитать кумулятивную функцию плотности распределения
std::pair<curve_t, curve_t> prepare_CDF(const curve_t &v1)
{
    //---------------------------------------------------------------------------------------------
    //                          рассчитаем функцию распределения
    //---------------------------------------------------------------------------------------------
    const unsigned nmid = 8; // число промежуточных точек
    double indent1 = std::max(0.5 * (v1.back()[0] - v1.front()[0]), double(v1.front()[0]));

    curve_t data;
    double prevx = v1.front()[0] - indent1, cnt = 0, dx = 0;
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
    dx = indent1 / nmid;
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

    //---------------------------------------------------------------------------------------------
    //                      расчет гистограммы цветения
    //---------------------------------------------------------------------------------------------
    return std::make_pair(data, wslm_simple_smooth<CDIM>(0, power, data, 0, 1, wx1, hx));
}

/// отклонение в заданной точке
double rmsd(const curve_t &f)
{
    double s2 = 0;
    for (auto e : f) s2 += sqr(e[2]);
    return sqrt(s2);
}

/// отклонение в заданной точке
double variance_change(const std::pair<curve_t, curve_t> &f, double x0, double x1)
{
    auto it = std::lower_bound(f.second.begin(), f.second.end(), x0, 
        [](const point_t &a, double x) { return a[0] < x; }
    );
    if (it == f.second.end()) --it;
    double y0 = (*it)[1]; // значение в кривой апроксимации f(x0)

    it = std::lower_bound(f.second.begin(), f.second.end(), x1, 
        [](const point_t &a, double x) { return a[0] < x; }
    );
    if (it == f.second.end()) --it;
    double y1 = (*it)[1]; // значение в кривой апроксимации f(x1)
    

    it = std::lower_bound(f.first.begin(), f.first.end(), x0, 
        [](const point_t &a, double x) { return a[0] < x; }
    );
    if (it == f.first.end()) --it;
    double q0 = (*it)[1]; // значение в данных g(x0)
    assert(q0 >= 1); 

    it = std::lower_bound(f.first.begin(), f.first.end(), x1, 
        [](const point_t &a, double x) { return a[0] < x; }
    );
    if (it == f.first.end()) --it;
    double q1 = (*it)[1]; // значение в данных g(x1)
    
    
    double before = sqr(y0 - q0) + sqr(y1 - q1);
    double after  = sqr(y0 - (q0 - 1)) + sqr(y1 - (q1 + 1));

    return after - before;
}

/*
 * \brief восстановление "точных" значений границ фаз методом улучшенного максимального правдоподобия
 * \param ne число экземляров во множестве экземпляров
 * \param ni число интервалов в каждом экземпляре
 * \param ps указатель на массив всех интервалов
 * \param num_iterations число итераций поиска
 */
std::pair<double, std::vector<double>> intrainterval_localization_by_XMLM(Unsigned ne, Unsigned ni,
    const interval<double> *ranges, Unsigned num_iterations)
{
    assert(ni == 3); // ограничение текущей версии
    
    Unsigned np = ne * ni; // полное число оптимизирумых параметров
    std::vector<int> x(np, 0); // настраиваемый для возвращения объект

    intraranges_date_generator generator;
    
    //----------------------------------------------------------------------------------------------
    //              рандомизация стартовых значений оптимизируемых параметров
    //----------------------------------------------------------------------------------------------
    for (Unsigned i=0; i<ne; i++)
    {
        std::array<int, 3> r = generator(&ranges[3*i]);
        x[3*i  ] = r[0];
        x[3*i+1] = r[1];
        x[3*i+2] = r[2];
    }
    
    auto CDF1 = prepare_CDF(make_histogram(ne, x, 0)); // BUD {CDF data, curve} 
    auto CDF2 = prepare_CDF(make_histogram(ne, x, 1)); // FLOWER {CDF data, curve} 
    
    double prevrmsd1 = rmsd(CDF1.second);
    double prevrmsd2 = rmsd(CDF2.second);
    double prevrmsd = prevrmsd1 + prevrmsd2;
    
    double newrmsd1 = 0;
    double newrmsd2 = 0;
    double newrmsd = newrmsd1 + newrmsd2;
    
//prg::printf(CDF1.first);
//SEPARATOR_LINE(80);
//prg::printf(CDF1.second);
//SEPARATOR_LINE(80);
//SEPARATOR_LINE(80);
//prg::printf(CDF2.first);
//SEPARATOR_LINE(80);
//prg::printf(CDF2.second);
//SEPARATOR_LINE(80);
//SEPARATOR_LINE(80);
//SEPARATOR_LINE(80);
    
    //---------------------------------------------------------------------------------------------
    //                              оптимизация x[i]
    //---------------------------------------------------------------------------------------------
    std::uniform_int_distribution<Unsigned> distribution(0, np-1);
    std::bernoulli_distribution leftdist(0.5); // вероятность сдвига вперед или назад
    double eps = 0.3;
    
    for (Unsigned iteration=0; iteration<num_iterations && prevrmsd > eps; iteration++)  
    {
        bool flag[2] = {false, false};
        Unsigned i = distribution(rangen);
    
        int a = int(ranges[i][0]); // правая дата
        int b = int(ranges[i][1]); // разрешенная левая дата
        assert(a <= b);
        
        int xx = x[i]; 
        assert(a <= xx);
        assert(xx <= b);
        
        if (i % 3 == 0) // меняет BUD
        {
            b = std::min(b, x[i+1]);
            flag[0] = true;
        }
        else if (i % 3 == 1) // меняет BUD, FLOWER
        { 
            a = std::max(x[i-1], a); 
            b = std::min(b, x[i+1]); 
            flag[0] = true;
            flag[1] = true;
        }
        else if (i % 3 == 2) // меняет FLOWER
        { 
            a = std::max(x[i-1], a); 
            flag[1] = true;
        }
        assert(a <= b);
        if (a == b) continue; 
                
        bool leftdir = leftdist(rangen);
        if (leftdir && xx == a) continue;
        if (!leftdir && xx == b) continue;
        int t = leftdir ? xx - 1 : xx + 1;

        x[i] = t; // пробуем новое значение
        
        
//        double x0 = xx - x[i-1]; // прежняя длина фазы BUD
//        double z0 = x[i+1] - xx; // прежняя длина фазы FLOWER
//        double x1 = xx - t; // новая длина фазы BUD
//        double z1 = x[i+1] - t; // новая длина фазы FLOWER
//            double prevrmsd = rmsd(FLW_CDF.second);
//            double d1 = variance_change(BUD_CDF, x0, x1);
//            double d2 = variance_change(FLW_CDF, z0, z1);
        
        
        std::pair<curve_t, curve_t> tmp1, tmp2;
        newrmsd1 = 0; 
        newrmsd2 = 0;

        if (flag[0])
        {
            tmp1 = prepare_CDF(make_histogram(ne, x, 0));
            newrmsd1 = rmsd(tmp1.second);
        }
        if (flag[1])
        {
            tmp2 = prepare_CDF(make_histogram(ne, x, 1));
            newrmsd2 = rmsd(tmp2.second);
        }
        newrmsd = newrmsd1 + newrmsd2;
        printf("%ld %lf -> %lf (%d)\n", iteration, prevrmsd, newrmsd, prevrmsd > newrmsd); FLUSH;
        
        if (newrmsd >= prevrmsd) x[i] = xx;  // восстанавливаем старое значение
        else 
        {
            CDF1 = tmp1;
            CDF2 = tmp2;
            prevrmsd1 = newrmsd1; 
            prevrmsd2 = newrmsd2; 
            prevrmsd = newrmsd; 
        }
    }

    
//prg::printf(CDF1.first);
//SEPARATOR_LINE(80);
//prg::printf(CDF1.second);
//SEPARATOR_LINE(80);
//SEPARATOR_LINE(80);
//prg::printf(CDF2.first);
//SEPARATOR_LINE(80);
//prg::printf(CDF2.second);
//SEPARATOR_LINE(80);

    std::pair<double, std::vector<double>> tmp;
    tmp.first = prevrmsd;
    for (auto e : x) tmp.second.push_back(e);
    return tmp;
}

//std::vector<double> intrainterval_localization_by_XMLM(Unsigned ne, Unsigned ni,
//    const interval<double> *ranges, Unsigned num_iterations)
//{
//    assert(ni == 3); // ограничение текущей версии
    
//    Unsigned np = ne * ni; // полное число оптимизирумых параметров
//    std::vector<int> x(np, 0); // настраиваемый для возвращения объект

//    intraranges_date_generator generator;
    
//    //----------------------------------------------------------------------------------------------
//    //              рандомизация стартовых значений оптимизируемых параметров
//    //----------------------------------------------------------------------------------------------
//    for (Unsigned i=0; i<ne; i++)
//    {
//        std::array<int, 3> r = generator(&ranges[3*i]);
//        x[3*i  ] = r[0];
//        x[3*i+1] = r[1];
//        x[3*i+2] = r[2];
//    }
    
//    auto BUD_CDF = prepare_CDF(make_histogram(ne, x, 0)); // BUD {CDF data, curve} 
//    auto FLW_CDF = prepare_CDF(make_histogram(ne, x, 1)); // FLOWER {CDF data, curve} 
    
//    //---------------------------------------------------------------------------------------------
//    //                              оптимизация x[i]
//    //---------------------------------------------------------------------------------------------
//    std::uniform_int_distribution<Unsigned> distribution(0, np-1);
//    std::uniform_int_distribution<int> dist;
//    typedef std::uniform_int_distribution<int>::param_type param_type;
    
//    for (Unsigned iteration=0; iteration<num_iterations; iteration++)  
//    {
//        Unsigned i = distribution(rangen);
    
//        int a = int(ranges[i][0]);
//        int b = int(ranges[i][1]);
//            // границы изменения x[i]
        
//        int xx = x[i]; 
//        if (i % 3 == 0) // есть влияние только на BUD
//        {
//            b = std::min(b, x[i+1]);
//            if (a == b) continue; 
//            dist.param(param_type{a, b});
            
//            int t = dist(rangen);
//            if (t == xx) continue;

//            double x0 = x[i+1] - x[i]; // прежняя длина фазы BUD
//            double x1 = x[i+1] - t; // новая длина фазы BUD
//            double d = variance_change(BUD_CDF, x0, x1);
//            if (d >= 0) continue;
//            x[i] = t;
//            BUD_CDF = prepare_CDF(make_histogram(ne, x, 0));
//        }
//        else if (i % 3 == 1) 
//        { 
//            a = std::max(x[i-1], a); 
//            b = std::min(b, x[i+1]); 
//            if (a == b) continue; 
//            dist.param(param_type{a, b});
            
//            int t = dist(rangen);
//            if (t == xx) continue;
            
//            double x0 = x[i] - x[i-1]; // прежняя длина фазы BUD
//            double z0 = x[i+1] - x[i]; // прежняя длина фазы FLOWER
//            double x1 = x[i] - t; // новая длина фазы BUD
//            double z1 = x[i+1] - t; // новая длина фазы FLOWER
            
//            double prevrmsd = rmsd(FLW_CDF.second);

//            double d1 = variance_change(BUD_CDF, x0, x1);
//            double d2 = variance_change(FLW_CDF, z0, z1);
//            if (d1 + d2 >= 0) continue;
            
//            x[i] = t;
//            BUD_CDF = prepare_CDF(make_histogram(ne, x, 0));
//            FLW_CDF = prepare_CDF(make_histogram(ne, x, 1));

//            double newrmsd = rmsd(FLW_CDF.second);
//            printf("(1) %ld %lf -> %lf\n", iteration, prevrmsd, newrmsd); FLUSH;
//        }
//        else if (i % 3 == 2) 
//        {
//            a = std::max(x[i-1], a);
//            if (a == b) continue; 
//            dist.param(param_type{a, b});
            
//            int t = dist(rangen);
//            if (t == xx) continue;

//            double z0 = x[i] - x[i-1]; // прежняя длина фазы FLOWER
//            double z1 = t - x[i-1]; // новая длина фазы FLOWER
            
//            double prevrmsd = rmsd(FLW_CDF.second);
//            double d = variance_change(FLW_CDF, z0, z1);
//            if (d >= 0) continue;
            
//            x[i] = t;
//            FLW_CDF = prepare_CDF(make_histogram(ne, x, 1));

//            double newrmsd = rmsd(FLW_CDF.second);
//            printf("(2) %ld %lf -> %lf\n", iteration, prevrmsd, newrmsd); FLUSH;
//        } 
//    }
//exit(0);
//    std::vector<double> tmp;
//    for (auto e : x) tmp.push_back(e);
//    return tmp;
//}
 
//std::vector<double> intrainterval_localization_by_XMLM(Unsigned ne, Unsigned ni,
//    const interval<double> *ranges, Unsigned num_iterations)
//{
//    assert(ni == 3); // ограничение текущей версии
    
//    Unsigned np = ne * ni; // полное число оптимизирумых параметров
//    std::vector<int> x(np, 0); // настраиваемый для возвращения объект

//    intraranges_date_generator generator;
//    //----------------------------------------------------------------------------------------------
//    //              вычисление стартовых значений оптимизируемых параметров
//    //----------------------------------------------------------------------------------------------

//    std::map<int, int> hist[2]; // гистограммы фаз (BUD, FLOWER)
//    for (Unsigned i=0; i<ne; i++)
//    {
//        std::array<int, 3> r = generator(&ranges[3*i]);
//        x[3*i  ] = r[0];
//        x[3*i+1] = r[1];
//        x[3*i+2] = r[2];
        
//        auto m = r[1] - r[0];
//        if (hist[0].count(m) == 0) hist[0][m] = 0;
//        hist[0][m]++;

//        m = r[2] - r[1];
//        if (hist[1].count(m) == 0) hist[1][m] = 0;
//        hist[1][m]++;

//    //#define DEBUG_PRINT
//    #ifdef DEBUG_PRINT
//        Timestamp ta1 = ps[3*i][0];
//        Timestamp ta2 = ps[3*i][1];

//        Timestamp tb1 = ps[3*i+1][0];
//        Timestamp tb2 = ps[3*i+1][1];

//        Timestamp tc1 = ps[3*i+2][0];
//        Timestamp tc2 = ps[3*i+2][1];

//        Timestamp t1 = x[3*i    ];
//        Timestamp t2 = x[3*i + 1];
//        Timestamp t3 = x[3*i + 2];

//        printf("   : [%6s %6s %6s] [%6s %6s %6s] [%6s %6s %6s] \n",
//            ta1.to_string().c_str(), t1.to_string().c_str(), ta2.to_string().c_str(),
//            tb1.to_string().c_str(), t2.to_string().c_str(), tb2.to_string().c_str(),
//            tc1.to_string().c_str(), t3.to_string().c_str(), tc2.to_string().c_str()
//        );
//        FLUSH;
//    #endif
//    #undef DEBUG_PRINT
//    }
    
//    //---------------------------------------------------------------------------------------------
//    //              сосчитаем кумулятивные функции плотности распределения
//    //---------------------------------------------------------------------------------------------
//    curve_t data[2]; 
//    for (auto e : hist[0]) data[0].push_back(point_t{double(e.first), double(e.second)});
//    for (auto e : hist[1]) data[1].push_back(point_t{double(e.first), double(e.second)});
    
//    auto BUD_CDF = prepare_CDF(data[0]); // BUD {data, curve} 
//    auto FLW_CDF = prepare_CDF(data[1]); // FLOWER {data, curve} 
    
//    //---------------------------------------------------------------------------------------------
//    //                              оптимизация x[i]
//    //---------------------------------------------------------------------------------------------
    
//    struct note
//    {
//        Unsigned i; // номер записи
//        int dx; // изменение начала фазы
//        double var2; // изменение варианты
//    };
//    std::vector<note> varchange; // массив объектов уменьшающих варианту
    
//    for (Unsigned i=0; i<np; i++)
//    {
//        int a = int(ranges[i][0]);
//        int b = int(ranges[i][1]);
//            // границы изменения x[i]
        
//        if (i % 3 == 0) // есть влияние только на BUD
//        {
//            b = std::min(b, x[i+1]);
//            if (a == b) continue; 

//            double x0 = x[i+1] - x[i]; // прежняя длина фазы
//            for (int k=a; k<b; k++)
//            {
//                double x1 = x[i+1] - k; // новая длина фазы
//                double d = variance_change(BUD_CDF, x0, x1);
//                if (d >= 0) continue;
                
//                varchange.push_back({i, k - x[i], d});
//            }
//        }
//        else if (i % 3 == 1) 
//        { 
//            a = std::max(x[i-1], a); 
//            b = std::min(b, x[i+1]); 
//            if (a == b) continue; 
            
//            double x0 = x[i] - x[i-1]; // прежняя длина фазы BUD
//            double z0 = x[i+1] - x[i]; // прежняя длина фазы FLOWER
//            for (int k=a; k<b; k++)
//            {
//                double x1 = x[i] - k; // новая длина фазы BUD
//                double z1 = x[i+1] - k; // новая длина фазы FLOWER
                
//                double d1 = variance_change(BUD_CDF, x0, x1);
//                double d2 = variance_change(FLW_CDF, z0, z1);
//                if (d1 + d2 >= 0) continue;
                
//                varchange.push_back({i, k - x[i], d1 + d2});
//            }
//        }
//        else if (i % 3 == 2) 
//        {
//            a = std::max(x[i-1], a);
//            if (a == b) continue; 
            
//            double z0 = x[i+1] - x[i]; // прежняя длина фазы FLOWER
//            for (int k=a; k<b; k++)
//            {
//                double z1 = x[i+1] - k; // новая длина фазы FLOWER
//                double d = variance_change(FLW_CDF, z0, z1);
//                if (d >= 0) continue;
                
//                varchange.push_back({i, k - x[i], d});
//            }
//        }
//    }

//    std::sort(varchange.begin(), varchange.end(), 
//        [](const note &a, const note &b){ return a.var2 < b.var2; }
//    );

//    std::vector<double> tmp;
//    return tmp;
//}


std::pair<double,double> intrainterval_localization_by_directMLM(
    std::vector<interval<double>> &restored_values,
    const std::vector<interval<double>> &lb,
    const std::vector<interval<double>> &rb)
{
    size_t N = lb.size();

    //----------------------------------------------------------------------------------------------
    //                  нахождение границ возможных значений среднего значения
    //----------------------------------------------------------------------------------------------
    double min_dt = 100, max_dt = 0;
    for (size_t i=0; i<N; i++)
    {
        double min_range = rb[i][0] < lb[i][1] ? 1 : rb[i][0] - lb[i][1] + 2;
            // превращение интервала "отк, закр" -> "закр, отк", т.е.") [" -> "[)" дает +2
        double max_range = rb[i][1] - lb[i][0];
        if (min_dt > min_range) min_dt = min_range;
        if (max_dt < max_range) max_dt = max_range;
    }

    //----------------------------------------------------------------------------------------------
    //                          получение среднего значения и дисперсии
    //----------------------------------------------------------------------------------------------

    // расчет убывающей функции f(x), прохождение которой через 0 указывает на ML среднее значение
    // и поиск для нее последнего (и единственного) интервала прохода через 0 и расчет среднего

    std::pair<point<2>, point<2>> zero_interval;
    size_t zeros = 0; // интервал перехода и число найденных таких интервалов
    {
        std::vector<double> tmp;
            // множество промежутков dt_i, которые получатся оптимизацией методом ML при попытке
            // установить среднее значение <dt>

        for (size_t i=0; i<N; i++)
        {
            double min_range = rb[i][0] < lb[i][1] ? 1 : rb[i][0] - lb[i][1] + 2;
            double max_range = rb[i][1] - lb[i][0];

            if (min_dt <= min_range) tmp.push_back(min_range);
            else if (min_dt >= max_range) tmp.push_back(max_range);
            else tmp.push_back(min_dt);
        }
        double prevf = average(tmp) - min_dt;

        for (double dt=min_dt+1; dt<=max_dt; dt++)
        {
            tmp.clear();
            for (size_t i=0; i<N; i++)
            {
                double min_range = rb[i][0] < lb[i][1] ? 1 : rb[i][0] - lb[i][1] + 2;
                double max_range = rb[i][1] - lb[i][0];

                if (dt <= min_range) tmp.push_back(min_range);
                else if (dt >= max_range) tmp.push_back(max_range);
                else tmp.push_back(dt);
            }
            double f = average(tmp) - dt;
            if (positive(prevf) && !positive(f))
            {
                zero_interval = std::make_pair(point<2>{dt-1, prevf}, point<2>{dt, f});
                zeros++;
            }
            prevf = f;
        }
    }
    assert(zeros == 1);

    point<2> a = zero_interval.first, b = zero_interval.second;
    double mean = prg::linear_interpolation(a[0], a[1], b[0], b[1]);
        // среднее значение для ML метода

    // дисперсия может быть расчитана как с учетом полностью подстаивающихся элементов, что
    // дает минимально возможное значение дисперсии, так и с исключением этих элементов, что
    // дает максимально возможную дисперсию. Расчет с исключением элементов требует
    // обоснования, так как приводит к увеличению среднего.

    std::vector<double> tmp;
    for (size_t i=0; i<N; i++)
    {
        double min_range = rb[i][0] < lb[i][1] ? 1 : rb[i][0] - lb[i][1] + 2;
        double max_range = rb[i][1] - lb[i][0];

        if (mean <= min_range) tmp.push_back(min_range);
        else if (mean >= max_range) tmp.push_back(max_range);
        tmp.push_back(mean);
    }
    double min_var = variance(tmp, mean);
    double min_sigma = sqrt(min_var);

    restored_values.clear();
    for (size_t i=0; i<N; i++)
    {
        double min_range = rb[i][0] < lb[i][1] ? 1 : rb[i][0] - lb[i][1] + 2;
        double max_range = rb[i][1] - lb[i][0];

        double center = 0.5 * (lb[i][1] - 1 + rb[i][0]);
        interval<double> range(center - 0.5 * mean, center + 0.5 * mean);
            // разместим по центру между левым и правым интервалом, что гарантирует
            // попадание в оба крайних граничных интервала

        if (mean <= min_range) range = interval<double>(lb[i][1]-1, rb[i][0]+1); // ") [" -> "[)"
        if (mean >= max_range) range = interval<double>(lb[i][0], rb[i][1]);
            // сделаем поправку размещения, если натыкаемся на границы

        restored_values.push_back(range);
    }
    return std::make_pair(mean, min_sigma);
}

std::pair<double,double> intrainterval_localization_by_MLM(
    std::vector<interval<double>> &restored_values,
    const std::vector<interval<double>> &left_bounds,
    const std::vector<interval<double>> &right_bounds,
    size_t num_iterations, bool print_flag)
{
    //----------------------------------------------------------------------------------------------
    //              задание функций для метода максимального правдоподобия
    //----------------------------------------------------------------------------------------------
    struct MLMab
    {
        const std::vector<interval<double>> *ta; // левые границы интервалов
        const std::vector<interval<double>> *tb; // правые границы интервалов
        double reg_coef; // коэффициент регуляризации
        double *mean, *sigma;
        size_t N;

        double f(const double p[])
        {
            assert(ta->size() == tb->size() && N == ta->size());
            double Gs = 1. / sqr(*sigma);
            double sum = 0;
            for (size_t i=0; i<N; i++)
            {
                double tta = p[2*i];
                double ttb = p[2*i + 1];
                sum += Gs * sqr(ttb - tta - *mean);
                    // основной член

                double offa1 = tta < (*ta)[i][0] ? tta - (*ta)[i][0] : 0;
                double offa2 = tta > (*ta)[i][1] ? tta - (*ta)[i][1] : 0;
                double offb1 = ttb < (*tb)[i][0] ? ttb - (*tb)[i][0] : 0;
                double offb2 = ttb > (*tb)[i][1] ? ttb - (*tb)[i][1] : 0;
                sum += reg_coef * (sqr(offa1) + sqr(offa2) + sqr(offb1) + sqr(offb2));
                    // регуляционный член
            }
            return N * log(*sigma) + 0.5 * sum;
        }

        double fdf(double g[], const double p[])
        {
            double Gs = 1. / sqr(*sigma);

            for (size_t i=0; i<2*N; i++) g[i] = 0;
            double sum = 0;
            for (size_t i=0; i<N; i++)
            {
                double tta = p[2*i];
                double ttb = p[2*i + 1];

                double f = ttb - tta - *mean;
                sum += Gs * sqr(f);

                g[2*i] -= Gs * f;
                g[2*i+1] += Gs * f;

                double offa1 = tta < (*ta)[i][0] ? tta - (*ta)[i][0] : 0;
                double offa2 = tta > (*ta)[i][1] ? tta - (*ta)[i][1] : 0;
                double offb1 = ttb < (*tb)[i][0] ? ttb - (*tb)[i][0] : 0;
                double offb2 = ttb > (*tb)[i][1] ? ttb - (*tb)[i][1] : 0;
                sum += reg_coef * (sqr(offa1) + sqr(offa2) + sqr(offb1) + sqr(offb2));

                g[2*i] += reg_coef * (offa1 + offa2);
                g[2*i+1] += reg_coef * (offb1 + offb2);
            }
            return N * log(*sigma) + 0.5 * sum;
        }
    } mlmab;

    size_t N = left_bounds.size();
    std::vector<double> x(2*N+2);

    mlmab.ta = &left_bounds;
    mlmab.tb = &right_bounds;
    mlmab.mean = &x[2*N];
    mlmab.sigma = &x[2*N+1];
    mlmab.reg_coef = 1000;
    mlmab.N = N;

    //----------------------------------------------------------------------------------------------
    //              вычисление стартовых значений оптимизируемых параметров
    //----------------------------------------------------------------------------------------------
    double sum = 0;
    for (size_t i=0; i<N; i++)
    {
        x[2*i  ] = 0.5 * (left_bounds[i][1] + left_bounds[i][0]);
        x[2*i+1] = 0.5 * (right_bounds[i][1] + right_bounds[i][0]);
        sum += x[2*i+1] - x[2*i];
    }

    double mean = sum / N;
    double sigma = 0;
    for (size_t i=0; i<N; i++) sigma += sqr(x[2*i+1] - x[2*i] - mean);
    sigma = sqrt(sigma / N);

    x[2*N] = mean;
    x[2*N+1] = sigma;

    //---------------------------------------------------------------------------------------------
    //                              оптимизация параметров
    //---------------------------------------------------------------------------------------------
    prg::Minimizer<MLMab> minimizerab(mlmab);
    auto &param = minimizerab.get_params();
    param.algorithm = prg::Minimizer<MLMab>::bfgs2;
    param.print_flag = -1; // строгий запрет на печать

    for (size_t iteration=0; iteration<num_iterations; iteration++)
    {
        std::vector<double> values;
        for (size_t i=0; i<N; i++) values.push_back(x[2*i+1] - x[2*i]);

        moment_t moment = statistics(values);
        x[2*N] = moment[0]; x[2*N+1] = moment[1];

        auto result = minimizerab.run(40, &x[0], &x[2*N]);
        if (print_flag) printf("Iteration %4ld L=%.9e mean = %.3lf sigma=%.3lf \n",
            iteration, result.first, x[2*N], x[2*N+1]);
            std::cout.flush();
        if (result.second) break;
    }

    restored_values.clear();
    for (size_t i=0; i<N; i++)
        restored_values.push_back(interval<double>(x[2*i], x[2*i+1]));

    return std::make_pair(x[2*N], x[2*N+1]);
}

/*!
 * \brief восстановление "точных" значений границ фаз методом максимального правдоподобия
 * \param ne число экземляров во множестве экземпляров
 * \param ni число интервалов в каждом экземпляре
 * \param ps указатель на массив всех интервалов
 * \param num_iterations число итераций поиска
 * \param print_flag флаг печати результатов
 * \return
 */
std::vector<double> intrainterval_localization_by_MLM(Unsigned ne, Unsigned ni,
    const interval<double> *ps, Unsigned num_iterations, bool print_flag)
{
    Unsigned np = ne * ni; // полное число оптимизирумых параметров

    //----------------------------------------------------------------------------------------------
    //              задание функций для метода максимального правдоподобия
    //----------------------------------------------------------------------------------------------
    struct MLM
    {
        Unsigned ne; // число экземпляров
        const interval<double> *ps; // указатель на 1-й элемент множества экземпляров
        Unsigned num_phase; // число фаз (едино для всех экземпляров) = ni -1
        Unsigned ni; //
        Unsigned np; // общее число переменных оптимизации = ne * ni
        double reg_coef; // коэффициент регуляризации

        std::vector<double> mean, sigma2; // [num_phase]
            // значения mu & sigma не оптимизируются, а вычисляются для каждой итерации

        /// вычисление средних значений длин всех фаз и дисперсий = sqr(sigma)
        void averages_and_dispersions(const double p[])
        {
            for (Unsigned j=0; j<num_phase; j++)
            {
                mean[j] = 0;
                for (Unsigned i=j; i<np; i+=ni) { mean[j] += p[i+1] - p[i]; }
                mean[j] /= ne;

                sigma2[j] = 0;
                for (Unsigned i=j; i<np; i+=ni) sigma2[j] += sqr(p[i+1] - p[i] - mean[j]);
                sigma2[j] /= ne;
            }
        }

        double f(const double p[])
        {
            averages_and_dispersions(p);

            double L = 0;
            for (Unsigned j=0; j<num_phase; j++) L += 1 + log(sigma2[j]);
            L *= 0.5 * ne;

            double R = 0;
            for (Unsigned i=0; i<np; i++)
            {
                double r0 = p[i] < ps[i][0] ? p[i] - ps[i][0] : 0;
                double r1 = p[i] > ps[i][1] ? p[i] - ps[i][1] : 0;
                R += sqr(r0) + sqr(r1);
            }
            R *= reg_coef;

            return L + 0.5 * R;
        }

        double fdf(double g[], const double p[])
        {
            averages_and_dispersions(p);
            //printf("%lf %lf %lf %lf  \n", mean[0], sqrt(sigma2[0]), mean[1], sqrt(sigma2[1]));

            double L = 0;
            for (Unsigned j=0; j<num_phase; j++) L += 1 + log(sigma2[j]);
            L *= 0.5 * ne;

            double R = 0;
            for (Unsigned i=0; i<np; i++)
            {
                double r0 = p[i] < ps[i][0] ? p[i] - ps[i][0] : 0;
                double r1 = p[i] > ps[i][1] ? p[i] - ps[i][1] : 0;
                R += sqr(r0) + sqr(r1);
                g[i] = reg_coef * (r0 + r1);
            }
            R *= reg_coef;

            double gprev = 0;
            for (Unsigned i=0; i<np; i++)
            {
                Unsigned j = i % ni;
                double gcur = (i+1) % ni == 0 ? 0 : (p[i+1] - p[i] - mean[j]) / sigma2[j];
                g[i] += gprev - gcur;
                gprev = gcur;
            }

            return L + 0.5 * R;
        }
    } mlm;

    mlm.ne = ne;
    mlm.ni = ni;
    mlm.ps = ps;
    mlm.num_phase = ni - 1;
    mlm.mean.resize(ni - 1, 0);
    mlm.sigma2.resize(ni - 1, 0);
    mlm.reg_coef = 100000;
    mlm.np = np;

    std::vector<double> x(np, 0.);

    //----------------------------------------------------------------------------------------------
    //              вычисление стартовых значений оптимизируемых параметров
    //----------------------------------------------------------------------------------------------

    for (size_t i=0; i<ne; i++)
    {
        x[3*i  ] = ps[3*i  ][0];
        x[3*i+1] = 0.5 * (ps[3*i+1][0] + ps[3*i+1][1]);
        x[3*i+2] = ps[3*i+2][1];


        if (print_flag)
        {
            Timestamp ta1 = ps[3*i][0];
            Timestamp ta2 = ps[3*i][1];

            Timestamp tb1 = ps[3*i+1][0];
            Timestamp tb2 = ps[3*i+1][1];

            Timestamp tc1 = ps[3*i+2][0];
            Timestamp tc2 = ps[3*i+2][1];

            Timestamp t1 = x[3*i    ];
            Timestamp t2 = x[3*i + 1];
            Timestamp t3 = x[3*i + 2];

            if (print_flag)
            {
                printf("   : [%6s %6s %6s] [%6s %6s %6s] [%6s %6s %6s] \n",
                    ta1.to_string().c_str(), t1.to_string().c_str(), ta2.to_string().c_str(),
                    tb1.to_string().c_str(), t2.to_string().c_str(), tb2.to_string().c_str(),
                    tc1.to_string().c_str(), t3.to_string().c_str(), tc2.to_string().c_str()
                );
                FLUSH;
            }
        }
    }

    //---------------------------------------------------------------------------------------------
    //                              оптимизация параметров
    //---------------------------------------------------------------------------------------------
    prg::Minimizer<MLM> minimizer(mlm);
    auto &param = minimizer.get_params();
    param.algorithm = prg::Minimizer<MLM>::bfgs;
    //param.algorithm = prg::Minimizer<MLM>::conjugate_pr;
    param.print_flag = -1;
    minimizer.run(num_iterations, &x[0], &x[np]);

    if (print_flag)
    {
        mlm.averages_and_dispersions(&x[0]);

        std::printf(" %8s %8s %8s\n", "phase", "mean", "sigma");
        for (Unsigned i=0; i<mlm.num_phase; i++)
            std::printf(" %8ld %8.5lf %8.5lf\n", i, mlm.mean[i], sqrt(mlm.sigma2[i]));
        std::cout.flush();
    }

    return x;
}
