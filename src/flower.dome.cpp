#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "prg.os.h"
#include "prg.string.h"
#include "prg.random.h"

#include "test.defs.h"
#include "math.defs.h"
#include "math.statistics.h"
#include "flower.dome.h"
#include "data.recovery.h"

// извлечение числа колонок и числа дат в файле csv
std::pair<unsigned, unsigned> file_columns(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file) no_file_exit(filename);

    char buffer[BUF_LEN], buf[BUF_LEN];

    // Узнаем число колонок в файле num_columns и число колонок дат num_dates подсчетом. Число
    // колонок для идентификатора образца num_names меняется от 1 до 3, они зачастую не имеют имени.
    // Данные в них отсутствуют, что означает "они такие же как предыдущие".
    // num_names = num_columns - num_dates

    unsigned num_columns = 0, num_dates = 0;
    file.getline(buffer, BUF_LEN);
    std::stringstream ss(buffer);
    ss.getline(buf, BUF_LEN, ',');
    while (ss)
    {
        const char *p = prg::trim(buf, remove_symbols);
        if (is_date(p)) num_dates++;
        num_columns++;
        ss.getline(buf, BUF_LEN, ',');
    }
    return std::make_pair(num_columns - num_dates, num_dates);
}

std::vector<flowerline2_t> Flowerdome::extract_phase_bounds(interval<Flowerstamp> fs,
    phase_enum_t phase, bool print_flag) const
{
    std::vector<flowerline2_t> result;

    //----------------------------------------------------------------------------------------------
    //            фильтрируем записи только для заданного интервала номеров цветков
    //----------------------------------------------------------------------------------------------
    auto first = std::lower_bound(begin(), end(), fs[0],
        [](const flower_t &f, const Flowerstamp &fs) { return f.fsp < fs; }
    );
    if (first == end()) return result; // возврат пустого массива

    auto last = std::upper_bound(begin(), end(), fs[1],
        [](const Flowerstamp &fs, const flower_t &f) { return fs < f.fsp; }
    );
    if (first == last) return result; // возврат пустого массива

    //----------------------------------------------------------------------------------------------
    //                  сделаем оценку максимальной продолжительности фазы
    //----------------------------------------------------------------------------------------------
    int max_time = 0, min_time = 100;
        // максимальный и минимальный интервалы времени, который выставляются если нет фазы
        // (до и после заданной и между пропущенной)
        // это оценочное значение, сделанное только по тем цветкам, где фаза зарегистрирована

    auto it = first;
    while (true)
    {
        //------------------------------------------------------
        // поиск первой записи, где появляется заданная фаза
        //------------------------------------------------------
        while (it != last && it->psp != phase) ++it;
        if (it == last) break;

        Timestamp ta2 = it->tsp;
        Flowerstamp fsa = it->fsp;

        //------------------------------------------------------
        // поиск последней записи, где сохраняется заданная фаза
        //------------------------------------------------------
        while (it != last && it->psp == phase && it->fsp == fsa) ++it;
        Timestamp tb1 = (it-1)->tsp;

        int time = tb1 - ta2 + 1;
        if (max_time < time) max_time = time;
        if (min_time > time) min_time = time;
    }
    if (print_flag)
    {
        printf("limits of phase \'%s\' = [%d, %d] days\n", Phasestamp(phase).to_string().c_str(),
            min_time, max_time); FLUSH
    }

    while (first != last)
    {
        Flowerstamp current_flower = first->fsp;

        auto right = first;
        while (right != last && right->fsp == current_flower) ++right;
            // находим границы записей [first, right), принадлежащх одному цветку

        while (first != right && first->psp < phase) ++first;
            // поиск первой записи, где появляется заданная (или следующая) фаза

        if (first != last && first->psp != DEATH)
        {
            Timestamp ta1, ta2, tb1, tb2;
            if (first->psp == phase) // фаза найдена
            {
                ta2 = first->tsp;
                std::pair<Timestamp, bool> prev = prev_date(ta2);
                ++ta2; // внеинтервальность после поиска prev (иначе prev укажет на ta2)
                ta1 = prev.second ? prev.first + 1 : ta2 - max_time;

                while (first != right && first->psp == phase) ++first;
                    // поиск записи, где фаза переключается

                tb1 = (first-1)->tsp; // внутриинтервальность
                std::pair<Timestamp, bool> next = next_date(tb1);
                tb2 = next.second ? next.first : tb1 + max_time;
                    // next.first без -1 для внеинтервальности

            }
            else // фаза отсутствует в данных
            {
                tb2 = first->tsp; // можно сразу внеинтервальность, так как фазы нет
                std::pair<Timestamp, bool> prev = prev_date(tb2);
                ta1 = prev.second ? prev.first + 1 : tb2 - max_time;

                int delta = std::min(min_time, tb2 - ta1) - 1;
                ta2 = tb2 - delta;
                tb1 = ta1 + delta;
            }

            if (print_flag)
            {
                printf("%s : [%s, %s] - [%s %s]\n", current_flower.to_string().c_str(),
                    ta1.to_string().c_str(), (ta2-1).to_string().c_str(),
                    tb1.to_string().c_str(), (tb2-1).to_string().c_str()
                ); // вывод делаем, включая границу
                FLUSH
            }

            interval<Timestamp> ta(ta1, ta2);
            interval<Timestamp> tb(tb1, tb2);
            result.push_back(flowerline2_t{current_flower, {ta, tb}});
        }

        while (first != last && first->fsp == current_flower) ++first;
            // перескочим на другой цветок для поиска нужной фазы
    }
    return result;
}

moment_t Flowerdome::control_statistics(interval<Flowerstamp> fs,
    phase_enum_t phase, bool print_flag) const
{
    auto intervals = extract_phase_bounds(fs, phase, print_flag);
    std::vector<double> days; // число дней фазы для каждого цветка

    for (auto e : intervals)
    {
        interval<Timestamp> ta = e[0];
        interval<Timestamp> tb = e[1];
        double a = 0.5 * (ta[0].to_double() + ta[1].to_double());
        double b = 0.5 * (tb[0].to_double() + tb[1].to_double());

        double duration = b - a;
        days.push_back(duration);
    }
    moment_t moments = statistics(days);
    double mean = moments[0];
    double sigma = moments[1];

    if (print_flag)
    {
        printf("control of phase \'%s\': mean=%4.2lf sigma=%4.2lf\n",
            Phasestamp(phase).to_string().c_str(), mean, sigma);

        size_t errors = 0;
        size_t uncorrectable_errors = 0;

        double N = double(days.size());
        double criteria2 = -2 * sqr(sigma) * log(sqrt(2*M_PI) * sigma / N);
            // критерий маловероятных объектов

        for (auto e : intervals)
        {
            Flowerstamp flower = e.fsp;
            interval<Timestamp> ta = e[0];
            interval<Timestamp> tb = e[1];
            double a = 0.5 * (ta[0].to_double() + ta[1].to_double());
            double b = 0.5 * (tb[0].to_double() + tb[1].to_double());

            double t = b - a;
            double ac = 0, bc = 0;
            if (t > mean)
            {
                ac = ta[1].to_double();
                bc = tb[0].to_double();
            }
            else if (t < mean)
            {
                ac = ta[0].to_double();
                bc = tb[1].to_double();
            }

            double s2 = sqr(t - mean);
            if (s2 > criteria2)
            {
                errors++;
                double p = exp(-s2 / (2 * sqr(sigma))) / (sqrt(2*M_PI) * sigma);

                printf("\twarning: \'%6s\' doubtful duration=%2d with probability=%12le ",
                    flower.to_string().c_str(), int(t), p);

                s2 = sqr(bc - ac - mean);
                if (s2 > criteria2) { printf("- UNCORRECTABLE\n"); uncorrectable_errors++; }
                else printf("\n");
            }
        }
        if (errors) { printf("\tuncorrectable errors=%ld all errors=%ld\n",
            uncorrectable_errors, errors); FLUSH }
    }
    return moments;
}

bool Flowerdome::load_file(const std::string &filename, const std::vector<unsigned> &sprouts,
    bool print_flag)
{
    if (filename.empty()) return false;
    species_name = prg::basename(filename + ".csv"); // сохраним имя вида

    bool valid_flag = true; // сигнал отсутствия ошибок в файле
    char buffer[BUF_LEN], buf[BUF_LEN];

    std::pair<size_t, size_t> n = file_columns(filename + ".csv");
    size_t num_names = n.first, num_dates = n.second;

    std::ifstream file(filename + ".csv");
    if (!file) no_file_exit(filename);

    file.getline(buffer, BUF_LEN); // разбор имен колонок для дат
    {
        std::stringstream ss(buffer);

        for (size_t i=0; i<num_names; i++) ss.getline(buf, BUF_LEN, ',');
            // удалим из потока первые num_names значений

        Timestamp prev(0);
        for (size_t i=0; i<num_dates; i++)
        {
            ss.getline(buf, BUF_LEN, ',');
            const char *p = prg::trim(buf, remove_symbols);
            bool flag = is_date(p);
            if (!flag)
            {
                printf("error: invalid date \'%s\'\n", buf);
                valid_flag = false;
            }

            Timestamp cur(p);
            if (cur < prev)
            {
                printf("error: wrong date order \'%s\'\n", buf);
                valid_flag = false;
            }

            dates.push_back(cur);
            prev = cur;
        }
    }

    flower_t record; // текущая загружаемая запись
    Flowerstamp fs_prev, fs_cur;

    std::set<unsigned> sproutset;
    for (auto e : sprouts) sproutset.insert(e);

    file.getline(buffer, BUF_LEN);
    while (!file.eof())
    {
        std::stringstream ss(buffer);
        for (size_t i=0; i<num_names; i++)
        {
            ss.getline(buf, BUF_LEN, ',');
            char *p = prg::trim(buf, remove_symbols);
            fs_cur.set(unsigned(i), unsigned(atoi(p)));
        }
        for (size_t i=0; i<num_names; i++) // наложили на пустые места прежние имена
        {
            if (fs_cur.get(unsigned(i)) != 0) break;
            fs_cur.set(unsigned(i), fs_prev.get(unsigned(i)));
        }

        if (!sproutset.empty() && sproutset.count(fs_cur.sprout()) == 0)
        {
            fs_prev = fs_cur;
            file.getline(buffer, BUF_LEN);
            continue; // игнорируем побеги, которых нет в списке
        }

        record.fsp = fs_cur;
        flowers.insert(fs_cur); // собираем все цветки файла

        Phasestamp prev_psp(UNDEF), cur_psp(UNDEF);
        for (size_t i=0; i<num_dates; i++)
        {
            ss.getline(buf, BUF_LEN, ',');
            char *p = prg::trim(buf, remove_symbols);
            cur_psp = make_phasestamp(prev_psp, std::string(p));

            if (cur_psp.phase == INVALID)
            {
                printf("error: invalid phase \'%s\' at (%s, %s)\n", p,
                    fs_cur.to_string().c_str(),
                    dates[i].to_string().c_str());
                valid_flag = false;
            }

            record.psp = cur_psp;
            record.tsp = dates[i];

            if (!record.psp.empty()) // не вставляем не выросшие образцы в начале лета
            {
                push_back(record);
                _DEBUG_PRINT(record);
            }
            prev_psp = cur_psp;
        }
        if (print_flag && (cur_psp == Phasestamp(BUD) || cur_psp == Phasestamp(FLOWER)))
        {
            printf("warning: unclosed record \'%s\'\n",
                fs_cur.to_string().c_str());
        }

        fs_prev = fs_cur;
        file.getline(buffer, BUF_LEN);
    }

    // обеспечение гарантии упорядочения записей если нет ошибок в файле
    if (valid_flag) std::sort(begin(), end(),
        [](const flower_t &a, const flower_t &b)
        {
            if (a.fsp < b.fsp) return true;
            if (b.fsp < a.fsp) return false;
            return a.tsp < b.tsp;
        });

    //---------------------------------------------------------------------------------------------
    //         контроль выхода данных за пределы доверительных интервалов
    //---------------------------------------------------------------------------------------------
    if (print_flag)
    {
        control_statistics(interval<Flowerstamp>(0, 10), BUD, print_flag);
        control_statistics(interval<Flowerstamp>(0, 10), FLOWER, print_flag);
        //control_statistics(OVARY, print_flag);
        //control_statistics(SEED, print_flag);
    }

    if (!valid_flag) printf("file \'%s\' contains errors\n", filename.c_str());
    //printf("%ld items are loaded from file \'%s\'\n", size(), filename.c_str());


    //---------------------------------------------------------------------------------------------
    //         сосчитаем число цветков на осях 1-го и 2-го порядков
    //---------------------------------------------------------------------------------------------
    for (auto e : flowers)
    {
        unsigned sprout = e.sprout();
        unsigned axis = e.axis_order();

        if (flowers_count.count(sprout) == 0)
            flowers_count[sprout] = std::array<unsigned, 3>({{0, 0, 0}});

        flowers_count[sprout][axis]++;
    }


    print_information();
    return valid_flag;
}

void Flowerdome::print_information() const
{
    PRINT_SEPARATOR_LINE(80, '-');
    printf("%s : \n", species_name.c_str());
    for (auto e : flowers_count)
        printf("sprout %d : #flowers on axes {%d %d}\n", e.first, e.second[1], e.second[2]);
    PRINT_SEPARATOR_LINE(80, '-');
}

std::vector<flowerline5_t> Flowerdome::extract_phases_intervals(interval<Flowerstamp> fs,
    bool print_flag) const
{
    std::vector<flowerline5_t> result;

    //----------------------------------------------------------------------------------------------
    //            фильтрируем записи только для заданного интервала номеров цветков
    //----------------------------------------------------------------------------------------------
    auto first = std::lower_bound(begin(), end(), fs[0],
        [](const flower_t &f, const Flowerstamp &fs) { return f.fsp < fs; }
    );
    if (first == end()) return result; // возврат пустого массива

    auto last = std::upper_bound(begin(), end(), fs[1],
        [](const Flowerstamp &fs, const flower_t &f) { return fs < f.fsp; }
    );
    if (first == last) return result; // возврат пустого массива

    auto it = first;
    while (it != last)
    {
        Flowerstamp current_flower = it->fsp;

        auto right = it;
        while (right != last && right->fsp == current_flower) ++right;
            // находим границы записей [left, right), принадлежащх одному цветку

        decltype(it) start[5] = {right, right, right, right, right}; // bud, flower, ovary, seed, death
        decltype(it) end[5] = {right, right, right};

        while (it != right)
        {
            if (start[0] == right && it->psp == BUD)
            {
                start[0] = it;
                while (it != right && it->psp == BUD) ++it;
                --it; end[0] = it;
            }
            if (start[1] == right && it->psp == FLOWER)
            {
                start[1] = it;
                while (it != right && it->psp == FLOWER) ++it;
                --it; end[1] = it;
            }
            if (start[2] == right && it->psp == OVARY)
            {
                start[2] = it;
                while (it != right && it->psp == OVARY) ++it;
                --it; end[2] = it;
            }
            if (start[3] == right && it->psp == SEED)
            {
                start[3] = it;
                while (it != right && it->psp == SEED) ++it;
                --it; end[3] = it;
            }
            if (start[4] == right && it->psp == DEATH)
            {
                start[4] = it;
                while (it != right && it->psp == DEATH) ++it;
                --it; end[4] = it;
            }
            ++it;
        }

        interval<Timestamp> range[5] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
        if (start[0] != right) range[0] = {start[0]->tsp, end[0]->tsp};
        if (start[1] != right) range[1] = {start[1]->tsp, end[1]->tsp};
        if (start[2] != right) range[2] = {start[2]->tsp, end[2]->tsp};
        if (start[3] != right) range[3] = {start[3]->tsp, end[3]->tsp};
        if (start[4] != right) range[4] = {start[4]->tsp, end[4]->tsp};

        flowerline5_t flowerline{current_flower, {range[0], range[1], range[2], range[3], range[4]}};
        result.push_back(flowerline);
        if (print_flag) { print(flowerline); FLUSH; }

        while (it != right) ++it;
    }
    return result;
}

std::vector<flowerline3_t> Flowerdome::extract_phase_bounds(interval<Flowerstamp> fs,
    bool skip_death, bool print_flag) const
{
    std::vector<flowerline5_t> phases = extract_phases_intervals(fs, false);

    //----------------------------------------------------------------------------------------------
    //              сделаем оценку максимальной продолжительности фаз (BUD, FLOWER)
    //----------------------------------------------------------------------------------------------
    int mintime[2] = {100, 100}, maxtime[2] = {0, 0};
        // минимальный и максимальный интервалы времени, который выставляются если фазы
        // зарегистрированы в неполном объеме (до BUD и после FLOWER)
        // это оценочное значение, сделанное только по тем цветкам, где обе фазы зарегистрирована
        // с обеими концами

    for (auto e : phases)
    {
        Flowerstamp current_flower = e.fsp;
        if (print_flag) { printf("%s : ", current_flower.to_string().c_str()); FLUSH; }

        for (auto iphase : {0, 1})
        {
            auto range = e.intervals[iphase];
            auto t2 = range[0];
            auto t3 = range[1];

            if (t2.valid())
            {
                int dur = t3.to_int() - t2.to_int() + 1;
                assert(dur > 0);
                if (mintime[iphase] > dur) mintime[iphase] = dur;
                if (maxtime[iphase] < dur) maxtime[iphase] = dur;

                if (print_flag)
                {
                    printf(" %c[%s %s] min=%d ", (iphase == 0 ? 'B' : 'F'),
                        t2.to_string().c_str(), t3.to_string().c_str(),dur
                    );
                }
            }
        }
        if (print_flag) printf("\n"); FLUSH;
    }
    if (print_flag)
    {
        printf("limits of phase \'bud   \' = [%d, %d] days \n", mintime[0], maxtime[0]);
        printf("limits of phase \'flower\' = [%d, %d] days \n", mintime[1], maxtime[1]);
        FLUSH;
    }

    std::vector<flowerline3_t> result;
    typedef interval<Timestamp> interval_type;

    interval_type bound[3] = {{0,0}, {0,0}, {0,0}};
    for (auto e : phases)
    {
        Flowerstamp current_flower = e.fsp;

        unsigned key = 0; // битовое число наличия фаз ([4] DEATH [2] OVARY, [1] FLOWER, [0] BUD)
        for (auto iphase : {0, 1, 2})
        {
            auto range = e.intervals[iphase];
            if (range[0].valid()) key += 1 << iphase;
        }
        if (!(key & 0x2) && skip_death && e.intervals[4][0].valid()) continue;
            // если нет FLOWER (бит 1) и задано удаление при DEATH, то эти записи не учитываются

        interval_type bud = e.intervals[0];
        interval_type flower = e.intervals[1];
        interval_type ovary = e.intervals[2];

        //---------------------------------------------------------------------------
        // макросы для установления границ:
        //---------------------------------------------------------------------------
        // (t1.. - гипотеза о возможном начале фазы
        // [t2.. - известное начало фазы
        // ..t3] - известный конец фазы
        // ..t4) - гипотеза о возможном конце фазы
        //---------------------------------------------------------------------------
        #define SET_LEFT(t1, t2, t3, t4, phaselen) { \
            prev = prev_date(t2); \
            t1 = prev.second ? prev.first + 1 : std::min(t2, t3 - phaselen); \
        }
        #define SET_RIGHT(t1, t2, t3, t4, phaselen) { \
            next = next_date(t3); \
            t4 = next.second ? next.first - 1 : std::max(t3, t2 + phaselen); \
        }
        //---------------------------------------------------------------------------

        Timestamp t1, t2, t3, t4; // (1, [2, 3], 4)
        std::pair<Timestamp, bool> prev, next;

        if (key == 0) continue; // пропускаем запись с полностью неопределенными фазами
        switch (key)
        {
        case 1 : // BUD
            t2 = bud[0]; t3 = bud[1];
            SET_LEFT(t1, t2, t3, t4, maxtime[0]);
            bound[0] = interval_type(t1, t2);
            //------------------------- прогноз -------------------------
            SET_RIGHT(t1, t2, t3, t4, maxtime[0]);
            t3 = t3 + 1;
            t4 = t4 + 1;
            bound[1] = interval_type(t3, t4);
            bound[2] = interval_type(t3 + mintime[1], t4 + maxtime[1]);
            break;

        case 2 : // FLOWER
            t2 = flower[0]; t3 = flower[1];
            SET_LEFT(t1, t2, t3, t4, maxtime[1]);
            bound[1] = interval_type(t1, t2);
            //------------------------- прогноз -------------------------
            bound[0] = interval_type(t1 - maxtime[0], t2);
            bound[2] = interval_type(t3 + 1, t2 + maxtime[1]);
            break;

        case 3 : // BUD + FLOWER
            t2 = bud[0]; t3 = bud[1];
            SET_LEFT(t1, t2, t3, t4, maxtime[0]);
            bound[0] = interval_type(t1, t2);

            bound[1] = interval_type(bud[1]+1, flower[0]);

            t2 = flower[0]; t3 = flower[1];
            SET_RIGHT(t1, t2, t3, t4, maxtime[1]);
            bound[2] = interval_type(t3 + 1, t4 + 1);
            break;

        case 4 : // OVARY
            t2 = ovary[0];
            bound[2] = interval_type(t2, t2); // полагаем это точный конец
            //------------------------- прогноз -------------------------
            {
//                std::pair<Timestamp, bool> prev = prev_date(t2);
//                if (prev.second)
//                {
//                    t1 = prev.first + 1;
//                    bound[1] = interval_type(t1 + 1, t2);
//                    bound[0] = interval_type(t1, t2 - 1);
//                }
//                else
//           дает ошибку длины фазы = 0

                {
                    t1 = t2 - maxtime[1]; t2 = t2 - mintime[1];
                    bound[1] = interval_type(t1, t2);
                    bound[0] = interval_type(t1 - maxtime[0], t2 - mintime[0]);
                }
            }
            break;

        case 5 : // BUD + OVARY
            t2 = bud[0]; t3 = bud[1];
            SET_LEFT(t1, t2, t3, t4, maxtime[0]);
            bound[0] = interval_type(t1, t2);

            t4 = ovary[0];
            bound[1] = interval_type(t3 + 1, t4 - 1);
            bound[2] = interval_type(t3 + 2, t4);
            break;

        case 6 : // FLOWER + OVARY
            bound[2] = interval_type(flower[1] + 1, ovary[0]);

            t2 = flower[0]; t3 = flower[1];
            SET_LEFT(t1, t2, t3, t4, maxtime[1]);
            //------------------------- прогноз -------------------------
            if (!prev.second)
            {
                t1 = std::min(t2, t3 - maxtime[1]);
                bound[1] = interval_type(t1, t2);
                bound[0] = interval_type(t1-1-maxtime[0], t2-1);
            }
            else
            {
                bound[0] = interval_type(t1, t2-1);
                t1 = std::min(t1 + mintime[0], t2);
                bound[1] = interval_type(t1, t2);
            }
            break;

        case 7 : // BUD + FLOWER + OVARY
            t2 = bud[0]; t3 = bud[1];
            SET_LEFT(t1, t2, t3, t4, maxtime[0]);
            bound[0] = interval_type(t1, t2);
            bound[1] = interval_type(t3 + 1, flower[0]);
            bound[2] = interval_type(flower[1] + 1, ovary[0]);
            break;
        }
        #undef SET_LEFT
        #undef SET_RIGHT

        if (!bound[0][0].valid()) bound[0] = {1, bound[0][1]};

        flowerline3_t flowerline{current_flower, {bound[0], bound[1], bound[2]}};
        if (print_flag) { print(flowerline); FLUSH; }

        result.push_back(flowerline);
    }
    return result;
}

std::vector<flowerline3_t> Flowerdome::interval_clipping(unsigned num_iterations,
    interval<Flowerstamp> fs, bool skip_death, bool print_flag) const
{
    std::vector<flowerline3_t> bounds = extract_phase_bounds(fs, skip_death, print_flag);
    if (bounds.empty()) return bounds;

    //---------------------------------------------------------------------------------------------
    // конвертируем данные к формату функции intrainterval_localization_by_MLM()
    //---------------------------------------------------------------------------------------------
    Unsigned ne = bounds.size();

    std::vector<interval<double>> intervals;
    for (Unsigned i=0; i<ne; i++)
    {
        for (Unsigned j=0; j<3; j++)
        {
            interval<Timestamp> range = bounds[i].intervals[j];
            intervals.push_back(interval<double>{range[0].to_double(), range[1].to_double()});
        }
    }

    //---------------------------------------------------------------------------------------------
    //                          рассчитаем суженные интервалы
    //---------------------------------------------------------------------------------------------

    auto restored_values = intrainterval_localization_by_MLM(ne, 3,
        &intervals[0], num_iterations, print_flag);

    //---------------------------------------------------------------------------------------------
    //                  вернем данные к формату возвращаемого значения
    //---------------------------------------------------------------------------------------------
    for (Unsigned i=0; i<ne; i++)
    {
        Flowerstamp flower = bounds[i].fsp;

        if (print_flag) printf("%8s : ", flower.to_string().c_str());

        auto r0 = bounds[i].intervals[0];
        Timestamp ta1 = r0[0];
        Timestamp ta2 = r0[1];

        auto r1 = bounds[i].intervals[1];
        Timestamp tb1 = r1[0];
        Timestamp tb2 = r1[1];

        auto r2 = bounds[i].intervals[2];
        Timestamp tc1 = r2[0];
        Timestamp tc2 = r2[1];

        Timestamp t1 = restored_values[3*i    ];
        Timestamp t2 = restored_values[3*i + 1];
        Timestamp t3 = restored_values[3*i + 2];

        int days[2] = { t2.to_int() - t1.to_int(), t3.to_int() - t2.to_int() };

        if (print_flag)
        {
            printf("%8s : [%6s %6s %6s] [%6s %6s %6s] [%6s %6s %6s] BUD=%d FLOWER=%d \n",
                flower.to_string().c_str(),
                ta1.to_string().c_str(), t1.to_string().c_str(), ta2.to_string().c_str(),
                tb1.to_string().c_str(), t2.to_string().c_str(), tb2.to_string().c_str(),
                tc1.to_string().c_str(), t3.to_string().c_str(), tc2.to_string().c_str(),
                days[0], days[1]
            );
            FLUSH;
        }
        assert(t1.to_int() <= t2.to_int() && t2.to_int() <= t3.to_int());
        assert(ta1.to_int() <= t1.to_int() && t1.to_int() <= ta2.to_int());
        assert(tb1.to_int() <= t2.to_int() && t2.to_int() <= tb2.to_int());
        assert(tc1.to_int() <= t3.to_int() && t3.to_int() <= tc2.to_int());
    }

    std::vector<flowerline5_t> phases = extract_phases_intervals(fs, print_flag);

    unsigned indent = 25; // число выдаваемых дней до и после дат файла

    unsigned startpos = dates.front().to_unsigned() - indent; // стартовая дата файла (со смещением)
    if (int(startpos) <= 0) startpos = 0;

    unsigned lastpos = dates.back().to_unsigned() + indent; // последняя дата файла (со смещением)
    unsigned len = lastpos + 1; // полное число отрисовываемых дней

    std::set<unsigned> fixdates; // все даты файла
    for (auto e : dates) fixdates.insert(e.to_unsigned());

    for (Unsigned i=0; i<ne; i++)
    {
        Flowerstamp flower = phases[i].fsp;
        std::string draw_phases(len, ' ');

        auto r = phases[i].intervals[0];
        Timestamp t1 = r[0], t2 = r[1], t3;
        if (t1.valid())
        for (unsigned i=t1.to_unsigned(); i<=t2.to_unsigned() && i<len; i++)
            if (fixdates.count(i)) draw_phases[i] = 'B';

        r = phases[i].intervals[1];
        t1 = r[0]; t2 = r[1];
        if (t1.valid())
        for (unsigned i=t1.to_unsigned(); i<=t2.to_unsigned() && i<len; i++)
            if (fixdates.count(i)) draw_phases[i] = 'F';

        r = phases[i].intervals[2];
        t1 = r[0]; t2 = r[1];
        if (t1.valid())
        for (unsigned i=t1.to_unsigned(); i<=t2.to_unsigned() && i<len; i++)
            if (fixdates.count(i)) draw_phases[i] = 'O';

        r = phases[i].intervals[3];
        t1 = r[0]; t2 = r[1];
        if (t1.valid())
        for (unsigned i=t1.to_unsigned(); i<=t2.to_unsigned() && i<len; i++)
            if (fixdates.count(i)) draw_phases[i] = 'S';

        r = phases[i].intervals[4];
        t1 = r[0]; t2 = r[1];
        if (t1.valid())
        for (unsigned i=t1.to_unsigned(); i<=t2.to_unsigned() && i<len; i++)
            if (fixdates.count(i)) draw_phases[i] = 'D';


        t1 = restored_values[3*i    ];
        t2 = restored_values[3*i + 1];
        t3 = restored_values[3*i + 2];

        std::string draw_bounds(len, ' ');
        for (auto e : fixdates) draw_bounds[e] = '.';
        for (unsigned i=t1.to_unsigned(); i<t2.to_unsigned() && i<len; i++) draw_bounds[i] = 'b';
        for (unsigned i=t2.to_unsigned(); i<t3.to_unsigned() && i<len; i++) draw_bounds[i] = 'f';
        //for (unsigned i=t3.to_unsigned(); i<len; i++) draw_bounds[i] = 'o';

        if (print_flag)
        {
            printf("%8s : \n", flower.to_string().c_str());
            printf("%s \n", draw_phases.c_str());
            printf("%s \n", draw_bounds.c_str());
            FLUSH;
        }
    }
    return bounds;
}

std::vector<flowerline2_t> Flowerdome::interval_clipping(unsigned,
    interval<Flowerstamp> fs, phase_enum_t phase, bool print_flag) const
{
    std::vector<flowerline2_t> bounds = extract_phase_bounds(fs, phase, false);
    if (bounds.empty()) return bounds;

    //---------------------------------------------------------------------------------------------
    // конвертируем данные к формату функции intrainterval_localization_by_MLM()
    //---------------------------------------------------------------------------------------------
    std::vector<interval<double>> restored_values;
    std::vector<interval<double>> left_bounds, right_bounds;

    for (auto e : bounds)
    {
        interval<Timestamp> left = e[0];
        left_bounds.push_back(interval<double>(left[0].to_double(), left[1].to_double()));

        interval<Timestamp> right = e[1];
        right_bounds.push_back(interval<double>(right[0].to_double(), right[1].to_double()));
    }

    //---------------------------------------------------------------------------------------------
    //                          рассчитаем суженные интервалы
    //---------------------------------------------------------------------------------------------
    std::pair<double,double> moments = intrainterval_localization_by_directMLM(
        restored_values, left_bounds, right_bounds);

    //---------------------------------------------------------------------------------------------
    //                  вернем данные к формату возвращаемого значения
    //---------------------------------------------------------------------------------------------
    size_t i = 0;
    for (auto &e : bounds)
    {
        Flowerstamp flower = e.fsp;

        interval<Timestamp> left = e[0];
        Timestamp ta1 = left[0];
        Timestamp ta2 = left[1];

        interval<Timestamp> right = e[1];
        Timestamp tb1 = right[0];
        Timestamp tb2 = right[1];

        interval<double> &r = restored_values[i];
        Timestamp ta(r[0]);
        Timestamp tb(r[1]);

        Timestamp tta1(ta);
        Timestamp ttb1(tb - 1);

        int ma = std::min(tta1 - ta1, ttb1 - tb1);
        tta1 -= ma;
        ttb1 -= ma;

        Timestamp tta2(tta1 + 1);
        Timestamp ttb2(ttb1 + 1);
        int mb = std::min(ta2 - tta2, tb2 - ttb2);
        tta2 += mb;
        ttb2 += mb;

        assert(tta1 < tta2);
        assert(ttb1 < ttb2);

        left = interval<Timestamp>(tta1, tta2);
        right = interval<Timestamp>(ttb1, ttb2);

        if (print_flag)
        {
            printf("%8s : [%6s, %6s] - [%6s %6s] --> ", flower.to_string().c_str(),
                ta1.to_string().c_str(), (ta2-1).to_string().c_str(),
                tb1.to_string().c_str(), (tb2-1).to_string().c_str()
            );

            if (tta1.to_int() == tta2.to_int()-1)
            {
                if (tta1.to_int() == ta1.to_int()) printf("[%6s        ] - ", tta1.to_string().c_str());
                else printf("[        %6s] - ", tta1.to_string().c_str());
            }
            else printf("[%6s, %6s] - ", tta1.to_string().c_str(), (tta2-1).to_string().c_str());

            if (ttb1.to_int() == ttb2.to_int()-1)
            {
                if (ttb1.to_int() == tb1.to_int()) printf("[%6s        ]", ttb1.to_string().c_str());
                else printf("[        %6s]", ttb1.to_string().c_str());
            }
            else printf("[%6s, %6s]", ttb1.to_string().c_str(), (ttb2-1).to_string().c_str());

            int case_number_before = (ta2 - ta1) * (tb2 - tb1);
            int case_number_after = tta2 - tta1;
            printf(" : %3d -> %3d ", case_number_before, case_number_after);

            int days = ttb1 - tta1 + 1;
            if (case_number_after == 1) printf("[%2d]\n", days);
            else printf("     (%2d)\n", days);
            FLUSH
        }
        i++;
    }
    if (print_flag)
    {
        printf("\t average=%lf sigma=%lf\n", moments.first, moments.second);
        FLUSH
    }

    return bounds;
}

std::vector<flowerline2_t> Flowerdome::maximum_likelihood_flowering_ensemble(
    Unsigned num_iterations, interval<Flowerstamp> fs, bool skip_death, bool print_flag) const
{
    std::vector<flowerline3_t> bounds = extract_phase_bounds(fs, skip_death, print_flag);
    if (bounds.empty()) return std::vector<flowerline2_t>();

    //---------------------------------------------------------------------------------------------
    // конвертируем данные к формату функции intrainterval_localization_by_MLM()
    //---------------------------------------------------------------------------------------------
    Unsigned ne = bounds.size();

    std::vector<interval<double>> intervals;
    for (Unsigned i=0; i<ne; i++)
    {
        for (Unsigned j=0; j<3; j++)
        {
            interval<Timestamp> range = bounds[i].intervals[j];
            intervals.push_back(interval<double>{range[0].to_double(), range[1].to_double()});
        }
    }

    //---------------------------------------------------------------------------------------------
    //                          рассчитаем суженные интервалы
    //---------------------------------------------------------------------------------------------

    auto restored_values = intrainterval_localization_by_MLM(ne, 3,
        &intervals[0], num_iterations, print_flag);

    //---------------------------------------------------------------------------------------------
    //                  вернем данные к формату возвращаемого значения
    //---------------------------------------------------------------------------------------------
    std::vector<flowerline_t> return_values;

    for (Unsigned i=0; i<ne; i++)
    {
        Flowerstamp flower = bounds[i].fsp;

        if (print_flag) printf("%8s : ", flower.to_string().c_str());

        auto r0 = bounds[i].intervals[0];
        Timestamp ta1 = r0[0];
        Timestamp ta2 = r0[1];

        auto r1 = bounds[i].intervals[1];
        Timestamp tb1 = r1[0];
        Timestamp tb2 = r1[1];

        auto r2 = bounds[i].intervals[2];
        Timestamp tc1 = r2[0];
        Timestamp tc2 = r2[1];

        Timestamp t1 = restored_values[3*i    ];
        Timestamp t2 = restored_values[3*i + 1];
        Timestamp t3 = restored_values[3*i + 2];

        int days[2] = { t2.to_int() - t1.to_int(), t3.to_int() - t2.to_int() };

        if (print_flag)
        {
            printf("%8s : [%6s %6s %6s] [%6s %6s %6s] [%6s %6s %6s] BUD=%d FLOWER=%d \n",
                flower.to_string().c_str(),
                ta1.to_string().c_str(), t1.to_string().c_str(), ta2.to_string().c_str(),
                tb1.to_string().c_str(), t2.to_string().c_str(), tb2.to_string().c_str(),
                tc1.to_string().c_str(), t3.to_string().c_str(), tc2.to_string().c_str(),
                days[0], days[1]
            );
            FLUSH;
        }
        assert(t1.to_int() <= t2.to_int() && t2.to_int() <= t3.to_int());
        assert(ta1.to_int() <= t1.to_int() && t1.to_int() <= ta2.to_int());
        assert(tb1.to_int() <= t2.to_int() && t2.to_int() <= tb2.to_int());
        assert(tc1.to_int() <= t3.to_int() && t3.to_int() <= tc2.to_int());

        return_values.push_back(flowerline_t{flower, {{t1, t2-1}, {t2, t3-1}}});
                // начало фазы цветения t2+1 на день позже конца предыдущей фазы
    }

    return return_values;
}

moment_t Flowerdome::number_of_flowers(const std::vector<unsigned> &sprouts, unsigned axis) const
{
    std::map<unsigned, unsigned> parents; // побег и максимальное число цветков на нем
    for (auto e : flowers)
    {
        unsigned sprout = e.sprout();
        if (parents.count(sprout) == 0) parents[sprout] = 0;
        if (e.axis_order() == unsigned(axis)) parents[sprout]++;
    }

    std::vector<double> tmp;
    for (auto e : parents)
    {
        for (unsigned i=0; i<sprouts.size(); i++)
            if (e.first == sprouts[i]) tmp.push_back(e.second);
    }
    return statistics(tmp);
}




