#ifndef FLOWER_DOME_H
#define FLOWER_DOME_H

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
#include "orn/recovery.h"
#include "orn/sprout.h"

namespace prg
{
    auto is_word = [](std::string s) -> bool
    {
        if (!s.empty() && std::isalpha(s[0])) return true;
        return false;
    };

    /// получить число колонок имени и число дат в файле
    inline std::pair<unsigned, unsigned> csv_file_columns(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file)
        {
            std::cout << "File opening failed\n";
            return {0, 0};
        }

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

    inline prg::vector<Time> orn_file_dates(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file)
        {
            std::cout << "File opening failed\n";
            return {0, 0};
        }

        std::set<Time> times;
        char buffer[BUF_LEN];

        file.getline(buffer, BUF_LEN);
        for (; !file.eof(); file.getline(buffer, BUF_LEN))
        {
            std::vector<std::string> tokens = split(buffer);
            if (tokens.empty()) continue;
            if (is_word(tokens[0])) continue;

            for (Unsigned i=1; i<tokens.size(); i++)
            {
                assert(is_date(tokens[i]));
                std::string ss = tokens[i];
                times.insert(Time(tokens[i].c_str()));
            }
        }

        // найдем значение date[-1]
        std::vector<double> tmp;
        Time prev; // invalid by default
        for (auto e : times)
        {
            if (prev.valid())
            {
                double dt = e - prev;
                tmp.push_back(dt);
            }
            prev = e;
        }
        double aver = std::floor(median(tmp));

        prg::vector<Time> dates;
        dates.push_back(*times.begin() - aver);
        for (auto e : times) dates.push_back(e);
        dates.minimal_index() = -1;

        return dates;
    }

    /*!
     * \brief база данных по динамике развития экземпляров растений
     * \note хранение цветков по возрастанию номеров для быстрого поиска
     */
    class Flowerdome : public std::map<Unsigned, Sprout>
    {
        typedef std::map<Unsigned, Sprout> Base;

        prg::vector<Time> dates; // упорядоченные даты наблюдений в файле,
            // включая dates[-1] дату до начала записывания данных
            // prg::vector vs std::vector

    public:

        std::string species_name; // имя вида

        /// создать пустую базу данных
        Flowerdome(){}

        /*!
         * \brief загрузка файла с контролем корректности его данных
         * \param filename полное имя файла с путем
         * \param print_flag печатать/не печатать диагностику потенциальных ошибок в файле
         * \note Функция не генерирует исключение, чтобы можно было использовать ее при множественном
         *  тестировании файлов. Реальные ошибки, обнаруженные в файле, печатаются независимо от
         *  параметра print_flag, и при этом происходит прерывание загрузки файла. Потенциальные
         *  ошибки (например, выход значений за достоверные интервалы) печатаются только если
         *  установлен флаг print_flag.
         */
        std::string load(const std::string &filename, bool print_flag=false)
        {
            std::string ext = extension(filename);
            if (ext == ".csv") return load_csv(filename, print_flag);
            else if (ext == ".orn") return load_orn(filename, print_flag);
            return "";
        }

        std::string load_csv(const std::string &filename, bool print_flag=false);
        std::string load_orn(const std::string &filename, bool print_flag=false);

        void save(const std::string &filename);

        /// немного о загруженном файле
        std::string get_info() const
        {
            std::string info = "";

            std::set<Unsigned> ids;
            for (auto sprout : *this)
            {
                if (sprout.second.size() > 1)
                    ids.insert(sprout.first);
            }

            Unsigned i1 = *ids.begin(), i2 = i1;
            for (auto i : ids)
            {
                if (i == i2) continue;
                if (i == i2 + 1) i2 = i;
                else
                {
                    info += std::to_string(i1);
                    if (i1 != i2) info += "-" + std::to_string(i2);
                    info += ",";
                    i1 = i; i2 = i;
                }
            }
            info += std::to_string(i1);
            if (i1 != i2) info += "-" + std::to_string(i2);

            return info;
        }

        /// возврат номеров побегов
        std::set<Unsigned> get_sprouts() const
            { std::set<Unsigned> tmp; for (auto key : *this) tmp.insert(key.first); return tmp; }

        const Sprout &sprout(Unsigned i) const { return at(i); }
        Sprout &sprout(Unsigned i) { return at(i); }
        bool has_sprout(Unsigned i) const { return count(i) > 0; }

        Unsigned num_sprouts() const { return size(); }
        Unsigned num_dynamic_sprouts() const
        {
            Unsigned num = 0;
            for (auto sprout : *this) if (sprout.second.size() > 1) num++;
            return num;
        }

        /// возврат предыдущей даты наблюдения
        std::pair<Time, bool> prev_date(Time t) const
        {
            auto it = std::upper_bound(dates.rbegin(), dates.rend(), t,
                [](Time a, Time b){ return b < a; });
            return (it != dates.rend()) ? std::make_pair(*it, true) : std::make_pair(t, false);
        }

        /// возврат следущей даты наблюдения
        std::pair<Time, bool> next_date(Time t) const
        {
            auto it = std::upper_bound(dates.begin(), dates.end(), t,
                [](Time a, Time b){ return a < b; });
            return (it != dates.end()) ? std::make_pair(*it, true) : std::make_pair(t, false);
        }

        //----------------------------------------------------------------------------------------
        //                                  STATISTICS
        //----------------------------------------------------------------------------------------

        /// число цветков на главной оси (для получения относительной высоты)
        Unsigned number_of_axis1_flowers(Unsigned sprout) const { return at(sprout).number_of_axis1_flowers(); }

        /// день начала фазы цветения побега
        Time day_of_fist_start(Unsigned sprout, Phase phase) const { return at(sprout).day_of_first_start(phase); }

        /// статистика по длине фазы для цветков указанных побегах и осях заданного порядка
        moment_t duration_of_phase(Phase phase, const std::set<Unsigned> &sprouts) const;

        /// проверка наличия сильных отклонений от нормальной статистики и их печать
        bool control_statistics(double pvalue=0.05) const;

        /*!
         * \brief установить максимально вероятные границы фазы (в предположении нормальности длины фазы)
         * \note при расчете выбрасываются цветки, для которых границы не определены
         * \return средняя продолжительность фазы и ее дисперсия
         */
        std::pair<double,double> recover_missed_data_by_direct_MLM(Phase phase);
        std::pair<double,double> recover_missed_data_by_direct_MLM();
        std::pair<double,double> recover_missed_data_by_random();

        /*!
         * \brief набирает информацию о любых событиях для построения функции распределения
         * \note событиями могут быть: продолжительность фазы цветка, число цветков на побеге
         * \note каждому записи соответствует одно событие в данных (не по порядку номеров записей)
         *
         * \param sprouts номера побегов, с которых собирается информация
         * \param op операция, выполняемая над записью
         * \return одномерный НЕупорядоченный массив наличия событий, типа {2, 12, 32, .., 45},
         */
        template <typename Oper>
        shared<curve<Oper::n>> foreach(const std::set<Unsigned> &sprouts, Oper op) const
        {
            shared<curve<Oper::n>> p(new curve<Oper::n>);
            for (Unsigned isprout : sprouts)
            {
                for (auto e : at(isprout))
                {
                    auto r = op(e);
                    auto flag = r.flag();
                    auto v = r.value();
                    if (flag) p->push_back(v);
                }
            }
            std::sort(p->begin(), p->end(), [](auto a, auto b){ return a[0] < b[0]; });
            return p;
        }

    };

    inline moment_t Flowerdome::duration_of_phase(Phase phase,
        const std::set<Unsigned> &sprouts) const
    {
        auto f = foreach(sprouts, phase_duration_functor(phase));

        std::vector<double> tmp;
        for (auto e : *f) tmp.push_back(e[1]);
        return statistics(tmp);
    }

    inline bool Flowerdome::control_statistics(double pvalue) const
    {
        bool flag = true;
        for (auto phase : { BUD, FLOWER })
        {
            // получение нормальной статистики
            moment_t moments = duration_of_phase(phase, get_sprouts());
            double mean = moments[0];
            double sigma = moments[1];

            std::cout << to_string("\n[test] of normal statistics for phase \'%s\': mean=%4.2lf sigma=%4.2lf\n",
                to_string(phase).c_str(), mean, sigma);

            Unsigned errors = 0;
            // печать маловероятных значений с pvalue < заданного
            for (auto sprout : *this)
            {
                for (auto line : sprout.second)
                {
                    auto dt = line.duration_of_phase(phase);
                    if (!dt.flag()) continue;

                    double pv = prg::pvalue(dt, mean, sigma);
                    if (pv < pvalue)
                    {
                        std::cout << to_string("\t%s pvalue=%lf\n", line.fsp.to_string().c_str(), pv);
                        errors++;
                    }
                }
            }

            if (errors > Unsigned(std::ceil(size() * pvalue)))
            {
                std::cout << std::printf("[warning] the number of errors allowed [%ld] for phase %s "
                    "is higher than the threshold %lf%%\n",
                    errors, to_string(phase).c_str(), pvalue);
                std::cout << to_string("[test] something is bad\n");
                std::cout.flush();
                flag = false;
            }
            else if (errors > 0)
            {
                std::cout << to_string("[test] statistics is more or less acceptable\n");
                std::cout.flush();
            }
            else
            {
                std::cout << to_string("[test] things are good\n");
                std::cout.flush();
            }
        }
        return flag;
    }

    inline std::string Flowerdome::load_csv(const std::string &filename, bool print_flag)
    {
        if (filename.empty()) return "filename error";
        species_name = prg::basename(filename); // сохраним имя вида

        bool valid_flag = true; // сигнал отсутствия ошибок в файле
        char buffer[BUF_LEN], buf[BUF_LEN];

        std::pair<Unsigned, Unsigned> n = csv_file_columns(filename);
        if (n.first == 0 || n.second == 0) return "file error";

        Unsigned num_names = n.first, num_dates = n.second;
        dates.resize(-1, 1 + num_dates);
            // установим смещенную нумерацию элементов массива
            // добавив день до наблюдений

        std::ifstream file(filename);
        if (!file) return "file opening error";

        file.getline(buffer, BUF_LEN); // разбор имен колонок для дат
        {
            std::stringstream ss(buffer);

            for (Unsigned i=0; i<num_names; i++) ss.getline(buf, BUF_LEN, ',');
                // удалим из потока первые num_names значений

            Time prev;
            for (Unsigned i=0; i<num_dates; i++)
            {
                ss.getline(buf, BUF_LEN, ',');
                const char *p = prg::trim(buf, remove_symbols);
                bool flag = is_date(p);
                if (!flag)
                {
                    std::cout << std::printf("error: %s invalid date \'%s\'\n",
                        filename.c_str(), buf);
                    valid_flag = false;
                }

                Time cur(p);
                if (cur < prev)
                {
                    std::cout << std::printf("error: %s wrong date order \'%s\'\n",
                        filename.c_str(), buf);
                    valid_flag = false;
                }

                dates[i] = cur;
                prev = cur;
            }

            // оценка даты до начала наблюдений, как отстоящую от первой на средний
            // интервал между датами наблюдений
            std::vector<double> tmp;
            for (Unsigned i=1; i<num_dates; i++) tmp.push_back(dates[i] - dates[i-1]);
            double aver = std::floor(average(tmp));
            dates[-1] = dates[0] - aver;
        }

        while (!file.eof())
        {
            auto filepos = file.tellg();
            file.getline(buffer, BUF_LEN);

            if (file.eof()) break;
            if (empty_line(buffer, num_names)) continue;

            file.seekg(filepos);

            Sprout sprout;
            bool flag = sprout.load_csv(file, n, &dates, print_flag);
            insert({sprout.sprout_index(), std::move(sprout)});

            valid_flag = valid_flag && flag;
        }

        std::sprintf(buffer, "%ld sprouts are loaded from file \'%s\'", size(), filename.c_str());
        if (print_flag) control_statistics();

        return buffer;
    }

    inline std::string Flowerdome::load_orn(const std::string &filename, bool print_flag)
    {
        if (filename.empty()) return "filename error";
        dates = orn_file_dates(filename);

        bool valid_flag = true; // сигнал отсутствия ошибок в файле
        char buffer[BUF_LEN];

        std::ifstream file(filename);
        if (!file) return "file opening error";

        auto file_pos = file.tellg();
        file.getline(buffer, BUF_LEN);
        for ( ; !file.eof(); file_pos = file.tellg(), file.getline(buffer, BUF_LEN))
        {
            auto tokens = split(buffer);
            if (tokens.empty()) continue;

            if (is_word(tokens[0]) && tokens[0] == "SPECIES")
                species_name = tokens.size() > 1 ? tokens[1] : "unknown";

            if (is_word(tokens[0]) && tokens[0] == "START")
            {
                file.seekg(file_pos);

                Sprout sprout;
                bool flag = sprout.load_orn(file, &dates, print_flag);
                insert({sprout.sprout_index(), std::move(sprout)});

                valid_flag = valid_flag && flag;
            }
        }
        std::sprintf(buffer, "%ld sprouts are loaded from file \'%s\'", size(), filename.c_str());
        if (print_flag) control_statistics();

        return buffer;
    }

    inline void Flowerdome::save(const std::string &filename)
    {
        std::ofstream file(filename);
        if (!file) no_file_exit(filename);

        species_name = prg::basename(filename); // сохраним имя вида
        file << "SPECIES " << species_name << "\n";

        for (auto key_value : *this)
        {
            Unsigned id = key_value.first;
            const Sprout &sprout = key_value.second;
            sprout.save(file, id);
        }
        std::printf("%ld sprouts are saved to file \'%s\'\n", size(), filename.c_str());
    }

    inline std::pair<double,double> Flowerdome::recover_missed_data_by_direct_MLM(Phase phase)
    {
        std::vector<rough_value<Time>*> lb;
        std::vector<rough_value<Time>*> rb;

        for (auto &mapval : *this)
        {
            for (auto &e : mapval.second)
            {
                if (!e.dynamic_flag()) continue;

                auto l = e.start_interval(phase);
                auto r = e.start_interval(phase+1);
                if (l->range().empty() || r->range().empty()) continue;

                lb.push_back(l);
                rb.push_back(r);
            }
        }
        return intrainterval_localization_by_directMLM(lb, rb);
    }

    inline std::pair<double,double> Flowerdome::recover_missed_data_by_direct_MLM()
    {
        std::vector<rough_value<Time>*> b1;
        std::vector<rough_value<Time>*> b2;
        std::vector<rough_value<Time>*> b3;

        for (auto &mapval : *this)
        {
            for (auto &e : mapval.second)
            {
                if (!e.dynamic_flag()) continue;

                auto l1 = e.start_interval(BUD);
                if (l1->range().empty()) continue;

                auto l2 = e.start_interval(FLOWER);
                if (l2->range().empty()) continue;

                auto l3 = e.start_interval(OVARY);
                if (l3->range().empty()) continue;

                b1.push_back(l1);
                b2.push_back(l2);
                b3.push_back(l3);
            }
        }
        return intrainterval_localization_by_directMLM(b1, b2, b3);
    }

    inline std::pair<double,double> Flowerdome::recover_missed_data_by_random()
    {
        std::vector<rough_value<Time>*> b1;
        std::vector<rough_value<Time>*> b2;
        std::vector<rough_value<Time>*> b3;

        for (auto &mapval : *this)
        {
            for (auto &e : mapval.second)
            {
                if (!e.dynamic_flag()) continue;

                auto l1 = e.start_interval(BUD);
                if (l1->range().empty()) continue;

                auto l2 = e.start_interval(FLOWER);
                if (l2->range().empty()) continue;

                auto l3 = e.start_interval(OVARY);
                if (l3->range().empty()) continue;

                b1.push_back(l1);
                b2.push_back(l2);
                b3.push_back(l3);
            }
        }
        intrainterval_localization_by_random(b1, b2, b3);
        return std::pair<double,double>(0, 0);
    }
}
#endif // FLOWER_DOME_H
