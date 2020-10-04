#ifndef SPROUT_H
#define SPROUT_H

#include <string>
#include <vector>
#include <string>
#include <set>
#include <tuple>
#include <memory>
#include <chrono>
#include <thread>

#include "container/vector.h"
#include "math/math.h"
#include "math/statistics.h"

#include "os/os.h"
#include "svg/primitive.h"

#include "orn/orn.defs.h"
#include "orn/time.h"
#include "orn/phase.h"
#include "orn/flower.h"
#include "orn/flowerline.h"
#include "orn/svg.h"
#include "orn/scheme.h"

using namespace std::chrono_literals;
namespace prg
{
    /// контроль на пустую запись в файле csv (проверка строки типа ,,,,,,,)
    inline bool empty_line(const char *buffer, Unsigned num_names)
    {
        std::stringstream ss(buffer);
        char buf[128];
        for (Unsigned i=0; i<num_names; i++)
        {
            ss.getline(buf, BUF_LEN, ',');
            char *p = prg::trim(buf, remove_symbols);
            if (std::strlen(p) != 0) return false;
        }
        return true;
    }

    /*!
     * \brief база данных по динамике развития экземпляров растений
     * \note хранение цветков по возрастанию номеров для быстрого поиска
     */
    class Sprout : public std::vector<flowerline>
    {
        typedef std::vector<flowerline> Base;

    public:

        /// создать пустой объект
        Sprout(){}

        /*!
         * \brief загрузка побег с контролем корректности его данных
         * \param file файловый поток загрузки
         * \param dim {1+число осей, число дат} размеры таблицы данных
         * \param print_flag печатать/не печатать диагностику потенциальных ошибок в файле
         * \note Функция не генерирует исключение, чтобы можно было использовать ее при множественном
         *  тестировании файлов. Реальные ошибки, обнаруженные в файле, печатаются независимо от
         *  параметра print_flag, и при этом происходит прерывание загрузки файла. Потенциальные
         *  ошибки (например, выход значений за достоверные интервалы) печатаются только если
         *  установлен флаг print_flag.
         */
        bool load_csv(std::ifstream &file, const paired<Unsigned> &dim, prg::vector<Time> *dates, bool print_flag=false);
        bool load_orn(std::ifstream &file, prg::vector<Time> *dates, bool print_flag=false);

        /*!
         * \brief сохранение побега в файле
         * \param file файловый поток
         * \note формат записи отличается от csv
         *      START SPROUT id n1,n2,n3,..
         *      1-1 30-Jun bud 02-Jul flower ..
         *      END SPROUT
         * \param sproutId идентификатор, под которым записывается побег
         * \param start дата старта, с которой идет запись в общей базе
         */
        void save(std::ofstream &file, Unsigned sproutId) const;

        /// идентификатор побега
        Unsigned sprout_index() const { return Base::at(0).fsp.sprout(); }

        /// оператор доступа к записям базы по номеру записи (или nullptr)
        const flowerline *operator[](Unsigned i) const { return i<size() ? &Base::at(i) : nullptr; }
        flowerline *operator[](Unsigned i) { return i<size() ? &Base::at(i) : nullptr; }

        /// оператор доступа к записям базы по номеру цветка (или nullptr)
        const flowerline *operator[](const FlowerIndex &fs) const
        {
            Unsigned i = std::distance(begin(),
                std::lower_bound(begin(), end(), fs, [](auto a, auto b) { return a.fsp < b; })
            );
            if (i == size()) return nullptr;

            const flowerline *line = (*this)[i];
            return line->fsp == fs ? line : nullptr;
        }
        flowerline *operator[](const FlowerIndex &fs)
        {
            Unsigned i = std::distance(begin(),
                std::lower_bound(begin(), end(), fs, [](auto a, auto b) { return a.fsp < b; })
            );
            if (i == size()) return nullptr;

            flowerline *line = (*this)[i];
            return line->fsp == fs ? line : nullptr;
        }

        interval<Unsigned> range(const FlowerIndex &fs1, const FlowerIndex &fs2) const
        {
            Unsigned i1 = std::distance(begin(),
                std::lower_bound(begin(), end(), fs1, [](auto a, auto b) { return a.fsp < b; })
            );
            Unsigned i2 = std::distance(begin(),
                std::lower_bound(begin(), end(), fs2, [](auto a, auto b) { return a.fsp < b; })
            );
            return {i1, i2};
        }

        //----------------------------------------------------------------------------------------
        //                                  STATISTICS
        //----------------------------------------------------------------------------------------

        /// число цветков на главной оси (для получения относительной высоты)
        Unsigned number_of_axis1_flowers() const { return (*this)[0]->fsp[1]; }
        Unsigned number_of_axis2_flowers() const { return (*this)[0]->fsp[2]; }

        /// день начала фазы цветения побега
        Time day_of_first_start(Phase phase) const;
        Time day_of_last_start(Phase phase) const;

        /// печать результата загрузки
        void print_information() const;

        /*!
         * \brief набирает информацию о любых событиях для построения функции распределения
         * \note событиями могут быть: продолжительность фазы цветка, число цветков на побеге
         * \note каждому записи соответствует одно событие в данных (не по порядку номеров записей)
         *
         * \param sprouts номера побегов, с которых собирается информация
         * \param axes номера осей, с которых собирается информация
         * \param op операция, выполняемая над записью
         * \return одномерный НЕупорядоченный массив наличия событий, типа {2, 12, 32, .., 45},
         */
        template <typename Oper>
        shared<curve<Oper::n>> foreach(Oper op) const
        {
            shared<curve<Oper::n>> p(new curve<Oper::n>);
            for (auto e : *this)
            {
                auto r = op(e);
                auto flag = r.flag();
                auto v = r.value();

                if (flag) p->push_back(v);
            }
            std::sort(p->begin(), p->end(), [](auto a, auto b){ return a[0] < b[0]; });
            return p;
        }

        /// генерация схемы побега для отрисовки
        SproutScheme scheme_generate(Time t=0.) const;

        /// вывод схемы побега
        void scheme_print(const std::string &filename, const SproutScheme &scheme,
            const std::string &caption="") const;

        void scheme_video(const std::string &filename, Time start, Time finish, std::chrono::seconds delay=2s) const;
    };

    inline void Sprout::scheme_video(const std::string &filename, Time start, Time finish,
        std::chrono::seconds delay) const
    {
        for (Time t=start; t<finish; ++t)
        {
            SproutScheme scheme = scheme_generate(t);
            std::string caption = t.to_string();
            scheme_print(filename, scheme, caption);
            std::this_thread::sleep_for(delay);
        }
    }

    inline Time Sprout::day_of_first_start(Phase phase) const
    {
        Time ts = 365.;
        for (Unsigned i=1; i<size(); i++)
        {
            const auto &e = (*this)[i];
            const FlowerIndex &f = e->fsp;
            if (f.axis() != 1) continue;
            Time t = e->value(phase);
            if (ts > t) ts = t;
        }
        return ts;
    }

    inline Time Sprout::day_of_last_start(Phase phase) const
    {
        Time ts = 0.;
        for (Unsigned i=1; i<size(); i++)
        {
            const auto &e = (*this)[i];
            Time t = e->value(phase);
            if (ts < t) ts = t;
        }
        return ts;
    }

    inline bool Sprout::load_csv(std::ifstream &file, const paired<Unsigned> &n, prg::vector<Time> *dates, bool print_flag)
    {
        bool valid_flag = true; // сигнал отсутствия ошибок в файле
        char buffer[BUF_LEN], buf[BUF_LEN];

        Unsigned num_names = n.first, num_dates = n.second;

        flowerline record; // текущая загружаемая запись, по умалчиванию динамическая
        FlowerIndex fcur = 0;

        file.getline(buffer, BUF_LEN);
        Unsigned sprout = 0;
        while (!file.eof())
        {
            // контроль на пустую запись
            if (empty_line(buffer, num_names)) break;

            std::stringstream ss(buffer);
            // не обнуляем предыдущее имя fcur, чтобы заполнить пустые поля в начале

            ss.getline(buf, BUF_LEN, ',');
            char *p = prg::trim(buf, remove_symbols);
            if (std::strlen(p) != 0)
            {
                sprout = atoi(p);
                fcur.sprout() = sprout;
                fcur.flag() = false;
            }
            else fcur.flag() = true;

            // считывание/установка отдельных индексов
            for (Unsigned i=1; i<num_names; i++)
            {
                ss.getline(buf, BUF_LEN, ',');
                char *p = prg::trim(buf, remove_symbols);
                if (std::strlen(p) != 0)
                {
                    Unsigned val = atoi(p);
                    fcur[i] = val; // индексы цветков | число цветков

                    for (Unsigned j=i+1; j<num_names; j++) fcur[j] = 0;
                        // обнуление вышележащих порядков
                }
            }
            record.fsp = fcur;

            if (fcur.flag())
            {
                Phase phprev(UNDEF), phcur(UNDEF);
                Time tprev = (*dates)[-1], tcur = (*dates)[-1];
                for (Unsigned i=0; i<num_dates; i++)
                {
                    tcur = (*dates)[i];

                    ss.getline(buf, BUF_LEN, ',');
                    char *p = prg::trim(buf, remove_symbols);
                    phcur = make_phase(phprev, std::string(p));

                    if (phcur == DEATH) break;
                        // следующая граница фазы отсутствует

                    if (phcur == INVALID)
                    {
                        std::cout << std::printf("error: invalid phase \'%s\' at (%s, %s)\n",
                            p, fcur.to_string().c_str(), tcur.to_string().c_str());
                        std::cout.flush();
                        valid_flag = false;
                        continue;
                    }

                    /// вставка всех пропущенных фаз при phprev != phcur
                    for (Phase p(phprev+1); p<=phcur; p=next_enum_value(p))
                        *record.start_interval(p) = {{tprev, tcur}, tcur};

                    phprev = phcur;
                    tprev = tcur;
                }
                if (phcur != INVALID && phcur != DEATH)
                {
                    phcur = next_enum_value(phprev);
                    *record.start_interval(phcur) = {{tprev, Time("1-Oct")}, Time("1-Oct")};
                        // изменение фазы за пределами наблюдений
                }

                if (print_flag && (phcur!= DEATH) && (phcur < Phase(OVARY)))
                {
                    std::cout << std::printf("[warning]: unclosed record \'%s\'\n",
                        fcur.to_string().c_str());
                    std::cout.flush();
                }
            }
            push_back(record);
            file.getline(buffer, BUF_LEN);
        }

        // обеспечение гарантии упорядочения записей если нет ошибок в файле
        if (valid_flag) std::sort(begin(), end(), [](auto a, auto b) { return a.fsp < b.fsp; });

        //---------------------------------------------------------------------------------------------
        //         сосчитаем число цветков на осях 1-го и 2-го порядков, если их нет
        //---------------------------------------------------------------------------------------------
        flowerline *p = (*this)[0]; // структурная запись является начальной
        if (p->fsp.axis() == 0)  // данные для структурной записи не установлены
        {
            std::array<unsigned, 4> sum = {0, 0, 0, 0};
            for (Unsigned i=1; i<size(); i++)
            {
                FlowerIndex f = (*this)[i]->fsp;
                Unsigned axis = f.axis();
                sum[axis-1]++;
            }
            p->fsp = sum;
        }

        if (print_flag) print_information();
        return valid_flag;
    }

    inline bool Sprout::load_orn(std::ifstream &file, prg::vector<Time> *dates, bool print_flag)
    {
        bool valid_flag = true; // сигнал отсутствия ошибок в файле
        char buffer[BUF_LEN];
        flowerline record; // текущая загружаемая запись, по умалчиванию динамическая

        // load structural record
        std::string buf; Unsigned sprout; std::array<unsigned, 4> sum = {0, 0, 0, 0};
        file.getline(buffer, BUF_LEN);
        std::stringstream ss(buffer);
        ss >> buf >> buf; // START SPROUT 1 11 3
        ss >> sprout >> sum[0] >> sum[1] >> sum[2] >> sum[3];
        record.fsp = FlowerIndex(sprout, false, sum);
        push_back(record);

        file.getline(buffer, BUF_LEN);
        for(; !file.eof(); file.getline(buffer, BUF_LEN))
        {
            if (std::strncmp("END", buffer, 3) == 0) break;

            auto tokens = split(buffer);
            record.fsp = FlowerIndex(sprout, true, tokens[0]);

            Phase phase = BUD;
            for (Unsigned i=1; i<tokens.size(); i++)
            {
                Time t(tokens[i].c_str());
                auto it = std::lower_bound(dates->begin(), dates->end(), t);
                Signed ndx = dates->minimal_index() + std::distance(dates->begin(), it);
                Time tprev = (*dates)[ndx - 1];
                 *record.start_interval(phase) = {{tprev, t}, t};
                phase = phase + 1;
            }
            push_back(record);
        }

        // обеспечение гарантии упорядочения записей если нет ошибок в файле
        if (valid_flag) std::sort(begin(), end(), [](auto a, auto b) { return a.fsp < b.fsp; });

        //---------------------------------------------------------------------------------------------
        //         сосчитаем число цветков на осях 1-го и 2-го порядков, если их нет
        //---------------------------------------------------------------------------------------------
        flowerline *p = (*this)[0]; // структурная запись является начальной
        if (p->fsp.axis() == 0)  // данные для структурной записи не установлены
        {
            std::array<unsigned, 4> sum = {0, 0, 0, 0};
            for (Unsigned i=1; i<size(); i++)
            {
                FlowerIndex f = (*this)[i]->fsp;
                Unsigned axis = f.axis();
                sum[axis-1]++;
            }
            p->fsp = sum;
        }

        if (print_flag) print_information();
        return valid_flag;
    }

    inline void Sprout::save(std::ofstream &file, Unsigned sproutId) const
    {
        // START SPROUT id n1,n2,n3,..
        //     1-1 30-Jun[bud] 02-Jul[flower] ..
        // END SPROUT

        const flowerline &line = Base::at(0);
        FlowerIndex fsp = line.fsp;
        const std::array<unsigned,4> &f = fsp.axes();
        file << "START SPROUT " << std::to_string(sproutId);
        if (f[0] != 0) file << " " << std::to_string(f[0]);
        if (f[1] != 0) file << " " << std::to_string(f[1]);
        if (f[2] != 0) file << " " << std::to_string(f[2]);
        if (f[3] != 0) file << " " << std::to_string(f[3]);
        file << "\n";

        for (Unsigned i=1; i<size(); i++)
        {
            const flowerline &line = Base::at(i);
            FlowerIndex fsp = line.fsp;

            file << fsp.to_short_string() << "\t";
            Time t1 = line.value(BUD);
            if (t1.valid()) file << t1.to_string() << "\t";

            Time t2 = line.value(FLOWER);
            if (t2.valid()) file << t2.to_string() << "\t";

            Time t3 = line.value(OVARY);
            if (t3.valid()) file << t3.to_string() << "\t";

            Time t4 = line.value(SEED);
            if (t4.valid()) file << t4.to_string() << "\t";

            Time t5 = line.value(DEATH);
            if (t5.valid()) file << t5.to_string() << "\t";

            file << "\n";
        }
        file << "END SPROUT \n";
    }

    inline void Sprout::print_information() const
    {
        const flowerline *line = (*this)[0];

        FlowerIndex f = line->fsp;
        Unsigned sprout = f.sprout();
        std::printf("sprout %ld : #flowers on axes {%d %d %d %d}\n",
            sprout, f[1], f[2], f[3], f[4]);
    }

    /*!
     * \brief генерация схемы соцветия побега на заданную дату
     * \note при t=0 генерация общей схемы соцветия без указания фаз развития
     */
    inline SproutScheme Sprout::scheme_generate(Time t) const
    {
        struct Node
        {
            Phase phase = UNDEF;
            std::vector<Node> children;

            void insert(const FlowerIndex &fs, Phase ph, Unsigned level=0)
            {
                Unsigned axis = fs.axis();
                if (axis == level) { phase = ph; return; }

                Unsigned i = fs[level+1]; // номер ребенка для текущего перехода
                if (i > children.size()) children.resize(i);
                children[i-1].insert(fs, ph, level+1);
            }

            SproutScheme generate(Signed level=0)
            {
                std::vector<SproutScheme> v;
                if (children.empty()) return SproutScheme(FLOWER_WIDTH, FLOWER_HEIGHT, phase);
                for (auto e : children) v.push_back(e.generate(level+1));
                return SproutScheme::construct(SproutScheme(FLOWER_WIDTH, FLOWER_HEIGHT, phase), v);
            }
        };

        Node tree;
        for (auto e : *this)
            if (e.dynamic_flag())
            {
                FlowerIndex f = e.fsp;
                if (f.axis() == 1 && f.position() == 1)
                {
                    // терминал главной оси вставляем отдельно из-за нарушения нумерации
                    tree.phase = t <= 0. ? UNDEF : e.phase(t);
                    continue;
                }
                f[1] = f[1] - 1; // корректируем нумерацию на 1 оси, чтобы не выступал терминал

                tree.insert(f, t <= 0. ? UNDEF : e.phase(t));
            }

        return tree.generate();
    }

    inline void Sprout::scheme_print(const std::string &filename, const SproutScheme &scheme,
        const std::string &caption) const
    {
        std::vector<ShemeRecord> scheme_records = scheme.data();

        Dimensions dimensions(400, 600);
        Document doc(filename, Layout(dimensions, Layout::BottomLeft));

        svgScheme svg_scheme(scheme_records);
        svg_scheme.scale(doc.layout);
        svg_scheme.center(doc.layout);
        doc << svg_scheme;

        doc << Text(Point(20,20), caption, Fill(Color::Black),
             Font(12, "Verdana"), Stroke(0.1, Color::Black));

        doc.save();
    }

}
#endif // SPROUT_H
