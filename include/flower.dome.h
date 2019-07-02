#ifndef FLOWER_DOME_H
#define FLOWER_DOME_H

#include <string>
#include <vector>
#include <string>
#include <set>
#include <tuple>
#include <memory>

#include "prg.os.h"

#include "defs.h"
#include "math.defs.h"
#include "math.statistics.h"
#include "record.time.h"
#include "record.phase.h"
#include "record.flower.h"

/// основной тип данных хранения начала и конца фазы для отдельного цветка
template <unsigned n>
struct flowerline
{
    Flowerstamp fsp;
    interval<Timestamp> intervals[n];
    interval<Timestamp> operator[](unsigned i) const { return intervals[i]; }
    void print() const;
};

template <unsigned n>
inline bool operator< (const flowerline<n> &a, const flowerline<n> &b) { return a.fsp < b.fsp; }

template <unsigned n>
inline void flowerline<n>::print() const
{
    printf("%8s : ", fsp.to_string().c_str());
    for (unsigned i=0; i<n; i++)
    {
        interval<Timestamp> t = intervals[i];
        std::printf("[%6s, %6s] ", t[0].to_string().c_str(), t[1].to_string().c_str());
    }
    std::printf("\n"); FLUSH;
}

typedef flowerline<1> flowerline1_t;
typedef flowerline<2> flowerline2_t;
typedef flowerline<3> flowerline3_t;
typedef flowerline<4> flowerline4_t;
typedef flowerline<5> flowerline5_t;
typedef flowerline2_t flowerline_t;

template <unsigned n> inline void print(const flowerline<n> &f) { f.print(); }

/// тестовая печать массива записей
inline void print(const std::vector<flowerline_t> &f) { for (auto &e : f) print(e); }

/// цветок с возможными датами его цветения
typedef std::pair<Flowerstamp, std::vector<interval<Timestamp>>> flower_ensemble;

/// Множество возможных дат цветения для заданного вида. Объекты множества упорядочены по номерам
/// цветков (требуется дополнительный вектор с этими номерами)
typedef std::vector<point<3>> species_ensemble;

/// получить число колонок имени и число дат в файле
std::pair<unsigned, unsigned> file_columns(const std::string &filename);

/// база данных по всем экземплярам растений
class Flowerdome : public std::vector<flower_t>
{
    // Реализация гарантирует сохранение следующего порядка записей: номера цветков по возрастанию,
    // (для заданного цветка) даты по возрастанию. Это позволяет делать быстрый поиск.

    typedef std::vector<flower_t> base_t;
    std::vector<Timestamp> dates; // упорядоченные даты в файле
    std::set<Flowerstamp> flowers; // все цветки файла
    std::string species_name; // имя вида

    std::map<unsigned, std::array<unsigned, 3>> flowers_count;
        // суммарное число цветов на осях (n1, n2) для заданного побега

public:

    /// создать пустую базу данных
    Flowerdome(){}

    /// загрузить базу данных из файла (выйти по исключению, если были ошибки)
    Flowerdome(const std::string &filename, const std::vector<unsigned> &sprouts={})
    {
        bool flag = load_file(filename, sprouts);
        if (!flag) throw std::exception();
    }

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
    bool load_file(const std::string &filename, const std::vector<unsigned> &sprouts={},
        bool print_flag=false);

    /// немного о загруженном файле
    void print_information() const;

    /*!
     * \brief извлечь интервалы границ одной(!) заданной фазы для всех цветков
     * \note используется функцией load_file для контроля файла
     * \note игнорирование DEATH концов, если фаза проявляется до DEATH
     * \param fs интервал цветков [первый включаемый, первый невключаемый[
     * \param phase фаза развития цветков
     * \param print_flag отладочная печать данных
     * \return интервалы границ начала фазы и ее конца для каждого цветка fs, [a1, a2[, [b1, b2[
     */
    std::vector<flowerline2_t> extract_phase_bounds(interval<Flowerstamp> fs, phase_enum_t phase,
        bool print_flag) const;

    /*!
     * \brief рассчитать статистику нормального распределения для заданной фазы
     * \note  используется функцией load_file для контроля файла
     * \param fs интервал цветков
     * \param phase контролируемая фаза
     */
    moment_t control_statistics(interval<Flowerstamp> fs, phase_enum_t phase,
        bool print_flag=false) const;


    /*!
     * \brief извлечь интервалы продолжительности всех фаз (BUD, FLOWER, OVERY, SEED, DEATH)
     * \note функция используется далее, чтобы выбрать междуфазные интервалы
     * \param fs интервал цветков
     * \param print_flag тестовая печать
     * \return интервалы фаз [фиксированное начало, фиксированный конец]
     */
    std::vector<flowerline5_t> extract_phases_intervals(interval<Flowerstamp> fs,
        bool print_flag) const;

    /*!
     * \brief извлечь интервалы границ фаз (BUD, FLOWER) для всех цветков
     * \param fs интервал цветков [первый включаемый, первый невключаемый[
     * \param skip_death удаление цветков, если DEATH произошло ранее FLOWER
     * \param phase фаза развития цветков
     * \param print_flag отладочная печать данных
     * \return интервалы границ начала фазы и ее конца для каждого цветка fs, [a1, a2[, [b1, b2[
     */
    std::vector<flowerline3_t> extract_phase_bounds(interval<Flowerstamp> fs, bool skip_death,
        bool print_flag) const;

    /*!
     * \brief уменьшение интервалов неопределенностей начала и конца фазы методом
     *  максимального правдоподобия
     * \param [empty] число итераций оптимизации интервалов (при прямом поиске игнорируется)
     * \param phase фаза
     * \return <номер цветка, интервалы начала и конца фазы>
     *
     * \note Общее число возможностей выбора начала и конца фазы в интервалах сужается.
     *  До вызова функции начало и конец в интервалах выбираются произвольно, что дает na * nb
     *  вариантов, где na и nb - число дней в начальном и завершающем интервалах, соответственно.
     *  После вызова функции известно число дней в фазе (это разность между стартовыми датами
     *  завершающего и начального интервалов), потому число вариантов равно na ([na1, nb1],
     *  [na1+1, nb1+1], ..
     */
    std::vector<flowerline2_t> interval_clipping(unsigned,
        interval<Flowerstamp> fs, phase_enum_t phase, bool print_flag) const;

    /*!
     * \brief уменьшение интервалов неопределенностей начала и конца фаз (BUD, FLOWER) методом
     *  максимального правдоподобия
     * \param num_iterations число итераций оптимизации интервалов
     * \param skip_death игнрировать записи, где DEATH произошло ранее FLOWER
     * \return <номер цветка, интервалы начала и конца фаз BUD + FLOWER>
     */
    std::vector<flowerline3_t> interval_clipping(unsigned num_iterations,
        interval<Flowerstamp> fs, bool skip_death, bool print_flag) const;

    /*!
     * \brief получить максимально вероятные границы FLOWER фазы
     * \param num_iterations число итераций подгона
     * \param fs интервал цветков
     * \return vector<fsp, BUD[начало фазы,конец фазы], FLOWER[начало фазы,конец фазы]>
     */
    std::vector<flowerline2_t> maximum_likelihood_flowering_ensemble(Unsigned num_iterations,
        interval<Flowerstamp> fs, bool skip_death, bool print_flag) const;

    /// возврат предыдущей даты наблюдения
    std::pair<Timestamp, bool> prev_date(Timestamp t) const
    {
        auto it = std::upper_bound(dates.rbegin(), dates.rend(), t,
            [](Timestamp a, Timestamp b){ return b < a; });
        return (it != dates.rend()) ? std::make_pair(*it, true) : std::make_pair(t, false);
    }

    /// возврат следущей даты наблюдения
    std::pair<Timestamp, bool> next_date(Timestamp t) const
    {
        auto it = std::upper_bound(dates.begin(), dates.end(), t,
            [](Timestamp a, Timestamp b){ return a < b; });
        return (it != dates.end()) ? std::make_pair(*it, true) : std::make_pair(t, false);
    }

    /// оператор доступа к записям базы
    const flower_t &operator[](size_t i) const { return base_t::at(i); }

    /// оператор доступа к записям базы
    const flower_t &operator()(const Flowerstamp &fs, const Phasestamp &ps) const
    {
        auto first = std::lower_bound(begin(), end(), flower_t(fs, Timestamp(), ps));
        assert(first != end());
        return *first;
    }
    flower_t &operator()(const Flowerstamp &fs, const Phasestamp &ps)
    {
        auto first = std::lower_bound(begin(), end(), flower_t(fs, Timestamp(), ps));
        assert(first != end());
        return *first;
    }

    /// число среднее цветков на осях заданного порядка для заданных побегов
    moment_t number_of_flowers(const std::vector<unsigned> &sprouts, unsigned axis) const;

    /// возврат массива номеров побегов
    std::vector<unsigned> get_sprouts() const
    {
        std::vector<unsigned> tmp;
        for (auto e : flowers_count) tmp.push_back(e.first);
        return tmp;
    }

    /// число цветков на оси заданного порядка для заданного побега
    unsigned number_of_flowers(unsigned sprout, unsigned axis) const
    { return flowers_count.at(sprout)[axis]; }

};

void print_flowering_bounds(const std::vector<flowerline3_t> &bounds, Unsigned i);
void print_flowering_phase(const std::pair<double, std::shared_ptr<species_ensemble>> &ens,
    Unsigned i);

#endif // FLOWER_DOME_H
