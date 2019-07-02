#ifndef FLOWER_STRUCT_H
#define FLOWER_STRUCT_H

#include <string>
#include <vector>
#include <string>
#include <set>

#include "prg.os.h"
#include "defs.h"
#include "math.defs.h"

struct Sproutstamp
{
    unsigned year; // год
    unsigned sprout; // побег
    unsigned n1; // полное число цветков 1-й оси
    unsigned n2; // полное число цветков 2-й оси
};

/// база данных по структуре растений
class Flowerstruct : public std::vector<Sproutstamp>
{
    std::string species_name; // имя вида

public:

    /// загрузить базу данных из файла (выйти по исключению, если были ошибки)
    Flowerstruct(const std::string &filename)
    {
        bool flag = load_file(filename);
        if (!flag) throw std::exception();
    }

    /*!
     * \brief загрузка файла с контролем корректности его данных
     * \param filename полное имя файла с путем
     */
    bool load_file(const std::string &filename, bool print_flag=false);

    /*!
     * \brief извлечь данные по заданным годам для всех побегов
     * \return интервалы границ начала фазы и ее конца для каждого цветка fs, [a1, a2[, [b1, b2[
     */
    curve_t count_select(const std::vector<unsigned> &years, unsigned axe) const;
    curve_t count_select(const std::vector<unsigned> &years) const;

    curve_t rate_select(const std::vector<unsigned> &years) const;
};

#endif // FLOWER_STRUCT_H
