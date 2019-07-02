#ifndef RECORD_FLOWER_H
#define RECORD_FLOWER_H

#include <cstdint>
#include <string>
#include <vector>

#include "defs.h"
#include "record.time.h"
#include "record.phase.h"

class Flowerstamp
{
    union
    {
        unsigned u;
        char b[4]; // цифровое представление имени в виде последовательности сегментов,
            // включающее (номер побега[3], номер цветка[2, 1, 0]), нумерация с 1.
            // Разделение полей делается умножением на 16, что позволяет просмотр как 0xFF.. в отладке.
            // Более удобное разделение степенью 10 не используется, номера цветков > 10.
            // Номер цветка представляет собой ряд сегментов, разделенных символом "-".
            // Символы разделения "-" не хранятся. Строковое представление "1-10-1".
            // Порядок хранения : начальный сегмент в самом верхнем байте, остальные ниже в памяти.
            // Это позволяет имеет стандартное упорядочивание имен.
    };

public:

    Flowerstamp(unsigned sprout=0, unsigned f1=0, unsigned f2=0, unsigned f3=0)
    : u(0) { b[3] = char(sprout); b[2] = char(f1); b[1] = char(f2); b[0] = char(f3); }

    Flowerstamp &operator=(unsigned u) { this->u = u; return *this; }

    /// установка значения заданного сегмента имени цветка
    void set(unsigned i, unsigned v) { assert(i<4); b[3-i] = char(v); }

    /// взятие значения заданного сегмента имени цветка
    unsigned get(unsigned i) const { assert(i<4); return unsigned(b[3-i]);}

    unsigned sprout() const { return unsigned(b[3]); }
    std::string flower() const
    {
        std::string s = std::to_string(unsigned(b[2]));
        if (b[1] != 0) { s += '-'; s += std::to_string(unsigned(b[1])); }
        if (b[0] != 0) { s += '-'; s += std::to_string(unsigned(b[0])); }
        return s;
    }

    /// строковое представление имени, типа "2.1-10-1"
    std::string to_string() const
    {
        std::string s = std::to_string(sprout());
        s += '.';
        s += flower();
        return s;
    }

    /// пустое имя "", которое содержит одни 0
    bool empty() const { return u == 0; }

    /// представление имени в виде интегрального (целого) типа
    unsigned to_integral() const { return u; }

    /// возвращает порядок оси, на которой располагается цветок
    unsigned axis_order() const
    {
        if (b[0]) return 3;
        if (b[1]) return 2;
        if (b[2]) return 1;
        return 0;
    }

    /// возвращает номер объекта на своей оси
    unsigned axis_position() const
    {
        if (b[0]) return unsigned(b[0]);
        if (b[1]) return unsigned(b[1]);
        if (b[2]) return unsigned(b[2]);
        return 0;
    }

    /// получение номера верхушки для данного цветка
    Flowerstamp parent(size_t offset=0) const
    {
        if (offset == 0)
        {
            if (b[0] != 0) return Flowerstamp(unsigned(b[3]), unsigned(b[2]), unsigned(b[1]), 0);
            if (b[1] != 0) return Flowerstamp(unsigned(b[3]), unsigned(b[2]), 0, 0);
            if (b[2] != 0) return Flowerstamp(unsigned(b[3]), 0, 0, 0);
        }

        if (offset == 1)
        {
            if (b[0] != 0) return Flowerstamp(unsigned(b[3]), unsigned(b[2]), 0, 0);
            if (b[1] != 0) return Flowerstamp(unsigned(b[3]), 0, 0, 0);
        }

        if (offset == 2)
        {
            if (b[0] != 0) return Flowerstamp(unsigned(b[3]), 0, 0, 0);
        }
        return Flowerstamp(0, 0, 0, 0);
    }

    /// получение номера верхушки заданного порядка для данного цветка
    Flowerstamp up(signed order=0) const
    {
        if (order == -1) return Flowerstamp(0, 0, 0, 0);
        if (order == 0) return Flowerstamp(unsigned(b[3]), 0, 0, 0);
        if (order == 1) return Flowerstamp(unsigned(b[3]), unsigned(b[2]), 0, 0);
        if (order == 2) return Flowerstamp(unsigned(b[3]), unsigned(b[2]), unsigned(b[1]), 0);
        return *this;
    }

    /// используется для получения верхней недействительной границы пары
    Flowerstamp operator++(int) const { return Flowerstamp(u + 0x01010101); }
};

/*!
 * \brief дает уровень цветка
 *   1 - цветoк 1-го порядка
 *   2 - цветoк 2-го порядка
 *   3 - цветoк 3-го порядка
 */
inline int rank(const Flowerstamp &a)
{
    if (a.get(3) != 0) return 3;
    if (a.get(2) != 0) return 2;
    if (a.get(1) != 0) return 1;
    return 0;
}

/*!
 * \brief дает уровень совпадения цветков
 *  -1 - цветки на разных главных побегах или ошибка
 *   0 - оба цветка лежат на одном главном побеге (совпадают sprouts)
 *   1 - оба цветка лежат на одном побеге 1-го порядка
 *   2 - оба цветка лежат на одном побеге 2-го порядка
 */
inline int rank(const Flowerstamp &a, const Flowerstamp &b)
{
    if (a.get(0) != b.get(0)) return -1;
    if (a.get(1) != b.get(1)) return 0;
    if (a.get(2) != b.get(2)) return 1;
    if (a.get(3) != b.get(3)) return 2;
    return -1;
}


inline bool operator<(const Flowerstamp &a, const Flowerstamp &b)
{ return a.to_integral() < b.to_integral(); }

inline bool operator<=(const Flowerstamp &a, const Flowerstamp &b)
{ return a.to_integral() <= b.to_integral(); }

inline bool operator!=(const Flowerstamp &a, const Flowerstamp &b)
{ return a.to_integral() != b.to_integral(); }

inline bool operator==(const Flowerstamp &a, const Flowerstamp &b)
{ return a.to_integral() == b.to_integral(); }

inline std::ostream &operator <<(std::ostream & os, const Flowerstamp &s)
{
    os << s.to_string().c_str();
    return os;
}

/*!
 * \brief структура хранения записи БД (dome) в памяти
 *  База данных хранится в стандартной для ML форме, то есть матрицы SAMPLES x ATRIBUTES.
 *  Цель такого хранения - обеспечить выборку данных по любой совокупности атрибутов.
 */
struct flower_t
{
    Flowerstamp fsp; // flowerstamp цветка
    Timestamp tsp; // timestamp текущей даты
    Phasestamp psp; // phasestamp фазы развития

    flower_t(const Flowerstamp &fsp=Flowerstamp(), const Timestamp &tsp=Timestamp(),
        const Phasestamp &psp=Phasestamp())
    : fsp(fsp), tsp(tsp), psp(psp) {}
};

/// оператор упорядочивания записей в БД для поиска по запросу
inline bool operator<(const flower_t &a, const flower_t &b)
{
    if (a.fsp != b.fsp) return a.fsp < b.fsp;
    return a.tsp < b.tsp;
}

inline std::ostream &operator <<(std::ostream & os, const flower_t &s)
{
    printf("{ %s   %s   %s   }", s.fsp.to_string().c_str(), s.tsp.to_string().c_str(),
        s.psp.to_string().c_str());
    return os;
}

#endif // RECORD_FLOWER_H
