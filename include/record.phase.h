#ifndef RECORD_PHASE_H
#define RECORD_PHASE_H

#include <cstdint>
#include <string>
#include <vector>

#include "prg.string.h"
#include "defs.h"

/*!
 * \file Основные структуры и функции для работы с фазами цветков
 */

enum phase_enum_t { UNDEF = 0, BUD, FLOWER, OVARY, SEED, DEATH, INVALID};
const char * const phase_name[] = {"", "bud", "flower", "ovary", "seed", "death", "?"};
constexpr unsigned PHASE_CNT = sizeof(phase_name) / sizeof(char *);

/// получить полное имя фазы
inline std::string phasename(phase_enum_t e) { return std::string(phase_name[e]); }

/*!
 * \brief класс является номером записи в масиве описателей фаз PHASES
 */
struct Phasestamp
{
    phase_enum_t phase;

    Phasestamp(phase_enum_t p=UNDEF) : phase(p) {}

    /// Получить фазу по имени или сокращению имени
    Phasestamp(const std::string &s) : phase(INVALID)
    {
        for (size_t i=0; i<PHASE_CNT; i++)
            if (strncmp(s.c_str(), phase_name[i], s.length()) == 0)
                // использование s.length() и однозначности имен
                // позволяет передавать в функцию сокращения имен
                { phase = static_cast<phase_enum_t>(i); break; }
    }
    Phasestamp &operator++()
    {
        if (phase == UNDEF || phase == INVALID) EXIT;
        auto p = static_cast<typename std::underlying_type<phase_enum_t>::type>(phase);
        phase = static_cast<phase_enum_t>(++p);
        return *this;
    }

    /// пустая фаза означает, что такая же как предыдущая
    bool empty() const { return phase == UNDEF; }

    std::string to_string() const { return std::string(phasename(phase)); }
    bool valid() const { return phase != INVALID; }
};

inline bool operator<(const Phasestamp &a, const Phasestamp &b)
    { if (a.phase == UNDEF || b.phase == INVALID) EXIT; return a.phase < b.phase; }

inline bool operator<=(const Phasestamp &a, const Phasestamp &b)
    { if (a.phase == UNDEF || b.phase == INVALID) EXIT; return a.phase <= b.phase; }

inline bool operator!=(const Phasestamp &a, const Phasestamp &b)
    { if (a.phase == UNDEF || b.phase == INVALID)
    EXIT;
    return a.phase != b.phase; }

inline bool operator==(const Phasestamp &a, const Phasestamp &b)
{
    if (a.phase == UNDEF || b.phase == INVALID)
        EXIT;
    return a.phase == b.phase;
}

inline std::ostream &operator <<(std::ostream & os, const Phasestamp &s)
    { return os << s.to_string(); }

/*!
 * \brief нахождение числового идентификтора фазы по строковому идентификатору фазы
 * \param s идентификатор фазы (пустой - возврат предыдущей фазы)
 */
inline Phasestamp make_phasestamp(Phasestamp prev, const std::string &s)
{
    std::string q(prg::trim(s));
    if (q.empty()) return prev;
    return Phasestamp(q);
}

#endif // RECORD_PHASE_H
