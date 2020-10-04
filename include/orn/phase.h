#ifndef ORN_PHASE_H
#define ORN_PHASE_H

#include <cstdint>
#include <string>
#include <vector>

#include "container/vector.h"
#include "orn/orn.defs.h"
#include "svg/primitive.h"

/*!
 * \file Основные структуры и функции для работы с фазами цветков
 */

namespace prg
{

    enum Phase { UNDEF=-1, BUD=0, FLOWER, OVARY, SEED, DEATH, INVALID };
        // UNDEF можно считать фазой отрастания побега (но лучше исключить)
        // BUD=0 для отсутствия пустого места в массивах по началам фаз

    const char * const phase_name[] = {"", "bud", "flower", "ovary", "seed", "death", "?"};
    constexpr Unsigned PHASE_CNT = 5; // "bud", "flower", "ovary", "seed" , "death"

    /// получить полное имя фазы
    inline std::string to_string(Phase e) { return phase_name[value(e) + 1]; } // список имен смещен (+1)

    /// представление фазы в виде индекса (есть отрицательный индекс)
    inline Signed to_signed(Phase e) { return value(e); }

    /// сконструировать объект из индекса
    inline Phase make_phase(Signed i) { Phase p; p = static_cast<Phase>(i); return p; }

    /// сконструировать объект из (короткого или полного) имени
    inline Phase make_phase(const std::string &s)
    {
        Phase phase(INVALID);  // если неопознан, то INVALID
        Unsigned cnt = sizeof(phase_name) / sizeof(char*); // учет всех значений
        for (Unsigned i=0; i<cnt; i++)
        {
            // s.length() и однозначность имен допускает сокращения имен
            if (std::strncmp(s.c_str(), phase_name[i], s.length()) == 0)
            {
                phase = static_cast<Phase>(i-1); // i-1, так как нумерация смещена
                break;
            }
        }
        return phase;
    }

    /// конструирование фазы по строке с учетом прежней фазы
    inline Phase make_phase(Phase prev, const std::string &s)
    {
        std::string q(prg::trim(s));
        return (q.empty()) ?  prev : make_phase(q);
    }

    #define MAKE_OPER(op) inline bool operator op(Phase a, Phase b) { return to_signed(a) op to_signed(b); }
        MAKE_OPER(<)
        MAKE_OPER(<=)
        MAKE_OPER(==)
        MAKE_OPER(!=)
    #undef MAKE_OPER

    /// "расстояние" между фазами, фазы пропущены, если > 1
    inline Unsigned distance(Phase a, Phase b) { return to_unsigned(b) - to_unsigned(a); }

    /// конструирование смещенной фазы
    template <typename T>
    inline Phase operator+(Phase a, T i) { return make_phase(to_unsigned(a) + i); }

    template <typename T>
    inline Phase operator-(Phase a, T i) { return operator+(a, -i); }

    /// вывод на печать
    inline std::ostream &operator <<(std::ostream & os, const Phase &s) { return os << to_string(s); }


    /*!
     * \brief create scheme of flower in different phase of development
     * \param phase phase of development (BUD, FLOWER, OVARY)
     * \return element of svg picture
     */
    inline shared<Path> flower_scheme(Phase phase, double scale=1.)
    {
        shared<Path> pvg(nullptr);
        if (phase == BUD)
        {
            pvg = shared<Path>(new Path(Color::Transparent, Stroke(.5, Color::Blue)));
            pvg->push_back('M',   0,   0);
            pvg->push_back('Q', -20,  -20);
            pvg->push_back(' ',   0,  -60);
            pvg->push_back('L',   0, -30);
            pvg->push_back('L',   0, -60);
            pvg->push_back('Q',  20,  -20);
            pvg->push_back(' ',   0,   0);
            pvg->scale(scale);
            pvg->rotate(180);
        }
        else if (phase == FLOWER)
        {
            pvg = shared<Path>(new Path(Color::Blue, Stroke(.5, Color::Red)));
            pvg->push_back('M',   0,   0);
            pvg->push_back('C', -25,-60);
            pvg->push_back(' ', -10,-110);
            pvg->push_back(' ', -40,-125);
            pvg->push_back('Q', -25,-125);
            pvg->push_back(' ', -10,-100);
            pvg->push_back('Q', -10,-120);
            pvg->push_back(' ',  0,-130);
            pvg->push_back('Q', 10,-120);
            pvg->push_back(' ', 10,-100);
            pvg->push_back('Q', 25,-125);
            pvg->push_back(' ', 40,-125);
            pvg->push_back('C',  10, -110);
            pvg->push_back(' ',  25, -60);
            pvg->push_back(' ', 0, 0);
            pvg->scale(scale);
            pvg->rotate(180);
        }
        else if (phase == OVARY)
        {
            pvg = shared<Path>(new Path(Color::Green, Stroke(.5, Color::Blue)));
            pvg->push_back('M',   0,   0);
            pvg->push_back('Q', -20,   0);
            pvg->push_back(' ', -20, -60);
            pvg->push_back('Q',   0, -50);
            pvg->push_back(' ',  20, -60);
            pvg->push_back('Q',  20,   0);
            pvg->push_back(' ',   0,   0);
            pvg->scale(scale);
            pvg->rotate(180);
        }
        else if (phase == SEED)
        {
            pvg = shared<Path>(new Path(Color::Brown, Stroke(.5, Color::Blue)));
            pvg->push_back('M',   0,   0);
            pvg->push_back('Q', -20,   0);
            pvg->push_back(' ', -20, -60);
            pvg->push_back('Q',   0, -50);
            pvg->push_back(' ',  20, -60);
            pvg->push_back('Q',  20,   0);
            pvg->push_back(' ',   0,   0);
            pvg->scale(scale);
            pvg->rotate(180);
        }
        return pvg;
    }

}
#endif // ORN_PHASE_H
