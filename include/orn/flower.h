#ifndef ORN_FLOWER_H
#define ORN_FLOWER_H

#include <cstdint>
#include <string>
#include <vector>
#include "orn/orn.defs.h"

namespace prg
{

    //-------------------------------------------------------------------------------------------------------
    //                                      FlowerIndex
    //-------------------------------------------------------------------------------------------------------
    class FlowerIndex
    {
        // Порядок отнесения цветков к осям определяется не нумерацией цветков (она дефектна
        // по основному терминалу побега, который не совпадает с принципом нумерации остальных
        // терминалов), этот порядок определяется только наличием/отсутствием дополнительных
        // индексов. То есть, все цветки типа {k, 0} есть цветки 1-го порядка, независимо
        // от значения k!=0, а цветки типа {k, l, 0} - 2-го и т.д. Сама нумерация в пределах
        // дополнительного индекса БЕЗРАЗЛИЧНА, она должна быть только упорядоченной.

        enum { n = sizeof(unsigned) }; // размерность мультииндекса (с исключением номера побега)
        union
        {
            Unsigned U;
            struct
            {
                unsigned char b[4]; // мультииндекс цветка
                char a; // тип записи (0 [структура], 1 [динамика])
                char _; // reserved
                unsigned short u; // номер побега
            };
        };

    public:

        FlowerIndex() : U(0) {}

        FlowerIndex(Unsigned U) : U(U) {}
        FlowerIndex(unsigned sprout, bool flag, const std::string &s) : U(0)
        {
            u = sprout; a = flag; if (s.empty()) return;

            Unsigned pos = 0;
            b[3] = char(atoi(&s[pos]));
            for (Signed i=2; i>=0; i--)
            {
                pos = s.find_first_of('-', pos);
                if (pos == std::string::npos) break;
                b[i] = char(atoi(&s[++pos]));
            }
        }

        /*!
         * \brief FlowerIndex
         * \param sprout номер побега
         * \param flag тип записи (0 - структура, 1- динамика)
         * \param f (заселенность осей [структура], номера цветка[динамика])
         */
        FlowerIndex(unsigned sprout, bool flag, const std::array<unsigned, 4> &f={0, 0, 0, 0})
        : U(0) { u = sprout; a = flag; b[3] = char(f[0]); b[2] = char(f[1]); b[1] = char(f[2]); b[0] = char(f[3]); }

        /// установка побега и типа записи
        FlowerIndex &operator=(const std::pair<unsigned, bool> &v) { u = v.first; a = v.second; return *this; }

        /// установка номеров осей цветка
        FlowerIndex &operator=(const std::array<unsigned, 4> &f)
            { b[3] = char(f[0]); b[2] = char(f[1]); b[1] = char(f[2]); b[0] = char(f[3]); return *this; }

        /// получение номера побега
        unsigned sprout() const { return u; }
        unsigned short &sprout() { return u; }

        /// получение типа записи
        bool dynamic() const { return a == 1; }
        bool flag() const { return a; }
        char &flag() { return a; }

        /// получение номеров осей
        std::array<unsigned, 4> axes() const { return {(unsigned)b[3], (unsigned)b[2], (unsigned)b[1], (unsigned)b[0]}; }

        /// пустое имя "", которое содержит одни 0
        bool empty() const { return U == 0; }

        /// представление имени в виде интегрального (целого) типа для упорядочения
        Unsigned to_integral() const { return U; }

        /// имя, содержащее только номера осей
        std::string to_short_string() const
        {
            std::string s;
            {
                if (b[n-1]) s += std::to_string(b[n-1]);
                for (Unsigned i=n-2; (Signed)i>=0 && b[i]!=0; i--)
                    s += '-' + std::to_string(b[i]);
            }
            return s;
        }

        /// строковое представление имени, типа "2.1-10-1"
        std::string to_string() const
        {
            std::string s = std::to_string(sprout());
            s += '.';
            if (!flag()) s+= '[';
            {
                if (b[n-1]) s += std::to_string(b[n-1]);
                for (Unsigned i=n-2; (Signed)i>=0 && b[i]!=0; i--)
                    s += '-' + std::to_string(b[i]);
            }
            if (!flag()) s+= ']';
            return s;
        }

        /// доступ к отдельным индексам [1..4]
        unsigned char operator[](Unsigned i) const { return b[n-i]; }
        unsigned char &operator[](Unsigned i) { return b[n-i]; }

        /// возвращает порядок оси, на которой располагается цветок
        Unsigned axis() const
            { Unsigned i = 0; while (i < 4 && b[i] == 0) i++; return n - i; }

        /// возвращает номер объекта на своей оси
        Unsigned position() const { Unsigned i = 0; while (b[i] == 0) i++; return b[i]; }

        /*!
         * \brief получение идентификатора номера верхушки для данного цветка
         * \note чтобы получить мультииндекс родителя нужно сделать присваивание FlowerIndex f = id;
         */
        FlowerIndex parent(Unsigned offset=0) const
        {
            FlowerIndex f = this->U; // копия
            Unsigned i = axis() - offset; // родитель это тот, у кого номер оси обнулен
            for (; i<n; i++) f.b[i] = 0; // обнуляем и более высокие индексы (на это рассчитывают другие функции)
            return f;
        }
    };

    inline bool operator<(const FlowerIndex &a, const FlowerIndex &b)
        { return a.to_integral() < b.to_integral(); }

    inline bool operator<=(const FlowerIndex &a, const FlowerIndex &b)
        { return a.to_integral() <= b.to_integral(); }

    inline bool operator!=(const FlowerIndex &a, const FlowerIndex &b)
        { return a.to_integral() != b.to_integral(); }

    inline bool operator==(const FlowerIndex &a, const FlowerIndex &b)
        { return a.to_integral() == b.to_integral(); }

    inline std::ostream &operator <<(std::ostream & os, const FlowerIndex &s)
        { os << s.to_string(); return os; }


    /*!
     * \brief relative height of the flower on sprout
     * \param f flower index
     * \param n number of flowers on the 1st axis
     * \return
     */
    inline flagged_value<double> relative_height(const FlowerIndex &f, Unsigned n)
    {
        if (n <= 1 || f.flag()==false) return {0, false};
        Unsigned k = f[1]; // номер оси 1-го порядка, где расположен цветок
        return {double(n - k) / (n - 1), true};
    }

}
#endif // ORN_FLOWER_H
