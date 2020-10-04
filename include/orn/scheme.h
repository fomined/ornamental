#ifndef ORN_SCHEME_H
#define ORN_SCHEME_H

#include <string>
#include <vector>
#include <string>
#include <set>
#include <tuple>
#include <memory>

#include "defs.h"
#include "container/vector.h"
#include "math/math.h"
#include "math/statistics.h"

#include "os/os.h"
#include "svg/primitive.h"

#include "orn/orn.defs.h"
#include "orn/time.h"
#include "orn/phase.h"
#include "orn/flower.h"

namespace prg
{
    const double sqrt_2 = std::sqrt(2);
    const double _1sqrt_2 = std::sqrt(0.5);
    const double sqrt_pi = std::sqrt(M_PI);
    const double sqrt_2pi = std::sqrt(2 * M_PI);

    const double RADIUS = 0.2;
    const double RADIUS_GAP = 0.05;

    const double SCALE_COEF = 30;
    const double PHASE_RADIUS = SCALE_COEF * (RADIUS - RADIUS_GAP);
    const double FLOWER_WIDTH = SCALE_COEF * (2 * RADIUS);
    const double FLOWER_HEIGHT = FLOWER_WIDTH * sqrt_2;

    /// запись для пересылки в latex для отрисовки части объекта
    struct ShemeRecord
    {
        double x1, y1, x2, y2; // координаты линии
        Signed phase; // фаза "цветка"
        Signed angle; // угол поворота цветка в градусах
        bool flag; // рисовать или нет в (x2,y2) цветок

        /// вывод в объекта строку
        std::string to_string() const;
    };

    inline std::string ShemeRecord::to_string() const
    {
        const char *format = "%5.2lf, %5.2lf, %5.2lf, %5.2lf, %1ld, %3ld, %1ld";
        char buffer[BUF_LEN];
        std::sprintf(buffer, format, x1, y1, x2, y2, phase, angle, flag);
        return std::string(buffer);
    }

    inline std::ostream &operator <<(std::ostream & os, const ShemeRecord &s)
        { return os << s.to_string(); }

    /*!
     * \brief Структура, описывающая схему соцветия побега для его визуализации
     */
    struct SproutScheme
    {
        double width = 0;       ///< ширина области, занимаемой "побегом"
        double heigth = 0;      ///< высота области, занимаемой "побегом"
        Phase phase = UNDEF;    ///< фаза для терминала (иначе UNDEF)
        Signed angle = 0;          ///< угол поворота (в числе 45 градусных сегментов)

        std::vector<SproutScheme> children; // побочные побеги, выросшие внутри данного
            // если 0, то это завершающий побег без потомков
        std::vector<double> children_d; // координата привязки побочных побегов вдоль линии основного
            // прочие требования на размещение исходят из следующего соглашений
            // 1) нечетные побочные побеги размещаются слева от корневого, четные - справа,
            // 2) угол размещения всегда 45 градусов
            // 3) координаты каждого элемента вычисляются после полного построения главного побега
            //  со всеми его потомками, при этом вычисление идет по пути от корня к потомкам

        SproutScheme(double width=0, double heigth=0, Phase phase=UNDEF,
            Signed angle=0) : width(width), heigth(heigth), phase(phase), angle(angle) {}

        const SproutScheme &child(Unsigned i) const { return children[i]; }
        double get_distance(Unsigned i) const { return children_d[i]; }
        Unsigned size() const { return children.size(); }
        bool empty() const { return children.empty(); }

        /// получить данные по сгенерированной схеме
        std::vector<ShemeRecord> data() const;

        /*!
         * \brief конструирование побега из его побочных побегов (составные элементы разные)
         * \param s набор сортированных побочных побегов (с 1 до конца)
         */
        static SproutScheme construct(const SproutScheme &main, const std::vector<SproutScheme> &s);
    };

    inline SproutScheme SproutScheme::construct(const SproutScheme &main, const std::vector<SproutScheme> &s)
    {
        SproutScheme sprout;

        sprout.phase = main.phase;
        sprout.angle = main.angle;

        Unsigned n = s.size();
        sprout.children = s;
        sprout.children_d.resize(n);

        double hcur = 0; // текущая высота побега
        double dh_next = sprout.children[n - 1].width / sqrt_2; // прогноз на следующую высоту побега
            // поскольку ширины внутренних побегов могут различаться, что их высоты крепления
            // будут случайным образом меняться. Чтобы не произошло перепутывания побегов по порядку,
            // прикрепление нового побега делается на высоту не меньше, чем прогноз (обусловленный
            // шириной последнего прикрепленного побега)

        sprout.width = 0;
        sprout.heigth = 0;

        for (Unsigned i=0; i<n; i++)
        {
            double w = sprout.children[n - 1 - i].width;
            double h = sprout.children[n - 1 - i].heigth;

            double dh_cur = w / sqrt_2; // изменение высоты крепления
            hcur += std::max(dh_cur, dh_next);
            sprout.children_d[n - 1 - i] = hcur;
            dh_next = dh_cur;

            sprout.width = std::max(sprout.width, h * sqrt_2);
            sprout.heigth =  std::max(sprout.heigth, hcur + h / sqrt_2);
        }
        return sprout;
    }

    inline std::vector<ShemeRecord> SproutScheme::data() const
    {
        struct record_t
        {
            const SproutScheme *sprout;
            double x, y; // начало побега в абсолютных координатах
            Signed iangle; // угол поворота побега в количестве 45 градусных поворотов
            record_t(const SproutScheme *sprout, double x=0, double y=0, Signed iangle=0)
                : sprout(sprout), x(x), y(y), iangle(iangle) {}
        };

        std::vector<ShemeRecord> flower_draw;

        std::stack<record_t> stack;
        stack.push(record_t(this, 0, 0, 0));

        double xangle_coef[8] = { 0, _1sqrt_2, 1, _1sqrt_2, 0, -_1sqrt_2, -1, -_1sqrt_2 };
        double yangle_coef[8] = { 1, _1sqrt_2, 0, -_1sqrt_2, -1, -_1sqrt_2, 0, _1sqrt_2 };

        while (!stack.empty())
        {
            record_t rec = stack.top();
            stack.pop();

            ShemeRecord f;

            //-------------------------------------------------------------------------------
            //                          обработка нетерминальных узлов
            //-------------------------------------------------------------------------------
            const SproutScheme *sprout = rec.sprout;
            if (!sprout->empty())
            {
                f.x1 = rec.x;
                f.y1 = rec.y;

                double xcoef = xangle_coef[mod(rec.iangle, Signed(8))];
                double ycoef = yangle_coef[mod(rec.iangle, Signed(8))];

                double h = sprout->heigth;
                f.x2 = rec.x + h * xcoef;
                f.y2 = rec.y + h * ycoef;

                f.flag = (sprout->phase == UNDEF) ? false : true;
                f.phase = to_signed(sprout->phase);
                f.angle = mod(2 - rec.iangle, Signed(8)) * 45;

                flower_draw.push_back(f);
                    // вставка побегов без "цветков"

                Unsigned n = sprout->size();
                for (Unsigned i=0; i<n; i++)
                {
                    double xy = sprout->get_distance(i);
                    double x = rec.x + xy * xcoef;
                    double y = rec.y + xy * ycoef;
                    Signed kangle = rec.iangle + ((i % 2 == 0) ? -1 : 1);

                    stack.push(record_t(&sprout->child(i), x, y, kangle));
                }
                continue;
            }
            //-------------------------------------------------------------------------------
            //                          обработка терминальных узлов
            //-------------------------------------------------------------------------------
            double xcoef = xangle_coef[mod(rec.iangle, Signed(8))];
            double ycoef = yangle_coef[mod(rec.iangle, Signed(8))];

            f.x1 = rec.x;
            f.y1 = rec.y;

            double h = 0.2 * sprout->heigth;
                // коэффициент у sprout->heigth определяет длину терминального побега до органа

            f.x2 = rec.x + h * xcoef;
            f.y2 = rec.y + h * ycoef;

            f.flag = true;
            f.angle = mod(2 - rec.iangle, Signed(8)) * 45;
            f.phase = to_signed(sprout->phase);

            flower_draw.push_back(f);
        }
        return flower_draw;
    }

}
#endif // ORN_SCHEME_H
