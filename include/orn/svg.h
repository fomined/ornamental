#ifndef ORN_SVG_H
#define ORN_SVG_H

#include <string>
#include <vector>
#include <string>
#include <set>
#include <tuple>
#include <memory>

#include "svg/primitive.h"
#include "orn/orn.defs.h"
#include "orn/time.h"
#include "orn/phase.h"
#include "orn/flower.h"
#include "orn/scheme.h"

namespace prg
{
    class svgScheme : public Shape
    {
        std::vector<shared<Line>> lines;
        std::vector<shared<Path>> flowers;

    public:

        svgScheme() : Shape() {}
        svgScheme(const std::vector<ShemeRecord> &s) : Shape() { *this << s; }

        svgScheme &operator<< (const std::vector<ShemeRecord> &s)
        {
            for (auto e : s)
            {
                shared<Line> line(new Line(Point(e.x1, e.y1), Point(e.x2, e.y2), Stroke(1, Color::Black)));
                lines.push_back(line);

                Phase phase = make_phase(e.phase);
                if (e.phase == BUD || e.phase == FLOWER || e.phase == OVARY || e.phase == SEED)
                {
                    auto f = flower_scheme(phase, 0.08);
                    f->rotate(e.angle-90);
                    f->mirror();
                    f->move(Point(e.x2, e.y2));
                    flowers.push_back(f);
                }
            }
            return *this;
        }

        Rectangle rectange() const
        {
            Rectangle r = lines[0]->rectange();
            for (auto e : lines) r = r + e->rectange();
                // for (auto e : flowers) r = r + e->rectange();
                // ограничим размеры только "ветками", чтобы видео "не прыгало"
            return r;
        }

        std::string toString(Layout const & layout) const
        {
            std::string ss;
            for (auto e : lines) ss += e->toString(layout);
            for (auto e : flowers) ss += e->toString(layout);
            return ss;
        }
        void offset(Point const &) {}

        void move(const Point &p)
        {
            for (auto &e : lines) e->move(p);
            for (auto &e : flowers) e->move(p);
        }

        /// центрировать по области полотна
        void center(Layout const & layout)
        {
            Dimensions dim = layout.dimensions;
            double lxc = 0.5 * dim.width;
            double lyc = 0.5 * dim.height;

            double xc = 0.5 * (lines[0]->start_point.x + lines[0]->end_point.x);
            double yc = 0.5 * (lines[0]->start_point.y + lines[0]->end_point.y);

            Point translation(lxc - xc, lyc - yc);
            move(translation);
        }

        /// масштабировать по размеру полотна
        void scale(Layout const & layout)
        {
            Dimensions dim = layout.dimensions;
            double width = dim.width;
            double height = dim.height;

            Rectangle ra = this->rectange();
            double xscale = 0.8 * width / ra.width;
            double yscale = 0.8 * height / ra.height;
            double scale = std::min(xscale, yscale);

            for (auto e : lines) e->scale(scale);
            for (auto e : flowers) e->scale(scale);
        }
    };

}
#endif // ORN_SVG_H
