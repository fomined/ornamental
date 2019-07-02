#ifndef PRG_CURVE
#define PRG_CURVE

#include <array>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <random>

#include "prg.defs.h"
#include "prg.array.h"
#include "prg.maxmin.h"
#include "prg.point.h"

namespace prg
{

    /// базисная структура "график"
    template <unsigned n> using curve = std::vector<point<n>>;

    /// сортировка данных в графике по заданной колонке (x, y, ..)
    template <unsigned power>
    inline void sort_by_column(curve<power> &g, unsigned col)
    {
        std::sort(g.begin(), g.end(),
            [col](const point<power> &a, const point<power> &b){ return a[col] < b[col]; }
        );
    }

    /// проверка наличия сортировки данных в графике по заданной колонке (x, y, ..)
    template <unsigned power>
    inline bool is_sorted_by_column(const curve<power> &g, unsigned col)
    {
        return std::is_sorted(g.begin(), g.end(),
            [col](const point<power> &a, const point<power> &b){ return a[col] < b[col]; }
        );
    }

    template <unsigned n>
    inline void printf(const curve<n> &v, const char *format="%lf ")
    {
        for (auto e : v)
        {
            for (unsigned i=0; i<n; i++) std::printf(format, e[i]);
            std::printf("\n");
        }
        std::cout.flush();
    }

    inline bool equal(double a, double b, double eps=0.001)
		{ return std::abs(a - b) < eps ? true : false; }

    inline bool less(double a, double b, double eps=0.001)
		{ return equal(a, b, eps) ? false : a < b; }

    template <unsigned n>
    inline curve<n> summarize(const curve<n> &a, const curve<n> &b, double eps=0.001)
    {
		curve<n> tmp;

		unsigned i=0, j=0;
		unsigned na = a.size(), nb = b.size();
		for (; i<na && j<nb; )
		{
			double xa = a[i][0], xb = b[j][0];
			if (less(xa, xb, eps)) { tmp.push_back(a[i++]); continue; }
			if (less(xb, xa, eps)) { tmp.push_back(b[j++]); continue; }

			point<n> p; p[0] = xa;
			for (unsigned k=1; k<n; k++) p[k] = a[i][k] + b[j][k];
			tmp.push_back(p);
			i++; j++;
		}
		if (i != na) for (; i<na; i++) tmp.push_back(a[i]);
		if (j != nb) for (; j<nb; j++) tmp.push_back(b[j]);

		return tmp;
    }

    template <unsigned n>
    inline curve<n> subtract(const curve<n> &a, const curve<n> &b, double eps=0.001)
    {
		curve<n> tmp;

		unsigned i=0, j=0;
		unsigned na = a.size(), nb = b.size();
		for (; i<na && j<nb; )
		{
			double xa = a[i][0], xb = b[j][0];
			if (less(xa, xb, eps)) { tmp.push_back(a[i++]); continue; }
			if (less(xb, xa, eps)) { tmp.push_back(negate<n>(b[j++])); continue; }

			point<n> p; p[0] = xa;
			for (unsigned k=1; k<n; k++) p[k] = a[i][k] - b[j][k];
			tmp.push_back(p);
			i++; j++;
		}
		if (i != na) for (; i<na; i++) tmp.push_back(a[i]);
		if (j != nb) for (; j<nb; j++) tmp.push_back(negate<n>(b[j]));

		return tmp;
    }

    template <unsigned n>
    inline curve<n> merge(const curve<n> &a, const curve<n> &b, unsigned icol=0)
    {
        curve<n> tmp(a);
        tmp.insert(tmp.end(), b.begin(), b.end());
        std::sort(tmp.begin(), tmp.end(),
            [icol](const point<n> &a, const point<n> &b)
            { return a[icol] < b[icol]; }
        );
        return tmp;
    }

    template <unsigned n>
    inline curve<n> lower_bound(const curve<n> &a, unsigned ix, unsigned iy, double eps=0.001)
    {
		curve<n> tmp(a);
		std::sort(tmp.begin(), tmp.end(),
            [ix,iy,eps](const point<n> &a, const point<n> &b)
            {
                if (equal(a[ix], b[ix], eps)) return a[iy] < b[iy];
                return a[ix] < b[ix];
            }
		);

        auto it = std::unique(tmp.begin(), tmp.end(),
            [ix,eps](const point<n> &a, const point<n> &b){ return equal(a[ix], b[ix], eps); }
        );
        tmp.erase(it, tmp.end());

 		return tmp;
    }

    template <unsigned n>
    inline curve<n> upper_bound(const curve<n> &a, unsigned ix, unsigned iy, double eps=0.001)
    {
		curve<n> tmp(a);
		std::sort(tmp.begin(), tmp.end(),
            [ix,iy,eps](const point<n> &a, const point<n> &b)
            {
                if (equal(a[ix], b[ix], eps)) return a[iy] > b[iy];
                return a[ix] < b[ix];
            }
		);

        auto it = std::unique(tmp.begin(), tmp.end(),
            [ix,eps](const point<n> &a, const point<n> &b){ return equal(a[ix], b[ix], eps); }
        );
        tmp.erase(it, tmp.end());

 		return tmp;
    }

    template <unsigned n, typename Iterator>
    inline curve<n> select_ynear(double ytarget, unsigned count, Iterator first, Iterator last,
        unsigned ix, unsigned iy, double eps=0.001)
    {
        point<n> concentrator_target;
        concentrator_target[iy] = ytarget;

        prg::concentrator<point<n>, void> sum(concentrator_target, count);
        for (auto it=first; it!=last; ++it) sum.insert(*it,
            [iy](const point<n> &a, const point<n> &b){ return std::abs(a[iy] - b[iy]); }
        );

        curve<n> tmp;
        for (auto e : sum) tmp.push_back(e.second);

		std::sort(tmp.begin(), tmp.end(),
            [ix, eps](const point<n> &a, const point<n> &b) { return a[ix] < b[ix]; }
		);

 		return tmp;
    }

    template <unsigned n, typename Iterator>
    inline curve<n> select_ynear(double ytarget, double portion, Iterator first, Iterator last,
        unsigned ix, unsigned iy, double eps=0.001)
    {
        unsigned count = unsigned(std::ceil(portion * std::distance(first, last)));
        return select_ynear<n>(ytarget, count, first, last, ix, iy, eps);
    }

    template <unsigned n, typename Iterator>
    inline curve<n> lower_bound(int direction, Iterator first, Iterator last, double portion=0.5,
        unsigned ix=0, unsigned iy=1, double eps=0.001)
    {
        if (first == last) return curve<n>();

        double ymin = infinity<double>();
        double xmin = infinity<double>(), xmax = -infinity<double>();

        for (auto it=first; it!=last; ++it)
        {
            if (ymin > (*it)[iy]) ymin = (*it)[iy];
            if (xmin > (*it)[ix]) xmin = (*it)[ix];
            if (xmax < (*it)[ix]) xmax = (*it)[ix];
        }
        auto f = select_ynear<n>(ymin, portion, first, last, ix, iy, eps);

        // вставляем пропущенные отбором точки начала кривой
        if (direction == -1 && xmin < f.front()[ix])
        {
            double xlim = f.front()[ix]; // граница вставки данных
            auto it = std::lower_bound(first, last, xlim,
                [ix](const point<n> &b, double a){ return b[ix] < a; }
            );

            curve<n> tmp = prg::lower_bound<n>(direction, first, it, portion, ix, iy, eps);
            tmp.insert(tmp.end(), f.begin(), f.end());
            f.swap(tmp);
        }

        if (direction == 1 && xmax > f.back()[ix])
        {
            double xlim = f.back()[ix]; // граница вставки данных
            auto it = std::upper_bound(first, last, xlim,
                [ix](double a, const point<n> &b){ return a < b[ix]; }
            );

            curve<n> tmp = prg::lower_bound<n>(direction, it, last, portion, ix, iy, eps);
            f.insert(f.end(), tmp.begin(), tmp.end());
        }
        assert(prg::is_sorted_by_column<n>(f, 0));

        return f;
    }

    template <unsigned n, typename Iterator>
    inline curve<n> lower_bound(Iterator first, Iterator last, double portion=0.5,
        unsigned ix=0, unsigned iy=1, double eps=0.001)
    {
        if (first == last) return curve<n>();

        double ymin = infinity<double>();
        double xmin = infinity<double>(), xmax = -infinity<double>();

        for (auto it=first; it!=last; ++it)
        {
            if (ymin > (*it)[iy]) ymin = (*it)[iy];
            if (xmin > (*it)[ix]) xmin = (*it)[ix];
            if (xmax < (*it)[ix]) xmax = (*it)[ix];
        }
        auto f = select_ynear<n>(ymin, portion, first, last, ix, iy, eps);

        // вставляем пропущенные отбором точки начала кривой
        if (xmin < f.front()[ix])
        {
            double xlim = f.front()[ix]; // граница вставки данных
            auto it = std::lower_bound(first, last, xlim,
                [ix](const point<n> &b, double a){ return b[ix] < a; }
            );

            curve<n> tmp = prg::lower_bound<n>(-1, first, it, portion, ix, iy, eps);
            tmp.insert(tmp.end(), f.begin(), f.end());
            f.swap(tmp);
        }

        if (xmax > f.back()[ix])
        {
            double xlim = f.back()[ix]; // граница вставки данных
            auto it = std::upper_bound(first, last, xlim,
                [ix](double a, const point<n> &b){ return a < b[ix]; }
            );

            curve<n> tmp = prg::lower_bound<n>(1, it, last, portion, ix, iy, eps);
            f.insert(f.end(), tmp.begin(), tmp.end());
        }
        assert(prg::is_sorted_by_column<n>(f, 0));

        return f;
    }

    template <unsigned n, typename Iterator>
    inline curve<n> upper_bound(int direction, Iterator first, Iterator last, double portion=0.5,
        unsigned ix=0, unsigned iy=1, double eps=0.001)
    {
        if (first == last) return curve<n>();

        double ymax = -infinity<double>();
        double xmin = infinity<double>(), xmax = -infinity<double>();

        for (auto it=first; it!=last; ++it)
        {
            if (ymax < (*it)[iy]) ymax = (*it)[iy];
            if (xmin > (*it)[ix]) xmin = (*it)[ix];
            if (xmax < (*it)[ix]) xmax = (*it)[ix];
        }
        auto f = select_ynear<n>(ymax, portion, first, last, ix, iy, eps);

        // вставляем пропущенные отбором точки начала кривой
        if (direction == -1 && xmin < f.front()[ix])
        {
            double xlim = f.front()[ix]; // граница вставки данных
            auto it = std::lower_bound(first, last, xlim,
                [ix](const point<n> &b, double a){ return b[ix] < a; }
            );

            curve<n> tmp = prg::upper_bound<n>(direction, first, it, portion, ix, iy, eps);
            tmp.insert(tmp.end(), f.begin(), f.end());
            f.swap(tmp);
        }

        if (direction == 1 && xmax > f.back()[ix])
        {
            double xlim = f.back()[ix]; // граница вставки данных
            auto it = std::upper_bound(first, last, xlim,
                [ix](double a, const point<n> &b){ return a < b[ix]; }
            );

            curve<n> tmp = prg::upper_bound<n>(direction, it, last, portion, ix, iy, eps);
            f.insert(f.end(), tmp.begin(), tmp.end());
        }
        assert(prg::is_sorted_by_column<n>(f, 0));

        return f;
    }

    template <unsigned n, typename Iterator>
    inline curve<n> upper_bound(Iterator first, Iterator last, double portion=0.5,
        unsigned ix=0, unsigned iy=1, double eps=0.001)
    {
        if (first == last) return curve<n>();

        double ymax = -infinity<double>();
        double xmin = infinity<double>(), xmax = -infinity<double>();

        for (auto it=first; it!=last; ++it)
        {
            if (ymax < (*it)[iy]) ymax = (*it)[iy];
            if (xmin > (*it)[ix]) xmin = (*it)[ix];
            if (xmax < (*it)[ix]) xmax = (*it)[ix];
        }
        auto f = select_ynear<n>(ymax, portion, first, last, ix, iy, eps);

        // вставляем пропущенные отбором точки начала кривой
        if (xmin < f.front()[ix])
        {
            double xlim = f.front()[ix]; // граница вставки данных
            auto it = std::lower_bound(first, last, xlim,
                [ix](const point<n> &b, double a){ return b[ix] < a; }
            );

            curve<n> tmp = prg::upper_bound<n>(-1, first, it, portion, ix, iy, eps);
            tmp.insert(tmp.end(), f.begin(), f.end());
            f.swap(tmp);
        }

        if (xmax > f.back()[ix])
        {
            double xlim = f.back()[ix]; // граница вставки данных
            auto it = std::upper_bound(first, last, xlim,
                [ix](double a, const point<n> &b){ return a < b[ix]; }
            );

            curve<n> tmp = prg::upper_bound<n>(1, it, last, portion, ix, iy, eps);
            f.insert(f.end(), tmp.begin(), tmp.end());
        }
        assert(prg::is_sorted_by_column<n>(f, 0));

        return f;
    }

    template <unsigned n>
    inline curve<n> closed_bound(const curve<n> &a, unsigned ix, unsigned iy, double eps=0.001)
    {
		curve<n> fl = prg::lower_bound<n>(a, ix, iy, eps);
		curve<n> fu = prg::upper_bound<n>(a, ix, iy, eps);

		for (auto it=fu.rbegin(); it!=fu.rend(); ++it) fl.push_back(*it);
 		return fl;
    }

    template <unsigned n>
    inline curve<n> closed_bound(const curve<n> &a, const curve<n> &b)
    {
		curve<n> f(a);
		for (auto it=b.rbegin(); it!=b.rend(); ++it) f.push_back(*it);
 		return f;
    }


    template <unsigned n>
    inline curve<n> closed_bound(const curve<n> &a, unsigned=0, unsigned iy=1, unsigned ierr=2)
    {
		curve<n> f;

		for (auto it=a.begin(); it!=a.end(); ++it)
            { point<n> p(*it); p[iy] -= p[ierr]; f.push_back(p); }

		for (auto it=a.rbegin(); it!=a.rend(); ++it)
            { point<n> p(*it); p[iy] += p[ierr]; f.push_back(p); }

 		return f;
    }

    template <unsigned n>
    inline curve<n> select(const curve<n> &a, unsigned ix=0, unsigned iy=1)
    {
        point<n> p; for (unsigned i=0; i<n; i++) p[i] = 0;

		curve<n> f;
		for (auto e : a)
		{
            p[0] = e[ix];
            p[1] = e[iy];
            f.push_back(p);
        }
 		return f;
    }

    template <unsigned n>
    inline curve<n> bounds(const curve<n> &a, unsigned ix=0, unsigned iy=1, unsigned ierr=2)
    {
		curve<n> f;

		for (auto e : a)
        {
            point<n> p(e);
            p[0] = e[ix];
            p[1] = e[iy] - e[ierr];
            p[2] = e[iy] + e[ierr];
            f.push_back(p);
        }

 		return f;
    }

    template <unsigned n>
    inline curve<n> to_int(const curve<n> &a, std::vector<unsigned> indices)
    {
		curve<n> f(a);
		for (auto &e : f) { for (auto i : indices) e[i] = std::round(e[i]); }
		return f;
    }

    /*!
     * \brief сосчитать число элементов с одним ключом
     * \param ix позиция ключа
     * \return массив из пар <key, count>
     */
    template <unsigned n>
    inline curve<n> count(const curve<n> &a, unsigned ix=0, bool norm=false)
    {
        assert(prg::is_sorted_by_column<n>(a, ix));

        std::map<double, unsigned> sum; unsigned cnt = 0;
        for (auto e : a)
        {
            if (sum.count(e[ix]) == 0) sum[e[ix]] = 0;
            sum[e[ix]]++; cnt++;
        }

		point<n> p = {0, 0, 0};

		curve<n> tmp;
		for (auto e : sum)
        {
            p[0] = e.first;
            p[1] = norm ? p[1] = double(e.second) / cnt : e.second;
            tmp.push_back(p);
        }
 		return tmp;
    }

}
#endif // PRG_CURVE
