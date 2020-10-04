#ifndef PRG_CURVE
#define PRG_CURVE

#include <array>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <random>

#include "math/math.h"
#include "math/point.h"

namespace prg
{

    /// базисная структура "график"
    template <unsigned n> using curve = std::vector<point<n>>;

    template <unsigned n, typename Unary>
    inline void foreach(curve<n> &curv, Unary op) { for (auto &e : curv) op(e); }

    template <unsigned n, typename Binary>
    double accumulate(curve<n> &curv, double init, Binary op)
        {  double sum = init; for (auto e : curv) sum = op(sum, e); return sum; }


    /// обеспечить, чтобы во всех точках функция была >= 0
    template <unsigned n>
    inline void positize(curve<n> &curv, Unsigned iy=1)
        { foreach(curv, [iy](auto &a){ a[iy] = a[iy] < 0 ? 0 : a[iy]; }); }

    /// нормализовать кривую ([NORM_MAXIMUM] на максимум, [NORM_POWER] на площадь)
    template <unsigned n>
    inline void normalize(curve<n> &curv, Unsigned ix=0, Unsigned iy=1, Unsigned norm=NORM_POWER)
    {
        double norm_coef;
        if (norm == NORM_POWER)
        {
            double sum = 0;
            double xprev = curv.front()[ix];
            for (auto e : curv)
            {
                double hx = e[ix] - xprev;
                sum += e[iy] * hx;
                xprev = e[ix];
            }
            norm_coef = 1. / sum;
        }
        else if (norm == NORM_MAXIMUN)
        {
            double max = -inf<double>(); for (auto e : curv) if (max < e[iy]) max = e[iy];
            norm_coef = 1. / max;
        }
        for (auto &e : curv) e[iy] *= norm_coef;
    }


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

    inline bool greater(double a, double b, double eps=0.001)
        { return equal(a, b, eps) ? false : a > b; }

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
    inline shared<curve<n>> summarize(shared<curve<n>> a, shared<curve<n>> b, double eps=0.1)
    {
        if (b->empty()) return a;
        if (a->empty()) return b;
        shared<curve<n>> out(new curve<n>(summarize<n>(*a, *b, eps)));
        return out;
    }

    template <unsigned n>
    inline curve<n> subtract(const curve<n> &a, const curve<n> &b, double eps=0.1)
    {
        curve<n> tmp;

		unsigned i=0, j=0;
		unsigned na = a.size(), nb = b.size();
		for (; i<na && j<nb; )
		{
			double xa = a[i][0], xb = b[j][0];
            if (less(xa, xb, eps))
            {
                tmp.push_back(a[i]);
                i++;
                continue;
            }
            if (less(xb, xa, eps))
            {
                point<n> p = negate<n>(b[j]);
                p[0] = xb; // negative all component but not p[0]
                tmp.push_back(p);
                j++;
                continue;
            }

			point<n> p; p[0] = xa;
			for (unsigned k=1; k<n; k++) p[k] = a[i][k] - b[j][k];
			tmp.push_back(p);
			i++; j++;
		}
		if (i != na) for (; i<na; i++) tmp.push_back(a[i]);
        if (j != nb) for (; j<nb; j++)
        {
            double x = b[j][0];
            point<n> p = negate<n>(b[j]);
            p[0] = x;
            tmp.push_back(p);
        }

		return tmp;
    }

    template <unsigned n>
    inline shared<curve<n>> subtract(shared<curve<n>> a, shared<curve<n>> b, double eps=0.001)
    {
        if (b->empty()) return a;
        if (a->empty())
        {
            shared<curve<n>> out(new curve<n>(*b));
            for (auto &e : *out) e[1] = -e[1];
            return out;
        }
        shared<curve<n>> out(new curve<n>(subtract<n>(*a, *b, eps)));
        return out;
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
