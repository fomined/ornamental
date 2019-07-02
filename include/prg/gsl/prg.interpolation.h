#ifndef MATH_INPERPOLATION_H
#define MATH_INPERPOLATION_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>
#include <array>
#include <assert.h>
#include <limits>

#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include "prg.defs.h"
#include "prg.curve.h"

namespace prg
{

    /// найти линейной интерполяцией значение x, где f(x) = 0
    inline double linear_interpolation(double xa, double ya, double xb, double yb)
    {
        assert(std::abs(yb - ya) > std::numeric_limits<double>::min());
        return (yb * xa - ya * xb) / (yb - ya);
    }

    /// найти значение f(x) линейной интерполяцией
    inline double linear_interpolation(double xa, double ya, double xb, double yb, double x)
    {
        assert(std::abs(xb - xa) > std::numeric_limits<double>::min());
        return ya + (yb - ya) * (x - xa) / (xb - xa);
    }


    enum { linear, polynomial, cspline, cspline_periodic, akima, akima_periodic, steffen };

    /*!
     * \brief нахождение интерполяцией точки x, где f(x) = 0
     * \param x[] массив x
     * \param y[] массив f(x)
     * \return значение x, где f(x) = 0
     */
    template <size_t> double find_zero_by_interpolation(size_t n, double x[], double y[]);

    /*!
     * \brief нахождение линйной интерполяцией точки x, где f(x) = 0
     * \param x[] массив x
     * \param y[] массив f(x)
     * \return значение x, где f(x) = 0
     */
    template <> inline double find_zero_by_interpolation<linear>(size_t n, double x[], double y[])
    {
        assert(n >= 2);
        if (n == 2) return (y[1] * x[0] - y[0] * x[1]) / (y[1] - y[0]);
        TODO("требуется реализация метода наименьших квадратов")
    }


    const std::map<size_t, const gsl_interp_type *> interp_method =
    {
        { linear, gsl_interp_linear },
        { polynomial, gsl_interp_polynomial },
        { cspline, gsl_interp_cspline },
        { cspline_periodic, gsl_interp_cspline_periodic },
        { akima, gsl_interp_akima },
        { akima_periodic, gsl_interp_akima_periodic },
        { steffen, gsl_interp_steffen }
    };

    /*!
     * \brief Объект, обеспечивающий интерполяцию функцию между заданными точками
     * \note тип интерполяции задается параметром шаблона, например Interpolator<steffen>
     */
    template <unsigned ID = steffen>
    class Interpolator
    {
        size_t N;
        gsl_interp_accel *acc;
        gsl_spline *spline;
        double xmin, xmax;

    public:

        /// создать интерполятор по заданому набору точек
        Interpolator(size_t n, double *x, double *y);
        Interpolator(const curve<2> &a);

        ~Interpolator();

        /// минимальное число точек для возможности интерполяции
        unsigned min_size() { return gsl_spline_min_size(spline); }

        double min() const { return xmin; }
        double max() const { return xmax; }

        /// интерполировать функцию или ее производные порядка k=[0..2] в заданной точке
        double operator()(double x) const;
        double operator()(unsigned k, double x) const;

        /// интерполировать функцию или ее производные порядка k=[0..2] во всех заданных точках
        template <typename Iterator, typename F>
        curve<2> operator()(unsigned k, Iterator first, Iterator last, F f)
        {
            curve<2> tmp;
            for (; first!=last; ++first)
            {
                double x = f(*first);
                tmp.push_back(point<2>{x, (*this)(k, x)});
            }
            return tmp;
        }

        /// найти минимум/максимум функции f/производной f в заданном интервале
        template <typename F>
        point<2> minimum(unsigned k, double xa, double xb, Unsigned num_ticks, F f) const
        {
            assert(k < 3);
            assert(xa >= xmin);
            assert(xb <= xmax);

            double hx = (xb - xa) / (num_ticks - 1);
            if (xa > xb) std::swap(xa, xb);

            double xmin = xa; double fmin = std::numeric_limits<double>::max();
            double x = xa;
            for (Unsigned i=0; i<num_ticks; i++, x+=hx)
            {
                double y = this->operator()(k, x);
                double fy = f(point<2>{x, y});
                if (fmin > fy) { xmin = x; fmin = fy; }
            }
            return point<2>{xmin, fmin};
        }

        /// найти минимум/максимум функции/производной в заданном интервале
        point<2> minimum(unsigned k, double xa, double xb, unsigned num_ticks);
        point<2> maximum(unsigned k, double xa, double xb, unsigned num_ticks);

    };

    #define THEADER template <unsigned ID>
    #define TARG ID

    THEADER
    inline Interpolator<TARG>::Interpolator(size_t n, double *x, double *y) : N(n)
    {
        acc = gsl_interp_accel_alloc();
        spline = gsl_spline_alloc(interp_method.at(ID), N);
        gsl_spline_init(spline, x, y, N);
        xmin = x[0]; xmax = std::nextafter(x[N-1], x[0]);
            // gsl не определяет поведение в правой точке диапазона (в текущей реализации возврат 0)
            // потому верхную границу чуть снизим, чтобы можно было обращаться к крайней точке
    }

    THEADER
    inline Interpolator<TARG>::Interpolator(const curve<2> &a) : N(a.size())
    {
        acc = gsl_interp_accel_alloc();
        spline = gsl_spline_alloc(interp_method.at(ID), N);

        std::vector<double> x, y; x.reserve(N); y.reserve(N);
        for (auto e : a) { x.push_back(e[0]); y.push_back(e[1]); }

        gsl_spline_init(spline, &x[0], &y[0], N);
        xmin = x[0]; xmax = std::nextafter(x[N-1], x[0]);
    }

    THEADER
    inline Interpolator<TARG>::~Interpolator()
    {
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }

    THEADER
    inline double Interpolator<TARG>::operator()(double x) const
    {
        if (x < xmin) x = xmin;
        if (x > xmax) x = xmax;
        return gsl_spline_eval(spline, x, acc);
    }

    THEADER
    inline double Interpolator<TARG>::operator()(unsigned k, double x) const
    {
        if (x < xmin) x = xmin;
        if (x > xmax) x = xmax;
        switch (k)
        {
        case 0 : return gsl_spline_eval(spline, x, acc);
        case 1 : return gsl_spline_eval_deriv(spline, x, acc);
        case 2 : return gsl_spline_eval_deriv2(spline, x, acc);
        }
        return 0;
    }

    THEADER
    inline point<2> Interpolator<TARG>::minimum(unsigned k, double xa, double xb, unsigned num_ticks)
    {
        return minimum(k, xa, xb, num_ticks, [](const point<2> &a){ return a[1]; });
    }

    THEADER
    inline point<2> Interpolator<TARG>::maximum(unsigned k, double xa, double xb, unsigned num_ticks)
    {
        point<2> p = minimum(k, xa, xb, num_ticks, [](const point<2> &a){ return -a[1]; });
        return point<2>{p[0], -p[1]};
    }

    #undef THEADER
    #undef TARG

    double test_interp1(size_t N, double *x, double *y, double xi);
    double test_interp2(size_t N, double *x, double *y, double xi);


//	/*!
//	 * \brief переинтерполирует заданный набор данных на точки новой кривой
//	 * \note необходимость переинтерполяции обусловлена тем, что все кривые задаются на наборах точек,
//	 *  и эти наборы могут не совпадать
//	 * \param f кривая, в точках которой переинтерполируются данные кривой data
//	 * \param data кривая, точки которой используются для аппроксимации функции в точках f
//	 */
//	inline void approximate_at_points_of(curve<2> &f, const curve<2> &data)
//	{
//		assert(prg::is_sorted_by_column<2>(data, 0));
//		prg::Interpolator<> interpolator(data); // создание интерполятора на точках data
//		interpolator(f); // переинтерполяция точек f
//	}

} // namespace prg
#endif // MATH_INPERPOLATION_H
