#ifndef PRG_MINIMIZER_H
#define PRG_MINIMIZER_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>

#include "test/logger.h"

namespace prg
{
    /**
     * @brief Образец обертки для вычисления минимума квадратичной функции.
     * @note Любая обертка обязана предоставить реализацию функции f.
     */
    struct minimizer1D_fn_sample
    {
        double f(double x) { return x * x; }
    };

    /**
     * \brief Класс обеспечивает интерфейс к одномерным gsl минимизаторам, позволяя скрыть
     *  от пользователя детали работы с ними и их интерфейсы.
     */
    template <class Wrapper>
    class Minimizer1D : public Wrapper // копия Wrapper в Minimizer дает безопасность в MT режиме
    {
        typedef Minimizer1D *this_ptr;
        static double f_(double x, void *p) { return this_ptr(p)->f(x); }

    public:

        enum algorithm_enum_t { goldensection = 0, brent, quad_golden };

        struct params_t
        {
            algorithm_enum_t algorithm = brent;
            double epsabs = 1e-4;
            double epsrel = 0;
            bool print_flag = false;
            std::string filename = "";

            /*!
             * \note условие выхода |a-b| < epsabs + epsrel*min(|a|,|b|)
             */
//            params_t(algorithm_enum_t algorithm=brent, double epsabs=1e-4, double epsrel=0,
//                bool print_flag=false, const std::string &filename="")
//            : algorithm(algorithm), epsabs(epsabs), epsrel(epsrel),
//                print_flag(print_flag), filename(filename) {}

        } params;

        typedef params_t params_type;

        /// Создание минимизатора для объектов, имеющих реализацию основных функций f(), fdf()
        Minimizer1D(const Wrapper &w): Wrapper(w) {}

        /// Получить параметры для переустановки при желании
        params_type &get_params() { return params; }

        /**
         * @brief найти координаты объекта, которые минимизируют функцию, заданную оберткой Wrapper
         * @note найденные координаты минимума НЕ сохраняются в объекте
         *
         * @param max_iteration верхний предел итераций минимизации
         * @param x оценка положения минимума
         * @param xa,xb интервал, в котором гарантированно есть хотя бы один минимум
         * @return значение функции в минимуме и значение вектора координат, обеспечивающего минимум
         */
        flagged_value<double> run(size_t max_iteration, double x, double xa, double xb);

        flagged_value<double> run(size_t max_iteration, double xa, double xb)
            { return run(max_iteration, 0.5 * (xa + xb), xa, xb); }
    };

    template <class Wrapper>
    flagged_value<double> Minimizer1D<Wrapper>::run(size_t max_iteration, double x, double xa, double xb)
    {
        prg::Logger logger(params.print_flag);

        gsl_function F;
        F.function = &f_;
        F.params = static_cast<void *>(this);

        const gsl_min_fminimizer_type *gsl_T = nullptr;
        switch (params.algorithm)
        {
        case goldensection  : gsl_T = gsl_min_fminimizer_goldensection; break;
        case brent          : gsl_T = gsl_min_fminimizer_brent; break;
        case quad_golden    : gsl_T = gsl_min_fminimizer_quad_golden; break;
        }

        gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(gsl_T);
        gsl_min_fminimizer_set(s, &F, x, xa, xb);

        Unsigned iteration = 0;
        logger.insert("%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "min", "err");
        logger.insert("%5ld [%.7lf, %.7lf] %.7lf %.7lf\n", iteration, xa, xb, x, xb - xa);

        int status = GSL_CONTINUE;
        while (status == GSL_CONTINUE && iteration < max_iteration)
        {
            iteration++;
            status = gsl_min_fminimizer_iterate(s);

            x = gsl_min_fminimizer_x_minimum(s);
            xa = gsl_min_fminimizer_x_lower(s);
            xb = gsl_min_fminimizer_x_upper(s);

            status = gsl_min_test_interval(xa, xb, params.epsabs, params.epsrel);

            if (status == GSL_SUCCESS) logger.insert("Converged:\n");
            logger.insert("%5ld [%.7lf, %.7lf] %.7lf %.7lf\n", iteration, xa, xb, x, xb - xa);
        }
        if (status == GSL_CONTINUE) logger.insert("*** Don't converged ***\n");
        else if (status != GSL_CONTINUE) logger.insert("*** Unknown error ***\n");
//        if (status != GSL_SUCCESS) logger.print_flag = true;

        gsl_min_fminimizer_free(s);
        return flagged_value<double>(x, status == GSL_SUCCESS);
    }

    /**
     * @brief Образец обертки для вычисления минимума квадратичной функции.
     * @note Любая обертка обязана предоставить реализацию функций f, fdf.
     */
    struct minimizer_fn_sample
    {
        double f(const double x[]) { return (*x) * (*x); }
        double fdf(double g[], const double x[]) { *g = 2 * (*x); return f(x); }
    };

    /**
     * \brief Класс обеспечивает интерфейс к gsl минимизаторам, позволяя скрыть от пользователя
     *  детали работы с ними и их интерфейсы.
     *  При реализации
     * \note детали работы с объектом скрываются в классе Wrapper
     */
    template <class Wrapper>
    class Minimizer : public Wrapper // копия Wrapper в Minimizer дает безопасность в MT режиме
    {
        typedef Minimizer *this_ptr;
        unsigned n; // размерности массивов
        double *x, *g; // координаты и производные (копии указателей), память захвачена в run()

		static double f_(const gsl_vector *gsl_x, void *p)
		{
			unsigned n = this_ptr(p)->n;
            double *x = this_ptr(p)->x;
			for (unsigned i=0; i<n; i++) x[i] = gsl_vector_get(gsl_x, i);
            return this_ptr(p)->f(x);
        }

        static void fdf_(const gsl_vector *gsl_x, void *p, double *value, gsl_vector *gsl_g)
        {
            unsigned n = this_ptr(p)->n;
            double *x = this_ptr(p)->x;
            double *g = this_ptr(p)->g;

            for (unsigned i=0; i<n; i++) x[i] = gsl_vector_get(gsl_x, i);
            *value = this_ptr(p)->fdf(g, x);
            for (unsigned i=0; i<n; i++) gsl_vector_set(gsl_g, i, g[i]);
        }

		static void df_(const gsl_vector *gsl_x, void *p, gsl_vector *gsl_g)
        {
            double value; // неиспользуемая переменная, нужная только для вызова fdf_
            fdf_(gsl_x, p, &value, gsl_g);
        }

    public:

        enum algorithm_enum_t { conjugate_fr = 0, conjugate_pr, bfgs, bfgs2, steepest_descent };

        struct params_t
        {
            algorithm_enum_t algorithm;
            double ftol;
            double gtol;
            int print_flag;
            const char *format;
            unsigned print_limit_left;
            unsigned print_limit_right;

            params_t(algorithm_enum_t algorithm=bfgs, double ftol=1e-4, double gtol=1e-3,
                int print_flag=0, const char *format=" %7.4lf",
                unsigned print_limit_left=3u, unsigned print_limit_right=3u)
            : algorithm(algorithm), ftol(ftol), gtol(gtol), print_flag(print_flag), format(format),
                print_limit_left(print_limit_left), print_limit_right(print_limit_right) {}
        } params;

        typedef params_t params_type;

        /// Создание минимизатора для объектов, имеющих реализацию основных функций f(), fdf()
        Minimizer(const Wrapper &w) : Wrapper(w) {}

        /// Получить параметры для переустановки при желании
        params_type &get_params() { return params; }

        /**
         * @brief найти координаты объекта, которые минимизируют функцию, заданную оберткой Wrapper
         * @note найденные координаты минимума НЕ сохраняются в объекте
         *
         * @param max_iteration верхний предел итераций минимизации
         * @param first стартовый итератор, определяющий начальные координаты объекта
         * @param last конечный итератор, определяющий начальные координаты объекта
         * @return значение функции в минимуме и значение вектора координат, обеспечивающего минимум
         */
        template <typename Iterator>
        flagged_value<double> run(size_t max_iteration, Iterator first, Iterator last);
    };

    template <class Wrapper>
    template <typename Iterator>
    flagged_value<double> Minimizer<Wrapper>::run(size_t max_iteration,
        Iterator first, Iterator last)
    {
        prg::Logger logger(params.print_flag);
        Iterator start = first; // сохраняем, чтобы записать результат по данному адресу

        n = std::distance(first, last);
        if (n == 0) return flagged_value<double>(0., false);

        std::vector<double> tmpx(n, 0.);
        std::vector<double> tmpg(n, 0.);
        x = tmpx.data();
        g = tmpg.data();

        for (size_t i=0; i<n; ++first, i++) x[i] = *first;

        gsl_vector *gsl_x = gsl_vector_alloc(n);
        for (unsigned i=0; i<n; i++) gsl_vector_set(gsl_x, i, x[i]);

        gsl_multimin_function_fdf my_func;
        my_func.n = n;
        my_func.f = f_;
        my_func.df = df_;
        my_func.fdf = fdf_;
        my_func.params = static_cast<void *>(this);

        const gsl_multimin_fdfminimizer_type *gsl_T = nullptr;
        switch (params.algorithm)
        {
        case conjugate_fr   : gsl_T = gsl_multimin_fdfminimizer_conjugate_fr; break;
        case conjugate_pr   : gsl_T = gsl_multimin_fdfminimizer_conjugate_pr; break;
        case bfgs           : gsl_T = gsl_multimin_fdfminimizer_vector_bfgs; break;
        case bfgs2          : gsl_T = gsl_multimin_fdfminimizer_vector_bfgs2; break;
        case steepest_descent : gsl_T = gsl_multimin_fdfminimizer_steepest_descent; break;
        }

        gsl_multimin_fdfminimizer *gsl_s = gsl_multimin_fdfminimizer_alloc(gsl_T, n);
        gsl_multimin_fdfminimizer_set(gsl_s, &my_func, gsl_x, 0.01, params.ftol);

        unsigned print_limit_left = params.print_limit_left;
        unsigned print_limit_right = params.print_limit_right;
        if (print_limit_left + print_limit_right >= n)
            { print_limit_left = n; print_limit_right = 0; }

        int status = GSL_CONTINUE; unsigned iter = 0;
        while (status == GSL_CONTINUE && iter < max_iteration)
        {
            status = gsl_multimin_fdfminimizer_iterate (gsl_s);
            status = gsl_multimin_test_gradient (gsl_s->gradient, params.gtol);

            if (status == GSL_SUCCESS) logger.insert("Converged:\n");
            double g2 = 0, dx2 = 0;
            for (size_t i=0; i<n; i++)
            {
                g2 += sqr(gsl_vector_get(gsl_s->gradient, i));
                dx2 += sqr(gsl_vector_get(gsl_s->dx, i));
            }

            logger.insert("Iteration %4d f=%.9e |g|=%.6e |dx|=%.6e [ ", iter, gsl_s->f, sqrt(g2), sqrt(dx2));

            for (size_t i=0; i<print_limit_left; i++)
                logger.insert(params.format, gsl_vector_get(gsl_s->x, i));

            if (print_limit_left + print_limit_right < n) logger.insert(" ... ");

            for (size_t i=0; i<print_limit_right; i++)
                logger.insert(params.format, gsl_vector_get(gsl_s->x, n-print_limit_right+i));

            logger.insert("]\n");

            iter++;
        }
        if (status == GSL_CONTINUE) logger.insert("*** Don't converged ***\n");
        else if (status != GSL_CONTINUE) logger.insert("*** Unknown error ***\n");
//        if (status != GSL_SUCCESS) logger.print_flag = true;

        if (params.print_flag < 0) logger.print_flag = false;
			// строгий запрет на печать

        double fmin = gsl_s->f;
        for (size_t i=0; i<n; ++start, i++) *start = gsl_vector_get(gsl_s->x, i);

        gsl_multimin_fdfminimizer_free(gsl_s);
        gsl_vector_free(gsl_x);

        return flagged_value<double>(fmin, status == GSL_SUCCESS);
    }

} // namespace prg
#endif // PRG_MINIMIZER_H

