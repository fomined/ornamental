#ifndef PRG_ROOTFINDER_H
#define PRG_ROOTFINDER_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "prg.defs.h"
#include "prg.logger.h"

namespace prg
{
    /**
     * @brief Образец обертки для вычисления минимума квадратичной функции.
     * @note Любая обертка обязана предоставить реализацию функции f.
     */
    struct rootfinder1D_fn_sample
    {
        double f(double x) { return x * x; }
    };

    /**
     * \brief Класс обеспечивает интерфейс к одномерным gsl минимизаторам, позволяя скрыть
     *  от пользователя детали работы с ними и их интерфейсы.
     */
    template <class Wrapper>
    class Rootfinder1D : public Wrapper // копия Wrapper в Minimizer дает безопасность в MT режиме
    {
        typedef Rootfinder1D *this_ptr;
        static double f_(double x, void *p) { return this_ptr(p)->f(x); }

    public:

        enum algorithm_enum_t { bisection = 0, falsepos, brent };

        struct params_t
        {
            algorithm_enum_t algorithm;
            double epsabs;
            double epsrel;
            bool print_flag;

            /*!
             * \note условие выхода |a-b| < epsabs + epsrel*min(|a|,|b|)
             */
            params_t(algorithm_enum_t algorithm=brent, double epsabs=1e-4, double epsrel=0,
                bool print_flag=false) : algorithm(algorithm), epsabs(epsabs), epsrel(epsrel),
                print_flag(print_flag) {}

        } params;

        typedef params_t params_type;

        /// Создание минимизатора для объектов, имеющих реализацию основных функций f(), fdf()
        Rootfinder1D(const Wrapper &w): Wrapper(w) {}

        /// Получить параметры для переустановки при желании
        params_type &get_params() { return params; }

        /**
         * @brief найти координаты объекта, которые минимизируют функцию, заданную оберткой Wrapper
         * @note найденные координаты минимума НЕ сохраняются в объекте
         *
         * @param max_iteration верхний предел итераций минимизации
         * @param xa,xb интервал, в котором гарантированно есть хотя бы один минимум
         * @return значение функции в минимуме и значение вектора координат, обеспечивающего минимум
         */
        flagged_interval<double> run(size_t max_iteration, double xa, double xb);
    };

    template <class Wrapper>
    flagged_interval<double> Rootfinder1D<Wrapper>::run(Unsigned max_iteration, double xa, double xb)
    {
//        char buffer[1024];
//        std::string debug_info;
        prg::Logger logger(params.print_flag);
            // буфер сброса информации о ходе поиска, при неудаче буфер выводится на экран

        gsl_function F;
        F.function = &f_;
        F.params = static_cast<void *>(this);

        const gsl_root_fsolver_type *gsl_T = nullptr;
        switch (params.algorithm)
        {
        case bisection  : gsl_T = gsl_root_fsolver_bisection; break;
        case falsepos   : gsl_T = gsl_root_fsolver_falsepos; break;
        case brent      : gsl_T = gsl_root_fsolver_brent; break;
        }

        gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_T);
        gsl_root_fsolver_set(s, &F, xa, xb);

        Unsigned iteration = 0;

//        std::sprintf(&buffer[0], "%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "root", "err");
//        debug_info.append(buffer);
        logger.insert("%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "root", "err");

        double x = 0;
        int status = GSL_CONTINUE;
        while (status == GSL_CONTINUE && iteration < max_iteration)
        {
            iteration++;
            status = gsl_root_fsolver_iterate(s);

            x = gsl_root_fsolver_root(s);
            xa = gsl_root_fsolver_x_lower(s);
            xb = gsl_root_fsolver_x_upper(s);

            status = gsl_root_test_interval(xa, xb, params.epsabs, params.epsrel);

//            if (status == GSL_SUCCESS)
//            {
//                std::sprintf(&buffer[0], "Converged:\n");
//                debug_info.append(buffer);
//            }
//            std::sprintf(&buffer[0], "%5ld [%.7lf, %.7lf] %.7lf %.7lf\n", iteration, xa, xb, x, xb - xa);
//            debug_info.append(buffer);
            if (status == GSL_SUCCESS) logger.insert("Converged:\n");
            logger.insert("%5ld [%.7lf, %.7lf] %.7lf %.7lf\n", iteration, xa, xb, x, xb - xa);
        }
        gsl_root_fsolver_free(s);

        if (status == GSL_CONTINUE) logger.insert("*** Don't converged ***\n");
        else if (status != GSL_CONTINUE) logger.insert("*** Unknown error ***\n");

        if (status != GSL_SUCCESS) logger.print_flag = true;

//        if (params.print_flag || status != GSL_SUCCESS)
//        {
//            if (status == GSL_CONTINUE)
//            {
//                std::sprintf(&buffer[0], "*** Don't converged ***\n");
//                debug_info.append(buffer);
//            }

//            std::cout << debug_info;
//            std::cout.flush();
//        }

        return flagged_interval<double>({xa, xb}, status == GSL_SUCCESS);
    }

} // namespace prg
#endif // PRG_ROOTFINDER_H
