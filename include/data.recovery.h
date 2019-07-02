#ifndef DATA_RECOVERY_H
#define DATA_RECOVERY_H

#include <tuple>
#include "gsl/prg.minimizer.h"

#include "defs.h"
#include "math.defs.h"
#include "record.flower.h"

/*!
 * \brief восстановление пропущенных данных методом максимального правдоподобия (минимизацией)
 * \param recovered_intervals восстановленные интервалы
 * \param left_bounds левые границы интервалов
 * \param right_intervals правые границы интервалов
 * \param num_iterations число итераций оптимизации при поиске
 * \param print_flag печать информации о процессе восстановления
 * \return (среднее значение, сигма) для восстанавливаемой фазы
 * \note значение регуляционого коэффициента подобрано вручную и "вшито" в код
 */
std::pair<double,double> intrainterval_localization_by_MLM
(
    std::vector<interval<double>> &recovered_intervals,
    const std::vector<interval<double>> &left_bounds,
    const std::vector<interval<double>> &right_bounds,
    size_t num_iterations, bool print_flag
);

/*!
 * \brief восстановление пропущенных данных методом максимального правдоподобия (точным решением)
 * \param recovered_intervals восстановленные интервалы
 * \param left_bounds левые границы интервалов
 * \param right_intervals правые границы интервалов
 * \return (среднее значение, сигма) для восстанавливаемой фазы
 */
std::pair<double,double> intrainterval_localization_by_directMLM(
    std::vector<interval<double>> &restored_values,
    const std::vector<interval<double>> &lb,
    const std::vector<interval<double>> &rb);

/*!
 * \brief восстановление "точных" значений границ фаз методом максимального правдоподобия
 * \param ne число экземляров во множестве экземпляров
 * \param ni число интервалов в каждом экземпляре
 * \param ps указатель на массив всех интервалов
 * \param num_iterations число итераций поиска
 * \param print_flag флаг печати результатов
 * \return
 */
std::vector<double> intrainterval_localization_by_MLM(Unsigned ne, Unsigned ni,
    const interval<double> *ps, Unsigned num_iterations, bool print_flag);

template <class Generator>
/*!
 * \brief восстановление "точных" значений границ фаз методом максимального правдоподобия
 * \param ne число экземляров во множестве экземпляров
 * \param ni число интервалов в каждом экземпляре
 * \param ps указатель на массив всех интервалов
 * \param num_iterations число итераций поиска
 * \param print_flag флаг печати результатов
 * \param generator объект для генерации дат внутри заданных интервалов
 * \return
 */
std::vector<double> intrainterval_localization_by_MLM(Unsigned ne, Unsigned ni,
    const interval<double> *ps, Unsigned num_iterations, Generator generator)
{
    Unsigned np = ne * ni; // полное число оптимизирумых параметров

    //----------------------------------------------------------------------------------------------
    //              задание функций для метода максимального правдоподобия
    //----------------------------------------------------------------------------------------------
    struct MLM
    {
        Unsigned ne; // число экземпляров
        const interval<double> *ps; // указатель на 1-й элемент множества экземпляров
        Unsigned num_phase; // число фаз (едино для всех экземпляров) = ni -1
        Unsigned ni; //
        Unsigned np; // общее число переменных оптимизации = ne * ni
        double reg_coef; // коэффициент регуляризации

        std::vector<double> mean, sigma2; // [num_phase]
            // значения mu & sigma не оптимизируются, а вычисляются для каждой итерации

        /// вычисление средних значений длин всех фаз и дисперсий = sqr(sigma)
        void averages_and_dispersions(const double p[])
        {
            for (Unsigned j=0; j<num_phase; j++)
            {
                mean[j] = 0;
                for (Unsigned i=j; i<np; i+=ni) { mean[j] += p[i+1] - p[i]; }
                mean[j] /= ne;

                sigma2[j] = 0;
                for (Unsigned i=j; i<np; i+=ni) sigma2[j] += sqr(p[i+1] - p[i] - mean[j]);
                sigma2[j] /= ne;
            }
        }

        double f(const double p[])
        {
            averages_and_dispersions(p);

            double L = 0;
            for (Unsigned j=0; j<num_phase; j++) L += 1 + log(sigma2[j]);
            L *= 0.5 * ne;

            double R = 0;
            for (Unsigned i=0; i<np; i++)
            {
                double r0 = p[i] < ps[i][0] ? p[i] - ps[i][0] : 0;
                double r1 = p[i] > ps[i][1] ? p[i] - ps[i][1] : 0;
                R += sqr(r0) + sqr(r1);
            }
            R *= reg_coef;

            return L + 0.5 * R;
        }

        double fdf(double g[], const double p[])
        {
            averages_and_dispersions(p);

            double L = 0;
            for (Unsigned j=0; j<num_phase; j++) L += 1 + log(sigma2[j]);
            L *= 0.5 * ne;

            double R = 0;
            for (Unsigned i=0; i<np; i++)
            {
                double r0 = p[i] < ps[i][0] ? p[i] - ps[i][0] : 0;
                double r1 = p[i] > ps[i][1] ? p[i] - ps[i][1] : 0;
                R += sqr(r0) + sqr(r1);
                g[i] = reg_coef * (r0 + r1);
            }
            R *= reg_coef;

            double gprev = 0;
            for (Unsigned i=0; i<np; i++)
            {
                Unsigned j = i % ni;
                double gcur = (i+1) % ni == 0 ? 0 : (p[i+1] - p[i] - mean[j]) / sigma2[j];
                g[i] += gprev - gcur;
                gprev = gcur;
            }

            return L + 0.5 * R;
        }
    } mlm;

    mlm.ne = ne;
    mlm.ni = ni;
    mlm.ps = ps;
    mlm.num_phase = ni - 1;
    mlm.mean.resize(ni - 1, 0);
    mlm.sigma2.resize(ni - 1, 0);
    mlm.reg_coef = 100000;
    mlm.np = np;

    std::vector<double> x(np, 0.);

    //----------------------------------------------------------------------------------------------
    //              получение оценочных средних значений и дисперсий
    //----------------------------------------------------------------------------------------------
    for (size_t i=0; i<ne; i++)
    {
        interval<double> ra = ps[3*i    ];
        interval<double> rb = ps[3*i + 1];
        interval<double> rc = ps[3*i + 2];
        x[3*i  ] = 0.5 * (ra[0] + ra[1]);
        x[3*i+1] = 0.5 * (rb[0] + rb[1]);
        x[3*i+2] = 0.5 * (rc[0] + rc[1]);
    }
    mlm.averages_and_dispersions(&x[0]);
    double *mean = &mlm.mean[0];
    double *sigma2 = &mlm.sigma2[0];

    //----------------------------------------------------------------------------------------------
    //              вычисление стартовых значений оптимизируемых параметров
    //----------------------------------------------------------------------------------------------
    generator.init(mean, sigma2); // иницилизация генератора средними и дисперсиями
        // для intraranges_uniform_generator это пустая операция
        // для intraranges_dispersion_generator это операция позволяет более точно получать
        // стартовые значения, чтобы уменьшить число итераций для MLM

    for (size_t i=0; i<ne; i++)
    {
        std::array<double, 3> r = generator(&ps[3*i]);
        x[3*i  ] = r[0];
        x[3*i+1] = r[1];
        x[3*i+2] = r[2];
    }

    //---------------------------------------------------------------------------------------------
    //                              оптимизация параметров
    //---------------------------------------------------------------------------------------------
    prg::Minimizer<MLM> minimizer(mlm);
    auto &param = minimizer.get_params();
    param.algorithm = prg::Minimizer<MLM>::bfgs;
    param.print_flag = -1; // строгий запрет на печать
    minimizer.run(num_iterations, &x[0], &x[np]);

    return x;
}

/*!
 * \brief восстановление "точных" значений границ фаз методом улучшенного максимального правдоподобия
 * \param ne число экземляров во множестве экземпляров
 * \param ni число интервалов в каждом экземпляре
 * \param ps указатель на массив всех интервалов
 * \param num_iterations число итераций поиска
 * \param print_flag флаг печати результатов
 * \param generator объект для генерации дат внутри заданных интервалов
 * \return
 */
std::pair<double, std::vector<double>> intrainterval_localization_by_XMLM(Unsigned ne, Unsigned ni,
    const interval<double> *ps, Unsigned num_iterations=20);

#endif // DATA_RECOVERY_H
