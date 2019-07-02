#ifndef FLOWER_ENSEMBLE_H
#define FLOWER_ENSEMBLE_H

#include <string>
#include <vector>
#include <string>
#include <set>
#include <memory>
#include <tuple>
#include <random>

#include "prg.maxmin.h"
#include "defs.h"
#include "math.defs.h"
#include "math.statistics.h"
#include "record.time.h"
#include "record.phase.h"
#include "record.flower.h"
#include "flower.dome.h"

const unsigned char NORM = 0x1; // нормировка на [0, 1] по высоте
const unsigned char ALIGN_SPROUT = 0x2; // выравнивание старта цветения по побегу
const unsigned char ALIGN_1AXIS = 0x4; // выравнивание старта цветения по 1й оси
const unsigned char DAY_ALIGN_CODE = 0x8; // выравнивание по границе дня

/// "скрытые" параметры являются моделью, построенной по учетным побегам
struct hidden_parameters_type
{
    typedef std::normal_distribution<>::param_type param_type;

    curve_t axis[3]; // начало относительного зацветания 1й и 2й оси vs высоты {null, axis1, axis2}
        // для 1й оси - относительно старта побега, для 2й оси - относительно 1й оси
    point<2> budspan; // продолжительность бутонизации
    curve_t flowerspan[3]; // продолжительность цветения 1й и 2й оси vs высоты {null, axis1, axis2}
    curve_t density; // плотность цветения 2й оси vs высоты
    curve_t between; // среднее расстояние между цветением цветков 2й оси vs высоты
    point<2> num_flowers[3]; // число цветков (с погрешностью) на осях {null, axis1, axis2}
    bool ready_flag = false;

    /// контроль готовности
    bool empty() const { return !ready_flag; }

    flagged_value<param_type> flowernum_params(unsigned iaxis) const
    { return {param_type{num_flowers[iaxis][0], num_flowers[iaxis][1]}, true}; }

    /// старт зацветения заданой оси на заданной высоте
    flagged_value<param_type> axis_params(double h, unsigned iaxis) const
    {
        auto it = std::lower_bound(axis[iaxis].begin(), axis[iaxis].end(), h,
            [](const point<CDIM> &a, double h) { return a[0] < h; }
        );
        if (it == axis[iaxis].end()) --it;
        return {param_type{(*it)[1], (*it)[2]}, true};
    }

    /// длительность цветения заданой оси на заданной высоте
    flagged_value<param_type> budspan_params() const
    { return {param_type{budspan[0], budspan[1]}, true}; }

    /// длительность цветения заданой оси на заданной высоте
    flagged_value<param_type> flowerspan_params(double h, unsigned iaxis) const
    {
        auto it = std::lower_bound(flowerspan[iaxis].begin(), flowerspan[iaxis].end(), h,
            [](const point<CDIM> &a, double h) { return a[0] < h; }
        );
        if (it == flowerspan[iaxis].end()) --it;
        return {param_type{(*it)[1], (*it)[2]}, true};
    }


    /// плотность цветения (2й оси) на заданной высоте
    flagged_value<param_type> density_params(double h, unsigned=0) const
    {
        const double eps = 0.025; // небольшое расширение диапазона для узких по h паракладиев
        auto it = std::lower_bound(density.begin(), density.end(), h,
            [](const point<CDIM> &a, double h) { return a[0] < h; }
        );
        if (it == density.end() && std::abs((*(it-1))[0] - h) < eps ) --it;

        bool flag = true;
        if (it == density.begin() && h + eps < (*it)[0]) flag = false; // при малых h нет 2й оси
        if (it == density.end()) flag = false; // при h > 1 нет 2й оси

        return {param_type{(*it)[1], (*it)[2]}, flag};
    }

    /// плотность цветения (2й оси) на заданной высоте
    flagged_value<param_type> between_params(double h, unsigned=0) const
    {
        auto it = std::lower_bound(between.begin(), between.end(), h,
            [](const point<CDIM> &a, double h) { return a[0] < h; }
        );
        if (it == between.end()) --it;

        bool flag = true;
        if (it == between.begin() && h < (*it)[0]) flag = false;

        return {param_type{(*it)[1], (*it)[2]}, flag};
    }
};



/*!
 * \brief генератор случайных дат внутри заданных диапазонов
 * \note даты генерятся однородным случайным образом, так что длины фаз не согласуются со средними
 */
struct intraranges_uniform_generator
{
    void init(double *, double *) {}

    std::array<double, 3> operator()(const interval<double> *r)
    {
        std::array<double, 3> tmp;
        std::uniform_real_distribution<> d1(r[0][0], r[0][1]);
        tmp[0] = d1(rangen);

        double a = std::max(tmp[0] + 1, r[1][0]);
        std::uniform_real_distribution<> d2(a, r[1][1]);
        tmp[1] = d2(rangen);

        a = std::max(tmp[1] + 1, r[2][0]);
        std::uniform_real_distribution<> d3(a, r[2][1]);
        tmp[2] = d3(rangen);

        return tmp;
    }
};

/*!
 * \brief генератор случайных дат с учетом средних значений внутри заданных диапазонов
 * \note даты генерятся случайным образом в диапазоне BUD-FLOWER и далее генерятся согласно
 * нормальным распределениям длины фаз BUD и FLOWER как отступы от этой даты
 */
struct intraranges_dispersion_generator
{
    typedef std::normal_distribution<> distribution_type;
    typedef std::shared_ptr<distribution_type> pointer_type;

    pointer_type norm_distribution[2]; // генераторы нормальных распределений

    void init(double *mean, double *sigma2)
    {
        norm_distribution[0] = pointer_type(new distribution_type(mean[0], sqrt(sigma2[0])));
        norm_distribution[1] = pointer_type(new distribution_type(mean[1], sqrt(sigma2[1])));
    }

    std::array<double, 3> operator()(const interval<double> *r)
    {
        // генерим значение в интервале смены BUD-FLOWER и от него в обе стороны, используя
        // нормальные распределения длин фаз, пытаемся попасть в интервалы *-BUD и FLOWER-OVARY
        std::array<double, 3> t;

        interval<double> ma = r[0]; // интервал *-BUD
        interval<double> mb = r[1]; // интервал BUD-FLOWER
        interval<double> mc = r[2]; // интервал FLOWER-OVARY

        std::uniform_real_distribution<> uni_distribution(mb[0], mb[1]);
        unsigned num_iteration = 300; // для избегания зацикливания

        bool flag = false; // флаг, что генерация прошла успешно
        while (!flag)
        {
            t[1] = uni_distribution(rangen);

            double mean = norm_distribution[0]->mean();
            double p = t[1] - mean;
            if (p < ma[0]) t[0] = ma[0];
            else if (p > ma[1]) t[0] = ma[1];
            else t[0] = t[1] - (*norm_distribution[0])(rangen);

            mean = norm_distribution[1]->mean();
            p = t[1] + mean;
            if (p < mc[0]) t[2] = mc[0];
            else if (p > mc[1]) t[2] = mc[1];
            else t[2] = t[1] + (*norm_distribution[1])(rangen);


            flag = (ma[0] <= t[0] && t[0] <= ma[1]) && (mc[0] <= t[2] && t[2] <= mc[1]);
            num_iteration--;

        //#define DEBUG_PRINT
        #ifdef DEBUG_PRINT
            if (flag) { printf("OK\n"); FLUSH; }
            if (num_iteration == 0)
            {
                printf("{%lf %lf %lf} {%lf %lf %lf} {%lf %lf %lf}\n",
                   ma[0], t[0], ma[1], mb[0], t[1], mb[1], mc[0], t[2], mc[1]);
                FLUSH;
            }
        #endif

            assert(num_iteration > 0);
        };
        return t;
    }
};

/// ансамбль из разных вариантов восстановления пропущенных дат цветения
class Flowerensemble
{
    Flowerdome dome; // база данных вида
    std::vector<flowerline3_t> bounds; // номер цветка и границы фаз (BUD, FLOWER)
    std::vector<interval<double>> intervals; // (служебное поле для MLM) то же, но в форме для MLM
    Timestamp tstart, tfinish; // предельные оценочные даты начала и конца цветения
        // для цветения предельные даты близки к датам, полученным MLM
        // реальный старт всегда позже, а реальный финиш всегда раньше предельных дат

    typedef std::vector<point<3>> exemple_type;
        // Множество возможных дат цветения для заданного вида.
        // Объекты множества упорядочены по номерам цветков (номера цветков в bounds)

    typedef std::pair<double, std::shared_ptr<exemple_type>> value_type;
        // тип хранения элемента ансамбля (дисперсия ансамблья, ссылка на ансамбль)

    prg::concentrator<value_type, void> concentrator;
        // накопитель наиболее "лучших" представителей ансамбля

protected:

    /*!
     * \brief генерация случайного ансамбля с датами цветения его цветков
     * \param num_iterations число итераций подгона
     * \return логарифм вероятности ансамбля и вектор дат цветения его цветков
     */
    std::pair<double, std::shared_ptr<species_ensemble>> generate_ensemble(
        Unsigned num_iterations, bool print_flag=false) const;

public:

    typedef prg::concentrator<value_type, void> concentrator_type;

    /*!
     * \brief загрузить базу данных из файла и извлечь границы фаз
     * \param n число объектов в ансамбле
     * \param sprouts список загружаемых побегов
     * \param skip_death удаление цветков, если DEATH произошло ранее FLOWER
     * \note чем больше загружено побегов, тем точнее определяются длины фаз, но хуже разбираться
     *  с выравниванием побегов друг относительно друга
     */
    Flowerensemble(Unsigned n, const std::string &filename,
        const std::vector<unsigned> &sprouts={}, bool skip_death=false)
    : dome(filename, sprouts), concentrator(value_type(0, nullptr), n)
    {
        // извлечь интервалы границ фаз (BUD, FLOWER) для всех цветков
        interval<Flowerstamp> fs(0, 255);
            // максимальный номер побега 255, так как используется char для sprout

        bounds = dome.extract_phase_bounds(fs, skip_death, false);

        //----------------------------------------------------------------------------------------
        //      конвертируем данные к формату intrainterval_localization_by_MLM()
        //----------------------------------------------------------------------------------------
        tstart = Timestamp(1000);
        tfinish = Timestamp(0);

        Unsigned ne = bounds.size();
        for (Unsigned i=0; i<ne; i++)
        {
            interval<Timestamp> r0 = bounds[i].intervals[0];
            interval<Timestamp> r1 = bounds[i].intervals[1];
            interval<Timestamp> r2 = bounds[i].intervals[2];

            intervals.push_back(interval<double>{r0[0].to_double(), r0[1].to_double()});
            intervals.push_back(interval<double>{r1[0].to_double(), r1[1].to_double()});
            intervals.push_back(interval<double>{r2[0].to_double(), r2[1].to_double()});

            if (tstart > r1[0]) tstart = r1[0]; // первая дата в переходе BUD-FLOWER
            if (tfinish < r2[1]) tfinish = r2[1]; // последняя дата в переходе FLOWER-OVARY
        }
        assert(tstart.valid());
        assert(tfinish.valid());
    }

    /*!
     * \brief загрузить базу данных из файла и извлечь границы фаз
     * \param runs число запусков генерации
     * \param num_iterations число итераций подгонки данных к более оптимальному ансамблю
     *  (чем больше итераций, тем более узким становится разброс данных в ансамбле,
     *  но тем медленней расчет)
     */
    void generate(Unsigned runs, Unsigned num_iterations, bool print_flag=false);

    /// оператор доступа к записям базы
    concentrator_type::iterator begin() { return concentrator.begin(); }
    concentrator_type::iterator end() { return concentrator.end(); }
    Unsigned size() { return concentrator.size(); }

    /*!
     * \brief получить реальные даты цветения в ансамбле (поиск начала по 1 оси, конца по всем)
     * \note реальная дата может на 1-2 дня быть раньше (зацветание цветка 2-го порядка),
     *  что требует при отрисовке брать отступ в 5 дней
     */
    interval<Timestamp> flowering_dates() const;

    /// определение старта цветения побега (по ансамблю) <mean, stddev>
    std::pair<double, double> flowering_start(unsigned sprout) const;

    /// определение старта цветения побега на заданной высоте (по ансамблю) <mean, stddev>
    std::pair<double, double> flowering_start(unsigned sprout, unsigned height) const;

    /// "высота" побега для нормировки
    unsigned flowering_height(unsigned sprout) const;

    /// получить функцию распределения по числу цветков на оси
    curve_t count_select(const std::vector<unsigned> &sprouts, unsigned axis) const;

    /*!
     * \brief получить усредненные по ансамблю точки стартов и концов цветения для заданных побега и оси
     * \note (y, <ts>, <dt>, sigma-ts, sigma-dt) !!! порядок записи параметров точек кривой
     * \param sprout номер побега
     * \param axis порядок оси
     * \param pos номер цветка на оси заданного порядка (0 - все номера)
     * \return множество записей <абсолютная позиция по вертикали, средний день цветения по ансамблю,
     *  средняя продолжительность цветения,
     *  среднеквадратичные отклонения дня начала и продолжительности цветения>
     */
    curve_t phase_select(phase_enum_t phase, unsigned sprout, unsigned axis,
        unsigned char flags, unsigned pos=0) const;
    curve_t phase_select(phase_enum_t phase, const std::vector<unsigned> &sprouts,
        unsigned axis, unsigned char flags, unsigned pos=0) const;

    curve_t phase_select(Unsigned sample, phase_enum_t phase, unsigned sprout, unsigned axis,
        unsigned char flags, unsigned pos=0) const;
    curve_t phase_select(Unsigned sample, phase_enum_t phase, const std::vector<unsigned> &sprouts,
        unsigned axis, unsigned char flags, unsigned pos=0) const;

    /// полный набор параметров для модели цветения (использовать только после запуска generate)
    hidden_parameters_type get_hidden_parameters(const std::vector<unsigned> &sprouts) const;
};

/// построение без сглаживания гистограммы длительности цветения (с нормировкой на полное число цветков)
curve_t construct_flowerspan(const curve_t &v1, unsigned num_flowers);

/// (отладка) визуализация границ интервала
void print_flowering_bounds(const std::vector<flowerline3_t> &bounds, Unsigned i);

/// (отладка) визуализация границ фазы
void print_flowering_phase(const std::pair<double, std::shared_ptr<species_ensemble>> &ens,
    Unsigned i);

#endif // FLOWER_ENSEMBLE_H
