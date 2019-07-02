#ifndef MODEL_AXIS2ND_H
#define MODEL_AXIS2ND_H

#include <tuple>
#include "gsl/prg.minimizer.h"

#include "defs.h"
#include "math.defs.h"
#include "record.flower.h"
#include "flower.dome.h"
#include "model.axis1st.h"

curve<2> make_2nd_axis_flowering_curve_pattern(const std::vector<flowerline_t> &flowerlines,
    unsigned indent=5);

flowering_pattern<2> make_2nd_axis_efflorescence_pattern(const std::string &filename,
    const std::vector<unsigned> &sprouts, bool ynorm = true, unsigned align=0);


/*!
 * \brief Построить картину плотности цветения оси 2-го порядка
 * \param filename файл данных
 * \param sprouts массив заданных побегов
 * \param ynorm выполнять или нет нормировку на интервал [0,1]
 */
flowering_pattern<2> make_2nd_axis_density_flowering_pattern(const std::string &filename,
    const std::vector<unsigned> &sprouts, bool ynorm=true);

/*!
 * \brief сгладить картину плотности цветения оси 2-го порядка
 * \param f вектор данных плотности цветения
 * \param sigma величина размытия данных
 * \param hy шаг дискретизации по вертикали
 */
flowering_pattern<2> smooth_2nd_axis_density_flowering_pattern(const flowering_pattern<2> &f,
    double sigma, double hy);

/*!
 * \brief _2nd_axis_density_flowering_match
 * \param a
 * \param b
 * \return
 */
double _2nd_axis_density_flowering_match(const flowering_pattern<2> &a, const flowering_pattern<2> &b);

/*!
 * \brief compare_density_flowering_patterns
 * \param ispecies
 * \param isprouts
 * \param jspecies
 * \param jsprouts
 * \param sigma
 * \param hy
 */
void compare_density_flowering_patterns(unsigned ispecies, const std::vector<unsigned> &isprouts,
    unsigned jspecies, const std::vector<unsigned> &jsprouts, double sigma, double hy);




//curve<2> make_sum_axis_flowering_curve_pattern(const std::vector<flowerline_t> &flowerlines,
//    unsigned indent=5);

//curve<2> make_diff_axis_flowering_curve_pattern(const std::vector<flowerline_t> &flowerlines,
//    unsigned indent=5);



///*!
// * \brief получить паттерн цветения 2-й оси заданного побега
// * \param filename файл данных
// * \param sprouts массив заданных побегов
// * \param ynorm выполнять или нет нормировку на интервал [0,1]
// * \return массив нормированных данных
// */
//flowering_pattern<2> make_2nd_axis_flowering_pattern(const std::string &filename,
//    const std::vector<unsigned> &sprouts, bool ynorm = true);

//inline flowering_pattern<2> make_2nd_axis_flowering_pattern(const std::string &filename,
//    unsigned sprout, bool ynorm = true)
//{
//    return make_2nd_axis_flowering_pattern(filename, std::vector<unsigned>(1, sprout), ynorm);
//}


///*!
// * \brief smooth_as_density
// * \param f
// * \param sigma
// * \param hy
// * \return
// */
////flowering_pattern<2> smooth_as_density(const flowering_pattern<2> &f, double sigma, double hy=0.025);
//flowering_pattern<2> interpolate_as_density(const flowering_pattern<2> &f, double sigma, double hy=0.025);

//flowering_pattern<2> make_2nd_axis_time_flowering_pattern(const std::string &filename,
//    const std::vector<unsigned> &sprouts, bool ynorm);

//double density_flowering_match(const flowering_pattern<2> &a, const flowering_pattern<2> &b);

//void print_density_flowering_pattern_match_matrix(unsigned ntrain, double sigma, double hy=0.025);

//void compare_density_flowering_patterns(unsigned ispecies, unsigned isprout,
//    unsigned jspecies, unsigned jsprout, double sigma, double hy);


//flowering_pattern<2> make_2nd_axis_time_flowering_curve(const std::string &filename,
//    const std::vector<unsigned> &sprouts, bool ynorm);

//inline flowering_pattern<2> make_2nd_axis_time_flowering_curve(const std::string &filename,
//    unsigned sprout, bool ynorm=true)
//{
//    return make_2nd_axis_time_flowering_curve(filename, std::vector<unsigned>(1, sprout), ynorm);
//}


//flowering_pattern<2> smooth_2nd_axis_time_flowering_curve(const flowering_pattern<2> &f,
//    double sigma, double hy);

///*!
// * \brief печатает таблицу совпадений по паттернам цветения главной оси
// * \note модель создается по каждому побегу, и проверяется на других побегах (исключая тот,
// *  который использовался для модели)
// */
//void print_2nd_axis_flowering_pattern_match_matrix(unsigned ntrain, double sigma, double hy);

#endif // MODEL_AXIS2ND_H
