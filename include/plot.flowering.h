#ifndef PLOT_FLOWERING_H
#define PLOT_FLOWERING_H

#include "defs.h"
#include "math.defs.h"
#include "record.flower.h"
#include "flower.dome.h"

//--------------------------------------------------------------------------------------------
//          функции отрисовки данных (общие для учетных и модельных данных
//--------------------------------------------------------------------------------------------

/// отрисовать функцию плотности вероятности (PDF) по числу цветков 1-го и 2-го порядков
void plot_countPDF(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2, const curve_t &vs);

/// отрисовать регрессию отношения цветков 1-го и 2-го порядков
void plot_countRF(const std::string &output, const std::string &species, const curve_t &data1);

/// отрисовать регрессию плотности цветения 1-го и 2-го порядков
void plot_densityRF(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2);

/// отрисовать функцию плотности вероятности (PDF) по длительности фазы цветков 1-го и 2-го порядков
void plot_spanPDF(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2);
void plot_spanPDF(const std::string &output, const std::string &species,
    const std::vector<curve_t> &v1, const std::vector<curve_t> &v2);

/// отрисовать регрессию длительности фазы цветков 1-го и 2-го порядков
void plot_spanRF(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2);
void plot_spanRF(const std::string &output, const std::string &species,
    const std::vector<curve_t> &v1, const std::vector<curve_t> &v2);

/// отрисовать регрессию начала фазы цветков 1-го и 2-го порядков
void plot_axeRF(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2, const std::string &xrange = "[-5:*]");

void plot_phase_range(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2,
    const std::string &xrange = "[-5:*]");
void plot_phase_range(const std::string &output, const std::string &species,
    const std::vector<curve_t> &v1, const std::vector<curve_t> &v2,
    const std::string &xrange = "[-5:*]");

void plot_phase_curve(const std::string &output, const std::string &species,
    const curve_t &v1, const curve_t &v2,
    const std::string &xrange = "[-5:*]");
void plot_phase_curve(const std::string &output, const std::string &species,
    const std::vector<curve_t> &v1, const std::vector<curve_t> &v2,
    const std::string &xrange = "[-5:*]");


#endif // PLOT_FLOWERING_H
