#include <iostream>
#include <map>
#include <limits>

#include "prg.defs.h"
#include "math/prg.distribution.h"

#include "defs.h"
#include "math.defs.h"
#include "math.WLSM.h"
#include "math.gnuplot.h"
#include "data.recovery.h"
#include "record.phase.h"
#include "flower.dome.h"
#include "test.input.file.h"
#include "flower.ensemble.h"
#include "flower.generator.h"
#include "flower.structure.h"
#include "plot.flowering.h"

const char *species_name[] = 
{
   "C.alliariifolia",   // 0
   "C.bononiensis",     // 1
   "C.carpatica",       // 2
   "C.lactiflora",      // 3
   "C.punctata",        // 4
   "C.rapunculoides",   // 5
   "C.sarmatica",       // 6
   "C.trachelium",      // 7
   "P.grandiflorus",    // 8
   "C.test"             // 9
};

std::random_device randevice;
//std::mt19937 rangen(randevice());
std::mt19937 rangen(4);

const unsigned runs = 100; // максимальное число генераций ансамбля
const unsigned used_runs = 10; // число используемых ансамблей с максимальной вероятностью
const unsigned num_iterations = 250; // число итераций MLM

#define PLOT_EXPER 
#define PLOT_MODEL 

#define USE_DAY_ALIGN
#ifdef USE_DAY_ALIGN
    #define DAY_ALIGN DAY_ALIGN_CODE  
#else
    #define DAY_ALIGN 0
#endif

int main(int argc, char * argv[])
{

    if (argc <= 1)
    {
        std::cout << "using: ornamental DATA_DIR OUTPUT_DIR\n";
        exit(0);
    }
    std::string DATA_DIR = argv[1]; 
    std::string OUTPUT_DIR = argv[2];


    SELECTED_RUN(false, "проверка корректности всех файлов данных")
    {
        test_directory(DATA_DIR);
    }

    std::vector<unsigned> years{1999, 2000, 2009, 2010, 2016};
    std::vector<unsigned> species{6};
    std::string xrange("[-5:30]");

    bool skip_death = false; // пропуск умерших цветков
    bool print_flag = false; // отладочная печать

    SELECTED_RUN(true, "демонстрация паттернов цветения учетных и модельных побегов файла")
    {
        for (auto ispecies : species)
        {
            auto filename = make_ifilename(DATA_DIR, ispecies);
            auto speciesname = species_name[ispecies];
            std::vector<unsigned> years{2009, 2010, 2016};
            std::vector<unsigned> sprouts{1, 2, 3};

            //----------------------------------------------------------------------------------
            //          демонстрация плотности распределения цветков учетных побегов файла
            //----------------------------------------------------------------------------------
            Flowerstruct flower_structure(filename);
            for (auto i : years)
            {
                std::string filename = OUTPUT_DIR + std::to_string(i) + "-";

                auto g1 = flower_structure.count_select({i}, 1);
                auto g2 = flower_structure.count_select({i}, 2);
                auto gs = flower_structure.count_select({i});
                auto gr = flower_structure.rate_select({i});
                plot_countPDF(filename + "countPDF[E].pdf", speciesname, g1, g2, gs);
                plot_countRF(filename + "countRF[E].pdf", speciesname, gr);
            }
            auto g1 = flower_structure.count_select(years, 1);
            auto g2 = flower_structure.count_select(years, 2);
            auto gs = flower_structure.count_select(years);
            auto gr = flower_structure.rate_select(years);
            plot_countPDF(OUTPUT_DIR + "countPDF[E].pdf", speciesname, g1, g2, gs);
            plot_countRF(OUTPUT_DIR + "countRF[E].pdf", speciesname, gr);

            Flowerensemble ensemble(used_runs, filename, sprouts, skip_death);
            ensemble.generate(runs, num_iterations, print_flag);

            for (auto i : sprouts)
            {
            #ifdef PLOT_MODEL
                Flowergenerator generator(ensemble, {i});
                auto model = generator.generate_sprout();
            #endif

                for (unsigned k=2; k<=2; k++)
                {
                    phase_enum_t phase = (k==1 ? BUD : FLOWER);
                    std::string filename = OUTPUT_DIR + std::to_string(i) + (k==1 ? "b-" : "f-");

            #ifdef PLOT_EXPER
                    //----------------------------------------------------------------------
                    //      демонстрация паттернов цветения учетных побегов файла
                    //----------------------------------------------------------------------
                    auto f1 = ensemble.phase_select(phase, i, 1, DAY_ALIGN | NORM | ALIGN_SPROUT);
                    auto f2 = ensemble.phase_select(phase, i, 2, DAY_ALIGN | NORM | ALIGN_1AXIS);
                    plot_axeRF(filename + "axe[E].pdf", speciesname, f1, f2, xrange);
                    plot_densityRF(filename + "density[E].pdf", speciesname, f1, f2);
                    plot_spanRF(filename + "spanRF[E].pdf", speciesname, f1, f2);
                    plot_spanPDF(filename + "spanPDF[E].pdf", speciesname, f1, f2);

                    f2 = ensemble.phase_select(phase, i, 2, DAY_ALIGN | NORM | ALIGN_SPROUT);
                    plot_phase_range(filename + "range[E].pdf", speciesname, f1, f2, xrange);
                    plot_phase_curve(filename + "curve[E].pdf", speciesname, f1, f2, xrange);
            #endif
            #ifdef PLOT_MODEL
                    //----------------------------------------------------------------------
                    //      демонстрация паттернов цветения модельных побегов файла
                    //----------------------------------------------------------------------
                    f1 = model->flowering_select(1, DAY_ALIGN | NORM | ALIGN_SPROUT);
                    f2 = model->flowering_select(2, DAY_ALIGN | NORM | ALIGN_1AXIS);
                    plot_axeRF(filename + "axe[M].pdf", speciesname, f1, f2, xrange);
                    plot_densityRF(filename + "density[M].pdf", speciesname, f1, f2);
                    plot_spanRF(filename + "spanRF[M].pdf", speciesname, f1, f2);
                    plot_spanPDF(filename + "spanPDF[M].pdf", speciesname, f1, f2);

                    f2 = model->flowering_select(2, DAY_ALIGN | NORM | ALIGN_SPROUT);
                    plot_phase_range(filename + "range[M].pdf", speciesname, f1, f2, xrange);
                    plot_phase_curve(filename + "curve[M].pdf", speciesname, f1, f2, xrange);
            #endif
                }
            }
        }
    }

    return 0;
}
