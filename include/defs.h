#ifndef DEFS_H
#define DEFS_H

#include <cstddef>
#include <numeric>
#include <array>
#include <type_traits>

#include "prg.defs.h"
#include "prg.curve.h"

using prg::point;
using prg::curve;

typedef std::ptrdiff_t long_t;
typedef std::ptrdiff_t diff_t;
typedef std::size_t size_t;

#define CDIM 6 // размерность обобщенной кривой отрисовки

// содержит следующие данные {h, ts, dt, sigma(ts), sigma(dt), sprout}
// необходимость номера побега нужна, чтобы при получении числа цветков на 2 оси не перепутались
// данные разных побегов
typedef curve<CDIM> curve_t;
typedef curve<CDIM>::value_type point_t;

extern const char *species[];

//const bool PRINT_YES = true;
//const bool PRINT_NO = false;

/// максимальное число побегов в одном файле
const unsigned NUM_SPROUT = 3;
const unsigned NUM_SPECIES = 9;
const unsigned MAX_FLOWERING_TIME = 90;

extern const char *species_name[];
extern const std::string PROJECT_DIR;
extern const std::string DATA_DIR;
extern const std::string FLOWER_DIR;
extern const std::string OUTPUT_DIR;
extern const std::string TMP_DIR;
extern const std::string RAWDATA_DIR;
extern const std::string MAP_FILE;

extern std::random_device randevice;
extern std::mt19937 rangen;

/*!
 * \file файл стандартных определений используемых типов проекта и функции отладки
 */

const size_t BUF_LEN = 2048; // длина записи в файле
const double TT_RATIO = 0.8; // стандартное отношение числа обучающих/тестирующих примеров
const int ITERATIONS = 150; // стандартное число итераций оптимизации
const double THICKNESS = 0.1; // доля объектов, принадлежащих границе

/// получить имя файла для заданного вида
inline std::string make_ifilename(const std::string &dirname, unsigned species)
{ return dirname + species_name[species]; }

/// получить имя файла для заданного вида
inline std::string make_ifilename(const std::string &dirname, const std::string &species)
{ return dirname + species; }

/// получить имя файла для заданного вида
inline std::string make_ofilename(const std::string &dirname, unsigned species,
    const std::vector<unsigned> &sprouts)
{
    std::string postfix("{" + std::to_string(sprouts[0]));
    for (unsigned i=1; i<sprouts.size(); i++)
        postfix += "," + std::to_string(sprouts[i]);
    postfix += "}";
    return dirname + std::string(species_name[species]) + postfix + ".pdf";
}

/// получить имя файла для заданного вида
inline std::string make_ofilename(const std::string &dirname, const std::string &species,
    const std::vector<unsigned> &sprouts)
{
    std::string postfix("{" + std::to_string(sprouts[0]));
    for (unsigned i=1; i<sprouts.size(); i++)
        postfix += "," + std::to_string(sprouts[i]);
    postfix += "}";
    return dirname + species + postfix + ".pdf";
}

/// конвертирование enum в целые числа
template <typename E> inline constexpr size_t to_underlying(E e) noexcept
{ return static_cast<typename std::underlying_type<E>::type>(e); }

inline std::ostream &operator <<(std::ostream & os, const point<2> &s)
{ return os << "{" << s[0] << "," << s[1] << "}"; }

/// литеральная строка со всеми пустыми символами и кавычками, удаляемая trim
const char * const remove_symbols = " \"\t\n\r\f\v";

/// выход из программы при отсутствии необходимого файла
inline void no_file_exit [[noreturn]] (const std::string &filename)
{
    printf("[ERR] no file \'%s\'\n", filename.c_str());
    throw std::exception();
}

/// перевод доли [0,1] в проценты
inline double to_percent(double v) { return v * 100; }

/// выход из программы при отсутствии необходимого файла
#define OPEN_FILE(file, filename) \
    std::ifstream file(filename); \
    if (!file) no_file_exit(filename); \

#define THIS_FILE std::string(std::string(__FILE__), std::string(__FILE__).find_last_of('/') + 1)
#define THIS_LINE std::to_string(__LINE__)
#define FILE_LINE std::string(THIS_FILE + std::string("(") + THIS_LINE + std::string(")"))
#define EXIT { printf("Error at %s\n", FILE_LINE.c_str()); exit(1); }

#define RUN(state, msg) if (state) printf("RUN : %s\n", msg); if (state)
#define SELECTED_RUN(state, msg) if (state) printf("RUN : %s\n", msg); if (state)

template <typename T, typename U>
inline std::ostream &operator <<(std::ostream & os, const std::pair<T, U> &s)
{ return os << "{" << s.first << "," << s.second << "}"; }

template <typename T>
void DEBUG_PRINT(T t) { std::cout << t << "\n"; std::cout.flush(); }

template <typename T, typename U>
void DEBUG_PRINT(T t, U u) { std::cout << t << " " << u << "\n"; std::cout.flush(); }

template <typename T>
void DEBUG_PRINT(const std::vector<T> &t)
{ for (auto e:t) std::cout << e << "\n"; std::cout.flush(); }

template <typename T> void _DEBUG_PRINT(const T &) {}

/// предупреждение о неиспользовании незавершенного или ошибочного кода
#define TODO(msg) { std::cout << msg << "\n"; assert(false); }

template <unsigned n>
inline curve<n> lower_bound(const curve<n> &f, double portion=0.5, unsigned ix=0, unsigned iy=1)
{
    if (f.empty()) return curve<n>();
    curve<n> tmp;

    std::map<double, std::vector<double>> m;
    for (auto e : f)
    {
        double h = e[ix];
        double v = e[iy];
        if (m.count(h) == 0) m[h] = std::vector<double>();
        m[h].push_back(v);
    }

    for (auto &e : m)
    {
        auto h = e.first;
        auto &vec = e.second;
        std::sort(vec.begin(), vec.end());

        double ts = vec.front();
        double te = vec.back();
        double tcut = ts + portion * (te - ts);
        auto it = std::upper_bound(vec.begin(), vec.end(), tcut,
            [](double tcut, double a){ return tcut < a; }
        );
        vec.erase(it, vec.end());
        for (auto e : vec) tmp.push_back({h, e});
    }
    assert(prg::is_sorted_by_column<n>(tmp, 0));

    return tmp;
}


#endif // DEFS_H
