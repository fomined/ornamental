#ifndef PRG_DEFS
#define PRG_DEFS

/**
 * @file хранение стандартных сокращений и макросов
 */

#include <assert.h>
#include <string>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <iostream>
#include <limits>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <random>
#include <regex>
#include <algorithm>

//--------------------------------------------------------------------------------------------------
// Следующие прагмы удаляют назойливые сообщения, например "format string is not a string literal",
// но потенциально могут приводить к ошибкам, которые не проверяются средствами контроля clang.
// Определение этих прагм будут действовать по всему тексту программы, так как данный файл
// включается везде. Чтобы ограничить действие прагм локально, они должны быть закоментированы
// в этом файле, но определены в нужном месте в виде:
//      #pragma clang diagnostic push
//      #pragma clang diagnostic ignored "-Wformat-nonliteral"
//          program text
//      #pragma clang diagnostic pop
//
#ifdef __clang__
    #pragma clang diagnostic ignored "-Wformat-nonliteral"
    #pragma clang diagnostic ignored "-Wdisabled-macro-expansion"
    #pragma clang diagnostic ignored "-Wmissing-braces"
#endif
//--------------------------------------------------------------------------------------------------

typedef ptrdiff_t Signed; ///< синоним ptrdiff_t
typedef size_t Unsigned; ///< синоним size_t
typedef double Real;
typedef Signed Int;

constexpr Unsigned nil = static_cast<Unsigned>(-1); ///< неопределенное значение индекса массива
constexpr Unsigned npos = nil; ///< не позиция (для обработки массивов)

template <typename T>
struct integral { typedef typename std::conditional<std::is_integral<T>::value, T, void>::type type; };
    // эта конструкция не позволяет инстанцировать объект нецелого типа на этапе компиляции.
    // точнее она дает неопределенный тип, с которым работать невозможно, и это видно сразу

const bool PRINT_YES = true;
const bool PRINT_NO = false;
const bool TOTAL_COPY = true;
const bool PARTIAL_COPY = false;
const bool TOTAL_OUTPUT = true;
const bool PARTIAL_OUTPUT = false;
const Unsigned NORM_MAXIMUN = 0x1;
const Unsigned NORM_POWER = 0x2;

const double RADIAN_FROM_DEGREE = M_PI / 180;
const double DEGREE_FROM_RADIAN = 180 / M_PI;
inline double to_radian(double v) { return v * RADIAN_FROM_DEGREE; }
inline double to_degree(double v) { return v * DEGREE_FROM_RADIAN; }

/// тип указатель на вектор с элементами любого типа
template <typename T> using pvector = std::shared_ptr<std::vector<T>>;

/// вектор с элементами, являющимися указателями на объекты типа T
template <typename T> using vector_of_shared = std::vector<std::shared_ptr<T>>;

/// короткое написание для разделяемого объекта
template <typename T> using shared = std::shared_ptr<T>;
template <typename T> inline shared<T> make_shared() { return shared<T>(new T); }
template <typename T> inline shared<T> clone_shared(const T &x) { return shared<T>(new T(x)); }

/// короткое написание для парного объекта
template <typename T> using paired = std::pair<T, T>;

#define bitsizeof(U) (sizeof(U) << 3)
#define STRING(expr) #expr
#define STR std::string
#define CSTR(s) std::to_string(s).c_str()

/// сброс потока
#ifdef _DEBUG
    #define FLUSH std::cout.flush();
#else
    #define FLUSH
#endif

#define THIS_FILE std::string(std::string(__FILE__), std::string(__FILE__).find_last_of('/') + 1)
#define THIS_LINE std::to_string(__LINE__)
#define FILE_LINE std::string(THIS_FILE + std::string("(") + THIS_LINE + std::string(")"))
#define EXIT { std::printf(" Error at %s\n", FILE_LINE.c_str()); exit(EXIT_FAILURE); }

#define PRINT(msg) { std::cout << msg << "\n"; std::cout.flush(); }
#define TOTEST { std::cout << "function must be tested\n";  throw std::exception(); }

/// предупреждение о неиспользовании незавершенного или ошибочного кода с прерыванием выполнения
#define TODO_STRONG(msg) { std::cout << "TODO : " << msg << " " << "\n"; EXIT; }
#define TODO(msg) TODO_STRONG(msg)

/// предупреждение о неиспользовании незавершенного или ошибочного кода без прерывания выполнения
#define TODO_WEAK(msg) { std::cout << "\n\tTODO : " << msg << " " << FILE_LINE << "\n\n"; }

#define ERROR_CODE(msg) { std::cout << "ERROR CODE : " << msg << "\n"; std::cout << FILE_LINE << "\n"; assert(false); }
#define ERROR_CODE_WARNING(msg) { std::cout << "ERROR CODE WARNING : " << msg << "\n"; std::cout << FILE_LINE << "\n"; }
#define INEFFICIENT_CODE_WARNING(msg) { std::cout << "INEFFICIENT CODE WARNING : " << msg << "\n"; std::cout << FILE_LINE << "\n"; }

#define MSG_ASSERT(condition, msg) \
    if (!(condition)) { std::printf("failed at %s(%d) : %s\n", __FILE__, __LINE__, msg); FLUSH; exit(1); }
#define PRINT_IF(condition, msg) \
    if ((condition)) { std::printf("failed at %s(%d) : %s\n", __FILE__, __LINE__, msg); }

/// печатать разделительную линию с заданным символом в output
#define PRINT_SEPARATOR_LINE(n, ch) std::printf("%s\n", std::string(n, ch).c_str()); std::cout.flush();
#define SEPARATOR_LINE PRINT_SEPARATOR_LINE(80, '-')
#define LINE(n) PRINT_SEPARATOR_LINE(n, '-')
#define COMMENT(msg) std::printf("%s\n", msg); std::cout.flush();

#define RUN(state, msg) if (state) std::printf("RUN : %s\n", msg); if (state)
#define SELECTED_RUN(state, msg) if (state) std::printf("RUN : %s\n", msg); if (state)

/*!
 * \brief вывести значение переменной (массива ..) вместе с ее именем
 * \note значение простых переменных выводятся на той же строке, а сложных - на последующих
 *  при этом детализация вывода выполняется соответствующим оператором<< сложного типа
 */
#define COUT(f) { std::cout << STRING(f) << (std::is_fundamental<decltype (f)>::value ? " : ": " ") << (f) << "\n"; FLUSH; }

/// появление в тексте выражения TEST(+) означает включение строки в код, а выражения TEST(-) исключение строки из кода,
/// что используется для тестирования в случае необходимости
#define TEST(op) if (1 op 1 == 2)

/// добавим в std специализацию шаблонной функции std::to_string для строк, чтобы унифицировать код
namespace std { inline string to_string(string s) { return s; }}

namespace prg
{
    using std::to_string;

    /// гарантирование, что значение элемента окажется в указанных границах
    template <typename T>
    inline void range_assert(T &a, T lower, T upper)
    {
        if (a < lower) a = lower;
        if (a > upper) a = upper;
    }


    template <typename T>
    inline std::ostream &operator<<(std::ostream &os, shared<T> p) { os << *p; return os; }

    /// "бесконечное" значение величины
    template <typename T> constexpr T inf() { return std::numeric_limits<T>::max(); }
    template <typename T> constexpr T inf(T) { return std::numeric_limits<T>::max(); }
    template <typename T> constexpr T infinity() { return std::numeric_limits<T>::max(); }
    template <typename T> constexpr T infinity(T) { return std::numeric_limits<T>::max(); }
        // функция infinity(T) более удобна, чем infinity(), поскольку мы не всегда помним тип
        //  T a = infinity(a); vs. T a = infinity<decltype(a)>();


    /*!
     * \brief Описание объекта (значение переменной + флаг корректности).
     * \note Объект используется чаще всего функциями из обертки gsl для возврата результата и
     *  флага того, корректен ли этот результат.
     */
    template <typename T>
    struct flagged_value
    {
        T first;
        bool second;

        flagged_value() : first(), second(false) {}
        flagged_value(T v=0, bool b=false) : first(v), second(b) {}

        T value() const { return first; }
        bool flag() const { return second; }

        T *operator->() { return &first; }
        T &operator*() { return first; }

        operator bool() const { return second; }
        operator T&() const { return first; }
        operator T&() { return first; }
    };

    template <typename T>
    inline std::string to_string(const flagged_value<T> &f)
    {
        std::string s;
        s += "{" + std::to_string(f.value()) + "," + std::to_string(f.flag()) + "}";
        return s;
    }

    /*!
     * \brief Описание объекта (значение переменной + число включений[+]/исключений[-] объекта).
     * \note Класс используется для описания есть или нет объекты и сколько их
     */
    template <typename T>
    struct counted_value : public std::pair<T, Signed>
    {
        typedef std::pair<T, Signed> base_type;
        using base_type::base_type;

        T value() const { return this->first; }
        Signed count() const { return this->second; }

        operator bool() const { return this->second > 0; }
        operator T&() const { return this->first; }
        operator T&() { return this->first; }
    };

    /// определяет минимальную разность между значениями заданного типа
    /// нужна для определения границ сумм и разностей интервалов
    template <typename T> constexpr typename T::difference_type minimal_difference(T);

    inline constexpr Signed minimal_difference(Signed) { return 1; }
    inline constexpr Unsigned minimal_difference(Unsigned) { return 1; }
    inline constexpr double minimal_difference(double) { return 0; }

    /*!
     * \brief итератор по диапазону целых значений
     * \note нужен, что удобства использования в for (auto j : interval)
     */
    template <typename T> class interval_iterator
    {
        typename integral<T>::type i;

    public:
        interval_iterator(T i) : i(i) {}
        interval_iterator *operator++() { i++; return this; }
        T operator*() { return i; }
        bool operator!=(interval_iterator v) { return i!=v.i; }
    };

    /*!
     * \brief пара чисел для описания интервала значений
     * \note интервал локализации числа, верхнее значение НЕРАЗРЕШЕНО
     * \note Флагированный интервал используется функциями из обертки gsl для возврата результата и
     *  флага того, корректен ли этот результат.
     */
    template <typename T>
    struct interval : public std::pair<T, T>
    {
        interval(T a=T(), T b=T()) { this->first = a; this->second = b; }

        bool invalid() const { return this->first > this->second; }
        bool empty() const { return this->first == this->second; }

        T &left() { return this->first; }
        T &right() { return this->second; }

        const T &left() const { return this->first; }
        const T &right() const { return this->second; }

        /// длина интервала для вещественного типа
        T length() const { return this->second - this->first; }

        /// длина интервала для целочисленного типа
        Unsigned size() const { return (Unsigned)(this->second - this->first); }

        // введены для удобства использования в for (auto j : interval)
        interval_iterator<T> begin() const { return interval_iterator<T>(this->first); }
        interval_iterator<T> end() const { return interval_iterator<T>(this->second); }

        // охватывает ли интервал заданное значение
        bool include(T v) const { return this->first <= v && v < this->second; }
    };

    template <typename T>
    inline std::string to_string(const interval<T> &f)
    {
        std::string s;
        s += "{" + std::to_string(f.left()) + "," + std::to_string(f.right()) + "}";
        return s;
    }

    template <typename T>
    inline interval<T> operator+(const interval<T> &a, const interval<T> &b)
        { return {a.left() + b.left(), a.right() + b.right()}; }

    template <typename T>
    inline interval<T> operator-(const interval<T> &b, const interval<T> &a)
        { return {b.left() - a.right(), b.right() - a.left()}; }

    /*!
     * \brief получить пересечение интервалов
     * \note если полученный интервал недействительный, то пересечения нет
     */
    template <typename T>
    inline interval<T> intersect(const interval<T> &a, const interval<T> &b)
        { return {std::max(a.first, b.first), std::min(a.second, b.second)}; }


    /// интервал с разными соглашениями по границам ([true] стандартные, [false] реверсионные
    /// [0] реверсионные [1] стандартные [2] обе включены
    template <typename T, Unsigned bound> struct Interval;

    template <typename T>
    struct Interval<T, 0> : public interval<T> // реверсионный интервал )a, b]
    {
        typedef interval<T> Base;
        typedef decltype(minimal_difference(T())) minimal_difference_type;

        using Base::invalid;
        using Base::empty;
        using Base::left;
        using Base::right;
        using Base::length;
        using Base::size;
        using Base::interval;

        /// минимальное и максимальное разрешенные значения
        T minimal_value() const { return left() + min_diff(); }
        T maximal_value() const { return right(); }
        minimal_difference_type min_diff() const { return minimal_difference(T()); }

        // охватывает ли интервал заданное значение
        bool include(T v) const { return minimal_value() <= v && v <= maximal_value(); }

        // введены для удобства использования в for (auto j : interval)
        interval_iterator<T> begin() const { return interval_iterator<T>(this->first + min_diff()); }
        interval_iterator<T> end() const { return interval_iterator<T>(this->second + min_diff()); }
    };

    template <typename T>
    inline Interval<T,0> operator-(const Interval<T,0> &b, const Interval<T,0> &a)
        { return {b.left() - a.right(), b.right() - a.left() - minimal_difference(T())}; }

    template <typename T>
    inline std::string to_string(const Interval<T, 0> &f)
    {
        std::string s;
        s += "(" + to_string(f.left()) + "," + to_string(f.right()) + "]";
        return s;
    }

    template <typename T>
    struct Interval<T, 1> : public interval<T> // стандартный интервал [a, b(
    {
        typedef interval<T> Base;
        typedef decltype(minimal_difference(T())) minimal_difference_type;

        using Base::invalid;
        using Base::empty;
        using Base::left;
        using Base::right;
        using Base::length;
        using Base::size;
        using Base::interval;

        /// минимальное и максимальное разрешенные значения
        T minimal_value() const { return left(); }
        T maximal_value() const { return right() - min_diff(); }
        minimal_difference_type min_diff() const { return minimal_difference(T()); }

        // охватывает ли интервал заданное значение
        bool include(T v) const { return minimal_value() <= v && v <= maximal_value(); }

        // введены для удобства использования в for (auto j : interval)
        interval_iterator<T> begin() const { return interval_iterator<T>(this->first); }
        interval_iterator<T> end() const { return interval_iterator<T>(this->second); }
    };

    template <typename T>
    inline std::string to_string(const Interval<T, 1> &f)
    {
        std::string s;
        s += "[" + to_string(f.left()) + "," + to_string(f.right()) + ")";
        return s;
    }

    template <typename T>
    struct Interval<T, 2> : public interval<T> // стандартный интервал [a, b] для распределений
    {
        typedef interval<T> Base;
        typedef decltype(minimal_difference(T())) minimal_difference_type;

        using Base::invalid;
        using Base::empty;
        using Base::left;
        using Base::right;
        using Base::length;
        using Base::size;
        using Base::interval;

        /// минимальное и максимальное разрешенные значения
        T minimal_value() const { return left(); }
        T maximal_value() const { return right(); }
        minimal_difference_type min_diff() const { return minimal_difference(T()); }

        // охватывает ли интервал заданное значение
        bool include(T v) const { return minimal_value() <= v && v <= maximal_value(); }

        // введены для удобства использования в for (auto j : interval)
        interval_iterator<T> begin() const { return interval_iterator<T>(this->first); }
        interval_iterator<T> end() const { return interval_iterator<T>(this->second + min_diff()); }
    };

    template <typename T>
    inline std::string to_string(const Interval<T, 2> &f)
    {
        std::string s;
        s += "[" + to_string(f.left()) + "," + to_string(f.right()) + ")";
        return s;
    }

    /*!
     * \brief значение, определенное достаточно неточно
     * \note оно задается предположительным значением и интервалом, в котором оно локализовано
     */
    template <typename T>
    class rough_value
    {
        typedef Interval<T, 0> interval_type;
        interval_type r_; // интервал локализации числа, оба значения разрешены как для ran0(..)
            // ВНИМАНИЕ, соглашения по границам для interval и rough_value РАЗНЫЕ

        T x_; // предположительное значение

    public:

        /// установка границ интервала и точки смены фазы в нем
        rough_value() {}
        rough_value(const interval_type &r) : r_(r) { x_ = r_.left() + (r_.right() - r_.left()) * 0.5; }
        rough_value(T x, const interval_type &r) : r_(r), x_(x) { assert (r_.include(x_)); }
        rough_value(const interval_type &r, T x) : r_(r), x_(x) { assert (r_.include(x_)); }

        /// получение/установка значения
        T value() const { return x_; }
        T &value() { return x_; }

        rough_value &operator=(const rough_value &r) { r_ = r.r_; x_ = r.x_; return *this; }

        /// получение границ интервала
        interval_type range() const { return r_; }
        interval_type &range() { return r_; }

        T left() const { return r_.first; }
        T right() const { return r_.second; }

        T &left() { return r_.first; }
        T &right() { return r_.second; }

        T minimal_value() const { return r_.minimal_value(); }
        T maximal_value() const { return r_.maximal_value(); }

        // охватывает ли интервал заданное значение
        bool include(T v) const { return r_.include(v); }

    };

    template <typename T>
    inline rough_value<T> operator+(const rough_value<T> &b, const rough_value<T> &a)
        { return {b.value() + a.value(), b.range() + a.range()}; }

    template <typename T>
    inline rough_value<T> operator-(const rough_value<T> &b, const rough_value<T> &a)
        { return {T(b.value() - a.value()), b.range() - a.range()}; }

    template <typename T>
    inline std::ostream &operator<< (std::ostream &os, const rough_value<T> &r)
        { os << ']' << r.left() << ',' << r.value() << ',' << r.right() << ']'; return os; }

    template <typename T>
    inline std::string to_string(const rough_value<T> &r)
        { std::string s; s += to_string(r.value()) + to_string(r.range()); return s; }

    /// синоним std::pair<T, T>
    template <typename T> struct twain : public std::pair<T, T>
    {
        twain(const T &a, const T &b) { this->first = a; this->second = b; }
    };

    template <typename T>
    inline twain<T> operator+(const twain<T> &a, const twain<T> &b)
        { return {a.first + b.first, a.second + b.second}; }

    inline void FILE_SAVED(const std::string &filename)
        { std::cout << filename << " is saved\n"; std::cout.flush(); }
    inline void FILE_OPEN(const std::string &filename)
        { std::printf("\'%s\' is open\n", filename.c_str()); std::cout.flush(); }
    inline void FILE_OPEN(const std::string &filename, Unsigned cnt=1)
        { std::printf("\'%s\' [%ld] is open\n", filename.c_str(), cnt); std::cout.flush(); }

    template <typename T> inline T sqr(T a) { return a * a; }
    template <typename T> inline T cube(T a) { return a * a * a; }
    template <typename T> inline T floor(T x, T m)
    {
        std::intmax_t quot = static_cast<std::intmax_t>(x / m);
        return static_cast<T>(x >= 0 ? quot * m : (quot - 1) * m);
    }
    template <typename T> inline T ceil(T x, T m) { return floor(x, m) + m; }
    template <typename T> inline T mod(T x, T m) { return x - floor(x, m); }
    template <typename T> inline T round(T x, T m)
        { return floor(x, m) + ((mod(x, m) <= 0.5 * m) ? 0 : m); }

    /// загрузка с заданного адреса числа или строки
    template <typename T> inline T load(const char *p);
    template <> inline int load<int>(const char *p) { return std::atoi(p); }
    template <> inline long load<long>(const char *p) { return std::atol(p); }
    template <> inline Unsigned load<Unsigned>(const char *p) { return std::atol(p); }
    template <> inline float load<float>(const char *p) { return float(std::atof(p)); }
    template <> inline double load<double>(const char *p) { return std::atof(p); }
    template <> inline std::string load<std::string>(const char *p) { return std::string(p); }


    /// общая операция явного преобразования типа
    template <typename T, typename F> T explicit_cast(F);

    /// конвертирование enum в целые числа
    template <typename E> inline constexpr Signed to_signed(E e) noexcept
    { return static_cast<typename std::underlying_type<E>::type>(e); }

    template <typename E> inline constexpr Unsigned to_unsigned(E e) noexcept
    { return static_cast<typename std::underlying_type<E>::type>(e); }

    template <typename E>
    inline constexpr typename std::underlying_type<E>::type value(E e) noexcept
    { return static_cast<typename std::underlying_type<E>::type>(e); }

    /// получение следующего значения enum (без контроля на допустимость)
    template <typename E> inline E next_enum_value(E e) noexcept
    {
        auto p = static_cast<typename std::underlying_type<E>::type>(e);
        return static_cast<E>(++p);
    }

    /// получение предыдущего значения enum (без контроля на допустимость)
    template <typename E> inline E prev_enum_value(E e) noexcept
    {
        auto p = static_cast<typename std::underlying_type<E>::type>(e);
        return static_cast<E>(--p);
    }

    // шаблонные функции для enum классов, такие как:
    // template <typename E> inline E operator+(E e, Signed i);
    // template <typename E> inline E operator-(E e, Signed i);
    // недопустимы, так как такие классы не отличаются от простых типов,
    // что приводит к ошибкам компиляции простых выражений, 1+2 ...

    /// преобразование по умалчиванию для стандарных числовых типов
    template <typename T> inline double to_double(T a) { return a; }

    template <typename T>
    inline void print(const T &v) { std::cout << v; }

    template <typename I, typename T>
    inline void print(const std::pair<I, T> &v)
    { std::cout << "{"; print(v.first); std::cout << ","; print(v.second); std::cout << "}"; }

    /// форматный вывод в строку, аналог sprintf(char *, const char *, ..)
    inline void printf(const char *s) { std::cout << s; std::cout.flush(); }

    template<typename T, typename... Args>
    void printf(const char *s, T value, Args... args)
    {
        while (*s)
        {
            if (*s == '%' && *(++s) != '%')
            {
                char argf[64];

                const char *p = s; int i = 0;
                argf[i++] = *p; p++;
                while (!isspace(*p) && *p!=0 && *p!='%')
                {
                    argf[i++] = *p;
                    ++p;
                }
                argf[i] = '\0';

                char argv[1024];
                std::sprintf(argv, argf, value);
                std::cout << argv;

                s = p;
                printf(s, args...);
                return;
            }
            std::cout << *s++;
        }
    }

    /*!
     * \brief печать уровня завершения работы в %
     * \param i текущий индекс прохождения массива
     * \param n максимальный индекс массива
     * \param h доля отработанных записей (от 1) при которой происходит вывод
     */
    inline void print_loop_completion(Unsigned i, Unsigned n, double h=0.1)
    {
        static Unsigned target_bound = ceil(i, Unsigned(n*h));

        if (i > target_bound)
        {
            std::printf("... completed by %lf%%\n", 100*double(i)/n);
            std::cout.flush();
            target_bound = ceil(i, Unsigned(n*h));
        }
    }

    //---------------------------------------------------------------------------------------------
    //                              вывод в поток сложных типов
    //---------------------------------------------------------------------------------------------
    template <typename T, typename Q>
    inline std::ostream &operator<<(std::ostream &os, const std::pair<T, Q> &p)
    {
        os << "[ " << p.first << " " << p.second << "]";
        return os;
    }

    #define OUTPUT(os, p) { os << "[\n"; for (auto e : p) os << e << "\n"; os << "]\n";  \
        std::cout.flush(); return os; }

    template <typename T, Unsigned n>
    inline std::ostream &operator<<(std::ostream &os, const std::array<T, n> &p) OUTPUT(os, p)

    template <typename T>
    inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &p) OUTPUT(os, p)

    template <typename T>
    inline std::ostream &operator<<(std::ostream &os, const std::set<T> &p) OUTPUT(os, p)

    template <typename K, typename T>
    inline std::ostream &operator<<(std::ostream &os, const std::map<K, T> &p) OUTPUT(os, p)

    #undef OUTPUT


    /// для упрощения поиска регулярных выражений
    inline flagged_value<std::string> regex_search(const std::string &target, const std::string &pattern)
    {
        std::smatch m;
        std::regex_search(target, m, std::regex(pattern));
        std::string s = m.empty() ? "" : std::string(m[0]);
        return flagged_value<std::string>(s, !s.empty());
    }

}
#endif // PRG_DEFS
