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
#include <algorithm>

typedef double Real;
typedef long Int;

typedef ptrdiff_t Signed; ///< синоним ptrdiff_t
typedef size_t Unsigned; ///< синоним size_t

constexpr Unsigned nil = static_cast<Unsigned>(-1); ///< неопределенное значение индекса массива
constexpr Unsigned npos = nil; ///< не позиция (для обработки массивов)

const bool PRINT_YES = true;
const bool PRINT_NO = false;

/// "бесконечное" значение величины
template <typename T> constexpr T infinity() { return std::numeric_limits<T>::max(); }

/*!
 * \brief Описание объекта (значение переменной + флаг корректности).
 * \note Объект используется чаще всего функциями из обертки gsl для возврата результата и
 *  флага того, корректен ли этот результат.
 */
template <typename T>
struct flagged_value : public std::pair<T, bool>
{
    typedef std::pair<T, bool> base_type;
    using base_type::base_type;
    operator bool() const { return this->second; }
    operator T() const { return this->first; }
};

/*!
 * \brief пара чисел для описания интервала значений
 * \note Флагированный интервал используется функциями из обертки gsl для возврата результата и
 *  флага того, корректен ли этот результат.
 */
template <typename T>
struct interval : public std::array<T, 2>
{
    interval(T a, T b)
    {
        if (b < a) std::swap(a, b);
        std::array<T, 2>::operator[](0) = a;
        std::array<T, 2>::operator[](1) = b;
    }
    T operator[](size_t i) const { return std::array<T, 2>::operator[](i); }
        // закрыли доступ по [] для изменения, чтобы изменить нужно i = interval(a, b)

	// охватывает ли интервал заданное значение
	bool include(T v) const { return (*this)[0] <= v && v < (*this)[1]; }
};

template <typename T>
inline bool intersect(const interval<T> &a, const interval<T> &b)
{
    std::pair<T, T> u = a[0] < b[0] ? std::make_pair(b[0], a[1]) : std::make_pair(a[0], b[1]);
        // образуем пересечение интервалов, если оно не образуется, то пересечения нет

    return u.first <= u.second;
}

template <typename T> using flagged_interval = flagged_value<interval<T>>;

#define bitsizeof(U) (sizeof(U) << 3)
#define STRING(expr) #expr
#define STR std::string
#define CSTR(s) std::to_string(s).c_str()
#define PRINT(msg) { std::cout << msg << "\n"; std::cout.flush(); }
#define TOTEST { std::cout << "function must be tested\n";  throw std::exception(); }

/// предупреждение о неиспользовании незавершенного или ошибочного кода
#define TODO(msg) { std::cout << msg << "\n"; assert(false); }

/// печатать разделительную линию с заданным символом в output
#define PRINT_SEPARATOR_LINE(n, ch) printf("%s\n", std::string(n, ch).c_str());
#define SEPARATOR_LINE(n) PRINT_SEPARATOR_LINE(n, '-')
#define LINE(n) PRINT_SEPARATOR_LINE(n, '-')

/// сброс потока
#define FLUSH std::cout.flush();

#define THIS_FILE std::string(std::string(__FILE__), std::string(__FILE__).find_last_of('/') + 1)
#define THIS_LINE std::to_string(__LINE__)
#define FILE_LINE std::string(THIS_FILE + std::string("(") + THIS_LINE + std::string(")"))
#define EXIT { printf("Error at %s\n", FILE_LINE.c_str()); exit(1); }

#define RUN(state, msg) if (state) printf("RUN : %s\n", msg); if (state)
#define SELECTED_RUN(state, msg) if (state) printf("RUN : %s\n", msg); if (state)


template <typename T> inline T sqr(T a) { return a * a; }
template <typename T> inline T cube(T a) { return a * a * a; }

template <typename T> inline T mod(T x, T m)
{
    std::intmax_t quot = static_cast<std::intmax_t>(x / m);
    return static_cast<T>(x >= 0 ? x - quot * m : x + (quot + 1) * m);
}

template <typename T> inline T floor(T x, T m)
{
    std::intmax_t quot = static_cast<std::intmax_t>(x / m);
    return static_cast<T>(x >= 0 ? quot * m : (quot + 1) * m);
}

template <typename T> inline T round(T x, T m)
{
	return floor(x, m) + ((mod(x, m) <= 0.5 * m) ? 0 : m);
}

/// печать уровня завершения работы в %
inline void print_loop_completion(Unsigned i, Unsigned cnt, Unsigned output_step=1)
{
    Unsigned val = Unsigned(std::ceil(output_step * cnt / 100.));
    if (i % val == 0)
    {
        std::string s("... completed by ");
        s += std::to_string(100 * double(i)/cnt) + "%";
        std::printf("%s \n", s.c_str()); FLUSH;
    }
}


/// определяет ширину доверительного интервала в числе сигм
constexpr double confidence_interval_multiplier = 4; // охватываем > 99.99% площади под куполом
constexpr double confidence_interval_sigma = 0.001; // значение по умалчиванию

/*!
 * \brief попадает ли число в доверительный интервал, где равенство его нулю возможно
 * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
 */
inline bool is_zero(double a, double sigma=confidence_interval_sigma)
{
    assert(sigma > 0);
    return std::abs(a) <= sigma * confidence_interval_multiplier;
}

/*!
 * \brief попадает ли число в доверительный интервал, где равенство его нулю возможно
 * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
 */
inline bool equal(double a, double b, double sigma=confidence_interval_sigma)
{ return is_zero(b - a, sigma); }

/*!
 * \brief меньше ли первое число левой границы доверительного интервала второго числа,
 *  где числа с некой вероятностью могут совпадать
 * \param a,b сравниваемые значения
 * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
 */
inline bool less(double a, double b, double sigma=confidence_interval_sigma)
{
    assert(sigma > 0);
    return a < b - sigma * confidence_interval_multiplier;
}

/*!
 * \brief меньше ли первое число левой границы доверительного интервала второго числа,
 *  где числа с некой вероятностью могут совпадать
 * \param a,b сравниваемые значения
 * \param sigma критерий сравнения (чем меньше, тем строже сравнение)
 */
inline bool lessequal(double a, double b, double sigma=confidence_interval_sigma)
{
    assert(sigma > 0);
    return a <= b + sigma * confidence_interval_multiplier;
}


namespace prg
{
    using std::to_string;
    using std::numeric_limits;

    /// преобразование по умалчиванию для стандарных числовых типов
    template <typename T> inline double to_double(T a) { return a; }

    template <typename T>
    inline void print(const T &v) { std::cout << v; }

    template <typename I, typename T>
    inline void print(const std::pair<I, T> &v)
    { std::cout << "{"; print(v.first); std::cout << ","; print(v.second); std::cout << "}"; }

    template <typename T>
    /*!
     * \brief создать разделяемый указатель на объект заданного типа
     * \return
     */
    inline std::shared_ptr<T> make() { return std::shared_ptr<T>(new T()); }

    template <typename T>
    /*!
     * \brief создать разделяемый указатель на клон объекта заданного типа
     * \return
     */
    inline std::shared_ptr<T> clone(const T &x) { return std::shared_ptr<T>(new T(x)); }

    template <typename T>
    /*!
     * \brief создать разделяемый указатель на клон объекта заданного типа
     * \return
     */
    inline std::shared_ptr<T> clone(const std::shared_ptr<T> &x)
        { return std::shared_ptr<T>(new T(*x.get())); }

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

}
#endif // PRG_DEFS
