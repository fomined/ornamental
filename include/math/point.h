#ifndef PRG_POINT
#define PRG_POINT

#include <array>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <random>

#include "defs.h"
#include "container/string.h"

namespace prg
{
    template <Unsigned n, typename T=double>
    class point
    {
        T v[n];

    public:

        /// инициализация компонент одним значением
        point(T value=0) { for (Unsigned i=0; i<n; i++) v[i] = value; }

        /// инициализация компонент из списка {}
        point(std::initializer_list<T> init) { Unsigned i=0; for (auto e : init) v[i++] = e; }

        /// инициализация компонент из другого объекта
        point(const point &other) { for (Unsigned i=0; i<n; i++) v[i] = other[i]; }

        /// присваивание единого значения всем компонентам
        point& operator=(T value) { for (Unsigned i=0; i<n; i++) v[i] = value; return *this; }

        /// присваивание значения объекту
        point& operator=(const point &other) { for (Unsigned i=0; i<n; i++) v[i] = other[i]; return *this; }

        /// умножение компонент на число
        point<n, T> &operator*=(T b) { for (Unsigned i=0; i<n; i++) (*this)[i] *= b; return *this; }

        /// прибавление другого вектора
        point<n, T> &operator+=(const point<n, T> &b)
        { for (Unsigned i=0; i<n; i++) (*this)[i] += b[i]; return *this; }

        /// отнимание другого вектора
        point<n, T> &operator-=(const point<n, T> &b)
        { for (Unsigned i=0; i<n; i++) (*this)[i] -= b[i]; return *this; }

        /// доступ к компонентам
        const T &operator[](Unsigned i) const { return v[i]; }
        T &operator[](Unsigned i) { return v[i]; }
        
        /// оператор смены знака
        point operator-() const { point p(*this); for (Unsigned i=0; i<n; i++) p[i] = -p[i]; return p; }

        /// число компонент
        Unsigned size() const { return n; }

        /// функции для совместимости с оператором for (auto e : point)
        const T *begin() const { return &v[0]; }
        const T *end() const { return &v[n]; }

        /// функции для совместимости с оператором for (auto &e : point)
        T *begin() { return &v[0]; }
        T *end() { return &v[n]; }
    };

	template <Unsigned  n, typename T>
    inline void print(const point<n, T> &a) 
    { 
        std::printf("[%ld]{ ", n);
        for (auto e : a) std::printf("%lf ", e);
        std::printf("}\n");
    }

    /// обменять точки друг с другом
	template <Unsigned  n, typename T>
	inline void swap(point<n, T> &a, point<n, T> &b)
        { for (Unsigned i=0; i<n; i++) std::swap(a[i], b[i]); }

    /// обменять заданные колонки внутри точки
	template <Unsigned  n, typename T>
	inline point<n, T> swap(const point<n, T> &data, Unsigned col1, Unsigned col2)
	{
		point<n, T> tmp(data);
		std::swap(tmp[col1], tmp[col2]);
		return tmp;
	}

    /// объект с противоположным знаком
	template <Unsigned  n, typename T>
    inline point<n, T> negate(const point<n, T> &p) { return -p; }

    /// сумму объектов
	template <Unsigned  n, typename T>
	inline point<n, T> operator+(const point<n, T> &a, const point<n, T> &b)
        { point<n, T> tmp(a); for (Unsigned i=0; i<n; i++) tmp[i] += b[i]; return tmp; }

    /// разность объектов
	template <Unsigned  n, typename T>
	inline point<n, T> operator-(const point<n, T> &a, const point<n, T> &b)
        { point<n, T> tmp(a); for (Unsigned i=0; i<n; i++) tmp[i] -= b[i]; return tmp; }

    /// умножение компонент на число
    template <Unsigned  n, typename T>
    inline point<n, T> operator*(const point<n, T> &a, T b)
        { point<n, T> tmp(a); for (Unsigned i=0; i<n; i++) tmp[i] *= b; return tmp; }

    /// умножение компонент на число
    template <Unsigned  n, typename T>
    inline point<n, T> operator*(T b, const point<n, T> &a) { return a * b; }

    /// скалярное произведение компонент объекта
	template <Unsigned  n, typename T>
    inline double scalar_product(const point<n, T> &a, const point<n, T> &b)
        { double sum = 0; for (Unsigned i=0; i<n; i++) sum += a[i] * b[i]; return sum; }

    /// векторное произведение компонент объекта
    inline point<3> vector_product(const point<3> &a, const point<3> &b)
    {
        point<3> out = { 0, 0, 0 };
        out[0] = a[1] * b[2] - a[2] * b[1];
        out[1] = a[2] * b[0] - a[0] * b[2];
        out[2] = a[0] * b[1] - a[1] * b[0];
        return out;
    }

    /// нормализовать длину вектора
    inline void normalize(point<3> &a) { a *= 1. / sqrt(scalar_product(a, a)); }
    inline point<3> normalize(const point<3> &a) { point<3> b(a); normalize(b); return b; }

    /// равенство объектов с заданной точностью
	template <Unsigned  n, typename T>
    inline bool equal(const point<n, T> &a, const point<n, T> &b, T eps=T(0.001))
		{ return scalar_product(a - b) < sqr(eps); }

    /// форматный вывод в строку с заданной шириной поля и числом знаков после запятой
	template <Unsigned  n, typename T>
    inline std::string to_string(const point<n, T> &a, unsigned W=6, unsigned w=3) 
    { 
        std::string s;
        for (auto e : a) s += prg::to_string(e, W, w);
        return s;
    }

	template <Unsigned  n, typename T>
    inline std::ostream &operator<<(std::ostream &os, const point<n, T> &p)
    {
        os << "[ ";
        for (auto e : p) os << e << " ";
        os << "]";
        return os;
    }


}
#endif // PRG_POINT
