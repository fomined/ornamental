#ifndef PRG_POINT
#define PRG_POINT

#include <array>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <random>

#include "prg.defs.h"
#include "prg.string.h"

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

        /// доступ к компонентам
        const T &operator[](Unsigned i) const { return v[i]; }
        T &operator[](Unsigned i) { return v[i]; }
        
        /// оператор смены знака
        point operator-() const { point p(*this); for (Unsigned i=0; i<n; i++) p[i] = -p[i]; return p; }

        /// число компонент
        Unsigned size() const { return n; }

        /// функции для совместимости с оператором for (auto e : point)
        const double *begin() const { return &v[0]; }
        const double *end() const { return &v[n]; }

        /// функции для совместимости с оператором for (auto &e : point)
        double *begin() { return &v[0]; }
        double *end() { return &v[n]; }
    };

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

    /// скалярное произведение компонент объекта
	template <Unsigned  n, typename T>
	inline point<n, T> scalar_product(const point<n, T> &a, const point<n, T> &b)
        { double sum = 0; for (Unsigned i=0; i<n; i++) sum += a[i] * b[i]; return sum; }

    /// равенство объектов с заданной точностью
	template <Unsigned  n, typename T>
    inline bool equal(const point<n, T> &a, const point<n, T> &b, T eps=T(0.001))
		{ return scalar_product(a - b) < sqr(eps); }

//    /// сравнение объектов с заданной точностью
//	template <Unsigned  n, typename T>
//    inline bool less(const point<n, T> &a, const point<n, T> &b, T eps=T(0.001))
//		{ return equal(a, b, eps) ? false : a < b; }

    /// форматный вывод в строку с заданной шириной поля и числом знаков после запятой
	template <Unsigned  n, typename T>
    inline std::string to_string(const point<n, T> &a, unsigned W=6, unsigned w=3) 
    { 
        std::string s;
        for (auto e : a) s += prg::to_string(e, W, w);
        return s;
    }


}
#endif // PRG_POINT
