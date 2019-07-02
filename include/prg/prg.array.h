#ifndef PRG_ARRAY_H
#define PRG_ARRAY_H

#include <array>
#include <exception>
#include <assert.h>

namespace prg
{

    /*!
     * \brief имитация класса std::array<T, N> с более удобным интерфейсом как в std::vector
     *
     * Целью разработки данного класса являлось избежать неудобств std::array<T, N> связанных
     * с тем, что тот имеет всегда фиксированное число элементов =N, а не разное число элементов
     * в диапазоне [0,N], и при этом сохранить преимущество этого класса, связанное с размещением
     * его в памяти.
     *
     * \note версия с фиксированным числом элементов в массиве (рабочая для проекта ornament)
     * \note позволяет (в отличие от std::array) передачу как значения в std::map
     *
     */
    template <typename T, size_t N>
    class array : public std::array<T, N>
    {
        typedef std::array<T, N> base_t;
    public:

        typedef typename base_t::const_iterator const_iterator;
        typedef typename base_t::iterator iterator;

        using base_t::front;
        using base_t::back;
        using base_t::operator[];
        using base_t::begin;
        using base_t::end;

        array(const T &value = T()) { base_t::fill(value); }
        array(const array &other) : array(other.begin(), other.end()) {}
        array(array &&other) { swap(other); }
        array(std::initializer_list<T> init) : array(init.begin(), init.end()) {}
        template<class Iter> array(Iter first, Iter last) { assign(first, last); }

        array& operator=(const array &other) { assign(other.begin(), other.end()); return *this; }
        array& operator=(array &&other) { swap(other); return *this; }
        array& operator=(std::initializer_list<T> init) { assign(init.begin(), init.end()); return *this; }

        template <class Iter> void assign(Iter first, Iter last)
        {
            for (size_t i=0; i<N; i++) (*this)[i] = 0;
            size_t i = 0; for (; first != last && i < N; ++first, i++) (*this)[i] = *first;
        }

        void swap(array &other) noexcept { base_t::swap(other); }

        size_t size() const { return N; }
        void clear() { base_t::fill(0); }
    };



//    template <typename T, size_t N>
//    class array : public std::array<T, N>
//    {
//        typedef std::array<T, N> base_t;
//        size_t n;

//    public:

//        typedef typename base_t::const_iterator const_iterator;
//        typedef typename base_t::iterator iterator;

//        using base_t::front;
//        //using base_t::operator[];
//        using base_t::data;
//        using base_t::begin;

//        array() : n(0) {}
//        array(size_t count, const T &value = T()) : n(count) { assert(count <= N); base_t::fill(value); }
//        array(const array &other) : array(other.begin(), other.end()) {}
//        array(array &&other) { swap(other); }
//        array(std::initializer_list<T> init) : array(init.begin(), init.end()) {}
//        template<class Iter> array(Iter first, Iter last) { assign(first, last); }

//        array& operator=(const array &other) { assign(other.begin(), other.end()); return *this; }
//        array& operator=(array &&other) { swap(other); return *this; }
//        array& operator=(std::initializer_list<T> init) { assign(init.begin(), init.end()); return *this; }

//        T operator[](size_t i) const { assert(i<n); return base_t::operator[](i); }
//        T& operator[](size_t i) { assert(i<n); return base_t::operator[](i); }

//        void assign(size_t count, const T &value) { n = count; assert(n <= N); base_t::fill(value); }
//        void assign(std::initializer_list<T> init) { assign(init.begin(), init.end()); }

//        // через эту функцию назначаются значения
//        template <class Iter> void assign(Iter first, Iter last)
//        {
//            size_t i = 0;
//            for (; first != last && i < N; ++first, i++)
//                base_t::operator[](i) = *first;
//            n = i;
//        }

//        constexpr T& at(size_t i) const { if (i >= size()) throw std::out_of_range(""); return (*this)[i]; }
//        T& at(size_t i) { if (i >= size()) throw std::out_of_range(""); return (*this)[i]; }
//        constexpr T& back() const { return (*this)[n-1]; }
//        T& back() { return (*this)[n-1]; }

//        iterator end() noexcept { return begin() + n; }
//        const_iterator end() const noexcept { return begin() + n; }

//        bool empty() const { return n == 0; }
//        size_t size() const { return n; }
//        size_t max_size() const { return N; }

//        void swap(array &other) noexcept { base_t::swap(other); std::swap(n, other.n); }
//        void clear() { n = 0; }

//        iterator insert(const_iterator pos, const T &value);
//        iterator insert(const_iterator pos, T &&value);
//        iterator insert(const_iterator pos, size_t count, const T& value );
//        template<class Iter> iterator insert(const_iterator pos, Iter first, Iter last);
//        iterator insert(const_iterator pos, std::initializer_list<T> init);

//        template<class... Args> iterator emplace(const_iterator pos, Args&&... args);

//        iterator erase(const_iterator pos);
//        iterator erase(const_iterator first, const_iterator last);

//        void push_back(const T &value)
//        {
//        assert(n<N);
//        base_t::operator[](n++) = value;
//        }
//        //void push_back(T &&value) { assert(n<N); base_t::operator[](n++) = value;  }

//        template<class... Args> void emplace_back(Args&&... args);
//        void pop_back() { n--; }

//        void resize(size_t count, const T &value = T())
//        {
//            assert(count <= N);
//            while (n < count) (*this)[n++] = value;
//            n = count; // сработает для случая, если исходно count < n
//        }
//    };

}
#endif // PRG_ARRAY_H

