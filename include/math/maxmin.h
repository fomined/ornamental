#ifndef PRG_MAXMIN_H
#define PRG_MAXMIN_H

/**
 * @file prg_maxmin.h
 * @brief операции накопления максимальных/минимальных чисел
*/

#include <string>
#include <ostream>
#include <utility>
#include <map>
#include <deque>
#include <vector>
#include <functional>
#include <utility>
#include <algorithm>
#include <fstream>

#include "defs.h"

namespace prg
{

    /**
     * @brief Объект, который "накапливает" заданное число значений
     */
    template <typename T>
    class accumulator : public std::deque<T>
    {
        typedef std::deque<T> Base;
        Unsigned N; ///< capacity

    public:

        accumulator (Unsigned N) : N(N) {}
        Unsigned capacity() const { return N; }

        void push(const T &v)
        {
            Base::push_front(v);
            if (Base::size() > N) Base::pop_back();
        }

        T top() const { return operator[](0); }
        std::string operator[](Signed k) const
        {
            if (k < 0) k = 0;
            if (k >= (Signed)Base::size()) k = Base::size() - 1;
            return Base::operator[](k);
        }

        std::string &operator[](Signed k)
        {
            if (k < 0) k = 0;
            if (k >= (Signed)Base::size()) k = Base::size() - 1;
            return Base::operator[](k);
        }
    };

    /**
     * @brief Объект, который "накапливает" заданное число значений (без последовательных повторений)
     */
    template <typename T>
    class History : public accumulator<T>
    {
        typedef accumulator<T> Base;

    public:
        using Base::top;
        using Base::operator[];

        History(Unsigned N) : Base(N) {}

        void push(const T &v)
        {
            if (v != Base::front()) Base::push_front(v);
            if (Base::size() > Base::capacity()) Base::pop_back();
        }

        T top() const { return operator[](0); }
    };

    /*!
     * \brief накопитель заданного числа значений, которые максимально близко совпадают с заданным
     * \note для поиска минимума нужно задать -inf, максимума +inf, заданного значения - некое число
     * \param T тип значения
     * \param D операция, вычисляющая расстояние до заданного значения
     */
    template <typename T, class D>
    class concentrator : public std::multimap<double, T>
    {
        typedef std::multimap<double, T> storage_type;
            // обязательно multimap, поскольку могут прийти значения с одинаковыми ключами, но
            // ссылающиеся на разные объекты

        T sample; ///< образец для сравнения
        Unsigned N; ///< максимальное число элементов в концентраторе
        double maxd; ///< текущее максимальное расстояние

    public:
        typedef std::pair<double, T> value_type;
        typedef T mapped_type;

        using storage_type::begin;
        using storage_type::end;
        using storage_type::size;

        /// создание пустого концентратора
        concentrator(T sample=T(), Unsigned N=1) : sample(sample), N(N),
            maxd(std::numeric_limits<double>::max())
        {
            assert(N > 0); // нет смысла хранить и делать обработку 0 объектов в функциях вставки
        }

        /// максимальная емкость концентратора
        Unsigned capacity() const { return N; }

        /// изменить число хранимых объектов
        void resize(Unsigned n)
        {
            assert(n > 0); // нет смысла хранить и делать обработку 0 объектов в функциях вставки
            Unsigned num = storage_type::size();
            if (num <= n) N = n;
            else
            {
                auto it = storage_type::end();
                Unsigned dn = num - n;
                for (Unsigned i=0; i<dn; i++) --it;

                storage_type::erase(it);
                maxd = (--storage_type::end())->first;
            }
        }

        /// операция вставки нового значения
        void insert(const T &v)
        {
            double d = D()(sample, v); // текущее расстояние

            if (storage_type::size() < N)
            {
                if (d > maxd) maxd = d;
                storage_type::insert(value_type(d, v));
            }
            else if (d < maxd)
            {
                storage_type::insert(value_type(d, v));
                storage_type::erase(--storage_type::end());
                maxd = (--storage_type::end())->first;
            }
        }

        template <typename Operation>
        void insert(const T &v, Operation op)
        {
            double d = op(sample, v); // текущее расстояние

            if (storage_type::size() < N)
            {
                if (d > maxd) maxd = d;
                storage_type::insert(value_type(d, v));
            }
            else if (d < maxd)
            {
                storage_type::insert(value_type(d, v));
                storage_type::erase(--storage_type::end());
                maxd = (--storage_type::end())->first;
            }
        }

        /*!
         * \brief сохранить данные в файл
         * \param filename имя файла
         * \note для возможности сброса данных у mapped_type должна быть операция сброса в поток
         */
        void save(const std::string &filename, std::ios::openmode mode=std::ios::app);
    };

    template <typename T, class D>
    void concentrator<T,D>::save(const std::string &filename, std::ios::openmode mode)
    {
        std::ofstream file(filename, mode);
        if (!file) { std::printf("can't open file \'%s\'", filename.c_str()); return; }
        for (auto e : *this)
        {
            mapped_type m = e.second;
            m.save(file); // требуется операция сброса в поток
            //file << m << "\n"; // требуется операция сброса в поток
        }
    }

    /**
     * @brief Объект, который "накапливает" заданное число максимальных значений
     */
    template <typename T, class Operation>
    class beam_accumalator
    {
        std::vector<T> data; ///<
        Unsigned N; ///< capacity
        Unsigned n; ///< number of elements inside
        Signed p; ///< insert position
        T min, max; ///< текущие минимальное и максимальное значения в beam векторе

        void ordering_insert(const T &v)
        {
            assert(n >= 1);

            data[p] = v; Signed k = n-1;
            Unsigned i = p, j = (p + 1) % N;
            while (k-- > 0 && Operation()(data[j], data[i]))
            {
                std::swap(data[i], data[j]);
                i = j; j = (j + 1) % N;
            }
            if (--p < 0) p = N-1;
        }

    public:

        beam_accumalator(Unsigned N) : data(N, T()), N(N), n(0), p(N-1) {}
        Unsigned size() const { return n; }
        Unsigned capacity() const { return N; }

        T minimal() const { return min; }
        T maximal() const { return max; }

        void insert(const T &v)
        {
            if (n == 0) { min = v; max = v; data[p] = v; if (--p < 0) p = N-1; n++; }
            else if (n < N)
            {
                if (Operation()(max, v)) { data[(p+n+1) % N] = v; max = v; n++; }
                else if (Operation()(v, min)) { data[p] = v; if (--p < 0) p = N-1; min = v; n++; }
                else { ordering_insert(v); n++; }
            }
            else
            {
                if (Operation()(max, v)) ;
                else if (Operation()(v, min))
                {
                    data[p] = v; if (--p < 0) p = N-1;
                    min = v; max = data[p];
                }
                else
                {
                    ordering_insert(v);
                    max = data[p];
                }
            }
        }
        const T &operator[](Unsigned i) const { return data[(i+p+1) % N]; }
    };

    /**
     * @brief Объект, который "накапливает" заданное число максимальных значений
     */
    template <typename T>
    class beam_min
    {
        std::vector<T> data; ///<
        Unsigned N; ///< capacity
        Unsigned n; ///< number of elements inside
        Signed p; ///< insert position
        T min, max; ///< текущие минимальное и максимальное значения в beam векторе

        void ordering_insert(const T &v)
        {
            assert(n >= 1);

            data[p] = v; Signed k = n-1;
            Unsigned i = p, j = (p + 1) % N;
            while (k-- > 0 && data[j] < data[i])
            {
                std::swap(data[i], data[j]);
                i = j; j = (j + 1) % N;
            }
            if (--p < 0) p = N-1;
        }

    public:

        beam_min(Unsigned N) : data(N, T()), N(N), n(0), p(N-1) {}
        Unsigned size() const { return n; }
        Unsigned capacity() const { return N; }

        T minimal() const { return min; }
        T maximal() const { return max; }

        void insert(const T &v)
        {
            if (n == 0) { min = v; max = v; data[p] = v; if (--p < 0) p = N-1; n++; }
            else if (n < N)
            {
                if (max < v) { data[(p+n+1) % N] = v; max = v; n++; }
                else if (v < min) { data[p] = v; if (--p < 0) p = N-1; min = v; n++; }
                else { ordering_insert(v); n++; }
            }
            else
            {
                if (max < v) ;
                else if (v < min) { data[p] = v; if (--p < 0) p = N-1; min = v; max = data[p]; }
                else
                {
                    ordering_insert(v);
                    max = data[p];
                }
            }
        }
        const T &operator[](Unsigned i) const { return data[(i+p+1) % N]; }
    };

    /**
     * @brief Объект, который "накапливает" экстремальное значение
     *
     * @tparam U тип значения
     * @tparam Info дополнительная информация об экстремуме
     * @tparam Operation операция подмены (если внутрений объект удовлетворяет требованию,
     *  то он заменяется на новый).
     *
     * @note Операция подмены может быть сложной, и включать дополнительные проверки например,
     *  на вхождение в интервал.
     */
    template <typename U, typename Info, template <typename> class Operation=std::less>
    class extremum : public std::pair<U, Info>
    {
        typedef std::pair<U, Info> base_type;

    public:

        typedef std::pair<U, Info> value_typeype;

        extremum(U u, Info s) : base_type(u, s){}
        extremum(const extremum &) = default;

        /// вставить новое число в конвейер
        extremum &operator <<(const std::pair<U,Info> &u)
        {
            if (Operation<U>()(base_type::first, u.first)) (base_type &)*this = u;
            return *this;
        }

        /// вставить новое число в конвейер
        extremum &operator <<(const extremum &u)
        {
            if (Operation<U>()(base_type::first, u.first)) *this = u;
            return *this;
        }

        /// вывод содержания конвейера в строку
        std::string to_string() const;
    };

    template <typename U, typename Info, template <typename> class Operation>
    inline std::string extremum<U, Info, Operation>::to_string() const
    {
        std::string s;
//        s += boost::lexical_cast<std::string>(base_type::first) + std::string("[");
//        s += boost::lexical_cast<std::string>(base_type::second) + std::string("]");
        return s;
    }

    template <typename U, typename Info, template <typename> class Operation>
    inline std::ostream &operator <<(std::ostream &os, const extremum<U, Info, Operation> &s)
        { return os << s.to_string(); }

    /// специализация extremum для случая отсутствия дополнительной информации
    template <typename U, template <typename> class Operation>
    class extremum<U, void, Operation>
    {
        U first;

    public:

        extremum(U u) : first(u){}
        extremum(const extremum &) = default;

        /// вставить новое число в конвейер
        extremum &operator <<(const U &u)
            { if (Operation<U>()(first, u.first)) first = u; return *this; }

        /// вставить новое число в конвейер
        extremum &operator <<(const extremum &u)
            { if (Operation<U>()(first, u.first)) first = u; return *this; }

        /// приведение к значению
        operator U() const { return first; }

        /// вывод содержания конвейера в строку
        std::string to_string() const;
    };

    template <typename U, template <typename> class Operation>
    inline std::string extremum<U, void, Operation>::to_string() const
    {
        std::string s;//(boost::lexical_cast<std::string>(first));
        return s;
    }

    template <typename U, template <typename> class Operation>
    inline std::ostream &operator <<(std::ostream &os, const extremum<U, void, Operation> &s)
        { return os << s.to_string(); }

    /**
     * @brief Объект, который "накапливает" максимальное значение
     * @tparam U тип значения
     * @tparam Info дополнительная информация об максимуме
     */
    template <typename U, typename Info>
    using maximum = extremum<U, Info, std::less>;

    /**
     * @brief Объект, который "накапливает" минимальное значение
     * @tparam U тип значения
     * @tparam Info дополнительная информация об минимуме
     */
    template <typename U, typename Info>
    using minimum = extremum<U, Info, std::greater>;

}
#endif // PRG_MAXMIN_H
