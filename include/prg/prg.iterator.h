/**
* \file _iterator.h
*
* \brief Итераторы для виртуальных контейнеров (вычисление элемента вместо его извлечения)
*	и для контейнеров, не входящих в STL (включая эти контейнеры).
*
*	При написании итераторов для указанных целей следует понять, что требование полной
*	совместимости итераторов с STL является плохой идеей. Уровень совместимости должен
*	определяться поставленными перед контейнерами и их итераторами задачами.
*
*	Совместимость с STL можно понимать следующим образом:
*	- совместимость со всеми контейнерами и алгоритмами STL (наивный подход),
*	- совместимость со всеми или отдельными алгоритмами STL (образцовый подход),
*	- совместимости на уровне схожести интерфейса (ленивый подход).

*	Наивность первого варианта обусловлена забыванием того факта, что любой итератор пишется
*	для возможности прохождения своего типа конвейера (наприсер, std::vector<_T>::iterator для
*	вектора или std::list<_T>::iterator для списка). Вновь написанный итератор не должен быть
*	совместим с уже написанными контейнерами, так у них уже есть свои итераторы.
*
*	Образцовый подход следует использовать только если новый контейнер действительно нуждается
*	в алгоритмах, предложенных STL, что включает сортировки, поиск, трансформации и т.д. Если
*	контейнер является виртуальным, как например "полная последовательность целых чисел"
*	(с бесконечным числом элементов), то есть который не реализуется как программная сущность,
*	а всего лишь "представляется в уме", то никакой STL алгоритм к нему неприменим. Здесь
*	использование образцового подхода является бессмысленной тратой сил и времени.
*
*	Итого, разумный подход связан с анализом необходимости уровня совместимости и написанием
*	кода под полученный результат анализа.
*
*	Вот что (для справки) должно быть описано в совместимом с STL итераторе:
*	\code
*	template <class Iterator> struct iterator_traits {
*		typedef typename Iterator::iterator_category iterator_category;
*		typedef typename Iterator::value_type        value_type;
*		typedef typename Iterator::difference_type   difference_type;
*		typedef typename Iterator::pointer           pointer;
*		typedef typename Iterator::reference         reference;
*	};
*	\endcode
*
*	либо при наследовании от std::iterator:
*	\code
*	template <typename _Category, typename _Tp, typename _Distance = ptrdiff_t,
*		typename _Pointer = _Tp*, typename _Reference = _Tp&>
*	\endcode
*
*	Заметим, что наследование практически бесполезно. Оно работает только для классов, но
*	не для шаблонов (требуется заново определять типы).
*/
#ifndef ITERATOR_H
#define ITERATOR_H

#include <iterator>
#include <vector>
#include <type_traits>
#include <cstddef>

/**
* \def MEMBER_ITERATOR_IS_FULL_STL_COMPATIBILE
* \brief Полная совместимость с STL для _member_iterator
* \note Полная совместимость с STL для _member_iterator нужна для экзотических случаев,
*	например, для сортировки отдельных член-данных в контейнерах объектов (что, скорее всего,
*	нарушает согласованность объектов). По умалчиванию значение равно 0 (false).
*/
#define MEMBER_ITERATOR_IS_FULL_STL_COMPATIBILE 0

namespace prg
{

	/**
	* \brief Итератор, проходящий всю последовательность целых чисел.
	*
	*	Хранить бесконечную последовательность целых чисел бессмысленно. Так что такой
	*	контейнер является виртуальным: реализации его нет, но извлекать данные из него возможно.
	*	Для такого контейнера уровень совместимости с STL на уровне схожести интерфейса и,
	*	возможно, совместимость с функциями расчета разности между итераторами, типа distance.
	*/
    class integer_iterator
	{
        long pos_; // текущая позиция итератора

	public:

		// Типов указателя и ссылки нет, так как отсутствующий контейнер недоступен для записи
		typedef std::random_access_iterator_tag iterator_category;
        typedef long value_type;
        typedef long difference_type;
        typedef long *pointer;
        typedef long reference;

		/// Конструктор
        integer_iterator(long pos=0) : pos_(pos) {}

		/// Оператор разыменования
		value_type operator*() const { return pos_; }

		/// Оператор сдвига итератора
        integer_iterator &operator++ () { ++pos_; return *this; }

		/// Оператор сдвига итератора
        integer_iterator &operator-- () { --pos_; return *this; }

		/// Оператор сдвига итератора
        integer_iterator &operator+= (difference_type n) { pos_ += n; return *this; }

		/// Оператор сдвига итератора
        integer_iterator &operator-= (difference_type n) { pos_ -= n; return *this; }

        /// Оператор сдвига итератора
        reference operator[](difference_type n) { return pos_ + n; }

        /// Оператор сдвига итератора
        integer_iterator operator+(difference_type n) { return integer_iterator(pos_ + n); }

        /// Оператор сдвига итератора
        integer_iterator operator-(difference_type n) { return integer_iterator(pos_ - n); }

        /// Оператор сравнения для использования с операторами ++ и --
        bool operator!= (const integer_iterator &t) const { return pos_ != t.pos_; }

        /// Оператор сравнения для использования с операторами ++ и --
        bool operator== (const integer_iterator &t) const { return pos_ == t.pos_; }

        /**
		* \brief Оператор сравнения для использования с операторами += и -=
		* \note Использование оператора != при инкрементировании через +=, -= чревато
		*	"проскоком" заданной другим итератором верхней границы с неверным поведением далее.
		*/
        bool operator<  (const integer_iterator &t) const { return pos_ < t.pos_; }
	};

    inline long operator-(integer_iterator a, integer_iterator b) { return *b - *a; }

	/// расстояние между элементами
    inline integer_iterator::difference_type
    distance(integer_iterator start, integer_iterator end) { return *end - *start; }


    /**
    * \brief Итератор, проходящий последовательность вещественных чисел по заданому шагу.
    *
    *	Хранить бесконечную последовательность целых чисел бессмысленно. Так что такой
    *	контейнер является виртуальным: реализации его нет, но извлекать данные из него возможно.
    *	Для такого контейнера уровень совместимости с STL на уровне схожести интерфейса и,
    *	возможно, совместимость с функциями расчета разности между итераторами, типа distance.
    */
    class real_iterator
    {
    public:

        // Типов указателя и ссылки нет, так как отсутствующий контейнер недоступен для записи
        typedef std::random_access_iterator_tag iterator_category;
        typedef double value_type;
        typedef double difference_type;
        typedef double *pointer;
        typedef double reference;

        /// конструктор
        real_iterator(value_type pos=0, value_type step=1.) : pos_(pos), step(step) {}

        /// разыменование
        value_type operator*() const { return pos_; }

        /// сдвиг итератора
        real_iterator &operator++ () { pos_ = pos_ + step; return *this; }

        /// сравнение ("!=" - не дошел до end(), "==" - перешел через end())
        bool operator== (const real_iterator &t) const { return pos_ >= t.pos_; }
        bool operator!= (const real_iterator &t) const { return pos_ < t.pos_; }

    private:

        value_type pos_; // текущая позиция итератора
        value_type step; // шаг сдвига
    };

    /// расстояние между элементами
    inline real_iterator::difference_type distance(real_iterator start, real_iterator end)
        { return *end - *start; }


	/**
	* \brief Итератор, оперирующий с любым подмножеством заданного вектора в произвольном порядке.
	*
	*	Данный итератор является альтернативным итератором доступа к вектору. Для него могут
	*	быть полезны алгоритмы STL (упорядочение, поиск, трансформация), которые что-то правят
	*	в подмножестве вектора. Требуется образцовая совместимость с STL.
	*
	*	Порядок обхода задается индексным итератором Iterator, который возвращает порядок
	*	обработки массива (последовательность его индексов). Если индексным итератор возвращает
	*	индексы, отсутствующие в заданной массиве, то результат неопределен.
	*
	*	Итератор позволяет не только пройти массив или его подмножество в любом порядке,
	*	включая повторные обращения к одному и тому же элементу, но и также он может быть
	*	использован для сортировки любого подмножества в порядке задаваемом индексным
	*	итератором Iterator. Если индексным итератор дает последовательность в которой есть
	*	совпадающие элементы (индексы), то результат сортировки неопределен. STL никто не
	*	проектировал на обработку таких случаев, которые невозможны для обычного итератора.
	*
	*	Испозовать данный итератор нужно следующим образом:
	*	- задать массив A (std::vector<_T> или _T*), элементы которого нужно перебирать,
	*	- задать массив индексов I (любой длины), перечислить в нем индексы массива A
	*		в произвольном порядке (можно с дублями),
	*	- создать 2 итератора на массив индексов (start, to), которые выбирают последовательность
	*		индексов для перебора (std::vector<int>::iterator или int*),
	*	- создать два итератора _array_iterator инициализированные массивом A и
	*		итераторами start и to,
	*	- и, наконец, перебирать (или сортировать) элементы A в порядке, заданном
	*		индексным итератором .
	*
	* \note Не следует путать обычный порядок индексов, при котором массив A рассматривается
	*	в порядке A[0], A[1], ... и порядок задаваемый индексным итератором A[I[0]], A[I[1]], ...
	*	Результат сортировки выглядит необычно. Например: для A={a, b, c, d} и I={2, 1, 0, 3}
	*	упорядоченый массив есть A'={c, b, a, d}. Действительно, A'[I[0]] = A'[2] = a и т.д.

	* \tparam _T тип элемента заданного std::vector<>
	* \tparam _Iterator внешний итератор, дающий индексы массива
	*/
	template <typename _Iterator, typename _T>
	class _array_iterator
	{
		_Iterator it_; // индексный итератор
		_T *v_; // опорный массив

	public:

		typedef _Iterator iterator_type;
		typedef std::random_access_iterator_tag iterator_category;
		typedef _T value_type;
		typedef ptrdiff_t difference_type;
		typedef _T* pointer;
		typedef _T& reference;

		_array_iterator() : it_(), v_(0) {}

		/**
		* \brief Конструктор итератора.
		*
		*	Пример конструирования: _array_iterator<float, std::vector<int>::iterator>
		*	it(&a[0], i.begin());
		*
		* \param a указатель на элемент массива, чьи элементы будут возвращаться
		* \param it итератор по индексам опорного массива
		*/
		_array_iterator(const _Iterator &it, _T *v) : it_(it), v_(v) {}

		/// "Скрытый" итератор для выполнения реальной работы.
		iterator_type base() const { return it_; }

		_array_iterator &operator++() { ++it_; return *this; }
		_array_iterator &operator--() { --it_; return *this; }
		_array_iterator &operator+= (difference_type n) { it_ += n; return *this; }
		_array_iterator &operator-= (difference_type n) { it_ -= n; return *this; }
		_array_iterator operator+ (difference_type n) const { return _array_iterator(it_ + n, v_); }
		_array_iterator operator- (difference_type n) const { return _array_iterator(it_ - n, v_); }

		reference operator* () const { return v_[*it_]; }
		pointer   operator->() const { return &v_[*it_]; }

		bool operator!=(const _array_iterator &it) const { return it_ != it.it_; }
	};

	template <typename _Iterator, typename _T>
	inline bool
	operator==(const _array_iterator<_Iterator, _T> &__x,
			   const _array_iterator<_Iterator, _T> &__y)
	{ return __x.base() == __y.base(); }

	template <typename _Iterator, typename _T>
	inline bool
	operator< (const _array_iterator<_Iterator, _T> &__x,
			   const _array_iterator<_Iterator, _T> &__y)
	{ return __y.base() < __x.base(); }

	template <typename _Iterator, typename _T>
	inline bool
	operator!=(const _array_iterator<_Iterator, _T> &__x,
			   const _array_iterator<_Iterator, _T> &__y)
	{ return !(__x == __y); }

	template <typename _Iterator, typename _T>
	inline typename _array_iterator<_Iterator, _T>::difference_type
	operator- (const _array_iterator<_Iterator, _T> &__y,
			   const _array_iterator<_Iterator, _T> &__x)
	{ return __y.base() - __x.base(); }


	/**
	* \brief Функтор, поддерживающий сортировку индексов массива вместо сортировки его элементов.
	*
	*	Для задач, в которых массив содержит "тяжелые" для копирования элементы, задачу сортировки
	*	можно решить через создание дополнительного массива индексов, которые указывают на
	*	элементы сортируемого массива. Сортировка в этом случае выполняется над индексами и
	*	для получения правильного порядка элементов следует обращаться к ним через промежуточный
	*	массив индексов. Чтобы выполнить индексную сортировку массива A[] нужно сделать следующее:
	*	- создать массив индексов I[] тех элементов A, которые нужно сортировать
	*	- вызвать std::sort() по массиву индексов std::sort(I, I + n, _index_less<_T>(A))
	* \note Функтор является дополнением к _array_iterator, но позволяет сортировать не элементы
	*	массива, а его индексы.
	*/
	template <typename _T>
	struct _index_less
	{
		const _T *a_; ///< адрес начального элемента массива для сортировки

		/// Инициализация начальным элементом массива для сортировки
		_index_less(const _T *a) : a_(a) {}

		/**
		* \brief Косвенное сравнение элементов.
		* \param i,j индексы массива
		*/
		bool operator()(size_t i, size_t j) const { return a_[i] < a_[j]; }
	};

	/**
	* \brief Итератор, проходящий любой контейнер, содержащий объекты структур, и дающий ссылку
	*	на член-данное объекта. Порядок обхода задается базовым итератором Iterator.
	*
	*	Очевидно, что данный итератор не должен быть полностью совместим с STL. Например,
	*	недопустима сортировка элементов. При сортировке значения член-данных будут "оторваны"
	*	от объектов, в которых они находятся и перетасованы между объектами. Вряд ли это имеет
	*	смысл, поскольку нарушает принцип сохранения корректного состояния объекта. Хотя
	*	с другой стороны, если дан прямой доступ к полю структуры, значит его модификация не
	*	разрушает корректность состояния объекта. Чтобы обеспечить неполную совместимость,
	*	достаточно не определять макрос MEMBER_ITERATOR_WITH_FULL_STL_COMPATIBILITY.
	*
	* \tparam _Iterator базовый итератор по объектам класса
	* \tparam _T тип объекта класса
	* \tparam _R тип член-данного объекта класса
	*/
	template <typename _Iterator, typename _T, typename _R>
	class _member_iterator
	{
		_Iterator it_;
		_R _T::* mp_;

	public:

		typedef _Iterator iterator_type;
		typedef typename std::iterator_traits<_Iterator>::iterator_category iterator_category;
		typedef _R value_type;
		typedef typename std::iterator_traits<_Iterator>::difference_type difference_type;
		typedef _R* pointer;
		typedef _R& reference;

		_member_iterator(_Iterator from, _R _T::* memptr) : it_(from), mp_(memptr) {}
		iterator_type base() const { return it_; }

		reference operator*() const { return (*it_).*mp_; }
		pointer operator->() const { return &((*it_).*mp_); }
		reference operator[](difference_type n) const { return it_[n].*mp_; }

		_member_iterator &operator++ () { ++it_; return *this; }
		_member_iterator &operator-- () { --it_; return *this; }
		_member_iterator &operator+= (difference_type n) { it_ += n; return *this; }
		_member_iterator &operator-= (difference_type n) { it_ -= n; return *this; }
		_member_iterator operator+   (difference_type n) const { return _member_iterator(it_ + n, mp_); }
		_member_iterator operator-   (difference_type n) const { return _member_iterator(it_ - n, mp_); }

		bool operator== (const _member_iterator &mi) const { return it_ == mi.it_; }
		bool operator!= (const _member_iterator &mi) const { return it_ != mi.it_; }
	};

#if MEMBER_ITERATOR_IS_FULL_STL_COMPATIBILE == 1

	template<typename _Iterator, typename _T, typename _R>
	inline bool
	operator< (const _member_iterator<_Iterator, _T, _R> &__x,
			   const _member_iterator<_Iterator, _T, _R> &__y)
	{ return __y.base() < __x.base(); }

	template<typename _Iterator, typename _T, typename _R>
	inline typename _member_iterator<_Iterator, _T, _R>::difference_type
	operator- (const _member_iterator<_Iterator, _T, _R> &__y,
			   const _member_iterator<_Iterator, _T, _R> &__x)
	{ return __y.base() - __x.base(); }

#endif

	/**
	* \brief Итератор по массиву указателей на объекты класса.
	*
	*	Разыменование указателя дает не содержимое ячейки массива (указатель на объект),
	*	а ссылку на сам объект. Итератор полностью совместим с STL.
	*
	* \note Класс имеет защищенные функции __next__ и __prev__. Данные функции перегружены
	*	по последнему параметру, имеющему категорию iterator_category, и определены только
	*	для типа random_access_iterator_tag из этой категории. Таким образом реализуется
	*	ограничение, которое не позволяет создавать функции-операторы +=, -=, +, - для
	*	итераторов других категорий. Произойдет ошибка компиляции error: no match for ‘operator+=’
	*	in ‘it += n’.
	*
	* \tparam _Iterator итератор любого STL контейнера, содержащего указатели, а не объекты
	*/
	template <typename _Iterator>
	struct _pointer_iterator : public _Iterator
	{
	public:

		typedef _Iterator iterator_type;
		typedef typename std::iterator_traits<_Iterator>::iterator_category iterator_category;
		typedef typename std::remove_pointer<typename std::iterator_traits<_Iterator>::value_type>::type value_type;
		typedef typename std::iterator_traits<_Iterator>::difference_type difference_type;
		typedef value_type* pointer;
		typedef value_type& reference;

		using iterator_type::operator++;
		using iterator_type::operator--;

		_pointer_iterator(const iterator_type &it) : iterator_type(it) {}

		iterator_type base() const { return (iterator_type) *this; }
		iterator_type &base() { return (iterator_type&) *this; }

		reference operator* () const { return ** (iterator_type)*this; }
		pointer   operator->() const { return * (iterator_type)*this; }

		_pointer_iterator &operator+= (difference_type n)
		{ base() = __advance__(base(), n, iterator_category()); return *this; }

		_pointer_iterator &operator-= (difference_type n)
		{ base() = __advance__(base(), -n, iterator_category()); return *this; }

		_pointer_iterator operator+ (difference_type n)
		{ return _pointer_iterator(__advance__(base(), n, iterator_category())); }

		_pointer_iterator operator- (difference_type n)
		{ return _pointer_iterator(__advance__(base(), -n, iterator_category())); }

	protected:

		iterator_type __advance__(iterator_type it, difference_type n, std::random_access_iterator_tag)
		{ it += n; return it; }

	};

	template <typename _Iterator>
	inline bool
	operator==(const _pointer_iterator<_Iterator> &__x,
			   const _pointer_iterator<_Iterator> &__y)
	{ return __x.base() == __y.base(); }

	template <typename _Iterator>
	inline bool
	operator< (const _pointer_iterator<_Iterator> &__x,
			   const _pointer_iterator<_Iterator> &__y)
	{ return __y.base() < __x.base(); }

	template <typename _Iterator>
	inline bool
	operator!=(const _pointer_iterator<_Iterator> &__x,
			   const _pointer_iterator<_Iterator> &__y)
	{ return !(__x == __y); }

	template <typename _Iterator>
	inline typename _pointer_iterator<_Iterator>::difference_type
	operator- (const _pointer_iterator<_Iterator> &__y,
			   const _pointer_iterator<_Iterator> &__x)
	{ return __y.base() - __x.base(); }


	// предварительное объявление для _piecewise_container
	template <typename> class _piecewise_iterator;


    /**
     * @brief область итерирования
     * @note область итерирования состоит из двух итераторов (начал и конца). Данный объект
     *  позволяет проходить область минуя чать объектов (фильрование), причем не выходя при
     *  таком фильтровании за пределы области. Для обычного подхода при фильтровании возможен
     *  пропуск итератора конца области с дальнейшим неопределенным поведением. Другим плюсом
     *  такого объекта является экономия ресурсов, так как ссылки на внешние объекты храняться в
     *  единственном экземпляре. Ну и третий плюс - большая простота алгоритмов, чем для пары
     *  итераторов (см Александреску).
     */
    template <typename _Iterator> struct _iterator_range : public std::pair<_Iterator, _Iterator>
    {
        typedef std::pair<_Iterator, _Iterator> _base;
        typedef std::forward_iterator_tag iterator_category;
        typedef typename std::iterator_traits<_Iterator>::value_type value_type;
        typedef typename std::iterator_traits<_Iterator>::difference_type difference_type;
        typedef typename std::iterator_traits<_Iterator>::pointer pointer;
        typedef typename std::iterator_traits<_Iterator>::reference reference;

        using _base::first; ///< доступ к полным возможностям вложенного итератора
        using _base::second;

        /// инициализация области итерирования
        _iterator_range(_Iterator t1, _Iterator t2) : _base(t1, t2) {}

        /// сужение области итерирования на один (или несколько при фильтрации) объектов
        _iterator_range &operator--() { ++_base::first; return *this; }

        /// доступ к объекту области
        reference operator*() const { return (reference) *_base::first; }

        /// доступ к объекту области
        pointer   operator->() const { return &(*_base::first); }

        /// проверка на завершение итерирования (сужение невозможно)
        bool empty() const { return _base::first == _base::second; }
    };

	/**
	* \brief Виртуальный контейнер, состоящий из нескольких стандартных STL контейнеров одного типа.
	* \tparam _Iterator тип итератора, характеризующий тип используемого стандарного STL контейнера
	* \note	Контейнер может быть связан как на обычным, так и на константным итератором, например:
	*	_piecewise_container<std::vector<T>::const_iterator>. При связи с константным итератором
	*	нет возможности модифицировать элементы контейнера.
	*/
	template <typename _Iterator>
    class _piecewise_container : public std::vector<_iterator_range<_Iterator> >
	{
        typedef std::vector<_iterator_range<_Iterator> > _base;

    public:

        /// Область просмотра используемого подконтейнера.
        typedef _iterator_range<_Iterator> iterator_range_type;

        /// Внутренний итератор по подконтейнерам
        typedef typename _base::const_iterator _container_iterator;

		/// Итератор по контейнеру
        typedef _piecewise_iterator<_Iterator> iterator;

	public:

        using _base::empty;
        using _base::size;
        using _base::push_back;
        using _base::pop_back;

		/**
		* \brief Оператор вставки области просмотра подконтейнера.
		* \param range пара итераторов начала и конца области просмотра контейнера
		* \note Можно многократно вставлять одну и ту же область, или несколько областей из
		*	одного и того же контейнера или из разных контейнеров одного типа.
		*/
        void push_back(_Iterator it, _Iterator ite) { push_back(iterator_range_type(it, ite)); }

		/// Итератор на первый элемент контейнера.
        iterator begin() const
        {
            if (empty()) return iterator(_base::begin(), _Iterator());
            return iterator(_base::begin(), _base::begin()->first);
        }

		/// Итератор на на последний (недействительный) элемент.
        iterator end() const
        {
            if (empty()) return iterator(_base::begin(), _Iterator());
            return iterator(_base::end(), _base::end()->first);
                // ссылка на _base::end()->first связана с переходом it_ = cit_->first
                // в _piecewise_iterator::__advance__(). Это неправильно, но без этого
                // не работает
        }

        _iterator_range<iterator> iterator_range() const
        { return _iterator_range<iterator>(begin(), end()); }

    };

	/**
	* \brief Итератор элементов во нескольких связанных контейнерах (_piecewise_container).
	*
	*	В текущей реализации итератор ограничен категорией forward_iterator_tag, что позволяет
	*	создать достаточно простую реализацию без излишних внутренних членов и полностью
	*	покрывает текущие задачи проекта. Кроме того, это позволяет использовть в качестве
	*	подконтейнеров любые STL контейнеры.
	*
	*	Структура может быть объявлена как на обычные, так и на константные итераторы, например:
	*	_piecewise_iterator<std::vector<T>::const_iterator>. Во втором случае нет возможности
	*	менять элементы.
	*/
	template <typename _Iterator>
	class _piecewise_iterator
	{
		typedef _piecewise_container<_Iterator> _container;
		typedef typename _container::_container_iterator _container_iterator;

		friend class _piecewise_container<_Iterator>;

		_container_iterator cit_; // текущая позиция итератора (подконтейнер)
		_Iterator it_; // текущая позиция итератора (в подконтейнере)

		_piecewise_iterator(_container_iterator cit, _Iterator it) : cit_(cit), it_(it) {}

	public:

		typedef std::forward_iterator_tag iterator_category;
		typedef typename std::iterator_traits<_Iterator>::value_type value_type;
		typedef typename std::iterator_traits<_Iterator>::difference_type difference_type;
		typedef typename std::iterator_traits<_Iterator>::pointer pointer;
		typedef typename std::iterator_traits<_Iterator>::reference reference;

		_piecewise_iterator() {}

		reference operator* () const { return (reference) *it_; }
		pointer   operator->() const { return &(*it_); }
		_piecewise_iterator &operator++ () { __advance__(); return *this; }

		bool operator==(const _piecewise_iterator<_Iterator> &__y) const
		{ return cit_ == __y.cit_ && it_ == __y.it_; }

		bool operator!=(const _piecewise_iterator<_Iterator> &__y) const
		{ return cit_ != __y.cit_ || it_ != __y.it_; }

	protected:

		/// реализация операции инкремента
		void __advance__()
		{
			if ( ++it_ == cit_->second) { ++cit_ ; it_ = cit_->first; }
		}
	};

    /**
    * \brief Расширенный виртуальный контейнер, состоящий из ряда последовательностей элементов
    *   типа T, входящих в тот или иной вектор.
    * \note	Контейнер включает подпоследовательности только из векторов типа (T*), посколько только
    *   для них возможно обеспечить гарантии безопасности использования, путем установки T* = 0.
    *   Итерирование по контейнеру защищено от пустых областей тем, что пустые области,
    *   в который begin() == end() не вставляются вообще. Защита от разыменования *end() и пустого
    *   итератора сделана путем добавления пустого игнорируемого элемента в начало контейнера.
    *
    * \tparam T тип элементов контейнера
    * \tparam S тип дополнительного объекта (информация о вставленной последовательности)
    *
    * \note	Итератор позволяет доступ к дополнительной информации о последовательности через
    *   вызов функции operator().
    */
    template <typename T, typename __S=int>
    class _piecewise_vector
    : public std::vector< std::pair< std::pair<T*, T*>, __S> >
    {
        typedef std::pair<T*, T*> __Range;
        typedef std::pair<__Range, __S> __Value;
        typedef std::vector<__Value> __Cont;
        typedef __Cont _base;

        /// собственный итератор, дающий наиболее эффективный обход контейнера
        class __Iterator
        {
            typedef typename __Cont::const_iterator container_iterator;
            friend class _piecewise_vector;

            container_iterator cit_; // текущая позиция итератора (указывает на подконтейнер)
            T *p_; // текущая позиция итератора (указывает на элемент в подконтейнере)

            __Iterator(container_iterator cit, T *p=nullptr) : cit_(cit), p_(p) {}

            void __advance__()
            { if ( ++p_ == cit_->first.second) { --cit_; p_ = cit_->first.first; } }

        public:

            typedef std::forward_iterator_tag iterator_category;
            typedef typename std::iterator_traits<T*>::value_type value_type;
            typedef typename std::iterator_traits<T*>::difference_type difference_type;
            typedef typename std::iterator_traits<T*>::pointer pointer;
            typedef typename std::iterator_traits<T*>::reference reference;

            __Iterator() {}
            __S operator()() const { return cit_->second; }
            reference operator* () const { return (reference)*p_; }
            pointer   operator->() const { return &(*p_); }
            __Iterator &operator++ () { __advance__(); return *this; }

            bool operator!=(const __Iterator &x) const { return cit_ != x.cit_ || p_ != x.p_; }
            bool operator==(const __Iterator &x) const { return !(*this != x); }
        };

    public:

        /// Внутренний итератор по подконтейнерам
        typedef typename __Cont::const_iterator container_iterator;

        /// Итератор по контейнеру
        typedef __Iterator iterator;

        /// Область просмотра используемого подконтейнера.
        typedef _iterator_range<iterator> iterator_range_type;

    public:

        using _base::empty;
        using _base::size;
        using _base::push_back;
        using _base::pop_back;

        _piecewise_vector() { _base::push_back(__Value(__Range(nullptr, nullptr), __S())); }

        /**
        * \brief Оператор вставки области просмотра подконтейнера.
        * \param range пара итераторов начала и конца области просмотра контейнера
        * \note Можно многократно вставлять одну и ту же область, или несколько областей из
        *	одного и того же контейнера или из разных контейнеров одного типа. Но невозможно
        *   вставить пустую область.
        */
        void push_back(T *ps, T *pe, __S s=__S())
        { if (ps != pe) _base::push_back(__Value(__Range(ps, pe), s)); }

        /// Итератор на первый элемент контейнера.
        iterator begin() const
        {
            container_iterator it = --_base::end();
            return iterator(it, (it->first).first);
        }

        /// Итератор на на последний (недействительный) элемент.
        iterator end() const
        {
            container_iterator it = _base::begin();
            return iterator(it, (it->first).first);
        }

        iterator_range_type iterator_range() const
        { return iterator_range_type(begin(), end()); }

    };

    /// Итератор для вставки элементов в массив с одновременным его расширением
	template <class _container>
	class sequental_write_iterator
	{
	protected:
		_container* container;
	public:
		explicit sequental_write_iterator(_container &x) : container(&x) {}
		sequental_write_iterator<_container>& operator=(typename _container::const_reference value)
		{ container->push_back(value); return *this; }
	};

	/// Итератор для чтения элементов массива с автоматическим его сдвигом на следующий элемент
	template <class _container>
	class sequental_read_iterator
	{
	protected:
		typename _container::const_iterator iterator;
	public:
		explicit sequental_read_iterator(const _container &x) : iterator(x.begin()) {}
		typename _container::const_reference operator*() { return *iterator++; }
	};


}
#endif // ITERATOR_H
