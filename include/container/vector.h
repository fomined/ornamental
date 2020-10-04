#ifndef PRG_VECTOR_H
#define PRG_VECTOR_H

#include <vector>
#include "defs.h"

namespace prg
{
    /*!
     * \brief вектор со смещенной нумерацией элементов
     * \note !ВНИМАНИЕ для std::vector нет необходимости узнавать значение верхнего и нижнего
     *  индексов массива, они известны заранее [0, size()), по этой причине широко используется
     *  следующий способ прохождения индексов как (Unsigned i=0; i<size(); i++).
     *  Вектор со смещенной нумерацией индексов показывает дефектность такой схемы, так как в нем,
     *  верхний индекс не равен size(), и такой цикл выходит за пределы индексов, ибо operator[]
     *  делает смещение. По этой причине правильным способом прохождения цикла является
     *      for (Signed i=v.lower_index(); i<v.upper_index(); i++)
     *      for (Signed i=v.minimal_index(); i<=v.maximal_index(); i++)
     */
    template <typename T>
    class vector : public std::vector<T>
    {
        typedef std::vector<T> Base;
        Signed start_pos;

    public:
        vector() : Base(0, T()), start_pos(0) {}
        vector(Signed pos, Unsigned count, T init=0) : Base(count, init), start_pos(pos) {}

        const T &operator[](Signed i) const { return Base::operator[](i - start_pos); }
        T &operator[](Signed i) { return Base::operator[](i - start_pos); }

        void resize(Signed pos, Unsigned count) { start_pos = pos; Base::resize(count, T()); }

        Signed lower_index() const { return start_pos; }
        Signed upper_index() const { return start_pos + Base::size(); }

        Signed minimal_index() const { return start_pos; }
        Signed &minimal_index() { return start_pos; }
        Signed maximal_index() const { return start_pos + Base::size() - 1; }
    };

}
#endif // PRG_VECTOR_H

