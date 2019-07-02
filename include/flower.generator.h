#ifndef FLOWER_GENERATOR_H
#define FLOWER_GENERATOR_H

#include <tuple>
#include <memory>

#include "defs.h"
#include "math.defs.h"
#include "flower.dome.h"
#include "flower.ensemble.h"

/// модель побега
class Flowermodel
{
    // модель побега состоит из совокупности цветков с нумерацией Flowerstamp и
    // массива[3], в котором лежат дни начала каждой фазы (0 - BUD, 1 - FLOWER, 2 - OVARY)

    std::vector<Flowerstamp> flowers;
    std::vector<point<3>> phases;

    Unsigned num_flowers1; // полное число цветков 1-го порядка
    Unsigned num_flowers2; // полное число цветков 2-го порядка

public:

    typedef std::map<Flowerstamp, point<3>> base_type;
    typedef std::pair<Flowerstamp, point<3>> value_type;

    /// создать пустой объект (создание происходит путем добавления цветков insert_flower())
    Flowermodel() : num_flowers1(0), num_flowers2(0) {}

    /// вставка цветка с датами фаз и учет вставленных цветков
    void insert_flower(const Flowerstamp &fs, const point<3> &v);

    /// полное число цветков на обеих осях вместе
    Unsigned num_flowers() const
    {
        assert(num_flowers1 + num_flowers2 == flowers.size());
        return flowers.size();
    }

    /// полное число осей 2-го порядка (паракладиев)
    Unsigned num_axis2() const { return 0; }

    /// получить данные по цветению {начало, конец} оси
    curve_t flowering_select(unsigned axis, unsigned char flags) const;

};
typedef std::shared_ptr<Flowermodel> pFlowermodel;

class Flowergenerator
{
    // "скрытые" параметры являются моделью, построенной по учетным побегам, они включают:
    //  axis1 - кривую распускания 1 порядка (высота, день распускания)
    //  axis2 - [0] кривая начала цветения 2 порядка (высота, день распускания)
    //          [1] кривая завершения цветения 2 порядка (высота, день завершения) (для демонстрации)
    //  dens - плотность цветения (высота, число цветков 2 порядка)
    hidden_parameters_type params;

public:

    /*!
     * \brief загрузить параметры генерации насчитанные в ансамбле для заданных побегов
     * \param sprouts список загружаемых побегов
     */
    Flowergenerator(const Flowerensemble &ensemble, const std::vector<unsigned> &sprouts)
    : params(ensemble.get_hidden_parameters(sprouts)) {}

    /*!
     * \brief создать модельный побег с заданным числом цветков на нем
     * \param num_flowers1 число цветов в побеге на 1й оси
     * \param num_flowers2 ПРИБЛИЗИТЕЛЬНОЕ число цветов в побеге на 2й осях
     * \note установка отрицательного числа цветков означает генерацию по умалчиванию
     */
    pFlowermodel generate_sprout(unsigned sprout=1, int num_flowers1=-1, int num_flowers2=-1);
};

#endif // FLOWER_GENERATOR_H
