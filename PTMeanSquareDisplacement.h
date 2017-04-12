#pragma once

#include "PTLattice.h"
#include "PTPotentials.h"
#include "PTRecordFile.h"
#include "PTHeader.h"
using std::cout;
using std::endl;


/*!
    \brief Унаследованный класс
    \author Комаров К.А.
    \version 2.0
    \date Апрель 2017 года

    Класс, для вычисления среднеквадратичного отклонения атомов в заданой системе.
    Здесь же вычисляются ангармоническая поправка для частоты \f$\Delta_k\f$ и обратное время жизни фотона \f$\Gamma_k\f$.

    Виртуально унаследован от класса Lattice и DynamicalMatrix, где хранятся все переменные, задающие кристалл.

    Цели класса: расчет MSD, ангармонических поправок и MSD с учетом ангармонизма.

*/
class MeanSquareDisplacement : public DynamicalMatrix, public Lattice {
    double msd;
	double omega_k, delta_k;
    double shift_1;
    complexDoubles shift_2;
    double smearing; ///< Параметр для расчета суммы в смысле главного значения.
public:
	MeanSquareDisplacement();
	~MeanSquareDisplacement();
    /*!
    Функция расчитывает ангармоническую поправку для частоты, а имеено величину \f$\Delta_k(\omega)\f$.
    \param [out] moda_0 - номер акустической ветки для которой вычисляется поправка.
    \param [out] k - коллекция обратных векторов вдоль которой необходимо вычислить поправку.
    \param [out] omega - фононный спектор (для всех веток) вдоль направления k, для которого ищется поправка по частоте.
    \param [out] omega_BZ - фононный спектр по первой зоне Бриллюэна для всех веток.
    \return [out] возвращает коллекцию поправок по частотам, \f$\Delta_k\f$, для заданой ветки.
    */
    cellDoubles getFrequencyShift(int moda_0, cellVectors& k, cellVectors& omega, cellVectors& omega_BZ);

    /*!
    This function calculate life time of the phonon by the way of reciprocal vector k, \f$\Delta_k(\omega)\f$.
    \param [out] moda_0 - number of the acoustic wave mode.
    \param [out] k - collection of reciprocal vector for calculation.
    \param [out] omega - фононный спектор (для всех веток) вдоль направления k, для которого ищется поправка по частоте.
    \param [out] omega_BZ - фононный спектр по первой зоне Бриллюэна для всех веток.
    \return [out] возвращает коллекцию поправок по частотам, \f$\Delta_k\f$, для заданой ветки.
    */
    cellDoubles getFrequencyAttenuation(cellVectors&, cellDoubles&, cellDoubles&);

    /*!
    This function MSD in BvK approch. Sum go to the first Brillouin zone (fBZ).
    \param [out] omega - phonon frequency on fBZ.
    \return [out] output the MSD for this cristall.
    */
    double getMSD(cellVectors& omega);

    /*!
    This function MSD in BvK approch. Sum go to the first Brillouin zone (fBZ).
    \param [out] omega - phonon frequency on fBZ.
    \return [out] output the MSD for this cristall.
    */
    double getRelativeMSD(cellVectors& k, cellVectors& omega);

    /*!
    This function MSD in BvK approch. Sum go to the first Brillouin zone (fBZ).
    \param [out] omega - phonon frequency by fBZ.
    \return [out] output the MSD for this cristall.
    */
    double getMSDwithAnharmonicCorrection(cellVectors&, cellVectors&);

    /*!
    This function MSD in BvK approch. Sum go to the first Brillouin zone (fBZ).
    \param [out] omega - phonon frequency by fBZ.
    \return [out] output the MSD for this cristall.
    */
    double getRelativeMSDwithAnharmonicCorrection(cellVectors&, cellVectors&);


    cellDoubles computerMSD_AE(cellVectors&, cellDoubles&, cellDoubles&);
	cellDoubles computerMSD(cellDoubles&, cellDoubles&);
};
