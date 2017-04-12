#pragma once
#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <Eigen/Eigen>

#include "PTMath.h"
#include "PTRecordFile.h"
#include "PTHeader.h"
#include "PTInput.h"

using  namespace math;
using std::cout;
using std::endl;


/*!
    \brief Псевдобазовый класс
    \author Комаров К.А.
    \version 2.0
    \date Апрель 2017 года

    Основной класс, создающий решетку системы частиц и обратную решетку волновых векторов.

    Виртуально унаследован от класса PTInput, где хранятся все переменные, задающие кристалл.

    Цели класса: залание типа решетки и создание кристалла и обратной решетки волновых векторов

*/

struct Lattice : virtual public MDPInput {
    double vol; ///< Указывает, что элемент недоступен для использования
    vector a_0, a_1, a_2, a_3; ///< Вектора, образующие элементарную ячейку, [vector].
    vector b_1, b_2, b_3; ///< Вектора, образующие обратную ячейку, [vector].
    cellVectors k_BZ; ///< коллекция векторов образующая первую зону Бриллюэна, [cellVector]
    cellVectors atoms; ///< коллекция векторов задающая кристалл, [cellVector]

    /*!
    Функция создающая кристалл
    \param [in] latticeType - тип кристаллической решетки.
    \return [in] atoms - коллекция векторов задающая кристалл.
    */
    void createLattice();

    /*!
    Функция создающая квадратную элементарную ячейку.
    \param [out] a - вектор трансляции.
    \return [out] a_1, a_2, a_3 - вектора, образующие элементарную ячейку.
    */
    void createSquareLattice2D(double a);


    void createSquareLattice3D(double);
    void createHexagonalLattice2D(double);
    void createHexagonalLattice3D(double);
    void setReciprocalVectors3D();
    void setReciprocalVectors2D();

	void setFirstBrillouinZone2D();

    void createHexagonalCristall2D();
    void createCristall2D(double, double);
    void createCristall3D(double, double, double);

    cellVectors getTranslationVectors();
    cellVectors getNearestNeighbors2D();
	cellVectors getFirstBrillouinZone2D();
    cellVectors getDirectionGK(int);
	cellVectors getDirectionMG(int);
	cellVectors getDirectionGKMG(int);
	bool kroneker_G(Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&);
	bool kroneker_G(Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&);
};


#endif /* LATTICE_H */
