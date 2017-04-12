#pragma once

#include "PTInput.h"
#include "PTMath.h"
#include "PTHeader.h"
#include <Eigen/Eigen>
using namespace math;


/**
    \brief Базовый класс, виртульно унаследованный от PTInput
    \author Комаров К.А.
    \version 2.0
    \date Апрель 2017 года
    \warning Следует наследовать только виртуально!

    Основной класс, который задает тип потенциала для всего дерева, а именно определяет зачение потенциала в заданой точке, а также его производные.

    Цель класса: задать тип потенцила Леннарда-Джонса
    \f{eqnarray*}{
       U(r) = 4\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right],
    \f}

    \param[in] timeStep Шаг по времени, при котором проводилость MD моделирование
    \param[in] numTimeStep Количество шагов по времени (известо из MD моделирования)
    \param[in] numParticles Количество частиц в системе (чаще определяется стуктурой и размерами решетки)
    \param[in] massParticles Масса частиц системы

    \param[in] energyScale безрамерный параметр потенциала отвечающий за номировку энергии
    \param[in] lengthScale безрамерный параметр потенциала отвечающий за номировку длины

*/
class LennardJones : virtual public MDPInput {
public:
	double getPotential(double);
	double get1DirivatePotential(double);
	double get2DirivatePotential(double);
	double get3DirivatePotential(double);
	double get4DirivatePotential(double);
};



/**
    \brief Базовый класс, виртульно унаследованный от PTInput
    \author Комаров К.А.
    \version 2.0
    \date Апрель 2017 года
    \warning Следует наследовать только виртуально!

    Основной класс, который задает тип потенциала для всего дерева, а именно определяет зачение потенциала в заданой точке, а также его производные.

    Цель класса: задать тип потенцила Леннарда-Джонса с смещенным минимумом
    \f{eqnarray*}{
       U(r) = \varepsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - 2\left(\frac{\sigma}{r}\right)^{6} \right],
    \f}

    \param[in] timeStep Шаг по времени, при котором проводилость MD моделирование
    \param[in] numTimeStep Количество шагов по времени (известо из MD моделирования)
    \param[in] numParticles Количество частиц в системе (чаще определяется стуктурой и размерами решетки)
    \param[in] massParticles Масса частиц системы

    \param[in] energyScale безрамерный параметр потенциала отвечающий за номировку энергии
    \param[in] lengthScale безрамерный параметр потенциала отвечающий за номировку длины

*/
class LennardJones2 : virtual public MDPInput {
public:
    double getPotential(double);
    double get1DirivatePotential(double);
    double get2DirivatePotential(double);
    double get3DirivatePotential(double);
    double get4DirivatePotential(double);
};



/**
    \brief Базовый класс, виртульно унаследованный от PTInput
    \author Комаров К.А.
    \version 2.0
    \date Апрель 2017 года
    \warning Следует наследовать только виртуально!

    Основной класс, который задает тип потенциала для всего дерева, а именно определяет зачение потенциала в заданой точке, а также его производные.

    Цель класса: задать тип потенцила WCA50/49
    \f{eqnarray*}{
       U(r) = \frac{50}{49}\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{50} - \left(\frac{\sigma}{r}\right)^{49} \right],
    \f}

    \param[in] timeStep Шаг по времени, при котором проводилость MD моделирование
    \param[in] numTimeStep Количество шагов по времени (известо из MD моделирования)
    \param[in] numParticles Количество частиц в системе (чаще определяется стуктурой и размерами решетки)
    \param[in] massParticles Масса частиц системы

    \param[in] energyScale безрамерный параметр потенциала отвечающий за номировку энергии
    \param[in] lengthScale безрамерный параметр потенциала отвечающий за номировку длины

*/
class WCA50_49 : virtual public MDPInput {
public:
	double getPotential(double);
	double get1DirivatePotential(double);
	double get2DirivatePotential(double);
	double get3DirivatePotential(double);
	double get4DirivatePotential(double);
};



/**
    \brief Базовый класс, виртульно унаследованный от PTInput
    \author Комаров К.А.
    \version 2.0
    \date Апрель 2017 года
    \warning Следует наследовать только виртуально!

    Основной класс, являющийся прородителм всего дерева, за исключением базовых классов,
    отвечающих за константы(синглоты) и стандартные алгоритмы

    Цель класса: задать тип потенцила Леннарда-Джонса
    \f{eqnarray*}{
       U(r) = 4\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right],
    \f}

    \param[in] timeStep Шаг по времени, при котором проводилость MD моделирование
    \param[in] numTimeStep Количество шагов по времени (известо из MD моделирования)
    \param[in] numParticles Количество частиц в системе (чаще определяется стуктурой и размерами решетки)
    \param[in] massParticles Масса частиц системы

    \param[in] energyScale безрамерный параметр потенциала отвечающий за номировку энергии
    \param[in] lengthScale безрамерный параметр потенциала отвечающий за номировку длины

*/
class IPL9 : virtual public MDPInput {
public:
    double getPotential(double);
    double get1DirivatePotential(double);
    double get2DirivatePotential(double);
    double get3DirivatePotential(double);
    double get4DirivatePotential(double);
};




class IPL12 : virtual public MDPInput {
public:
	double getPotential(double);
	double get1DirivatePotential(double);
	double get2DirivatePotential(double);
	double get3DirivatePotential(double);
	double get4DirivatePotential(double);
};


class Yukawa : virtual public MDPInput {
public:
	double getPotential(double);
	double get1DirivatePotential(double);
	double get2DirivatePotential(double);
};


/**
    \brief Псевдобазовый класс, виртульно унаследованный от PTInput
    \author Комаров К.А.
    \version 2.0
    \date Апрель 2017 года

    Основной класс, который задает тип потенциала для всего дерева, а именно определяет зачение потенциала в заданой точке, а также его производные.

    Цель класса: задать потенцил в зависимости от макро-команд, полученных из файла "GInputData.txt"

    Макро-команды:
    - LJ - получаем потенциал Леннарда-Джонса:
    \f{eqnarray*}{
       U(r) = 4\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right],
    \f}
    - IPL  n - получаем потенциал IPLn:
    \f{eqnarray*}{
       U(r) = \varepsilon \left(\frac{\sigma}{r}\right)^{n},
    \f}
    - WCA  n    m - получаем потенциал WCAn/m:
    \f{eqnarray*}{
       U(r) = \frac{n}{m}\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{n} - \left(\frac{\sigma}{r}\right)^{m} \right] + \varepsilon,
    \f}

    \param[in] timeStep Шаг по времени, при котором проводилость MD моделирование
    \param[in] numTimeStep Количество шагов по времени (известо из MD моделирования)
    \param[in] numParticles Количество частиц в системе (чаще определяется стуктурой и размерами решетки)
    \param[in] massParticles Масса частиц системы

    \param[in] energyScale безрамерный параметр потенциала отвечающий за номировку энергии
    \param[in] lengthScale безрамерный параметр потенциала отвечающий за номировку длины

*/
class Hybrid : virtual public MDPInput {
public:
    double getPotential(double);
    double get1DirivatePotential(double);
    double get2DirivatePotential(double);
    double get3DirivatePotential(double);
    double get4DirivatePotential(double);
};




class Potential : virtual public Hybrid {
	double f_1, f_2, f_3, f_4;
	double df_1, df_2, df_3, df_4;
    double ddf_1, ddf_2, ddf_3, ddf_4;
public:
    vector R;
	double getFull2DirivatePotential(int, int);
	double getFull3DirivatePotential(int, int, int);
	double getFull4DirivatePotential(int, int, int, int);
	double getSimmetricPartDirivate(int);
};

class DynamicalMatrix : public Potential {
	std::complex<double> sum_P;
public:
	int mod;
	cellVectors e_1;
	cellVectors e_2;
	cellVectors e_3;
	cellDoubles lambda_1;
	cellDoubles lambda_2;
	Eigen::MatrixXcd D;
	void createDynamicalMatrix2D(Eigen::VectorXd&, cellVectors&);
	void createDynamicalMatrix2D(cellVectors&, cellVectors&);
    void setEigenVecrtors();
    void setEigenValues();
    double getEigenValues(int moda);
    vector getEigenVecrtors(int moda);
	cell3ComplexVectors getPotencial_3(cellVectors&, cellVectors&);
    complexDoubles getPotencial_3(int, Eigen::VectorXd&, int, Eigen::VectorXd&, int, Eigen::VectorXd&, cellVectors&);
    complexDoubles getPotencial_4(int, Eigen::VectorXd&, int, Eigen::VectorXd&, int, Eigen::VectorXd&, int, Eigen::VectorXd&, cellVectors&);
};

