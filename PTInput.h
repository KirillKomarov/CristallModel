#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include "PTPhysicalList.h"

/*!
    \brief Базовый класс
    \author Комаров К.А.
    \version 2.0
    \date Апрель 2017 года
    \warning Следует наследовать только виртуально!

    Основной класс, являющийся прородителм всего дерева, за исключением базовых классов,
    отвечающих за константы(синглоты) и стандартные алгоритмы

    Данный класс(структура) необходим для считывания входных данных из файла с названием "GInputData.txt".
    Поэтому в папке сборки необходимо создать данный файл и записать все необходимые макро-команды.

    Цели класса: хранение и обработка входных данных(переменных)

    \param[in] timeStep Шаг по времени, при котором проводилость MD моделирование
    \param[in] numTimeStep Количество шагов по времени (известо из MD моделирования)
    \param[in] numParticles Количество частиц в системе (чаще определяется стуктурой и размерами решетки)
    \param[in] massParticles Масса частиц системы

    \param[in] energyScale безрамерный параметр потенциала отвечающий за номировку энергии
    \param[in] lengthScale безрамерный параметр потенциала отвечающий за номировку длины

*/


struct MDPInput {
	std::string fileName;
	double timeStep;
	int numTimeStep;
	int numParticles;
	double massParticles;

	double energyScale;
	double lengthScale;
	
	double Temp;
	std::string latticeDim;
	std::string latticeType;
    std::string cristallType;
    std::string potentialType;
    int numIPL;

    int cristallRow;
    int cristallLength;
    int cristallWidth;
	double translation;
	void getInputData();
};
