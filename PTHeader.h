#pragma once

#ifndef DHEADER_H_
#define DHEADER_H_

#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <map>


//typedef std::map<std::string, cellDoubles> mapCellDoubles;
//typedef std::map<std::string, double> mapDoubles;



typedef Eigen::VectorXd vector;
typedef Eigen::Vector3d vector3d;

typedef std::vector<double> cellDoubles;
typedef std::vector<Eigen::MatrixXd> cellMatrixs;
typedef std::vector<Eigen::VectorXd> cellVectors;
typedef std::vector<Eigen::Vector3d> cellVectors3d;



typedef std::complex<double> complexDoubles;
typedef std::vector<std::complex<double>> cellComplexVectors;
typedef std::vector<std::vector<std::complex<double>>> cell2ComplexVectors;
typedef std::vector<std::vector<std::vector<std::complex<double>>>> cell3ComplexVectors;


cellVectors sumCellVectors(cellVectors& A, cellVectors& B);

cellVectors sumCellVectors(cellVectors& A, cellDoubles& a_1, cellDoubles& a_2);

cellDoubles sumCellDoubles(cellDoubles& A, cellDoubles& B);


#endif // !DHEADER_H_

