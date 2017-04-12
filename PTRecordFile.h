#pragma once

#ifndef MYALGEBRA_H_
#define MYALGEBRA_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>  
#include <Eigen/Eigen>

//#include "ProgramConstant.h"

void recordVector(Eigen::VectorXd& vector, std::string adress, std::string name, std::string format);
void recordCollection(std::vector<Eigen::Vector3d>& collection, std::string adress, std::string name, std::string format, int s);
void recordCollection(std::vector<Eigen::VectorXd>& collection, std::string adress, std::string name, std::string format);
void recordCollection(std::vector<Eigen::Vector3d>& collection, std::string adress, std::string name, std::string format);
void recordVector(Eigen::VectorXd& vector, std::string adress, std::string name, std::string format, int s);
void recordMatrix(Eigen::MatrixXd& matrix, std::string adress, std::string name, std::string format);
void recordMatrix(Eigen::MatrixXd& matrix, std::string adress, std::string name, std::string format, int s);
void recordVector(std::vector<double> vector, std::string adress, std::string name, std::string format);


#endif


