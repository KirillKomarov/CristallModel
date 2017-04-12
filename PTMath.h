#pragma once

#ifndef MDPMATH_H_
#define MDPMATH_H_

#include <iostream>
#include <vector>
#include <complex>  
#include <Eigen/Eigen>

namespace math {
	double kronecker(int i, int j);
	double pi();
	std::complex<double> i();
}

Eigen::VectorXd getOrtogonalVector(double x, double y);
Eigen::VectorXd getUnitVector(Eigen::VectorXd vec);
Eigen::VectorXd cross2D(Eigen::VectorXd& vec_0);
Eigen::VectorXd cross3D(Eigen::VectorXd& a, Eigen::VectorXd& b);
double approx(double, double, double);

#endif // !MDPMATH_H_


