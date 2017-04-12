#include "stdafx.h"
#include "PTMath.h"

namespace math {
	double kronecker(int i, int j)
	{
		if (i == j)
		{
			return 1;
		}
		else
		{
			return 0.0;
		}
	}
	double pi()
	{
		return 3.1415926535;
	}
	std::complex<double> i()
	{
		std::complex<double> im(0, 1);
		return im;
	}
}

Eigen::VectorXd getOrtogonalVector(double x, double y)
{
	Eigen::VectorXd e_n(2);
	Eigen::VectorXd e_k(2); e_k << x, y;

	double norm = std::hypot(x, y);
	if (norm != 0)
	{
		e_k = e_k / norm;
	}
	else
	{
		e_k = e_k;
	}
	e_n(0) = -e_k(1), e_n(1) = e_k(0);
	return e_n;
}
Eigen::VectorXd getUnitVector(Eigen::VectorXd vec)
{
	double norm = vec.norm();
	if (norm != 0)
	{
		return vec / norm;
	}
	else
	{
		return vec;
	}
}
Eigen::VectorXd cross2D(Eigen::VectorXd& vec_0)
{
	Eigen::VectorXd vec(2);
	vec(0) = -vec_0(1);
	vec(1) = vec_0(0);
	return vec;
}
Eigen::VectorXd cross3D(Eigen::VectorXd& a, Eigen::VectorXd& b)
{
	Eigen::VectorXd vec(3);
	vec(0) = a(2)*b(1) - a(1)*b(2);
	vec(1) = a(0)*b(2) - a(2)*b(0);
	vec(2) = a(1)*b(0) - a(0)*b(1);
	return vec;
}

double approx(double x, double y, double e)
{
	if (abs((x - y)) < e) {
		return 0.0;
	}
	else { 
		return 1.0;
	}
}
