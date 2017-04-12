#pragma once

#ifndef MDPPHYSICALLIST_H_
#define MDPPHYSICALLIST_H_

#include "PTMath.h"
#include <string>

namespace phys {
	double epsilon();
	double e();
	double m_e();
	double k_B();
	double E_0();
	double c();
}

double sm();
double mm();
double um();
double nm();
double keV();

double sec();
double ms();
double us();
double ns();
double ps();
double fs();
double kg();

double getDimLength(std::string dim);
double getDimTime(std::string dim);
double getDimEnergy(std::string dim);
double getDimMass(std::string dim);


#endif // !MDPPHYSICALLIST_H_


