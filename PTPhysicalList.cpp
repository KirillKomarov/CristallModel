#include "stdafx.h"
#include "PTPhysicalList.h"

namespace phys {
	double epsilon()
	{
		return 8.85418781762039e-12;
	}
	double e()
	{
		return 1.6e-19;
	}
	double m_e()
	{
		return 9.10938356e-31;
	}
	double k_B()
	{
		return 1.38e-23;
	}
	double E_0()
	{
		return 511;	//[keV], rest energy for electron
	}
	double c()
	{
		return 3e8;	// [m/sec], speed of light
	}
}


double sm()
{
	return 1e-2;
}
double mm()
{
	return 1e-3;
}
double um()
{
	return 1e-6;
}
double nm()
{
	return 1e-9;
}

double keV()
{
	return 1e-3;
}

double sec()
{
	return 1;
}
double ms()
{
	return 10e-4;
}
double us()
{
	return 10e-7;
}
double ns()
{
	return 10e-10;
}
double ps()
{
	return 10e-13;
}
double fs()
{
	return 10e-16;
}
double kg()
{
	return 1;
}

double getDimLength(std::string dim)
{
	if (dim == "(sm)") return sm();
	if (dim == "(mm)") return mm();
	if (dim == "(um)") return um();
	if (dim == "()") return 1.0;
}
double getDimTime(std::string dim)
{
	if (dim == "(sec)") return sec();
	if (dim == "(ms)") return ms();
	if (dim == "(us)") return us();
	if (dim == "()") return 1.0;
}
double getDimEnergy(std::string dim)
{
	if (dim == "(keV)") return keV();
}
double getDimMass(std::string dim)
{
	if (dim == "(kg)") return kg();
	if (dim == "()") return 1.0;
}
