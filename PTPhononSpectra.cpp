#include "stdafx.h"
#include "PTPhononSpectra.h"


PhononSpectra::PhononSpectra(int _numPoints)
{
    getInputData();
	createLattice();
	numPoints = _numPoints;
}
PhononSpectra::~PhononSpectra()
{

}
void PhononSpectra::loadMDfile(int step)
{
	r.clear();
	v.clear();
	//std::cout << "Loading File, wait..." << std::endl;
	int s;
	std::string f;
	Eigen::VectorXd r_(2); r_ << 0.0, 0.0;
	Eigen::VectorXd v_(2); v_ << 0.0, 0.0;
	std::string fullAdress_1 = fileName + "/" +"dump.lammpstrj." + std::to_string(step);
	std::fstream F;
	F.open(fullAdress_1);
	if (F)
	{
		s = 0;
		F >> f; F >> f;
		F >> f;
		F >> f; F >> f; F >> f; F >> f;
		F >> f;
		F >> f; F >> f; F >> f; F >> f; F >> f; F >> f;
        F >> f;
        F >> f;
        F >> f;
        F >> f;
        F >> f; F >> f;
		F >> f; F >> f; F >> f; F >> f; F >> f; F >> f;	F >> f;	F >> f;	// ITEM: ATOMS id x y  !!z
		while (!F.eof())
		{
			F >> f;
			F >> f;
			F >> r_(0);
			F >> r_(1);
			F >> v_(0);
			F >> v_(1);
			r.push_back(r_);
			v.push_back(v_);
			//std::cout << s << std::endl;
			s++;
		}
		F.close();
	}
	//std::cout << "File is downloading" << std::endl;
}
void PhononSpectra::getFrequencyAxis()
{
	std::vector<double> freq;
	for (int i = 0; i < numTimeStep; i++)
	{
		freq.push_back(2 * pi() / (numTimeStep*timeStep)*i);
	}
	recordVector(freq, "", "freq", "txt");
}
cellVectors PhononSpectra::computerSpectraBvK(cellVectors& k)
{
    lambda_1.clear();
    lambda_2.clear();
    cellVectors psBvK;
    vector omega(2);
    for (int i = 0; i < k.size(); i++)
    {
        createDynamicalMatrix2D(k[i], getNearestNeighbors2D());
        setEigenValues();
        omega(0) = sqrt(lambda_1.back()/massParticles);
        omega(1) = sqrt(lambda_2.back()/massParticles);
        psBvK.push_back(omega);
    }
    return psBvK;
}
cellDoubles PhononSpectra::computerSpectraBvK(cellVectors& k, int moda)
{
    lambda_1.clear();
    lambda_2.clear();
    cellDoubles psBvK;
	for (int i = 0; i < k.size(); i++)
	{
		createDynamicalMatrix2D(k[i], getNearestNeighbors2D());
        setEigenValues();
        if (moda == 0)   {
            psBvK.push_back(sqrt(lambda_1.back()/massParticles));
        }
        if (moda == 1)   {
            psBvK.push_back(sqrt(lambda_2.back()/massParticles));
        }
    }
    return psBvK;

}
cellDoubles PhononSpectra::analiticSpectraBvK(cellVectors& k, int mod)
{
	double delta_1, delta_2, s_a, s_b, s_c;
	cellVectors nn;
	cellDoubles phsp;
	for (int i = 0; i < k.size(); i++)
	{
		nn = getNearestNeighbors2D();
		R = nn[0];
		delta_1 = getSimmetricPartDirivate(+1);
		delta_2 = getSimmetricPartDirivate(-1);

		s_a = pow(sin(k[i].dot(nn[0]) / 2), 2);
		s_b = pow(sin(k[i].dot(nn[1]) / 2), 2);
		s_c = pow(sin(k[i].dot(nn[2]) / 2), 2);
		phsp.push_back(sqrt(2 * delta_1*(s_a + s_b + s_c) + mod * 2 * delta_2*sqrt(s_a*s_a + s_b*s_b + s_c*s_c - s_a*s_b - s_b*s_c - s_c*s_a)));
	}
	return phsp;
}
