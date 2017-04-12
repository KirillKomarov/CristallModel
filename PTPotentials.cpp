#include "stdafx.h"
#include "PTPotentials.h"

using namespace math;


double LennardJones::getPotential(double r)
{
    return 4*energyScale*(pow(lengthScale / r, 12) - pow(lengthScale / r, 6));
}
double LennardJones::get1DirivatePotential(double r)
{
	//std::cout << -energyScale * 1 / lengthScale*(12 * pow(lengthScale / r, 13) - 6 * 2 * pow(lengthScale / r, 7)) << std::endl;
	//system("pause");
    return -4*energyScale * 1 / lengthScale*(12 * pow(lengthScale / r, 13) - 6 * pow(lengthScale / r, 7));
}
double LennardJones::get2DirivatePotential(double r)
{
    return 4*energyScale * 1 / (lengthScale*lengthScale)*(12 * 13 * pow(lengthScale / r, 14) - 6 * 7 * pow(lengthScale / r, 8));
}
double LennardJones::get3DirivatePotential(double r)
{
    return -4*energyScale * 1 / (lengthScale*lengthScale*lengthScale)*(12 * 13 * 14 * pow(lengthScale / r, 15) - 6 * 7 * 8 * pow(lengthScale / r, 9));
}
double LennardJones::get4DirivatePotential(double r)
{
    return 4*energyScale * 1 / (lengthScale*lengthScale*lengthScale*lengthScale)*(12 * 13 * 14 * 15 * pow(lengthScale / r, 16) - 6 * 7 * 8 * 9 * pow(lengthScale / r, 10));
}


double LennardJones2::getPotential(double r)
{
    return energyScale*(pow(lengthScale / r, 12) - 2*pow(lengthScale / r, 6));
}
double LennardJones2::get1DirivatePotential(double r)
{
    return -energyScale * 1 / lengthScale*(12 * pow(lengthScale / r, 13) - 2*6 * pow(lengthScale / r, 7));
}
double LennardJones2::get2DirivatePotential(double r)
{
    return energyScale * 1 / (lengthScale*lengthScale)*(12 * 13 * pow(lengthScale / r, 14) - 2*6 * 7 * pow(lengthScale / r, 8));
}
double LennardJones2::get3DirivatePotential(double r)
{
    return -energyScale * 1 / (lengthScale*lengthScale*lengthScale)*(12 * 13 * 14 * pow(lengthScale / r, 15) - 2*6 * 7 * 8 * pow(lengthScale / r, 9));
}
double LennardJones2::get4DirivatePotential(double r)
{
    return energyScale * 1 / (lengthScale*lengthScale*lengthScale*lengthScale)*(12 * 13 * 14 * 15 * pow(lengthScale / r, 16) - 2*6 * 7 * 8 * 9 * pow(lengthScale / r, 10));
}



double WCA50_49::getPotential(double r)
{
	return 50 * pow(50.0 / 49.0, 49)*energyScale*(pow(lengthScale / r, 50) - pow(lengthScale / r, 49));
}
double WCA50_49::get1DirivatePotential(double r)
{
	return -50 * energyScale * 1 / lengthScale*(50 * pow(lengthScale / r, 51) - 49 * pow(lengthScale / r, 50));
}
double WCA50_49::get2DirivatePotential(double r)
{
	return 50 * energyScale * 1 / (lengthScale*lengthScale)*(50 * 51 * pow(lengthScale / r, 52) - 49 * 50 * pow(lengthScale / r, 51));
}
double WCA50_49::get3DirivatePotential(double r)
{
	return 50 * pow(50 / 49, 49)*energyScale * 1 / (lengthScale*lengthScale)*(50 * 51 * pow(lengthScale / r, 52) - 49 * 50 * pow(lengthScale / r, 51));
}
double WCA50_49::get4DirivatePotential(double r)
{
	return 50 * pow(50 / 49, 49)*energyScale * 1 / (lengthScale*lengthScale)*(50 * 51 * pow(lengthScale / r, 52) - 49 * 50 * pow(lengthScale / r, 51));
}


double IPL9::getPotential(double r)
{
    return energyScale*pow(lengthScale / r, 9);
}
double IPL9::get1DirivatePotential(double r)
{
    return -energyScale * 9 * pow(lengthScale / r, 10)/lengthScale;
}
double IPL9::get2DirivatePotential(double r)
{
    return energyScale * 9 * 10 * pow(lengthScale / r, 11)/lengthScale/lengthScale;
}
double IPL9::get3DirivatePotential(double r)
{
    return -energyScale * 9*10*11 * pow(lengthScale / r, 12)/lengthScale/lengthScale/lengthScale;
}
double IPL9::get4DirivatePotential(double r)
{
    return energyScale * 9*10*11*12 * pow(lengthScale / r, 13)/lengthScale/lengthScale/lengthScale/lengthScale;
}



double IPL12::getPotential(double r)
{
    return energyScale*pow(1.0 / r, 12);
}
double IPL12::get1DirivatePotential(double r)
{
    return -energyScale * 12 * pow(1.0 / r, 13);
}
double IPL12::get2DirivatePotential(double r)
{
    return energyScale * 12 * 13 * pow(1.0 / r, 14);
}
double IPL12::get3DirivatePotential(double r)
{
    return -energyScale * 12*13*14 * pow(1.0 / r, 15);
}
double IPL12::get4DirivatePotential(double r)
{
    return energyScale * 12 * 13*14*15 * pow(1.0 / r, 16);
}



double Yukawa::getPotential(double r)
{
	return energyScale*lengthScale*exp(-r / lengthScale) / r;
}
double Yukawa::get1DirivatePotential(double r)
{
	return  -energyScale*lengthScale*exp(-r / lengthScale) / pow(r, 2)*(r / lengthScale + 1.0);
}
double Yukawa::get2DirivatePotential(double r)
{
	return energyScale*lengthScale*exp(-r / lengthScale) / pow(r, 3)*(2.0*(r / lengthScale + 1.0) + pow(r / lengthScale, 2));
}



double Hybrid::getPotential(double r)
{
    if (potentialType == "LJ")
    {
        return 4*energyScale*(pow(lengthScale / r, 12) - pow(lengthScale / r, 6));
    }
    if (potentialType == "IPL")
    {
        return energyScale*pow(lengthScale / r, numIPL);
    }
    if (potentialType == "Yukawa")
    {
        return energyScale*lengthScale*exp(-r / lengthScale) / r;
    }
}
double Hybrid::get1DirivatePotential(double r)
{
    if (potentialType == "LJ")
    {
        return -4*energyScale * 1 / lengthScale*(12 * pow(lengthScale / r, 13) - 6 * pow(lengthScale / r, 7));
    }
    if (potentialType == "IPL")
    {
        return -energyScale * numIPL * pow(lengthScale / r, numIPL + 1)/lengthScale;
    }
    if (potentialType == "Yukawa")
    {
        return  -energyScale*lengthScale*exp(-r / lengthScale) / pow(r, 2)*(r / lengthScale + 1.0);
    }
}
double Hybrid::get2DirivatePotential(double r)
{
    if (potentialType == "LJ")
    {
        return 4*energyScale * 1 / (lengthScale*lengthScale)*(12 * 13 * pow(lengthScale / r, 14) - 6 * 7 * pow(lengthScale / r, 8));
    }
    if (potentialType == "IPL")
    {
        return energyScale * numIPL * (numIPL + 1) * pow(lengthScale / r, (numIPL + 2))/lengthScale/lengthScale;
    }
    if (potentialType == "Yukawa")
    {
        return energyScale*lengthScale*exp(-r / lengthScale) / pow(r, 3)*(2.0*(r / lengthScale + 1.0) + pow(r / lengthScale, 2));
    }
}
double Hybrid::get3DirivatePotential(double r)
{
    if (potentialType == "LJ")
    {
        return -4*energyScale * 1 / (lengthScale*lengthScale*lengthScale)*(12 * 13 * 14 * pow(lengthScale / r, 15) - 6 * 7 * 8 * pow(lengthScale / r, 9));
    }
    if (potentialType == "IPL")
    {
        return -energyScale * numIPL*(numIPL + 1)*(numIPL + 2) * pow(lengthScale / r, numIPL + 3)/lengthScale/lengthScale/lengthScale;
    }
    if (potentialType == "Yukawa")
    {

    }
}
double Hybrid::get4DirivatePotential(double r)
{
    if (potentialType == "LJ")
    {
        return 4*energyScale * 1 / (lengthScale*lengthScale*lengthScale*lengthScale)*(12 * 13 * 14 * 15 * pow(lengthScale / r, 16) - 6 * 7 * 8 * 9 * pow(lengthScale / r, 10));
    }
    if (potentialType == "IPL")
    {
        return energyScale * numIPL*(numIPL + 1)* (numIPL + 2)*(numIPL + 3) * pow(lengthScale / r, numIPL + 4)/lengthScale/lengthScale/lengthScale/lengthScale;
    }
    if (potentialType == "Yukawa")
    {

    }
}



double  Potential::getFull2DirivatePotential(int alpha, int beta)
{
	double r = R.norm();
    return (get2DirivatePotential(r) - get1DirivatePotential(r) / r)*(R(alpha)*R(beta) / (r*r)) + get1DirivatePotential(r) / r*math::kronecker(alpha, beta);
}
double  Potential::getFull3DirivatePotential(int alpha, int beta, int gamma)
{
	double r = R.norm();
	/*
	f_1 = (get2DirivatePotential(r) - get1DirivatePotential(r) / r);
	df_1 = (get3DirivatePotential(r) - get2DirivatePotential(r) / r + get1DirivatePotential(r) / (r*r))*R(gamma) / r;
	f_2 = (R(alpha)*R(beta) / (r*r));
	df_2 = (kronecker(alpha, gamma)*R(beta) + R(alpha)*kronecker(beta, gamma)) / (r*r) - R(alpha)*R(beta)*R(gamma) / (r*r*r);
	//f_3 = get1DirivatePotential(r) / r*math::kronecker(alpha, beta);
	df_3 = (get2DirivatePotential(r)*R(gamma) / (r*r) - get1DirivatePotential(r) *R(gamma) / (r*r*r))*math::kronecker(alpha, beta);
	return f_1*df_2 + df_1*f_2 + df_3;
	*/
	/*f_1 = get2DirivatePotential(r) - get1DirivatePotential(r) / r;
	f_2 = get3DirivatePotential(r) - 3/r*f_1;
	f_3 = R(alpha)*R(beta)*R(gamma) / (r*r*r);
	f_4 = (R(alpha)*kronecker(beta, gamma) + R(beta)*kronecker(alpha, gamma) + R(gamma)*kronecker(alpha, beta)) / (r*r);
	return  f_2 * f_3 + f_1 * f_4;*/
    f_1 = get2DirivatePotential(r) - get1DirivatePotential(r) / r;
	f_2 = R(alpha)*R(beta) / (r*r);
    f_3 = get1DirivatePotential(r)*kronecker(alpha, beta) / r;
    df_1 = (get3DirivatePotential(r) - get2DirivatePotential(r) / r + get1DirivatePotential(r) / (r*r))*R(gamma) / r;
	df_2 = (R(alpha)*kronecker(gamma, beta) + R(beta)*kronecker(alpha, gamma)) / (r*r) - 2 * R(alpha)*R(beta)*R(gamma) / (r*r*r*r);
    df_3 = (get2DirivatePotential(r) - get1DirivatePotential(r) / r)*R(gamma)*kronecker(alpha, beta) / (r*r);
	return df_1*f_2 + f_1*df_2 + df_3;
}
double  Potential::getFull4DirivatePotential(int alpha, int beta, int gamma, int sigma)
{
	double r = R.norm();

    f_1 = get2DirivatePotential(r) - get1DirivatePotential(r) / r;
	f_2 = R(alpha)*R(beta) / (r*r);
    f_3 = get1DirivatePotential(r)*kronecker(alpha, beta) / r;
    df_1 = ((get3DirivatePotential(r) - get2DirivatePotential(r)) / r + get1DirivatePotential(r) / (r*r))*R(gamma) / r;
	df_2 = (R(alpha)*kronecker(gamma, beta) + R(beta)*kronecker(alpha, gamma)) / (r*r) - 2 * R(alpha)*R(beta)*R(gamma) / (r*r*r*r);
    df_3 = (get2DirivatePotential(r) / r - get1DirivatePotential(r) / r)*R(gamma)*kronecker(alpha, beta) / (r*r);
    ddf_1 = ((get4DirivatePotential(r) - get3DirivatePotential(r)) / r + 2 * get2DirivatePotential(r) / (r*r) - 2 * get1DirivatePotential(r) / (r*r*r))*R(gamma)*R(sigma) / (r*r)
		- df_1*R(sigma) / (r*r)
        + ((get3DirivatePotential(r) - get2DirivatePotential(r)) / r + get1DirivatePotential(r) / (r*r))*kronecker(gamma, sigma) / r;
	ddf_2 = (kronecker(alpha, sigma)*kronecker(gamma, beta) + kronecker(beta, sigma)*kronecker(alpha, gamma)) / (r*r)
		- 2 * (R(alpha)*kronecker(gamma, beta)*R(sigma) + R(beta)*kronecker(alpha, gamma)*R(sigma)) / (r*r*r*r)
		+ 8 * R(alpha)*R(beta)*R(gamma)*R(sigma) / (r*r*r*r*r*r)
		- 2 * (kronecker(alpha, sigma)*R(beta)*R(gamma) + R(alpha)*kronecker(beta, sigma)*R(gamma) + R(alpha)*R(beta)*kronecker(gamma, sigma)) / (r*r*r*r);
    ddf_3 = (get3DirivatePotential(r) - get2DirivatePotential(r) / r + get1DirivatePotential(r) / (r*r))*R(gamma)*R(sigma)*kronecker(alpha, beta) / (r*r*r)
        + (get2DirivatePotential(r) - get1DirivatePotential(r) / r)*kronecker(alpha, beta) / (r*r)*kronecker(gamma, sigma)
		- 2 * df_3*R(sigma) / (r*r);


	return ddf_1*f_2 + 2 * df_1*df_2 + f_1*ddf_2 + ddf_3;
}
double  Potential::getSimmetricPartDirivate(int coef)
{
	double r = R.norm();
	//return (get2DirivatePotential(r) + coef*get1DirivatePotential(r) / r);
	//std::cout << (get2DirivatePotential(r) + coef*get1DirivatePotential(r) / r) << std::endl;
	//system("pause");
	return 1.0;
}



double DynamicalMatrix::getEigenValues(int moda)
{
    double lambda;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(D);
    lambda = es.eigenvalues()[moda].real();
    return lambda;
}
vector DynamicalMatrix::getEigenVecrtors(int moda)
{
    vector e;
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(D);
    e = es.eigenvectors().col(moda).real();
    return e;
}
void DynamicalMatrix::setEigenVecrtors()
{
	// Finding of the eigen vectors: 
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(D);
	e_1.push_back(es.eigenvectors().col(0).real());
	e_2.push_back(es.eigenvectors().col(1).real());
}
void DynamicalMatrix::setEigenValues()
{
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(D);
    lambda_1.push_back(es.eigenvalues()[0].real());
    lambda_2.push_back(es.eigenvalues()[1].real());
}
void DynamicalMatrix::createDynamicalMatrix2D(vector& k, cellVectors& nn)
{
	D.resize(2, 2);
    for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			sum_P = 0;
			for (int s = 0; s < nn.size(); s++)
			{
				R = nn[s];
                sum_P += getFull2DirivatePotential(i, j)*(math::i().imag() - std::exp(math::i()*k.dot(R)));
                //sum_P += getFull2DirivatePotential(i, j)*(cos(k.dot(R)) - 1);
			}
			D(i, j) = sum_P;
		}
    }
}
void DynamicalMatrix::createDynamicalMatrix2D(cellVectors& k, cellVectors& nn)
{
    e_1.clear();
    e_2.clear();
    lambda_1.clear();
    lambda_2.clear();
    for (int s = 0; s < k.size(); s++)
	{
		createDynamicalMatrix2D(k[s], nn);
        setEigenValues();
        setEigenVecrtors();
	}
}
cell3ComplexVectors DynamicalMatrix::getPotencial_3(cellVectors& k, cellVectors& nn)
{
	createDynamicalMatrix2D(k, nn);
	cellComplexVectors F_1;
	cell2ComplexVectors F_2;
	cell3ComplexVectors F;
	sum_P = 0.0;
	for (int s_1 = 0; s_1 < k.size(); s_1++)
	{
		for (int s_2 = 0; s_2 < k.size(); s_2++)
		{
			for (int s_3 = 0; s_3 < k.size(); s_3++)
			{
				sum_P = 0;
				for (int j = 0; j < nn.size(); j++)
				{
					for (int alpha = 0; alpha < 2; alpha++)
					{
						for (int beta = 0; beta < 2; beta++)
						{
							for (int gamma = 0; gamma < 2; gamma++)
							{
								R = nn[j];
								sum_P += getFull3DirivatePotential(alpha, beta, gamma)*e_1[s_1](alpha)*e_1[s_2](beta)*e_1[s_3](gamma)
									*(math::i().imag() - std::exp(math::i()*k[s_1].dot(R)))*(math::i().imag() - std::exp(math::i()*k[s_2].dot(R)))*(math::i().imag() - std::exp(math::i()*k[s_3].dot(R)));
							}
						}
					}
				}
				system("pause");
				F_1.push_back((1.0 / 12.0 / sqrt(numParticles*pow(massParticles, 3)))*sum_P);

			}
			F_2.push_back(F_1);
		}
		F.push_back(F_2);
	}
	return F;
}
complexDoubles DynamicalMatrix::getPotencial_3(int moda_1, vector& k_1,
                                               int moda_2, vector& k_2,
                                               int moda_3, vector& k_3,    cellVectors& nn)
{
    vector e_k1, e_k2, e_k3;

	createDynamicalMatrix2D(k_1, nn);
    e_k1 = getEigenVecrtors(moda_1);
	createDynamicalMatrix2D(k_2, nn);
    e_k2 = getEigenVecrtors(moda_2);
	createDynamicalMatrix2D(k_3, nn);
    e_k3 = getEigenVecrtors(moda_3);
	complexDoubles F_1;
	sum_P = 0.0;
	int numNeigbor = nn.size();
	for (int j = 0; j < numNeigbor; j++)
	{
		for (int alpha = 0; alpha < 2; alpha++)
		{
			for (int beta = 0; beta < 2; beta++)
			{
				for (int gamma = 0; gamma < 2; gamma++)
				{
					R = nn[j];
                    sum_P += getFull3DirivatePotential(alpha, beta, gamma)*e_k1(alpha)*e_k2(beta)*e_k3(gamma)
                        *(math::i().imag() - std::exp(-math::i()*k_1.dot(R)))*(math::i().imag() - std::exp(-math::i()*k_2.dot(R)))*(math::i().imag() - std::exp(-math::i()*k_3.dot(R)));
				}
			}
		}
	}
    F_1 = sum_P;
	return F_1;
}
complexDoubles DynamicalMatrix::getPotencial_4(int moda_1, vector& k_1,
                                               int moda_2, vector& k_2,
                                               int moda_3, vector& k_3,
                                               int moda_4, vector& k_4,    cellVectors& nn)
{
    vector e_k1, e_k2, e_k3, e_k4;
    createDynamicalMatrix2D(k_1, nn);
    e_k1 = getEigenVecrtors(moda_1);
	createDynamicalMatrix2D(k_2, nn);
    e_k2 = getEigenVecrtors(moda_2);
	createDynamicalMatrix2D(k_3, nn);
    e_k3 = getEigenVecrtors(moda_3);
	createDynamicalMatrix2D(k_4, nn);
    e_k4 = getEigenVecrtors(moda_4);
	complexDoubles F_1;
	sum_P = 0.0;
	for (int j = 0; j < nn.size(); j++)
	{
		for (int alpha = 0; alpha < 2; alpha++)
		{
			for (int beta = 0; beta < 2; beta++)
			{
				for (int gamma = 0; gamma < 2; gamma++)
				{
					for (int sigma = 0; sigma < 2; sigma++)
					{
						R = nn[j];
                        sum_P += getFull4DirivatePotential(alpha, beta, gamma, sigma)*e_k1(alpha)*e_k2(beta)*e_k3(gamma)*e_k4(sigma)
                            *(math::i().imag() - std::exp(-math::i()*k_1.dot(R)))*(math::i().imag() - std::exp(-math::i()*k_2.dot(R)))*(math::i().imag() - std::exp(-math::i()*k_3.dot(R)))*(math::i().imag() - std::exp(math::i()*k_4.dot(R)));

					}
				}
			}
		}
	}
    F_1 = sum_P;
	return F_1;
}
