#include "stdafx.h"
#include "PTMeanSquareDisplacement.h"

MeanSquareDisplacement::MeanSquareDisplacement()
{
    smearing = 10;
	getInputData();
    createLattice();
    setFirstBrillouinZone2D();
}
MeanSquareDisplacement::~MeanSquareDisplacement()
{

}
cellDoubles MeanSquareDisplacement::getFrequencyShift(int moda_0, cellVectors& k, cellVectors& omega, cellVectors& omega_BZ)
{
    auto numWaveVectors = k.size();
    auto numWaveVectors_BZ = k_BZ.size();
	cellDoubles Delta_k;
    double F;
    complexDoubles F_2;
    double omega_0, omega_1, omega_2, omega_p;
    vector dk, k_neg_1, k_neg_2;
    auto f_lim = [](double om, double eps) { return om / (om*om + eps*eps); };
    for (int s = 0; s < numWaveVectors; s++)
	{
        omega_0 = omega[s](moda_0);
        if (omega_0 != 0)
		{
			shift_1 = 0;
			shift_2 = 0;
            for (int i = 0; i < numWaveVectors_BZ; i++)
			{
                for (int eta = 0; eta < 2; eta++)
                {
                    omega_1 = omega_BZ[i](eta);
                    if (omega_1 != 0.0)
                    {
                        int xi = 1;
                        //for (int xi = 0; xi < 2; xi++)
                        {

                            dk = k[s] + k_BZ[i];
                            createDynamicalMatrix2D(dk, getNearestNeighbors2D());
                            omega_2 = getEigenValues(xi);

                            if (omega_2 != 0.0)
                            {
                                if (kroneker_G(k[s], k_BZ[i], k_BZ[i]))
                                {
                                    F = fabs(abs(getPotencial_3(moda_0, k[s],
                                                                eta, k_BZ[i],
                                                                xi,  dk,
                                                                getNearestNeighbors2D()   )));
                                    shift_1 += pow(F,2)/(omega_1*omega_2)*
                                            (
                                                (1/omega_1 + 1/omega_2)*
                                                (   f_lim(omega_1 + omega_2 + omega_0, smearing) +
                                                    f_lim(omega_1 + omega_2 - omega_0, smearing)) +
                                                (1/omega_1 - 1/omega_2)*
                                                (   f_lim(omega_1 - omega_2 + omega_0, smearing) -
                                                    f_lim(omega_2 - omega_1 + omega_0, smearing))
                                            );
                                    /*omega_p = (pow(omega_1*omega_1 - omega_2*ome
                                     * ga_2 - omega_0*omega_0, 2) - 4*omega_1*omega_1*omega_2*omega_2);
                                    shift_1 += pow(F, 2)/pow(omega_1*omega_2, 2) *
                                            (
                                                pow(omega_1*omega_1 - omega_2*omega_2, 2) - omega_0*omega_0*(omega_1*omega_1 + omega_2*omega_2)
                                            )
                                                *
                                            (
                                                omega_p / (omega_p*omega_p + smearing*smearing)
                                            );*/
                                }
                            }
                        }

                        k_neg_1 = (-1)*k[s];
                        k_neg_2 = (-1)*k_BZ[i];
                        if (kroneker_G(k[s], k_neg_1, k_BZ[i], k_neg_2))
                        {
                            F_2 = (getPotencial_4(moda_0, k[s],
                                                   moda_0,  k_neg_1,
                                                   eta, k_BZ[i],
                                                   eta, k_neg_2, getNearestNeighbors2D()));
                            shift_2 += F_2 / pow(omega_1, 2);
                        }
                    }
                }
            }

        }
        cout << "s = " << s << " / " << numWaveVectors << endl;
        if (omega_0 != 0) {
            Delta_k.push_back((-shift_1/4 + abs(shift_2)*massParticles/4/4)*(Temp/(pow(massParticles, 3)* numWaveVectors_BZ*omega_0)));
        }
    }
    recordVector(Delta_k, "output/", "Delta_k" + std::to_string(moda_0), "txt");
    return Delta_k;
 }
cellDoubles MeanSquareDisplacement::getFrequencyAttenuation(cellVectors& k, cellDoubles& omega, cellDoubles& omega_BZ)
{
	cellDoubles Delta_k;
	double F;
	double omega_1, omega_2;
	double domega_pos_1, domega_neg_1;
	double domega_pos_2, domega_neg_2;
	Eigen::VectorXd k_neg_1, k_neg_2;
   //Вариант из статьи.//
    /*auto f_Lim = [](double om, double delta) { return delta / pi() / (om*om + delta*delta); };

	for (int s = 1; s < k.size() - 1; s++)
	{
		shift_1 = 0;
		for (int i = 0; i < k_BZ.size(); i++)
		{
			omega_1 = omega_BZ[i];
			if (omega_1 != 0.0)
			{
				for (int j = 0; j < k_BZ.size(); j++)
				{
					omega_2 = omega_BZ[j];
					if (omega_2 != 0.0)
					{
						domega_pos_1 = omega_1 + omega_2 + omega[s];
						domega_neg_1 = omega_1 + omega_2 - omega[s];
						domega_pos_2 = omega_1 - omega_2 + omega[s];
						domega_neg_2 = omega_2 - omega_1 + omega[s];

						if (kroneker_G(k[s], k_BZ[i], k_BZ[j]))
						{
							F = fabs(abs(getPotencial_3(k[s], k_BZ[i], k_BZ[j], getNearestNeighbors())));
						}
						else
						{
							F = 0.0;
						}

						shift_1 += pow(F, 2) / (omega_1*omega_2)*(

							(1 / omega_2 + 1 / omega_1)*(
								f_Lim(domega_neg_1, smearing) - f_Lim(domega_pos_1, smearing)) +

								(1 / omega_2 - 1 / omega_1)*(
									f_Lim(domega_neg_2, smearing) - f_Lim(domega_pos_2, smearing))

							);
					}
				}
			}
		}

		cout << "s = " << s << " / " << k.size() << endl;
		if (omega[s] != 0) {
			Delta_k.push_back((Temp*math::pi()*(-2.25 * shift_1) / omega[s]) / pow(k_BZ.size(), 1.5));
		}
    }*/
	return Delta_k;
}
double MeanSquareDisplacement::getMSD(cellVectors& omega)
{
    double msd_1 = 0;
    int s_1 = 0;
	for (int i = 0; i < omega.size(); i++)
	{
        omega_k = omega[i](0);
		if (omega_k != 0.0)
		{
            msd_1 += 1 / (omega_k*omega_k);
            s_1++;
		}
	}
    double msd_2 = 0;
    int s_2 = 0;
    for (int i = 0; i < omega.size(); i++)
    {
        omega_k = omega[i](1);
        if (omega_k != 0.0)
        {
            msd_2 += 1 / (omega_k*omega_k);
            s_2++;
        }
    }
    return Temp*(msd_1/s_1 + msd_2/s_2) / (massParticles);
}
double MeanSquareDisplacement::getRelativeMSD(cellVectors& k, cellVectors& omega)
{
    double msd_1 = 0;
    int s_1 = 0;
    for (int i = 0; i < omega.size(); i++)
    {
        omega_k = omega[i](0);
        if (omega_k != 0.0)
        {
            msd_1 += 1 / (omega_k*omega_k)*(1 - cos(k_BZ[i].dot(getNearestNeighbors2D()[1])));
            s_1++;
        }
    }
    double msd_2 = 0;
    int s_2 = 0;
    for (int i = 0; i < omega.size(); i++)
    {
        omega_k = omega[i](1);
        if (omega_k != 0.0)
        {
            msd_2 += 1 / (omega_k*omega_k)*(1 - cos(k_BZ[i].dot(getNearestNeighbors2D()[1])));;
            s_2++;
        }
    }
    cout << Temp*msd_1/s_1 / (massParticles) << std::endl;
    cout << Temp*msd_2/s_2 / (massParticles) << std::endl;
    return Temp*(msd_1/s_1 + msd_2/s_2) / (massParticles);
}
double MeanSquareDisplacement::getMSDwithAnharmonicCorrection(cellVectors& k, cellVectors& omega_BZ)
{
    cellDoubles Eps_1;
    cellDoubles Eps_2;
    Eps_1 = getFrequencyShift(0, k, omega_BZ, omega_BZ);
    Eps_2 = getFrequencyShift(1, k, omega_BZ, omega_BZ);

    double msd_1 = 0;
    int s = 0;
    for (int i = 0; i < Eps_1.size(); i++)
    {
        omega_k = omega_BZ[i](0);
        if (omega_k != 0.0)
        {
            delta_k = Eps_1[i];
            //msd_1 += 1 / ((omega_k + delta_k)*(omega_k + delta_k));
            msd_1 += 1 / (omega_k*omega_k)*(1 + 2*delta_k / (omega_k))*(1 - cos(k_BZ[i].dot(getNearestNeighbors2D()[1])));
            s++;
        }
    }
    msd_1 = msd_1 /s;
    double msd_2 = 0;
    s=0;
    for (int i = 0; i < Eps_2.size(); i++)
    {
        omega_k = omega_BZ[i](1);
        if (omega_k != 0.0)
        {
            delta_k = Eps_2[i];
            //msd_2 += 1 / ((omega_k + delta_k)*(omega_k + delta_k));
            msd_2 += 1 / (omega_k*omega_k)*(1 + 2*delta_k / (omega_k))*(1 - cos(k_BZ[i].dot(getNearestNeighbors2D()[1])));
            s++;
        }
    }
    msd_2 = msd_2/s;
    return Temp*(msd_1 + msd_2) / (massParticles);
}
double MeanSquareDisplacement::getRelativeMSDwithAnharmonicCorrection(cellVectors& k, cellVectors& omega_BZ)
{
    cellDoubles Eps_1;
    cellDoubles Eps_2;
    Eps_1 = getFrequencyShift(0, k, omega_BZ, omega_BZ);
    Eps_2 = getFrequencyShift(1, k, omega_BZ, omega_BZ);

    double msd_1 = 0;
	int s = 0;
    for (int i = 0; i < Eps_1.size(); i++)
	{
        omega_k = omega_BZ[i](0);
		if (omega_k != 0.0)
		{
            delta_k = Eps_1[i];
            //msd_1 += 1 / ((omega_k + delta_k)*(omega_k + delta_k));
            msd_1 += 1 / (omega_k*omega_k)*(1 + 2*delta_k / (omega_k))*(1 - cos(k_BZ[i].dot(getNearestNeighbors2D()[1])));
			s++;
		}
	}
    msd_1 = msd_1 /s;
    double msd_2 = 0;
    s=0;
    for (int i = 0; i < Eps_2.size(); i++)
    {
        omega_k = omega_BZ[i](1);
        if (omega_k != 0.0)
        {
            delta_k = Eps_2[i];
            //msd_2 += 1 / ((omega_k + delta_k)*(omega_k + delta_k));
            msd_2 += 1 / (omega_k*omega_k)*(1 + 2*delta_k / (omega_k))*(1 - cos(k_BZ[i].dot(getNearestNeighbors2D()[1])));
            s++;
        }
    }
    msd_2 = msd_2/s;
    return Temp*(msd_1 + msd_2) / (massParticles);
}
cellDoubles MeanSquareDisplacement::computerMSD(cellDoubles& omega_L, cellDoubles& omega_T)
{
	cellDoubles msd;
	for (int i = 0; i < 10; i++)
	{
		Temp = i*0.1;
        //msd.push_back(getMSD(omega_L) + getMSD(omega_T));
	}
	return msd;
}
cellDoubles MeanSquareDisplacement::computerMSD_AE(cellVectors& k, cellDoubles& omega, cellDoubles& omega_BZ)
{
	cellDoubles msd;
	for (int i = 0; i < 10; i++)
	{
		Temp = i*0.1;
        //msd.push_back(getMSD_AE(k, omega, omega_BZ));
	}
	return msd;
}


