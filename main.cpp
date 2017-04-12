
#include "stdafx.h"
#include "mainwindow.h"
#include <QApplication>

#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <Eigen/Eigen>

#include "PTRecordFile.h"
#include "PTPhononSpectra.h"
#include "PTMeanSquareDisplacement.h"
#include "PTPlotDate.h"

using namespace math;
using std::cout; using std::endl;


int main(int argc, char *argv[])
{
/*
    MeanSquareDisplacement msd;
    msd.createHexagonalCristall2D();
    plotTarget(msd.atoms, argc, argv);
*/



    PhononSpectra BvK(1.0);
    cellDoubles omega = BvK.computerSpectraBvK(BvK.getDirectionGKMG(20), 0);
    plotVector(omega, argc, argv);
    recordVector(omega, "output/", "omega_BvK_0", "txt");
    omega = BvK.computerSpectraBvK(BvK.getDirectionGKMG(20), 1);
    plotVector(omega, argc, argv);
    recordVector(omega, "output/", "omega_BvK_1", "txt");


/*
    MeanSquareDisplacement msd;
    PhononSpectra BvK(1.0);
    cellDoubles omega = BvK.computerSpectraBvK(msd.k_BZ, 0);
    recordVector(omega, "output/", "omega_BvK_0", "txt");
    omega = BvK.computerSpectraBvK(msd.k_BZ, 1);
    recordVector(omega, "output/", "omega_BvK_1", "txt");
*/


/*
    PhononSpectra BvK(1.0);
    MeanSquareDisplacement msd;
    cout << msd.k_BZ.size() << endl;
    cout << msd.getRelativeMSD(msd.k_BZ, BvK.computerSpectraBvK(msd.k_BZ)) << endl;
    cout << msd.getRelativeMSDwithAnharmonicCorrection(msd.k_BZ, BvK.computerSpectraBvK(msd.k_BZ)) << endl;
*/

/*
    PhononSpectra BvK(1.0);
    cellVectors k = BvK.getDirectionGKMG(50);
    cellDoubles Delta_k = msd.getFrequencyShift(1, k, BvK.computerSpectraBvK(k), BvK.computerSpectraBvK(msd.k_BZ));
    plotVector(Delta_k, argc, argv);
*/


    return EXIT_SUCCESS;
}

