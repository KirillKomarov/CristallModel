#pragma once

#ifndef PHONONSPECTRA_H_
#define PHONONSPECTRA_H_

#include "PTLattice.h"
#include "PTHeader.h"
#include "PTPotentials.h"
#include "PTRecordFile.h"

class PhononSpectra : public Lattice, public DynamicalMatrix {
	std::vector<Eigen::VectorXd> k;
	std::vector<Eigen::VectorXd> r;
	std::vector<Eigen::VectorXd> v;
public:
	int numPoints;
	std::vector<std::complex<double>> mode_L, mode_T;
	PhononSpectra(int);
	~PhononSpectra();
	void loadMDfile(int);
	void getFrequencyAxis();
	cellDoubles computerSpectraBvK(cellVectors&, int);
    cellVectors computerSpectraBvK(cellVectors&);
	cellDoubles analiticSpectraBvK(cellVectors&, int);
};

#endif // !PHONONSPECTRA_H_
