#include "stdafx.h"
#include "PTInput.h"

void MDPInput::getInputData()
{
	double f;
	std::string I;
	std::string dim;
    std::string fullAdress = "GInputData.txt";
    std::fstream input;
	input.open(fullAdress);
	int s = 0;
	if (input)
	{
		while (!input.eof())
        {
            input >> I;
			if (I == "fileName") { input >> I; fileName = I; }
			if (I == "timeStep") { input >> f; input >> dim; timeStep = f*getDimTime(dim); }
			if (I == "numTimeStep") { input >> f; input >> dim; numTimeStep = f; }
			if (I == "Particles") {
				input >> I; input >> I; input >> f; input >> dim; massParticles = f*getDimMass(dim);
				input >> I; input >> I; input >> f; numParticles = f;
				input >> I; input >> I; input >> f; Temp = f;
			}
			if (I == "Potential") {
                input >> I; potentialType = I;
                if (I == "IPL")	{input >> f; numIPL = f;}
				input >> I; input >> I; input >> f; input >> dim; energyScale = f;
				input >> I; input >> I; input >> f; input >> dim; lengthScale = f;
			}
			if (I == "Lattice") {
				input >> I; latticeType = I;	input >> I;	latticeDim = I;
				input >> I; 
                if (I == "a")	{input >> I; input >> f; input >> dim; translation = f*getDimLength(dim);}
                if (I == "rho") {input >> I; input >> f; input >> dim; translation = f*getDimLength(dim);}
                input >> I; cristallType = I;
                if (I == "hexahedron")	{input >> f; cristallRow = f;}
                if (cristallType == "rectangle")	{input >> f; cristallLength = f;  input >> f; cristallWidth = f;}
			}
			s++;
		}
	}
	input.close();
}
