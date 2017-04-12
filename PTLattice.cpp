#include "stdafx.h"
#include "PTLattice.h"


void Lattice::createLattice()
{
    if (latticeType == "hexagonal") {
        if (latticeDim == "2d") {
            createHexagonalLattice2D(translation);
        }
        if (latticeDim == "3d") {
            createHexagonalLattice3D(translation);
        }
    }
    if (latticeType == "square") {
        if (latticeDim == "2d") {
            createSquareLattice2D(translation);
        }
        if (latticeDim == "3d") {
            createSquareLattice3D(translation);
        }
    }
}
void Lattice::createHexagonalLattice2D(double a)
{
    a_0.resize(2);
    a_0 << 0.0, 0.0;
    a_1.resize(2);
    a_1 << 0.5*a, 0.5*a*sqrt(3);
    a_2.resize(2);
    a_2 << 0.5*a, -0.5*a*sqrt(3);
    a_3.resize(2);
    a_3 << -a, 0;
    vol = a*a*sqrt(3)/2;
}
void Lattice::createSquareLattice2D(double a)
{
    a_0.resize(2);
    a_0 << 0.0, 0.0;
    a_1.resize(2);
    a_1 << a, 0.0;
    a_2.resize(2);
    a_2 << 0.0, a;
    vol = 0.5*a_1.dot(a_2);
}
void Lattice::createHexagonalLattice3D(double a)
{
    a_1.resize(3);
    a_1 << a, 0.0, 0.0;
    a_2.resize(3);
    a_2 << a*0.5, a*0.5*sqrt(3), 0.0;
    a_3.resize(3);
    a_3 << 0.0, 0.0, a;
    vol = a_1.dot(cross3D(a_2, a_3));
}
void Lattice::createSquareLattice3D(double a)
{
    a_0.resize(3);
    a_0 << 0.0, 0.0, 0.0;
    a_1.resize(3);
    a_1 << a, 0.0, 0.0;
    a_2.resize(3);
    a_2 << 0.0, a, 0.0;
    a_3.resize(3);
    a_3 << 0.0, 0.0, a;

    vol = a_1.dot(cross3D(a_2, a_3));
}
void Lattice::createHexagonalCristall2D()
{

    if (cristallType == "hexahedron")
    {
        atoms.push_back(a_0);
        for (int m = 0; m <= cristallRow; m++)
        {
            for (int n = 1; n <= cristallRow; n++)
            {
                atoms.push_back(n*a_1 + m*a_2);
                //atoms.push_back(n*a_2 + m*a_3);
                //atoms.push_back(n*a_3 + m*a_1);
             }
        }
        for (int m = 0; m <= cristallRow; m++)
        {
            for (int n
                 = 1; n <= cristallRow; n++)
            {
                atoms.push_back(n*a_2 + m*a_3);
             }
        }
        for (int m = 0; m <= cristallRow; m++)
        {
            for (int n
                 = 1; n <= cristallRow; n++)
            {
                atoms.push_back(n*a_3 + m*a_1);
             }
        }
    }
    if (cristallType == "rectangle")    {
        vector a_12 = a_1 + a_2;
        vector a_23 = a_2 - a_1;
        for (int m = -cristallWidth; m <= cristallWidth; m++)
        {
            for (int n = -cristallLength; n <= cristallLength; n++)
            {
                atoms.push_back(n*a_12 + m*a_23);
            }
        }
        for (int m = -cristallWidth; m <= cristallWidth; m++)
        {
            for (int n = -cristallLength; n <= cristallLength; n++)
            {
                atoms.push_back(n*a_12 + m*a_23 + a_2);
            }
        }
    }
    cout << "Creat is " << atoms.size() << " atoms" << endl;
}
void Lattice::createCristall2D(double num_x, double num_y)
{
    for (int n = 0; n < num_y; n++)
    {
        for (int m = 0; m < num_x; m++)
        {
            atoms.push_back(n*a_1 + m*a_2);
        }
    }
}
void Lattice::createCristall3D(double num_x, double num_y, double num_z)
{
    for (int k = 0; k < num_z; k++)
    {
        for (int n = 0; n < num_y; n++)
        {
            for (int m = 0; m < num_x; m++)
            {
                atoms.push_back(n*a_1 + m*a_2 + k*a_3);
            }
        }
    }
}

void Lattice::setReciprocalVectors2D()
{
    double coef = 1.0;//1.1547;
    double norm = 4 * pi() / translation/sqrt(3);
    if (latticeType == "hexagonal")
    {
        b_1 = 2 * pi() / vol*cross2D(a_1)*coef;
        b_2 = 2 * pi() / vol*cross2D(a_2)*coef;
        b_3 = 2 * pi() / vol*cross2D(a_3)*coef;
        /*b_1 = cross2D(a_1)/cross2D(a_1).norm()*norm;
        b_2 = cross2D(a_2)/cross2D(a_2).norm()*norm;
        b_3 = cross2D(a_3)/cross2D(a_3).norm()*norm;*/
    }
    else
    {
        b_1 = 2 * pi() / vol*cross2D(a_1);
        b_2 = 2 * pi() / vol*cross2D(a_2);
    }
}
void Lattice::setReciprocalVectors3D()
{
    b_1 = 2 * pi() / vol*cross3D(a_2, a_3);
    b_2 = 2 * pi() / vol*cross3D(a_3, a_1);
    b_3 = 2 * pi() / vol*cross3D(a_1, a_2);
}
void Lattice::setFirstBrillouinZone2D()
{
    setReciprocalVectors2D();

    b_1 = getUnitVector(a_1)*b_1.norm()/sqrt(3);
    b_2 = getUnitVector(a_2)*b_2.norm()/sqrt(3);
    b_3 = getUnitVector(a_3)*b_3.norm()/sqrt(3);

    if (cristallType == "hexahedron")
    {
        if (latticeType == "hexagonal")
        {
            for (int m = 0; m <= cristallRow; m++)
            {
                for (int n = 1; n <= cristallRow; n++)
                {
                    k_BZ.push_back((n*b_1 + m*b_2)/(cristallRow));
                    //k_BZ.push_back((n*b_2 + m*b_3)/cristallRow);
                    //k_BZ.push_back((n*b_3 + m*b_1)/cristallRow);
                }
            }
            for (int m = 0; m <= cristallRow; m++)
            {
                for (int n = 1; n <= cristallRow; n++)
                {
                    k_BZ.push_back((n*b_2 + m*b_3)/(cristallRow));
                }
            }
            for (int m = 0; m <= cristallRow; m++)
            {
                for (int n = 1; n <= cristallRow; n++)
                {
                    k_BZ.push_back((n*b_3 + m*b_1)/(cristallRow));
                }
            }
        }
    }
    if (cristallType == "rectangle")
    {
        if (latticeType == "hexagonal")
        {
            for (int m = 0; m <= cristallRow; m++)
            {
                for (int n = 1; n <= cristallRow; n++)
                {
                    k_BZ.push_back((n*b_1 + m*b_2)/cristallRow*2/3);
                    k_BZ.push_back((n*b_2 + m*b_3)/cristallRow*2/3);
                    k_BZ.push_back((n*b_3 + m*b_1)/cristallRow*2/3);
                }
            }
        }
    }
    cout << "Creat is " << k_BZ.size() << " vectors" << endl;
    recordCollection(k_BZ, "output/", "BZ", "txt");
}


cellVectors Lattice::getTranslationVectors()
{
    cellVectors tv;
    tv.push_back(a_1);
    tv.push_back(a_2);
    tv.push_back(a_3);
    return tv;
}
cellVectors Lattice::getNearestNeighbors2D()
{
    cellVectors nn;
    if (latticeType == "hexagonal")
    {
        nn.push_back(a_1);
        nn.push_back(a_2);
        nn.push_back(a_3);

        nn.push_back(a_1 + a_2);
        nn.push_back(a_2 + a_3);
        nn.push_back(a_3 + a_1);
    }
    if (latticeType == "square")
    {
        nn.push_back(a_1);
        nn.push_back(a_2);
        nn.push_back(a_1 + a_2);
        nn.push_back(a_1 - a_2);

        nn.push_back(-a_1);
        nn.push_back(-a_2);
        nn.push_back(-a_1 - a_2);
        nn.push_back(-a_1 + a_2);
    }
    return nn;
}
cellVectors Lattice::getFirstBrillouinZone2D()
{
    setReciprocalVectors2D();

    b_1 = getUnitVector(a_1)*b_1.norm()*2/sqrt(3);
    b_2 = getUnitVector(a_2)*b_2.norm()*2/sqrt(3);
    b_3 = getUnitVector(a_3)*b_3.norm()*2/sqrt(3);


    cellVectors k_bz;
    if (latticeType == "hexagonal")
    {
        for (int m = 0; m <= cristallRow; m++)
        {
            for (int n = 1; n <= cristallRow; n++)
            {
                k_bz.push_back(n*b_1 + m*b_2);
                //k_bz.push_back(n*b_2 + m*b_3);
                //k_bz.push_back(n*b_3 + m*b_1);
             }
        }
        for (int m = 0; m <= cristallRow; m++)
        {
            for (int n = 1; n <= cristallRow; n++)
            {
                k_bz.push_back(n*b_2 + m*b_3);
             }
        }
        for (int m = 0; m <= cristallRow; m++)
        {
            for (int n = 1; n <= cristallRow; n++)
            {
                k_bz.push_back(n*b_3 + m*b_1);
             }
        }
    }
    cout << "Creat is " << k_bz.size() << " vectors" << endl;
    recordCollection(k_bz, "output/", "BZ", "txt");
    return k_bz;
}
cellVectors Lattice::getDirectionGK(int numPoints)
{
    double dk_x, dk_y;
    vector bz_1, bz_2, bz_3;
    setReciprocalVectors2D();

    bz_1 = getUnitVector(a_1)*b_1.norm();
    bz_2 = getUnitVector(a_2)*b_2.norm();
    bz_3 = getUnitVector(a_3)*b_3.norm();

/*
    bz_1 = b_1;
    bz_2 = b_2;
    bz_3 = b_3;
*/

    Eigen::VectorXd k_(2);
    cellVectors k_bz;
    //if ("GK" == dir)
    {
        //k_ = 0.5*(2*b_1 + b_2);
        //k_ = pi() / fabs(k_.dot(a_1))*k_*2/3;
        k_ = bz_1*2/sqrt(3);
        cout << k_.norm() << endl;
        dk_x = k_(0) / (numPoints - 1);
        dk_y = k_(1) / (numPoints - 1);
        for (int i = 0; i < numPoints; i++)
        {
            k_(0) = dk_x*i;
            k_(1) = dk_y*i;
            k_bz.push_back(k_);
        }
    }
    return k_bz;
}
cellVectors Lattice::getDirectionMG(int numPoints)
{
    setReciprocalVectors2D();

    Eigen::VectorXd k_(2); //k_ << 0.0, 0.0;
    cellVectors k_bz;

    //if (dir == "MG")
    {
        Eigen::VectorXd k_mg = b_2;
        k_mg = pi() / fabs(k_mg.dot(a_2))*k_mg;
        Eigen::VectorXd k_abs = k_mg;
        Eigen::VectorXd dk_mg = k_mg / (numPoints - 1);
        for (int i = 0; i < numPoints; i++)
        {
            k_mg = k_abs - dk_mg*i;
            k_bz.push_back(k_mg);
        }
    }

    return k_bz;
}
cellVectors Lattice::getDirectionGKMG(int numPoints)
{
    double dk_x, dk_y;
    vector bz_1, bz_2, bz_3;
    setReciprocalVectors2D();

    bz_1 = getUnitVector(a_1)*b_1.norm()/sqrt(3);
    bz_2 = getUnitVector(a_2)*b_2.norm()/sqrt(3);
    bz_3 = getUnitVector(a_3)*b_3.norm()/sqrt(3);

    Eigen::VectorXd k_(2);
    cellVectors k_bz;
    //if ("GK" == dir)
    {
        //k_ = 0.5*(2*b_1 + b_2);
        //k_ = pi() / fabs(k_.dot(a_1))*k_*2/3;
        k_ = bz_1;
        cout << k_.norm() << endl;
        dk_x = k_(0) / (numPoints - 1);
        dk_y = k_(1) / (numPoints - 1);
        for (int i = 0; i < numPoints; i++)
        {
            k_(0) = dk_x*i;
            k_(1) = dk_y*i;
            k_bz.push_back(k_);
        }
    }

    //if (dir == "KM")
    {
        Eigen::VectorXd k_gk = bz_1;
        Eigen::VectorXd k_mg = 0.5*((bz_1 + bz_2) + bz_1);
        //k_mg = pi() / fabs(k_mg.dot(a_2))*k_mg;

        Eigen::VectorXd k_km = k_mg - k_gk;
        Eigen::VectorXd dk_km = k_km / (numPoints - 1);
        for (int i = 1; i < numPoints; i++)
        {
            k_km = k_gk + dk_km*i;
            k_bz.push_back(k_km);
        }
    }


    //if (dir == "MG")
    {
        Eigen::VectorXd k_mg =  0.5*((bz_1 + bz_2) + bz_1);
        Eigen::VectorXd k_abs = k_mg;
        Eigen::VectorXd dk_mg = k_mg / (numPoints - 1);
        cout << k_mg.norm() << endl;
        for (int i = 1; i < numPoints; i++)
        {
            k_mg = k_abs - dk_mg*i;
            k_bz.push_back(k_mg);
        }
    }

    recordCollection(k_bz, "output/", "BZ", "txt");
    return k_bz;
}
bool Lattice::kroneker_G(Eigen::VectorXd& k_1, Eigen::VectorXd& k_2, Eigen::VectorXd& k_3)
{
    //if (approx(exp(math::i().imag()*(k_1 + k_2 + k_3).dot(a_1 + a_2)), 1.0, 5e-1) == 0	)
    if ((k_1 + k_2 + k_3).norm() < pi())
    {
        return true;
    }
    else
    {
        return false;
    }
}
bool Lattice::kroneker_G(Eigen::VectorXd& k_1, Eigen::VectorXd& k_2, Eigen::VectorXd& k_3, Eigen::VectorXd& k_4)
{
    if ((k_1 + k_2 + k_3 + k_4).norm() < pi())
    {
        return true;
    }
    else
    {
        return false;
    }
}
