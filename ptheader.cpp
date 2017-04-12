#include "PTHeader.h"

cellVectors sumCellVectors(cellVectors& A, cellVectors& B)
{
    cellVectors C;
    try {
        if (A.size() != B.size())   throw "don't equal";

    }
    catch(std::string error)
    {
        std::cout<< "Error! " << error << std::endl;
        system("pause");
    }
    for (int i = 0; i<A.size(); i++)
    {
        //C.push_back(A[i] + B[i]);
    }
    return C;
}

cellVectors sumCellVectors(cellVectors& A, cellDoubles& a_1, cellDoubles& a_2)
{
    vector D(2);
    cellVectors C;
    C = A;
    try {
        if (A.size() != a_1.size())   throw "don't equal";

    }
    catch(std::string error)
    {
        std::cout<< "Error! " << error << std::endl;
        system("pause");
    }
    for (int i = 0; i<A.size(); i++)
    {
        D(0) = a_1[i];
        D(1) = a_2[i];
        D+=A[i];
        C.push_back(D);
    }
    return C;
}


cellDoubles sumCellDoubles(cellDoubles& A, cellDoubles& B)
{
    cellDoubles C;
    try {
        if (A.size() != B.size())   throw "don't equal";

    }
    catch(std::string error)
    {
        std::cout<< "Error! " << error << std::endl;
        system("pause");
    }
    for (int i = 0; i<A.size(); i++)
    {
        C.push_back(A[i] + B[i]);
    }
    return C;
}
