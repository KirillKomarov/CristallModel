#include "stdafx.h"
#include "PTRecordFile.h"

void recordVector(Eigen::VectorXd& vector, std::string adress, std::string name, std::string format)
{
	std::string fullAdress = adress + name + "." + format;
	std::ofstream out(fullAdress, std::ios::out | std::ios::binary);

	for (int i = 0; i < vector.size(); i++)
	{
		out << vector(i) << "   ";
	}
	out.close();
}
void recordCollection(std::vector<Eigen::Vector3d>& collection, std::string adress, std::string name, std::string format, int s)
{
	std::string fullAdress = adress + name + "_" + std::to_string(s) + "." + format;
	std::ofstream out(fullAdress, std::ios::out | std::ios::binary);

	for (int i = 0; i < collection.size(); i++)
	{
		out << collection[i](0) << "   " << collection[i](1) << "   " << collection[i](2) << "\n";
	}
	out.close();
}
void recordCollection(std::vector<Eigen::VectorXd>& collection, std::string adress, std::string name, std::string format)
{
	std::string fullAdress = adress + name + "." + format;
	std::ofstream out(fullAdress, std::ios::out | std::ios::binary);

	for (int i = 0; i < collection.size(); i++)
	{
		for (int j = 0; j < collection[i].size(); j++)
		{
			out << collection[i](j) << "   ";
		}
		out << "\n";
	}

	out.close();
}
void recordCollection(std::vector<Eigen::Vector3d>& collection, std::string adress, std::string name, std::string format)
{
	std::string fullAdress = adress + name + "." + format;
	std::ofstream out(fullAdress, std::ios::out | std::ios::binary);

	for (int i = 0; i < collection.size(); i++)
	{
		out << collection[i](0) << "   " << collection[i](1) << "   " << collection[i](2) << "\n";
	}
	out.close();
}
void recordVector(Eigen::VectorXd& vector, std::string adress, std::string name, std::string format, int s)
{
	std::string fullAdress = adress + name + "_" + std::to_string(s) + "." + format;
	std::ofstream out(fullAdress, std::ios::out | std::ios::binary);

	for (int i = 0; i < vector.size(); i++)
	{
		out << vector(i) << "   ";
	}
	out.close();
}
void recordMatrix(Eigen::MatrixXd& matrix, std::string adress, std::string name, std::string format)
{
	std::string fullAdress = adress + name + "." + format;
	std::ofstream out(fullAdress, std::ios::out | std::ios::binary);

	for (int j = 0; j < matrix.cols(); j++)
	{
		for (int i = 0; i < matrix.rows(); i++)
		{
			out << matrix(i, j) << "   ";
		}
		out << "\n";
	}
	out.close();
}
void recordMatrix(Eigen::MatrixXd& matrix, std::string adress, std::string name, std::string format, int s)
{
	std::string fullAdress = adress + name + "_" + std::to_string(s) + "." + format;
	std::ofstream out(fullAdress, std::ios::out | std::ios::binary);
	for (int j = 0; j < matrix.cols(); j++)
	{
		for (int i = 0; i < matrix.rows(); i++)
		{
			out << matrix(i, j) << "   ";
		}
		out << "\n";
	}
	out.close();
}
void recordVector(std::vector<double> vector, std::string adress, std::string name, std::string format)
{
	std::string fullAdress = adress + name + "." + format;
	std::ofstream out(fullAdress, std::ios::out | std::ios::binary);
	for (int i = 1; i <= vector.size(); i++)
	{
		out << vector[i - 1] << "   ";
	}
	out.close();
}
