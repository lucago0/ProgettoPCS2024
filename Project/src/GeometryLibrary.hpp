#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <Eigen/Eigen>
#include <algorithm>
#include <map>


using namespace std;
using namespace Eigen;

namespace FracturesLib{

struct Fractures{
    unsigned int NumberFractures = 0;
    map<unsigned int, Matrix<double, 3, Dynamic>> Vertices;
    double tol = numeric_limits<double>::epsilon(); // precisione
    double tol_aree = tol*tol/2;
};

struct Traces{
    unsigned int NumberTraces = 0;
    map<unsigned int, array<unsigned int, 2>> FracturesId;
    map<unsigned int, Matrix<double,3,2>> Vertices;
    map<unsigned int, bool> Tips;
    map<unsigned int, double> Lengths;
};

struct Line{
    Vector3d point;
    Vector3d direction;
};


bool importFracture(const string& filename, Fractures& fracture);
double distanceSquared(const Vector3d& A,const Vector3d& B);
void OutputFile(Traces& TR, Fractures& FR);
bool areClose(Fractures& fracture, unsigned int& Id1, unsigned int& Id2);
Vector4d Piano(unsigned int& id, Fractures& FR);
Line Inter(const Vector4d &coeff1, const Vector4d &coeff2);
VectorXd PuntiIntersRetta(const Line& r,const Line& rj);
Vector4d intersection(const Vector4d& Q);
bool compareByValue(const std::pair<unsigned int, double> &a, const std::pair<unsigned int, double> &b);

