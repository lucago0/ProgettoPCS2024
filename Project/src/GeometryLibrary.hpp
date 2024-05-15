#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <Eigen/Eigen>
#include <map>


using namespace std;
using namespace Eigen;

namespace FracturesLib{

struct Fractures{
    unsigned int NumberFractures = 0;
    unsigned int Id;
    map<unsigned int, Matrix<double, 3, Dynamic>> Vertices;
    double tol = numeric_limits<double>::epsilon(); // precisione
    double tol_aree = tol*tol/2;
    map<unsigned int, array<double,4>> Coeff;
};

struct Traces{
    unsigned int Id;
    map<unsigned int, array<unsigned int, 2>> FracturesId;
    map<unsigned int, Vector3d> Vertices;
    bool Tips;
};

struct Line{
    array<double,3> point;
    array<double,3> direction;
};


bool importFracture(const string& filename, Fractures& fracture);
double distanceSquared(Vector3d& A, Vector3d& B);
void OutputFile(Traces& TR, Fractures& FR);
bool areClose(Fractures& mesh, unsigned int& Id1, unsigned int& Id2);
array<double,4> Piano(unsigned int& id, Fractures& FR);
array<double,6> Inter(array<double,4>& coeff1, array<double,4>& coeff2);
Matrix<double,4,4> PuntiIntersRetta(Fractures& fracture, unsigned int& Id1, unsigned int& Id2, array<double,6>& v);
array<double,4> intersection(const Matrix<double,4,4>&Q);

};
