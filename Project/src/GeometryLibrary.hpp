#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace FracturesLib{

struct Fractures{
    unsigned int NumberFractures = 0;
    unsigned int Id;
    map<unsigned int, Matrix<double, 3, Dynamic>> Vertices;
    double tol = numeric_limits<double>::epsilon(); // precisione
    double tol_aree = tol*tol/2;
};

struct Traces{
    unsigned int Id;
    map<unsigned int, array<unsigned int, 2>> FracturesId;
    map<unsigned int, array<double, 3>> Vertices;
    bool Tips;
};

bool importFracture(const string& filename, Fractures& fracture);

bool areClose(unsigned int Id1, unsigned int Id2);

double distanceSquared(Fractures& mesh, array<double,3> P1, array<double,3> P2);
};
