#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace FracturesLib{

struct Fractures{
    unsigned int NumberFractures = 0;
    map<unsigned int, vector<array<double, 3>>> Vertices;
    double tol = numeric_limits<double>::epsilon(); // precisione
    double tol_aree = tol*tol/2;
};

struct Traces{
    unsigned int Id;
    map<unsigned int, array<unsigned int, 2>> FracturesId;
    map<unsigned int, array<double, 3>> Vertices;
    bool Tips;
};
};
