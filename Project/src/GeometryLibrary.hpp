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
    map<unsigned int, array<double,4>> Coeff;
};

struct Traces{
    unsigned int Id;
    map<unsigned int, array<unsigned int, 2>> FracturesId;
    map<unsigned int, Vector3d> Vertices;
    bool Tips;
};

struct Line{
<<<<<<< Updated upstream
    unsigned int id;
    array<double,3> point;
    array<double,3> direction;
=======
    unsigned int Id;
    unsigned int Id1;
    unsigned int Id2;
    map<unsigned int, array<array<double,3>,array<double,3>>> EqLine;
>>>>>>> Stashed changes
};

bool importFracture(const string& filename, Fractures& fracture);

bool areClose(unsigned int Id1, unsigned int Id2);

void OutputFile(Traces& TR, Fractures& FR);



double distanceSquared(Fractures& mesh, array<double,3> P1, array<double,3> P2);
};
