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
    map<unsigned int, Matrix<double, 3, Dynamic>> Vertices;
    double tol = numeric_limits<double>::epsilon(); // precisione
    double tol_aree = tol*tol/2;
    VectorXi NumTracce;
    map<unsigned int, vector<tuple<unsigned int, bool, double>>> tracce;
    // map<unsigned int, Matrix<double, 3, Dynamic>> SubFractures = {};
};

struct Traces{
    unsigned int NumberTraces = 0;
    vector<array<unsigned int, 4>> FracturesId; //{id1,id2,tips1,tips2}
    map<unsigned int, Matrix<double,3,2>> Vertices;
    map<unsigned int, double> Lengths;
};

struct Line{
    Vector3d point;
    Vector3d direction;
};

struct PolygonalMesh{
    unsigned int NumberCell0Ds = 0;
    vector<unsigned int> IdCell0Ds = {};
    vector<Vector2d> CoordinateCell0Ds = {};

    unsigned int NumberCell1Ds = 0;
    vector<unsigned int> IdCell1Ds = {};
    vector<Vector2i> VerticesCell1Ds = {};

    unsigned int NumberCell2Ds = 0;
    vector<unsigned int> IdCell2Ds = {};
    vector<VectorXi> VerticesCell2Ds = {};
    vector<VectorXi> EdgesCell2Ds = {};

    /* double tau = numeric_limits<double>::epsilon();
    double tol = tau*tau/2; */
};
}

bool importFracture(const string& filename, Fractures& fracture);
double distanceSquared(const Vector3d& A,const Vector3d& B);
void OutputFile(Traces& TR, Fractures& FR);
bool areClose(Fractures& fracture, unsigned int& Id1, unsigned int& Id2);
Vector4d Piano(unsigned int& id, Fractures& FR);
Line Inter(const Vector4d &coeff1, const Vector4d &coeff2);
VectorXd PuntiIntersRetta(const Line& r,const Line& rj);
Vector4d intersection(const Vector4d& Q);
bool almostEqual(double a, double b, double tol);
bool arePlanesParallel(Vector4d v1, Vector4d v2, double tol);
bool compareByValue(const pair<unsigned int, double> &a, const pair<unsigned int, double> &b);
bool compareTuple(const tuple<unsigned int, bool, double>& a, const tuple<unsigned int, bool, double>& b);

};
