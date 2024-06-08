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
    vector<Matrix<double, 3, Dynamic>> Vertices;
    vector<Vector4d> CoeffPiano;
    vector<unsigned int> NumTracce = {};
    vector<vector<tuple<unsigned int, bool, double>>> tracce; //ad ogni id frattura associo vettore di idtracce con le info per ogni traccia
    vector<vector<Matrix<double, 3, Dynamic>>> SubFracVert; //vertici sotto frattura
};

struct SubFracture{
    Matrix<double, 3 ,Dynamic> Vertices;
    vector<unsigned int> VerticesId;
    vector<unsigned int> EdgesId;
    vector<unsigned int> traceId;
};

struct Traces{
    unsigned int NumberTraces = 0;
    vector<array<unsigned int, 4>> FracturesId; //{id1,id2,tips1,tips2} (id1 e id2 sono gli id delle fratture associate ad idtraccia
    vector<Matrix<double,3,2>> Vertices;
    //map<unsigned int, array<bool,2>> Tips;
    map<unsigned int, double> Lengths;
};

struct Line{
    Vector3d point;
    Vector3d direction;
};
struct PolygonalMesh{
    unsigned int NumberCell0Ds = 0;
    vector<unsigned int> IdCell0Ds = {};
    vector<Vector3d> CoordinateCell0Ds = {};

    unsigned int NumberCell1Ds = 0;
    vector<unsigned int> IdCell1Ds = {};
    vector<array<unsigned int,2>> VerticesCell1Ds = {};

    unsigned int NumberCell2Ds = 0;
    vector<unsigned int> IdCell2Ds = {};
    vector<vector<unsigned int>> VerticesCell2Ds = {};
    vector<vector<unsigned int>> EdgesCell2Ds = {};

    /* double tau = numeric_limits<double>::epsilon();
    double tol = tau*tau/2; */
};


bool importFracture(const string& filename, Fractures& fracture);
double distanceSquared(const Vector3d& A,const Vector3d& B);
void OutputFile(Traces& TR, Fractures& FR);
bool areClose(Fractures& fracture,const unsigned int& Id1, const unsigned int& Id2, const double& tol);
Vector4d Piano(const unsigned int& id, Fractures& FR);
Line Inter(const Vector4d &coeff1, const Vector4d &coeff2, const double &tol);
VectorXd PuntiIntersRetta(const Line& r,const Line& rj);
Vector4d intersection(const Vector4d& Q, const double& tol);
bool almostEqual(const double a,const double b,const double& tol);
bool arePlanesParallel(const Vector4d v1, const Vector4d v2, double& tol);
bool compareByValue(const pair<unsigned int, double> &a, const pair<unsigned int,const double> &b);
bool compareTuple(const tuple<unsigned int, bool, double>& a, const tuple<unsigned int, bool, double>& b);
double posizionePuntoPiano(const Vector4d& coeffPiano, const Vector3d& coordPunto);
void splitSubFractures(SubFracture& subFract, const Fractures& fractures, const Traces& traces, PolygonalMesh &mesh, const unsigned int &idFrac, double& tol);


};
