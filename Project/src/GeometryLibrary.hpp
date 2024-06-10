#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <Eigen/Eigen>
#include <map>


using namespace std;
using namespace Eigen;

namespace fracturesLib{

struct Fractures{
    unsigned int numberOfFractures = 0;
    vector<Matrix<double, 3, Dynamic>> vertices;
    vector<Vector4d> planeCoeffs;
    vector<unsigned int> numberOfTraces = {};
    vector<vector<tuple<unsigned int, bool, double>>> traces; //ad ogni id frattura associo vettore di idtracce con le info per ogni traccia
    vector<vector<Matrix<double, 3, Dynamic>>> SubFracVert; //vertici sotto frattura
};

struct SubFracture{
    Matrix<double, 3 ,Dynamic> Vertices;
    vector<unsigned int> VerticesId;
    vector<unsigned int> EdgesId;
    vector<unsigned int> traceId;
};

struct Traces{
    unsigned int numberOfTraces = 0;
    vector<array<unsigned int, 4>> fracturesId; //{id1,id2,tips1,tips2} (id1 e id2 sono gli id delle fratture associate ad idtraccia
    vector<Matrix<double,3,2>> vertices;
    //map<unsigned int, array<bool,2>> Tips;
    map<unsigned int, double> lengths;
};

struct line{
    Vector3d point;
    Vector3d direction;
};
struct PolygonalMesh{
    unsigned int numberCell0Ds = 0;
    vector<unsigned int> idCell0Ds = {};
    vector<Vector3d> coordinateCell0Ds = {};

    unsigned int numberCell1Ds = 0;
    vector<unsigned int> idCell1Ds = {};
    vector<array<unsigned int,2>> verticesCell1Ds = {};

    unsigned int numberCell2Ds = 0;
    vector<unsigned int> idCell2Ds = {};
    vector<vector<unsigned int>> verticesCell2Ds = {};
    vector<vector<unsigned int>> edgesCell2Ds = {};

    /* double tau = numeric_limits<double>::epsilon();
    double tol = tau*tau/2; */
};


bool importFractures(const string& filename, Fractures& fractures);
double distanceSquared(const Vector3d& a, const Vector3d& b);
void outputFile(Traces& traces, Fractures& fractures);
bool areClose(Fractures& fractures,const unsigned int& id1, const unsigned int& id2, const double& tol);
Vector4d plane(const unsigned int& id, Fractures& fractures);
line planesIntersection(const Vector4d& coeff1, const Vector4d& coeff2, const double& tol);
VectorXd linesIntersection(const line& r,const line& rj);
Vector4d intervalsIntersection(const Vector4d& q, const double& tol);
bool almostEqual(const double& a, const double& b,const double& tol);
bool arePlanesParallel(const Vector4d& v1, const Vector4d& v2, double& tol);
bool compareByValue(const pair<unsigned int, double>& a, const pair<unsigned int,const double>& b);
bool compareTuple(const tuple<unsigned int, bool, double>& a, const tuple<unsigned int, bool, double>& b);
// double posizionePuntoPiano(const Vector4d& coeffPiano, const Vector3d& coordPunto);
// void splitSubfractures(SubFracture& subFract, const fractures& fractures, const Traces& traces, PolygonalMesh &mesh, const unsigned int &idFrac, double& tol);


};
