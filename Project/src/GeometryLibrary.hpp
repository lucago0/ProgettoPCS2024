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
};

struct Traces{
    unsigned int numberOfTraces = 0;
    vector<array<unsigned int, 4>> fracturesId; //{id1,id2,tips1,tips2} id1 e id2 sono gli id delle fratture associate ad idtraccia
    vector<Matrix<double,3,2>> vertices;
    map<unsigned int, double> lengths;
};

struct line{
    Vector3d point;
    Vector3d direction;
};
struct PolygonalMesh{
    unsigned int numberCell0Ds = 0;
    vector<Vector3d> coordinateCell0Ds = {};

    unsigned int numberCell1Ds = 0;
    vector<array<unsigned int,2>> verticesCell1Ds = {};
    vector<vector<unsigned int>> neighCell1Ds = {};
    vector<bool> isOn1D;

    unsigned int numberCell2Ds = 0;
    vector<vector<unsigned int>> verticesCell2Ds = {};
    vector<vector<unsigned int>> edgesCell2Ds = {};
    vector<bool> isOn2D;

    /* double tau = numeric_limits<double>::epsilon();
    double tol = tau*tau/2; */
};


bool importFractures(const string& filename, Fractures& fractures);
void outputFile(Traces& traces, Fractures& fractures);
line planesIntersection(const Vector4d& coeff1, const Vector4d& coeff2, const double& tol);
VectorXd linesIntersection(const line& r,const line& rj);
void print(PolygonalMesh& mesh);
PolygonalMesh mergeMesh(vector<PolygonalMesh>& finalMesh);

inline double distanceSquared(const Vector3d& a,const Vector3d& b){
    return pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2);
}

inline bool areClose(Fractures& fractures,const unsigned int& id1, const unsigned int& id2, const double& tol){
    Vector3d c1 = {0,0,0};
    const unsigned int n1 = fractures.vertices[id1].cols();
    Vector3d c2 = {0,0,0};
    const unsigned int n2 = fractures.vertices[id2].cols();

    for (unsigned int i = 0; i < 3; i++){
        for (unsigned int j = 0; j < n1; j++){
            c1[i] += fractures.vertices[id1](i,j);
        }
        c1[i] /= n1;
        for (unsigned int j = 0; j < n2; j++){
            c2[i] += fractures.vertices[id2](i,j);
        }
        c2[i] /= n2;
    }

    Vector3d vertex;
    VectorXd raysSquared1;
    raysSquared1.resize(n1);
    for (unsigned int i = 0; i < n1; i++){
        vertex = fractures.vertices[id1].col(i);
        raysSquared1[i] = distanceSquared(c1,vertex);
    }
    VectorXd raysSquared2;
    raysSquared2.resize(n2);
    for (unsigned int i = 0; i < n2; i++){
        vertex = fractures.vertices[id2].col(i);
        raysSquared2[i] = distanceSquared(c2,vertex);
    }

    double r1 = *max_element(raysSquared1.begin(), raysSquared1.end());
    double r2 = *max_element(raysSquared2.begin(), raysSquared2.end());

    return distanceSquared(c1,c2) <= (r1 + r2 + 2*sqrt(r1)*sqrt(r2) + tol);
}

inline Vector4d plane(const unsigned int& id, Fractures& fractures)
{
    const Vector3d v0 = fractures.vertices[id].col(0);
    const Vector3d v1 = fractures.vertices[id].col(1);
    const Vector3d v2 = fractures.vertices[id].col(2);

    Vector4d coeff = {0,0,0,0};

    coeff[0] = (v1[1]-v0[1])*(v2[2]-v0[2])-(v1[2]-v0[2])*(v2[1]-v0[1]);
    coeff[1] = -((v1[0]-v0[0])*(v2[2]-v0[2])-(v2[0]-v0[0])*(v1[2]-v0[2]));
    coeff[2] = (v1[0]-v0[0])*(v2[1]-v0[1])-(v2[0]-v0[0])*(v1[1]-v0[1]);
    coeff[3] = -(coeff[0]*v0[0]+coeff[1]*v0[1]+coeff[2]*v0[2]); //da capire i segni

    return coeff;
}

inline bool compareByValue(const pair<unsigned int, double> &a, const pair<unsigned int, const double> &b) {
    return a.second > b.second;
}

inline bool almostEqual(const double& a, const double& b, const double& tol) {
    return fabs(a - b) < tol;
}

inline bool arePlanesParallel(const Vector4d& v1, const Vector4d& v2, double& tol) {

    Vector3d w1 = v1.head(3);
    Vector3d w2 = v2.head(3);
    Vector3d v = w1.cross(w2);

    bool par = true;

    if(abs(v[0]) > tol || abs(v[1]) > tol || abs(v[2]) > tol)
        par = false;

    return par;
}

inline bool compareTuple(const tuple<unsigned int, bool, double>& a, const tuple<unsigned int, bool, double>& b) {
    return get<1>(a) < get<1>(b); // Ordina in base al valore booleano all'interno delle tuple
}

inline Vector4d intervalsIntersection(const Vector4d& q,const double& tol){
    double a = min(q[0],q[1]);
    double b = max(q[0],q[1]);
    double c = min(q[2],q[3]);
    double d = max(q[2],q[3]);

    // Calcola l'estremo sinistro dell'intersezione
    double sx = max(a, c);
    // Calcola l'estremo destro dell'intersezione
    double dx = min(b, d);
    //Se gli intervalli non si sovrappongono, l'intersezione sarà vuota
    if (sx > (dx-tol)) {
        sx = dx = numeric_limits<double>::quiet_NaN(); // Non un numero
    }
    double other_sx = (a < (c+tol)) ? a : c; // other_sx è pari ad a se a<c, altrimenti è pari a c
    double other_dx = (d > (b-tol)) ? d : b; // other_dx è pari a d se d>b, altrimenti è pari a b
    Vector4d output = {sx,dx,other_sx,other_dx}; // in ordine restituiamo l'intervallo di intersezione e gli altri due estremi ordinati
    return output;
}

};
