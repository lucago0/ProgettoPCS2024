#include <GeometryLibrary.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <map>

using namespace std;
using namespace Eigen;

namespace fracturesLib{

bool importFractures(const string& filename, Fractures& fractures) {
    ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    string header;
    getline(file, header); // Read the header line to discard it

    string line;
    getline(file,line);
    fractures.numberOfFractures = stoi(line);
    int numberOfFractures = stoi(line);
    fractures.vertices.resize(numberOfFractures);
    char pv;
    int numberOfVertices;
    unsigned int k;

    while (numberOfFractures--) {
        getline(file, line);
        getline(file, line);
        stringstream ss(line);
        ss >> k >> pv >> numberOfVertices;

        Matrix<double, 3, Dynamic> actualVert(3, numberOfVertices);

        getline(file, line);
        string val;
        // Read Vertices
        for(int j = 0; j < 3; j++){
            for(int i = 0; i < numberOfVertices; i++){
                file >> val;
                actualVert(j, i) = stod(val);
            }
        }
        file >> pv;

        // Add fracture data to the fractures map
        fractures.vertices[k] = actualVert;
    }
    file.close();
    return true;
}

double distanceSquared(const Vector3d& a,const Vector3d& b){
    return pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2);
}

bool compareByValue(const pair<unsigned int, double> &a, const pair<unsigned int, const double> &b) {
    return a.second > b.second;
}

void outputFile(Traces& traces, Fractures& fractures)
{
    string outputFracts = "Traces.txt";
    string outputFractsTraces = "Fratture-Tracce.txt";
    ofstream ofs1(outputFracts);
    ofstream ofs2(outputFractsTraces);

    if (ofs1.fail() || ofs2.fail())
    {
        cerr << "Error in creating the output file" << endl;
        return;
    }

    ofs1 << "# Number of Traces" << endl;
    ofs1 << traces.numberOfTraces << endl;
    ofs1 << "# TraceId; FracturesId1; FracturesId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

    for(unsigned int i = 0; i < traces.numberOfTraces;i++)
    {
        ofs1 << i << ";" << traces.fracturesId[i][0] << ";" << traces.fracturesId[i][1] << ";" << traces.vertices[i](0,0) << ";" << traces.vertices[i](1,0) << ";" << traces.vertices[i](2,0) << ";" << traces.vertices[i](0,1) << ";" << traces.vertices[i](1,1) << ";" << traces.vertices[i](2,1) << endl;
    }

    for(unsigned int i = 0; i < fractures.numberOfFractures; i++)
    {
        ofs2 << "# FractureId; NumTraces" << endl;
        ofs2 << i << ";" << fractures.numberOfTraces[i] << endl;
        ofs2 << "# TraceId; Tips; Length" << endl;
        for (const auto& elem : fractures.traces[i]) {
            ofs2 << get<0>(elem) << ";" << get<1>(elem) << ";" << get<2>(elem) << endl;
        }
    }


}

bool areClose(Fractures& fractures,const unsigned int& id1, const unsigned int& id2, const double& tol){
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

Vector4d plane(const unsigned int& id, Fractures& fractures)
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

line planesIntersection(const Vector4d &coeff1, const Vector4d &coeff2, const double &tol)
{
    line output;
    output.direction = {0,0,0};
    output.point = {0,0,0};
    Vector3d v1;
    Vector3d v2;
    v1 = coeff1.head(3);
    v2 = coeff2.head(3);

    output.direction = v2.cross(v1);

    Matrix<double,2,3> M;
    M << v1[0],v1[1],v1[2],
        v2[0],v2[1],v2[2];
    Vector2d b = {-coeff1[3],-coeff2[3]};

    bool colZero = false;
    bool flag = false;
    for (unsigned int j = 0; j < 3; j++){
        if(almostEqual(M.col(j)[0],0,tol) && almostEqual(M.col(j)[1],0,tol)){
            colZero = true;
            Matrix2d A(M.rows(), M.cols() - 1);
            A << M.block(0, 0, M.rows(), j),
                M.block(0, j + 1, M.rows(), M.cols() - j - 1);
            Vector2d P = A.colPivHouseholderQr().solve(b);
            for (unsigned int k = 0; k < 3; k++){
                if (k != j){
                    output.point[k] = P[k-flag];
                }
                else{
                    output.point[k] = 0;
                    flag = 1;
                }
            }
        }
    }
    if(!colZero)
    {
        Matrix<double,2,2> M;
        M << v1[0], v1[2],
            v2[0],v2[2];
        Vector2d b = {-coeff1[3],-coeff2[3]};
        Vector2d P = M.colPivHouseholderQr().solve(b);
        output.point[0] = P[0];
        output.point[2] = P[1];
    }
    return output;
}

VectorXd linesIntersection(const line& r,const line& rj){  //primi 3 punto di inters. poi t e s
    VectorXd q;
    q.resize(5);
    q[4] = -2; //così se non c'è soluzione non salvo
    double t;
    double s;
    Matrix<double,3,2> A;
    A << r.direction[0], -rj.direction[0],
        r.direction[1], -rj.direction[1],
        r.direction[2], -rj.direction[2];

    Vector3d b = {rj.point[0] - r.point[0],
                  rj.point[1] - r.point[1],
                  rj.point[2] - r.point[2]};

    Vector2d k = A.colPivHouseholderQr().solve(b); //non è quadrata la matrice
    t = k[0];
    s = k[1];
    q << rj.point[0] + rj.direction[0]*s,
        rj.point[1] + rj.direction[1]*s,
        rj.point[2] + rj.direction[2]*s,
        t,
        s;
    return q;
}

Vector4d intervalsIntersection(const Vector4d& q,const double& tol){
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

bool almostEqual(const double& a, const double& b, const double& tol) {
    return fabs(a - b) < tol;
}

bool arePlanesParallel(const Vector4d& v1, const Vector4d& v2, double& tol) {

    Vector3d w1 = v1.head(3);
    Vector3d w2 = v2.head(3);
    Vector3d v = w1.cross(w2);

    bool par = true;

    if(abs(v[0]) > tol || abs(v[1]) > tol || abs(v[2]) > tol)
        par = false;

    return par;
}

bool compareTuple(const tuple<unsigned int, bool, double>& a, const tuple<unsigned int, bool, double>& b) {
    return get<1>(a) < get<1>(b); // Ordina in base al valore booleano all'interno delle tuple
}

void print(PolygonalMesh& mesh){

    string outputMesh = "PolygonalMesh.txt";
    ofstream ofs(outputMesh);

    if (ofs.fail())
    {
        cerr << "Error in creating the output file" << endl;
        return;
    }

    ofs << "# ID Cell0D; X; Y; Z" << endl;
    unsigned int n = 0;

    for(unsigned int i = 0; i < mesh.numberCell0Ds; i++)
    {
        ofs << n++ << "; " << mesh.coordinateCell0Ds[i][0] << "; " << mesh.coordinateCell0Ds[i][1] << "; " << mesh.coordinateCell0Ds[i][2] << endl;
    }


    ofs << endl;
    ofs << "# ID Cell1D; Origin; End" << endl;
    n = 0;

    for(unsigned int i = 0; i < mesh.numberCell1Ds; i++){
        if(mesh.isOn1D[i]){
            ofs << n++ << "; " << mesh.verticesCell1Ds[i][0] << "; " << mesh.verticesCell1Ds[i][1] << endl;
        }
    }


    ofs << endl;
    ofs << "# ID Cell2D; NumVertices; Vertices; NumEdges; Edges" << endl;
    n = 0;
    for(unsigned int i = 0; i < mesh.numberCell2Ds; i++){
        if (mesh.isOn2D[i]){
            ofs << n++ << "; " << mesh.verticesCell2Ds[i].size();
            for(unsigned int j = 0; j < mesh.verticesCell2Ds[i].size();j++)
                ofs << "; " << mesh.verticesCell2Ds[i][j];
        }
        ofs << i << "; " << mesh.edgesCell2Ds[i].size();
        for(unsigned int j = 0; j < mesh.edgesCell2Ds[i].size();j++)
            if (mesh.isOn1D[mesh.edgesCell2Ds[i][j]]){
                ofs << "; " << mesh.edgesCell2Ds[i][j];
            }
        ofs << endl;
    }
}


PolygonalMesh mergeMesh(vector<PolygonalMesh>& finalMesh){
    PolygonalMesh outputMesh;
    for (PolygonalMesh mesh : finalMesh){
        outputMesh.numberCell0Ds += mesh.numberCell0Ds;
        outputMesh.coordinateCell0Ds.insert(outputMesh.coordinateCell0Ds.end(), mesh.coordinateCell0Ds.begin(), mesh.coordinateCell0Ds.end());

        outputMesh.numberCell1Ds += mesh.numberCell1Ds;
        for (unsigned int i = 0; i < mesh.verticesCell1Ds.size(); i++){
            outputMesh.verticesCell1Ds.push_back({outputMesh.numberCell0Ds + mesh.verticesCell1Ds[i][0],outputMesh.numberCell0Ds + mesh.verticesCell1Ds[i][1]});
        }
        outputMesh.isOn1D.insert(outputMesh.isOn1D.end(), mesh.isOn1D.begin(), mesh.isOn1D.end());

        outputMesh.numberCell2Ds += mesh.numberCell2Ds;
        for (unsigned int i = 0; i < mesh.verticesCell2Ds.size(); i++){
            outputMesh.verticesCell2Ds.resize(outputMesh.verticesCell2Ds.size()+1);
            for (unsigned int j = 0; j < mesh.verticesCell2Ds[i].size(); j++){
                outputMesh.verticesCell2Ds[i][j] = outputMesh.numberCell0Ds + mesh.verticesCell2Ds[i][j];
            }
        }
        for (unsigned int i = 0; i < mesh.edgesCell2Ds.size(); i++){
            outputMesh.edgesCell2Ds.resize(outputMesh.edgesCell2Ds.size()+1);
            for (unsigned int j = 0; j < mesh.verticesCell2Ds[i].size(); j++){
                outputMesh.edgesCell2Ds[i][j] = outputMesh.numberCell1Ds + mesh.edgesCell2Ds[i][j];
            }
        }
        outputMesh.isOn2D.insert(outputMesh.isOn2D.end(), mesh.isOn2D.begin(), mesh.isOn2D.end());
    }
    return outputMesh;
}
}


