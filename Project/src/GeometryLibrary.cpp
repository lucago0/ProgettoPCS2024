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
        ofs << n++ << "; " << scientific << setprecision(8) << mesh.coordinateCell0Ds[i][0] << "; " << mesh.coordinateCell0Ds[i][1] << "; " << mesh.coordinateCell0Ds[i][2] << endl;
    }


    ofs << endl;
    ofs << "# ID Cell1D; Origin; End" << endl;
    n = 0;

    for(unsigned int i = 0; i < mesh.numberCell1Ds; i++){
        if(mesh.isOn1D[i]){
            ofs << fixed << setprecision(0) << i << "; " << mesh.verticesCell1Ds[i][0] << "; " << mesh.verticesCell1Ds[i][1] << endl;
        }
    }


    ofs << endl;
    ofs << "# ID Cell2D; NumVertices; Vertices; NumEdges; Edges" << endl;
    n = 0;
    for(unsigned int i = 0; i < mesh.numberCell2Ds; i++){
        if (mesh.isOn2D[i]){
            ofs << i << "; " << mesh.verticesCell2Ds[i].size();
            for(unsigned int j = 0; j < mesh.verticesCell2Ds[i].size();j++){
                ofs  << "; " << mesh.verticesCell2Ds[i][j];
            }

            ofs << "; " << mesh.edgesCell2Ds[i].size();
            for(unsigned int j = 0; j < mesh.edgesCell2Ds[i].size();j++){
                if (mesh.isOn1D[mesh.edgesCell2Ds[i][j]]){
                    ofs << "; " << mesh.edgesCell2Ds[i][j];
                }
            }
            ofs << endl;
        }
    }
}


PolygonalMesh mergeMesh(vector<PolygonalMesh>& finalMesh, const unsigned int& NCell0D, const unsigned int& NCell1D, const unsigned int& NCell2D){
    PolygonalMesh outputMesh;
    outputMesh.coordinateCell0Ds.resize(NCell0D);
    outputMesh.verticesCell1Ds.resize(NCell1D);
    outputMesh.isOn1D.resize(NCell1D);
    outputMesh.verticesCell2Ds.resize(NCell2D);
    outputMesh.edgesCell2Ds.resize(NCell2D);
    outputMesh.isOn2D.resize(NCell2D);

    for (PolygonalMesh mesh : finalMesh){

        for (unsigned int i = 0; i < mesh.numberCell0Ds; i++){
            outputMesh.coordinateCell0Ds[outputMesh.numberCell0Ds + i] = mesh.coordinateCell0Ds[i];
        }

        for (unsigned int i = 0; i < mesh.numberCell1Ds; i++){
            outputMesh.verticesCell1Ds[i+outputMesh.numberCell1Ds][0] = outputMesh.numberCell0Ds + mesh.verticesCell1Ds[i][0];
            outputMesh.verticesCell1Ds[i+outputMesh.numberCell1Ds][1] = outputMesh.numberCell0Ds + mesh.verticesCell1Ds[i][1];
            outputMesh.isOn1D[i+outputMesh.numberCell1Ds] = mesh.isOn1D[i];
        }

        for (unsigned int i = 0; i < mesh.numberCell2Ds; i++){
            outputMesh.verticesCell2Ds[i+outputMesh.numberCell2Ds].resize(mesh.verticesCell2Ds[i].size());
            outputMesh.edgesCell2Ds[i+outputMesh.numberCell2Ds].resize(mesh.edgesCell2Ds[i].size());
            for (unsigned int j = 0; j < mesh.verticesCell2Ds[i].size(); j++){
                outputMesh.verticesCell2Ds[i+outputMesh.numberCell2Ds][j] = (outputMesh.numberCell0Ds + mesh.verticesCell2Ds[i][j]);
                outputMesh.edgesCell2Ds[i+outputMesh.numberCell2Ds][j] = (outputMesh.numberCell1Ds + mesh.edgesCell2Ds[i][j]);
            }
            outputMesh.isOn2D[i+outputMesh.numberCell2Ds] = mesh.isOn2D[i];
        }

        outputMesh.numberCell0Ds += mesh.numberCell0Ds;
        outputMesh.numberCell1Ds += mesh.numberCell1Ds;
        outputMesh.numberCell2Ds += mesh.numberCell2Ds;
    }
    return outputMesh;
}
}


