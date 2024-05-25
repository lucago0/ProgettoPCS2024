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

namespace FracturesLib{

bool importFracture(const string& filename, Fractures& fracture) {
    ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    string header;
    getline(file, header); // Read the header line to discard it

    string line;
    getline(file,line);
    fracture.NumberFractures = stoi(line);
    int numFrac = stoi(line);
    char pv;
    int numVertices;
    unsigned int k;

    while (numFrac--) {
        getline(file, line);
        getline(file, line);
        stringstream ss(line);
        ss >> k >> pv >> numVertices;

        Matrix<double, 3, Dynamic> actualVert(3, numVertices);

        getline(file, line);
        string val;
        // Read Vertices
        for(int j = 0; j < 3; j++){
            for(int i = 0; i < numVertices; i++){
                file >> val;
                actualVert(j, i) = stod(val);
            }
        }
        file >> pv;

        // Add fracture data to the fractures map
        fracture.Vertices.insert(make_pair(k, actualVert));
    }
    file.close();
    return true;
}

double distanceSquared(const Vector3d& A,const Vector3d& B){
    return pow(A[0]-B[0],2) + pow(A[1]-B[1],2) + pow(A[2]-B[2],2);
}

bool compareByValue(const pair<unsigned int, double> &a, const pair<unsigned int, double> &b) {
    return a.second > b.second;
}

void OutputFile(Traces& TR, Fractures& FR)
{
    string Tracce = "Traces.txt";
    string FrattTracc = "Fratture-Tracce.txt";
    ofstream ofs(Tracce);
    ofstream ofs2(FrattTracc);

    if (ofs.fail())
    {
        cout << "Impossibile creare il file di output" << endl;
        return;
    }

    ofs << "# Number of Traces" << endl;
    ofs << TR.NumberTraces << endl;
    ofs << "# TraceId; FracturesId1; FracturesId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

    for(unsigned int i = 0; i < TR.NumberTraces;i++)
    {
        ofs << i << ";" << TR.FracturesId[i][0] << ";" << TR.FracturesId[i][1] << ";" << TR.Vertices[i](0,0) << ";" << TR.Vertices[i](1,0) << ";" << TR.Vertices[i](2,0) << ";" << TR.Vertices[i](0,1) << ";" << TR.Vertices[i](1,1) << ";" << TR.Vertices[i](2,1) << endl;
    }



    map<unsigned int, unsigned int> FracTrace;

    for(unsigned int i = 0; i < FR.NumberFractures; i++)
    {
        for(unsigned int j = 0; j < TR.NumberTraces; j++)
        {
            if(i == TR.FracturesId[j][0] || i == TR.FracturesId[j][1])
            {
                FracTrace[i] += 1;
            }
        }
    }

    // Copia gli elementi della mappa in un vettore di coppie
    vector<pair<unsigned int, double>> mapElements(TR.Lengths.begin(), TR.Lengths.end());

    // Ordina il vettore in base ai valori
    sort(mapElements.begin(), mapElements.end(), compareByValue);

    for(unsigned int i = 0; i < FR.NumberFractures; i++)
    {
        ofs2 << "# FractureId; NumTraces" << endl;
        ofs2 << i << ";" << FracTrace[i] << endl;
        ofs2 << "# TraceId; Tips; Length" << endl;
        for(auto& couple : mapElements){
            if(TR.FracturesId[couple.first][0] == i){
                ofs2 << couple.first << ";" << TR.Tips[couple.first][0] << ";" << couple.second << endl;
            }
            else if(TR.FracturesId[couple.first][1] == i){
                ofs2 << couple.first << ";" << TR.Tips[couple.first][1] << ";" << couple.second << endl;
            }
        }
    }
}

bool areClose(Fractures& fracture, unsigned int& Id1, unsigned int& Id2){
    Vector3d C1 = {0,0,0};
    const unsigned int n1 = fracture.Vertices[Id1].cols();
    Vector3d C2 = {0,0,0};
    const unsigned int n2 = fracture.Vertices[Id2].cols();

    for(unsigned int i=0; i<3; i++){
        for (unsigned int j=0; j<n1; j++){
            C1[i] += fracture.Vertices[Id1](i,j);
        }
        C1[i] /= n1;
        for (unsigned int j=0; j<n2; j++){
            C2[i] += fracture.Vertices[Id2](i,j);
        }
        C2[i] /= n2;
    }

    VectorXd rays1;
    for(unsigned int i=0; i<n1; i++){
        rays1.resize(rays1.size() + 1);
        Vector3d point = fracture.Vertices[Id1].col(i);
        rays1(rays1.size() - 1) = distanceSquared(C1,point);
    }
    VectorXd rays2;
    for(unsigned int i=0; i<n2; i++){
        rays2.resize(rays2.size() + 1);
        Vector3d point = fracture.Vertices[Id2].col(i);
        rays2(rays2.size() - 1) = distanceSquared(C2,point);
    }

    double R1 = *max_element(rays1.begin(), rays1.end());
    double R2 = *max_element(rays2.begin(), rays2.end());

    return distanceSquared(C1,C2) <= (R1+R2+(2*sqrt(R1)*sqrt(R2))); //capire bene
}

Vector4d Piano(unsigned int& id, Fractures& FR)
{

    Vector3d v0 = FR.Vertices[id].col(0);
    Vector3d v1 = FR.Vertices[id].col(1);
    Vector3d v2 = FR.Vertices[id].col(2);

    Vector4d coeff = {0,0,0,0};

    coeff[0] = (v1[1]-v0[1])*(v2[2]-v0[2]) - (v1[2]-v0[2])*(v2[1]-v0[1]);
    coeff[1] = -((v1[0]-v0[0])*(v2[2]-v0[2])-(v2[0]-v0[0])*(v1[2]-v0[2]));
    coeff[2] = (v1[0]-v0[0])*(v2[1]-v0[1])-(v2[0]-v0[0])*(v1[1]-v0[1]);
    coeff[3] = -(coeff[0]*v0[0]+coeff[1]*v0[1]+coeff[2]*v0[2]); //da capire i segni

    return coeff;
}

Line Inter(const Vector4d& coeff1, const Vector4d& coeff2)
{
    Line line;
    Vector3d v1;
    Vector3d v2;
    v1 = coeff1.head(3);
    v2 = coeff2.head(3);

    line.direction = v2.cross(v1);

    if(line.direction[2] != 0)
    {
        Matrix<double,2,2> M;
        M << v1[0], v1[1],
            v2[0],v2[1];
        Vector2d b = {-coeff1[3],-coeff2[3]};
        Vector2d P = M.colPivHouseholderQr().solve(b);
        line.point[0] = P[0];
        line.point[1] = P[1];
    }
    else
    {
        Matrix<double,2,2> M;
        M << v1[0], v1[2],
            v2[0],v2[2];
        Vector2d b = {-coeff1[3],-coeff2[3]};
        Vector2d P = M.colPivHouseholderQr().solve(b);
        line.point[0] = P[0];
        line.point[2] = P[1];
    }

    /*else{
     * ...inserire altri casi...
        } */

    return line;
}

VectorXd PuntiIntersRetta(const Line& r,const Line& rj){  //primi 3 punto di inters. poi t e s
    VectorXd Q;
    Q.resize(5);
    Q[4] = -2; //così se non c'è soluzione non salvo
    double t;
    double s;
    Matrix<double,3,2> A;
    A << r.direction[0], rj.direction[0],
        r.direction[1], rj.direction[1],
        r.direction[2], rj.direction[2];

    Vector3d b = {rj.point[0] - r.point[0],
                  rj.point[1] - r.point[1],
                  rj.point[2] - r.point[2]};

    // FullPivLU<MatrixXd> lu_decomp(A);
    // int rank_A = lu_decomp.rank();

    // // Creare la matrice aumentata [A | B]
    // MatrixXd augmented(A.rows(), A.cols() + 1);
    // augmented << A, b;

    // // Calcolare il rango della matrice aumentata [A | B] usando la decomposizione LU con pivotaggio completo
    // FullPivLU<MatrixXd> lu_decomp_aug(augmented);
    // int rank_augmented = lu_decomp_aug.rank();

    //if(rank_A == rank_augmented){
    Vector2d k = A.colPivHouseholderQr().solve(b); //non è quadrata la matrice
    t = k[0];
    s = -k[1]; //da capire perchè ci va il -
    Q << rj.point[0] + rj.direction[0]*s,
        rj.point[1] + rj.direction[1]*s,
        rj.point[2] + rj.direction[2]*s,
        t,
        s;
    //};
    return Q;
}

Vector4d intersection(const Vector4d& Q){
    double a = min(Q[0],Q[1]);
    double b = max(Q[0],Q[1]);
    double c = min(Q[2],Q[3]);
    double d = max(Q[2],Q[3]);

    // Calcola l'estremo sinistro dell'intersezione
    double sx = max(a, c);
    // Calcola l'estremo destro dell'intersezione
    double dx = min(b, d);
    //Se gli intervalli non si sovrappongono, l'intersezione sarà vuota
    if (sx > dx) {
        sx = dx = numeric_limits<double>::quiet_NaN(); // Non un numero
    }
    double other_sx = (a < c) ? a : c; // other_sx è pari ad a se a<c, altrimenti è pari a c
    double other_dx = (d > b) ? d : b; // other_dx è pari a d se d>b, altrimenti è pari a b
    Vector4d output = {sx,dx,other_sx,other_dx}; // in ordine restituiamo l'intervallo di intersezione e gli altri due estremi ordinati
    return output;
}

bool almostEqual(double a, double b, double tol) {
    return fabs(a - b) < tol;
}

bool arePlanesParallel(double A1, double B1, double C1, double A2, double B2, double C2, double tol) {
    // Compute the normal vectors of the planes
    double normal1_length = sqrt(A1 * A1 + B1 * B1 + C1 * C1);
    double normal2_length = sqrt(A2 * A2 + B2 * B2 + C2 * C2);

    // Check if the normal vectors are parallel (i.e., dot product is 1 or -1)
    double dot_product = (A1 * A2 + B1 * B2 + C1 * C2) / (normal1_length * normal2_length);

    // Check if the dot product is close to 1 or -1 within tolerance
    return almostEqual(dot_product, 1.0, tol) || almostEqual(dot_product, -1.0, tol);
}


};
