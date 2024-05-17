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

    while (numFrac--) {
        getline(file, line);
        getline(file, line);
        stringstream ss(line);
        ss >> fracture.Id >> pv >> numVertices;

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
        fracture.Vertices.insert(make_pair(fracture.Id, actualVert));
    }
    file.close();
    return true;
}

double distanceSquared(Vector3d& A, Vector3d& B){
    return pow(A[0]-B[0],2) + pow(A[1]-B[1],2) + pow(A[2]-B[2],2);
}

void OutputFile(Traces& TR, Fractures& FR)
{
    string nameFileO = "Traces.txt";
    ofstream ofs(nameFileO);

    if (ofs.fail())
    {
        cout << "Impossibile creare il file di output" << endl;
        return;
    }

    ofs << "# Number of Traces" << endl;
    ofs << TR.FracturesId.size() << endl;
    ofs << "# TraceId; FracturesId1; FracturesId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

    for(unsigned int i = 0; i < TR.FracturesId.size();i++)
    {
        ofs << i+1 << ";" << TR.FracturesId[i][0] << ";" << TR.FracturesId[i][1] << ";" << TR.Vertices[i][0] << ";" << TR.Vertices[i][1] << ";" << TR.Vertices[i][2] << endl;
    }

    map<unsigned int, unsigned int> FracTrace;

    for(unsigned int i = 0; i < FR.NumberFractures; i++)
    {
        for(unsigned int j = 0; j < TR.FracturesId.size(); j++)
        {
            if(i == TR.FracturesId[j][0] || i == TR.FracturesId[j][1])
            {
                FracTrace[i] += 1;
            }
        }
    }

    ofs << "# FractureId; NumTraces" << endl;
    for(unsigned int i = 0; i < FR.NumberFractures; i++)
    {
        ofs << i << ";" << FracTrace[i] << endl;
    }

    ofs << "# TraceId; Tips; Length" << endl;
    for(unsigned int i = 0; i < TR.FracturesId.size();i++)
    {
        ofs << i << ";" << TR.Tips << sqrt(distanceSquared(TR.Vertices[0],TR.Vertices[1])) << endl;
    }

}

bool areClose(Fractures& fracture, unsigned int& Id1, unsigned int& Id2){
    Vector3d C1;
    const unsigned int n1 = fracture.Vertices[Id1].cols();
    Vector3d C2;
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

    return distanceSquared(C1,C2) <= pow(R1+R2,2);
}

Vector4d Piano(unsigned int& id, Fractures& FR)
{

    Vector3d v0 = FR.Vertices[id].col(1);
    Vector3d v1 = FR.Vertices[id].col(2);
    Vector3d v2 = FR.Vertices[id].col(3);

    array<double, 4> coeff = {0,0,0,0};

    coeff[0] = (v1[1]-v0[1])*(v2[2]-v0[2]) - (v1[2]-v0[2])*(v2[1]-v0[1]);
    coeff[1] = -((v1[0]-v0[0])*(v2[2]-v0[2])-(v2[0]-v0[0])*(v1[2]-v0[2]));
    coeff[2] = (v1[0]-v0[0])*(v1[1]-v0[1])-(v2[0]-v0[0])*(v1[1]-v0[1]);
    coeff[3] = -coeff[0]*v0[0]-coeff[1]*v0[1]-coeff[2]*v0[2];

    return coeff;
}

Line Inter(Vector4d& coeff1, Vector4d& coeff2)
{
    Line line;
    Vector3d v1;
    Vector3d v2;
    for(unsigned int i = 0; i < 3; i++)
    {v1[i] = coeff1[i];
        v2[i] = coeff2[i];}

    line.direction[0] = v1[1]*v2[2] - v1[2]*v2[1];
    line.direction[1] = v1[2]*v2[0] - v1[0]*v2[2];
    line.direction[2] = v1[0]*v2[1] - v1[1]*v2[0];

    if(v2[2] != 0 || v1[2] != 0)
    {
        Matrix<double,2,2> M;
        M << v1[0], v1[1],
            v2[0],v2[1];
        Vector2d b = {coeff1[3],coeff2[3]};
        Vector2d P = M.lu().solve(b);
        line.point[0] = P[0];
        line.point[1] = P[1];
    }
    else if (v2[1] != 0 || v1[1] != 0)
    {
        Matrix<double,2,2> M;
        M << v1[0], v1[2],
            v2[0],v2[2];
        Vector2d b = {coeff1[3],coeff2[3]};
        Vector2d P = M.lu().solve(b);
        line.point[0] = P[0];
        line.point[2] = P[1];
    }

    /*else{
     * ...inserire altri casi...
        } */

    return line;
}

VectorXd PuntiIntersRetta(Line& r, Line& rj){  //primi 3 punto di inters. poi t e s
    VectorXd Q;
    Q.resize(5);
    double t;
    double s;
    Vector4d Qtemp;
    int counter = 0;
    Matrix<double,4,4> Q; //matrice dei punti di intersezione

    for(unsigned int i = 0; i < 2; i++) {
        unsigned int currentId = (i == 0) ? Id1 : Id2;
        for(unsigned int j = 0; j < fracture.Vertices[currentId].cols()-1; j++){
            Matrix<double,3,2> A;
            A << v[0], (fracture.Vertices[currentId](0,j+1)-fracture.Vertices[currentId](0,j)),
                v[1], (fracture.Vertices[currentId](1,j+1)-fracture.Vertices[currentId](1,j)),
                v[2], (fracture.Vertices[currentId](2,j+1)-fracture.Vertices[currentId](2,j));

            Vector3d b = {fracture.Vertices[currentId](0,j) - v[3],
                          fracture.Vertices[currentId](1,j) - v[4],
                          fracture.Vertices[currentId](2,j) - v[5]};

            FullPivLU<MatrixXd> lu_decomp(A);
            int rank_A = lu_decomp.rank();

            // Creare la matrice aumentata [A | B]
            MatrixXd augmented(A.rows(), A.cols() + 1);
            augmented << A, b;

            // Calcolare il rango della matrice aumentata [A | B] usando la decomposizione LU con pivotaggio completo
            FullPivLU<MatrixXd> lu_decomp_aug(augmented);
            int rank_augmented = lu_decomp_aug.rank();

            if(rank_A == rank_augmented){
                Vector2d r = A.colPivHouseholderQr().solve(b); //non è quadrata la matrice
                t = r[0];
                s = r[1];
                Q = (fracture.Vertices[currentId](0,j) + A(0,1)*s,
                     fracture.Vertices[currentId](1,j) + A(1,1)*s,
                     fracture.Vertices[currentId](2,j) + A(2,1)*s,
                     s
                         t);

            };
        };
    }
}

return Q;
}

Vector4d intersection(const Matrix<double,4,4>& Q){
    double a = Q(3,0);
    double b = Q(3,1);
    double c = Q(3,2);
    double d = Q(3,3);

    Vector4d v = {a,b,c,d};
    // Calcola l'estremo sinistro dell'intersezione
    double sx = max(a, c);
    // Calcola l'estremo destro dell'intersezione
    double dx = min(b, d);
    // Se gli intervalli non si sovrappongono, l'intersezione sarà vuota
    if (sx > dx) {
        sx = dx = numeric_limits<double>::quiet_NaN(); // Non un numero
    }
    double other_sx = (a < c) ? a : c; // other_sx è pari ad a se a<c, altrimenti è pari a c
    double other_dx = (d > b) ? d : b; // other_dx è pari a d se d>b, altrimenti è pari a b
    Vector4d output = {sx,dx,other_sx,other_dx}; // in ordine restituiamo l'intervallo di intersezione e gli altri due estremi ordinati
    return output;
}


};
