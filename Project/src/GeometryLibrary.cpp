#include <GeometryLibrary.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

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


Matrix<unsigned int, Dynamic,2> areClose(Fractures& frac){

    Matrix<unsigned int, Dynamic,2> Temp;

    for (auto it1 = frac.Vertices.begin(); it1 != prev(frac.Vertices.end()); ++it1) {
        for (auto it2 = next(it1); it2 != frac.Vertices.end(); ++it2) {
            unsigned int iD1 = it1->first;
            unsigned int iD2 = it2->first;
            Vector3d C1;
            const unsigned int n1 = frac.Vertices[Id1].cols();
            Vector3d C2;
            const unsigned int n2 = frac.Vertices[Id2].cols();

            for(unsigned int i=0; i<3; i++){
                for (unsigned int j=0; j<n1; j++){
                    C1[i] += frac.Vertices[Id1](j,i);
                }
                C1[i] /= n1;
                for (unsigned int j=0; j<n2; j++){
                    C2[i] += frac.Vertices[Id2](j,i);
                }
                C2[i] /= n2;
            }

            VectorXd rays1;
            for(unsigned int i=0; i<n1; i++){
                rays1.resize(rays1.size() + 1);
                Vector3d point = frac.Vertices[Id1].col(i);
                rays1(rays1.size() - 1) = distanceSquared(C1,point);
            }
            VectorXd rays2;
            for(unsigned int i=0; i<n2; i++){
                rays2.resize(rays2.size() + 1);
                Vector3d point = frac.Vertices[Id2].col(i);
                rays2(rays2.size() - 1) = distanceSquared(C2,point);
            }

            double R1 = *max_element(rays1.begin(), rays1.end());
            double R2 = *max_element(rays2.begin(), rays2.end());

            if distanceSquared(C1,C2) <= pow(R1+R2,2){
                    Temp = Temp.append(Vector2d(Id1,Id2));
                };
        }
        }
    return Temp;
}

map<unsigned int, array<double,4>> Piano(Fractures& FR)
{
    map<unsigned int, array<double, 4>> coeff;

    for(unsigned int i = 0; i < FR.NumberFractures; i++)
    {
        Vector3d v0 = FR.Vertices[i].col(1);
        Vector3d v1 = FR.Vertices[i].col(2);
        Vector3d v2 = FR.Vertices[i].col(3);

        coeff[i][0] = (v1[1]-v0[1])*(v2[2]-v0[2]) - (v1[2]-v0[2])*(v2[1]-v0[1]);
        coeff[i][1] = -((v1[0]-v0[0])*(v2[2]-v0[2])-(v2[0]-v0[0])*(v1[2]-v0[2]));
        coeff[i][2] = (v1[0]-v0[0])*(v1[1]-v0[1])-(v2[0]-v0[0])*(v1[1]-v0[1]);
        coeff[i][3] = -coeff[i][0]*v0[0]-coeff[i][1]*v0[1]-coeff[i][2]*v0[2];
        return coeff;
    }
}

array<double,6> Inter(Fractures& FR, unsigned int& id1, unsigned int& id2)
{
     Vector3d v1;
     Vector3d v2;
     for(unsigned int i = 0; i < 3; i++)
     {v1[i] = FR.Coeff[id1][i];
     v2[i] = FR.Coeff[id2][i];}

     array<double,6> vect;

     vect[0] = v1[1]*v2[2] - v1[2]*v2[1];
     vect[1] = v1[2]*v2[0] - v1[0]*v2[2];
     vect[2] = v1[0]*v2[1] - v1[1]*v2[0];

     if(v2[2] != 0 && v1[2] != 0)
     {
         Matrix<double,2,2> M;
         M << v1[0], v1[1],
                 v2[0],v2[1];
         Vector2d b = {FR.Coeff[id1][3],FR.Coeff[id2][3]};
         Vector2d P = M.lu().solve(b);
         vect[3] = P[0];
         vect[4] = P[1];
     }
}


//FRClose matrice con gli Id delle fratture che possono intersecarsi, restituita da areClose
Matrix<double,4,4> PuntoIntersRetta(const struct Fractures& fracture, Matrix<unsigned int, dynamic, 2>& FrClose){
    double t;
    double s;
    Vector3d Qtemp;
    Matrix<double,4,4> Q; //matrice dei punti di intersezione

    unsigned int NumRow = FRClose.rows();

    for (int i = 0; i < NumRow; ++i) {
        unsigned int Id1 = FRClose(i, 0);
        unsigned int Id2 = FRClose(i, 1);
        for(unsigned int k = 0; k < 2; k++) {
            unsigned int currentId = (k == 0) ? Id1 : Id2;
            for (unsigned int j = 0; j < fracture.Vertices[i].cols(); j++){

                Matrix<double,3,2> A;
                A << vx, (fracture.Vertices[currentId](0,j+1)-fracture.Vertices[currentId](0,j)),
                     vy, (fracture.Vertices[currentId](1,j+1)-fracture.Vertices[currentId](1,j)),
                     vz, (fracture.Vertices[currentId](2,j+1)-fracture.Vertices[currentId](2,j));

                Vector3d b = {fracture.Vertices[currentId](0,j) - xbar,
                              fracture.Vertices[currentId](1,j) - ybar,
                              fracture.Vertices[currentId](2,j) - zbar};

                bool soluzione_esiste = (A.fullPivLu().rank() == A.fullPivLu().append(b).rank()); //verifico esistenza soluzione

                if(soluzione_esiste){
                    Vector2d r = A.colPivHouseholderQr().solve(b); //non è quadrata la matrice
                    t = r[0];
                    s = r[1];
                    if(s >= 0 && s <= 1){
                        Qtemp = {
                            fracture.Vertices[currentId](0,j) + A(0,1)*s,
                            fracture.Vertices[currentId](1,j) + A(1,1)*s,
                            fracture.Vertices[currentId](2,j) + A(2,1)*s,
                            t
                        };
                        Q = Q.append(Qtemp);
                    };
                };
            }
        }
    }
    return Q;
}

array<double,4> intersection(const Matrix&Q){
    a = Q(3,0);
    b = Q(3,1);
    c = Q(3,2);
    d = Q(3,3);

    vector<double> v = {a,b,c,d};
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
    array<double,4> output = {sx,dx,other_sx,other_dx}; // in ordine restituiamo l'intervallo di intersezione e gli altri due estremi ordinati
    return output;
}


};

