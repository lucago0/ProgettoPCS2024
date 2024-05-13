#include <GeometryLibrary.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>

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
    cout << fracture.NumberFractures << endl;
    for (const auto& pair : fracture.Vertices) {
        std::cout << "Chiave: " << pair.first << ", Valore: " << pair.second << std::endl;
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

<<<<<<< Updated upstream
bool areClose(Fractures& mesh, unsigned int& Id1, unsigned int& Id2){

    Vector3d C1;
    const unsigned int n1 = mesh.Vertices[Id1].cols();
    Vector3d C2;
    const unsigned int n2 = mesh.Vertices[Id2].cols();

    for(unsigned int i=0; i<3; i++){
        for (unsigned int j=0; j<n1; j++){
            C1[i] += mesh.Vertices[Id1](j,i);
        }
        C1[i] /= n1;
        for (unsigned int j=0; j<n2; j++){
            C2[i] += mesh.Vertices[Id2](j,i);
        }
        C2[i] /= n2;
    }

    VectorXd rays1;
    for(unsigned int i=0; i<n1; i++){
        rays1.resize(rays1.size() + 1);
        Vector3d point = mesh.Vertices[Id1].col(i);
        rays1(rays1.size() - 1) = distanceSquared(C1,point);
    }
    VectorXd rays2;
    for(unsigned int i=0; i<n2; i++){
        rays2.resize(rays2.size() + 1);
        Vector3d point = mesh.Vertices[Id2].col(i);
        rays2(rays2.size() - 1) = distanceSquared(C2,point);
    }

    double R1 = *max_element(rays1.begin(), rays1.end());
    double R2 = *max_element(rays2.begin(), rays2.end());

    return distanceSquared(C1,C2) <= pow(R1+R2,2);
}
=======
map<unsigned int, array<double,4>> Piano(Fractures& FR)
{
    for(unsigned int i = 0; i < FR.NumberFractures; i++)
    {
        Vector3d v0 = FR.Vertices[i].col(1);
        Vector3d v1 = FR.Vertices[i].col(2);
        Vector3d v2 = FR.Vertices[i].col(3);

        map<unsigned int, array<double, 4>> coeff;

        coeff[i][0] = (v1[1]-v0[1])*(v2[2]-v0[2]) - (v1[2]-v0[2])*(v2[1]-v0[1]);
        coeff[i][1] = -((v1[0]-v0[0])*(v2[2]-v0[2])-(v2[0]-v0[0])*(v1[2]-v0[2]));
        coeff[i][2] = (v1[0]-v0[0])*(v1[1]-v0[1])-(v2[0]-v0[0])*(v1[1]-v0[1]);
        coeff[i][3] = -coeff[i][0]*v0[0]-coeff[i][1]*v0[1]-coeff[i][2]*v0[2];

        return coeff;

    }
>>>>>>> Stashed changes
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
}


