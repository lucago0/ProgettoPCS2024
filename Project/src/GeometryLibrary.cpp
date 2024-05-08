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

double distanceSquared(Vector3d& A, Vector3d& B){
    return pow(A[0]-B[0],2) + pow(A[1]-B[1],2) + pow(A[2]-B[2],2);
}

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
}
