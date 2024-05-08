#include <GeometryLibrary.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <array>
#include <vector>

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

}
