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
        file >> pv; //per gestire la fine con lo spazio vuoto

        // Add fracture data to the fractures map
        fracture.Vertices.insert(make_pair(fracture.Id, actualVert));
    }
    cout << fracture.NumberFractures << endl;
    file.close();
    return true;
}
}
