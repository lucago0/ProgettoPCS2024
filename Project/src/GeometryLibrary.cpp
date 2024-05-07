#include <GeometryLibrary.hpp>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

namespace FracturesLib{

bool importFracture(const string& filename, Fractures& fracture){

    ifstream file;
    file.open(filename);
    if(file.fail()){
        return false;
    }

    ifstream reading;

    while(file >> reading){
        string header;
        getline(file,header);

        string line;
        getline(file,line);

        fracture.NumberFractures = stoi(line);

        string header;
        getline(file,header);

        int numVertices;
        char pv;

        file >> fracture.Id >> pv >> numVertices;

        vector<array<double, 3>> vertices;

        for(j=0; j<3; j++)
        {
            for(i=0; i<numVertices; i++)
            {
                file >> vertices[i][j] >> pv;
            }
        }

        fracture.Vertices = vertices;
    }

}
