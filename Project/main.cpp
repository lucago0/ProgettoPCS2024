#include <iostream>
#include "GeometryLibrary.hpp"

using namespace std;
using namespace FracturesLib;

int main()
{
    Fracture fracture;
    string filepath = "DFN";
    if(!importFracture(filepath, fracture))
    {
        return 1;
    }
    return 0;
}
