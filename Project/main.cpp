#include <iostream>
#include <GeometryLibrary.hpp>
#include <Utils.hpp>

using namespace std;
using namespace FracturesLib;

int main()
{
    Fractures fracture;
    string filepath = "/Users/andreacigna/Desktop/PCS24/ProgettoPCS2024/Project/DFN/FR10_data.txt";
    if(!importFracture(filepath, fracture))
    {
        return 1;
    }
    return 0;
}
