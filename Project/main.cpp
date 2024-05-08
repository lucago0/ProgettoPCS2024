#include <iostream>
#include <GeometryLibrary.hpp>
#include <Utils.hpp>

using namespace std;
using namespace FracturesLib;

int main()
{
    Fractures fracture;
    string filepath = "./FR362_data.txt";
    if(!importFracture(filepath, fracture))
    {
        return 1;
    }
    return 0;
}
