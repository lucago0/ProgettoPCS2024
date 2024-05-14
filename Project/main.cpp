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


    /* Per la funzione che associa a una retta r e a uno
     * scalare t il valore P+tv (dove P è il punto iniziale
     * e v la sua direzione) scriviamo
     * array<double,3> P = r.point + t*r.direction */

    return 0;


}
