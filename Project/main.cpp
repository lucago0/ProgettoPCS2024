#include <iostream>
#include <GeometryLibrary.hpp>
#include <Utils.hpp>

using namespace std;
using namespace FracturesLib;

int main()
{
    Fractures fractures;
    string filepath = "./FR362_data.txt";
    if(!importFracture(filepath, fractures))
    {
        return 1;
    }


    for (unsigned int id1 = 0; id1<fractures.NumberFractures; id1++) {
        for (unsigned int id2 = id1+1; id2<fractures.NumberFractures; id2++) {
            if areClose(fractures,id1,id2){

            }

        }
    }



    /* Per la funzione che associa a una retta r e a uno
     * scalare t il valore P+tv (dove P è il punto iniziale
     * e v la sua direzione) scriviamo
     * array<double,3> P = r.point + t*r.direction */

    return 0;


}
