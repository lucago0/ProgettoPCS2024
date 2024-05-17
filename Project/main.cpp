#include <iostream>
#include <GeometryLibrary.hpp>
#include <Utils.hpp>

using namespace std;
using namespace FracturesLib;

int main()
{
    Fractures fractures;
    string filepath = "./FR3_data.txt";
    if(!importFracture(filepath, fractures))
    {
        return 1;
    }

<<<<<<< Updated upstream
    for (unsigned int i=0;i<fractures.NumberFractures;i++){
        array<double,4> coeff = Piano(i, fractures);
        fractures.Coeff.insert(make_pair(fractures.Id, coeff));
=======

    for (unsigned int id1 = 0; i<fractures.NumberFractures; id1++) {
        for (unsigned int id2 = id1+1; i<fractures.NumberFractures; id2++) {
            if (areClose(fractures,id1,id2)){
                array<double,4> coeff1 = Piano(id1,fractures);
                array<double,4> coeff2 = Piano(id2,fractures);
                Line r = Inter(coeff1,coeff2);

            }

        }
>>>>>>> Stashed changes
    }

    for (unsigned int id1 = 0; id1<fractures.NumberFractures; id1++) {
        for (unsigned int id2 = id1+1; id2<fractures.NumberFractures; id2++) {
            if (areClose(fractures,id1,id2)){
                array<double,6> v = Inter(fractures.Coeff[id1],fractures.Coeff[id2]); //trova i coefficienti della retta di intersezione
                                                                                      //in forma parametrica (vx,vy,vz, xbar, ybar, zbar)
                Matrix<double, 4,4> Q = PuntiIntersRetta(fractures,id1,id2,v); //trova i punti di intersezione
                array<double,4> inters = intersection(Q); //restuisce i punti ordinati
                };
            };
        }

    /* Per la funzione che associa a una retta r e a uno
     * scalare t il valore P+tv (dove P è il punto iniziale
     * e v la sua direzione) scriviamo
     * array<double,3> P = r.point + t*r.direction */

    return 0;
}
