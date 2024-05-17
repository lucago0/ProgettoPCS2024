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

    Traces traces;
    unsigned int numberTraces = 0;
    for (unsigned int id1 = 0; id1<fractures.NumberFractures; id1++) {
        for (unsigned int id2 = id1+1; id2<fractures.NumberFractures; id2++) {
            if (areClose(fractures,id1,id2)){
                Vector4d coeff1 = Piano(id1,fractures);
                Vector4d coeff2 = Piano(id2,fractures);
                // piani non paralleli
                Line r = Inter(coeff1,coeff2);
                Line r_j;
                Matrix<double,4,4> intersectionPoints;
                unsigned int points = 0;
                for(unsigned int i = 0; i < 2; i++){
                    unsigned int currentId = (i==0) ? id1:id2;
                    for (unsigned int j=0; j<fractures.Vertices[currentId].rows()-1; j++){
                        r_j.point(fractures.Vertices[currentId].col(j));
                        r_j.direction(fractures.Vertices[currentId].col(j+1)-fractures.Vertices[currentId].col(j));
                        // mi assicuro che ci sia intersezione tra r ed r_j con cross
                        VectorXd Q = PuntiIntersRetta(r,r_j); // Q,t,s
                        if (Q[4]>=0 && Q[4]<=1){ // Q[4] è s!!! e tau?
                            intersectionPoints.col(points) = Q.head(4);
                            points++;
                        }
                    }
                }
                points = 0;
                numberTraces++;
                Vector4d t = intersection(intersectionPoints);
                Vector4d t_star = intersectionPoints.row(3);
                array<unsigned int,2> v = {id1,id2};
                traces.FracturesId[numberTraces] = v;
                Matrix<double,3,2> vertices;
                vertices.col(0) = r.point + t[0]*r.direction;
                vertices.col(1) = r.point + t[1]*r.direction;
                traces.Vertices.insert(make_pair(numberTraces, vertices));

                // Tips
                double a = min(t_star[0],t_star[1]);
                double b = max(t_star[0],t_star[1]);
                double c = min(t_star[2],t_star[3]);
                double d = max(t_star[2],t_star[3]);
                if (t[0] == a && t[1] == c){
                    traces.Tips[numberTraces] = true;
                }
                else if (t[0] == a && t[1] == b){
                    traces.Tips[numberTraces] = false;
                }
                else if (t[0] == b && t[1] == d){
                    traces.Tips[numberTraces] = true;
                }
                else{
                    traces.Tips[numberTraces] = false;
                }
            }
        }
    }
    traces.NumberTraces = numberTraces;

    //OutputFile(traces,fractures);

    return 0;
}
