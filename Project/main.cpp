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

    Traces traces;
    unsigned int numberTraces = 0;
    bool flag;
    for (unsigned int id1 = 0; id1<fractures.NumberFractures; id1++) {
        for (unsigned int id2 = id1+1; id2<fractures.NumberFractures; id2++) {
            if (areClose(fractures,id1,id2)){
                flag = false;
                Vector4d coeff1 = Piano(id1,fractures);
                Vector4d coeff2 = Piano(id2,fractures);
                if(!arePlanesParallel(coeff1[0],coeff1[1],coeff1[2],coeff2[0],coeff2[1],coeff2[2],pow(10,-10))){
                    Line r = Inter(coeff1,coeff2);
                    Line r_j;
                    Matrix<double,4,4> intersectionPoints;
                    unsigned int points = 0;
                    for(unsigned int i = 0; i < 2; i++){
                        unsigned int currentId = (i==0) ? id1:id2;
                        for (unsigned int j=0; j<fractures.Vertices[currentId].cols(); j++){
                            if(j<fractures.Vertices[currentId].cols()-1){
                                r_j.point = fractures.Vertices[currentId].col(j);
                                r_j.direction = (fractures.Vertices[currentId].col(j+1)-fractures.Vertices[currentId].col(j));
                            }
                            else{
                                r_j.point = fractures.Vertices[currentId].col(j);
                                r_j.direction = (fractures.Vertices[currentId].col(0)-fractures.Vertices[currentId].col(j));
                            }
                            // mi assicuro che ci sia intersezione tra r ed r_j con cross
                            Vector3d test = (r.direction).cross(r_j.direction);
                            if(!almostEqual(test[0],0,pow(10,-10)) || !almostEqual(test[1],0,pow(10,-10)) || !almostEqual(test[2],0,pow(10,-10))){
                                VectorXd Q = PuntiIntersRetta(r,r_j); // Q,t,s
                                if (Q[4]>= 0 && Q[4]<=1){ // Q[4] è s!!! e tau?
                                    intersectionPoints.col(points) = Q.head(4);
                                    points++;
                                    flag = true;
                                };
                            };
                        }
                    }
                    if(flag && points == 4){
                        Vector4d t_star = intersectionPoints.row(3);
                        Vector4d t = intersection(t_star);
                        if(!isnan(t[0])){
                            array<unsigned int,2> v = {id1,id2};
                            traces.FracturesId[numberTraces] = v;
                            Matrix<double,3,2> vertices;
                            vertices.col(0) = r.point + t[0]*r.direction;
                            vertices.col(1) = r.point + t[1]*r.direction;
                            traces.Vertices.insert(make_pair(numberTraces, vertices));
                            numberTraces++;

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
                        };
                    };
                }
            }
        }
    }
    traces.NumberTraces = numberTraces;

    for (unsigned int id = 0; id < traces.NumberTraces; id++){
        traces.Lengths[id] = sqrt(distanceSquared(traces.Vertices[id].col(0),traces.Vertices[id].col(1)));
    }

    OutputFile(traces,fractures);

    return 0;
}
