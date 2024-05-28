#include <iostream>
#include <GeometryLibrary.hpp>
#include <Utils.hpp>

using namespace std;
using namespace FracturesLib;

int main()
{
    Fractures fractures;
    string filepath = "./FR10_data.txt";
    if(!importFracture(filepath, fractures))
    {
        return 1;
    }

    fractures.NumTracce.resize(fractures.NumberFractures);
    Traces traces;
    unsigned int numberTraces = 0;
    for (unsigned int id1 = 0; id1<fractures.NumberFractures; id1++) {
        for (unsigned int id2 = id1+1; id2<fractures.NumberFractures; id2++) {
            if (areClose(fractures,id1,id2)){
                Vector4d coeff1 = Piano(id1,fractures);
                Vector4d coeff2 = Piano(id2,fractures);
                if(!arePlanesParallel(coeff1,coeff2,pow(10,-10))){
                    Line r = Inter(coeff1,coeff2);
                    Line r_j;
                    Matrix<double,4,4> intersectionPoints;
                    unsigned int points = 0;
                    for(unsigned int i = 0; i < 2; i++){
                        unsigned int currentId = (i==0) ? id1:id2;
                        unsigned int n = fractures.Vertices[currentId].cols();
                        for (unsigned int j=0; j<n; j++){
                            r_j.point = fractures.Vertices[currentId].col(j);
                            r_j.direction = fractures.Vertices[currentId].col((j+1)%n)-fractures.Vertices[currentId].col(j);
                            // mi assicuro che ci sia intersezione tra r ed r_j con cross
                            Vector3d test = (r.direction).cross(r_j.direction);
                            if(!almostEqual(test[0],0,pow(10,-10)) || !almostEqual(test[1],0,pow(10,-10)) || !almostEqual(test[2],0,pow(10,-10))){
                                VectorXd Q = PuntiIntersRetta(r,r_j); // Q,t,s
                                if (Q[4]>= 0 && Q[4]<=1){ // Q[4] è s!!! e tau?
                                    intersectionPoints.col(points) = Q.head(4);
                                    points++;
                                };
                            };
                        }
                    }
                    if(points == 4){
                        Vector4d t_star = intersectionPoints.row(3);
                        Vector4d t = intersection(t_star);
                        if(!isnan(t[0])){
                            array<unsigned int,4> v = {id1,id2,0,0};
                            traces.FracturesId.push_back(v);
                            Matrix<double,3,2> vertices;
                            vertices.col(0) = r.point + t[0]*r.direction;
                            vertices.col(1) = r.point + t[1]*r.direction;
                            traces.Lengths[numberTraces] = sqrt(distanceSquared(vertices.col(0),vertices.col(1)));
                            traces.Vertices.insert(make_pair(numberTraces, vertices));

                            // Tips
                            double a = min(t_star[0],t_star[1]);
                            double b = max(t_star[0],t_star[1]);
                            double c = min(t_star[2],t_star[3]);
                            double d = max(t_star[2],t_star[3]);
                            if (t[0] == a && t[1] == b){
                                traces.FracturesId[numberTraces][2] = 0;
                                traces.FracturesId[numberTraces][3] = 1;
                            }
                            else if (t[0] == c && t[1] == d){
                                traces.FracturesId[numberTraces][2] = 1;
                                traces.FracturesId[numberTraces][3] = 0;
                            }
                            else{
                                traces.FracturesId[numberTraces][2] = 1;
                                traces.FracturesId[numberTraces][3] = 1;
                            }
                            numberTraces++;
                            fractures.NumTracce[id1]++;
                            fractures.NumTracce[id2]++;

                        };
                    };
                }
            }
        }
    }
    traces.NumberTraces = numberTraces;

    // Copia gli elementi della mappa in un vettore di coppie
    vector<pair<unsigned int, double>> mapElements(traces.Lengths.begin(), traces.Lengths.end());

    // Ordina il vettore in base ai valori
    sort(mapElements.begin(), mapElements.end(), compareByValue);


    for(unsigned int i = 0; i < fractures.NumberFractures; i++)
    {
        for(auto& couple : mapElements){
            if(traces.FracturesId[couple.first][0] == i){
                fractures.tracce[i].push_back(make_tuple(couple.first, traces.FracturesId[couple.first][2], couple.second));
            }
            else if(traces.FracturesId[couple.first][1] == i){
                fractures.tracce[i].push_back(make_tuple(couple.first, traces.FracturesId[couple.first][3], couple.second));
            }
        }
        stable_sort(fractures.tracce[i].begin(), fractures.tracce[i].end(), [](const auto& a, const auto& b) {
            return get<1>(a) < get<1>(b);
        });
    }

    OutputFile(traces,fractures);

    /* for (unsigned int id = 0; id < fractures.NumberFractures; id++) {
        for (unsigned int j = 0; j < fractures.NumTracce[id]; j++)
            for (unsigned int idTrace = 0; idTrace < fractures.tracce[id][j][0]; idTrace++){
                unsigned int n = fractures.Vertices[id].cols();
                for (unsigned int column = 0; column < n; column++){
                    Vector3d firstPoint = fractures.Vertices.col(column);
                    Vector3d secondPoint = fractures.Vertices.col((column+1)%n);
                    if ((traces.Vertices[idTrace].col(0)-firstPoint).cross(secondPoint-firstPoint)==Zeros(3)){
                        unsigned int k = 1;
                        while((traces.Vertices[idTrace].col(1)-currentPoint).cross(consecutivePoint-currentPoint)==Zeros(3)){
                            currentPoint = fractures.Vertices.col(column+k);
                            consecutivePoint = fractures.Vertices.col((column+k+1)%n);
                            k++;
                        }
                    }
                }
            }
        }
    }
    */



    return 0;
}
