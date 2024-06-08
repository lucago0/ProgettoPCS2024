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

    double tol = max(pow(10,-10), numeric_limits<double>::epsilon());

    fractures.NumTracce.resize(fractures.NumberFractures);
    fractures.NumTracce.assign(fractures.NumberFractures,0);
    fractures.CoeffPiano.resize(fractures.NumberFractures);
    fractures.CoeffPiano.assign(fractures.NumberFractures, Vector4d::Zero());
    Traces traces;
    unsigned int numberTraces = 0;
    for (unsigned int id1 = 0; id1<fractures.NumberFractures; id1++) {
        for (unsigned int id2 = id1+1; id2<fractures.NumberFractures; id2++) {
            if (areClose(fractures,id1,id2,tol)){
                if(fractures.CoeffPiano[id1].isZero()){
                    fractures.CoeffPiano[id1] = Piano(id1,fractures);
                }
                if(fractures.CoeffPiano[id2].isZero()){
                    fractures.CoeffPiano[id2] = Piano(id2,fractures);
                }
                if(!arePlanesParallel(fractures.CoeffPiano[id1],fractures.CoeffPiano[id2],tol)){
                    Line r = Inter(fractures.CoeffPiano[id1],fractures.CoeffPiano[id2],tol);
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
                            if(!almostEqual(test[0],0,tol) || !almostEqual(test[1],0,tol) || !almostEqual(test[2],0,tol)){
                                VectorXd Q = PuntiIntersRetta(r,r_j); // Q,t,s
                                if ((Q[4]>= (0-tol)) && (Q[4]<=(1+tol))){ // Q[4] è s
                                    intersectionPoints.col(points) = Q.head(4);
                                    points++;
                                };
                            };
                        }
                    }
                    if(points == 4){
                        Vector4d t_star = intersectionPoints.row(3);
                        Vector4d t = intersection(t_star,tol);
                        if(!isnan(t[0])){
                            array<unsigned int,4> v = {id1,id2,0,0};
                            traces.FracturesId.push_back(v);
                            Matrix<double,3,2> vertices;
                            vertices.col(0) = r.point + t[0]*r.direction;
                            vertices.col(1) = r.point + t[1]*r.direction;
                            traces.Lengths[numberTraces] = sqrt(distanceSquared(vertices.col(0),vertices.col(1)));
                            traces.Vertices.push_back(vertices);

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

    fractures.tracce.resize(fractures.NumberFractures);
    for (unsigned int i = 0; i < fractures.NumberFractures; i++) {
        fractures.tracce[i].resize(fractures.NumTracce[i]);
        unsigned int index = 0;
        for (const auto& couple : mapElements) {
            if (traces.FracturesId[couple.first][0] == i) {
                fractures.tracce[i][index++] = make_tuple(couple.first, traces.FracturesId[couple.first][2], couple.second);
            }
            else if (traces.FracturesId[couple.first][1] == i) {
                fractures.tracce[i][index++] = make_tuple(couple.first, traces.FracturesId[couple.first][3], couple.second);
            }
        }

        // Ordinamento stabile
        stable_sort(fractures.tracce[i].begin(), fractures.tracce[i].end(), [](const auto& a, const auto& b) {
            return get<1>(a) < get<1>(b);
        });
    }

    traces.Lengths.clear(); //per evitare il doppione di memoria

    OutputFile(traces,fractures);

    fractures.SubFracVert.resize(fractures.NumberFractures);

    PolygonalMesh mesh;
    for (unsigned int idFrac = 0; idFrac < fractures.NumberFractures; idFrac ++){
        SubFracture original;
        original.Vertices = fractures.Vertices[idFrac];
        for(unsigned int z = 0; z < original.Vertices.cols(); z++){
            mesh.CoordinateCell0Ds.push_back(original.Vertices.col(z));
            array<unsigned int, 2> idVertici;
            unsigned int a = mesh.NumberCell0Ds+z;
            unsigned int b = mesh.NumberCell0Ds+((z+1)%original.Vertices.cols());
            idVertici = {a,b};
            mesh.VerticesCell1Ds.push_back(idVertici);
            original.VerticesId.push_back(mesh.NumberCell0Ds+z);
            original.EdgesId.push_back(mesh.NumberCell0Ds+z);
        }
        mesh.NumberCell0Ds = mesh.NumberCell0Ds + original.Vertices.cols(); //inserisco i punti esterni originali
        mesh.NumberCell1Ds = mesh.NumberCell1Ds + original.Vertices.cols();

        for(unsigned int k = 0; k < fractures.tracce[idFrac].size(); k++){
            original.traceId.push_back(get<0>(fractures.tracce[idFrac][k]));
        }
        splitSubFractures(original,fractures,traces,mesh,idFrac,tol);
    }

    return 0;
}
