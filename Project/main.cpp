#include <iostream>
#include <GeometryLibrary.hpp>
#include <Utils.hpp>

using namespace std;
using namespace fracturesLib;

int main()
{
    Fractures fractures;
    string filepath = "./FR3_data.txt";
    if(!importFractures(filepath, fractures))
    {
        return 1;
    }

    double tol = max(pow(10,-10), numeric_limits<double>::epsilon());

    fractures.numberOfTraces.resize(fractures.numberOfFractures);
    fractures.numberOfTraces.assign(fractures.numberOfFractures,0);
    fractures.planeCoeffs.resize(fractures.numberOfFractures);
    fractures.planeCoeffs.assign(fractures.numberOfFractures, Vector4d::Zero());
    Traces traces;
    unsigned int numberOfTraces = 0;
    for (unsigned int id1 = 0; id1 < fractures.numberOfFractures; id1++) {
        for (unsigned int id2 = id1+1; id2 < fractures.numberOfFractures; id2++) {
            if (areClose(fractures,id1,id2,tol)){
                if(fractures.planeCoeffs[id1].isZero()){
                    fractures.planeCoeffs[id1] = plane(id1,fractures);
                }
                if(fractures.planeCoeffs[id2].isZero()){
                    fractures.planeCoeffs[id2] = plane(id2,fractures);
                }
                if(!arePlanesParallel(fractures.planeCoeffs[id1],fractures.planeCoeffs[id2],tol)){
                    line r = planesIntersection(fractures.planeCoeffs[id1],fractures.planeCoeffs[id2],tol);
                    line r_j;
                    Matrix<double,4,4> intersectionPoints;
                    unsigned int points = 0;
                    for (unsigned int i = 0; i < 2; i++){
                        unsigned int currentId = (i==0) ? id1:id2;
                        unsigned int n = fractures.vertices[currentId].cols();
                        for (unsigned int j=0; j < n; j++){
                            r_j.point = fractures.vertices[currentId].col(j);
                            r_j.direction = fractures.vertices[currentId].col((j+1)%n)-fractures.vertices[currentId].col(j);
                            // mi assicuro che ci sia intersezione tra r ed r_j con cross
                            Vector3d test = (r.direction).cross(r_j.direction);
                            if(!almostEqual(test[0],0,tol) || !almostEqual(test[1],0,tol) || !almostEqual(test[2],0,tol)){
                                VectorXd q = linesIntersection(r,r_j); // Q,t,s
                                if ((q[4]>= (0-tol)) && (q[4]<=(1+tol))){ // Q[4] è s
                                    intersectionPoints.col(points) = q.head(4);
                                    points++;
                                };
                            };
                        }
                    }
                    if (points == 4){
                        Vector4d t_star = intersectionPoints.row(3);
                        Vector4d t = intervalsIntersection(t_star,tol);
                        if (!isnan(t[0])){
                            array<unsigned int,4> v = {id1,id2,0,0};
                            traces.fracturesId.push_back(v);
                            Matrix<double,3,2> vertices;
                            vertices.col(0) = r.point + t[0]*r.direction;
                            vertices.col(1) = r.point + t[1]*r.direction;
                            traces.lengths[numberOfTraces] = sqrt(distanceSquared(vertices.col(0),vertices.col(1)));
                            traces.vertices.push_back(vertices);

                            // Tips
                            double a = min(t_star[0],t_star[1]);
                            double b = max(t_star[0],t_star[1]);
                            double c = min(t_star[2],t_star[3]);
                            double d = max(t_star[2],t_star[3]);
                            if (t[0] == a && t[1] == b){
                                traces.fracturesId[numberOfTraces][2] = 0;
                                traces.fracturesId[numberOfTraces][3] = 1;
                            }
                            else if (t[0] == c && t[1] == d){
                                traces.fracturesId[numberOfTraces][2] = 1;
                                traces.fracturesId[numberOfTraces][3] = 0;
                            }
                            else{
                                traces.fracturesId[numberOfTraces][2] = 1;
                                traces.fracturesId[numberOfTraces][3] = 1;
                            }
                            numberOfTraces++;
                            fractures.numberOfTraces[id1]++;
                            fractures.numberOfTraces[id2]++;

                        };
                    };
                }
            }
        }
    }
    traces.numberOfTraces = numberOfTraces;

    // Copia gli elementi della mappa in un vettore di coppie
    vector<pair<unsigned int, double>> mapElements(traces.lengths.begin(), traces.lengths.end());

    // Ordina il vettore in base ai valori
    sort(mapElements.begin(), mapElements.end(), compareByValue);

    fractures.traces.resize(fractures.numberOfFractures);
    for (unsigned int i = 0; i < fractures.numberOfFractures; i++) {
        fractures.traces[i].resize(fractures.numberOfTraces[i]);
        unsigned int index = 0;
        for (const auto& couple : mapElements) {
            if (traces.fracturesId[couple.first][0] == i) {
                fractures.traces[i][index++] = make_tuple(couple.first, traces.fracturesId[couple.first][2], couple.second);
            }
            else if (traces.fracturesId[couple.first][1] == i) {
                fractures.traces[i][index++] = make_tuple(couple.first, traces.fracturesId[couple.first][3], couple.second);
            }
        }

        // Ordinamento stabile
        stable_sort(fractures.traces[i].begin(), fractures.traces[i].end(), [](const auto& a, const auto& b) {
            return get<1>(a) < get<1>(b);
        });
    }

    traces.lengths.clear(); //per evitare il doppione di memoria

    outputFile(traces,fractures);

    // PARTE 2

    PolygonalMesh mesh;
    Vector3d coordinatesIntersectionPoint;
    unsigned int idIntersectionEdge;
    unsigned int idNewIntersectionEdge;
    unsigned int idIntersectionPoint;
    unsigned int idNewIntersectionPoint;
    unsigned int actualNeigh;
    unsigned int idSuccessivePoint;
    unsigned int idNewSuccessivePoint;
    unsigned int idPrevoiusPoint;
    unsigned int idNewPrevoiusPoint;
    unsigned int idInitialEdge0;
    unsigned int idNewInitialEdge0;
    unsigned int idInitialEdge1;
    unsigned int idNewInitialEdge1;
    unsigned int idInitialCell2D;
    bool out;
    double s;
    double t;

    for (unsigned int idFrac = 0; idFrac < fractures.numberOfFractures; idFrac++){ //ciclo su ogni frattura
        unsigned int numberOfVertices = fractures.vertices[idFrac].cols();
        mesh.verticesCell2Ds.resize(mesh.numberCell2Ds+1);
        mesh.edgesCell2Ds.resize(mesh.numberCell2Ds+1);
        for(unsigned int v = 0; v < numberOfVertices; v++){ //aggiungo la frattura su cui sono alla mesh
            mesh.coordinateCell0Ds.push_back(fractures.vertices[idFrac].col(v));
            mesh.verticesCell1Ds.push_back({mesh.numberCell0Ds+v,(mesh.numberCell0Ds+v+1)%numberOfVertices});
            mesh.isOn1D.push_back(true);
            mesh.neighCell1Ds.resize(mesh.neighCell1Ds.size()+1);
            mesh.neighCell1Ds[mesh.numberCell1Ds+v].push_back(mesh.numberCell2Ds);
            mesh.verticesCell2Ds[mesh.numberCell2Ds].push_back(mesh.numberCell0Ds+v);
            mesh.edgesCell2Ds[mesh.numberCell2Ds].push_back(mesh.numberCell1Ds+v);
        }
        mesh.numberCell0Ds += numberOfVertices;
        mesh.numberCell1Ds += numberOfVertices;
        mesh.numberCell2Ds++;
        mesh.isOn2D.push_back(true);

        unsigned int numberOfCell2DsNow = mesh.numberCell2Ds;

        for(unsigned int w = 0; w < fractures.traces[idFrac].size(); w++){ //ciclo su ogni traccia della frattura
            unsigned int idTrace = get<0>(fractures.traces[idFrac][w]); //prendo id traccia
            bool found = false;
            line trace;
            trace.point = traces.vertices[idTrace].col(0);
            trace.direction = traces.vertices[idTrace].col(1) - trace.point;
            double s1 = std::numeric_limits<double>::max();

            for (unsigned int idCell2D = numberOfCell2DsNow-1; idCell2D < mesh.numberCell2Ds; idCell2D++){ //cella 2D su cui lavoro
                if (mesh.isOn2D[idCell2D]){
                    unsigned int n = mesh.verticesCell2Ds[idCell2D].size();
                    for (unsigned int i = 0; i < n; i++){ //per ogni lato della cella
                        line edge;
                        edge.point = mesh.coordinateCell0Ds[mesh.verticesCell2Ds[idCell2D][i]];
                        edge.direction = mesh.coordinateCell0Ds[mesh.verticesCell2Ds[idCell2D][(i+1)%n]] - edge.point;
                        Vector3d test = (trace.direction).cross(edge.direction);
                        if (!almostEqual(test[0],0,tol) || !almostEqual(test[1],0,tol) || !almostEqual(test[2],0,tol)){
                            VectorXd q = linesIntersection(edge,trace); //Q, t, s
                            t = q[3];
                            s = q[4];
                            if ((t>=(0-tol)) && (t<=(1+tol))){ //se sono sul lato allora entro
                                if ((s>=(0-tol)) && (s<=(1+tol))){
                                    coordinatesIntersectionPoint = q.head(3);
                                    idIntersectionEdge = mesh.edgesCell2Ds[idCell2D][i];
                                    idPrevoiusPoint = mesh.verticesCell2Ds[idCell2D][(i-1)%n];
                                    idSuccessivePoint = mesh.verticesCell2Ds[idCell2D][(i+1)%n];
                                    found = true;
                                    break;
                                }
                                else if (s < (s1+tol) && s > (1-tol)){
                                    s1 = s;
                                    coordinatesIntersectionPoint = q.head(3);
                                    idIntersectionEdge = mesh.edgesCell2Ds[idCell2D][i];
                                    idPrevoiusPoint = mesh.verticesCell2Ds[idCell2D][(i-1)%n];
                                    idSuccessivePoint = mesh.verticesCell2Ds[idCell2D][(i+1)%n];
                                    idInitialCell2D = idCell2D;
                                }
                            }
                        }
                        if (found){
                            break;
                        }
                    }
                }
            }
            mesh.coordinateCell0Ds.push_back(coordinatesIntersectionPoint);
            idIntersectionPoint = mesh.numberCell0Ds++;
            mesh.isOn1D[idIntersectionEdge] = false; //spengo il lato da cui parto

            idInitialEdge0 = mesh.numberCell1Ds++;
            mesh.verticesCell1Ds.push_back({idIntersectionPoint,idSuccessivePoint});
            mesh.neighCell1Ds.resize(mesh.numberCell1Ds); mesh.neighCell1Ds[idInitialEdge0].push_back(mesh.numberCell2Ds);
            mesh.isOn1D.push_back(true);

            idInitialEdge1 = mesh.numberCell1Ds++;
            mesh.verticesCell1Ds.push_back({idPrevoiusPoint,idIntersectionPoint});
            mesh.neighCell1Ds.resize(mesh.numberCell1Ds); mesh.neighCell1Ds[idInitialEdge1].push_back(mesh.numberCell2Ds+1);
            mesh.isOn1D.push_back(true);

            bool pass = true;
            for (unsigned int neigh : mesh.neighCell1Ds[idIntersectionEdge]){ //per ogni vicino del lato di intersezione
                if (!found && neigh != idInitialCell2D){
                    pass = false;
                }
                if (mesh.isOn2D[neigh] && pass){
                    actualNeigh = neigh;
                    out = false;
                    while (!out){ //finchè non finisco la traccia
                        unsigned int polygon = 0;
                        //trovo la posizione di idSuccessivePoint in mesh.verticesCell2Ds[actualNeigh]
                        unsigned int indexOfInitialPoint = 0;
                        unsigned int numberOfVertices = mesh.verticesCell2Ds[actualNeigh].size();
                        for (unsigned int j = 0; j < numberOfVertices; j++){
                            if (mesh.verticesCell2Ds[actualNeigh][j] == idSuccessivePoint){
                                indexOfInitialPoint = j;//num vertici poligono precedente
                                break;
                            }
                        }

                        mesh.verticesCell2Ds.resize(mesh.verticesCell2Ds.size()+2);
                        mesh.edgesCell2Ds.resize(mesh.edgesCell2Ds.size()+2);

                        mesh.verticesCell2Ds[mesh.numberCell2Ds].push_back(idIntersectionPoint);
                        mesh.edgesCell2Ds[mesh.numberCell2Ds].push_back(idInitialEdge0);

                        mesh.neighCell1Ds.resize(mesh.neighCell1Ds.size()+2); // troveremo due nuovi lati
                        for (unsigned int v =  1; v < numberOfVertices; v++){
                            unsigned int idCell0D = mesh.verticesCell2Ds[actualNeigh][(indexOfInitialPoint+v)%numberOfVertices];
                            unsigned int idCell1D = mesh.edgesCell2Ds[actualNeigh][(indexOfInitialPoint+v)%numberOfVertices];
                            line edge;
                            edge.point = mesh.coordinateCell0Ds[idCell0D];
                            edge.direction = mesh.coordinateCell0Ds[mesh.verticesCell2Ds[actualNeigh][(indexOfInitialPoint+v+1)%numberOfVertices]] - edge.point;
                            Vector3d test = (trace.direction).cross(edge.direction);
                            VectorXd q = linesIntersection(edge,trace);
                            double t = q[3];
                            if (t>=(0-tol) && t<=(1+tol) && (!almostEqual(test[0],0,tol) || !almostEqual(test[1],0,tol) || !almostEqual(test[2],0,tol))){
                                s = q[4];

                                idNewIntersectionPoint = mesh.numberCell0Ds++;
                                mesh.coordinateCell0Ds.push_back(q.head(3));

                                idNewSuccessivePoint = idCell0D;
                                idNewPrevoiusPoint = mesh.verticesCell2Ds[actualNeigh][(indexOfInitialPoint+v+1)%numberOfVertices];

                                idNewIntersectionEdge = idCell1D;
                                mesh.isOn1D[idCell1D] = false;

                                idNewInitialEdge0 = mesh.numberCell1Ds++;
                                mesh.verticesCell1Ds.push_back({idCell0D,idNewIntersectionPoint});
                                mesh.neighCell1Ds[idNewInitialEdge0].push_back(mesh.numberCell2Ds);
                                mesh.isOn1D.push_back(true);

                                mesh.verticesCell2Ds[mesh.numberCell2Ds].push_back(idCell0D);
                                mesh.edgesCell2Ds[mesh.numberCell2Ds].push_back(idNewInitialEdge0);
                                mesh.verticesCell2Ds[mesh.numberCell2Ds].push_back(idIntersectionPoint);

                                polygon = 1;

                                idNewInitialEdge1 = mesh.numberCell1Ds++;
                                mesh.verticesCell1Ds.push_back({idIntersectionPoint,idNewPrevoiusPoint});
                                mesh.neighCell1Ds[idNewInitialEdge1].push_back(mesh.numberCell2Ds+1);
                                mesh.isOn1D.push_back(true);

                                mesh.verticesCell2Ds[mesh.numberCell2Ds+1].push_back(idNewIntersectionPoint);
                                mesh.edgesCell2Ds[mesh.numberCell2Ds+1].push_back(idNewInitialEdge1);
                            }
                            else {
                                mesh.verticesCell2Ds[mesh.numberCell2Ds+polygon].push_back(idCell0D);
                                mesh.edgesCell2Ds[mesh.numberCell2Ds+polygon].push_back(idCell1D);
                                mesh.neighCell1Ds[idCell1D].push_back(mesh.numberCell2Ds+polygon);
                            }
                        }
                        mesh.verticesCell2Ds[mesh.numberCell2Ds+1].push_back(idPrevoiusPoint);
                        mesh.edgesCell2Ds[mesh.numberCell2Ds].push_back(idInitialEdge1);

                        mesh.neighCell1Ds.resize(mesh.neighCell1Ds.size()+1); mesh.neighCell1Ds[mesh.numberCell1Ds].push_back(mesh.numberCell2Ds); mesh.neighCell1Ds[mesh.numberCell1Ds].push_back(mesh.numberCell2Ds+1);
                        mesh.verticesCell1Ds.push_back({idIntersectionPoint,idNewIntersectionPoint});
                        mesh.isOn1D.push_back(true);
                        mesh.edgesCell2Ds[mesh.numberCell2Ds].push_back(mesh.numberCell1Ds);
                        mesh.edgesCell2Ds[mesh.numberCell2Ds+1].push_back(mesh.numberCell1Ds);
                        mesh.numberCell1Ds++;

                        mesh.numberCell2Ds += 2;
                        mesh.isOn2D[actualNeigh] = false;
                        mesh.isOn2D.push_back(true); mesh.isOn2D.push_back(true);

                        idIntersectionPoint = idNewIntersectionPoint;
                        idPrevoiusPoint = idNewPrevoiusPoint;
                        idSuccessivePoint = idNewSuccessivePoint;
                        idInitialEdge0 = idNewInitialEdge0;
                        idInitialEdge1 = idNewInitialEdge1;
                        idIntersectionEdge = idNewIntersectionEdge;

                        if(almostEqual(s,0,tol) || almostEqual(s,1,tol) || (!(s>=(0-tol) && s<=(1+tol)))){
                            out = true;
                        }
                        for (unsigned int x : mesh.neighCell1Ds[idIntersectionEdge]){
                            if (actualNeigh != x && mesh.isOn2D[x]){
                                actualNeigh = x;
                            }
                        }
                    }
                    for (unsigned int i = 0; i < mesh.verticesCell2Ds[actualNeigh].size(); i++){
                        if (!mesh.isOn1D[mesh.edgesCell2Ds[actualNeigh][i]]){
                            mesh.verticesCell2Ds[actualNeigh].insert(mesh.verticesCell2Ds[actualNeigh].begin() + i+1, idIntersectionPoint);
                            mesh.edgesCell2Ds[actualNeigh].erase(mesh.edgesCell2Ds[actualNeigh].begin() + i);
                            mesh.edgesCell2Ds[actualNeigh].insert(mesh.edgesCell2Ds[actualNeigh].begin() + i, idInitialEdge1);
                            mesh.edgesCell2Ds[actualNeigh].insert(mesh.edgesCell2Ds[actualNeigh].begin() + i+1, idInitialEdge0);
                            break;
                        }
                    }
                    mesh.neighCell1Ds[idInitialEdge0].push_back(actualNeigh);
                    mesh.neighCell1Ds[idInitialEdge1].push_back(actualNeigh);
                }
            }
        }
    }
    return 0;
}
