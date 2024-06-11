#include <iostream>
#include <GeometryLibrary.hpp>
#include <Utils.hpp>

using namespace std;
using namespace fracturesLib;

int main()
{
    Fractures fractures;
    string filepath = "./FR362_data.txt";
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


    PolygonalMesh mesh;
    Vector3d initialPoint;
    unsigned int idInitialEdge;
    bool out;

    for (unsigned int idFrac = 0; idFrac < fractures.numberOfFractures; idFrac++){
        unsigned int numberOfVertices = fractures.vertices[idFrac].cols();
        mesh.verticesCell2Ds.resize(mesh.numberCell2Ds+1);
        mesh.edgesCell2Ds.resize(mesh.numberCell2Ds+1);
        for(unsigned int v = 0; v < numberOfVertices; v++){
            mesh.coordinateCell0Ds.push_back(fractures.vertices[idFrac].col(v));
            mesh.coordinateCell1Ds.push_back({mesh.numberCell0Ds+v,(mesh.numberCell0Ds+v+1)%numberOfVertices});
            mesh.isOn1D.push_back(true);
            mesh.hasNeigh.push_back(false);
            mesh.neighCell1Ds.push_back(mesh.numberCell2Ds);
            mesh.verticesCell2Ds[mesh.numberCell2Ds].push_back(mesh.numberCell0Ds+v);
            mesh.edgesCell2Ds[mesh.numberCell2Ds].push_back(mesh.numberCell1Ds+v);
        }
        mesh.numberCell0Ds += numberOfVertices;
        mesh.numberCell1Ds += numberOfVertices;
        mesh.numberCell2Ds++;
        mesh.isOn2D.push_back(true);

        unsigned int numberOfCell2DsNow = mesh.numberCell2Ds;

        for(unsigned int k = 0; k < fractures.traces[idFrac].size(); k++){
            unsigned int idTrace = fractures.traces[idFrac][k];
            bool found = false;
            line trace;
            trace.point = traces.vertices[idTrace].col(0);
            trace.direction = traces.vertices[idTrace].col(1) - trace.point;
            for(unsigned int idCell2D = numberOfCell2DsNow-1; idCell2D < mesh.numberCell2Ds; idCell2D++){
                for(unsigned int idCell1D : mesh.edgesCell2Ds[idCell2D]){
                    line edge;
                    edge.point = mesh.coordinateCell0Ds[mesh.verticesCell1Ds[idCell1D][0]];
                    edge.direction = mesh.coordinateCell0Ds[mesh.verticesCell1Ds[idCell1D][1]] - edge.point;
                    Vector3d test = (trace.direction).cross(edge.direction);
                    if (!almostEqual(test[0],0,tol) || !almostEqual(test[1],0,tol) || !almostEqual(test[2],0,tol)){
                        VectorXd q = linesIntersection(trace,edge); // Q,t,s
                        if ((q[4]>= (0-tol)) && (q[4]<=(1+tol))){ // Q[4] è s
                            found = true;
                            out = false;
                            initialPoint = q.head(3);
                            idInitialEdge = idCell1D;
                            mesh.numberCell0Ds++;
                            mesh.coordinateCell0Ds.push_back(initialPoint);
                            break;
                        };
                    };
                }
                if (found){
                    break;
                };
            }
            // non passante?
            // if (found){
            //     mesh.isOn[initialEdge] = false;
            //     mesh.numberCell1Ds += 2;
            //     mesh.verticesCell1Ds.push_back({mesh.verticesCell1Ds[initialEdge][0],mesh.numberCell0Ds});
            //     mesh.verticesCell1Ds.push_back({mesh.numberCell0Ds,mesh.verticesCell1Ds[initialEdge][1]});
            //     mesh.isOn1D.push_back(true); mesh.isOn1D.push_back(true);
            // }
            while (!out){
                for (unsigned int neigh : mesh.neighCell1Ds[idInitialEdge]){
                    bool polygon = 0;
                    unsigned int idInitialPoint = mesh.numberCell0Ds;
                    for (unsigned int idCell1D : mesh.edgesCell2Ds[neigh]){
                        line edge;
                        edge.point = mesh.coordinateCell0Ds[mesh.verticesCell1Ds[idCell1D][0]];
                        edge.direction = mesh.coordinateCell0Ds[mesh.verticesCell1Ds[idCell1D][1]] - edge.point;
                        Vector3d test = (trace.direction).cross(edge.direction);
                        if (!almostEqual(test[0],0,tol) || !almostEqual(test[1],0,tol) || !almostEqual(test[2],0,tol)){
                            VectorXd q = linesIntersection(trace,edge);
                            double t = q[3]; // Q,t,s
                            if (almostEqual(q.head(3),initialPoint,tol)){
                                mesh.isOn[idCell1D] = false;

                                mesh.numberCell1Ds++;
                                mesh.verticesCell1Ds.push_back({mesh.verticesCell1Ds[idCell1D][0],idInitialPoint});
                                mesh.neighCell1Ds.push_back(mesh.numberCell2Ds+polygon);
                                mesh.isOn1D.push_back(true);
                                mesh.verticesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.verticesCell1Ds[idCell1D][0]);
                                mesh.verticesCell2Ds[mesh.numberCell2Ds+polygon].push_back(idInitialPoint);
                                mesh.edgesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.numberCell1Ds);

                                polygon = !polygon;

                                mesh.numberCell1Ds++;
                                mesh.verticesCell1Ds.push_back({idInitialPoint,mesh.verticesCell1Ds[idCell1D][1]});
                                mesh.neighCell1Ds.push_back(mesh.numberCell2Ds+polygon);
                                mesh.isOn1D.push_back(true);
                                mesh.verticesCell2Ds[mesh.numberCell2Ds+polygon].push_back(idInitialPoint);
                                mesh.verticesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.verticesCell1Ds[idCell1D][0]);
                                mesh.edgesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.numberCell1Ds);
                            }
                            else if (t>=(0-tol) && t<=(1+tol)){
                                double s = q[4];

                                mesh.numberCell0Ds++;
                                mesh.coordinateCell0Ds.push_back(q.head(3));

                                mesh.isOn[idCell1D] = false;

                                mesh.numberCell1Ds++;
                                mesh.verticesCell1Ds.push_back({mesh.verticesCell1Ds[idCell1D][0],mesh.numberCell0Ds});
                                mesh.neighCell1Ds.push_back(mesh.numberCell2Ds+polygon);
                                mesh.isOn1D.push_back(true);
                                mesh.verticesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.verticesCell1Ds[idCell1D][0]);
                                mesh.verticesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.numberCell0Ds);
                                mesh.edgesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.numberCell1Ds);

                                polygon = !polygon;

                                mesh.numberCell1Ds++;
                                mesh.verticesCell1Ds.push_back({mesh.numberCell0Ds,mesh.verticesCell1Ds[idCell1D][1]});
                                mesh.neighCell1Ds.push_back(mesh.numberCell2Ds+polygon);
                                mesh.isOn1D.push_back(true);
                                mesh.verticesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.numberCell0Ds);
                                mesh.verticesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.verticesCell1Ds[idCell1D][1]);
                                mesh.edgesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.numberCell1Ds);
                            }
                            else {
                                mesh.verticesCell2Ds[mesh.numberCell2Ds+polygon].push_back(mesh.verticesCell1Ds[idCell1D][0]);
                                mesh.edgesCell2Ds[mesh.numberCell2Ds+polygon].push_back(idCell1D);
                            }

                            // Ed EF?
                        }
                    }
                    mesh.isOn2D[neigh] = false;
                    mesh.isOn2D.push_back(true); mesh.isOn2D.push_back(true);
                }
                if (!(s>=(0-tol) && s<=(1+tol))){
                    out = true;
                }
            }

        }

    }
}













/*
    for (unsigned int idFrac = 0; idFrac < fractures.numberOfFractures; idFrac ++){
        SubFracture original;
        original.Vertices = fractures.vertices[idFrac];
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
    */

return 0;
}
