#include <GeometryLibrary.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <map>

using namespace std;
using namespace Eigen;

namespace fracturesLib{

bool importFractures(const string& filename, Fractures& fractures) {
    ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    string header;
    getline(file, header); // Read the header line to discard it

    string line;
    getline(file,line);
    fractures.numberOfFractures = stoi(line);
    int numberOfFractures = stoi(line);
    fractures.vertices.resize(numberOfFractures);
    char pv;
    int numberOfVertices;
    unsigned int k;

    while (numberOfFractures--) {
        getline(file, line);
        getline(file, line);
        stringstream ss(line);
        ss >> k >> pv >> numberOfVertices;

        Matrix<double, 3, Dynamic> actualVert(3, numberOfVertices);

        getline(file, line);
        string val;
        // Read Vertices
        for(int j = 0; j < 3; j++){
            for(int i = 0; i < numberOfVertices; i++){
                file >> val;
                actualVert(j, i) = stod(val);
            }
        }
        file >> pv;

        // Add fracture data to the fractures map
        fractures.vertices[k] = actualVert;
    }
    file.close();
    return true;
}

double distanceSquared(const Vector3d& a,const Vector3d& b){
    return pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2);
}

bool compareByValue(const pair<unsigned int, double> &a, const pair<unsigned int, const double> &b) {
    return a.second > b.second;
}

void outputFile(Traces& traces, Fractures& fractures)
{
    string outputFracts = "Traces.txt";
    string outputFractsTraces = "Fratture-Tracce.txt";
    ofstream ofs1(outputFracts);
    ofstream ofs2(outputFractsTraces);

    if (ofs1.fail() || ofs2.fail())
    {
        cerr << "Error in creating the output file" << endl;
        return;
    }

    ofs1 << "# Number of Traces" << endl;
    ofs1 << traces.numberOfTraces << endl;
    ofs1 << "# TraceId; FracturesId1; FracturesId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

    for(unsigned int i = 0; i < traces.numberOfTraces;i++)
    {
        ofs1 << i << ";" << traces.fracturesId[i][0] << ";" << traces.fracturesId[i][1] << ";" << traces.vertices[i](0,0) << ";" << traces.vertices[i](1,0) << ";" << traces.vertices[i](2,0) << ";" << traces.vertices[i](0,1) << ";" << traces.vertices[i](1,1) << ";" << traces.vertices[i](2,1) << endl;
    }

    for(unsigned int i = 0; i < fractures.numberOfFractures; i++)
    {
        ofs2 << "# fracturesactureId; NumTraces" << endl;
        ofs2 << i << ";" << fractures.numberOfTraces[i] << endl;
        ofs2 << "# TraceId; Tips; Length" << endl;
        for (const auto& elem : fractures.traces[i]) {
            ofs2 << get<0>(elem) << ";" << get<1>(elem) << ";" << get<2>(elem) << endl;
        }
    }


}

bool areClose(Fractures& fractures,const unsigned int& id1, const unsigned int& id2, const double& tol){
    Vector3d c1 = {0,0,0};
    const unsigned int n1 = fractures.vertices[id1].cols();
    Vector3d c2 = {0,0,0};
    const unsigned int n2 = fractures.vertices[id2].cols();

    for (unsigned int i = 0; i < 3; i++){
        for (unsigned int j = 0; j < n1; j++){
            c1[i] += fractures.vertices[id1](i,j);
        }
        c1[i] /= n1;
        for (unsigned int j = 0; j < n2; j++){
            c2[i] += fractures.vertices[id2](i,j);
        }
        c2[i] /= n2;
    }

    Vector3d vertex;
    VectorXd raysSquared1;
    raysSquared1.resize(n1);
    for (unsigned int i = 0; i < n1; i++){
        vertex = fractures.vertices[id1].col(i);
        raysSquared1[i] = distanceSquared(c1,vertex);
    }
    VectorXd raysSquared2;
    raysSquared2.resize(n2);
    for (unsigned int i = 0; i < n2; i++){
        vertex = fractures.vertices[id2].col(i);
        raysSquared2[i] = distanceSquared(c2,vertex);
    }

    double r1 = *max_element(raysSquared1.begin(), raysSquared1.end());
    double r2 = *max_element(raysSquared2.begin(), raysSquared2.end());

    return distanceSquared(c1,c2) <= (r1 + r2 + 2*sqrt(r1)*sqrt(r2) + tol);
}

Vector4d plane(const unsigned int& id, Fractures& fractures)
{
    const Vector3d v0 = fractures.vertices[id].col(0);
    const Vector3d v1 = fractures.vertices[id].col(1);
    const Vector3d v2 = fractures.vertices[id].col(2);

    Vector4d coeff = {0,0,0,0};

    coeff[0] = (v1[1]-v0[1])*(v2[2]-v0[2])-(v1[2]-v0[2])*(v2[1]-v0[1]);
    coeff[1] = -((v1[0]-v0[0])*(v2[2]-v0[2])-(v2[0]-v0[0])*(v1[2]-v0[2]));
    coeff[2] = (v1[0]-v0[0])*(v2[1]-v0[1])-(v2[0]-v0[0])*(v1[1]-v0[1]);
    coeff[3] = -(coeff[0]*v0[0]+coeff[1]*v0[1]+coeff[2]*v0[2]); //da capire i segni

    return coeff;
}

line planesIntersection(const Vector4d &coeff1, const Vector4d &coeff2, const double &tol)
{
    line output;
    output.direction = {0,0,0};
    output.point = {0,0,0};
    Vector3d v1;
    Vector3d v2;
    v1 = coeff1.head(3);
    v2 = coeff2.head(3);

    output.direction = v2.cross(v1);

    if(!almostEqual(output.direction[2],0,tol))
    {
        Matrix<double,2,2> M;
        M << v1[0], v1[1],
            v2[0],v2[1];
        Vector2d b = {-coeff1[3],-coeff2[3]};
        Vector2d P = M.colPivHouseholderQr().solve(b);
        output.point[0] = P[0];
        output.point[1] = P[1];
    }
    else
    {
        Matrix<double,2,2> M;
        M << v1[0], v1[2],
            v2[0],v2[2];
        Vector2d b = {-coeff1[3],-coeff2[3]};
        Vector2d P = M.colPivHouseholderQr().solve(b);
        output.point[0] = P[0];
        output.point[2] = P[1];
    }

    return output;
}

VectorXd linesIntersection(const line& r,const line& rj){  //primi 3 punto di inters. poi t e s
    VectorXd q;
    q.resize(5);
    q[4] = -2; //così se non c'è soluzione non salvo
    double t;
    double s;
    Matrix<double,3,2> A;
    A << r.direction[0], rj.direction[0],
        r.direction[1], rj.direction[1],
        r.direction[2], rj.direction[2];

    Vector3d b = {rj.point[0] - r.point[0],
                  rj.point[1] - r.point[1],
                  rj.point[2] - r.point[2]};

    Vector2d k = A.colPivHouseholderQr().solve(b); //non è quadrata la matrice
    t = k[0];
    s = -k[1]; //da capire perchè ci va il -
    q << rj.point[0] + rj.direction[0]*s,
        rj.point[1] + rj.direction[1]*s,
        rj.point[2] + rj.direction[2]*s,
        t,
        s;
    return q;
}

Vector4d intervalsIntersection(const Vector4d& q,const double& tol){
    double a = min(q[0],q[1]);
    double b = max(q[0],q[1]);
    double c = min(q[2],q[3]);
    double d = max(q[2],q[3]);

    // Calcola l'estremo sinistro dell'intersezione
    double sx = max(a, c);
    // Calcola l'estremo destro dell'intersezione
    double dx = min(b, d);
    //Se gli intervalli non si sovrappongono, l'intersezione sarà vuota
    if (sx > (dx-tol)) {
        sx = dx = numeric_limits<double>::quiet_NaN(); // Non un numero
    }
    double other_sx = (a < (c+tol)) ? a : c; // other_sx è pari ad a se a<c, altrimenti è pari a c
    double other_dx = (d > (b-tol)) ? d : b; // other_dx è pari a d se d>b, altrimenti è pari a b
    Vector4d output = {sx,dx,other_sx,other_dx}; // in ordine restituiamo l'intervallo di intersezione e gli altri due estremi ordinati
    return output;
}

bool almostEqual(const double& a, const double& b, const double& tol) {
    return fabs(a - b) < tol;
}

bool arePlanesParallel(const Vector4d& v1, const Vector4d& v2, double& tol) {

    Vector3d w1 = v1.head(3);
    Vector3d w2 = v2.head(3);
    Vector3d v = w1.cross(w2);

    bool par = true;

    if(abs(v[0]) > tol || abs(v[1]) > tol || abs(v[2]) > tol)
        par = false;

    return par;
}

bool compareTuple(const tuple<unsigned int, bool, double>& a, const tuple<unsigned int, bool, double>& b) {
    return get<1>(a) < get<1>(b); // Ordina in base al valore booleano all'interno delle tuple
}

}

/*
double posizionePuntoPiano(const Vector4d& coeffPiano, const Vector3d& coordPunto){
    double a = coeffPiano[0];
    double b = coeffPiano[1];
    double c = coeffPiano[2];
    double d = coeffPiano[3];
    double x = coordPunto[0];
    double y = coordPunto[1];
    double z = coordPunto[2];

    return a*x + b*y +c*z + d;
}

void splitSubFractures(SubFracture& subFract, const Fractures& fractures,const Traces& traces,PolygonalMesh& mesh,const unsigned int& idFrac,double& tol){
    if (subFract.traceId.size()==0){
        //importSubFract
        //cancella SubFract
    }
    else{
        array<SubFracture,2> sottofratt;
        unsigned int idTraccia = subFract.traceId[0];
        line currTrac;
        currTrac.point = traces.Vertices[idTraccia].col(0); //punto iniziale
        Vector3d puntoFinTraccia = traces.Vertices[idTraccia].col(1);
        //estraggo direzione traccia
        currTrac.direction = puntoFinTraccia-currTrac.point;
        bool divisione = 0;
        for(unsigned int latoFratt = 0; latoFratt < subFract.Vertices.cols(); latoFratt ++){
            Line currLato;
            currLato.point = subFract.Vertices.col(latoFratt);
            currLato.direction = subFract.Vertices.col((latoFratt+1)%subFract.Vertices.cols())-currLato.point;
            // mi assicuro che ci sia intersezione tra traccia e lato con cross
            Vector3d test = (currTrac.direction).cross(currLato.direction);
            if(!almostEqual(test[0],0,tol) || !almostEqual(test[1],0,tol) || !almostEqual(test[2],0,tol)){
                VectorXd Q = PuntiIntersRetta(currTrac, currLato); // Q,t,s
                if ((Q[4]>= (0-tol)) && (Q[4]<=(1+tol))){ // Q[4] è s
                    sottofratt[divisione].Vertices.conservativeResize(3, sottofratt[divisione].Vertices.cols() + 2);
                    sottofratt[!divisione].Vertices.conservativeResize(3, sottofratt[!divisione].Vertices.cols() + 1);
                    sottofratt[divisione].Vertices.col(sottofratt[divisione].Vertices.cols() - 2) = currLato.point;
                    sottofratt[divisione].Vertices.col(sottofratt[divisione].Vertices.cols() - 1) = currLato.point+Q[4]*currLato.direction;
                    sottofratt[divisione].VerticesId.push_back(subFract.VerticesId[latoFratt]);
                    sottofratt[divisione].VerticesId.push_back(mesh.NumberCell0Ds++);
                    mesh.CoordinateCell0Ds.push_back(currLato.point+Q[4]*currLato.direction);
                    mesh.VerticesCell1Ds[subFract.EdgesId[latoFratt]] = {numeric_limits<unsigned int>::max(),numeric_limits<unsigned int>::max()};
                    mesh.VerticesCell1Ds.push_back({subFract.VerticesId[latoFratt], mesh.NumberCell0Ds});
                    mesh.NumberCell1Ds++;
                    sottofratt[divisione].EdgesId.push_back(mesh.NumberCell1Ds);

                    divisione = !divisione;

                    sottofratt[divisione].Vertices.col(sottofratt[divisione].Vertices.cols() - 1) = currLato.point+Q[4]*currLato.direction;
                    sottofratt[divisione].VerticesId.push_back(mesh.NumberCell0Ds);
                    mesh.VerticesCell1Ds[subFract.EdgesId[latoFratt]] = {numeric_limits<unsigned int>::max(),numeric_limits<unsigned int>::max()};



                    array<unsigned int,2> temp = {subFract.VerticesId[latoFratt], mesh.NumberCell0Ds};
                    mesh.VerticesCell1Ds.push_back(temp);
                    mesh.VerticesCell1Ds.push_back({mesh.NumberCell0Ds, subFract.VerticesId[(latoFratt+1)%subFract.Vertices.cols()]});
                }
                else{
                    sottofratt[divisione].Vertices.conservativeResize(3, sottofratt[divisione].Vertices.cols() + 1);
                    sottofratt[divisione].Vertices.col(sottofratt[divisione].Vertices.cols() - 1) = currLato.point;
                }
            }
        }
        for(unsigned int m = 1; m < subFract.traceId.size(); m++){
            unsigned int idTr = subFract.traceId[m];
            Vector3d puntoInTr = traces.Vertices[idTr].col(0);
            Vector3d puntoFinTr = traces.Vertices[idTr].col(1);
            Vector4d pianoTemp;
            if(traces.FracturesId[idTr][0] == idFrac){
                pianoTemp = fractures.CoeffPiano[traces.FracturesId[idTr][1]];
            }
            else if(traces.FracturesId[idTr][1] == idFrac){
                pianoTemp = fractures.CoeffPiano[traces.FracturesId[idTr][0]];
            };

            if (posizionePuntoPiano(pianoTemp, sottofratt[0].Vertices.col(0))>0){
                if(posizionePuntoPiano(pianoTemp,puntoInTr) > 0 && posizionePuntoPiano(pianoTemp,puntoFinTr) > 0){
                    sottofratt[0].traceId.push_back(idTr);
                }
                else if(posizionePuntoPiano(pianoTemp,puntoInTr) < 0 && posizionePuntoPiano(pianoTemp,puntoFinTr) < 0){
                    sottofratt[1].traceId.push_back(idTr);
                }
            }
            else if(posizionePuntoPiano(pianoTemp, sottofratt[1].Vertices.col(1))>0){
                if(posizionePuntoPiano(pianoTemp,puntoInTr) > 0 && posizionePuntoPiano(pianoTemp,puntoFinTr) > 0){
                    sottofratt[1].traceId.push_back(idTr);
                }
                else if(posizionePuntoPiano(pianoTemp,puntoInTr) < 0 && posizionePuntoPiano(pianoTemp,puntoFinTr) < 0){
                    sottofratt[0].traceId.push_back(idTr);
                }
            }
            else{
                sottofratt[0].traceId.push_back(idTr);
                sottofratt[1].traceId.push_back(idTr);
            };
        }
        splitSubFractures(sottofratt[0],fractures,traces,mesh,idFrac,tol);
        splitSubFractures(sottofratt[1],fractures,traces,mesh,idFrac,tol);
    }
}

void importSubFract(SubFracture& subFr, PolygonalMesh& mesh){
    mesh.VerticesCell2Ds[mesh.NumberCell2Ds++] = subFr.traceId;
    mesh.EdgesCell2Ds[mesh.NumberCell2Ds] = subFr.EdgesId;
}


};
*/
