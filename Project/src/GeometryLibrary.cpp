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

namespace FracturesLib{

bool importFracture(const string& filename, Fractures& fracture) {
    ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    string header;
    getline(file, header); // Read the header line to discard it

    string line;
    getline(file,line);
    fracture.NumberFractures = stoi(line);
    int numFrac = stoi(line);
    fracture.Vertices.resize(numFrac);
    char pv;
    int numVertices;
    unsigned int k;

    while (numFrac--) {
        getline(file, line);
        getline(file, line);
        stringstream ss(line);
        ss >> k >> pv >> numVertices;

        Matrix<double, 3, Dynamic> actualVert(3, numVertices);

        getline(file, line);
        string val;
        // Read Vertices
        for(int j = 0; j < 3; j++){
            for(int i = 0; i < numVertices; i++){
                file >> val;
                actualVert(j, i) = stod(val);
            }
        }
        file >> pv;

        // Add fracture data to the fractures map
        fracture.Vertices[k] = actualVert;
    }
    file.close();
    return true;
}

double distanceSquared(const Vector3d& A,const Vector3d& B){
    return pow(A[0]-B[0],2) + pow(A[1]-B[1],2) + pow(A[2]-B[2],2);
}

bool compareByValue(const pair<unsigned int, double> &a, const pair<unsigned int, const double> &b) {
    return a.second > b.second;
}

void OutputFile(Traces& TR, Fractures& FR)
{
    string Tracce = "Traces.txt";
    string FrattTracc = "Fratture-Tracce.txt";
    ofstream ofs(Tracce);
    ofstream ofs2(FrattTracc);

    if (ofs.fail() || ofs2.fail())
    {
        cout << "Impossibile creare il file di output" << endl;
        return;
    }

    ofs << "# Number of Traces" << endl;
    ofs << TR.NumberTraces << endl;
    ofs << "# TraceId; FracturesId1; FracturesId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

    for(unsigned int i = 0; i < TR.NumberTraces;i++)
    {
        ofs << i << ";" << TR.FracturesId[i][0] << ";" << TR.FracturesId[i][1] << ";" << TR.Vertices[i](0,0) << ";" << TR.Vertices[i](1,0) << ";" << TR.Vertices[i](2,0) << ";" << TR.Vertices[i](0,1) << ";" << TR.Vertices[i](1,1) << ";" << TR.Vertices[i](2,1) << endl;
    }

    for(unsigned int i = 0; i < FR.NumberFractures; i++)
    {
        ofs2 << "# FractureId; NumTraces" << endl;
        ofs2 << i << ";" << FR.NumTracce[i] << endl;
        ofs2 << "# TraceId; Tips; Length" << endl;
        for (const auto& elem : FR.tracce[i]) {
            ofs2 << get<0>(elem) << ";" << get<1>(elem) << ";" << get<2>(elem) << endl;
        }
    }


}

bool areClose(Fractures& fracture,const unsigned int& Id1,const unsigned int& Id2, const double& tol){
    Vector3d C1 = {0,0,0};
    const unsigned int n1 = fracture.Vertices[Id1].cols();
    Vector3d C2 = {0,0,0};
    const unsigned int n2 = fracture.Vertices[Id2].cols();

    for(unsigned int i=0; i<3; i++){
        for (unsigned int j=0; j<n1; j++){
            C1[i] += fracture.Vertices[Id1](i,j);
        }
        C1[i] /= n1;
        for (unsigned int j=0; j<n2; j++){
            C2[i] += fracture.Vertices[Id2](i,j);
        }
        C2[i] /= n2;
    }

    VectorXd rays1;
    for(unsigned int i=0; i<n1; i++){
        rays1.resize(rays1.size() + 1);
        Vector3d point = fracture.Vertices[Id1].col(i);
        rays1(rays1.size() - 1) = distanceSquared(C1,point);
    }
    VectorXd rays2;
    for(unsigned int i=0; i<n2; i++){
        rays2.resize(rays2.size() + 1);
        Vector3d point = fracture.Vertices[Id2].col(i);
        rays2(rays2.size() - 1) = distanceSquared(C2,point);
    }

    double R1 = *max_element(rays1.begin(), rays1.end());
    double R2 = *max_element(rays2.begin(), rays2.end());

    return distanceSquared(C1,C2) <= ((R1+R2+(2*sqrt(R1)*sqrt(R2)))+tol); //capire bene
}

Vector4d Piano(const unsigned int& id, Fractures& FR)
{

    Vector3d v0 = FR.Vertices[id].col(0);
    Vector3d v1 = FR.Vertices[id].col(1);
    Vector3d v2 = FR.Vertices[id].col(2);

    Vector4d coeff = {0,0,0,0};

    coeff[0] = (v1[1]-v0[1])*(v2[2]-v0[2]) - (v1[2]-v0[2])*(v2[1]-v0[1]);
    coeff[1] = -((v1[0]-v0[0])*(v2[2]-v0[2])-(v2[0]-v0[0])*(v1[2]-v0[2]));
    coeff[2] = (v1[0]-v0[0])*(v2[1]-v0[1])-(v2[0]-v0[0])*(v1[1]-v0[1]);
    coeff[3] = -(coeff[0]*v0[0]+coeff[1]*v0[1]+coeff[2]*v0[2]); //da capire i segni

    return coeff;
}

Line Inter(const Vector4d& coeff1, const Vector4d& coeff2,const double& tol)
{
    Line line;
    line.direction = {0,0,0};
    line.point = {0,0,0};
    Vector3d v1;
    Vector3d v2;
    v1 = coeff1.head(3);
    v2 = coeff2.head(3);

    line.direction = v2.cross(v1);

    if(!almostEqual(line.direction[2],0,tol))
    {
        Matrix<double,2,2> M;
        M << v1[0], v1[1],
            v2[0],v2[1];
        Vector2d b = {-coeff1[3],-coeff2[3]};
        Vector2d P = M.colPivHouseholderQr().solve(b);
        line.point[0] = P[0];
        line.point[1] = P[1];
    }
    else
    {
        Matrix<double,2,2> M;
        M << v1[0], v1[2],
            v2[0],v2[2];
        Vector2d b = {-coeff1[3],-coeff2[3]};
        Vector2d P = M.colPivHouseholderQr().solve(b);
        line.point[0] = P[0];
        line.point[2] = P[1];
    }

    /*else{
     * ...inserire altri casi...
        } */

    return line;
}

VectorXd PuntiIntersRetta(const Line& r,const Line& rj){  //primi 3 punto di inters. poi t e s
    VectorXd Q;
    Q.resize(5);
    Q[4] = -2; //così se non c'è soluzione non salvo
    double t;
    double s;
    Matrix<double,3,2> A;
    A << r.direction[0], rj.direction[0],
        r.direction[1], rj.direction[1],
        r.direction[2], rj.direction[2];

    Vector3d b = {rj.point[0] - r.point[0],
                  rj.point[1] - r.point[1],
                  rj.point[2] - r.point[2]};

    // FullPivLU<MatrixXd> lu_decomp(A);
    // int rank_A = lu_decomp.rank();

    // // Creare la matrice aumentata [A | B]
    // MatrixXd augmented(A.rows(), A.cols() + 1);
    // augmented << A, b;

    // // Calcolare il rango della matrice aumentata [A | B] usando la decomposizione LU con pivotaggio completo
    // FullPivLU<MatrixXd> lu_decomp_aug(augmented);
    // int rank_augmented = lu_decomp_aug.rank();

    //if(rank_A == rank_augmented){
    Vector2d k = A.colPivHouseholderQr().solve(b); //non è quadrata la matrice
    t = k[0];
    s = -k[1]; //da capire perchè ci va il -
    Q << rj.point[0] + rj.direction[0]*s,
        rj.point[1] + rj.direction[1]*s,
        rj.point[2] + rj.direction[2]*s,
        t,
        s;
    //};
    return Q;
}

Vector4d intersection(const Vector4d& Q,const double& tol){
    double a = min(Q[0],Q[1]);
    double b = max(Q[0],Q[1]);
    double c = min(Q[2],Q[3]);
    double d = max(Q[2],Q[3]);

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

bool almostEqual(const double a, const double b, const double &tol) {
    return fabs(a - b) < tol;
}

bool arePlanesParallel(const Vector4d v1, const Vector4d v2, double &tol) {

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
    Line currTrac;
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
                divisione = !divisione;
                sottofratt[divisione].Vertices.col(sottofratt[divisione].Vertices.cols() - 1) = currLato.point+Q[4]*currLato.direction;
                mesh.CoordinateCell0Ds.push_back(currLato.point+Q[4]*currLato.direction);
                mesh.NumberCell0Ds++;
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

// void importSubFract(SubFracture& subFr, PolygonalMesh& mesh){
//     for(unsigned int latoFrat = 0; latoFrat < subFr.Vertices.cols(); latoFrat ++){


//         if(//non esiste già){
//         mesh.NumberCell0Ds++;
//         //}
//     }
// }


};
