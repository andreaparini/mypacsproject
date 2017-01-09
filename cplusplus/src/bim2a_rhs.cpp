/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "../include/bim3a_rhs.hpp"
#include "quadmesh.hpp"
#include <vector>

using namespace std;


vector<double> bim2a_rhs
(Mesh mesh, vector<double> f, std::vector<double> g){
    
    vector<double> rhs(mesh.nelem,0);
    vector< vector<double> > x(4);
    vector< vector<double> > y(4);
    vector< vector<double> > F(4);
    for(int m = 0; m < 4; m++){
        x[m].resize(mesh.nelem,0);
        y[m].resize(mesh.nelem,0);
        F[m].resize(mesh.nelem,0);
        for (int k = 0; k < mesh.nelem; k++){
            x[m][k] = mesh.p[1][ mesh.t[m][k] ];
            y[m][k] = mesh.p[2][ mesh.t[m][k] ];
            F[m][k] = f[k]*g[mesh.t[m][k]];
        }
    }
    vector<vector<double>> rhs_loc(4);
    for(int i = 0; i < 4; i++){
        rhs_loc[i].resize(mesh.nelem);
    }
    
    for(int k = 0; k < mesh.nelem; k++)
    
}