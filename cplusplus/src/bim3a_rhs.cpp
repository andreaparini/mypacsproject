/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "../include/bim3a_rhs.hpp"

using namespace std;


void bim3a_rhs
(Mesh& mesh, const vector<double> f, const vector<double> g, vector<double>& rhs){
    
    if(int(rhs.size()) < mesh.nnodes){
        rhs.resize(mesh.nnodes,0);
    }
    int i;
    
    for(int k = 0; k < mesh.nelem; k++){
        for(int m = 0; m < 8; m++){
            
            i = mesh.t(m,k);
            rhs[i] += mesh.hx * mesh.hy * mesh.hz / 8 * f[k] * g[i];
            
        }
    }
}