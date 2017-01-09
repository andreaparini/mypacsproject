/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#include "../include/bim3a_advection_diffusion.hpp"

using namespace std;

void bim3a_advection_diffusion
(Mesh& mesh, vector<double> alpha, vector<double> beta, sparse_matrix &A){
    
    typedef vector<vector<vector<double>>> mat3d;
    mat3d Aloc(8);
    for(int i = 0; i < 8; i++){
        Aloc[i].resize(8);
        for(int j = 0; j < 8; j++){
            Aloc[i][j].resize(mesh.nelem);
        }
    }
    
    vector< vector<double> > x(mesh.nelem);
    vector< vector<double> > y(mesh.nelem);
    vector< vector<double> > z(mesh.nelem);
    
    for(int m = 0; m < mesh.nelem; m++){
        x[m].resize(8,0);
        y[m].resize(8,0);
        z[m].resize(8,0);
    }
    for(int k = 0; k < mesh.nelem; k++){
        for(int m = 0; m < 8; m++){
            x[m][k] = mesh.p(0,mesh.t(m,k));
            y[m][k] = mesh.p(1,mesh.t(m,k));
            z[m][k] = mesh.p(2,mesh.t(m,k));            
        }
    }
    
    vector <double> psi12(mesh.nelem);
    vector <double> psi23(mesh.nelem);
    vector <double> psi43(mesh.nelem);
    vector <double> psi14(mesh.nelem);
    vector <double> psi15(mesh.nelem);
    vector <double> psi26(mesh.nelem);
    vector <double> psi37(mesh.nelem);
    vector <double> psi48(mesh.nelem);
    vector <double> psi56(mesh.nelem);
    vector <double> psi67(mesh.nelem);
    vector <double> psi87(mesh.nelem);
    vector <double> psi58(mesh.nelem);
    
    for (int k = 0; k < mesh.nelem; k++){
        psi12[k] = beta[0]* (x[k][1] - x[k][0]) +
                   beta[1]* (y[k][1] - y[k][0]) +
                   beta[2]* (z[k][1] - z[k][0]);
        
        psi23[k] = beta[0]* (x[k][2] - x[k][1]) +
                   beta[1]* (y[k][2] - y[k][1]) +
                   beta[2]* (z[k][2] - z[k][1]);
        
        psi43[k] = beta[0]* (x[k][2] - x[k][3]) +
                   beta[1]* (y[k][2] - y[k][3]) +
                   beta[2]* (z[k][2] - z[k][3]);
        
        psi14[k] = beta[0]* (x[k][3] - x[k][0]) +
                   beta[1]* (y[k][3] - y[k][0]) +
                   beta[2]* (z[k][3] - z[k][0]);
        
        psi15[k] = beta[0]* (x[k][4] - x[k][0]) +
                   beta[1]* (y[k][4] - y[k][0]) +
                   beta[2]* (z[k][4] - z[k][0]);
        
        psi26[k] = beta[0]* (x[k][5] - x[k][1]) +
                   beta[1]* (y[k][5] - y[k][1]) +
                   beta[2]* (z[k][5] - z[k][1]);
        
        psi37[k] = beta[0]* (x[k][6] - x[k][2]) +
                   beta[1]* (y[k][6] - y[k][2]) +
                   beta[2]* (z[k][6] - z[k][2]);
        
        psi48[k] = beta[0]* (x[k][7] - x[k][3]) +
                   beta[1]* (y[k][7] - y[k][3]) +
                   beta[2]* (z[k][7] - z[k][3]);
        
        psi56[k] = beta[0]* (x[k][5] - x[k][4]) +
                   beta[1]* (y[k][5] - y[k][4]) +
                   beta[2]* (z[k][5] - z[k][4]);
        
        psi67[k] = beta[0]* (x[k][6] - x[k][5]) +
                   beta[1]* (y[k][6] - y[k][5]) +
                   beta[2]* (z[k][6] - z[k][5]);
        
        psi87[k] = beta[0]* (x[k][1] - x[k][0]) +
                   beta[1]* (y[k][1] - y[k][0]) +
                   beta[2]* (z[k][1] - z[k][0]);
        
        psi58[k] = beta[0]* (x[k][6] - x[k][7]) +
                   beta[1]* (y[k][6] - y[k][7]) +
                   beta[2]* (z[k][6] - z[k][7]);
    }
    
    vector <double> bp12;
    vector <double> bp23;
    vector <double> bp43;
    vector <double> bp14;
    vector <double> bp15;
    vector <double> bp26;
    vector <double> bp37;
    vector <double> bp48;
    vector <double> bp56;
    vector <double> bp67;
    vector <double> bp87;
    vector <double> bp58;
    
    vector <double> bm12;
    vector <double> bm23;
    vector <double> bm43;
    vector <double> bm14;
    vector <double> bm15;
    vector <double> bm26;
    vector <double> bm37;
    vector <double> bm48;
    vector <double> bm56;
    vector <double> bm67;
    vector <double> bm87;
    vector <double> bm58;
    
    bernoulli(bp12,bm12,psi12);
    bernoulli(bp23,bm23,psi23);
    bernoulli(bp43,bm43,psi43);
    bernoulli(bp14,bm14,psi14);
    bernoulli(bp15,bm15,psi15);
    bernoulli(bp26,bm26,psi26);
    bernoulli(bp37,bm37,psi37);
    bernoulli(bp48,bm48,psi48);
    bernoulli(bp56,bm56,psi56);
    bernoulli(bp67,bm67,psi67);
    bernoulli(bp87,bm87,psi87);
    bernoulli(bp58,bm58,psi58);
    
    for (int k = 0; k < mesh.nelem; k++){
        bp12[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
        bm12[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
        bp23[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
        bm23[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
        bp43[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
        bm43[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
        bp14[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
        bm14[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
        
        bp15[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
        bm15[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
        bp26[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
        bm26[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
        bp37[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
        bm37[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
        bp48[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
        bm48[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
        
        bp56[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
        bm56[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
        bp67[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
        bm67[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
        bp87[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
        bm87[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
        bp58[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
        bm58[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
        
        
        
        Aloc[0][0][k] = bm12[k] + bm14[k] + bm15[k];
        Aloc[0][1][k] = -bp12[k];
        Aloc[0][3][k] = -bp14[k];
        Aloc[0][4][k] = -bp15[k];
        
        Aloc[1][0][k] = -bm12[k];
        Aloc[1][1][k] = bp12[k] + bm23[k] + bm26[k];
        Aloc[1][2][k] = -bp23[k];
        Aloc[1][5][k] = -bp26[k];
        
        Aloc[2][1][k] = -bm23[k];
        Aloc[2][2][k] = bp23[k] + bp43[k] + bm37[k];
        Aloc[2][3][k] = -bm43[k];
        Aloc[2][6][k] = -bp37[k];
        
        Aloc[3][0][k] = -bm14[k];
        Aloc[3][2][k] = -bp43[k];
        Aloc[3][3][k] = bm43[k] + bp14[k] + bm48[k];
        Aloc[3][7][k] = -bp48[k];
        
        Aloc[4][0][k] = -bm15[k];
        Aloc[4][4][k] = bp15[k] + bm56[k] + bm58[k];
        Aloc[4][5][k] = -bp56[k];
        Aloc[4][7][k] = -bp58[k];
        
        Aloc[5][1][k] = -bm26[k];
        Aloc[5][4][k] = -bm56[k];
        Aloc[5][5][k] = bp56[k] + bp26[k] + bm67[k];
        Aloc[5][6][k] = -bp67[k];
        
        Aloc[6][2][k] = -bm37[k];
        Aloc[6][5][k] = -bm67[k];
        Aloc[6][6][k] = bp37[k] + bp67[k] + bp87[k];
        Aloc[6][7][k] = -bm87[k];
        
        Aloc[7][3][k] = -bm48[k];
        Aloc[7][4][k] = -bm58[k];
        Aloc[7][6][k] = -bp87[k];
        Aloc[7][7][k] = bp48[k] + bp58[k] + bm87[k];
        
        for(int ii = 0; ii < 8; ii++){
            for(int jj = 0; jj < 8; jj++){
                A[ mesh.t(ii,k) ][ mesh.t(jj,k) ] += Aloc[ii][jj][k];
            }
        }
        
    }
    
}