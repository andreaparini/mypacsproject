/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "../include/quadmesh.hpp"

using namespace std;

Mesh::Mesh(const double L1, const double L2, const double H1, const double H2,
           const double W1, const double W2, const unsigned int Nx_,
           const unsigned int Ny_, const unsigned int Nz_, const unsigned int region_,
           vector<unsigned int> sides_):
        L (L2 - L1),
        H (H2 - H1),
        W (W2 - W1),
        Nx (Nx_),
        Ny (Ny_),
        Nz (Nz_),
        hx ( L / double(Nx) ),
        hy ( H / double(Ny) ),
        hz ( W / double(Nz) ),
        nnodes ( (Nx+1)*(Ny+1)*(Nz+1) ),
        nelem ( Nx * Ny *Nz ),
        region (region_),
        sides(sides_)
{
    p_data = new double [3*nnodes];
    e_data = new unsigned int [11*2*(Nx*Ny + Ny*Nz + Nx*Nz)];
    t_data = new unsigned int [9*nelem];
    
    // generate points
    for (size_t k = 0; k < Nz + 1; k++){
        for (size_t i = 0; i < Nx + 1; i++){
            for (size_t j = 0; j < Ny + 1; j++){
                p(0, k*(Nx+1)*(Ny+1) + i*(Ny+1) + j) = i * hx;
                p(1, k*(Nx+1)*(Ny+1) + i*(Ny+1) + j) = j * hy;
                p(2, k*(Nx+1)*(Ny+1) + i*(Ny+1) + j) = k * hz;
            }
        }
    }
    
    // generate cube vertexes 
    unsigned int *iiv = new unsigned int[Nx*Ny*Nz];
    unsigned int count = 0;
    for(size_t k = 0; k < Nz + 1; k++){
        for(size_t i = 0; i < Nx + 1; i++ ){
            for(size_t j = 0; j< Ny + 1; j++ ){
                if(i != Nx+1 && j != Ny+1 && k != Nz+1){
                    iiv[k*Nx*Ny + i*Ny + j] = count;
                }
                count++;
            }
        }
    }
    for (size_t k = 0; k < nelem; k++){
        t(0,k) = iiv[k];
        t(1,k) = iiv[k] + Ny + 1;
        t(2,k) = iiv[k] + Ny + 2;
        t(3,k) = iiv[k] + 1;
        t(4,k) = iiv[k] + (Nx+1)*(Ny+1);
        t(5,k) = iiv[k] + (Nx+1)*(Ny+1) + Ny + 1;
        t(6,k) = iiv[k] + (Nx+1)*(Ny+1) + Ny + 2;
        t(7,k) = iiv[k] + (Nx+1)*(Ny+1) + 1;
        
    }
    delete [] iiv;
    
    //generate boundary face list
    
    //front
    
    vector< vector< unsigned int > > T(nelem);
    vector< vector< unsigned int > > order(nelem);
    
    for (size_t k = 0; k < nelem; k++){
        order[k].resize(8,0);
        T[k].resize(8,0);
    }
    
    for(size_t k = 0; k < nelem; k++){
        for (size_t m = 0; m < 8; m++){
            if( p(0 , t(m,k) ) == L1 ){
                T[k][m] = 1;
            }else{
                T[k][m] = 0;
            }
        }
    }
    vector<unsigned int> ii;
    vector < vector <unsigned int> > e1(11);
    for (size_t k = 0; k < 11; k++){
        e1[k].resize(Ny*Nz);
    }
    for (size_t k = 0; k < nelem; k++){
        order[k] = sort_indexes(T[k]);
        if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
         + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4){
            ii.push_back(k);
        }
    }    
    for (size_t jj = 0; jj < ii.size(); jj++){
        for (size_t m = 0; m < 4; m++){
            e1[m][jj] = t( order[ ii[jj] ][ m+4 ] , ii[jj] );
        }
    }    
    
    //back
    
    for(size_t k = 0; k < nelem; k++){
        for (size_t m = 0; m < 8; m++){
            if( p(0, t(m,k) ) == L2 ){
                T[k][m] = 1;
            }else{
                T[k][m] = 0;
            }
        }
    }
    ii.erase(ii.begin(),ii.end());
    vector < vector <unsigned int> > e2(11);
    for (size_t k = 0; k < 11; k++){
        e2[k].resize(Ny*Nz);
    }
    
    for (size_t k = 0; k < nelem; k++){
        order[k] = sort_indexes(T[k]);
        if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
         + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4){
            ii.push_back(k);
        }
    }    
    for (size_t jj = 0; jj < ii.size(); jj++){
        for (size_t m = 0; m < 4; m++){
            e2[m][jj] = t( order[ ii[jj] ][ m+4 ] , ii[jj] );
        }
    }
    
    //left
    
    for(size_t k = 0; k < nelem; k++){
        for (size_t m = 0; m < 8; m++){
            if( p(1, t(m,k) ) == H1 ){
                T[k][m] = 1;
            }else{
                T[k][m] = 0;
            }
        }
    }
    ii.erase(ii.begin(),ii.end());
    vector < vector <unsigned int> > e3(11);
    for (size_t k = 0; k < 11; k++){
        e3[k].resize(Nx*Nz);
    }
    
    for (size_t k = 0; k < nelem; k++){
        order[k] = sort_indexes(T[k]);
        if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
         + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4){
            ii.push_back(k);
        }
    }    
    for (size_t jj = 0; jj < ii.size(); jj++){
        for (size_t m = 0; m < 4; m++){
            e3[m][jj] = t( order[ ii[jj] ][ m+4 ] , ii[jj] );
        }
    }
    
    
    //right
    
    for(size_t k = 0; k < nelem; k++){
        for (size_t m = 0; m < 8; m++){
            if( p(1, t(m,k) ) == H2 ){
                T[k][m] = 1;
            }else{
                T[k][m] = 0;
            }
        }
    }
    ii.erase(ii.begin(),ii.end());
    vector < vector <unsigned int> > e4(11);
    for (size_t k = 0; k < 11; k++){
        e4[k].resize(Nx*Nz);
    }
    
    for (size_t k = 0; k < nelem; k++){
        order[k] = sort_indexes(T[k]);
        if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
         + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4){
            ii.push_back(k);
        }
    }    
    for (size_t jj = 0; jj < ii.size(); jj++){
        for (size_t m = 0; m < 4; m++){
            e4[m][jj] = t( order[ ii[jj] ][ m+4 ], ii[jj] );
        }
    }
    
    //bottom
    
    for(size_t k = 0; k < nelem; k++){
        for (size_t m = 0; m < 8; m++){
            if( p(2, t(m,k) ) == W1 ){
                T[k][m] = 1;
            }else{
                T[k][m] = 0;
            }
        }
    }
    ii.erase(ii.begin(),ii.end());
    vector < vector <unsigned int> > e5(11);
    for (size_t k = 0; k < 11; k++){
        e5[k].resize(Nx*Ny);
    }
    
    for (size_t k = 0; k < nelem; k++){
        order[k] = sort_indexes(T[k]);
        if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
         + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4){
            ii.push_back(k);
        }
    }    
    for (size_t jj = 0; jj < ii.size(); jj++){
        for (size_t m = 0; m < 4; m++){
            e5[m][jj] = t( order[ ii[jj] ][ m+4 ], ii[jj] );
        }
    }
    
    //top
    
    for(size_t k = 0; k < nelem; k++){
        for (size_t m = 0; m < 8; m++){
            if( p(2, t(m,k) ) == W2 ){
                T[k][m] = 1;
            }else{
                T[k][m] = 0;
            }
        }
    }
    ii.erase(ii.begin(),ii.end());
    vector < vector <unsigned int> > e6(11);
    for (size_t k = 0; k < 11; k++){
        e6[k].resize(Nx*Ny);
    }
    
    for (size_t k = 0; k < nelem; k++){
        order[k] = sort_indexes(T[k]);
        if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
         + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4){
            ii.push_back(k);
        }
    }    
    for (size_t jj = 0; jj < ii.size(); jj++){
        for (size_t m = 0; m < 4; m++){
            e6[m][jj] = t( order[ ii[jj] ][ m+4 ], ii[jj] );
        }
    }
    
    //building mesh.e field
    
    for(size_t k = 0; k < Ny*Nz; k++){
        for(size_t m = 0; m < 4; m++){
            e(m,k) = e1[m][k];
            e(m,k + Ny*Nz) = e2[m][k];
            e(m+4, k) = 0;
            e(m+4, k + Ny*Nz) = 0;
        }
        e(9,k) = sides[0];
        e(9,k + Ny*Nz) = sides[1];
        e(10,k) = region;
        e(10,k + Ny*Nz) = region;
    }
    for(size_t k = 2*Ny*Nz; k < 2*Ny*Nz + Nx*Nz; k++){
        for(size_t m = 0; m < 4; m++){
            e(m,k) = e3[m][k];
            e(m,k + Nx*Nz) = e4[m][k];
            e(m+4,k) = 0;
            e(m+4,k + Nx*Nz) = 0;
        }
        e(9,k) = sides[2];
        e(9,k + Nx*Nz) = sides[3];
        e(10,k) = region;
        e(10,k + Nx*Nz) = region;
    }
    for(size_t k = 2*(Ny*Nz + Nx*Nz); k < 2*(Ny*Nz + Nx*Nz) + Nx*Ny; k++){
        for(size_t m = 0; m < 4; m++){
            e(m,k) = e5[m][k];
            e(m,k + Nx*Ny) = e6[m][k];
        }
        e(9,k) = sides[4];
        e(9,k + Nx*Ny) = sides[5];
        e(10,k) = region;
        e(10,k + Nx*Ny) = region;
    }
    for(size_t k = 0; k < nelem; k++){
        t(8,k) = region;
    }
    
};

vector<unsigned int> Mesh::sort_indexes(const vector<unsigned int>& v){
    
    vector<unsigned int> idx;
    for(size_t k = 0; k < v.size(); k++){
        idx.push_back(k);
    }
    sort(idx.begin(),idx.end(), [&v](size_t i1, size_t i2){return v[i1] <= v[i2];} );
    
    return idx;
      
};


Mesh::~Mesh()
{
    delete [] p_data;
    delete [] e_data;
    delete [] t_data;
};

