/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   quadmesh.hpp
 * Author: pacs_student
 *
 * Created on December 13, 2016, 11:50 AM
 */

#ifndef QUADMESH_HPP
#define QUADMESH_HPP

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

class Mesh
{
public:
    
    const double L;
    const double H;
    const double W;
    const int Nx;
    const int Ny;
    const int Nz;
    const double hx;
    const double hy;
    const double hz;
    const int nnodes;
    const int nelem;
    double *p_data;
    int *e_data;
    int *t_data;
    const int region;
    vector<int> sides;
    
    Mesh
    (const double L1, const double L2, const double H1, const double H2,
     const double W1, const double W2, const int Nx, const int Ny,
     const int Nz,    const int region, vector<int> sides);
    
    ~Mesh();
    
    // accessing p data
    inline double&
    p(int idir, int inode){
        return (*(p_data + idir + 3*inode));
    };
    
    // accessing t data
    inline int&
    t(int inode, int iel){
        return (*(t_data + inode + 9*iel));
    };
    
    // accessing e data
    inline int&
    e(int i, int iface){
        return (*(e_data + i + 11*iface));
    };
    
    // accessing p data (const version)
    inline const double&
    p(int idir, int inode) const{
        return (*(p_data + idir + 3*inode));
    };
    
    // accessing t data (const version)
    inline const int&
    t(int inode, int iel) const{
        return (*(t_data + inode + 9*iel));
    };
    
    // accessing e data (const version
    inline const int&
    e(int i, int iface) const{
        return (*(e_data + i + 11*iface));
    };
    
    // takes a vector v as input and returns a vector w in which 
    // w[i] contains the index of the position j in which v[i] 
    // would be in the sorted v vector
    vector <int> sort_indexes(const vector<int> &v);
    
    
};
#endif /* QUADMESH_HPP */

