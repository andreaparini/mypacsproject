/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   msh3m_nodes_on_faces.hpp
 * Author: pacs_student
 *
 * Created on December 26, 2016, 11:55 AM
 */

#ifndef MSH3M_NODES_ON_FACES_HPP
#define MSH3M_NODES_ON_FACES_HPP

#include <algorithm>
#include <vector>
#include "quadmesh.hpp"
using namespace std;

vector<int> msh3m_nodes_on_faces(Mesh& M, const vector<int> facelist){
    
    
    vector<int> facefaces;
    for (size_t k = 0; k < facelist.size(); k++){
        for (int j = 0; j < 2*M.Nx*M.Ny + 2*M.Ny*M.Nz + 2*M.Nx*M.Nz; j++){
            if (M.e(9,j) == facelist[k]){
                facefaces.push_back(j);
            }
        }
    }
    vector<int> facenodes;
    for (int i = 0; i < 4; i++){
        for (size_t j = 0; j < facefaces.size(); j++){
            facenodes.push_back(M.e(i, facefaces[j]));
        }
    }
    sort(facenodes.begin(), facenodes.end());
    auto last = unique(facenodes.begin(), facenodes.end());
    facenodes.erase(last, facenodes.end());
    
    
    
    return facenodes;
}

#endif /* MSH3M_NODES_ON_FACES_HPP */

