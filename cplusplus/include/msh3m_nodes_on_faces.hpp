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
bool isequal(int i, int j){
    return(i==j);
}

vector<unsigned int> msh3m_nodes_on_faces(Mesh M, const vector<unsigned int> facelist){
    
    
    vector<unsigned int> facefaces;
    for (size_t k = 0; k < facelist.size(); k++){
        for (size_t j = 0; j < 2*M.Nx + 2*M.Ny + 2*M.Nz; j++){
            if (M.e[9][j] == facelist[k]){
                facelist.push_back(j);
            }
        }
    }
    vector<unsigned int> facenodes;
    for (size_t i = 0; i < 4; i++){
        for (size_t j = 0; j < facelist.size(); j++){
            facenodes.push_back(M.e[i][facelist(j)]);
        }
    }
    unique(facenodes.begin(), facenodes.end(), isequal);
    
    return facenodes;
}

#endif /* MSH3M_NODES_ON_FACES_HPP */
