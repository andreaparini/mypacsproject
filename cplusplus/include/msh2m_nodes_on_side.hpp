/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   msh2m_nodes_on_side.hpp
 * Author: pacs_student
 *
 * Created on December 13, 2016, 11:47 AM
 */

#ifndef MSH2M_NODES_ON_SIDE_HPP
#define MSH2M_NODES_ON_SIDE_HPP

#include <algorithm>
#include <vector>
#include "quadmesh.hpp"

bool isequal(int i, int j){
    return(i==j);
}

std::vector<int> msh2m_nodes_on_side(Mesh M, std::vector<int>sidelist){
    
    std::vector edgelist;
    for (size_t k = 0; k < sidelist.size(); k++){
        for (size_t j = 0; j < 2*M.Nx + 2*M.Ny; j++){
            if (M.e[5][j] == sidelist[k]){
                edgelist.push_back(j);
            }
        }
    }
    std::vector nodelist;
    for (size_t i = 0; i < 2; i++){
        for (size_t j = 0; j < edgelist.size(); j++){
            nodelist.push_back(M.e[i][edgelist(j)]);
        }
    }
    std::unique(nodelist.begin(), nodelist.end(), isequal);
    
    return nodelist;
}





#endif /* MSH2M_NODES_ON_SIDE_HPP */

