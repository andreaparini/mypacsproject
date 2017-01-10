/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bim3a.hpp
 * Author: pacs_student
 *
 * Created on December 28, 2016, 10:01 PM
 */

#ifndef BIM3A_HPP
#define BIM3A_HPP

#include "quadmesh.hpp"
#include "bim_sparse.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <limits>

using namespace std;

void
bim3a_structure (const Mesh& msh, sparse_matrix& SG);


void 
bernoulli (vector<double>& bp, vector<double>& bn, vector<double> x);

void
bim3a_dirichlet_bc (sparse_matrix& M,
                    std::vector<double>& b,
                    const std::vector<int>& bnodes,
                    const std::vector<double>& vnodes);



#endif /* BIM3A_HPP */

