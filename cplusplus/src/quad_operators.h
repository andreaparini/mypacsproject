/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   quad_operators.h
 * Author: andrea
 *
 * Created on 19 gennaio 2017, 10.45
 */

#ifndef QUAD_OPERATORS_H
#define QUAD_OPERATORS_H

#include "quad_mesh.h"
#include "operators.h"
#include "bim_sparse.h"
#include <vector>

using namespace std;

void bim3a_advection_diffusion (quadmesh& mesh, 
                                vector<double> alpha, 
                                vector<double> beta, 
                                sparse_matrix &A);

void bim3a_rhs (quadmesh& mesh,
                const vector<double> f,
                const vector<double> g,
                vector<double>& rhs);

void
bim3a_structure (const quadmesh& msh,
                 sparse_matrix& A);

#endif /* QUAD_OPERATORS_H */

