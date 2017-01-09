/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bim3a_advection_diffusion.hpp
 * Author: pacs_student
 *
 * Created on December 26, 2016, 12:15 PM
 */

#ifndef BIM3A_ADVECTION_DIFFUSION_HPP
#define BIM3A_ADVECTION_DIFFUSION_HPP

#include "quadmesh.hpp"
#include "bim3a.hpp"
#include "bim_sparse.hpp"
#include <vector>

using namespace std;

void bim3a_advection_diffusion
(Mesh& mesh, vector<double> alpha, vector<double> beta, sparse_matrix &A);;

#endif /* BIM3A_ADVECTION_DIFFUSION_HPP */

