/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bim3a_rhs.hpp
 * Author: pacs_student
 *
 * Created on December 26, 2016, 10:59 AM
 */

#ifndef BIM3A_RHS_HPP
#define BIM3A_RHS_HPP

#include "quadmesh.hpp"
#include "bim3a.hpp"
#include "../muparser-2.2.5/include/muParser.h"
#include <functional>
#include <vector>

using namespace std;


void bim3a_rhs
(Mesh& mesh, const vector<double> f , const vector<double> g, vector<double>& rhs);
    

#endif /* BIM3A_RHS_HPP */

