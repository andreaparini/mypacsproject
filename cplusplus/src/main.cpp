/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: pacs_student
 *
 * Created on December 13, 2016, 6:46 PM
 */

#include "../include/quadmesh.hpp"
#include "../include/msh3m_nodes_on_faces.hpp"
#include "../include/GetPot.hpp"
#include "../muparser-2.2.5/include/muParser.h"
#include "../include/bim3a.hpp"
#include "../include/bim_sparse.hpp"
#include "../include/mumps_class.h"
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <functional>

/*
 * 
 */
using namespace std;

int main(int argc, char** argv) {
        
    
    GetPot cl(argc, argv);
    
    // set mesh parameters
    
    const double     L1 = cl("-L1", 0.0);
    const double     L2 = cl(2, "-L", "-L2", 1.0);
    const double     H1 = cl("-H1", 0.0);
    const double     H2 = cl(2, "-H", "-H2", 1.0);    
    const double     W1 = cl("-W1", 0.0);
    const double     W2 = cl(2, "-W", "-W2", 1.0);    
    const unsigned   Nx = cl("-Nx", 20);
    const unsigned   Ny = cl("-Ny", 20);
    const unsigned   Nz = cl("-Nz", 20);
    
    // set equation parameters
    
    const string     alpha_ = cl(2, "-a", "--alpha", "1.0");
    
    const double     beta1 = cl(2, "-b", "--beta", 1.0, 0);
    const double     beta2 = cl(2, "-b", "--beta", 0.0, 1);
    const double     beta3 = cl(2, "-b", "--beta", 0.0, 2);
    
    const vector <double> beta = {beta1, beta2, beta3};
    
    const string     f_ = cl("-f", "1.0");
    const string     g_ = cl("-g", "cos(2*y) * cos(z) *(cos(x) + 6*sin(x)");
    
    coefficient      alpha( alpha_ );
    coefficient      f(f_);
    coefficient      g(g_);
    
    
    
    // set Dirichlet sidelist
    
    const string     dir_sides_input = cl(3,"-d", "-s", "--sidelist", "1 2");
    stringstream     iss(dir_sides_input);

    const unsigned   side;
    vector<unsigned> dir_sidelist;

    while(iss >> side){
        dir_sidelist.push_back(side);
    }
    
    const unsigned  region = cl(2, "-r", "--region", 1);
    
    // set default cube faces
    vector <unsigned int> sides;  //{left,right,front,back,bottom,top} // change to runtime defined
    sides.push_back(1);
    sides.push_back(2);
    sides.push_back(3);
    sides.push_back(4);
    sides.push_back(5);
    sides.push_back(6);
    
    // define mesh
    Mesh mesh(L1, L2, H1, H2, W1, W2, Nx, Ny, Nz, region, sides);
    
    
    vector<double> f_values (mesh.nelem);
    vector<double> g_values (mesh.nnodes);
    vector<double> alpha_values (mesh.nelem);
    
    f.get_element_values(mesh, f, f_values);
    g.get_nodes_values(mesh, g, g_values);
    alpha.get_element_values(mesh, alpha, alpha_values);    
    
    // set Dirichlet nodes indexes
    vector<unsigned int> Dnodes = msh3m_nodes_on_faces(mesh, dir_sidelist);
    sort(Dnodes.begin(), Dnodes.end);
    
    vector<unsigned int> Nnodes(mesh.nelem,0);
    iota(Nnodes.begin(),Nnodes.end(),0);
    
    // set the non border nodes
    vector<unsigned int> Varnodes();
    set_difference(Nnodes.begin(),Nnodes.end(), Dnodes.begin(), Dnodes.end(), Varnodes.begin());
    
    // define A and rhs
    sparse_matrix A;
    bim3a_structure(mesh, A);
    
    bim3a_advection_diffusion(mesh, alpha_values, beta, A);
    
    vector<double> rhs(mesh.nelem,0);
    bim3a_rhs(mesh, f_values, g_values, rhs);
    
    //converting sparse matrix to aij format
    
    vector <double> vals;
    vector <int>    irow;
    vector <int>    jcol;
    A.aij(vals, irow, jcol, 1);
    
    mumps mumps_solver;
    
    mumps_solver.set_lhs_structure(A.rows(), irow, jcol);
    
    mumps_solver.analyze();
    
    mumps_solver.set_lhs_data(vals);
    
    mumps_solver.factorize();
    
    mumps_solver.set_rhs(rhs);
    
    mumps_solver.solve();
    
    std::cout << "\nResult of Stationary Test"
              << std::endl << "will be written in solution.txt"
              << std::endl;
    std::ofstream fout ("solution.txt");
    fout << std::endl;

    for (unsigned int k = 0; k < rhs.size (); ++k){
        
        fout << rhs[k] << std::endl;
    }
    fout.close ();

    mumps_solver.cleanup ();

    return 0;
}

