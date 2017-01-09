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

#include "../include/bim3a.hpp"
#include "../include/bim_sparse.hpp"

#include "../include/bim3a_advection_diffusion.hpp"
#include "../include/bim3a_rhs.hpp"

#include "../include/mumps_class.h"

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <functional>
#include <iterator>

/*
 * 
 */
using namespace std;

int main(int argc, char** argv) {
        
    
    // set mesh parameters
    
    const double     L1 = 0;
    const double     L2 = 1;
    const double     H1 = 0;
    const double     H2 = 1;    
    const double     W1 = 0;
    const double     W2 = 1;    
    const int        Nx = 20;
    const int        Ny = 20;
    const int        Nz = 20;
    
    // set equation parameters
    
    const double     alpha = 1.0;
    
    const double     beta1 = 0.0;
    const double     beta2 = 0.0;
    const double     beta3 = 0.0;
    
    
    
    const vector <double> beta = {beta1, beta2, beta3};
    
    
    const int        region = 1;
    
    // set default cube faces
    vector <int> sides;  //{left,right,front,back,bottom,top} // change to runtime defined
    sides.push_back(1);
    sides.push_back(2);
    sides.push_back(3);
    sides.push_back(4);
    sides.push_back(5);
    sides.push_back(6);
    
    // define mesh
    Mesh mesh(L1, L2, H1, H2, W1, W2, Nx, Ny, Nz, region, sides);
    
    // define alpha vector
    vector<double> alpha_values(mesh.nelem);
    for(int k = 0; k < mesh.nelem; k++){
        alpha_values[k] = alpha;
    }
    
    // define nodes coordinates
    vector<double> xnodes(mesh.nnodes);
    vector<double> ynodes(mesh.nnodes);
    vector<double> znodes(mesh.nnodes);
    
    for(int j = 0; j < mesh.nnodes; j++){
        xnodes[j] = mesh.p(0,j);
        ynodes[j] = mesh.p(1,j);
        znodes[j] = mesh.p(2,j);
    }
    
    // define forcing
    vector<double> f(mesh.nelem, 1.0);
    vector<double> g(mesh.nnodes);
    
    for(int j = 0; j < mesh.nnodes; j++){
        g[j] = 1;
        //g[j] = cos( 2*ynodes[j] ) * cos( znodes[j] ) *( cos(xnodes[j]) + 6*sin(xnodes[j]) );
    }
    
    // define dirichlet side list
    vector<int> dir_sidelist;
    dir_sidelist.push_back(1);
    dir_sidelist.push_back(2);
    
    // set Dirichlet nodes indexes
    vector<int> Dnodes = msh3m_nodes_on_faces(mesh, dir_sidelist);
    sort(Dnodes.begin(), Dnodes.end());
    
    vector<int> Nnodes;
    for(int k = 0; k < mesh.nnodes; k++){
        Nnodes.push_back(k);
    }
    
    // set the non border nodes
    vector<int> Varnodes;
    set_difference(Nnodes.begin(),Nnodes.end(), Dnodes.begin(), Dnodes.end(),
                   std::inserter(Varnodes, Varnodes.begin()) );
    
    // define A and rhs
    sparse_matrix A;
    bim3a_structure(mesh, A);
    
    bim3a_advection_diffusion(mesh, alpha_values, beta, A);
    cout << "MATRICE A"<< endl;
    cout << A;
    
    vector<double> rhs(mesh.nnodes,0);
    bim3a_rhs(mesh, f, g, rhs);
    cout << "VETTORE RHS" << endl;
    
    for(size_t k = 0; k< rhs.size(); k++){
        cout << k << " " << rhs[k] << endl;
    }
    
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

    for (size_t k = 0; k < rhs.size (); ++k){
        
        fout << rhs[k] << std::endl;
    }
    fout.close ();

    mumps_solver.cleanup ();

    return 0;
}

