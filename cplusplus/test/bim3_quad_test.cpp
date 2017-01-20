/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   quad_main.cpp
 * Author: pacs_student
 *
 * Created on December 13, 2016, 6:46 PM
 */




#include <quad_mesh.h>
#include <quad_operators.h>
#include <mpi.h>
#include <mumps_class.h>
#include <bim_sparse.h>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <functional>
#include <iterator>
#include <bim_config.h>
/*
 * 
 */
using namespace std;

int main(int argc, char** argv) {
  
  
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  
  
  vector <double> vals;
  vector <int>    irow;
  vector <int>    jcol;
  sparse_matrix   A;
  vector<double>  rhs;
      
  if(rank == 0)
    {
      // set mesh parameters
    
      const double     L1 = 0;
      const double     L2 = 1;
      const double     H1 = 0;
      const double     H2 = 1;    
      const double     W1 = 0;
      const double     W2 = 1;    
      const int        Nx = 20;
      const int        Ny = 10;
      const int        Nz = 10;
    
      // set equation parameters
    
      const double     alpha = 1.0;
    
      const double     beta1 = 0.0;
      const double     beta2 = 0.0;
      const double     beta3 = 0.0;
    
      const vector <double> beta = {beta1, beta2, beta3};
    
      const int        region = 1;
    
      // set default cube faces
      vector <int> sides;  //{left,right,front,back,bottom,top} // change to runtime defined
      sides.push_back(1);  //left
      sides.push_back(2);  //right
      sides.push_back(3);  //front
      sides.push_back(4);  //back
      sides.push_back(5);  //bottom
      sides.push_back(6);  //top
    
      // define mesh
      quadmesh mesh(L1, L2, H1, H2, W1, W2, Nx, Ny, Nz, region, sides);
    
      // define alpha vector
      vector<double> alpha_values(mesh.nelem);
      for(int k = 0; k < mesh.nelem; k++)
        {
          alpha_values[k] = alpha;
        }
    
        // define nodes coordinates
      vector<double> xnodes(mesh.nnodes);
      vector<double> ynodes(mesh.nnodes);
      vector<double> znodes(mesh.nnodes);
    
      for(int j = 0; j < mesh.nnodes; j++)
        {
          xnodes[j] = mesh.p(0,j);
          ynodes[j] = mesh.p(1,j);
          znodes[j] = mesh.p(2,j);
        }
    
        // define forcing
      vector<double> f(mesh.nelem, 1.0);
      vector<double> g(mesh.nnodes);
    
      for(int j = 0; j < mesh.nnodes; j++)
        {
          g[j] = 0.0;
          //g[j] = cos( 2*ynodes[j] ) * cos( znodes[j] ) *( cos(xnodes[j]) + 6*sin(xnodes[j]) );
        }
    
      // define dirichlet side list
      vector<int> dir_sidelist;
      dir_sidelist.push_back(5);
      dir_sidelist.push_back(6);
    
      // set Dirichlet nodes indexes
      vector<int> Dnodes;
      mesh.nodes_on_faces(dir_sidelist, Dnodes);
    
      // set boundary conditions
      vector<double> Vnodes(Dnodes.size(),0);
      for (size_t j = 0; j < Vnodes.size(); j++)
        {
          if(mesh.p(2,Dnodes[j]) == W1)
            {
              Vnodes[j] = 0;
            }
          if(mesh.p(2,Dnodes[j]) == W2)
            {
              Vnodes[j] = 1;
            }
        }
        
        // define A and rhs
      bim3a_structure(mesh, A);
      
      // set matrix and rhs
      bim3a_advection_diffusion(mesh, alpha_values, beta, A);
        
      bim3a_rhs(mesh, f, g, rhs);
        
      bim3a_dirichlet_bc(A, rhs, Dnodes, Vnodes);
    }
  //converting sparse matrix to aij format
    
  
  A.aij(vals, irow, jcol, 1);
    
  // solving problem
    
  mumps mumps_solver;
  
  if(rank == 0)
    mumps_solver.set_lhs_structure(A.rows(), irow, jcol);
  
  mumps_solver.analyze();
  
  if(rank == 0)
    mumps_solver.set_lhs_data(vals);
  
  mumps_solver.factorize();
  
  if(rank == 0)
    mumps_solver.set_rhs(rhs);
  
  mumps_solver.solve();
  
  if(rank == 0)
    {
      std::cout << "\nResult of Stationary Test"
                << std::endl << "will be written in solution.txt"
                << std::endl;
      std::ofstream fout ("solution.txt");
      fout << std::endl;
        
      for (size_t k = 0; k < rhs.size (); ++k)
        {
          fout << rhs[k] << std::endl;
        }
        
      fout.close ();
    }
  
  mumps_solver.cleanup ();
  MPI_Finalize();
  return 0;
}

