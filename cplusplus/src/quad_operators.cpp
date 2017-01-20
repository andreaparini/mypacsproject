/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "quad_operators.h"


void 
bim3a_advection_diffusion(quadmesh& mesh, 
                          vector<double> alpha, 
                          vector<double> beta, 
                          sparse_matrix &A)
{
    
  typedef vector< vector< vector<double> > > mat3d;
  mat3d Aloc(8, vector< vector<double> >(8, vector<double>(mesh.nelem)));
  
  vector< vector<double> > x(mesh.nelem, vector<double>(8));
  vector< vector<double> > y(mesh.nelem, vector<double>(8));
  vector< vector<double> > z(mesh.nelem, vector<double>(8));
  
  for(int k = 0; k < mesh.nelem; k++)
    {
      for(int m = 0; m < 8; m++)
        {
          x[k][m] = mesh.p(0, mesh.t(m,k) );
          y[k][m] = mesh.p(1, mesh.t(m,k) );
          z[k][m] = mesh.p(2, mesh.t(m,k) );            
        }
    }
    
  vector <double> psi12(mesh.nelem);
  vector <double> psi23(mesh.nelem);
  vector <double> psi43(mesh.nelem);
  vector <double> psi14(mesh.nelem);
  vector <double> psi15(mesh.nelem);
  vector <double> psi26(mesh.nelem);
  vector <double> psi37(mesh.nelem);
  vector <double> psi48(mesh.nelem);
  vector <double> psi56(mesh.nelem);
  vector <double> psi67(mesh.nelem);
  vector <double> psi87(mesh.nelem);
  vector <double> psi58(mesh.nelem);
  
  for (int k = 0; k < mesh.nelem; k++)
    {
      psi12[k] = beta[0]* (x[k][1] - x[k][0]) +
                 beta[1]* (y[k][1] - y[k][0]) +
                 beta[2]* (z[k][1] - z[k][0]);
      
      psi23[k] = beta[0]* (x[k][2] - x[k][1]) +
                 beta[1]* (y[k][2] - y[k][1]) +
                 beta[2]* (z[k][2] - z[k][1]);
      
      psi43[k] = beta[0]* (x[k][2] - x[k][3]) +
                 beta[1]* (y[k][2] - y[k][3]) +
                 beta[2]* (z[k][2] - z[k][3]);
      
      psi14[k] = beta[0]* (x[k][3] - x[k][0]) +
                 beta[1]* (y[k][3] - y[k][0]) +
                 beta[2]* (z[k][3] - z[k][0]);
      
      psi15[k] = beta[0]* (x[k][4] - x[k][0]) +
                 beta[1]* (y[k][4] - y[k][0]) +
                 beta[2]* (z[k][4] - z[k][0]);
      
      psi26[k] = beta[0]* (x[k][5] - x[k][1]) +
                 beta[1]* (y[k][5] - y[k][1]) +
                 beta[2]* (z[k][5] - z[k][1]);
      
      psi37[k] = beta[0]* (x[k][6] - x[k][2]) +
                 beta[1]* (y[k][6] - y[k][2]) +
                 beta[2]* (z[k][6] - z[k][2]);
      
      psi48[k] = beta[0]* (x[k][7] - x[k][3]) +
                 beta[1]* (y[k][7] - y[k][3]) +
                 beta[2]* (z[k][7] - z[k][3]);
      
      psi56[k] = beta[0]* (x[k][5] - x[k][4]) +
                 beta[1]* (y[k][5] - y[k][4]) +
                 beta[2]* (z[k][5] - z[k][4]);
      
      psi67[k] = beta[0]* (x[k][6] - x[k][5]) +
                 beta[1]* (y[k][6] - y[k][5]) +
                 beta[2]* (z[k][6] - z[k][5]);
      
      psi87[k] = beta[0]* (x[k][1] - x[k][0]) +
                 beta[1]* (y[k][1] - y[k][0]) +
                 beta[2]* (z[k][1] - z[k][0]);
      
      psi58[k] = beta[0]* (x[k][6] - x[k][7]) +
                 beta[1]* (y[k][6] - y[k][7]) +
                 beta[2]* (z[k][6] - z[k][7]);
  }
    
  vector <double> bp12(mesh.nelem);
  vector <double> bp23(mesh.nelem);
  vector <double> bp43(mesh.nelem);
  vector <double> bp14(mesh.nelem);
  vector <double> bp15(mesh.nelem);
  vector <double> bp26(mesh.nelem);
  vector <double> bp37(mesh.nelem);
  vector <double> bp48(mesh.nelem);
  vector <double> bp56(mesh.nelem);
  vector <double> bp67(mesh.nelem);
  vector <double> bp87(mesh.nelem);
  vector <double> bp58(mesh.nelem);
  
  vector <double> bm12(mesh.nelem);
  vector <double> bm23(mesh.nelem);
  vector <double> bm43(mesh.nelem);
  vector <double> bm14(mesh.nelem);
  vector <double> bm15(mesh.nelem);
  vector <double> bm26(mesh.nelem);
  vector <double> bm37(mesh.nelem);
  vector <double> bm48(mesh.nelem);
  vector <double> bm56(mesh.nelem);
  vector <double> bm67(mesh.nelem);
  vector <double> bm87(mesh.nelem);
  vector <double> bm58(mesh.nelem);
  
  for (size_t k = 0; k < mesh.nelem; k++)
    {
      bimu_bernoulli(psi12[k], bp12[k], bm12[k]);
      bimu_bernoulli(psi23[k], bp23[k], bm23[k]);
      bimu_bernoulli(psi43[k], bp43[k], bm43[k]);
      bimu_bernoulli(psi14[k], bp14[k], bm14[k]);
      bimu_bernoulli(psi15[k], bp15[k], bm15[k]);
      bimu_bernoulli(psi26[k], bp26[k], bm26[k]);
      bimu_bernoulli(psi37[k], bp37[k], bm37[k]);
      bimu_bernoulli(psi48[k], bp48[k], bm48[k]);
      bimu_bernoulli(psi56[k], bp56[k], bm56[k]);
      bimu_bernoulli(psi67[k], bp67[k], bm67[k]);
      bimu_bernoulli(psi87[k], bp87[k], bm87[k]);
      bimu_bernoulli(psi58[k], bp58[k], bm58[k]);
    }

  for (int k = 0; k < mesh.nelem; k++)
    {
      bp12[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
      bm12[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
      bp23[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
      bm23[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
      bp43[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
      bm43[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
      bp14[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
      bm14[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
      
      bp15[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
      bm15[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
      bp26[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
      bm26[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
      bp37[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
      bm37[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
      bp48[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
      bm48[k] *= alpha[k] * mesh.hx * mesh.hy / (4 * mesh.hz);
        
      bp56[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
      bm56[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
      bp67[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
      bm67[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
      bp87[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
      bm87[k] *= alpha[k] * mesh.hy * mesh.hz / (4 * mesh.hx);
      bp58[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
      bm58[k] *= alpha[k] * mesh.hx * mesh.hz / (4 * mesh.hy);
        
        
        
      Aloc[0][0][k] = bm12[k] + bm14[k] + bm15[k];
      Aloc[0][1][k] = -bp12[k];
      Aloc[0][3][k] = -bp14[k];
      Aloc[0][4][k] = -bp15[k];
        
      Aloc[1][0][k] = -bm12[k];
      Aloc[1][1][k] = bp12[k] + bm23[k] + bm26[k];
      Aloc[1][2][k] = -bp23[k];
      Aloc[1][5][k] = -bp26[k];
        
      Aloc[2][1][k] = -bm23[k];
      Aloc[2][2][k] = bp23[k] + bp43[k] + bm37[k];
      Aloc[2][3][k] = -bm43[k];
      Aloc[2][6][k] = -bp37[k];
        
      Aloc[3][0][k] = -bm14[k];
      Aloc[3][2][k] = -bp43[k];
      Aloc[3][3][k] = bm43[k] + bp14[k] + bm48[k];
      Aloc[3][7][k] = -bp48[k];
        
      Aloc[4][0][k] = -bm15[k];
      Aloc[4][4][k] = bp15[k] + bm56[k] + bm58[k];
      Aloc[4][5][k] = -bp56[k];
      Aloc[4][7][k] = -bp58[k];
        
      Aloc[5][1][k] = -bm26[k];
      Aloc[5][4][k] = -bm56[k];
      Aloc[5][5][k] = bp56[k] + bp26[k] + bm67[k];
      Aloc[5][6][k] = -bp67[k];
        
      Aloc[6][2][k] = -bm37[k];
      Aloc[6][5][k] = -bm67[k];
      Aloc[6][6][k] = bp37[k] + bp67[k] + bp87[k];
      Aloc[6][7][k] = -bm87[k];
        
      Aloc[7][3][k] = -bm48[k];
      Aloc[7][4][k] = -bm58[k];
      Aloc[7][6][k] = -bp87[k];
      Aloc[7][7][k] = bp48[k] + bp58[k] + bm87[k];
        
      for(int ii = 0; ii < 8; ii++)
        {
          for(int jj = 0; jj < 8; jj++)
            {
              A[ mesh.t(ii,k) ][ mesh.t(jj,k) ] += Aloc[ii][jj][k];
            }
        }        
    }    
}

void 
bim3a_rhs (quadmesh& mesh, 
           const vector<double> f, 
           const vector<double> g, 
           vector<double>& rhs)
{
  if(int(rhs.size()) < mesh.nnodes)
    {
      rhs.resize(mesh.nnodes,0);
    }
  int i;
    
  for(int k = 0; k < mesh.nelem; k++)
    {
      for(int m = 0; m < 8; m++)
        {
          i = mesh.t(m,k);
          rhs[i] += mesh.hx * mesh.hy * mesh.hz / 8 * f[k] * g[i];
        }
  }
}

void
bim3a_structure (const quadmesh& msh, sparse_matrix& SG)
{
  SG.resize (msh.nnodes);
  for (int iel = 0; iel < msh.nelem; ++iel)
    for (int inode = 0; inode < 4; ++inode)
      for (int jnode = 0; jnode < 4; ++jnode)
        {
          int ig = msh.t (inode, iel);
          int jg = msh.t (jnode, iel);
          SG[ig][jg] = 0.0;
        }
  SG.set_properties ();
};