/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "quad_mesh.h"

using namespace std;

quadmesh::quadmesh
(const double L1, const double L2, const double H1, const double H2,
 const double W1, const double W2, const int Nx_, const int Ny_, const int Nz_,
 const int region_, vector<int> sides_):
 L (L2 - L1),
 H (H2 - H1),
 W (W2 - W1),
 Nx (Nx_),
 Ny (Ny_),
 Nz (Nz_),
 hx ( L / double(Nx) ),
 hy ( H / double(Ny) ),
 hz ( W / double(Nz) ),
 nnodes ( (Nx+1)*(Ny+1)*(Nz+1) ),
 nelem ( Nx * Ny *Nz ),
 region (region_),
 sides(sides_)
{
  p_data = new double [3*nnodes];
  e_data = new int [11*2*(Nx*Ny + Ny*Nz + Nx*Nz)];
  t_data = new int [9*nelem];
    
  // generate points
  for (int k = 0; k < Nz + 1; k++)
    {
    for (int i = 0; i < Nx + 1; i++)
      {
      for (int j = 0; j < Ny + 1; j++)
        {
          p(0, k*(Nx+1)*(Ny+1) + i*(Ny+1) + j) = i * hx;
          p(1, k*(Nx+1)*(Ny+1) + i*(Ny+1) + j) = j * hy;
          p(2, k*(Nx+1)*(Ny+1) + i*(Ny+1) + j) = k * hz;
        }
      }
    }
    
  int count = 0;
  // creating iiv
  vector < vector < vector < int > > > iiv_c(Nz+1);
  for( int k = 0; k < Nz + 1; k++ )
    {
      iiv_c[k].resize(Ny+1);
      for( int i = 0; i < Ny + 1; i++)
      {
        iiv_c[k][i].resize(Nx + 1);
      }
    }
    // initializing iivc
  for(int k = 0; k < Nz + 1; k++ )
    {
      for(int j = 0; j < Nx + 1; j++)
        {
          for(int i = 0; i < Ny +1; i++)
            {
              iiv_c[k][i][j] = count;
              count++;
            }
        }
    }
  
  for( int k = 0; k < Nz + 1; k++ )
    {
      for( int i  = 0; i < Ny + 1; i++)
        {
          iiv_c[k][i].pop_back();
        }
      iiv_c[k].pop_back();
    }
  iiv_c.pop_back();
  
  vector <int> iiv;
  count = 0;
  for (int k = 0; k < Nz; k++)
    {
      for (int j = 0; j < Nx; j++)
        {
          for (int i = 0; i < Ny; i++)
            {
              iiv.push_back(iiv_c[k][i][j]);
            }   
        }    
    }

  for (int k = 0; k < nelem; k++)
    {
      t(0,k) = iiv[k];
      t(1,k) = iiv[k] + Ny + 1;
      t(2,k) = iiv[k] + Ny + 2;
      t(3,k) = iiv[k] + 1;
      t(4,k) = iiv[k] + (Nx+1)*(Ny+1);
      t(5,k) = iiv[k] + (Nx+1)*(Ny+1) + Ny + 1;
      t(6,k) = iiv[k] + (Nx+1)*(Ny+1) + Ny + 2;
      t(7,k) = iiv[k] + (Nx+1)*(Ny+1) + 1;
    }
  
  //generate boundary face list
  
  //front
  
  vector< vector<int> > T(nelem);
  vector< vector<int> > order(nelem);
  
  for (int k = 0; k < nelem; k++)
    {
      order[k].resize(8,0);
      T[k].resize(8,0);
   }
  
  for(int k = 0; k < nelem; k++)
    {
      for (int m = 0; m < 8; m++)
        {
          if( p(0 , t(m,k) ) == L1 )
            {
              T[k][m] = 1;
            }
          else
            {
              T[k][m] = 0;
            }
        }
    }
  vector<int> ii;
  vector< vector<int> > e1(11);
  for (int k = 0; k < 11; k++)
    {
      e1[k].resize(Ny*Nz);
    }
  for (int k = 0; k < nelem; k++)
    {
      order[k] = sort_indexes(T[k]);
      if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
       + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4)
        {
          ii.push_back(k);
        }
    }    
  for (size_t jj = 0; jj < ii.size(); jj++)
    {
      for (int m = 0; m < 4; m++)
        {
          e1[m][jj] = t( order[ ii[jj] ][ m+4 ] , ii[jj] );
        }
    }    
   
  //back
  
  for(int k = 0; k < nelem; k++)
    {
      for (int m = 0; m < 8; m++)
        {
          if( p(0, t(m,k) ) == L2 )
            {
              T[k][m] = 1;
            }
          else
            {
              T[k][m] = 0;
            }
        }
    }
  ii.erase(ii.begin(),ii.end());
  vector< vector<int> > e2(11);
  for (int k = 0; k < 11; k++)
    {
      e2[k].resize(Ny*Nz);
    }
  
  for (int k = 0; k < nelem; k++)
    {
    order[k] = sort_indexes(T[k]);
      if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
       + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4)
        {
          ii.push_back(k);
        }
    }    
  for (size_t jj = 0; jj < ii.size(); jj++)
    {
      for (int m = 0; m < 4; m++)
        {
          e2[m][jj] = t( order[ ii[jj] ][ m+4 ] , ii[jj] );
        }
    }
  
  //left
   
  for(int k = 0; k < nelem; k++)
    {
      for (int m = 0; m < 8; m++)  
         {
          if( p(1, t(m,k) ) == H1 )
            {
              T[k][m] = 1;
            }
          else
            {
              T[k][m] = 0;
            }
        }
    }
  ii.erase(ii.begin(),ii.end());
  vector < vector <int> > e3(11);
  for (int k = 0; k < 11; k++)
    {
      e3[k].resize(Nx*Nz);
    }
  
  for (int k = 0; k < nelem; k++)
    {
      order[k] = sort_indexes(T[k]);
      if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
       + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4)
        {
          ii.push_back(k);
        }
    }    
  for (size_t jj = 0; jj < ii.size(); jj++)
    {
      for (int m = 0; m < 4; m++)
        {
          e3[m][jj] = t( order[ ii[jj] ][ m+4 ] , ii[jj] );
        }
    }
  
  
  //right
  
  for(int k = 0; k < nelem; k++)
    {
      for (int m = 0; m < 8; m++)
        {
          if( p(1, t(m,k) ) == H2 )
            {
              T[k][m] = 1;
            }
          else
            {
              T[k][m] = 0;
            }
        }
    }
  ii.erase(ii.begin(),ii.end());
  vector < vector <int> > e4(11);
  for (int k = 0; k < 11; k++)
    {
      e4[k].resize(Nx*Nz);
    }
  
  for (int k = 0; k < nelem; k++)
    {
      order[k] = sort_indexes(T[k]);
      if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
       + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4)
        {
          ii.push_back(k);
        }
    }     
 for (size_t jj = 0; jj < ii.size(); jj++)
    {
      for (int m = 0; m < 4; m++)
        {
          e4[m][jj] = t( order[ ii[jj] ][ m+4 ], ii[jj] );
        }
    }
  
  //bottom
  
  for(int k = 0; k < nelem; k++)
    {
      for (int m = 0; m < 8; m++)
        {
          if( p(2, t(m,k) ) == W1 )
            {
              T[k][m] = 1;
            }
          else
            {
              T[k][m] = 0;
            }
        }
    }
  ii.erase(ii.begin(),ii.end());
  vector < vector <int> > e5(11);
  for (int k = 0; k < 11; k++)
    {
      e5[k].resize(Nx*Ny);
    }
  
  for (int k = 0; k < nelem; k++)
    {
      order[k] = sort_indexes(T[k]);
      if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
       + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4)
        {
          ii.push_back(k);
        }
    }    
  for (size_t jj = 0; jj < ii.size(); jj++)
    {
      for (int m = 0; m < 4; m++)
        {
          e5[m][jj] = t( order[ ii[jj] ][ m+4 ], ii[jj] );
        }
    }
  
  //top
  
  for(int k = 0; k < nelem; k++)
    {
      for (int m = 0; m < 8; m++)
        {
          if( p(2, t(m,k) ) == W2 )
            {
              T[k][m] = 1;
            }
          else 
            {
              T[k][m] = 0;
            }
        }
    }
  ii.erase(ii.begin(),ii.end());
  vector < vector <int> > e6(11);
  for (int k = 0; k < 11; k++)
    {
      e6[k].resize(Nx*Ny);
    }
  
  for (int k = 0; k < nelem; k++)
    {
      order[k] = sort_indexes(T[k]);
      if(T[k][0] + T[k][1] + T[k][2] + T[k][3]
       + T[k][4] + T[k][5] + T[k][6] + T[k][7] == 4)
        {
          ii.push_back(k);
        }
    }    
  for (size_t jj = 0; jj < ii.size(); jj++)
    {
      for (int m = 0; m < 4; m++)
        {
          e6[m][jj] = t( order[ ii[jj] ][ m+4 ], ii[jj] );
        }
    }
  
  //building mesh.e field
  
  for(int k = 0; k < Ny*Nz; k++)
    {
      for(int m = 0; m < 4; m++)
        {
          e(m,k) = e1[m][k];
          e(m,k + Ny*Nz) = e2[m][k];
          e(m+4, k) = 0;
          e(m+4, k + Ny*Nz) = 0;
        }
      e(9,k) = sides[0];
      e(9,k + Ny*Nz) = sides[1];
      e(10,k) = region;
      e(10,k + Ny*Nz) = region;
    }
  for(int k = 0; k < Nx*Nz; k++)
    {
      for(int m = 0; m < 4; m++)  
        {
          e(m,k + 2*Ny*Nz) = e3[m][k];
          e(m,k + 2*Ny*Nz + Nx*Nz) = e4[m][k];
          e(m+4,k + 2*Ny*Nz) = 0;
          e(m+4,k + 2*Ny*Nz+ Nx*Nz) = 0;
        }
      e(9,k + 2*Ny*Nz) = sides[2];
      e(9,k + 2*Ny*Nz + Nx*Nz) = sides[3];
      e(10,k + 2*Ny*Nz) = region;
      e(10,k + 2*Ny*Nz + Nx*Nz) = region;
    }
  for(int k = 0; k < Nx*Ny; k++)
    {
      for(int m = 0; m < 4; m++)
        {
          e(m,k + 2*Ny*Nz + 2*Nx*Nz) = e5[m][k];
          e(m,k + 2*Ny*Nz + 2*Nx*Nz + Nx*Ny) = e6[m][k];
          e(m+4,k + 2*Ny*Nz + 2*Nx*Nz) = 0;
          e(m+4,k + 2*Ny*Nz  + 2*Nx*Nz+ Nx*Ny) = 0;
        }
      e(9,k + 2*Ny*Nz + 2*Nx*Nz) = sides[4];
      e(9,k + 2*Ny*Nz + 2*Nx*Nz + Nx*Ny) = sides[5];
      e(10,k + 2*Ny*Nz + 2*Nx*Nz) = region;
      e(10,k + 2*Ny*Nz + 2*Nx*Nz + Nx*Ny) = region;
    }
  for(int k = 0; k < nelem; k++)
    {
      t(8,k) = region;
    }
  
};

vector<int> quadmesh::sort_indexes
(const vector<int>& v)
{
    
  vector<int> idx;
    for(size_t k = 0; k < v.size(); k++)
      {
        idx.push_back(k);
      }
  sort(idx.begin(),idx.end(), [&v](int i1, int i2){return v[i1] <= v[i2];} );
  
  return idx;
      
};


void quadmesh::nodes_on_faces
(const vector<int> facelist, vector<int> &facenodes)
{
  vector<int> facefaces;
  for (size_t k = 0; k < facelist.size(); k++)
    {
      for (int j = 0; j < 2*this->Nx * this->Ny + 2 * this->Ny * this->Nz 
                        + 2*this->Nx * this->Nz; j++)
        {
          if (this->e(9,j) == facelist[k])
            {
              facefaces.push_back(j);
            }
        }
    }
  
  for (int i = 0; i < 4; i++)
    {
      for (size_t j = 0; j < facefaces.size(); j++)
        {
          facenodes.push_back(this->e(i, facefaces[j]));
        }
    }
  sort(facenodes.begin(), facenodes.end());
  auto last = unique(facenodes.begin(), facenodes.end());
  facenodes.erase(last, facenodes.end());
 
};


quadmesh::~quadmesh()
{
  delete [] p_data;
  delete [] e_data;
  delete [] t_data;
};

