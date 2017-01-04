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

#include "../muparser-2.2.5/include/muParser.h"
#include "quadmesh.hpp"
#include "bim_sparse.hpp"

#include <vector>
#include <string>


using namespace std;
using namespace mu;



class coefficient
{
private:
    string expr;
    double var_x;
    double var_y;
    double var_z;
    Parser p;
    
public:
    coefficient(const string &s) : expr(s){
        try{
            p.DefineVar("x", &var_x);
            p.DefineVar("y", &var_y);
            p.DefineVar("z", &var_z);
            p.SetExpr(expr.c_str());
        }catch(mu::Parser::exception_type &e){
            std::cerr << e.GetMsg() << std::endl;
        }
    };
    double operator() (double x, double y, double z){
        double retval;
        try{
            retval = p.Eval();
        }catch(mu::Parser::exception_type &e){
            std::cerr << e.GetMsg() << std::endl;
        }
        return retval;
    };
    void get_element_values (Mesh mesh, vector<double>& f_values){
        
        vector< vector <double> > xk(mesh.nelem);
        vector< vector <double> > yk(mesh.nelem);
        vector< vector <double> > zk(mesh.nelem);
    
        for(size_t m = 0; m < mesh.nelem; m++){
            xk[m].resize(8,0.0);
            yk[m].resize(8,0.0);
            zk[m].resize(8,0.0);
        }
        
        for(size_t k = 0; k < mesh.nelem; k++){
            for(size_t m = 0; m < 8; m++){
                xk[m][k] = mesh.p(0,mesh.t(m,k));
                yk[m][k] = mesh.p(1,mesh.t(m,k));
                zk[m][k] = mesh.p(2,mesh.t(m,k));            
            }
        }
        for(size_t k = 0; k < mesh.nelem; k++){
            f_values[k] = this->operator ()( (xk[0][k] + xk[1][k])/2,
                                             (yk[0][k] + yk[2][k])/2,
                                             (zk[0][k] + zk[4][k])/2  );
        }
       
    };
    
    void get_nodes_values (Mesh mesh, vector<double>& g_values){
        
        vector<double> x(mesh.nnodes);
        vector<double> y(mesh.nnodes);
        vector<double> z(mesh.nnodes);
    
        //evaluate f and g coefficients
        for (size_t j = 0; j < mesh.nnodes; j++){
            x[j] = mesh.p[0][j];
            y[j] = mesh.p[1][j];
            z[j] = mesh.p[2][j];
        }
        
        for(size_t j = 0; j < mesh.nnodes; j++){
            g_values[j] = this->operator ()( x[j], y[j], z[j]);
        }
    };

};

void
bim3a_structure (const Mesh& msh, sparse_matrix& SG)
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


#endif /* BIM3A_HPP */

