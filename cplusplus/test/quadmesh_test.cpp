/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   quadmesh_test.cpp
 * Author: pacs_student
 *
 * Created on January 3, 2017, 5:07 PM
 */

#include "../include/quadmesh.hpp"
#include "../include/GetPot.hpp"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

/*
 * Simple C++ Test Suite
 * Testing quadmesh.cpp function
 */
void test(double L1, double L2, double H1, double H2, double W1, double W2,
           int Nx, int Ny, int Nz, int region, vector <int> sides) {
    
    Mesh mesh(L1, L2, H1, H2, W1, W2, Nx, Ny, Nz, region, sides);
        
    std::ofstream fout ("mesh.txt");
    fout << std::endl;

    std::cout << "printing to mesh.txt file" << std::endl;
    std::cout << "printing mesh nodes coordinates" << std::endl;
        
    for (int j = 0; j < mesh.nnodes; j++){
        
        fout << mesh.p(0,j) << "\t" << mesh.p(1,j) << "\t" << mesh.p(2,j) << "\t" << std::endl;
    }
        
    std::cout << "printing cubes vertexes" << std::endl;
        
    for (int k = 0; k < mesh.nelem; k++){
    
        fout << k << "\t" << mesh.t(0,k) << "\t" << mesh.t(1,k) << "\t" << mesh.t(2,k) << "\t" << mesh.t(3,k) << "\t"
                          << mesh.t(4,k) << "\t" << mesh.t(5,k) << "\t" << mesh.t(6,k) << "\t" << mesh.t(7,k) << std::endl;
    }
        
    std::cout << "printing border faces nodes indexes" << std::endl;
        
    for (int i = 0; i < 2*(mesh.Nx*mesh.Ny + mesh.Ny*mesh.Nz + mesh.Nx*mesh.Nz); i++){
        
        fout << i << "\t" << mesh.e(0,i) << "\t" << mesh.e(1,i) << "\t" << mesh.e(2,i) << "\t" << mesh.e(3,i) 
                          << "\t side = " << mesh.e(9,i) << std::endl;
    }
        
    fout.close ();
    

}



int main(int argc, char** argv) {
    
    GetPot cl(argc, argv);
    
    
    std::cout << " -----quadmesh_test----- " << std::endl;
    std::cout << " -h, --help to get test information" << std::endl;
    
    
    
    if (cl.search(2, "-h", "--help") ){
        std::cout << "-----quadmesh_test----- " << std::endl;
        std::cout << "test modes:" << std::endl;
        std::cout << "-h, --help: prints this help message" << std::endl;
        std::cout << "-t1, --test1: creates a uniform cube mesh with unit edge and h = 0.1" << std::endl;
        std::cout << "-t2, --test2: creates a mesh with these sizes: L = 2, H = 1, W = 0.5; hx = 0.1; hy = 0.1; hz = 0.1" << std::endl;
        
    }
    
    if(cl.search(2, "-t1", "--test1") ){
        
        std::cout << "-----quadmesh_test----- " << std::endl;
        std::cout << "starting test 1" << std::endl;
        std::cout << "creating cube mesh, unit edge, 10 elements per side" << std::endl;
        
        const double     L1 = 0;
        const double     L2 = 1;
        const double     H1 = 0;
        const double     H2 = 1;
        const double     W1 = 0;
        const double     W2 = 1;
        const int        Nx = 10;
        const int        Ny = 10;
        const int        Nz = 10;
        const int        region = 1;
        vector<int> sides;
        sides.push_back(1);
        sides.push_back(2);
        sides.push_back(3);
        sides.push_back(4);
        sides.push_back(5);
        sides.push_back(6);
        
        test(L1, L2, H1, H2, W1, W2, Nx, Ny, Nz, region, sides);
        
        
    }
    
    if(cl.search(2, "-t2", "--test2") ){
        
        std::cout << "-----quadmesh_test----- " << std::endl;
        std::cout << "starting test 2" << std::endl;
        std::cout << "creating a mesh with these sizes: L = 2, H = 1, W = 0.5; hx = 0.1; hy = 0.1; hz = 0.1" << std::endl;
        
        const double     L1 = 0;
        const double     L2 = 2;
        const double     H1 = 0;
        const double     H2 = 1;
        const double     W1 = 0;
        const double     W2 = 0.5;
        const int        Nx = 20;
        const int        Ny = 10;
        const int        Nz = 5;
        const int        region = 1;
        vector<int> sides;
        sides.push_back(1);
        sides.push_back(2);
        sides.push_back(3);
        sides.push_back(4);
        sides.push_back(5);
        sides.push_back(6);
        
        test(L1, L2, H1, H2, W1, W2, Nx, Ny, Nz, region, sides);
        
    }
    
    return 0;
}

