/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bernoulli.hpp
 * Author: pacs_student
 *
 * Created on December 13, 2016, 12:21 PM
 */

#ifndef BERNOULLI_HPP
#define BERNOULLI_HPP

#include <vector>
#include <cmath>
#include <limits>
using namespace std;

void bernoulli (vector<double>& bp, vector<double>& bn, vector<double> x){
    
    bp.resize(x.size(),0);
    bn.resize(x.size(),0);
    double xlim = 0.01;
    vector<double> absx(x.size());
    for (size_t i = 0; i < absx.size(); i++){
        absx[i] = abs(x[i]);
    }
    
    double bp_tmp;
    double bn_tmp;
    for (size_t i = 0; i < x.size(); i++){
        if(x[i] == 0) {
            bp_tmp = 1;
            bn_tmp = 1;
        }
        if(absx[i] > 80 && x[i] > 0 ){
            bp_tmp = 0;
            bn_tmp = x[i];
        }
        if(absx[i] > 80 && x[i] < 0 ){
            bp_tmp = - x[i];
            bn_tmp = 0;
        }
        if(absx[i] <= 80 && absx[i] > xlim){
            bp_tmp = x[i]/( exp(x[i]) - 1 );
            bn_tmp = x[i] + bp_tmp;
        }
        if(absx[i] > 0 && absx[i] <= xlim){
            double jj = 1;
            double fp = 1;
            double fn = 1;
            double df = 1;
            double segno = 1;
            while(df > std::numeric_limits<double>::epsilon()){
                jj++;
                segno = -segno;
                df = df * x[i] / jj;
                fp += df;
                fn += segno * df;
            }
            bp_tmp = 1/fp;
            bn_tmp = 1/fn;
        }
        bp[i] = bp_tmp;
        bn[i] = bn_tmp;
    }
    
    
    
}

#endif /* BERNOULLI_HPP */

