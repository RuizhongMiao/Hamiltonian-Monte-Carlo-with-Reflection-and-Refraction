//  A simple C++ implementation of Hamiltonian MCMC
//  HMC.cpp
//
//  Created by Ruizhong Miao on 5/3/19.
//  Inspired by Radford Neal's "MCMC Using Hamiltonian Dynamics" url: https://arxiv.org/pdf/1206.1901.pdf
//  Copyright © 2019 Ruizhong Miao. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <random>
#include <vector>

using namespace std;

// Helper function: calculate the kinetic energy (squared norm of a numeric vector)
double kinetic(vector<double> p){
    double K = 0;
    for(int i=0; i<p.size(); i++){
        K = K + pow(p[i],2)/2;
    }
    return K;
}

// this function performs one iteration of HMC
vector<double> HMC(double (*U)(vector<double>), vector<double> (*grad_U)(vector<double>), double epsilon, int L, const vector<double> & current_q, default_random_engine & generator){
    
    uniform_real_distribution<double> u_dist(0.0,1.0);
    normal_distribution<double> n_dist(0.0, 1.0);
    
    vector<double> q = current_q;
    vector<double> p;
    vector<double> current_grad_U;
    
    unsigned long dim = q.size();
    int j,i;
    
    for(i=0; i<dim; i++){
        p.push_back(n_dist(generator));
    }
    vector<double> current_p = p;
    
    // Make a half step for momentum at the beginning
    current_grad_U = grad_U(q);
    for(i=0; i<dim; i++){
        p[i] = p[i] - epsilon * current_grad_U[i] / 2;
    }
    
    for(j=0; j<L; j++){
        for(i=0; i<dim; i++){
            // Make a full step for the position
            q[i] = q[i] + epsilon * p[i];
        }
        if(j != L-1){
            current_grad_U = grad_U(q);
            for(i=0; i<dim; i++){
                p[i] = p[i] - epsilon * current_grad_U[i];
            }
        }
    }
    
    // Make a half step for momentum at the end.
    current_grad_U = grad_U(q);
    for(i=0; i<dim; i++){
        p[i] = p[i] - epsilon * current_grad_U[i] / 2;
    }
    
    // Negate momentum at end of trajectory to make the proposal symmetric
    for(i=0; i<dim; i++){
        p[i] = -p[i];
    }
    
    double current_U = U(current_q);
    double proposed_U = U(q);
    double current_K = kinetic(current_p);
    double proposed_K = kinetic(p);
    
    if(u_dist(generator)<exp(current_U-proposed_U+current_K-proposed_K)){
        return q;
    }
    
    return current_q;
    
}


// this function performs one iteration of HMC
// No gradient of U is required
// Incorporated relection and refraction http://papers.nips.cc/paper/5801-reflection-refraction-and-hamiltonian-monte-carlo
vector<double> HMC(double (*U)(vector<double>), double epsilon, int L, const vector<double> & current_q, default_random_engine & generator){
    
    uniform_real_distribution<double> u_dist(0.0,1.0);
    normal_distribution<double> n_dist(0.0, 1.0);
    
    vector<double> q = current_q;
    vector<double> temp_q;
    vector<double> delta_q;
    vector<double> p;
    
    unsigned long dim = q.size();
    int j,i;
    double temp_U;
    double delta_U;
    
    for(i=0; i<dim; i++){
        p.push_back(n_dist(generator));
    }
    vector<double> current_p = p;
    
    for(j=0; j<L; j++){
        // Make a step for position and momentum
        temp_U = U(q);
        temp_q = q;
        for(i=0; i<dim; i++){
            delta_q = q;
            delta_q[i] = q[i] + epsilon * p[i];
            delta_U = U(delta_q) - temp_U;
            if( 2*delta_U > pow(p[i],2) ){
                p[i] = -p[i];
            }else{
                temp_q[i] = q[i] + epsilon * p[i] / 2;
                p[i] = sqrt( pow(p[i],2) - 2*delta_U ) * ((p[i]>0)-0.5)*2;
                temp_q[i] = q[i] + epsilon * p[i] / 2;
            }
        }
        q = temp_q;
    }
    
    // Negate momentum at end of trajectory to make the proposal symmetric
    for(i=0; i<dim; i++){
        p[i] = -p[i];
    }
    
    double current_U = U(current_q);
    double proposed_U = U(q);
    double current_K = kinetic(current_p);
    double proposed_K = kinetic(p);
    
    if(u_dist(generator)<exp(current_U-proposed_U+current_K-proposed_K)){
        return q;
    }
    
    return current_q;
    
}





