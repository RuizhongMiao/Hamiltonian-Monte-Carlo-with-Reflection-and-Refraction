//
//  main.cpp
//  BayesianProject
//
//  Created by Ruizhong Miao on 4/14/19.
//  Copyright © 2019 Ruizhong Miao. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "HMC.hpp"

using namespace std;

double U(vector<double> q){
    return pow(q[0],2)/2;
}

// Optional input to HMC function
vector<double> dU(vector<double> q){
    return q;
}


int main(int argc, const char * argv[]) {
    
    default_random_engine generator;
    generator.seed((unsigned)time(NULL));
    vector<double> q;
    ofstream myfile;
    
    q.push_back(0);
    cout<<q.size()<<endl;
    myfile.open("example.txt");
    for(int i=0; i<10000; i++){
        q = HMC(U, 0.001, 500, q, generator);
        myfile << q[0] << endl;
    }
    
    myfile.close();
    
    
    return 0;
}



