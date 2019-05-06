//
//  HMC.hpp
//  BayesianProject
//
//  Created by Ruizhong Miao on 4/15/19.
//  Copyright © 2019 Ruizhong Miao. All rights reserved.
//

#ifndef HMC_hpp
#define HMC_hpp

// #include <vector>
// #include <random>

using namespace std;

double kinetic(vector<double> p);

vector<double> HMC(double (*U)(vector<double>), vector<double> (*grad_U)(vector<double>), double epsilon, int L, const vector<double> & current_q, default_random_engine & generator);

vector<double> HMC(double (*U)(vector<double>), double epsilon, int L, const vector<double> & current_q, default_random_engine & generator);


#endif /* HMC_hpp */
