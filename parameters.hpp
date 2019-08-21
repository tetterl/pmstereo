#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

/// define common parameters, variables, ...

// max_disp is a global variable as it can be set from pm.cpp by 
// setting the cl argument --maxdisp
int max_disp = 64;
#define WINDOW_SIZE 35
#define ALPHA 0.9f
#define ONEMINUSALPHA 0.1f
#define GAMMA 10.f
#define GAMMA_INV 0.1f
#define TAUCOL 10.f
#define TAUGRAD 2.f

#endif // PARAMETERS_HPP
