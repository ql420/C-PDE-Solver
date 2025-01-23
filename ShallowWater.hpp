#ifndef SHALLOWWATER_HPP
#define SHALLOWWATER_HPP

#include <iostream>

using namespace std;

class ShallowWater
{
public:
    double dt;
    double T;
    int Nx;
    int Ny;
    int ic;
    int in;
    double* u;
    double* v;
    double* h;
    
    ShallowWater();
    ShallowWater(const ShallowWater& rhs);
    ~ShallowWater();
    
    int SetParameters(int argc, char **argv);
    void SetInitialConditions(double H, double* x , double* y);
    void TimeIntegrate(double dx, double dy);
    void MatrixTrans(double* A,double* Atem,int Nx,int Ny);
    void MatrixMultiply(double* A, double* B, double* C,int rA, int cA, int rB, int cB);

};

#endif // SHALLOWWATER_HPP
