#include "ShallowWater.hpp"
#include <iostream>
#include <boost/program_options.hpp>
#include <exception>
#include <cblas.h>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <omp.h>

using namespace std;
namespace po = boost::program_options;

//set the parameters of model from command line input
ShallowWater::ShallowWater() 
:dt(0.0),T(80.0),Nx(100),Ny(100)
{
}


//copy constructor
ShallowWater::ShallowWater(const ShallowWater& rhs) 
{
    dt = rhs.dt;
    T = rhs.T;
    Nx = rhs.Nx;
    Ny = rhs.Ny;
    ic = rhs.ic;
    in = rhs.in;
    u = rhs.u;
    v = rhs.v;
    h = rhs.h;
}

//destroy constructor
ShallowWater::~ShallowWater()
{
    //clean dynamic memory
    delete[] u;
    delete[] v;
    delete[] h;
}

//Get parameters from command line
int ShallowWater::SetParameters(int argc, char **argv)
{
    po::options_description desc("Allowed options");
        desc.add_options()
            ("help","produce help message")(
            "dt",po::value<double>()->default_value(0.1),"Time-step to use")(
            "T",po::value<double>()->default_value(80),"Total integration time")(
            "Nx",po::value<int>()->default_value(100),"Number of grid points in x")(
            "Ny",po::value<int>()->default_value(100),"Number of grid points in y")(
            "ic",po::value<int>()->default_value(1),"Index of the initial condition to use (1-4)")(
            "in",po::value<int>()->default_value(1),"Index of the method used for f evaluation (1:BLAS , 2:Loop)");
        po::variables_map vm;
        po::store(po::parse_command_line(argc,argv,desc),vm);
        po::notify(vm);
        //print help message
        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;//when changing this to function, if function return 0 then main.cpp return 0
        }
        //store value as variables
        dt = vm["dt"].as<double>();
        T = vm["T"].as<double>();
        Nx = vm["Nx"].as<int>();
        Ny = vm["Ny"].as<int>();
        ic = vm["ic"].as<int>();
        in = vm["in"].as<int>();
        if (ic<1 || ic >4) {
            throw ic;
        }
        if (in<1 || in >2) {
            throw in;
        }
        //dynamically allocate memory for matrix u,v,h
        u = new double[Nx*Ny];
        v = new double[Nx*Ny];
        h = new double[Nx*Ny];
        return 1;
}


//set initial condition of wave
void ShallowWater::SetInitialConditions(double H, double* x , double* y)
{
        for (int i = 0;i<Ny;i++){//coloum (y)
            for (int j = 0; j<Nx; j++){//row (x)
                u[i*Nx +j] = 0;
                v[i*Nx +j] = 0;
                //different test case
                if (ic == 1){//test 1
                    h[i*Nx + j] = H + exp(-pow((x[j] - 50.0),2)/25.0);
                }
                else if (ic == 2){//test 2
                    h[i*Nx + j] = H + exp(-pow((y[i] - 50.0),2)/25.0);
                }
                else if (ic == 3){//test 3
                    h[i*Nx + j] = H + exp(-(pow((x[j] - 50.0),2)+pow((y[i] - 50.0),2) )/25.0);
                }
                else{//test 4
                    h[i*Nx + j] = H + exp(-(pow((x[j] - 25.0),2) + pow((y[i] - 25.0),2))/25.0) +exp(-(pow((x[j] - 75.0),2) +pow((y[i] - 75.0),2))/25.0);
                }
            }
        }
}


void ShallowWater::MatrixTrans(double* A,double* Atem, int N1, int N2)
{
     int i;
     int j;
    for (i = 0;i<N2;i++){
        for (j = 0; j<N1; j++){
            Atem[N2*j+i] = A[N1*i+j];
        }
    }
}

void ShallowWater::MatrixMultiply(double* A, double* B, double* C, int rA, int cA, int rB, int cB)
{
     double sum1,sum2;
     int i,j,k;
    if (cA != rB){
        throw length_error("Matrix Multiplication not valid with given matrix size");
    }
    for ( i = 0; i < rA; i++){
        for (j = 0; j<cB-1 ; j+=2){
                sum1 = 0;
                sum2 = 0;
                for ( k = 0; k<cA;k++){
                    sum1+=A[i*cA + k] * B[k*cB + j];
                    sum2+=A[i*cA + k] * B[k*cB + j+1];
                }
                C[i*cB + j] = sum1;
                C[i*cB + j +1] = sum2;
        }
    }
}


void ShallowWater::TimeIntegrate(double dx, double dy)
{    
        //initialise memory
        double* Mtem = new double[Nx*Ny];
        double* M2tem = new double[Nx*Ny];
        double* utem = new double[Nx*Ny];
        double* vtem = new double[Nx*Ny];
        double* htem = new double[Nx*Ny];
        double* du_dx = new double[Nx*Ny];
        double* du_dy = new double[Nx*Ny];
        double* dv_dx = new double[Nx*Ny];
        double* dv_dy = new double[Nx*Ny];
        double* dh_dx = new double[Nx*Ny];
        double* dh_dy = new double[Nx*Ny];
        double* du_dt = new double[Nx*Ny];
        double* dv_dt = new double[Nx*Ny];
        double* dh_dt = new double[Nx*Ny];
        double* k1u = new double[Nx*Ny];
        double* k2u = new double[Nx*Ny];
        double* k3u = new double[Nx*Ny];
        double* k4u = new double[Nx*Ny];
        double* k1v = new double[Nx*Ny];
        double* k2v = new double[Nx*Ny];
        double* k3v = new double[Nx*Ny];
        double* k4v = new double[Nx*Ny];
        double* k1h = new double[Nx*Ny];
        double* k2h = new double[Nx*Ny];
        double* k3h = new double[Nx*Ny];
        double* k4h = new double[Nx*Ny];
        double* Ax = new double[Nx*Nx];
        double* Ay = new double[Ny*Ny];
        double g = 9.81;
        int i = 0;
        int j = 0;
        int t = 0;
        int k = 0;
        int threadid,nthread;
        //matrix Ax contains index used for finite difference for x
        //initialise parallel region
        #pragma omp parallel 
        {
            //split initialising of Ax and Ay matrices to different sections
            #pragma omp sections
            {
                #pragma omp section
                {
                    //matrix Ax contains index used for finite difference for x
                    #pragma omp parallel for 
                    for (int i = 0; i<Nx;i++){
                        Ax[i*Nx+((i-3+Nx)%Nx)] = -1.0/(60.0*dx);
                        Ax[i*Nx+((i-2+Nx)%Nx)] = 3.0/(20.0*dx);
                        Ax[i*Nx+((i-1+Nx)%Nx)] = -3.0/(4.0*dx);
                        Ax[i*Nx+((i+1)%Nx)] = 3.0/(4.0*dx);
                        Ax[i*Nx+((i+2)%Nx)] = -3.0/(20.0*dx);
                        Ax[i*Nx+((i+3)%Nx)] = 1.0/(60.0*dx);
                    }
                }
                
                #pragma omp section
                {
                    //matrix Ay contains index used for finite difference for y
                    #pragma omp parallel for 
                    for (int i = 0; i<Ny;i++){
                        Ay[i*Ny+((i-3+Ny)%Ny)] = -1.0/(60.0*dy);
                        Ay[i*Ny+((i-2+Ny)%Ny)] = 3.0/(20.0*dy);
                        Ay[i*Ny+((i-1+Ny)%Ny)] = -3.0/(4.0*dy);
                        Ay[i*Ny+((i+1)%Ny)] = 3.0/(4.0*dy);
                        Ay[i*Ny+((i+2)%Ny)] = -3.0/(20.0*dy);
                        Ay[i*Ny+((i+3)%Ny)] = 1.0/(60.0*dy);
                    }
                }

            }
        }
        
        
        if (in == 1){
            //use blas
            for (t = 0;t<T/dt;t++){
                //copy matrix u,v,h to utem, utem = u
                cblas_dcopy(Nx*Ny, u,1,utem,1);
                cblas_dcopy(Nx*Ny, v,1,vtem,1);
                cblas_dcopy(Nx*Ny, h,1,htem,1);
                for (k = 0; k<4;k++){//4th runge kutta index
                    switch(k){
                        case 0:
                            break;
                        case 1:
                            //u haven't be changed at this case
                            //yn+dt*k1/2
                            cblas_daxpy(Nx*Ny,dt/2.0,k1u,1,u,1);
                            cblas_daxpy(Nx*Ny,dt/2.0,k1v,1,v,1);                          
                            cblas_daxpy(Nx*Ny,dt/2.0,k1h,1,h,1);
                            break;
                        case 2:
                            //reset u = utem
                            cblas_dcopy(Nx*Ny, utem,1,u,1);
                            cblas_dcopy(Nx*Ny, vtem,1,v,1);
                            cblas_dcopy(Nx*Ny, htem,1,h,1);
                            //yn+dt*k2/2
                            cblas_daxpy(Nx*Ny,dt/2.0,k2u,1,u,1);
                            cblas_daxpy(Nx*Ny,dt/2.0,k2v,1,v,1);                          
                            cblas_daxpy(Nx*Ny,dt/2.0,k2h,1,h,1);
                            break;
                        case 3:
                            //reset u = utem
                            cblas_dcopy(Nx*Ny, utem,1,u,1);
                            cblas_dcopy(Nx*Ny, vtem,1,v,1);
                            cblas_dcopy(Nx*Ny, htem,1,h,1);
                            //yn+dt*k3/2
                            cblas_daxpy(Nx*Ny,dt,k3u,1,u,1);
                            cblas_daxpy(Nx*Ny,dt,k3v,1,v,1);                          
                            cblas_daxpy(Nx*Ny,dt,k3h,1,h,1);
                            break;
                    }
                    
                    //u derivatives
                    //calculate du_dx derivative
                    //apply matrix-vector multiplication for each row, therefore u,v,h need to be transposed
                    MatrixTrans(u,Mtem,Nx,Ny);
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,Nx,Ny,Nx,1.0,Ax,Nx,Mtem,Ny,0.0,M2tem,Ny);
                    MatrixTrans(M2tem,du_dx,Ny,Nx);//transpose M2tem to get du_dt
                    //calculate du_dy derivative
                    //no need for transpose
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,Ny,Nx,Ny,1.0,Ay,Ny,u,Nx,0.0,du_dy,Nx);
                    
                    //v derivatives
                    MatrixTrans(v,Mtem,Nx,Ny);//transpose v
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,Nx,Ny,Nx,1.0,Ax,Nx,Mtem,Ny,0.0,M2tem,Ny);
                    MatrixTrans(M2tem,dv_dx,Ny,Nx);
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,Ny,Nx,Ny,1.0,Ay,Ny,v,Nx,0.0,dv_dy,Nx);
                    
                    //h derivatives
                    MatrixTrans(h,Mtem,Nx,Ny);//transpose h
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,Nx,Ny,Nx,1.0,Ax,Nx,Mtem,Ny,0.0,M2tem,Ny);
                    MatrixTrans(M2tem,dh_dx,Ny,Nx); 
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,Ny,Nx,Ny,1.0,Ay,Ny,h,Nx,0.0,dh_dy,Nx);
                    
                    for (i = 0; i<Nx*Ny; i++){
                            du_dt[i] = -(u[i]*du_dx[i] +v[i]*du_dy[i] + g*dh_dx[i]);
                            dv_dt[i] = -(u[i]*dv_dx[i] +v[i]*dv_dy[i] + g*dh_dy[i]);
                            dh_dt[i] = -(h[i]*du_dx[i] +u[i]*dh_dx[i] + h[i]*dv_dy[i] +v[i]*dh_dy[i]);
                        }
                    switch(k){
                        case 0:
                            cblas_dcopy(Nx*Ny,du_dt,1,k1u,1);
                            cblas_dcopy(Nx*Ny,dv_dt,1,k1v,1);
                            cblas_dcopy(Nx*Ny,dh_dt,1,k1h,1);
                            break;
                        case 1:
                            cblas_dcopy(Nx*Ny,du_dt,1,k2u,1);
                            cblas_dcopy(Nx*Ny,dv_dt,1,k2v,1);
                            cblas_dcopy(Nx*Ny,dh_dt,1,k2h,1);
                            break;
                        case 2:
                            cblas_dcopy(Nx*Ny,du_dt,1,k3u,1);
                            cblas_dcopy(Nx*Ny,dv_dt,1,k3v,1);
                            cblas_dcopy(Nx*Ny,dh_dt,1,k3h,1);
                            break;
                        case 3:
                            cblas_dcopy(Nx*Ny,du_dt,1,k4u,1);
                            cblas_dcopy(Nx*Ny,dv_dt,1,k4v,1);
                            cblas_dcopy(Nx*Ny,dh_dt,1,k4h,1);
                            break;
                    }
                }
                for (int i= 0;i<Nx*Ny;i++){
                        u[i] = utem[i] + (1.0/6.0) * (k1u[i] +2.0*k2u[i] +2.0*k3u[i]+k4u[i])*dt;
                        v[i] = vtem[i] + (1.0/6.0) * (k1v[i] +2.0*k2v[i] +2.0*k3v[i]+k4v[i])*dt;
                        h[i] = htem[i] + (1.0/6.0) * (k1h[i] +2.0*k2h[i] +2.0*k3h[i]+k4h[i])*dt;
                        }  
            }
        }
        else{
            //use loop
            for (int t = 0;t<T/dt;t++){
                //store value of u,v,h in a temperory storage variable, so they will not be adjusted when evaluating k
                    for ( i = 0; i<Nx*Ny; i++){
                            utem[i] = u[i];
                            vtem[i] = v[i];
                            htem[i] = h[i];
                    }
                    
                for (int k = 0; k<4;k++){//4th runge kutta index
                    switch(k){
                        case 0://evaluate k1
                            break;//dont change u,v,h
                        case 1://use k1 to evaluate k2, k2 = f(yn + dt*k1/2)
                        for (i = 0; i<Nx*Ny ; i++){
                            u[i] = utem[i] + dt*k1u[i]/2.0;
                            v[i] = vtem[i] + dt*k1v[i]/2.0;
                            h[i] = htem[i] + dt*k1h[i]/2.0;
                        }
                        
                            break;
                        case 2://k3 = f(yn + dt*k2/2)
                        for (i = 0; i<Nx*Ny ; i++){
                            u[i] = utem[i] + dt*k2u[i]/2.0;
                            v[i] = vtem[i] + dt*k2v[i]/2.0;
                            h[i] = htem[i] + dt*k2h[i]/2.0;
                        }
                        
                            break;
                        case 3://k4 = f(yn+dt*k3)
                        for (i = 0; i<Nx*Ny ; i++){
                            u[i] = utem[i] + dt*k3u[i];
                            v[i] = vtem[i] + dt*k3v[i];
                            h[i] = htem[i] + dt*k3h[i];
                        }
                        
                            break;
                    }
                            
                    //calculate x derivative
                    //apply matrix-vector multiplication for each row
                    //therefore u,v,h need to be transposed
                    MatrixTrans(u,Mtem,Nx,Ny);
                    MatrixMultiply(Ax,Mtem,M2tem,Nx,Nx,Nx,Ny);
                    MatrixTrans(M2tem,du_dx,Ny,Nx);//transpose M2tem to get du_dx
                    //transpose v
                    MatrixTrans(v,Mtem,Nx,Ny);
                    MatrixMultiply(Ax,Mtem,M2tem,Nx,Nx,Nx,Ny);
                    MatrixTrans(M2tem,dv_dx,Ny,Nx);
                    //transpose h
                    MatrixTrans(h,Mtem,Nx,Ny);
                    MatrixMultiply(Ax,Mtem,M2tem,Nx,Nx,Nx,Ny);
                    MatrixTrans(M2tem,dh_dx,Ny,Nx); 
                    
                    //calculate y derivative
                    //apply matrix-vector multiplication for each column
                    //no need to transpose 
                    MatrixMultiply(Ay,u,du_dy,Ny,Ny,Ny,Nx);
                    MatrixMultiply(Ay,v,dv_dy,Ny,Ny,Ny,Nx);
                    MatrixMultiply(Ay,h,dh_dy,Ny,Ny,Ny,Nx);
                    
                    //calculate time derivative(f)
                    for (int i = 0; i < Nx*Ny; i++){
                        du_dt[i] = -(u[i]*du_dx[i] +v[i]*du_dy[i] + g*dh_dx[i]);
                        dv_dt[i] = -(u[i]*dv_dx[i] +v[i]*dv_dy[i] + g*dh_dy[i]);
                        dh_dt[i] = -(h[i]*du_dx[i] +u[i]*dh_dx[i] + h[i]*dv_dy[i] +v[i]*dh_dy[i]);
                        //store value into different k matrix
                        switch(k){
                            case 0:
                                k1u[i] = du_dt[i];
                                k1v[i] = dv_dt[i];
                                k1h[i] = dh_dt[i];
                                break;
                            case 1:
                                k2u[i] = du_dt[i];
                                k2v[i] = dv_dt[i];
                                k2h[i] = dh_dt[i];
                                break;
                            case 2:
                                k3u[i] = du_dt[i];
                                k3v[i] = dv_dt[i];
                                k3h[i] = dh_dt[i];
                                break;
                            case 3:
                                k4u[i] = du_dt[i];
                                k4v[i] = dv_dt[i];
                                k4h[i] = dh_dt[i];
                                break;
                        }
                    }
                }
             //assign value for u, v, h at the new time step (yn+1) using RK4
             //using OpenMP for here increase runing time
                for (int i= 0;i<Nx*Ny;i++){
                        u[i] = utem[i] + (1.0/6.0) * (k1u[i] +2.0*k2u[i] +2.0*k3u[i]+k4u[i])*dt;
                        v[i] = vtem[i] + (1.0/6.0) * (k1v[i] +2.0*k2v[i] +2.0*k3v[i]+k4v[i])*dt;
                        h[i] = htem[i] + (1.0/6.0) * (k1h[i] +2.0*k2h[i] +2.0*k3h[i]+k4h[i])*dt;
                }
          }
        }
        
        
        //free memory
        delete[] Mtem;
        delete[] M2tem;
        delete[] Ax;
        delete[] Ay;
        delete[] utem;
        delete[] vtem;
        delete[] htem;
        delete[] du_dx;
        delete[] du_dy;
        delete[] dv_dx;
        delete[] dv_dy;
        delete[] dh_dx;
        delete[] dh_dy;
        delete[] du_dt;
        delete[] dv_dt;
        delete[] dh_dt;
        delete[] k1u;
        delete[] k2u;
        delete[] k3u;
        delete[] k4u;
        delete[] k1v;
        delete[] k2v;
        delete[] k3v;
        delete[] k4v;
        delete[] k1h;
        delete[] k2h;
        delete[] k3h;
        delete[] k4h;
}
