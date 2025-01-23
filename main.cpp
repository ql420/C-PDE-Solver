#include <iostream>
#include <exception>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include "ShallowWater.hpp"

using namespace std;


int main(int argc, char **argv)
{
    try {
        
        //initialise class
        ShallowWater wave;
        //set parameters of model
        int help;
        help = wave.SetParameters(argc, argv);
        if (help==0){//return 0 if command line type "--help"
            return 0;
        }
        
        //initialise x,y coordinates
        double* x = new double[wave.Nx];
        double* y = new double[wave.Ny];
        double dx = 1.0;
        double dy = 1.0;
        double H = 10.0;//mean surface height
        x[0] = 0;
        y[0] = 0;
        //populte x and y, point x = 100 is not considered because of periodic boundary condition
        for (int i=1; i<wave.Nx;i++){
            x[i] = x[i-1] + dx;
        }
        for (int i=1; i<wave.Ny;i++){
            y[i] = y[i-1] + dy;
        }
        
        //set initial condition of wave
        wave.SetInitialConditions(H,x,y);
        
        //Time Integrate
        wave.TimeIntegrate(dx,dy);

        //output data to output.txt
        ofstream vOut("output.txt", ios::out | ios::trunc);
        vOut.precision(5);
        for (int i = 0;i<wave.Ny;i++){
            for (int j = 0; j<wave.Nx;j++){
                vOut << setw(20) << x[j] << setw(20) << y[i] << setw(20) << wave.u[wave.Nx*i+j] << setw(20) << wave.v[wave.Nx*i+j] << setw(20) << wave.h[wave.Nx*i+j] << endl;
            }
            vOut<<endl;
        }
        vOut.close();
        
        //free memory
        delete[] x;
        delete[] y;
    }
    
    catch (int ic) {
        cout << "Index error, index ic should be within 1-4, index in should be 1 or 2" << "\n";
        return 1;
    }
    catch (const length_error& e) {
        cout << "Error!" << e.what() << endl;
        return 2;
    }
    
	return 3;
}
