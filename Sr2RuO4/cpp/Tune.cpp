/*
Program to tune the bands of a Hamiltonian written in the file "Hamltn.cpp" to fit ARPES data in kPath.cpp

Author: Anirudh Chandrasekaran
Date: 15 June 2022
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>
#include <stdlib.h>
#include "zheev.h"
#include <sstream>
#include <cstring>



using namespace std;

const complex<double> I(0.0,1.0);
const double Pi(acos(-1.));
const double r3(sqrt(3.));

// We begin by (re)defining some functions to accommodate Mathematica's sort of unusual CForm convention.
inline complex<double> Complex(double a, double b)
{
    return (1. * a + 1. * b * I);
}

inline double Cos(double a)
{
    return cos(a);
}

inline double Sin(double a)
{
    return sin(a);
}

inline double Power(double a, double b)
{
    return pow(a, b);
}

// The file "Hamltn.cpp" contains the definition of the Hamiltonian matrix.
#include "HamltnTune.cpp"


int main()
{
    // cout<<"\n Hello World! \n";  // <- remove!!

    double *evals;
    evals = new double[hmlt_dim];

    complex<double> **hmlt;
    hmlt = new complex<double>*[hmlt_dim];

    for(int i = 0; i < hmlt_dim; i++)
    {
        hmlt[i] = new complex<double>[hmlt_dim];
    }

    ofstream f2;

    f2.open("TuningSolution.dat");

    bool found = false;

    // now we include the loops that are defined in the Loops.cpp file here:

    #include "ParamLoops.cpp"
    
    // Time to write the data in FS to the output file

    if(!found)
    {
        cout<<"\n No solution found! Errors were larger than trial starting trial error.\n";
    }
    else
    {
        cout<<"\n Successful execution. Writing data to file now.\n";
    }

    f2<<"";

    for(int i = 0; i < param_dim; i++)
    {
        f2<<t_soln[i]<<" ";
    }

    f2<<Z_soln<<" "<<b_soln<<" "<<err_soln;

        
    for(int i = 0; i < hmlt_dim; i++)
    {
        delete[] hmlt[i];
    }

    f2.close();

    delete[] evals;     delete[] hmlt;      

    return 0;
}
