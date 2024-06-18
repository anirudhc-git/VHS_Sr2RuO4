/*
Program to compute the Fermi surface of a Hamiltonian written in the file "Hamltn.cpp"
It is written in the file "FS.dat", as a list with the following format:
{{%band 1% {{k1,k2,...,kn}, E(k1,k2,...,kn)}, {{k1,k2,...,kn}, E(k1,k2,...,kn)},...}, {%band 2% {{k1,k2,...,kn}, E(k1,k2,...,kn)}, {{k1,k2,...,kn}, E(k1,k2,...,kn)},...},.., {%band N% {{k1,k2,...,kn}, E(k1,k2,...,kn)}, {{k1,k2,...,kn}, E(k1,k2,...,kn)},...}}

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
#include "Hamltn.cpp"


// We need our custom made functions to check if the numbers read off params.dat are legitimitate parameters for the calculation to follow.
bool is_number(const char *str)
{
    int dots = 0;
    
    if(strlen(str) == 0)
    {
        return false;
    }

    if(str[0] == '-')
    {
        for(int i = 1; i < strlen(str); i++)
        {
            if(!isdigit(str[i]))
            {
                if(str[i] == '.')
                {
                    if(dots == 0)
                    {
                        dots++;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    return false;
                }

            }
        }
    }
    else
    {
        for(int i = 0; i < strlen(str); i++)
        {
            if(!isdigit(str[i]))
            {
                if(str[i] == '.')
                {
                    if(dots == 0)
                    {
                        dots++;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }
        }
    }

    return true;
}

bool is_positive_integer(const char *str)
{
    if(strlen(str) == 0)
    {
        return false;
    }

    else
    {
        for(int i = 0; i < strlen(str); i++)
        {
            if(!isdigit(str[i]))
            {
                return false;
            }
        }

        return true;
    }
}


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

    ifstream f1("params.dat");
    ofstream f2;

    double kmin[vec_dim], kmax[vec_dim];
    int grid_pts[vec_dim], no_FS_pts[hmlt_dim];

    double EF, resln;
    string outfile, str_line, wrd;

    /*
    for(int i = 0; i < vec_dim; i++)
    {
        f1>>kmin[i]>>kmax[i];
    }
    */

    // The first few lines of params.dat contains the values of kmin and kmax.
    for(int i = 0; i < vec_dim; i++)
    {
        getline(f1, str_line);

        istringstream strm(str_line);

        if(strm>>wrd)
        {
            if(is_number(wrd.c_str()))
            {
                kmin[i] = atof(wrd.c_str());
            }
            else
            {
                return 1;
            }
        }
        else
        {
            return 1;
        }

        if(strm>>wrd)
        {
            if(is_number(wrd.c_str()))
            {
                kmax[i] = atof(wrd.c_str());
            }
            else
            {
                return 1;
            }
        }
        else
        {
            return 1;
        }

    }

    // Next we load the grid points
    for(int i = 0; i < vec_dim; i++)
    {
        getline(f1, str_line);

        if(is_positive_integer(str_line.c_str()))
        {
            grid_pts[i] = atoi(str_line.c_str());
        }
        else
        {
            return 1;
        }

    }

    getline(f1, str_line);
    if(is_number(str_line.c_str()))
    {
        EF = atof(str_line.c_str());
    }
    else
    {
        return 1;
    }

    getline(f1, str_line);
    if(is_number(str_line.c_str()))
    {
        resln = abs(atof(str_line.c_str()));
    }
    else
    {
        return 1;
    }

    getline(f1, outfile);


    // cout<<"\n space dimension = "<<vec_dim<<"\n"; 
    // cout<<"\n Fermi level = "<<EF;
    // cout<<"\n Resolution = "<<resln;
    // cout<<"\n kmin = "<<kmin[0]<<", kmax = "<<kmax[0];
    // cout<<"\n kmin = "<<kmin[1]<<", kmax = "<<kmax[1];
    // cout<<"\n grid_x = "<<grid_pts[0]<<", grid_y = "<<grid_pts[1];
    // cout<<"\n EF = "<<EF;
    // cout<<"\n resln = "<<resln;
    // cout<<"\n address = "<<outfile<<"\n";

    int no_pts[hmlt_dim], tot_pts;

    tot_pts = 1;

    for(int i = 0; i < vec_dim; i++)
    {
        tot_pts = tot_pts * grid_pts[i];
    }


    // Now we define the variable that will be used to store the Fermi surface data and allocate memory for it:
    double*** FS;

    FS = new double**[hmlt_dim];

    for(int i = 0; i < hmlt_dim; i++)
    {
        FS[i] = new double*[tot_pts];
    }

    for(int i = 0; i < hmlt_dim; i++)
    {
        for(int j = 0; j < tot_pts; j++)
        {
            FS[i][j] = new double[vec_dim + 1];
        }
    }

    // The array no_FS_pts will contain the total no of points in each branch of the Fermi surface:

    for(int i = 0; i < hmlt_dim; i++)
    {
        no_FS_pts[i] = -1;
    }

    // now we include the loops that are defined in the Loops.cpp file here:

    #include "Loops.cpp"
    
    // Time to write the data in FS to the output file

    
    f2.open(outfile.c_str());

    f2<<"{";

    for(int i = 0; i < hmlt_dim; i++)
    {
        f2<<"{";

        if(no_FS_pts[i] >= 0)
        {
            // cout<<"\n Number of points in band "<<i+1<<" = "<<no_FS_pts[i+1];
            
            for(int j = 0; j <= no_FS_pts[i]; j++)
            {
                f2<<"{{"<<FS[i][j][0];

                for(int k = 1; k < vec_dim; k++)
                {
                    f2<<","<<FS[i][j][k];
                }

                f2<<"},"<<FS[i][j][vec_dim]<<"}";

                if(j != no_FS_pts[i])
                {
                    f2<<",";
                }
            }

        }

        f2<<"}";

        if(i != hmlt_dim - 1)
        {
            f2<<",";
        }
    }

    f2<<"}";

    

    for(int i = 0; i < hmlt_dim; i++)
    {
        for(int j = 0; j < tot_pts; j++)
        {
            delete[] FS[i][j];
        }
    }
    
    for(int i = 0; i < hmlt_dim; i++)
    {
        delete[] hmlt[i];
        delete[] FS[i];
    }

    f2.close();
    f1.close();

    delete[] evals;     delete[] hmlt;      delete[] FS;

    return 0;
}
