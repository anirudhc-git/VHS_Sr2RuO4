#include "kPath.cpp"

int no_k_pts = 38;

double k1, k2;
double t1, t2;

int grid_pts[2] = {130, 50};

double tmin[2] = {0., 0.01};
double tmax[2] = {12., 0.4};
double t_soln[2];

double Z_min, Z_soln, b_min, b_soln, A_sum, B_sum, Asqr_sum, Bsqr_sum, AB_sum, Z_numer, Z_denom, err, err_soln;
err_soln = 1.5e8;

for(int i1 = 0; i1 < grid_pts[0]; i1++)
{
	t1 = tmin[0] + (1. * i1 / grid_pts[0]) * (tmax[0] - tmin[0]);

	for(int i2 = 0; i2 < grid_pts[1]; i2++)
	{
		t2 = tmin[1] + (1. * i2 / grid_pts[1]) * (tmax[1] - tmin[1]);

		Z_numer = 0.;
		Z_denom = 0.;
		Z_min = 1.;
		err = 0.;
		A_sum = 0.;
		B_sum = 0.;
		AB_sum = 0.;
		Bsqr_sum = 0.;
		Asqr_sum = 0.;

		#include "PathLoop.cpp"

		if(err < err_soln)
		{
			found = true;
			err_soln = err;
			Z_soln = Z_min;
			b_soln = b_min;		
			t_soln[0] = t1;		
			t_soln[1] = t2;
		}
	}
}
