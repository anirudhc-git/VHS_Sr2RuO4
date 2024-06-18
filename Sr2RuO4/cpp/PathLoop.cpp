int no_E_pts, no_bands;
no_E_pts = 0;

for(int l = 0; l < no_k_pts; l++)
{
	k1 = kPath[l][0];
	k2 = kPath[l][1];

	hmlt_def(hmlt, k1, k2, t1, t2);
	zheev(hmlt, hmlt_dim, evals);


	no_bands = (int) kPath[l][2];

	no_E_pts += no_bands;

	for(int j = 0; j < no_bands; j++)
	{

		AB_sum += kPath[l][3 + j] * evals[(int)(kPath[l][3 + no_bands + j])];

		A_sum += kPath[l][3 + j];

		B_sum += evals[(int)(kPath[l][3 + no_bands + j])];

		Asqr_sum += kPath[l][3 + j] * kPath[l][3 + j];

		Bsqr_sum += evals[(int)(kPath[l][3 + no_bands + j])] * evals[(int)(kPath[l][3 + no_bands + j])];

	}
}

Z_numer = AB_sum - A_sum * B_sum / (1. * no_E_pts);
Z_denom = Bsqr_sum - B_sum * B_sum / (1. * no_E_pts);

if(Z_denom == 0 || (Z_numer * Z_denom) < 0)
{
	cout<<"\n Z factor denominator vanished at iteration TODO!\n";
}
else
{
	Z_min = Z_numer / Z_denom;
	b_min = (Z_min * B_sum - A_sum)/(1. * no_k_pts);
}

err = Asqr_sum + Z_min * Z_min * Bsqr_sum + 2. * b_min * (A_sum - Z_min * B_sum) - 2. * Z_min * AB_sum + b_min * b_min * no_E_pts;