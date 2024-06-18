double k1, k2;

for(int i1 = 0; i1 < grid_pts[0]; i1++)
{
	k1 = kmin[0] + (1. * i1 / grid_pts[0]) * (kmax[0] - kmin[0]);

	for(int i2 = 0; i2 < grid_pts[1]; i2++)
	{
		k2 = kmin[1] + (1. * i2 / grid_pts[1]) * (kmax[1] - kmin[1]);

		hmlt_def(hmlt, k1, k2);
		zheev(hmlt, hmlt_dim, evals);

		for(int j = 0; j < hmlt_dim; j++)
		{
			if(abs(evals[j] - EF) < resln)
			{
				no_FS_pts[j]++;

				FS[j][no_FS_pts[j]][0] = k1;
				FS[j][no_FS_pts[j]][1] = k2;
				FS[j][no_FS_pts[j]][2] = evals[j];
			}
		}
	}
}
