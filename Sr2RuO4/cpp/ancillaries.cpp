const complex<double> I(0.0,1.0);

// We begin by (re)defining some functions to accommodate Mathematica's sort of unusual CForm convention.
inline complex<double> Complex(double a, double b)
{
	return (1. * a + 1. * b * I);
}

inline double Cos(double a)
{
	return cos(a);
}

inline double Power(double a, double b)
{
	return pow(a, b);
}

