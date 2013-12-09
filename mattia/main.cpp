//
//  main.cpp
//  spline
//
//  Created by Mattia on 11/27/13.
//  Copyright (c) 2013 Mattia. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <vector>
#include <assert.h>
#include <cstdio>

double _coxDeBoorRecursion(const int j, const int d, const double * t, const double u)
{
	const double thres = 1e-6;
	double fac1 = 0.0;
	double fac2 = 0.0;
    
	if(fabs((double)d-1.0)<thres )
		return (t[j]<=u && u<t[j+1])?1:0;
    
	if(fabs(t[j+d-1]-t[j])<thres)
		fac1 = 0.0;
	else
		fac1 = (u-t[j])/(t[j+d-1]-t[j]);
    
	if(fabs(t[j+d]-t[j+1])<thres)
		fac2 = 0.0;
	else
		fac2 = (t[j+d]-u)/(t[j+d]-t[j+1]);
    
	const double bOld = _coxDeBoorRecursion(j,d-1,t,u);
	const double bOldRight = _coxDeBoorRecursion(j+1,d-1,t,u);
    
	return fac1*bOld+fac2*bOldRight;
}


double _getWidthNormalized(const double & ss, std::vector<double> WIDTH)
{
	// Very inefficient implementation, but the filling of the W vector
	// is done once at the beginning and that s it, so no big deal!
    
	// Constants
    const int N = 500;
	const int Nw = 9;
    assert(WIDTH.size()==Nw-2); // ---> The size of the vector WIDTH should be Nw-2!
	const int Nx = 2*N;
	const int n = Nw - 1;
	const int ndegree = 3;
	const int d = ndegree + 1;
	const int m = n + d + 1;
	double t[m];
    
	// Allocate control points
	double Pwx[Nw];
	double Pwy[Nw];
    
	// Allocate x and y vectors and set them to zero
	double xr[Nx];
	double yr[Nx];
	for(int i=0; i<Nx; i++){ xr[i] = 0.0; yr[i] = 0.0; }
    
	// Define control points, first the ends
	Pwx[0] = 0.0;
	Pwy[0] = 0.0;
	Pwx[Nw-1] = 1.0;
	Pwy[Nw-1] = 0.0;
    
	// Now the interior
	int counter = 0;
	const double dxw = 1.0/((double)Nw-3.0);
	for(int i=1; i<Nw-1; i++)
	{
		Pwx[i] = (i-1)*dxw;
		Pwy[i] = WIDTH[counter];
		counter++;
	}
    
	// Build t vector
	for(int i=0; i<m; i++)
	{
		if(i<d)
			t[i] = 0.0;
		else if(i>n)
			t[i] = n - d + 2.0;
		else
			t[i] = i - d + 1.0;
	}
    
	for(int i=0; i<m; i++)
		t[i] /= t[m-1];
    
	// Build the spline
	const double du = t[m-1]/((double)Nx-1.0);
	for(int i=0; i<Nx; i++)
	{
		const double u = (double)i*du;
		for(int j=0; j<n+1; j++)
		{
			xr[i] += Pwx[j]*_coxDeBoorRecursion(j,d,t,u);
			yr[i] += Pwy[j]*_coxDeBoorRecursion(j,d,t,u);
		}
	}
    
	xr[Nx-1] = Pwx[Nw-1];
	yr[Nx-1] = Pwy[Nw-1];
    
	// Generate profile via linear interpolation
	const double end = 1.0;
	const double s = ss;
	double width = 0.0;
	if( (s>=0.0) && (s<=end) )
	{
		int idx = 0;
		for(int i=0; i<Nx; i++)
		{
			if(xr[i]>s)
			{
				idx = std::max(0,i-1);
				break;
			}
		}
		const int idxPlus = std::min(idx+1,Nx-1);
		width = yr[idx] + (yr[idxPlus]-yr[idx])*(s-xr[idx])/(xr[idxPlus]-xr[idx]);
	}
    
	return width;
}

void getPoints(const double scaling, const int N, std::vector<double> WIDTH, std::vector<std::pair<double, double> > & POINTS)
{
    std::vector<double> SYMMETRICWIDTH;
    
    // First half
    for(int i=0; i < WIDTH.size(); i++)
        SYMMETRICWIDTH.push_back(WIDTH[i]);
    
    // Second half
    for(int i=0; i < WIDTH.size()-1; i++)
    {
        const int revidx = (int)WIDTH.size() - 2 - i;
        SYMMETRICWIDTH.push_back(WIDTH[revidx]);
    }
    
    
    const double ds = 1.0/(N-1);
    double s[N];
    double w[N];
    for(int i=0; i<N; i++){ s[i] = 1.0; w[i] = 0.0; }
    
    for(int i = 0; i < N-1; i++)
	{
		s[i] = ds*(double)i;
		w[i] = _getWidthNormalized(s[i],SYMMETRICWIDTH);
	}
    
    for(int i = 0; i < N; i++)
    {
        std::pair<double, double> coord(scaling*s[i],scaling*w[i]);
        POINTS.push_back(coord);
    }
    
    for(int i = N-2; i > 0; i--)
    {
        std::pair<double, double> coord(scaling*s[i],-scaling*w[i]);
        POINTS.push_back(coord);
    }
}

int main(int argc, const char * argv[])
{
    // Coordinate container
    std::vector<std::pair<double, double> > POINTS;
    
    // Width (normailzed given length 1) at the control points (half+1 due to symmetry)
    std::vector<double> WIDTH;
    WIDTH.push_back(8.9e-2);
    WIDTH.push_back(1.7e-2);
    WIDTH.push_back(1.6e-2);
    WIDTH.push_back(1.3e-2);
    
    // Number of discretization points
    const int N = 200;
    
    // Scaling factor
    const double scaling = 2.0;
    
    // Get coordinates
    getPoints(scaling, N, WIDTH, POINTS);
    
    // Print coordinates
    for(int i = 0; i < POINTS.size(); i++)
        printf("%f %f\n", POINTS[i].first, POINTS[i].second);
    
    return 0;
}

