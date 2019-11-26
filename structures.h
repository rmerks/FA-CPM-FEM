// file: structures.h
#ifndef _STRUCTURE
#define _STRUCTURE
#include <stdlib.h>
#include "def.h"
#include <math.h>

struct ForceComponent {

  /// The model parameters
  double a, b;

  //extra variable we do not need now
  std::vector<double> x;

  // The observed points
  std::vector<double> obs;

  ForceComponent(const std::vector<double> &x):
    a(1.0),
    b(2.0),
    x(x)
  {
  }

  void operator() (std::vector<double> &res)
  {
    const size_t n=x.size();
    res.resize(n);
    res[0]=a*sin(b);
    res[1]=a*(cos(b)+sin(b))/sqrt(2);
    res[2]=a*cos(b);
    res[3]=a*(cos(b)-sin(b))/sqrt(2);
  }
};

typedef struct
{
   	//double x, y;			// minimum coordinates
	int ctag;				// id of occupying cell, 0 if no cell
} VOX;

typedef struct
{
	//double x, y;			// coordinates
	double fx, fy;			// predefined forces
	double tx, ty;			// predefined tensile forces
	double ux, uy;			// calculated displacements
	BOOL restrictx, restricty; // nodal restrictions
} NOD;


#endif
