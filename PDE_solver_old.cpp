#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
void upwind_dens(double *dens, NOD* pn, NOD* pnold)
{

//use to solve drho/dt = -grad(rho*du/dt)
  double *alt_dens;
  alt_dens = new double[NV];
  double dt = par.PDEdt;
  double dx = par.VOXSIZE;

  int P,N,E,S,W;

for(int v =0;v<NV;v++){alt_dens[v]=0;}

for(int vx=0;vx<par.NVX;vx++)
{
	for(int vy=0;vy<par.NVY;vy++)
	{

	int v = vx+vy*par.NVX;
	P=v;

	N=vx+(vy-1)*par.NVX;
	if(vy==1){N=vx+vy*par.NVX;} //no flux boundary condition
	E=vx+1+vy*par.NVX;
	if(vx==par.NVX-2){E=vx+vy*par.NVX;}
	S=vx+(vy+1)*par.NVX;
	if(vy==par.NVY-2){S=vx+vy*par.NVX;}
	W=vx-1+vy*par.NVX;
	if(vx==1){W=vx+vy*par.NVX;}

	//for gradx
	double gradx;
	double uE[2], uP[2], uW[2];
	double Fex, Fwx;
	get_deform(pn,E,uE);
	get_deform(pn,P,uP);
	get_deform(pn,W,uW);
	Fex = par.COLLAGENSPEED*(uE[0]+uP[0])/(2*dt*par.COLLAGENREPEAT); 
	Fwx = par.COLLAGENSPEED*(uW[0]+uP[0])/(2*dt*par.COLLAGENREPEAT);
	if(Fex>=0){gradx = (Fex*dens[P]-Fwx*dens[W])/dx;}
	if(Fex<0){gradx = (Fex*dens[E]-Fwx*dens[P])/dx;}

	//for grad y
	double grady;
	double uN[2], uS[2];
	double Fsy, Fny;
	get_deform(pn,S,uS);
	get_deform(pn,P,uP);
	get_deform(pn,N,uN);
	Fsy = par.COLLAGENSPEED*(uS[1]+uP[1])/(2*dt*par.COLLAGENREPEAT); 
	Fny = par.COLLAGENSPEED*(uN[1]+uP[1])/(2*dt*par.COLLAGENREPEAT);
	if(Fsy>=0){grady = (Fsy*dens[P]-Fny*dens[N])/dx;} 
	if(Fsy<0){grady = (Fsy*dens[S]-Fny*dens[P])/dx;}

	alt_dens[v]=dens[v]+dt*(-gradx-grady);

	//no flux boundary condition
	if(vx==0){alt_dens[v]=dens[vx+1+vy*par.NVX];}
	if(vx==par.NVX-1){alt_dens[v]=dens[vx-1+vy*par.NVX];}
	if(vy==0){alt_dens[v]=dens[vx+(vy+1)*par.NVX];}
	if(vy==par.NVY-1){alt_dens[v]=dens[vx+(vy-1)*par.NVX];}

	}

}
for(int v =0;v<NV;v++){dens[v]=alt_dens[v];}
delete [] alt_dens;
}


