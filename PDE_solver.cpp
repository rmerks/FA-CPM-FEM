#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
void upwind_dens(double *dens, NOD* pn, NOD* pnold,VOX* pv)
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
	double uEold[2], uPold[2], uWold[2];
	double Fex, Fwx;
	get_deform(pn,E,uE);
	get_deform(pn,P,uP);
	get_deform(pn,W,uW);
	get_deform(pnold,E,uEold);
	get_deform(pnold,P,uPold);
	get_deform(pnold,W,uWold);
	Fex = par.COLLAGENSPEED*(uE[0]-uEold[0]+uP[0]-uPold[0])/(2*dt); 
	Fwx = par.COLLAGENSPEED*(uW[0]-uWold[0]+uP[0]-uPold[0])/(2*dt);
	if(Fex>=0){gradx = (Fex*dens[P]-Fwx*dens[W])/dx;}
	if(Fex<0){gradx = (Fex*dens[E]-Fwx*dens[P])/dx;}

	//for grad y
	double grady;
	double uN[2], uS[2];
	double uNold[2], uSold[2];
	double Fsy, Fny;
	get_deform(pn,S,uS);
	get_deform(pn,P,uP);
	get_deform(pn,N,uN);
	get_deform(pnold,S,uSold);
	get_deform(pnold,P,uPold);
	get_deform(pnold,N,uNold);
//if(S==100){cout << "uS[1]" << uS[1] << endl; cout << "uSold[1]" << uSold[1] << endl;}
	Fsy = par.COLLAGENSPEED*(uS[1]-uSold[1]+uP[1]-uPold[1])/(2*dt); 
	Fny = par.COLLAGENSPEED*(uN[1]-uNold[1]+uP[1]-uPold[1])/(2*dt);
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
//for(int v =0;v<NV;v++){double diff = dens[v]-alt_dens[v];if(pv[v].ctag){Acell[pv[v].ctag-1]+=-diff;}}
//for(int v =0;v<NV;v++){dens[v]=alt_dens[v];}
//we should change Acell2 when using the flux
delete [] alt_dens;
}

/*
void forwardeuler_FA(double* FA, NOD* pn,NOD* pnold, VOX* pv, double* dens,int NRc,double* tensionold,double* sumFA,double** klocal)
{

	double etractionstress[2];
        double dt = par.PDEdt;

	for(int tt=0;tt<par.PDEREPEAT;tt++)
	{
	double *Nt = new double[NRc]; for(int c=0;c<NRc;c++){Nt[c]=par.CAPACITYFA-sumFA[c];}//20
	for(int n=1;n<NV;n++){if(pv[n].ctag>0){

		int cellnr=pv[n].ctag;
		double tension=0;
		double tension0=0; //later veranderen
		get_etractionstress(pn,n, etractionstress);
	        double strx=etractionstress[0]; double stry=etractionstress[1];
		tension0=sqrt(strx*strx + stry*stry)*par.VOXSIZE*par.VOXSIZE;
		double Km = par.YOUNGS*par.VOXSIZE*1/(1-par.POISSON*par.POISSON);
		if(par.DUROTAXIS){
			double v= par.YOUNGS*par.VOXSIZE*1/(1-par.POISSON*par.POISSON);
			double Km = v*(1+par.GRADIENT*double(n%par.NVX)/par.YOUNGS);
			Km=Km*par.THICKNESS;

		}
		double tk=tension0/(par.VISC*Km);
		tension = tensionold[n]+(tension0-tensionold[n])*(1-exp(-(tt*dt)/tk));
		if(tk==0){tension=tension0;}
		double g =0;
		if(!par.PATTERN){g=par.GROWTHFA;}
		if(par.PATTERN)
		{
			int nx=n%par.NVX;
			int ny=n/par.NVX;
			if((nx%par.PATTERNC==1) && (ny%par.PATTERNC==1)){g = par.GROWTHFA;}
		} 
		double sliptension = par.SLIPTENSION;  
		double catchtension = par.CATCHTENSION;
		double tension2 = tension*1e12; //because pN
		double rate = exp((tension2/FA[n]-sliptension))+exp(-(tension2/FA[n]-catchtension));
		double Na = Nt[cellnr-1]; if(Na<0){Na=0;}
		//if(Na+FA[n]>10000){Na=10000-FA[n];}

		double maxfa = 50*par.VOXSIZE*par.VOXSIZE/8e-15;
		double s=FA[n]-dt*rate*FA[n]+dt*g*Na*(1-FA[n]/maxfa); //logistic growth 
		if(FA[n]>=maxfa){s=FA[n]-dt*rate*FA[n];} //try



		if(s<0){s=0;}
		//try
		if(FA[n]==0){s=0;} //want rate is infinite
		double diff = FA[n]-s;
		FA[n]=s;


		sumFA[cellnr-1]=sumFA[cellnr-1]-diff;


		}}


	}
}
*/
	


void calctension(NOD* pn,double** klocal)
{

        for(int n=0;n<NN;n++)
        {



		double fx = pn[n].fx;
		double fy = pn[n].fy;
		double ux = pn[n].ux;
		double uy = pn[n].uy;

		double forceleft = 0;
		double forceright = 0;
		double forceup = 0;
		double forcedown = 0;

		//calculate tension (assume we have only one cell)
		
		//there are 4 forces from surrounding elements


		//check: since this is not ok for boundary nodes


		double kx,ky,fx1,fx2,fx3;
		kx=klocal[0][6]; ky=klocal[0][7];
		fx1=kx*(pn[n+NNX].ux-ux)+ky*(pn[n+NNX].uy-uy); 
		kx=klocal[0][4]; ky=klocal[0][5];
		fx2=kx*(pn[n+NNX+1].ux-ux)+ky*(pn[n+NNX+1].uy-uy); 
		kx=klocal[0][2]; ky=klocal[0][3];
		fx3=kx*(pn[n+1].ux-ux)+ky*(pn[n+1].uy-uy); 
		//for element 3 (top right)
		double fxe3 = fx1+fx2+fx3; 

		kx=klocal[6][4]; ky=klocal[6][5];
		fx1=kx*(pn[n+1].ux-ux)+ky*(pn[n+1].uy-uy); 
		kx=klocal[6][2]; ky=klocal[6][3];
		fx2=kx*(pn[n-NNX+1].ux-ux)+ky*(pn[n-NNX+1].uy-uy); 
		kx=klocal[6][0]; ky=klocal[6][1];
		fx3=kx*(pn[n-NNX].ux-ux)+ky*(pn[n-NNX].uy-uy); 
		//for element 1 (bottom right)
		double fxe1 = fx1+fx2+fx3; 

		kx=klocal[4][2]; ky=klocal[4][3];
		fx1=kx*(pn[n-NNX].ux-ux)+ky*(pn[n-NNX].uy-uy); 
		kx=klocal[4][0]; ky=klocal[4][1];
		fx2=kx*(pn[n-NNX-1].ux-ux)+ky*(pn[n-NNX-1].uy-uy);  
		kx=klocal[4][6]; ky=klocal[4][7];
		fx3=kx*(pn[n-1].ux-ux)+ky*(pn[n-1].uy-uy);  
		//for element 0 (bottom left)
		double fxe0 = fx1+fx2+fx3; 




		
		kx=klocal[2][0]; ky=klocal[2][1];
		fx1=kx*(pn[n-1].ux-ux)+ky*(pn[n-1].uy-uy);  
		kx=klocal[2][6]; ky=klocal[2][7];
		fx2=kx*(pn[n+NNX-1].ux-ux)+ky*(pn[n+NNX-1].uy-uy); 
		kx=klocal[2][4]; ky=klocal[2][5];
		fx3=kx*(pn[n+NNX].ux-ux)+ky*(pn[n+NNX].uy-uy); 
		//for element 2 (top left)
		double fxe2 = fx1+fx2+fx3;

		if(fx<0){forceleft += -fx; }
		if(fx>0){forceright += -fx;}

		if(fxe0<0){forceright += fxe0;}
		if(fxe0>0){forceleft += fxe0;}
		if(fxe1<0){forceright += fxe1;}
		if(fxe1>0){forceleft += fxe1;}
		if(fxe2<0){forceright += fxe2;}
		if(fxe2>0){forceleft += fxe2;}
		if(fxe3<0){forceright += fxe3;}
		if(fxe3>0){forceleft += fxe3;}


		//check: since this is not ok for boundary nodes
		double fy1,fy2,fy3;
		kx=klocal[1][6]; ky=klocal[1][7];
		fy1=kx*(pn[n+NNX].ux-ux)+ky*(pn[n+NNX].uy-uy); 
		kx=klocal[1][4]; ky=klocal[1][5];
		fy2=kx*(pn[n+NNX+1].ux-ux)+ky*(pn[n+NNX+1].uy-uy); 
		kx=klocal[1][2]; ky=klocal[1][3];
		fy3=kx*(pn[n+1].ux-ux)+ky*(pn[n+1].uy-uy); 
		//for element 3 (top right)
		double fye3 = fy1+fy2+fy3; 

		kx=klocal[7][4]; ky=klocal[7][5];
		fy1=kx*(pn[n+1].ux-ux)+ky*(pn[n+1].uy-uy); 
		kx=klocal[7][2]; ky=klocal[7][3];
		fy2=kx*(pn[n-NNX+1].ux-ux)+ky*(pn[n-NNX+1].uy-uy); 
		kx=klocal[7][0]; ky=klocal[7][1];
		fy3=kx*(pn[n-NNX].ux-ux)+ky*(pn[n-NNX].uy-uy); 
		//for element 1 (bottom right)
		double fye1 = fy1+fy2+fy3; 

		kx=klocal[5][2]; ky=klocal[5][3];
		fy1=kx*(pn[n-NNX].ux-ux)+ky*(pn[n-NNX].uy-uy); 
		kx=klocal[5][0]; ky=klocal[5][1];
		fy2=kx*(pn[n-NNX-1].ux-ux)+ky*(pn[n-NNX-1].uy-uy);  
		kx=klocal[5][6]; ky=klocal[5][7];
		fy3=kx*(pn[n-1].ux-ux)+ky*(pn[n-1].uy-uy);  
		//for element 0 (bottom left)
		double fye0 = fy1+fy2+fy3; 
		
		kx=klocal[3][0]; ky=klocal[3][1];
		fy1=kx*(pn[n-1].ux-ux)+ky*(pn[n-1].uy-uy);  
		kx=klocal[3][6]; ky=klocal[3][7];
		fy2=kx*(pn[n+NNX-1].ux-ux)+ky*(pn[n+NNX-1].uy-uy); 
		kx=klocal[3][4]; ky=klocal[3][5];
		fy3=kx*(pn[n+NNX].ux-ux)+ky*(pn[n+NNX].uy-uy); 
		//for element 2 (top left)
		double fye2 = fy1+fy2+fy3;

		if(fy<0){forcedown += -fy; }
		if(fy>0){forceup += -fy;}

		if(fye0<0){forceup += fye0;}
		if(fye0>0){forcedown += fye0;}
		if(fye1<0){forceup += fye1;}
		if(fye1>0){forcedown += fye1;}
		if(fye2<0){forceup += fye2;}
		if(fye2>0){forcedown += fye2;}
		if(fye3<0){forceup += fye3;}
		if(fye3>0){forcedown += fye3;}

if(fx==0){
if(ux<0 & forceleft<forceright){pn[n].tx=forceleft;}
if(ux<0 & forceright <= forceleft){pn[n].tx=forceright;}
if(ux>=0 & forceright < forceleft){pn[n].tx=forceleft;}
if(ux>=0 & forceleft<=forceright){pn[n].tx=forceright;}
}

if(fy==0){
if(uy<0 & forceup>forcedown){pn[n].ty=forcedown;}
if(uy<0 & forcedown>=forceup){pn[n].ty=forceup;}
if(uy>=0 & forceup>forcedown){pn[n].ty=forceup;}
if(uy>=0 & forcedown>=forceup){pn[n].ty=forcedown;}
}

if(!fx==0)
{
if(fx<0 & forceleft<forceright){pn[n].tx=forceleft;}
if(fx<0 & forceright <= forceleft){pn[n].tx=forceright;}
if(fx>=0 & forceright < forceleft){pn[n].tx=forceleft;}
if(fx>=0 & forceleft<=forceright){pn[n].tx=forceright;}
}

if(!fy==0){
if(fy<0 & forceup>forcedown){pn[n].ty=forcedown;}
if(fy<0 & forcedown>=forceup){pn[n].ty=forceup;}
if(fy>=0 & forceup>forcedown){pn[n].ty=forceup;}
if(fy>=0 & forcedown>=forceup){pn[n].ty=forcedown;}
}



if(abs(fye0+fye1-fy-forceup)<1e-15 & !fy==0){pn[n].ty=0; }
if(abs(fxe0+fxe2-fx-forceleft)<1e-15 & !fx==0){pn[n].tx=0; }
if(fye0+fye1<0 & fy==0){pn[n].ty=0;}
if(fxe1+fxe3<0 & fx==0){pn[n].tx=0;}


/* if(n==NN/2+2+4*NNX-2-1-NNX-2*NNX+3+2-3*NNX){
cout << "n " << n << endl;
cout << "ux " << ux << endl;
cout << "fx " << fx << endl;
cout << "fxe0 " << fxe0 << endl;
cout << "fxe1 " << fxe1 << endl;
cout << "fxe2 " << fxe2 << endl;
cout << "fxe3 " << fxe3 << endl;
cout << "forceleft " << forceleft << endl;
cout << "forceright " << forceright << endl;
cout << "pn[n].tx " << pn[n].tx << endl;
cout << "uy " << uy << endl;
cout << "fy " << fy << endl;
cout << "fye0 " << fye0 << endl;
cout << "fye1 " << fye1 << endl;
cout << "fye2 " << fye2 << endl;
cout << "fye3 " << fye3 << endl;
cout << "forceup " << forceup << endl;
cout << "forcedown " << forcedown << endl;
cout << "pn[n].ty " << pn[n].ty << endl;
cout << "angle " << atan(pn[n].ty/pn[n].tx)*180/3.14 << endl;
}
*/


}

}

	





