#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
void calc_Kdotx(int* kcol, double* kval, double* diag, double* x, double* b, int nrrdof)
{
	int r, a, lim;

	for(r=0;r<nrrdof;r++)
		b[r] = diag[r]*x[r];

	for(r=0;r<nrrdof;r++)
	{
		lim = 10*r+kcol[10*r];
		for(a=10*r+1;a<lim;a++)
		{
			b[r]+=kval[a]*x[kcol[a]];
			b[kcol[a]]+=kval[a]*x[r];
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
void solvePCG(int* kcol, double* kval, double* u, double* f, int nrrdof)
{
	int i, a, iter;
	double *ui, *ri, *diag, *invC, *zi, *pi, *qi;
	double rhoi, rhoinew, initrho;
	double beti, alfi, pq;

	//double rhoreal = 1E6;
	//initrho = 1E6;

	//while(rhoreal > par.ACCURACY*initrho){
	ui = new double[nrrdof];//calloc(nrrdof, sizeof(double)); 
	for(i=0;i<nrrdof;i++){ui[i]=u[i];} //u[i];
	ri   = new double[nrrdof];
	diag = new double[nrrdof];
	invC = new double[nrrdof];
	zi   = new double[nrrdof];
	pi   = new double[nrrdof];
	qi   = new double[nrrdof];

	for(i=0;i<nrrdof;i++)
		diag[i] = kval[10*i];

	for(i=0;i<nrrdof;i++) // for each row in K
	{
		if(diag[i] != 0.0) // if Kii not zero
			invC[i] = 1.0/diag[i];  					// invC = inv(diag(K))
		else
			invC[i] = 0.0;
	}

	calc_Kdotx(kcol,kval,diag,ui,qi,nrrdof);
	for(i=0;i<nrrdof;i++)			// r0 = f-K*u0
		ri[i]=f[i]-qi[i];

	for(i=0;i<nrrdof;i++)
		zi[i] = invC[i]*ri[i]; 							// z0 = inv(C)*r0

	for(i=0;i<nrrdof;i++)
		pi[i]=zi[i]; 									// p0 = z0

	for(i=0,rhoinew=.0;i<nrrdof;i++)
		rhoinew += ri[i]*zi[i]; 						// rhoi = zi*ri

	for(i=0,initrho=.0;i<nrrdof;i++)
		initrho += invC[i]*f[i]*f[i]; 				// FOR (par.ACCURACY)

	// start iterative solve
	for(iter=0; (rhoinew>(par.ACCURACY)*initrho); iter++) // error due to rounding errors?
	{

if(iter>par.MAXNRITER){break;}

		rhoi = rhoinew;

		calc_Kdotx(kcol,kval,diag, pi, qi, nrrdof); 	// qi = K*pi


		for(i=0,pq=0;i<nrrdof;i++)
			pq+=pi[i]*qi[i];
		alfi = rhoi/pq;									// alfi = rhoi/(pi*qi)

		for(i=0;i<nrrdof;i++)
	   		ui[i]+=alfi*pi[i];							// ui+1 = ui+alfi*pi

		for(i=0;i<nrrdof;i++){
		 	ri[i]-=alfi*qi[i];}			// ri+1 = ri-alfi*qi


		for(i=0;i<nrrdof;i++)
			zi[i] = invC[i]*ri[i]; 						// zi+1 = inv(C)*ri+1

		for(i=0,rhoinew=.0;i<nrrdof;i++)
			rhoinew += ri[i]*zi[i]; 					// rhoi+1 = ri+1*zi+1

		beti = rhoinew/rhoi;							// beti = rhoinew/rhoi

		for(i=0;i<nrrdof;i++)
			pi[i] = zi[i] + beti*pi[i]; 				// pi+1 = zi+1 + betai*pi

		//if(iter%10==0)
		//	printf("\ni %4d, rhoinew/initrho=%18.11lf",iter, rhoinew/initrho);
	}

	printf("\n Stop iterating at iter %d",iter);
	for(i=0;i<nrrdof;i++)
		u[i]=ui[i];

	//double *Kureal; 
	//Kureal=new double[nrrdof];
	//calc_Kdotx(kcol,kval,diag,u,Kureal,nrrdof);
	//double *rreal;
	//rreal=new double[nrrdof];
	//for(i=0;i<nrrdof;i++)
	//  	rreal[i]=f[i]-Kureal[i];
	//for(i=0,rhoreal=0;i<nrrdof;i++)
	//   	rhoreal+=rreal[i]*rreal[i]*invC[i];

	//}

	//free(ui);free(ri);free(diag);free(invC);free(zi);free(pi);free(qi);
	delete [] ui; delete [] ri; delete [] diag; delete [] invC; delete [] zi; delete [] pi; delete [] qi;
}
