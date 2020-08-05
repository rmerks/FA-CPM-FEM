#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
/*
int** set_topology(void)
{
	int **top;
	int v, vx, vy;
	int n00, n10, n11, n01;

	top = calloc(NV,sizeof(int*));
	for(v=0;v<NV;v++)
		top[v] = calloc(8,sizeof(int));

	// set topology
   	for(vy=0; vy<par.NVY; vy++)
   	for(vx=0; vx<(par.NVX); vx++)
   	{
   		v = vx + vy*(par.NVX);

		// determine corner node numbers of this voxel
		n00 = (vx  ) + (vy  )*NNX;
		n10 = (vx+1) + (vy  )*NNX;
		n11 = (vx+1) + (vy+1)*NNX;
		n01 = (vx  ) + (vy+1)*NNX;

		top[v][0] = 2*n00;
		top[v][1] = 2*n00+1;
   		top[v][2] = 2*n10;
		top[v][3] = 2*n10+1;
		top[v][4] = 2*n11;
		top[v][5] = 2*n11+1;
		top[v][6] = 2*n01;
		top[v][7] = 2*n01+1;
	}

	return top;
}
*/

////////////////////////////////////////////////////////////////////////////////
void assembly(int* kcol, double* kval, double** klocal, VOX* pv)
{
	int v, vx, vy;
	double value;
	int il, jl, ig, jg; // row i & column j in klocal and K
	int d, a, b, lim;
	BOOL alreadynonzero;
	int n00, n10, n11, n01;
	int topv[8];
	double Ef; // multiply klocal with Ef depending on local E

	for(d=0;d<NDOF;d++)
	{
		kcol[10*d] = 1;
		kval[10*d] = .0;
	}

	for(vy=0; vy<par.NVY; vy++)
   	for(vx=0; vx<(par.NVX); vx++)
	{
		Ef = 1;
		if(par.DUROTAXIS)
		{


		Ef=(1+par.GRADIENT*(vx-100)/par.YOUNGS);
		if(Ef<1/par.YOUNGS){Ef=1/par.YOUNGS;}

		}



		if(par.TESTLOCALSINE)
		{

			Ef=Ef*(1+par.SINEAMP*cos((vx-par.NVX/2)*2*M_PI/par.SINEPER)*cos((vy-par.NVY/2)*2*M_PI/par.SINEPER)/par.YOUNGS);
			if(Ef<1/par.YOUNGS){Ef=1/par.YOUNGS;}

		}

		if(par.YOUNGSNOISE)
		{

			double r = rand()/(double)RAND_MAX; //between 0 and 1

			r=-1+2*r; //between -1 and 1

			r =1+(par.YOUNGSNOISE/par.YOUNGS)*r;
			Ef=Ef*r;
			if(Ef<0){Ef=0;}



		}
		// determine corner node numbers of this element
		n00 = (vx  ) + (vy  )*NNX;
		n10 = (vx+1) + (vy  )*NNX;
		n11 = (vx+1) + (vy+1)*NNX;
		n01 = (vx  ) + (vy+1)*NNX;

		topv[0] = 2*n00;
		topv[1] = 2*n00+1;
   		topv[2] = 2*n10;
		topv[3] = 2*n10+1;
		topv[4] = 2*n11;
		topv[5] = 2*n11+1;
		topv[6] = 2*n01;
		topv[7] = 2*n01+1;


		// place klocal in K matrix
		for(il=0;il<8;il++) // go through rows in klocal
		for(jl=0;jl<8;jl++) // go through columns in klocal
		{
			value = Ef*klocal[il][jl];
			ig = topv[il]; // row in K
			jg = topv[jl]; // column in K

			if(jg==ig) // if on diagonal
				kval[10*ig] += value;

			if(jg>ig) // if right of diagonal
			{
				lim = 10*ig+kcol[10*ig];

				// check if there was already a nonzero on K(ig,jg)
				alreadynonzero = FALSE;
				for(a=10*ig+1;a<lim;a++) // go over ig-th row of K
				{
					if(kcol[a]==jg) // if storage for jg-th column of K
					{
						alreadynonzero = TRUE;
						b = a;
					}
				}

				if(alreadynonzero) // if already a nonzero on K(ig,jg)
					kval[b] += value; // add klocal(il,jl)

				else // if nothing on K(ig,jg)
				{
					b = lim;
					kcol[b] = jg; // make storage for jg-th column of K
					kval[b] = value; // add klocal(il,jl)
					kcol[10*ig]++;
				}
			} //endfor go through klocal
		} //endif relevant element
	} // endfor go though elements
	printf("\nASSEMBLY COMPLETED");
}

////////////////////////////////////////////////////////////////////////////////
int arrange_dofpos(int* dofpos, NOD* pn)
{
	int n, cnt;

	for(n=0,cnt=0;n<NN;n++)
	{
		if(pn[n].restrictx)
			dofpos[2*n] = -1;
		else
			dofpos[2*n] = cnt++;

		if(pn[n].restricty)
			dofpos[2*n+1] = -1;
		else
			dofpos[2*n+1] = cnt++;
	}
	return cnt;
}


////////////////////////////////////////////////////////////////////////////////
void reduce_K(int* kcol, double* kval, int* dofpos)
{
	int ro, co; // old row and column in K
	int rn, cn; // new row and column in K
	int a, lim, shift;

	for(ro=0;ro<NDOF;ro++)
	{
		rn = dofpos[ro];
		if(rn>-1) // if this row is not to be removed
		{
			lim = 10*ro+kcol[10*ro];

			// change column indices for this row:
			for(a=10*ro+1;a<lim;a++)
			{
				co = kcol[a]; // old column index
				cn = dofpos[co];         // new column index
				kcol[a] = cn; // give new column index (some get -1)
			}
			// remove columns with -1 index
			shift = 0;
			for(a=10*ro+1;a<lim;a++)
			{
				kcol[a-shift] = kcol[a];
				kval[a-shift] = kval[a];
			   	if(kcol[a]==-1)
					shift++;
			}
			kcol[10*ro] = kcol[10*ro]-shift;

			// shift row itself
			for(a=0;a<10;a++)
			{
				kcol[10*rn+a] = kcol[10*ro+a];
				kval[10*rn+a] = kval[10*ro+a];
			}
		}
	}
}

