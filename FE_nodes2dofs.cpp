// file: nodelabels.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
void disp_to_nodes(NOD* pn, double* u)
{
	int n, cnt;

	for(n=0,cnt=0; n<NN; n++)
	{

		pn[n].ux=.0;
		pn[n].uy=.0;
		if(!pn[n].restrictx)
		{
			pn[n].ux=u[cnt++];
		}
		if(!pn[n].restricty)
		{
			pn[n].uy=u[cnt++];
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
void kval_to_nodes(double* Km, double* kval,NOD* pn)
{
	int n, cnt;

	for(n=0,cnt=0; n<NN; n++)
	{

		if(!pn[n].restrictx)
		{
			Km[2*n]=kval[10*cnt++];
		}
		if(!pn[n].restricty)
		{
			Km[2*n+1]=kval[10*cnt++];
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
void set_disp_of_prev_incr(NOD* pn, double* u)
{
	int n, cnt;

	for(n=0,cnt=0; n<NN; n++)
	{
		if(!pn[n].restrictx)
			u[cnt++]=pn[n].ux;
		if(!pn[n].restricty)
			u[cnt++]=pn[n].uy;
	}
}

////////////////////////////////////////////////////////////////////////////////
void place_node_forces_in_f(NOD* pn, double* f)
{
	int n, cnt;

	for(n=0,cnt=0; n<NN; n++)
	{
		if(!pn[n].restrictx)
			f[cnt++]=pn[n].fx;
		if(!pn[n].restricty)
			f[cnt++]=pn[n].fy;
	}
}

