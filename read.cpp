// file: read.c
// If the simulation was interrupted, the functions in load.c can retreive
// necessary output data to restart the simulation.
// The structure of these functions is often similar to those in write.c.

#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
int read_increment(void)
{
	int incr;
	FILE *ifp;

	ifp = fopen("increment_number.out","r");

	if(ifp)
	{
		fscanf(ifp,"%d",&incr);
		incr++; // because this increment was finished when written down
		fflush(ifp); fclose(ifp);
	}
	else  // if file does not exist fopen returns NULL and incr becomes 0
		incr = 0;

	return incr;
}

////////////////////////////////////////////////////////////////////////////////
int read_cells(VOX* pv, int increment,string sreadcells)
{
   	int v = 0;
   	int anint, NRc;
   	FILE *ifp;

	const char * filename = sreadcells.c_str();

   	ifp = fopen(filename,"r");
	while((v<NV) && (fscanf(ifp,"%d",&anint)==1))
   		pv[v++].ctag = anint;
   	if(v != NV)
   		printf("\nERROR while loading densities");
   	fflush(ifp); fclose(ifp);

	for(v=0,NRc=0;v<NV;v++)
		NRc = (pv[v].ctag>NRc)?pv[v].ctag:NRc;

	return NRc;
}
