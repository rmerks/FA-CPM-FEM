// file: write.c
#include "functions.h"
#include <fstream>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
void write_increment(int increment)
{
	FILE *ofp;

	ofp = fopen("increment_number.out","w");
	fprintf(ofp,"%d",increment);
	fflush(ofp); fclose(ofp);
}


////////////////////////////////////////////////////////////////////////////////
void write_cells(VOX* pv,int increment,string sigdir)
{
	int v;
    int vx,vy;
	stringstream str;
	string str1 = sigdir;
	string str2;
	//str2="sigma%05d.out";
	str2="sigmas/sigma%05d.txt";
	str<<str1<<str2;
	string s = str.str();
	const char * fcstr = s.c_str();
        char fname[200];
        sprintf(fname,fcstr,increment);

stringstream ssfname;
string sfname;
ssfname << fname;
ssfname >> sfname;

const char * cfname = sfname.c_str();

	FILE *ofp;

	ofp = fopen(cfname,"a");
        fprintf(ofp,"\n");
    for(vx=0; vx<(par.NVX); vx++) {
        for (vy=0; vy<par.NVY; vy++) {
            v = vx + vy * (par.NVX);
     fprintf(ofp ,"%d ", pv[v].ctag);
        }
        fprintf(ofp, "\n");
    }

   	//fflush(ofp);  
	fclose(ofp);

}


void write_fa(double* FA, int increment, string fadir)
{
	int v;
    int vx,vy;
	stringstream str;
	string str1 = fadir;
	string str2;
	//str2="sigma%05d.out";
	str2="fa/fa%05d.txt";
	str<<str1<<str2;
	string s = str.str();
	const char * fcstr = s.c_str();
        char fname[200];
        sprintf(fname,fcstr,increment);



stringstream ssfname;
string sfname;
ssfname << fname;
ssfname >> sfname;

const char * cfname = sfname.c_str();

	FILE *ofp;





	ofp = fopen(cfname,"a"); 

        fprintf(ofp,"\n");  
    for(vx=0; vx<(par.NVX); vx++) {
        for (vy=0; vy<par.NVY; vy++) {
            v = vx + vy * (par.NVX);
     fprintf(ofp ,"%e ", FA[v]);
        }
        fprintf(ofp, "\n");
    }

   	//fflush(ofp);  
	fclose(ofp);
}
/*
////////////////////////////////////////////////////////////////////////////////
void write_sed(VOX* pv, NOD* pn, int increment)
{
	int v, vx, vy;
	double sed;
   	char filename[20];
   	char astring[20];
   	FILE *ofp;

   	myitostr(increment, astring);
	strcpy(filename, "sed");
   	strcat(filename, astring);
   	strcat(filename, ".out");

	ofp = fopen(filename,"w");

	for(vy=0; vy<par.NVY; vy++)
   	for(vx=0; vx<(par.NVX); vx++)
   	{
   		v = vx + vy*(par.NVX);
		sed = get_sed(vx,vy,pn);
		fprintf(ofp ,"%lf\n",sed);
	}

   	fflush(ofp); fclose(ofp);
}
*/

////////////////////////////////////////////////////////////////////////////////
void write_pstrain(VOX* pv, NOD* pn, int increment)
{
	int v;
	double estrains[3],L1,L2,v1[2],v2[2];
	char filename[20],filename2[20];
   	char astring[20];
   	FILE *ofp;

	myitostr(increment, astring);
	strcpy(filename, "pstrain");
   	strcat(filename, astring);
   	strcat(filename, ".out");

	ofp = fopen(filename,"w");
	for(v=0;v<NV;v++)
	{
		get_estrains(pn,v,estrains);
		L1=L2=.0; get_princs(estrains,&L1,&L2,v1,v2,1);
		if(L1>L2)
		{
			fprintf(ofp,"%d ", (int)(1000000*L1));
			fprintf(ofp,"%d ", (int)(1000*v1[0]));
			fprintf(ofp,"%d ", (int)(1000*v1[1]));
			fprintf(ofp,"%d\n",(int)(1000000*L2));
		}
		else
		{
			fprintf(ofp,"%d ", (int)(1000000*L2));
			fprintf(ofp,"%d ", (int)(1000*v2[0]));
			fprintf(ofp,"%d ", (int)(1000*v2[1]));
			fprintf(ofp,"%d\n",(int)(1000000*L1));
		}
	}
	fflush(ofp); fclose(ofp);
}


/*
////////////////////////////////////////////////////////////////////////////////
void write_pstress(VOX* pv, NOD* pn, int increment)
{
	int v;
	double estrains[3],estress[3],L1,L2,v1[2],v2[2];
	char filename[20];
   	char astring[20];
   	FILE *ofp;

	myitostr(increment, astring);
	strcpy(filename, "pstress");
   	strcat(filename, astring);
   	strcat(filename, ".out");

	ofp = fopen(filename,"w");
	for(v=0;v<NV;v++)
	{
		get_estrains(pn,v,estrains);
		get_estress(v,estrains,estress);
		L1=L2=.0; get_princs(estress,&L1,&L2,v1,v2,0);
		if(L1>L2)
			fprintf(ofp ,"%e\n",L1);
		else
			fprintf(ofp ,"%e\n",L2);
	}
	fflush(ofp); fclose(ofp);
}
*/

////////////////////////////////////////////////////////////////////////////////
void write_forces(NOD* pn, int increment)
{
	int n;
	char filename[20];
   	char astring[20];
	FILE *ofp;

	myitostr(increment, astring);
	strcpy(filename, "forces");
   	strcat(filename, astring);
   	strcat(filename, ".out");

	ofp = fopen(filename,"w");
	for(n=0;n<NN;n++)
		fprintf(ofp ,"%e\n",pn[n].fx);
	for(n=0;n<NN;n++)
		fprintf(ofp ,"%e\n",pn[n].fy);
	fflush(ofp); fclose(ofp);
}


////////////////////////////////////////////////////////////////////////////////
void write_disps(NOD* pn, int increment)
{
	int n;
	char filename[20];
   	char astring[20];
	FILE *ofp;

	myitostr(increment, astring);
	strcpy(filename, "disps");
   	strcat(filename, astring);
   	strcat(filename, ".out");

	ofp = fopen(filename,"w");
	for(n=0;n<NN;n++)
		fprintf(ofp ,"%e\n",pn[n].ux);
	for(n=0;n<NN;n++)
		fprintf(ofp ,"%e\n",pn[n].uy);
	fflush(ofp); fclose(ofp);
}

void write_ratiopa(int* cper,int* csize,int NRc,int increment,FILE* ofp){
	fprintf(ofp,"%d ",increment);
	for(int n=0;n<NRc;n++)
	{
		double ratio = (double)cper[n]/(double)csize[n];
		fprintf(ofp ,"%e ",ratio);
	}
        fprintf(ofp,"\n");
	//fflush(ofp); 

}

void write_length(double* length,int NRc,int increment,FILE* ofp){

	fprintf(ofp,"%d ",increment);
	for(int n=0;n<NRc;n++)
	{
		fprintf(ofp ,"%e ",length[n]);

	}
        fprintf(ofp,"\n");
	//fflush(ofp);
}

void write_area(int* area,int NRc,int increment,FILE* ofp){

	fprintf(ofp,"%d ",increment);
	for(int n=0;n<NRc;n++)
	{
		fprintf(ofp ,"%e ",(double) area[n]);
	}
        fprintf(ofp,"\n");
	//fflush(ofp); 

}

void write_cangle(double* cangle,int NRc,int increment,FILE* ofp){
	if(!par.WTWOCELLANGLECM)
	{
	fprintf(ofp,"%d ",increment);
	for(int n=0;n<NRc;n++)
	{
		fprintf(ofp ,"%e ",cangle[n]);
	}
	//fflush(ofp);
        fprintf(ofp,"\n");
	}
	if(par.TWOCELL & par.WTWOCELLANGLECM)
	{
	fprintf(ofp,"%d ",increment);
	for(int n=0;n<3;n++)
	{
		fprintf(ofp ,"%e ",cangle[n]);
	}
	//fflush(ofp);
        fprintf(ofp,"\n");
	}
	
}


void write_sqdis(double* sqdis,int NRc,int increment,FILE* ofp){

	fprintf(ofp,"%d ",increment);
	for(int n=0;n<NRc;n++)
	{
		fprintf(ofp ,"%e ",sqdis[n]);
	}
        fprintf(ofp,"\n");
	//fflush(ofp);

}

void write_eccentricity(double* ecc,int NRc,int increment,FILE* ofp){

	fprintf(ofp,"%d ",increment);
	for(int n=0;n<NRc;n++)
	{
		fprintf(ofp ,"%e ",ecc[n]);
	}
	//fflush(ofp);      
        fprintf(ofp,"\n");

}

void write_twocellcontact(int contact,int increment,FILE* ofp){

	fprintf(ofp,"%d ",increment);
	fprintf(ofp ,"%e ",(double) contact);
        fprintf(ofp,"\n");
	//fflush(ofp); 

}

void write_totshape(int* totshape,int increment,string totshapedir){

	int v;
    	int vx,vy;
	stringstream str;
	string str1 = totshapedir;
	string str2;
	//str2="sigma%05d.out";
	str2="totshape/totshape%05d.txt";
	str<<str1<<str2;
	string s = str.str();
	const char * fcstr = s.c_str();
        char fname[200];
        sprintf(fname,fcstr,increment);

stringstream ssfname;
string sfname;
ssfname << fname;
ssfname >> sfname;

const char * cfname = sfname.c_str();

	FILE *ofp;

	ofp = fopen(cfname,"a");
        fprintf(ofp,"\n");
    for(vx=0; vx<2*(par.NVX)+1; vx++) {
        for (vy=0; vy<2*par.NVY+1; vy++) {
            v = vx + vy*(2*(par.NVX)+1); 
     fprintf(ofp ,"%d ", totshape[v]); 
        }
        fprintf(ofp, "\n");
    }

   	//fflush(ofp);  
	fclose(ofp);

}
