// file: cpm_dh.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
double calcdH(VOX* pv, NOD* pn, int* csize, int xt, int xs, int pick, int ttag, int stag, int incr, double* dens, bool matrix,int NRc,int* celltypes)
{


	double dH, dHcontact, dHvol, dHadh;
	dH=0;


		if(!matrix) 
		{


			dHcontact = calcdHcontact(pv,xt,ttag,stag,celltypes);
			dHvol = calcdHvol(csize,ttag,stag,celltypes);
			dHadh = calcdHadh(pn,xt,xs,pick,ttag,stag,dens,csize,pv); 

			//let cells grow first
			if(incr<par.RELAXTIME)
			{
				double dHvolnews=0; double dHvololds=0;
				double dHvolnewt=0; double dHvololdt=0;
				double Atarget = 300;
				if(stag){dHvolnews=(csize[stag-1]+1-Atarget)*(csize[stag-1]+1-Atarget);}
				if(stag){dHvololds=(csize[stag-1]-Atarget)*(csize[stag-1]-Atarget);}
				if(ttag){dHvolnewt=(csize[stag-1]-1-Atarget)*(csize[stag-1]-1-Atarget);}
				if(ttag){dHvololdt=(csize[ttag-1]-Atarget)*(csize[ttag-1]-Atarget);}

				dH=+dHvolnews-dHvololds+dHvolnewt-dHvololdt;
				dH=dH/80;
//cout << "dHcontact " << dHcontact << endl;
//cout << "dH " << dH << endl;

			//dH+=dHcontact;

			}
			

			if(incr>=par.RELAXTIME)
			{
			dH+=dHvol;
			dH+=dHcontact;
			dH+=dHadh;

			}


	//cout << "stag " << stag << endl;
	//cout << "ttag " << ttag << endl;
	if(ttag)
	{
	if(pick==8)
	{

	dH=0; //TEST


	//cout << "middlecell " << endl;	cout << "dHcontact " << dHcontact << endl;
	//cout << "dHvol " << dHvol << endl;
	//cout << "dHadh " << dHadh << endl;
	//cout << "dH " << dH << endl;


	}

 	double estrains[3];
	double estress[3];
      	get_estrains(pn,xt,estrains);
	get_estress(xt, estrains,estress);
	double hstr = (estress[0]+estress[1])/2;
	if(hstr<0){hstr=0;}
	if(hstr>0)
	{
	//cout << endl;
	//cout << "dHcontact " << dHcontact << endl;
	//cout << "dHvol " << dHvol << endl;
	//cout << "dHadh " << dHadh << endl;
	//cout << "dH " << dH << endl;

	}
	}
	if(stag)
	{
	//cout << endl;
	//cout << "extension " << endl;
	//cout << "dHcontact " << dHcontact << endl;
	//cout << "dHvol " << dHvol << endl;
	//cout << "dHadh " << dHadh << endl;
	//cout << "dH " << dH << endl;
	}
	//cout << endl;
		}



		if(matrix) 
		{

				
				dHcontact=0;
				//dHcontact=calcdHcontact(pv,xt,ttag,stag,celltypes);
				dHvol=0;
				if(!par.LEMMONROMER){dHvol = calcdHvol(csize,ttag,stag,celltypes); }

				dH += dHcontact+dHvol; 

			


		}

	return dH;


}




////////////////////////////////////////////////////////////////////////////////
double calcdHcontact(VOX* pv, int xt, int ttag, int stag,int* celltypes)
{


	double dHcontact, Hcontact, Hcontactn;
	int nbh=par.NBHRAD;
	int rad; if(nbh>2){rad=par.NBHRAD;}
	if(nbh<=2){rad=0;}



	if(rad==0)
	{
		if(nbh==2){

		int nbs[20],n,nbtag;

		Hcontact=0; Hcontactn=0;
		nbs[0]=xt-1+(par.NVX); nbs[1]=xt+(par.NVX); nbs[2]=xt+1+(par.NVX);
		nbs[7]=xt-1;                    nbs[3]=xt+1;
		nbs[6]=xt-1-(par.NVX); nbs[5]=xt-(par.NVX); nbs[4]=xt+1-(par.NVX);

		nbs[8]=nbs[6]-par.NVX;
		nbs[9]=nbs[5]-par.NVX;
		nbs[10]=nbs[4]-par.NVX;
		nbs[11]=nbs[4]+1;
		nbs[12]=nbs[3]+1;
		nbs[13]=nbs[2]+1;
		nbs[14]=nbs[2]+par.NVX;
		nbs[15]=nbs[1]+par.NVX;
		nbs[16]=nbs[0]+par.NVX;
		nbs[17]=nbs[0]-1;
		nbs[18]=nbs[7]-1;
		nbs[19]=nbs[6]-1;




		for(n=0;n<20;n++)
		{
			nbtag = pv[nbs[n]].ctag; 
			if(nbs[n]>=0&nbs[n]<=par.NVX*par.NVY)
			{
				Hcontact += contactenergy(ttag,nbtag,celltypes);
				Hcontactn += contactenergy(stag,nbtag,celltypes);
			}
		}
		dHcontact = Hcontactn-Hcontact;}


		if(nbh==1){
		int nbs[8],n,nbtag;

		Hcontact=0; Hcontactn=0;
		nbs[0]=xt-1+(par.NVX); nbs[1]=xt+(par.NVX); nbs[2]=xt+1+(par.NVX);
		nbs[7]=xt-1;                    nbs[3]=xt+1;
		nbs[6]=xt-1-(par.NVX); nbs[5]=xt-(par.NVX); nbs[4]=xt+1-(par.NVX);

		for(n=0;n<8;n++)
		{
			nbtag = pv[nbs[n]].ctag; 
			if(nbs[n]>=0&nbs[n]<=par.NVX*par.NVY)
			{
				Hcontact += contactenergy(ttag,nbtag,celltypes);
				Hcontactn += contactenergy(stag,nbtag,celltypes);
			}
		}

		dHcontact = Hcontactn-Hcontact;}
	}
	
	if(rad>0)
	{

		Hcontact=0; Hcontactn=0;

		int y = xt/(par.NVX); int x = xt%(par.NVX);




		for(int vy=y-rad; vy<=y+rad; vy++) 
		{
			for(int vx=x-rad; vx<=x+rad; vx++) 
			{
					int v=vx+ vy*par.NVX; 

				if((x-vx)*(x-vx)+(y-vy)*(y-vy)<=rad*rad)
				{


						if(vx>0 & vx<par.NVX & vy>0 & vy<par.NVY){int nbtag = pv[v].ctag;
						Hcontact += contactenergy(ttag,nbtag,celltypes);
						Hcontactn += contactenergy(stag,nbtag,celltypes);}
						if(vx<0 | vx>par.NVX | vy<0 | vy>par.NVY){int nbtag = 0;
						Hcontact += contactenergy(ttag,nbtag,celltypes);
						Hcontactn += contactenergy(stag,nbtag,celltypes);}
				
				}			
			}
		}


	dHcontact = Hcontactn-Hcontact;


	}


	return dHcontact;
}

////////////////////////////////////////////////////////////////////////////////
double contactenergy(int tag1, int tag2,int* celltypes)
{

	double J;

	J = 0;

	if(tag1!=tag2)
	{
    	if((tag1==0)||(tag2==0))
	{
        	J = JCM;
		if(tag1!=0 && celltypes[tag1-1]==2)
        		J = JCM2;
		if(tag2!=0 && celltypes[tag2-1]==2)
        		J = JCM2;
	} 
    	else
	{

        	J = JCC;

		if(celltypes[tag1-1]==2)
        		{J = JCC2;}
		if(celltypes[tag1-1]!=celltypes[tag2-1])
			{J = JCCb;}

	}
	} 



	return J;


}





////////////////////////////////////////////////////////////////////////////////
double calcdHvol(int* csize, int ttag, int stag,int* celltypes)
{
double dHvol=0;
	if(!par.CLASSICCPM)
	{
		double dHvolA,dHvolB,V0,eV,eVn;
		int V;
		dHvolA=0; dHvolB=0;

		// assume 2 cells, A (ttag) and B (stag) are involved
		if(ttag) // cell ttag retracts
		{
			dHvolA = -par.INELASTICITY;		
				if(celltypes[ttag-1]==2){dHvolA = -par.INELASTICITY2;}
		}
		if(stag) // cell stag expands
		{
			dHvolB = par.INELASTICITY;	
				if(celltypes[stag-1]==2){dHvolB = par.INELASTICITY2;}
		}
		dHvol = dHvolA+dHvolB;



	}
	if(par.CLASSICCPM)
	{
		double dHvolA,dHvolB,V0,eV,eVn;
		int V;

		// assume 2 cells, A (ttag) and B (stag) are involved
		dHvolA=0; dHvolB=0; V0=par.TARGETVOLUME;
		if(ttag) // cell ttag retracts
		{
			V=csize[ttag-1]; eV=(V-V0); eVn=(V-1-V0);
			dHvolA = (par.INELASTICITY)*(eVn*eVn-eV*eV);
		}
		if(stag) // cell stag expands
		{
			V=csize[stag-1]; eV=(V-V0); eVn=(V+1-V0);
			dHvolB = (par.INELASTICITY)*(eVn*eVn-eV*eV);

		}
		dHvol = dHvolA+dHvolB;


	}
	
		return dHvol;

}


double calcdHadh(NOD* pn, int xt, int xs, int pick, int ttag, int stag, double *dens, int* csize, VOX* pv)
{
double dHadh=0;
	//double Aref=2000;//200
	double Aref=par.MAXAREA;
	//nog niet dens gebruiken
	double dHadh1 =0;
	double adh1=0;
	double adh1new=0;
	double dHadh2 = 0;
	double adh2 = 0;
	double adh2new =0;
	if(stag)
	{
		adh1=csize[stag-1];
		adh1new=csize[stag-1]+1;
	}
	if(ttag)
	{
		adh2=csize[ttag-1];
		adh2new=csize[ttag-1]-1;
	}
	dHadh1=par.LAMBDAADHESION*adh1/(Aref+adh1)-par.LAMBDAADHESION*adh1new/(Aref+adh1new);
	dHadh2=par.LAMBDAADHESION*adh2/(Aref+adh2)-par.LAMBDAADHESION*adh2new/(Aref+adh2new);
	dHadh = dHadh1+dHadh2;



	return dHadh;

}




