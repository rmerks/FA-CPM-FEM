// file: cellmoves.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
void CPM_moves(VOX* pv, NOD* pn, int* csize, int *csumx, int *csumy, int incr,double *FA,int NRc,int* celltypes,double* sumFA,NOD* pnold)
// cellular potts model: one Monte Carlo step
{
	int i,NRsteps = NV;
	int xs, xt; // source and target pixel
	int xtx,xty; // x and y position of target pixel
	int xsx,xsy; // x and y position of source pixel
	int ttag, stag; // target and source label
	int nbs[8],pick; // neighbors of target pixel
	BOOL go_on;
	double dH, prob;
	int fc;



	for(i=0;i<NRsteps;i++)
	{

		xt = mt_random()%NV; // pick random element
		xty = xt/(par.NVX); xtx = xt%(par.NVX);


		if((xtx>par.FORBIDDENZONE)&&(xtx<(par.NVX)-par.FORBIDDENZONE-1)&&(xty>par.FORBIDDENZONE)&&(xty<par.NVY-par.FORBIDDENZONE-1)) // exclude outer rim
		{


			nbs[0]=xt-1+(par.NVX); nbs[1]=xt+(par.NVX); nbs[2]=xt+1+(par.NVX);
			nbs[7]=xt-1;                    nbs[3]=xt+1;
			nbs[6]=xt-1-(par.NVX); nbs[5]=xt-(par.NVX); nbs[4]=xt+1-(par.NVX);

			ttag = pv[xt].ctag;


			int middlecell=0;

			pick = mt_random()%8;




			if(pick<8){
			xs = nbs[pick]; stag = pv[xs].ctag; }// pick random neighbor}
			if(pick==8){
			stag=0;
			pick = mt_random()%8;
			xs=nbs[pick]; 
			pick=8;	//test
			}

			xsy = xs/(par.NVX); xsx = xs%(par.NVX);


			go_on = 0; 
				
			if(ttag!=stag) //don't bother if no difference
			{
				go_on = 1; 




				if(ttag) // if a cell in xt (retracting)
				{
			    		if (splitcheckCCR(pv,csize,xt,ttag))
						go_on = 0;
			
	
			    		if(csize[ttag-1]==1) // cell cannot disappear (constraint may be removed)
						go_on = 0;
				}


				for(fc = 1;fc<=par.NRcf;fc++) //fixed cells
				{
					if(ttag==fc){go_on=0;}	
					if(stag==fc){go_on=0;}	
				}

			}


			if(go_on)
			{



				dH = calcdH(pv,pn,csize,xt,xs,pick,ttag,stag,incr,FALSE,NRc,celltypes);



				if(ttag>0 & FA[xt]>0)//
				{

				    	double estrains[3];
					double estress[3];
				      	get_estrains(pnold,xt,estrains);
					get_estress(xt, estrains,estress);
					double hstr = (estress[0]+estress[1])/2;
					if(hstr<0){hstr=0;}

					double dHdiss=par.INELASTICITY/2;
					dHdiss = 0;

					if(FA[xt]>par.BASEFA) 
					{
					double maxfa = 50*par.VOXSIZE*par.VOXSIZE/8e-15; 

					double p = par.LAMBDAPLAQUE;
					double k = par.CONFSTRESS;
					dHdiss = par.LAMBDAFA*FA[xt]*(1+p*hstr/(k+hstr)); 
					dHdiss = par.LAMBDAFA*(FA[xt]-par.BASEFA)*(1+p*hstr/(k+hstr)); 
					dHdiss = par.LAMBDAFA*(FA[xt]-par.BASEFA)/(par.FAH-par.BASEFA+FA[xt]-par.BASEFA)*(1+p*hstr/(k+hstr)); 

					
					dH+=dHdiss;

					}


	
				}
				


				prob = exp(-(dH)/(par.MOTILITY));
				if(incr<par.RELAXTIME){prob=exp(-(dH)/(0.5));}

				if (prob>(rand()/(double)RAND_MAX))
				{

					if(stag==ttag){cout << "BEW HIER " << xt<< " FA[xt] " << FA[xt] << endl;}

					if(pick<8){
			    		if(ttag) {
				        csize[ttag-1]--;
				        csumx[ttag-1]-=xtx;
				        csumy[ttag-1]-=xty;
				    	}
			    		if(stag) {
				        csize[stag-1]++;
				        csumx[stag-1]+=xtx;
				        csumy[stag-1]+=xty;
					
				   	}}



					if(pick<8){pv[xt].ctag = stag;} // a move is made
					if(ttag)
					{	


						sumFA[ttag-1]=sumFA[ttag-1]-FA[xt];
						FA[xt]=0; 


					}
					if(stag)
					{

						if(sumFA[stag-1]+par.BASEFA<=par.CAPACITYFA)
						{
							FA[xt]=par.BASEFA;
							sumFA[stag-1]=sumFA[stag-1]+par.BASEFA;
						}
						else
						{
							FA[xt]=par.CAPACITYFA-sumFA[stag-1];
							sumFA[stag-1]=sumFA[stag-1]+FA[xt];
						}
					

					}


                    
				}
	
			} 
		}
	}

}

////////////////////////////////////////////////////////////////////////////////
BOOL splitcheckCCR(VOX* pv, int* csize, int xt, int ttag)
{
	BOOL split;
	int nbs[8],n,nb,prev,curr,in;
	int v, nrblue, nrgrey, startnb;
	int greys[csize[ttag-1]];
	short CCAlabels[NV];
	int i, nrgrey0, g, nbsg[8];

	nbs[0]=xt-1+(par.NVX); nbs[1]=xt+(par.NVX); nbs[2]=xt+1+(par.NVX);
	nbs[7]=xt-1;                    nbs[3]=xt+1;
	nbs[6]=xt-1-(par.NVX); nbs[5]=xt-(par.NVX); nbs[4]=xt+1-(par.NVX);

	prev = pv[nbs[7]].ctag; in = 0;
	for(n=0;n<8;n++)
	{
		curr = pv[nbs[n]].ctag;
		if((prev!=ttag)&&(curr==ttag))
			in++;
		prev = curr;
	}

	split = FALSE;
	if(in>1)
	{
	// CONNECTED COMPONENT ALGORITHM Rene(van Oers)-style (CCR)
    	// connected checking "label":
    	// 0: blue;    neighbors of retracted element
    	// 1: white;   undiscovered
    	// 2: grey;    discovered but not finished processing
    	// 3: black;   finished processing

		for(v=0;v<NV;v++) {CCAlabels[v] = 1;}
		CCAlabels[xt] = 3;

		nrblue = -1;
		for(n=0;n<8;n++)
		{
			nb = nbs[n];
			if(pv[nb].ctag==ttag)
			{
				CCAlabels[nb]=0; nrblue++;
				startnb = nb;
			}
		}
		CCAlabels[startnb]=2; nrgrey=1; greys[0]=startnb;

		while(nrgrey&&nrblue)
		{
			nrgrey0 = nrgrey;
			// make neighbors of discovered greys grey
			for(i=0;i<nrgrey0;i++)
			{
				g = greys[i];
				nbsg[0]=g-1+(par.NVX); nbsg[1]=g+(par.NVX); nbsg[2]=g+1+(par.NVX);
				nbsg[7]=g-1;                    nbsg[3]=g+1;
				nbsg[6]=g-1-(par.NVX); nbsg[5]=g-(par.NVX); nbsg[4]=g+1-(par.NVX);
				for(n=0;n<8;n++)
				{
					nb = nbsg[n];
					if((pv[nb].ctag==ttag)&&(CCAlabels[nb]<2))
					{
						if(CCAlabels[nb]==0) {nrblue--;}
						CCAlabels[nb]=2; nrgrey++; greys[nrgrey-1]=nb;
					}
				}
			}

			// make processed greys black
			for(i=0;i<nrgrey0;i++)
			{
				g = greys[i];
				CCAlabels[g]=3;
				greys[i]=greys[nrgrey-1]; nrgrey--;
			}

		}
		if(nrblue) {split = TRUE;}

	}
	return split;
}

void CalcPerimeters(VOX* pv, int* cper, int NRc){

	bool already;
	for(int c=0;c<NRc;c++) {cper[c]=0;}

	for(int v = 0;v<NV;v++)
	{
		already=0;
		int ct = pv[v].ctag;
		if(ct)
		{
			int vy = v/(par.NVX); int vx = v%(par.NVX);
			//determine neighbouring voxels
			//int vxn[8] = {0,1,1,1,0,-1,-1,-1};
			//int vyn[8] = {1,1,0,-1,-1,-1,0,1};
			int vxn[4] = {0,1,0,-1};
			int vyn[4] = {1,0,-1,0};
			//check if one of these neighbours has a different cell type
			//if yes, this voxel is on the perimeter of the cell
			for(int i = 0;i<4;i++)
			{
			
				int vn = vx+vxn[i] + vy*par.NVX + vyn[i]*par.NVX;
				if((vx+vxn[i])>=0 && (vy+vyn[i])>=0 && (vy+vyn[i]) < par.NVY && (vx+vxn[i]) < par.NVX)
				{
					int ctn = pv[vn].ctag;
					if(!(ct==ctn) && !already)
					{	
						cper[ct-1]++;already=1;
					}
				}
			}
		}
	}

}


void CalcLengths(VOX* pv,NOD* pn,double* clength,double* ecc,double* cangle,int NRc,int* csize,int* csumx, int* csumy){
	for(int c=0;c<NRc;c++)
	{
		double Ixx = 0;    
		double Iyy = 0;
		double Ixy = 0;
		double lambda = 0;
		double lambda_small = 0;
		double length = 0;
		double e = 0;
		double shortlength =0;
		double cmx = double(csumx[c])/csize[c];
		double cmy = double(csumy[c])/csize[c];
		double inertia[3],v1[2],v2[2];


		for(int v = 0;v<NV;v++)
		{

			if(pv[v].ctag == c+1)
			{
				int vy = v/(par.NVX); int vx = v%(par.NVX);
				Ixx += (vx-cmx)*(vx-cmx);
				Iyy += (vy-cmy)*(vy-cmy);
				Ixy -= (vy-cmy)*(vx-cmx);
			}

		}

	inertia[0] =  Ixx;
	inertia[1] =  Iyy;
	inertia[2] =  Ixy;



	lambda=lambda_small=.0; get_princs(inertia,&lambda,&lambda_small,v1,v2,0);
			//lambda = 0.5*(Ixx+Iyy)+0.5*sqrt((Ixx+Iyy)*(Ixx+Iyy)+4*(Ixy*Ixy-Ixx*Iyy));
			//lambda_small = 0.5*(Ixx+Iyy)-0.5*sqrt((Ixx+Iyy)*(Ixx+Iyy)+4*(Ixy*Ixy-Ixx*Iyy));
			length = 4*sqrt(lambda/(double)csize[c]);

			shortlength = 4*sqrt(lambda_small/(double)csize[c]);
			e=sqrt(1-(shortlength/length)*(shortlength/length));
		clength[c]=length;	
		ecc[c]=e;
		double PI = 3.14159265;
		cangle[c] = atan2(v1[1],v1[0])*180/PI;


	}


	if(par.TWOCELL & par.WTWOCELLANGLECM)
	{
		double cmx1 = double(csumx[0])/csize[0];
		double cmy1 = double(csumy[0])/csize[0];
		double cmx2 = double(csumx[1])/csize[1];
		double cmy2 = double(csumy[1])/csize[1];
		double PI =3.14159265;
		if(cmx1==cmx2 & cmy1==cmy2){cangle[2]=0;}
		if(cmx1<cmx2){cangle[2]=-atan2(cmy2-cmy1,cmx2-cmx1)*180/PI;}
		if(cmx2<cmx1){cangle[2]=-atan2(cmy1-cmy2,cmx1-cmx2)*180/PI;}

	}
}

int check_contact(VOX* pv)
{
	int contact;
	int v,vx,vy;
	int nbs[8],n,nb;

	contact = 0;

	for(vy=1; vy<par.NVY-1; vy++)
   	for(vx=1; vx<par.NVX-1; vx++)
   	{
   		v = vx + vy*par.NVX;

		if(pv[v].ctag==1) // first cell
		{
			nbs[0]=v-1+par.NVX; nbs[1]=v+par.NVX; nbs[2]=v+1+par.NVX;
			nbs[7]=v-1;                   nbs[3]=v+1;
			nbs[6]=v-1-par.NVX; nbs[5]=v-par.NVX; nbs[4]=v+1-par.NVX;

			for(n=0;n<8;n++)
			{
				nb = nbs[n];
				if(pv[nb].ctag==2) // second cell
					contact = 1;
					
			}
		}
	}

	return contact;
}




