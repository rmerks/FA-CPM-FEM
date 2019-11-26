// file: cellmoves.c
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
void CPM_moves(VOX* pv, NOD* pn, int* csize, int *csumx, int *csumy, int incr, double *dens, double *FA,double* FAcell,double* Acell,int NRc)
// cellular potts model: one Monte Carlo step
{
	int i,NRsteps = NV;
	int xs, xt; // source and target pixel
	int xtx,xty; // x and y position of target pixel
	int ttag, stag; // target and source label
	int nbs[8],pick; // neighbors of target pixel
	BOOL go_on;
	double dH, prob;
	int fc;



	for(i=0;i<NRsteps;i++)
	{
		//xt = (rand()*NV/RAND_MAX); // pick random element
		xt = mt_random()%NV; // pick random element
		xty = xt/(par.NVX); xtx = xt%(par.NVX);


		if((xtx>0)&&(xtx<(par.NVX)-1)&&(xty>0)&&(xty<par.NVY-1)) // exclude outer rim: not anymore
		{
			nbs[0]=xt-1+(par.NVX); nbs[1]=xt+(par.NVX); nbs[2]=xt+1+(par.NVX);
			nbs[7]=xt-1;                    nbs[3]=xt+1;
			nbs[6]=xt-1-(par.NVX); nbs[5]=xt-(par.NVX); nbs[4]=xt+1-(par.NVX);

			ttag = pv[xt].ctag;

			if(par.INSERTMEDIUM)
			{
				pick = mt_random()%8;
				//check if in middle of cell
				int middlecell=0;
				int allcell=0;
				for(int pi=0;pi<8;pi++){if(pv[nbs[pi]].ctag>0){allcell++;}}
				for(int pi=0;pi<8;pi++){if(pv[nbs[pi]].ctag==ttag){middlecell++;}}
				if(middlecell<8 & allcell==8){pick=mt_random()%9;}

			}
			if(!par.INSERTMEDIUM){pick = mt_random()%8;}




			if(pick<8){
			xs = nbs[pick]; stag = pv[xs].ctag; }// pick random neighbor}
			if(pick==8){stag=0;pick = mt_random()%8;xs=nbs[pick];}

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


				dH = calcdH(pv,pn,csize,xt,xs,pick,ttag,stag,incr,dens,FALSE,FA,FAcell,Acell);
			
				//Energy that needs to be overcome

				prob = exp(-(dH)/(par.MOTILITY));

				if (prob>(rand()/(double)RAND_MAX))
				{




		double fa=0;
		for(int n=0;n<NN;n++)
		{
		//four surrounding pixels of node, lowerright voxel vx1 vy1
		int vx1=n%NNX; int vy1 = n/NNX;int v1 = vx1 + vy1*(par.NVX);
		int vx2=vx1;int vy2=vy1-1;int v2 = vx2 + vy2*(par.NVX);
		int vx3=vx1-1;int vy3=vy1-1;int v3 = vx3 + vy3*(par.NVX);
		int vx4=vx1-1;int vy4=vy1;int v4 = vx4 + vy4*(par.NVX);


			int nrstag = 0;
			if(v1>=0&v1<NV){if(pv[v1].ctag==1){nrstag++;}}
			if(v2>=0&v2<NV){if(pv[v2].ctag==1){nrstag++;}}		
			if(v3>=0&v3<NV){if(pv[v3].ctag==1){nrstag++;}}		
			if(v4>=0&v4<NV){if(pv[v4].ctag==1){nrstag++;}}
			if(nrstag)
			{
			fa+=FA[n];
			}	

		}
if(FAcell[0] < fa-0.000001 | FAcell[0]>fa+0.000001 ){ break;}

			    		if(ttag) {
				        csize[ttag-1]--;
				        csumx[ttag-1]-=xtx;
				        csumy[ttag-1]-=xty;
				    	}
			    		if(stag) {
				        csize[stag-1]++;
				        csumx[stag-1]+=xtx;
				        csumy[stag-1]+=xty;
				   	}


					//four surrounding nodes of target pixel
					int n00 = (xtx  ) + (xty  )*NNX;
					int n10 = (xtx+1) + (xty  )*NNX;
					int n11 = (xtx+1) + (xty+1)*NNX;
					int n01 = (xtx  ) + (xty+1)*NNX;


					//update FAcell and FA for retraction
					if(ttag>0 & stag==0)
					{
						//check if nodes will be retracted by evaluating four pixels 
						int n2;
						for(int k = 1; k<5;k++)
						{
							if(k==1){n2=n00;}
							if(k==2){n2=n10;}
							if(k==3){n2=n11;}
							if(k==4){n2=n01;}

							//four surrounding voxels of this node
							int vx1=n2%NNX; int vy1 = n2/NNX; int v1 = vx1 + vy1*(par.NVX);
							int vx2=vx1;int vy2=vy1-1;int v2 = vx2 + vy2*(par.NVX);
							int vx3=vx1-1;int vy3=vy1-1;int v3 = vx3 + vy3*(par.NVX);
							int vx4=vx1-1;int vy4=vy1;int v4 = vx4 + vy4*(par.NVX);

							int count = 0; //count how many of these voxels belong to cell ttag
							if(v1>=0&v1<NV){if(pv[v1].ctag==ttag){count++;}}
							if(v2>=0&v2<NV){if(pv[v2].ctag==ttag){count++;}}
							if(v3>=0&v3<NV){if(pv[v3].ctag==ttag){count++;}}
							if(v4>=0&v4<NV){if(pv[v4].ctag==ttag){count++;}}

							//count how many belong to a cell
							int nr = 0;
							if(v1>=0&v1<NV){if(pv[v1].ctag>0){nr++;}}
							if(v2>=0&v2<NV){if(pv[v2].ctag>0){nr++;}}		
							if(v3>=0&v3<NV){if(pv[v3].ctag>0){nr++;}}		
							if(v4>=0&v4<NV){if(pv[v4].ctag>0){nr++;}}

							//if it is just one voxel, this one will be removed, so the focal adhesion will break down
							if(count<=1){FAcell[ttag-1]+= -FA[n2];} //there is a retraction of this node
							if(nr==1){FA[n2]=0;} //whole focal adhesion breaks down because is not attached to other cells


						}
	
	
					}


					//update FAcell and FA for extension
					if(stag>0&ttag==0)
					{
						//check if nodes will be added by evaluating four pixels 
						int n2;
						for(int k = 1; k<5;k++)
						{


							if(k==1){n2=n00;}
							if(k==2){n2=n10;}
							if(k==3){n2=n11;}
							if(k==4){n2=n01;}

							//four surrounding voxels of this node
							int vx1=n2%NNX; int vy1 = n2/NNX; int v1 = vx1 + vy1*(par.NVX);
							int vx2=vx1;int vy2=vy1-1;int v2 = vx2 + vy2*(par.NVX);
							int vx3=vx1-1;int vy3=vy1-1;int v3 = vx3 + vy3*(par.NVX);
							int vx4=vx1-1;int vy4=vy1;int v4 = vx4 + vy4*(par.NVX);

							int count = 0; //count how many of these voxels belong to stag
							if(v1>=0&v1<NV){if(pv[v1].ctag==stag){count++;}}
							if(v2>=0&v2<NV){if(pv[v2].ctag==stag){count++;}}
							if(v3>=0&v3<NV){if(pv[v3].ctag==stag){count++;}}
							if(v4>=0&v4<NV){if(pv[v4].ctag==stag){count++;}}

							//does FA attach to another cell
							int nr =0; int othercells[4];
							if(v1>=0&v1<NV & v1!=xt){if(pv[v1].ctag>0 & pv[v1].ctag!=stag){othercells[nr]=pv[v1].ctag;nr++;}}
							if(v2>=0&v2<NV & v2!=xt){if(pv[v2].ctag>0 & pv[v2].ctag!=stag){othercells[nr]=pv[v2].ctag;nr++;}}
							if(v3>=0&v1<NV & v3!=xt){if(pv[v3].ctag>0 & pv[v3].ctag!=stag){othercells[nr]=pv[v3].ctag;nr++;}}
							if(v4>=0&v4<NV & v4!=xt){if(pv[v4].ctag>0 & pv[v4].ctag!=stag){othercells[nr]=pv[v4].ctag;nr++;}}


//done
							//if all belong to another cell, than it is a new focal adhesion for this cell
							if(count==0)
							{
								if(FA[n2]==0)
								{
									FA[n2]=par.BASEFA; //attached to medium and cell and was zero -> set to baseline
									//also add to the other cells
									bool done[4];done[0]=FALSE;done[1]=FALSE;done[2]=FALSE;done[3]=FALSE;
									for(int oc=0;oc<nr;oc++)
									{
										if(!done[othercells[oc]-1])
										{
											FAcell[othercells[oc]-1]+=FA[n2];
											done[othercells[oc]-1]=TRUE;
										}
									}
								}
							
								FAcell[stag-1]+= FA[n2]; //cell may stick to existing FA, make sure this does not give weird bias
							
							}


							
							
			
						}

					}



			    		pv[xt].ctag = stag; // a move is made


					//update adhesion cells
					if(stag){Acell[stag-1]+= dens[xt];}
					if(ttag){Acell[ttag-1]+= - dens[xt];}





                    
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
		// CONNECTED COMPONENT ALGORITHM Rene-style (CCR)
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




