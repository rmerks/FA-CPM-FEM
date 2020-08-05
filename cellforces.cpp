// file: par.CELLFORCEs.c

#include <functions.h>
#include <math.h>
//#include <boost/math/tools/minima.hpp>

#include <vector>
#include <cmath>
#include <iostream>

void cell_forces(VOX* pv, NOD* pn, int* csize, int NRc,int* csumx,int* csumy, double* FA, NOD* pnold)
{


	int c;
	int n,nx,ny;
	int v,vx,vy, cnttag;
	int NRcelln,cellnodes[NN];
	int i,j, n2;
	double dnx,dny,forcex,forcey;




		for(c=0;c<NRc;c++)//NRc
		{


		double cmx = double(csumx[c])/csize[c];
		double cmy = double(csumy[c])/csize[c];

			// determine which nodes belong to cell c
			NRcelln = 0;
			for(ny=1; ny<NNY-1; ny++)
	   		for(nx=1; nx<NNX-1; nx++)
	   		{
	   			n = nx + ny*NNX;
				cnttag = 0;
				for(vy=ny-1; vy<ny+1; vy++)
				for(vx=nx-1; vx<nx+1; vx++)
				{
					v = vx + vy*(par.NVX);
					if(pv[v].ctag == c+1)
						cnttag++;
				}
				if(cnttag>0) // all cell nodes
				{
					cellnodes[NRcelln] = n;
					NRcelln++;
				}
			}


			// forces between cellnodes
			for(i=0;i<NRcelln;i++)
			{
				n = cellnodes[i];
				int count =0;


				for(j=0;j<NRcelln;j++)
				{



					n2 = cellnodes[j];
					ny=n/NNX; nx=n%NNX;
					dny=(n2/NNX-ny)*(par.VOXSIZE); // y distance between n and n2
					dnx=(n2%NNX-nx)*(par.VOXSIZE); // x distance between n and n2
					//check if cell nodes n and n2 are connected by a straight line that stays within the cell c

					if(par.NODECONNECTION)
					{

						if(CheckCellNodeConnection(pv,pn,n,n2,c))
						{


							count++;


							if(par.CLASSICCPM)
							{
								forcex = par.LRTENSION*dnx/double(csize[c])/par.THICKNESS;
								forcey = par.LRTENSION*dny/double(csize[c])/par.THICKNESS;

		
							}

							if(!par.CLASSICCPM)
							{
								forcex = par.LRTENSION*dnx/double(csize[c]*csize[c])/par.THICKNESS;
								forcey = par.LRTENSION*dny/double(csize[c]*csize[c])/par.THICKNESS;


							}



						pn[n].fx += forcex;
						pn[n].fy += forcey;



						}


					}
					else
					{



							if(par.CLASSICCPM)
							{
								forcex = par.LRTENSION*dnx/double(csize[c])/par.THICKNESS;
								forcey = par.LRTENSION*dny/double(csize[c])/par.THICKNESS;

							}

							if(!par.CLASSICCPM)
							{
								forcex = par.LRTENSION*dnx/double(csize[c]*csize[c])/par.THICKNESS;
								forcey = par.LRTENSION*dny/double(csize[c]*csize[c])/par.THICKNESS;
							}


						pn[n].fx += forcex;
						pn[n].fy += forcey;
					}

				
				}



				if(par.FORCEFA)
				{

					//four surrounding voxels of this node
					int vx1=n%NNX; int vy1 = n/NNX; int v1 = vx1 + vy1*(par.NVX);
					int vx2=vx1;int vy2=vy1-1;int v2 = vx2 + vy2*(par.NVX);
					int vx3=vx1-1;int vy3=vy1-1;int v3 = vx3 + vy3*(par.NVX);
					int vx4=vx1-1;int vy4=vy1;int v4 = vx4 + vy4*(par.NVX);


						double k = par.CONFSTRESS;
						double hstr=0;

						double sumfactor=0;
						double estrains[3];
						double estress[3];

					      	get_estrains(pnold,v1,estrains);
						get_estress(v1, estrains,estress);
						hstr = (estress[0]+estress[1])/2;

						if(hstr>0){sumfactor+=hstr/(k+hstr); }

					      	get_estrains(pnold,v2,estrains);
						get_estress(v2, estrains,estress);
						hstr = (estress[0]+estress[1])/2;

						if(hstr>0){sumfactor+=hstr/(k+hstr); }
					      	get_estrains(pnold,v3,estrains);
						get_estress(v3, estrains,estress);
						hstr = (estress[0]+estress[1])/2;

						if(hstr>0){sumfactor+=hstr/(k+hstr); }
					      	get_estrains(pnold,v4,estrains);
						get_estress(v4, estrains,estress);
						hstr = (estress[0]+estress[1])/2;

						if(hstr>0){sumfactor+=hstr/(k+hstr); }

						pn[n].fx=pn[n].fx*(1+par.LAMBDAFORCEFA*(sumfactor/4));
						pn[n].fy=pn[n].fy*(1+par.LAMBDAFORCEFA*(sumfactor/4));



				}



			}



		}

}

////////////////////////////////////////////////////////////////////////////////


BOOL CheckCellNodeConnection(VOX* pv, NOD* pn, int n1, int n2, int c){


	int ny1 = n1/NNX; int nx1=n1%NNX;
	int ny2 = n2/NNX; int nx2 = n2%NNX;
	int counter = 0;
	double diffx = nx2-nx1; double diffy = ny2-ny1;
	int vx1, vy1, vx, vy, v, n11, n22, nx, ny, n;
	double yideal, xideal;
	BOOL nodesconnected=TRUE;
	int bound = 0;
	int nocellvox = 0;
	double slope;
	int bn;
	int steps;
	int vxend1; int vyend1; int vxend2; int vyend2;

	if(!(diffx==0)){slope = diffy/diffx;}


//We begin with dividing into two cases:
//1: abs(diffx)>=abs(diffy) and diffy can not be zero
//2: abs(diffy)>abs(diffx) and the case diffy==0 is dealt with as well
//these two cases are equally dealth with


//If !diffx==0 and !diffy==0 we use Bresenham's line algorithm see: http://en.wikipedia.org/wiki/Bresenham's_line_algorithm
//we check if the pixels that denote the line from node1 to node2 are cell pixels, if not, we break


//if diffx==0 or diffy==0 we do something different
//in this case
//we check all nodes between node1 and node2, and say that of their neighbouring pixels there should be at least two cellular pixels
//so the number of neighbouring "other cellular" pixels, should be smaller than three

//we consider nodes on the boundary differently
//if there is no node in between node1 and node2, we check the neighbouring pixels of node1 and node2 themselves

	if(((abs(diffx)>=abs(diffy))&(!diffy==0))||(diffx==0))
	{
			n11=n1;n22=n2;
			if(diffx<0){n22=n1;n11=n2;} //swap the nodes, NEW
			ny1 = n11/NNX; nx1=n11%NNX;
			ny2 = n22/NNX; nx2 = n22%NNX;
			diffx = nx2-nx1; diffy = ny2-ny1;
			if(diffy>0&diffx>0){vx1=n11%NNX; vy1 = n11/NNX;}
			if(diffy<0&diffx>0){vx1=n11%NNX; vy1 = n11/NNX-1;}
			if(diffy>0&diffx<0){vx1=n11%NNX-1; vy1 = n11/NNX;}
			if(diffy<0&diffx<0){vx1=n11%NNX-1; vy1 = n11/NNX-1;}
			if(diffy>0&diffx>0){vyend2=n22/NNX-1; vyend1=0;}
			if(diffy<0&diffx>0){vyend1=n22/NNX; vyend2=par.NVY-1;}
			if(diffy>0&diffx<0){vyend1=0; vyend2=n22/NNX-1;}
			if(diffy<0&diffx<0){vyend1=n22/NNX; vyend2=par.NVY-1;}
			vx=vx1;
			vy=vy1;
			ny = ny1;
			nx=nx1;
			if(diffx==0){

			if(abs(diffy)==1)
			{		
					bound=1;
					if(ny2>ny1){n=nx1+ny1*NNX;}
					else{n=nx2+ny2*NNX;}
						vy = n/NNX; 
						vx = n%NNX;
				                v = vx + vy*(par.NVX);
						if(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1){bound=0;}
						if(!(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1)){if(!(pv[v].ctag==c+1)){nocellvox++;}}
						vy = n/NNX; 
						vx = n%NNX-1;
				                v = vx + vy*(par.NVX);
						if(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1){bound=0;}
						if(!(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1)){if(!(pv[v].ctag==c+1)){nocellvox++;}}
						if(nocellvox>bound){nodesconnected=FALSE;}
			}
			else
			{

				if(ny1<ny2){ny=ny1; bn = ny2;}else{ny=ny2; bn=ny1;}
				ny++;

				while(ny<=bn){
					nocellvox=0;
					bound = 1;
					n = nx + ny*(NNX);
					//this node should not be on a corner, so we check all four voxels
					//first lower right corner
					vy = n/NNX-1; //-1
					vx = n%NNX;
		                        v = vx + vy*(par.NVX);
						if(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1){bound=0;}
						if(!(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1)){if(!(pv[v].ctag==c+1)){nocellvox++;}}

					vy = n/NNX-1; //-1
					vx = n%NNX-1;
		                        v = vx + vy*(par.NVX);
						if(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1){bound=0;}
						if(!(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1)){if(!(pv[v].ctag==c+1)){nocellvox++;}}
					if(nocellvox>bound){nodesconnected=FALSE; }
					ny++;
				}}
			}
				else{
		                v = vx + vy*(par.NVX);
				if(!(pv[v].ctag==c+1)){nodesconnected=FALSE; }
				steps=1;
				while(steps<abs(diffx)){
					if(diffx>0){vx++;}
					if(diffx<0){vx--;}
					yideal= fabs(slope*steps);
					if(diffy>0){vy = vy1+round(yideal);}
					if(diffy<0){vy = vy1-round(yideal);}
					if(vy<vyend1){vy=vyend1;}
					if(vy>vyend2){vy=vyend2;}
		                        v = vx + vy*(par.NVX);
					if(!(pv[v].ctag==c+1)){ nodesconnected=FALSE; }
					steps++;
				}}

	}
	else
	{

			n11=n1;n22=n2;
			if(diffy<0){n22=n1;n11=n2;} //swap the nodes
			ny1 = n11/NNX; nx1=n11%NNX;
			ny2 = n22/NNX; nx2 = n22%NNX;
			diffx = nx2-nx1; diffy = ny2-ny1;
			if(diffx>0&diffy>0){vx1 = n11%NNX; vy1 = n11/NNX; }
			if(diffx<0&diffy>0){vx1 = n11%NNX-1; vy1 = n11/NNX;}
			if(diffx<0&diffy<0){vx1 = n11%NNX-1; vy1 = n11/NNX-1; }
			if(diffx>0&diffy<0){vx1 = n11%NNX; vy1 = n11/NNX-1; }
			if(diffx>0&diffy>0){vxend1=0; vxend2=n22%NNX-1;}
			if(diffx<0&diffy>0){vxend2=par.NVX-1; vxend1=n22%NNX;}
			if(diffx<0&diffy<0){vxend1=n22%NNX; vxend2=par.NVX-1;}
			if(diffx>0&diffy<0){ vxend2=n22%NNX-1; vxend1=0;}
			vx=vx1;
			vy=vy1;
			nx=nx1;
			ny=ny1;
			if(diffy==0){
			if(abs(diffx)==1)
			{	
					bound=1;
					if(nx2>nx1){n=nx1+ny1*NNX;}
					else{n=nx2+ny2*NNX;}
						vy = n/NNX;
						vx = n%NNX;
				                v = vx + vy*(par.NVX);
						if(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1){bound=0;}
						if(!(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1)){if(!(pv[v].ctag==c+1)){nocellvox++;}}
						vy = n/NNX-1;
						vx = n%NNX;
				                v = vx + vy*(par.NVX);
						if(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1){bound=0;}
						if(!(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1)){if(!(pv[v].ctag==c+1)){nocellvox++;}}
						if(nocellvox>bound){nodesconnected=FALSE;}
			}
			else
			{
				if(nx1<nx2){nx=nx1; bn = nx2;}else{nx=nx2; bn=nx1;}
				nx++;
				while(nx<=bn){
					bound = 1;
					n = nx + ny*(NNX);
					nocellvox=0;
					//this node should not be on a corner, so we check all four voxels
					vy = n/NNX-1;
					vx = n%NNX-1;
		                        v = vx + vy*(par.NVX);
						if(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1){bound=0;}
						if(!(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1)){if(!(pv[v].ctag==c+1)){nocellvox++;}}
					vy = n/NNX;
					vx = n%NNX-1;
		                        v = vx + vy*(par.NVX);
						if(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1){bound=0;}
						if(!(vx < 0 || vy <0 || vx > par.NVX-1 || vy >par.NVX-1)){if(!(pv[v].ctag==c+1)){nocellvox++;}}

					if(nocellvox>bound){nodesconnected=FALSE; }
					nx++;
				}}
			}
				else{
		                v = vx + vy*(par.NVX);
				if(!(pv[v].ctag==c+1)){ nodesconnected=FALSE; }
				steps=1;
				while(steps<abs(diffy)){
					if(diffy>0){vy++;}
					if(diffy<0){vy--;}
					xideal = fabs((1/slope)*steps);
					if(diffx>0){vx=vx1+round(xideal);}
					if(diffx<0){vx=vx1-round(xideal);}
					if(vx<vxend1){vx=vxend1;}
					if(vx>vxend2){vx=vxend2;}
		                        v = vx + vy*(par.NVX);
					if(!(pv[v].ctag==c+1)){ nodesconnected=FALSE; }
					steps++;
				}}

	}

	return nodesconnected;
}
