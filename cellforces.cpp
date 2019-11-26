// file: par.CELLFORCEs.c

#include <functions.h>
#include <math.h>
//#include <boost/math/tools/minima.hpp>

#include <vector>
#include <cmath>
#include <iostream>

//#include "/ufs/rens/opt/inmin/include/inmin_lm.h"

//#include <boost/assign/list_of.hpp>


//res[1]=obs[1]-x[1]*sin(x[2]);
//res[2]=obs[2]-x[1]*(cos(x[2])+sin(x[2]))/sqrt(2);
//res[3]=obs[3]-x[1]*cos(x[2]);
//res[4]=obs[4]-x[1]*(cos(x[2])-sin(x[2]))/sqrt(2);



/*
/// This function interfaces between the imnin library and simple C++
/// code. It is hard coded for this particular problem.
void LMHelper(const double *x,double *res,void *data)
{
  ForceComponent *pmodel=reinterpret_cast<ForceComponent*>(data);
  pmodel->a=x[0];
  pmodel->b=x[1];

  //std::cout<<">> Interim point: "<<x[0]<<","<<x[1]<<std::endl;


  std::vector<double> scratch;

  (*pmodel)(scratch);

  for (size_t i=0; i<pmodel->obs.size(); ++i)
  {
    res[i]=scratch[i]-pmodel->obs[i];
    //cout << "i " << i << " " << "scratch[i] " << scratch[i] << endl;
    //cout << "i " << i << " " << "pmodel->obs[i] " << pmodel->obs[i] << endl;
    //cout << "i " << i << " " << "res[i] " << res[i] << endl;
  }
}
*/

/*
std::vector<double> minimise(ForceComponent &model,std::vector<double> &obs,std::vector<double> &lowerbounds,std::vector<double> &upperbounds, double* ic)
{

  inmin_lm_in in;
  memset(&in, 0, sizeof(inmin_lm_in));
  in.N=2;
  in.M=obs.size();
  in.f=LMHelper;

  in.ftol=1e-8;
  in.xtol=1e-8;
  in.gtol=1e-8;
  in.maxfev=100000;

  model.obs=obs;

	//cout << "obs[0] " << obs[0] << endl;
	//cout << "obs[1] " << obs[1] << endl;
	//cout << "obs[2] " << obs[2] << endl;
	//cout << "obs[3] " << obs[3] << endl;




  double box[4];
  for (size_t i=0; i<2; ++i)
  {
     box[i*2]=lowerbounds[i];
     box[i*2+1]=upperbounds[i];
  }
  in.box=box;

  std::vector<double> res(2);
  inmin_lm_out out;
  out.x=&res[0];
  out.covar=NULL;


  inmin_lm_run(&in, ic,&out,&model);

  //std::cout<<"The initial is:"<<ic[0]<<" , "<<ic[1]<<std::endl;
  //std::cout<<"The final result is:"<<res[0]<<" , "<<res[1]<<std::endl;

	return res;

};
*/

/*
void cell_forces_ham(VOX* pv, NOD* pn, int NRc, int* csize,int incr, double* dens,int* celltypes)
{


	int c;
	int n,nx,ny;
	int v,vy,vx;
	int i,j;
	int count=0;

			for(ny=1; ny<NNY-1; ny++)
	   		for(nx=1; nx<NNX-1; nx++)
	   		{
	   			n = nx + ny*NNX;

				//lower right voxel is 4
				int vx4=n%NNX; int vy4 = n/NNX;int v4 = vx4 + vy4*(par.NVX);
				int vx2=vx4;int vy2=vy4-1;int v2 = vx2 + vy2*(par.NVX);
				int vx1=vx4-1;int vy1=vy4-1;int v1 = vx1 + vy1*(par.NVX);
				int vx3=vx4-1;int vy3=vy4;int v3 = vx3 + vy3*(par.NVX);

					int pixels[4];



					//determine if this node is on the boundary of a cell
					bool cellboundary=FALSE;
					int nrp=0;
					if(v1>=0&v1<NV){pixels[nrp]=v1;nrp++;}
					if(v2>=0&v2<NV){pixels[nrp]=v2;nrp++;}
					if(v3>=0&v3<NV){pixels[nrp]=v3;nrp++;}
					if(v4>=0&v4<NV){pixels[nrp]=v4;nrp++;}
					for(int pix1=0;pix1<nrp;pix1++)
					{
						for(int pix2=0;pix2<nrp;pix2++)
						{
								if(pv[pixels[pix1]].ctag!=pv[pixels[pix2]].ctag){cellboundary=TRUE; break;}
						}	
					}


				if(cellboundary) //if boundary node then we calculate a force
				{count++;



					double f1,f2,f3,f4;
					int xt,xs,ttag,stag,pick;
					xt = v3;
					xs = v1;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=5;
					double dH13=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc, celltypes);



					xt = v1;
					xs = v3;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=1;
					double dH31=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc, celltypes);

					xt = v2;
					xs = v4;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=1;
					double dH42=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc, celltypes);

					xt = v4;
					xs = v2;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=5;
					double dH24=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc, celltypes);

					xt = v3;
					xs = v2;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=6;
					double dH23=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc, celltypes);



					xt = v2;
					xs = v3;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=2;
					double dH32=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc, celltypes);

					xt = v1;
					xs = v2;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=7;
					double dH21=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc, celltypes);

					xt = v2;
					xs = v1;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=3;
					double dH12=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc,celltypes);

					xt = v3;
					xs = v4;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=7;
					double dH43=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc, celltypes);

					xt = v4;
					xs = v3;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=3;
					double dH34=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc, celltypes);

					xt = v4;
					xs = v1;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=4;
	
					double dH14=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc,celltypes);
			



					xt = v1;
					xs = v4;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=0;
					double dH41=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,NRc,celltypes);




					f1=-(dH13-dH31)/4-(dH24-dH42)/4;
					f2=(dH32-dH23)/(2*1.4142);
					f3=-(dH21-dH12)/4-(dH43-dH34)/4;
					f4=(dH14-dH41)/(2*1.4142);




  					std::vector<double> obs=boost::assign::list_of
					    (f1)(f2)(f3)(f4);

					  std::vector<double> x = boost::assign::list_of
					    (-10)(-8)(-5)(1)(4)(10)(40);	//extra variable we do not need now

				  	ForceComponent model(x);

					  std::vector<double> lower = boost::assign::list_of
					    (0)(-10000*3.1416);
					  std::vector<double> upper = boost::assign::list_of
					    (1000000)(10000*3.1416);

					std::vector<double> res;

					double *ic = new double[2];
					ic[0]=2; ic[1]=0;

					double af1=3.1416/2;
					double af2=3.1416/4;
					double af3=0;
					double af4=2*3.1416-3.1416/4;

					if(abs(f1)>=abs(f2) & abs(f1)>=abs(f3) & abs(f1)>=abs(f4))
					{
						ic[0]=f1;ic[1]=af1;
						if(f1==f2){ic[1]=(af1+af2)/2;}
						if(f1==f3){ic[1]=(af1+af3)/2;}
						if(f1==f4){ic[1]=(af1+af4)/2;}
					}
					if(abs(f2)>=abs(f1) & abs(f2)>=abs(f3) & abs(f2)>=abs(f4))
					{
						ic[0]=f2;ic[1]=af2;
						if(f2==f1){ic[1]=(af2+af1)/2;}
						if(f2==f3){ic[1]=(af2+af3)/2;}
						if(f2==f4){ic[1]=(af2+af4)/2;}
					}
					if(abs(f3)>=abs(f1) & abs(f3)>=abs(f2) & abs(f3)>=abs(f4))
					{
						ic[0]=f3;ic[1]=af3;
						if(f3==f1){ic[1]=(af3+af1)/2;}
						if(f3==f2){ic[1]=(af3+af2)/2;}
						if(f3==f4){ic[1]=(af3+af4)/2;}
					}
					if(abs(f4)>=abs(f1) & abs(f4)>=abs(f2) & abs(f4)>=abs(f3))
					{
						ic[0]=f4;ic[1]=af4;
						if(f4==f1){ic[1]=(af4+af1)/2;}
						if(f4==f2){ic[1]=(af4+af2)/2;}
						if(f4==f3){ic[1]=(af4+af3)/2;}
					}

		
					if(ic[0]<0){ic[0]=-ic[0];ic[1]=ic[1]+3.1416;}
					if(ic[1]>2*3.1416){ic[1]=ic[1]-2*3.1416;}



					

					  // Run the minimiser

					res=minimise(model,obs,lower,upper,ic);
					//res=boost::assign::list_of(ic[0])(ic[1]);

					double V=res[0]; double alpha=res[1];
					//test
					//V=ic[0];alpha=ic[1];





					//minus because f=-dH
					//extra minus for fy, because of plot coordinates
					pn[n].fx=-V*cos(alpha);
					pn[n].fy=V*sin(alpha);



//cout << endl;
//cout << "n " << n << endl;
//cout << "ic[0] " << ic[0] << endl;
//cout << "V " << V << endl;
//cout << "ic[1] " << ic[1] << endl;
//cout << "alpha " << alpha << endl;
//cout << "pn[n].fx " << pn[n].fx << endl;
//cout << "pn[n].fy " << pn[n].fy << endl;


		
					delete [] ic;

pn[n].fx=pn[n].fx/par.THICKNESS;
pn[n].fy=pn[n].fy/par.THICKNESS;




				}

		}
cout << "count " << count << endl;


}
*/
////////////////////////////////////////////////////////////////////////////////
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

		//calc inertia
	double Ixx=0;double Iyy=0;
		double cmx = double(csumx[c])/csize[c];
		double cmy = double(csumy[c])/csize[c];
		for(int v = 0;v<NV;v++)
		{

			if(pv[v].ctag == c+1)
			{
				int vy = v/(par.NVX); int vx = v%(par.NVX);
				Ixx += (vx-cmx)*(vx-cmx);
				Iyy += (vy-cmy)*(vy-cmy);

			}

		}

if(Ixx==0){Ixx=0.5;}
if(Iyy==0){Iyy=0.5;}

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

//cout << "NRcelln " << NRcelln << endl;

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

//cout << "n " << n << " , n2 " << n2 << " , connected " << CheckCellNodeConnection(pv,pn,n,n2,c) << endl;

					if(par.NODECONNECTION)
					{

						if(CheckCellNodeConnection(pv,pn,n,n2,c))
						{


							count++;

							//forcex = par.LRTENSION*dnx/double(csize[c]*csize[c])/par.THICKNESS;
							//forcey = par.LRTENSION*dny/double(csize[c]*csize[c])/par.THICKNESS;
							//forcex = 0.5*par.LRTENSION*dnx/(Ixx+Iyy)/par.THICKNESS;
							//forcey =  0.5*par.LRTENSION*dny/(Ixx+Iyy)/par.THICKNESS;

							//forcex = par.LRTENSION*dnx/par.THICKNESS;
							//forcey =  par.LRTENSION*dny/par.THICKNESS;

							if(par.CLASSICCPM)
							{
								forcex = par.LRTENSION*dnx/double(csize[c])/par.THICKNESS;
								forcey = par.LRTENSION*dny/double(csize[c])/par.THICKNESS;

			

								//forcex = par.LRTENSION*dnx/par.THICKNESS;
								//forcey = par.LRTENSION*dny/par.THICKNESS;
							}

							if(!par.CLASSICCPM)
							{
								forcex = par.LRTENSION*dnx/double(csize[c]*csize[c])/par.THICKNESS;
								forcey = par.LRTENSION*dny/double(csize[c]*csize[c])/par.THICKNESS;

								//forcex = par.LRTENSION*dnx/par.THICKNESS;
								//forcey = par.LRTENSION*dny/par.THICKNESS;
							}



						pn[n].fx += forcex;
						pn[n].fy += forcey;



						}


					}
					else
					{


							//forcex = par.LRTENSION*dnx/double(csize[c]*csize[c])/par.THICKNESS;
							//forcey = par.LRTENSION*dny/double(csize[c]*csize[c])/par.THICKNESS;
							//forcex =  0.5*par.LRTENSION*dnx/(Ixx+Iyy)/par.THICKNESS;
							//forcey =  0.5*par.LRTENSION*dny/(Ixx+Iyy)/par.THICKNESS;


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
//cout << "count " << count << endl;
//cout << "csize " << csize[c] << endl;
			//pn[n].fx=pn[n].fx/(count*count); //try to average
			//pn[n].fy=pn[n].fy/(count*count);


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
				/*if(v1==18504-1*10| v1 ==18504+9*200+5*200 | v1==18504 |v1 ==18504+9*200+5*200-10 )
				{
				cout << "hstr " << hstr << endl;
				}*/
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


	if(sumfactor>0){
	//cout << "n " << n << endl;
	//cout << "vx1 " << vx1 << endl;
	//cout << "vy1 " << vy1 << endl;
	//cout << "par.LAMBDAFORCEFA*sumfactor/4 " << par.LAMBDAFORCEFA*sumfactor/4 << endl;
	}

				}


//cout << "pn[n].fx " << pn[n].fx << endl;
//cout << "pn[n].fy " << pn[n].fy << endl;

			}




		}

	//for(int n=1;n<NV;n++)
	//{
	//	if(pn[n].fx*pn[n].fx+pn[n].fy*pn[n].fy<2e-5){pn[n].fx=0;pn[n].fy=0;}
	//}

	




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
