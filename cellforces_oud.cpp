// file: par.CELLFORCEs.c

#include <functions.h>
#include <math.h>
#include <boost/math/tools/minima.hpp>

#include <vector>
#include <cmath>
#include <iostream>

#include "/ufs/rens/opt/inmin/include/inmin_lm.h"

#include <boost/assign/list_of.hpp>


//res[1]=obs[1]-x[1]*sin(x[2]);
//res[2]=obs[2]-x[1]*(cos(x[2])+sin(x[2]))/sqrt(2);
//res[3]=obs[3]-x[1]*cos(x[2]);
//res[4]=obs[4]-x[1]*(cos(x[2])-sin(x[2]))/sqrt(2);




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

std::vector<double> minimise(ForceComponent &model,std::vector<double> &obs,std::vector<double> &lowerbounds,std::vector<double> &upperbounds, double* ic)
{

  inmin_lm_in in;
  memset(&in, 0, sizeof(inmin_lm_in));
  in.N=2;
  in.M=obs.size();
  in.f=LMHelper;

  in.ftol=1e-3;
  in.xtol=1e-3;
  in.gtol=1e-3;
  in.maxfev=1000;

  model.obs=obs;

	//cout << "obs[0] " << obs[0] << endl;
	//cout << "obs[1] " << obs[1] << endl;
	//cout << "obs[2] " << obs[2] << endl;
	//cout << "obs[3] " << obs[3] << endl;

	//cout << "ic[0] " << ic[0] << endl;
	//cout << "ic[1] " << ic[1] << endl;


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

  //double ic[] = {0.0005,3.1416/2};

  inmin_lm_run(&in, ic,&out,&model);

  //std::cout<<"The final result is:"<<res[0]<<","<<res[1]<<std::endl;

	return res;

};


void cell_forces_ham(VOX* pv, NOD* pn, int NRc, int* csize,int incr, double* dens,double* FA)
{


	int c;
	int n,nx,ny;
	int v,vy,vx,cnttag;
	int NRcelln,cellnodes[NN];
	int i,j;
	
	for(ny=1; ny<NNY-1; ny++)
   	for(nx=1; nx<NNX-1; nx++)
   	{
   		n = nx + ny*NNX;
		pn[n].fx = 0;
		pn[n].fy = 0;
	}


		for(c=0;c<NRc;c++)
		{

			// determine which nodes belong to cell c and are on the boundary
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



			// forces at node
			for(i=0;i<NRcelln;i++)
			{
				n = cellnodes[i];

				int ny=n/NNX;
				int nx=n%NNX;

				//lower right voxel is 4
				int vx4=n%NNX; int vy4 = n/NNX;int v4 = vx4 + vy4*(par.NVX);
				int vx2=vx4;int vy2=vy4-1;int v2 = vx2 + vy2*(par.NVX);
				int vx1=vx4-1;int vy1=vy4-1;int v1 = vx1 + vy1*(par.NVX);
				int vx3=vx4-1;int vy3=vy4;int v3 = vx3 + vy3*(par.NVX);

					int nr = 0;
					if(v1>=0&v1<NV){if(pv[v1].ctag>0){nr++;}}
					if(v2>=0&v2<NV){if(pv[v2].ctag>0){nr++;}}		
					if(v3>=0&v3<NV){if(pv[v3].ctag>0){nr++;}}		
					if(v4>=0&v4<NV){if(pv[v4].ctag>0){nr++;}}

				if(nr<4) //if boundary node
				{

					double f1,f2,f3,f4;
					int xt,xs,ttag,stag,pick;
				
					xt = v3;
					xs = v1;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=5;
					double dH13=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v1;
					xs = v3;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=1;
					double dH31=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v2;
					xs = v4;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=1;
					double dH42=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v4;
					xs = v2;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=5;
					double dH24=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v3;
					xs = v2;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=6;
					double dH23=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v2;
					xs = v3;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=2;
					double dH32=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v1;
					xs = v2;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=7;
					double dH21=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v2;
					xs = v1;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=3;
					double dH12=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v3;
					xs = v4;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=7;
					double dH43=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v4;
					xs = v3;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=3;
					double dH34=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v4;
					xs = v1;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=4;
					double dH14=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);

					xt = v1;
					xs = v4;
					ttag=pv[xt].ctag;
					stag=pv[xs].ctag;
					pick=0;
					double dH41=calcdH(pv, pn, csize, xt, xs, pick, ttag, stag, incr,dens,TRUE,FA);


					f1=(dH13-dH31)/4+(dH24-dH42)/4;
					f2=-(dH32-dH23)/(2*1.4142);
					f3=(dH21-dH12)/4+(dH43-dH34)/4;
					f4=-(dH14-dH41)/(2*1.4142);


  					std::vector<double> obs=boost::assign::list_of
					    (f1)(f2)(f3)(f4);

					  std::vector<double> x = boost::assign::list_of
					    (-10)(-8)(-5)(1)(4)(10)(40);	//extra variable we do not need now

				  	ForceComponent model(x);

					  std::vector<double> lower = boost::assign::list_of
					    (0)(0);
					  std::vector<double> upper = boost::assign::list_of
					    (100)(2*3.1416);

					std::vector<double> res;

					double *ic = new double[2];

					if(abs(f1)>=abs(f2) & abs(f1)>=abs(f2) & abs(f1)>=abs(f3)){ic[0]=f1;ic[1]=3.1416/2;}
					if(abs(f2)>=abs(f1) & abs(f2)>=abs(f3) & abs(f2)>=abs(f4)){ic[0]=f2;ic[1]=3.1416/4;}
					if(abs(f3)>=abs(f1) & abs(f3)>=abs(f2) & abs(f3)>=abs(f4)){ic[0]=f3;ic[1]=0;}
					if(abs(f4)>=abs(f1) & abs(f4)>=abs(f2) & abs(f4)>=abs(f3)){ic[0]=f4;ic[1]=2*3.1416-3.1416/4;}
		
					if(ic[0]<=0 & ic[0]==f1){ic[0]=-ic[0];ic[1]=ic[1]+3.1416;}
					if(ic[0]<=0 & ic[0]==f2){ic[0]=-ic[0];ic[1]=ic[1]+3.1416;}
					if(ic[0]<=0 & ic[0]==f3){ic[0]=-ic[0];ic[1]=ic[1]+3.1416;}
					if(ic[0]<=0 & ic[0]==f4){ic[0]=-ic[0];ic[1]=ic[1]-3.1416;}


					  // Run the minimiser

			
					  //res=minimise(model,obs,lower,upper,ic);
					res=boost::assign::list_of(ic[0])(ic[1]);


					double V=res[0]; double alpha=res[1];
					pn[n].fx=V*cos(-alpha);
					pn[n].fy=V*sin(-alpha);
					delete [] ic;


					/* if(n==25){
					cout << "v1 " << v1 << endl;
					cout << "v2 " << v2 << endl;
					cout << "v3 " << v3 << endl;
					cout << "v4 " << v4 << endl;
					cout << "dH13 " << dH13 << endl;
					cout << "dH31 " << dH31 << endl;
					cout << "dH24 " << dH24 << endl;
					cout << "dH42 " << dH42 << endl;
					cout << "dH32 " << dH32 << endl;
					cout << "dH23 " << dH23 << endl;
					cout << "dH21 " << dH21 << endl;
					cout << "dH12 " << dH12 << endl;
					cout << "dH43 " << dH43 << endl;
					cout << "dH34 " << dH34 << endl;
					cout << "dH14 " << dH14 << endl;
					cout << "dH41 " << dH41 << endl;
					cout << "f1 " << f1 << endl;
					cout << "f2 " << f2 << endl;
					cout << "f3 " << f3 << endl;
					cout << "f4 " << f4 << endl;
					}  */


	

if(n==54){
//cout << "n " << n << endl;				
//cout << "V " << V << endl;
//cout << "f4 " << f4 << endl;
//cout << "alpha " << alpha << endl;
}


//pn[n].fx=pn[n].fx*(1+FA[n]/par.CAPACITYFA)/2;
//pn[n].fy=pn[n].fy*(1+FA[n]/par.CAPACITYFA)/2;

	//int xty = xs/(par.NVX); int xtx = xs%(par.NVX);
	//int n00 = (xtx  ) + (xty  )*NNX;
	//int n10 = (xtx+1) + (xty  )*NNX;
	//int n11 = (xtx+1) + (xty+1)*NNX;
	//int n01 = (xtx  ) + (xty+1)*NNX;
	//double FAstick = FA[n00]+FA[n10]+FA[n11]+FA[n01];
	//w=FAstick/(4*par.CAPACITYFA);

	

//cout << "V " << V << endl;




				}


			}
		}

}
