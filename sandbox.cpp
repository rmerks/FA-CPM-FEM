#include "myparameters.h"
extern "C" {
#include <stdlib.h>
}
#include "functions.h"
#include "plotcpm.h"
#include <QApplication>
#include <QMainWindow>
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <string>
#include <fstream>
#include <sstream>


Parameter par;


int main(int argc, char *argv[])
{
       cerr << "Hello 1\n";
	int nrrdof, d;
	int *dofpos;
	double *f, *u, *uold;
	double **klocal;
	int *kcol;
	double *kval;
	NOD *pn, *pnold,*pnstall; 	
	VOX *pv;
	int NRc,c,v;
	int *csize,*csumx,*csumy;
	int *cper;
	double *clength, *ecc, *cangle;
	int incr, startincr;
	int *cmxi, *cmyi;
	double* sqdis;
	int contact = 0;
	int* celltypes;



	
	
	//FOR SAVING FILES			
	stringstream currentdir;
	string thisdir;
	thisdir=argv[1];
	currentdir << thisdir;
	string scurrentdir = currentdir.str();
	const char * ccurrentdir = scurrentdir.c_str();

	stringstream lstr;
	string lstr1 = scurrentdir;
	string lstr2;
	lstr2 = "/length.txt";
	lstr<<lstr1<<lstr2;
	string ls = lstr.str();
	const char * lfcstr = ls.c_str();
	FILE *ofpl;
	ofpl = fopen(lfcstr,"w");

	stringstream pastr;
	string pastr1 = scurrentdir;
	string pastr2;
	pastr2 = "/ratiopa.txt";
	pastr<<pastr1<<pastr2;
	string pas = pastr.str();
	const char * pafcstr = pas.c_str();
	FILE *ofppa;
	ofppa = fopen(pafcstr,"w");

	stringstream astr;
	string astr1 = scurrentdir;
	string astr2;
	astr2 = "/area.txt";
	astr<<astr1<<astr2;
	string as = astr.str();
	const char * afcstr = as.c_str();
	FILE *ofpa;
	ofpa = fopen(afcstr,"w");

	stringstream sdstr;
	string sdstr1 = scurrentdir;
	string sdstr2;
	sdstr2 = "/sqdis.txt";
	sdstr<<sdstr1<<sdstr2;
	string sds = sdstr.str();
	const char * sdfcstr = sds.c_str();
	FILE *ofpsd;
	ofpsd = fopen(sdfcstr,"w");

	stringstream estr;
	string estr1 = scurrentdir;
	string estr2;
	estr2 = "/ecc.txt";
	estr<<estr1<<estr2;
	string es = estr.str();
	const char * efcstr = es.c_str();
	FILE *ofpe;
	ofpe = fopen(efcstr,"w");

	stringstream anstr;
	string anstr1 = scurrentdir;
	string anstr2;
	anstr2 = "/angle.txt";
	anstr<<anstr1<<anstr2;
	string ans = anstr.str();
	const char * anfcstr = ans.c_str();
	FILE *ofpan;
	ofpan = fopen(anfcstr,"w");

	stringstream tccstr;
	string tccstr1 = scurrentdir;
	string tccstr2;
	tccstr2 = "/twocellcontact.txt";
	tccstr<<tccstr1<<tccstr2;
	string tccs = tccstr.str();
	const char * tccfcstr = tccs.c_str();
	FILE *ofptcc;
	ofptcc = fopen(tccfcstr,"w");

	stringstream ssreadcells;
	string sreadcells1 = scurrentdir;
	string sreadcells2;
	sreadcells2 = "/ctags0.out";
	ssreadcells<<sreadcells1<<sreadcells2;
	string sreadcells = ssreadcells.str();

	stringstream sssigmas;
	string ssigmas1 = scurrentdir;
	string ssigmas2;
	ssigmas2 = "/sigmas/";
	sssigmas<<ssigmas1<<ssigmas2;
	string ssigmas = sssigmas.str();
	const char * sigmasdir = ssigmas.c_str();
	mkdir(sigmasdir,0777);

	stringstream sstotshape;
	string stotshape1 = scurrentdir;
	string stotshape2;
	stotshape2 = "/totshape/";
	sstotshape<<stotshape1<<stotshape2;
	string stotshape = sstotshape.str();
	const char * totshapedir = stotshape.c_str();
	mkdir(totshapedir,0777);  
        
        stringstream ssfa;	
	string sfa1 = scurrentdir;
        string sfa2;
        sfa2 = "/fa/";
        ssfa<<sfa1<<sfa2;
        string sfa = ssfa.str();
        const char * fadir = sfa.c_str();
        mkdir(fadir,0777);




    //to read parameters
    if (argc <= 1)
    {
        cout << "Usage: " << argv[0] << " <ParametersDirectory>" << endl;
        exit(1);
    }

   cerr << "Hello 1\n";

	stringstream str;
	string str1, str2;
	str1=argv[1];
	str2="/parameters.txt";
	str<<str1<<str2;
	string s = str.str();
	const char * cstr = s.c_str();
	//First argument is directory
	par.Read(cstr);


    /// GRAPHICS //
    QApplication a(argc, argv);
    setlocale(LC_NUMERIC,"C");
    if(!par.WIDTHCOLORBAR){par.WIDTHCOLORBAR=par.NVX*par.PIXPERVOX/10;}
    if(!par.COLORBAR){par.WIDTHCOLORBAR=0;}

	/// INITIALIZE ///
    if ((par.SEED)) { srand((par.SEED)); 
    }


    mt_init();
    pv = init_voxels();
    pn = init_nodes();
    pnold = init_nodes();
    pnstall = init_nodes();



    startincr = 0;
    if(startincr==0) {
    
	if(par.CELLCOL)
	{
	        NRc = init_cells(pv);
		celltypes = new int[NRc];
		for(int i=0;i<NRc;i++)
		{
			celltypes[i]=1;
			double r01 = rand()/(double)RAND_MAX;
			if(r01<1/double(2)) //i%2==0
				{

					if(par.TWOCELLTYPES){celltypes[i]=2;}
				}

		}
	}


        if(par.TWOCELL)
	{
      	// two cells in the middle
        int vx1=(par.NVX)/2-par.DISTWOCELLS, vx2=(par.NVX)/2+par.DISTWOCELLS;
        int vy=par.NVY/2;
        
        v = vx1 + vy*(par.NVX);
        pv[v].ctag=1;
        v = vx2 + vy*(par.NVX);
        pv[v].ctag=2;
        NRc=2; 
	celltypes = new int[NRc];
	celltypes[0]=1;celltypes[1]=1;if(par.TWOCELLTYPES){celltypes[1]=2;}
	}
      

	if(par.ONECELL)
	{
	 // one cell in the middle
		int vx=(par.NVX)/2;
		int vy=par.NVY/2;

		v = vx + vy*(par.NVX);
		pv[v].ctag=1;
		NRc=1; 
		celltypes = new int[NRc];
		celltypes[0]=1;
	}
	

	//USE /and* ....  *and/  to make a big section into commentary

        //write_cells(pv,0);
    }
	if(par.READCELLS) 
	{
	NRc = read_cells(pv,startincr,sreadcells);
		celltypes = new int[NRc];
		if(NRc==1){celltypes[0]=1;}
		if(NRc>=2){
				for(int i=0;i<NRc;i++)
				{
					celltypes[i]=1;
					double r01 = rand()/(double)RAND_MAX;
					if(r01<1/double(2)) //i%2==0
					{
						if(par.TWOCELLTYPES){celltypes[i]=2;}
					}

				}
		}

	}




	if(NRc==2){par.TWOCELL=true;}
	cout << "NRc" << NRc << endl;



	//write NRC
	stringstream nrcstr;
	string nrcstr1 = scurrentdir;
	string nrcstr2;
	nrcstr2 = "/nrc.txt";
	nrcstr<<nrcstr1<<nrcstr2;
	string nrcs = nrcstr.str();
	const char * nrcfcstr = nrcs.c_str();
	FILE *ofpnrc;
	ofpnrc = fopen(nrcfcstr,"w");
	fprintf(ofpnrc ,"%d ",NRc);
	fclose(ofpnrc);


	if(par.WRATIOPA){cper = new int[NRc];
		for(c=0;c<NRc;c++) {cper[c]=0;}}
	if(par.WLENGTH | par.WECC | par.WANGLE)
	{
		clength = new double[NRc];
		for(c=0;c<NRc;c++) {clength[c]=0;}
		ecc = new double[NRc];
		for(c=0;c<NRc;c++) {ecc[c]=0;}
		if(!par.WTWOCELLANGLECM){cangle = new double[NRc]; for(c=0;c<NRc;c++) {cangle[c]=0;}}
		if(par.TWOCELL & par.WTWOCELLANGLECM){cangle = new double[3];for(c=0;c<3;c++) {cangle[c]=0;}}
		
	}

	if(par.WSQDIS){sqdis = new double[NRc];
		for(c=0;c<NRc;c++) {sqdis[c]=0;}}
	csize = new int[NRc];
		for(c=0;c<NRc;c++) {csize[c]=0;}
	csumx = new int[NRc];
		for(c=0;c<NRc;c++) {csumx[c]=0;}
	csumy = new int[NRc];
		for(c=0;c<NRc;c++) {csumy[c]=0;}
	cmxi = new int[NRc];
	cmyi = new int[NRc];

	int* totshape = new int[(par.NVX*2+1)*(par.NVY*2+1)];
	for(int vx=0; vx<2*(par.NVX)+1; vx++) {
        for (int vy=0; vy<2*par.NVY+1; vy++) {
            int v = vx + vy * (2*(par.NVX)+1);totshape[v]=0; }}


    
    // set initial volumes and center of mass tags
    for (v=0;v<NV;v++) {
            int y = v/(par.NVX); int x = v%(par.NVX);
            if (pv[v].ctag) {
                csize[pv[v].ctag-1]++;
                csumx[pv[v].ctag-1]+=x;
                csumy[pv[v].ctag-1]+=y;
	
            }
    }



	//initial center of mass
    for(c=0;c<NRc;c++)
	{
		cmxi[c]=csumx[c]/csize[c];
		cmyi[c]=csumy[c]/csize[c];
	}

    
    // set initial volumes
    //for(v=0;v<NV;v++) {if(pv[v].ctag) {csize[pv[v].ctag-1]++;}}
    


	set_forces(pn,0);
	set_forces(pnstall,0);
	set_forces(pnold,0);
        

	//dont set restrictions if static stress 
	if(!par.GLOBALSTRAIN){set_restrictions(pn);set_restrictions(pnold);set_restrictions(pnstall);}
	// local K matrix
	klocal = set_klocal();

	// global K matrix:
	kcol = new int[10*NDOF];
	kval = new double[10*NDOF];
	assembly(kcol,kval,klocal,pv);
	dofpos = new int[NDOF];
	nrrdof = arrange_dofpos(dofpos,pn);


	reduce_K(kcol,kval,dofpos);


	double* Kms;
	Kms= new double[NDOF];
	for(int n=0;n<NDOF;n++){Kms[n]=0;}
	kval_to_nodes(Kms, kval,pn);

	
	uold = new double[nrrdof];
	for(int nn=0;nn<NN;nn++){uold[nn]=0;}

	//NEW: FA now only one dimensional (circle), later maybe growth in direction (ellipse)
	double *FA;
	FA = new double[NV];
	for(int n=0;n<NV;n++){FA[n]=0;} 
	//for(int n=0;n<NV;n++)
	//{
	//	if(pv[n].ctag){FA[n]=par.BASEFA;}
	//}

	double* sumFA2;
	sumFA2 = new double[NRc];


	double *sumFA;
	sumFA = new double[NRc];
	for(int c=0;c<NRc;c++){sumFA[c]=0;}
	for(int n=0;n<NV;n++)
	{
		if(FA[n]>0)
		{
			sumFA[pv[n].ctag-1]=sumFA[pv[n].ctag-1]+FA[n];
		}
	}


	for(int vx=0; vx<par.NVX; vx++) {
        for (int vy=0; vy<par.NVY; vy++) {
		
	    int v = vx+vy*par.NVX;
	    if(pv[v].ctag==1)
		{
		    int cmx = round(csumx[0]/csize[0]); 
		    int cmy = round(csumy[0]/csize[0]); 
		    int disx = cmx-vx; 
		    int disy =cmx-vy; 
		    int x = par.NVX+1+disx; 
		    int y = par.NVY+1+disy; 
		    int v2 = x+y*(2*(par.NVX)+1); 
		    totshape[v2]=totshape[v2]+1;
		}
	

	}}



    	FILE *of=fopen("com.dat","w");




	double *Nt = new double[NRc]; 
	/// START SIMULATION ///
	for(incr=startincr; incr<(par.NRINC); incr++)
	{
		printf("\nSTART INCREMENT %d",incr);

	set_forces(pn,incr); //set all to zero 
	set_forces(pnstall,incr);
	if(par.LEMMONROMER){cell_forces(pv, pnstall, csize, NRc,csumx,csumy,FA, pnold);}


	for(int c=0;c<NRc;c++){sumFA2[c]=0;}
	for(int n=0;n<NV;n++)
	{
		if(FA[n]>0)
		{
			sumFA2[pv[n].ctag-1]=sumFA2[pv[n].ctag-1]+FA[n];
		}
	}
	cout << endl;
	cout << "sumFA2 " << sumFA2[0] << endl;
	cout << "sumFA " << sumFA[0] << endl;
	cout << "faode " << endl;


		double dt = par.PDEdt;
	//FORCE BUILD UP AND FOCAL ADHESION GROWTH


	for(int tt=0;tt<par.PDEREPEAT;tt++)
	{


		for(int n=0;n<NN;n++)
		{
		if(pnstall[n].fx==0&pnstall[n].fy==0){pnold[n].fx=0;pnold[n].fy=0;} //if here no cell at the moment, build up from zero next time there is a cell

		//surrounding voxels of this node
		int vx1=n%NNX; int vy1 = n/NNX; int v1 = vx1 + vy1*(par.NVX);
		int vx2=vx1;int vy2=vy1-1;int v2 = vx2 + vy2*(par.NVX);
		int vx3=vx1-1;int vy3=vy1-1;int v3 = vx3 + vy3*(par.NVX);
		int vx4=vx1-1;int vy4=vy1;int v4 = vx4 + vy4*(par.NVX);

		double nrFA=0; 
		if(vx1>=0 && vx1<par.NVX && vy1>=0 && vy1<par.NVY){if(FA[v1]>par.BASEFA){nrFA++;}}
		if(vx2>=0 && vx2<par.NVX && vy2>=0 && vy2<par.NVY){if(FA[v2]>par.BASEFA){nrFA++;}}
		if(vx3>=0 && vx3<par.NVX && vy3>=0 && vy3<par.NVY){if(FA[v3]>par.BASEFA){nrFA++;}}
		if(vx4>=0 && vx4<par.NVX && vy4>=0 && vy4<par.NVY){if(FA[v4]>par.BASEFA){nrFA++;}}
		double factor = 1;
		if(nrFA==0){factor=0;}//if no focal adhesion surrounding the node, then start with no force
		factor = nrFA/4;

		if(tt==0){pnold[n].fx=pnold[n].fx*factor;pnold[n].fy=pnold[n].fy*factor;}


		if(!(pnstall[n].fx==0)&!(pnstall[n].fy==0))//there is a cell here
		{





			double Km=Kms[2*n]; //the x

		
			double forcemag=sqrt(pnstall[n].fx*pnstall[n].fx+pnstall[n].fy*pnstall[n].fy);
			double tk=forcemag/(par.VISC*Km);

			pn[n].fx=pnstall[n].fx+(pnold[n].fx-pnstall[n].fx)*exp(-dt*tt/tk);
			pn[n].fy=pnstall[n].fy+(pnold[n].fy-pnstall[n].fy)*exp(-dt*tt/tk);
			if(tt==par.PDEREPEAT-1){pnold[n].fx=pn[n].fx;pnold[n].fy=pn[n].fy;}
			
			}
		}
		for(int c=0;c<NRc;c++){Nt[c]=par.CAPACITYFA-sumFA[c];}//20
		for(int n=1;n<NV;n++)
		{
			if(pv[n].ctag>0)
			{

				double etractionstress[2];
				int cellnr=pv[n].ctag;
				double tension=0;
				get_etractionstress(pn,n, etractionstress);
				double strx=etractionstress[0]; double stry=etractionstress[1];
				tension=sqrt(strx*strx + stry*stry)*par.VOXSIZE*par.VOXSIZE;
				double g =0;
				g=par.GROWTHFA;
				
				double sliptension = par.SLIPTENSION;  
				double catchtension = par.CATCHTENSION;
				double tension2 = tension*1e12; //because picoNewton

				double rate = exp((tension2/FA[n]-sliptension))+exp(-(tension2/FA[n]-catchtension));
				rate=rate*par.DECAYFA;

				double Na = Nt[cellnr-1]; if(Na<0){Na=0;}

				double maxfa = par.MAXFAPP;

				double s=FA[n]-dt*rate*FA[n]+dt*g*Na*(1-par.LOGISTICPAR*FA[n]/maxfa); //logistic growth 


				if(FA[n]>=maxfa){s=FA[n]-dt*rate*FA[n];} 
				if(s<0){s=0;}
				
				if(FA[n]==0){s=0;} //rate is infinite
				double diff = FA[n]-s;
			
				FA[n]=s;

				sumFA[cellnr-1]=sumFA[cellnr-1]-diff;

			}
		}
	}
	


	for(int c=0;c<NRc;c++){sumFA2[c]=0;}
	for(int n=0;n<NV;n++)
	{
		if(FA[n]>0)
		{
			sumFA2[pv[n].ctag-1]=sumFA2[pv[n].ctag-1]+FA[n];
		}
	}
	cout << "sumFA2 " << sumFA2[0] << endl;
	cout << "sumFA " << sumFA[0] << endl;	
	cout << endl;

	int totfa=0;
	for(int n=0;n<NV;n++) 
        {
		if(FA[n]>0 & pv[n].ctag>0){totfa++; }
	}
	cout << "totfa " << totfa << endl;
	cout << "csize " << csize[0] << endl;




	// FEA part // parts of this can go out the loop depending on what changes
	f = new double[nrrdof];
	place_node_forces_in_f(pn,f);
	//place_node_forces_in_f(pnold,f);
	u = new double[nrrdof];

	set_disp_of_prev_incr(pn,u);
	solvePCG(kcol,kval,u,f,nrrdof); 
	disp_to_nodes(pn,u);
	disp_to_nodes(pnold,u);

		//free(u); free(f);
		delete [] u; delete [] f;

	stringstream fstr;
	string fstr2;
	fstr2="cpmfem%05d.png";
	fstr<<str1<<fstr2;
	string fs = fstr.str();
	const char * fcstr = fs.c_str();

	//calctension(pn,klocal);
        if (!(incr%(par.STRIDE))) {
    PlotCPM plot((par.NVX)*(par.PIXPERVOX)+3*par.WIDTHCOLORBAR,par.NVY*(par.PIXPERVOX));

             char fname[200];
             sprintf(fname,fcstr,incr);

		if(par.PRINCFIELD | par.STRAINFIELD){
		plot.CalculateStrainField(pn);}

		if(par.PRINCFIELD | par.STRESSFIELD){
		plot.CalculateStressField(pn);}



	        if(!par.TESTNEWFAPLOT){plot.PlotHueStressField(par.STRESSFIELD);}
	        if(par.TESTNEWFAPLOT){
plot.PlotHueStressField_CTB(par.STRESSFIELD);
// plot.PlotHueStressField(par.STRESSFIELD);
}
	        //plot.PlotHueStrainField(par.STRAINFIELD);


                if(par.PRINCFIELD){ 		
	     	plot.PlotPrincipleStressField(pv);}

		plot.Plot(pv,celltypes);

		

		if(par.FAFIELD)
		{

		if(par.TESTNEWFAPLOT==false) {plot.PlotNodalFA(pn,FA);}
		else
			plot.PlotNodalFA2(pv,FA);
		}

		if(par.STRESSTENSOR){plot.PlotStressTensor(pn,pv);}
		if(par.DEFORMFIELD){plot.PlotNodalDeform(pn);}



		if(par.FORCEFIELD){
		plot.CalculateForceField(pn); 
		plot.PlotNodalForces(pn,pv);}

	
		if(par.COLORBAR){
		if(par.STRAINFIELD){plot.StrainColorBar();}
		if(par.STRESSFIELD && !par.TESTNEWFAPLOT){plot.StressColorBar();}
		if(par.STRESSFIELD && par.TESTNEWFAPLOT){
plot.StressColorBar_CTB();
//plot.StressColorBar();
}
		}


		cout << fname << endl;

                plot.Write(fname);

		plot.DeleteStrainForce();



           
			//write_forces(pn,incr);
			//write_disps(pn,incr);
            		//write_pstrain(pv,pn,incr);



		}




        for (int c=0;c<NRc;c++) {
            fprintf(of, "%d %lf %lf\n", c, (double)csumx[c]/(double)csize[c], (double)csumy[c]/(double)csize[c]);
        }


  if (!(incr%(par.WSTRIDE)))
	{



		if(par.WRATIOPA)
		{
			CalcPerimeters(pv,cper,NRc);
			write_ratiopa(cper,csize,NRc,incr,ofppa);

		}

		if(par.WLENGTH)
		{
			CalcLengths(pv,pn,clength,ecc,cangle,NRc,csize,csumx,csumy);
			write_length(clength,NRc,incr,ofpl);
		}
		if(par.WECC)
		{
			CalcLengths(pv,pn,clength,ecc,cangle,NRc,csize,csumx,csumy);
			write_eccentricity(ecc,NRc,incr,ofpe);
		}

		if(par.WANGLE)
		{
			CalcLengths(pv,pn,clength,ecc,cangle,NRc,csize,csumx,csumy);
			write_cangle(cangle,NRc,incr,ofpan);
		}
		if(par.WAREA)
		{
			write_area(csize,NRc,incr,ofpa);
		}
		if(par.WSIGMA)
		{
			write_cells(pv,incr,scurrentdir);
		}




		if(par.WTOTSHAPE & NRc==1)
		{
			write_totshape(totshape,incr,scurrentdir);
		}

	        if(par.WFA)
		{
		write_fa(FA,incr,scurrentdir);
		}



		if(par.WSQDIS)
		{
			for(int c=0;c<NRc;c++)
			{
				sqdis[c] = (cmxi[c]-csumx[c]/csize[c])*(cmxi[c]-csumx[c]/csize[c])+(cmyi[c]-csumy[c]/csize[c])*(cmyi[c]-csumy[c]/csize[c]);
			}
			write_sqdis(sqdis,NRc,incr,ofpsd);
		}
		if(par.TWOCELL & par.WTWOCELLCONTACT)
		{	
			contact=check_contact(pv);
			write_twocellcontact(contact,incr,ofptcc);
		}
	}


	for(int c=0;c<NRc;c++){sumFA2[c]=0;}
	for(int n=0;n<NV;n++)
	{
		if(FA[n]>0)
		{
			sumFA2[pv[n].ctag-1]=sumFA2[pv[n].ctag-1]+FA[n];
		}
	}
	cout << endl;

	cout << "sumFA2 " << sumFA2[0] << endl;
	cout << "sumFA " << sumFA[0] << endl;
	cout << "cellmoves " << endl;

		//cout<<endl;
		//cout << "move1 " << endl;
		//now = time(0);
		//cout << now << endl;
				CPM_moves(pv,pn,csize, csumx, csumy, incr,FA,NRc,celltypes,sumFA,pnold);


		//cout << "move2 " << endl;
		//now = time(0);
		//cout << now << endl;
	for(int c=0;c<NRc;c++){sumFA2[c]=0;}
	for(int n=0;n<NV;n++)
	{
		if(FA[n]>0)
		{
			sumFA2[pv[n].ctag-1]=sumFA2[pv[n].ctag-1]+FA[n];
		}
	}
	cout << "sumFA2 " << sumFA2[0] << endl;
	cout << "sumFA " << sumFA[0] << endl;
	cout << endl;



	int nrfa=0; int nrupdated=0;
	for(int n=0;n<NV;n++) 
        {
	if(FA[n]<1 & pv[n].ctag>0){nrfa++;}
	}
        while(nrupdated<nrfa)
        {

	int n = mt_random()%NV;
		//lowerright voxel

			if(pv[n].ctag)
			{
				if(FA[n]<par.BASEFA) 
				{
						sumFA[pv[n].ctag-1]=sumFA[pv[n].ctag-1]-FA[n];

						
						if(sumFA[pv[n].ctag-1]+par.BASEFA<=par.CAPACITYFA)
						{
							FA[n]=par.BASEFA;
							sumFA[pv[n].ctag-1]=sumFA[pv[n].ctag-1]+par.BASEFA;
						}
						else
						{
							FA[n]=par.CAPACITYFA-sumFA[pv[n].ctag-1];
							sumFA[pv[n].ctag-1]=sumFA[pv[n].ctag-1]+FA[n];
						}
						nrupdated++;


				}
			}

		
	}



	for(int vx=0; vx<par.NVX; vx++) {
        for (int vy=0; vy<par.NVY; vy++) {
		
	    int v = vx+vy*par.NVX;
	    if(pv[v].ctag==1)
		{
		    int cmx = round(csumx[0]/csize[0]); 
		    int cmy = round(csumy[0]/csize[0]); 
		    int disx = cmx-vx; 
		    int disy =cmx-vy; 
		    int x = par.NVX+1+disx; 
		    int y = par.NVY+1+disy; 
		    int v2 = x+y*(2*(par.NVX)+1); 
		    totshape[v2]=totshape[v2]+1;
		}
	

	}}

		set_disp_of_prev_incr(pn,uold);
		disp_to_nodes(pnold,uold);


      

		}
    fclose(of);
	/// END ///
	printf("\nSIMULATION FINISHED!");
	//free(pv); free(pn); free(klocal); free(kcol); free(kval); free(dofpos);
	delete [] pv; delete [] pn; delete [] klocal; delete [] kcol; delete [] kval; delete [] dofpos;
	delete [] csize; delete [] csumx; delete [] csumy; delete [] cmxi; delete [] cmyi;
	if(par.WRATIOPA){delete [] cper;} if(par.WLENGTH){delete [] clength;} 
	if(par.WSQDIS){delete [] sqdis;}
	if(par.WANGLE){delete [] cangle;}
	if(par.WECC){delete [] ecc;}
	delete [] pnold; delete [] uold; delete [] FA; delete [] pnstall;
	delete [] totshape; delete [] sumFA; delete [] sumFA2; //delete [] tensionold;
	delete [] Kms;
	delete [] celltypes;
	delete [] Nt;
	fclose(ofpl); fclose(ofppa); fclose(ofpe); fclose(ofpan); fclose(ofpa); fclose(ofpsd); fclose(ofptcc);
	return 0;
}
 

