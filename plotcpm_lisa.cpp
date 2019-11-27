#include "plotcpm.h"
//#include <iostream>
//using namespace std;

#include "functions.h"
#include <string>
#include <sstream>


PlotCPM::PlotCPM(int xfield, int yfield, char *moviefile) :
    QtGraphics(xfield, yfield, moviefile)
{

    int vecsize=3;
    // allocate field for strain magnitude visualization
    strain = new double**[par.NVX];
    strain[0] = new double*[par.NVX*par.NVY];
    strain[0][0] = new double[vecsize*par.NVX*par.NVY];

    for (int i=1;i<par.NVX;i++) {
        strain[i]=strain[i-1]+par.NVY;
    }
    for (int i=1;i<par.NVX*par.NVY;i++) {
        strain[0][i]=strain[0][i-1]+vecsize;
    }

    for (int i=0;i<vecsize*par.NVX*par.NVY;i++) {
        strain[0][0][i]=0.;
    }

    force = new double**[NNX];
    force[0] = new double*[NNX*NNY];
    force[0][0] = new double[vecsize*NNX*NNY];

    for (int i=1;i<NNX;i++) {
        force[i]=force[i-1]+NNY;
    }
    for (int i=1;i<NNX*NNY;i++) {
        force[0][i]=force[0][i-1]+vecsize;
    }

    for (int i=0;i<vecsize*NNX*NNY;i++) {
        force[0][0][i]=0.;
    }

    tension = new double**[NNX];
    tension[0] = new double*[NNX*NNY];
    tension[0][0] = new double[vecsize*NNX*NNY];

    for (int i=1;i<NNX;i++) {
        tension[i]=tension[i-1]+NNY;
    }
    for (int i=1;i<NNX*NNY;i++) {
        tension[0][i]=tension[0][i-1]+vecsize;
    }

    for (int i=0;i<vecsize*NNX*NNY;i++) {
        tension[0][0][i]=0.;
    }

    // allocate field for stress magnitude visualization
    stress = new double**[par.NVX];
    stress[0] = new double*[par.NVX*par.NVY];
    stress[0][0] = new double[vecsize*par.NVX*par.NVY];

    for (int i=1;i<par.NVX;i++) {
        stress[i]=stress[i-1]+par.NVY;
    }
    for (int i=1;i<par.NVX*par.NVY;i++) {
        stress[0][i]=stress[0][i-1]+vecsize;
    }

    for (int i=0;i<vecsize*par.NVX*par.NVY;i++) {
        stress[0][0][i]=0.;
    }




}

void PlotCPM::Plot(VOX *pv, int* celltypes) {
    BeginScene();
    for (int vx = 0; vx < par.NVX-1; vx++ ) {
        for (int vy = 0; vy < par.NVY-1; vy++ ) {
            int v=vx + vy*par.NVX;

int celltype = celltypes[pv[v].ctag-1];

	int colour;

            if (pv[v].ctag<=0) {
                   colour=-1; //-1
            } else {
	      if(par.CELLCOLOUR){colour =Qt::green;if(celltype==2){colour =15;}} // make cells grey (yes...Qt::green is grey?)
		else{colour=-1;}
            }
            
	    for(int pixx = 0; pixx < (par.PIXPERVOX)-1; pixx++){ 	
		for(int pixy = 0; pixy < (par.PIXPERVOX)-1; pixy++){ 	
            		Point( colour, (par.PIXPERVOX)*vx+pixx, (par.PIXPERVOX)*vy+pixy);
		}
	    }
		

            if ( pv[v].ctag != pv[v+1].ctag )  // if cellborder  etc. etc.
            {
		for(int pixy = 0; pixy < (par.PIXPERVOX)-1; pixy++){ 	
		for(int l = 0; l < (par.LINEWIDTH); l++){
                Point( 1, (par.PIXPERVOX)*vx+(par.PIXPERVOX)-1-l, (par.PIXPERVOX)*vy +pixy -1*l);
		}}

            }
            else{
		for(int pixy = 0; pixy < (par.PIXPERVOX)-1; pixy++){ 	
                Point( colour, (par.PIXPERVOX)*vx+(par.PIXPERVOX)-1, (par.PIXPERVOX)*vy +pixy );
		}
	   }


            if ( pv[v].ctag != pv[v+par.NVX].ctag ) {
		for(int pixx = 0; pixx < (par.PIXPERVOX)-1; pixx++){ 	
		for(int l = 0; l < (par.LINEWIDTH); l++){
                Point( 1, (par.PIXPERVOX)*vx+pixx-1*l, (par.PIXPERVOX)*vy+(par.PIXPERVOX)-1-l );
		}}
            } else
		{
		for(int pixx = 0; pixx < (par.PIXPERVOX)-1; pixx++){ 	
                Point( colour, (par.PIXPERVOX)*vx+pixx, (par.PIXPERVOX)*vy+(par.PIXPERVOX)-1 );
		}
		}

            // Cells that touch eachother's corners are NO neighbours //

            if (pv[v].ctag!=pv[v+par.NVX+1].ctag
                    || pv[v+1].ctag!=pv[v+par.NVX].ctag ) {
		for(int l = 0; l < (par.LINEWIDTH); l++){
                Point( 1, (par.PIXPERVOX)*vx+(par.PIXPERVOX)-1-l, (par.PIXPERVOX)*vy+(par.PIXPERVOX)-1-l );
            }}
            else
                Point( colour, (par.PIXPERVOX)*vx+(par.PIXPERVOX)-1, (par.PIXPERVOX)*vy+(par.PIXPERVOX)-1 );

        }
    }
    EndScene();
}
/*
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
*/


void PlotCPM::PlotHueStrainField(bool STRAINFIELD) {


    BeginScene();
QColor c;
	    double max_strain;
	    //if(par.MAXCOLORBAR<=0.01){par.MAXCOLORBAR = MaxStrainMagnitude();}
	    //if(par.MAXCOLORBAR<=0.01){par.MAXCOLORBAR=0.01;}
	    if(par.COLORBAR)
	    {
		max_strain=par.MAXCOLORBARSTRAIN;
		if(par.MAXCOLORBARSTRAIN==0){max_strain=MaxStrainMagnitude();}

	    }
	    else{max_strain=MaxStrainMagnitude();}

if(max_strain==0){ max_strain=0.0001;}

c=Qt::white;

    // Plot vector field

    // Draw strain magnitude

    for (int vx=0;vx<par.NVX;vx++) {
        for (int vy=0;vy<par.NVY;vy++) {

	if(par.STRAINFIELD)
	{	
            int hue=240-240*strain[vx][vy][2]/max_strain;
	    if(strain[vx][vy][2]>max_strain){hue=0;}
            if (hue<0) {
                cerr << "Panic. Hue is: " << hue << endl;
                cerr << "Strain at ( " << vx << ", " << vy << ") is " << strain[vx][vy][2] << ", max strain: " << max_strain << endl;
                exit(0);
            }
            c.setHsv(hue,255,255);
	}
        picture->setPen( c );

	
	    for(int pixx = 0; pixx < (par.PIXPERVOX); pixx++){ 	
		for(int pixy = 0; pixy < (par.PIXPERVOX); pixy++){ 	
            		picture->drawPoint((par.PIXPERVOX)*vx+pixx, (par.PIXPERVOX)*vy+pixy);
		}
	    }

	    
        }
    }
EndScene();
}

void PlotCPM::PlotHueChemField(double *chem) {


    BeginScene();
QColor c;


	    double max_chem;
	    //if(par.MAXCOLORBAR<=0.01){par.MAXCOLORBAR = MaxStrainMagnitude();}
	    //if(par.MAXCOLORBAR<=0.01){par.MAXCOLORBAR=0.01;}
	    if(par.COLORBAR)
	    {
		max_chem=par.MAXCOLORBARDENS;
		if(par.MAXCOLORBARDENS==0){for (int v=0;v<NV;v++){if(chem[v]>max_chem){max_chem=chem[v];}}}

	    }
	    else{for (int v=0;v<NV;v++){if(chem[v]>max_chem){max_chem=chem[v];}}}

if(max_chem==0){ max_chem=0.0001;}

cout << "max_chem " << max_chem << endl;


c=Qt::white;

    // Plot vector field

    // Draw strain magnitude

    for (int vx=0;vx<par.NVX;vx++) {
        for (int vy=0;vy<par.NVY;vy++) {

	int v=vx+par.NVX*vy;
            int hue=240-240*chem[v]/max_chem;
		if(chem[v]>max_chem){hue=0;}

            if (hue<0) {
                cerr << "Panic. Hue is: " << hue << endl;
                cerr << "Chem at ( " << vx << ", " << vy << ") is " << chem[v] << ", max chem: " << max_chem << endl;
                exit(0);
            }
            c.setHsv(hue,255,255);
	
        picture->setPen( c );



	
	    for(int pixx = 0; pixx < (par.PIXPERVOX); pixx++){ 	
		for(int pixy = 0; pixy < (par.PIXPERVOX); pixy++){ 	
            		picture->drawPoint((par.PIXPERVOX)*vx+pixx, (par.PIXPERVOX)*vy+pixy);
		}
	    }

	    
        }
    }
EndScene();
}


void PlotCPM::PlotHueStressField(bool STRESSFIELD) {


    BeginScene();
QColor c;
	    double max_stress;
	    //if(par.MAXCOLORBAR<=0.01){par.MAXCOLORBAR = MaxStrainMagnitude();}
	    //if(par.MAXCOLORBAR<=0.01){par.MAXCOLORBAR=0.01;}
	    if(par.COLORBAR)
	    {
		max_stress=par.MAXCOLORBARSTRESS;
		if(par.MAXCOLORBARSTRESS==0){max_stress=MaxStressMagnitude();}

	    }
	    else{max_stress=MaxStressMagnitude();}

//if(max_stress==0){ max_stress=0.0001;}



c=Qt::white;

    // Plot vector field

    // Draw strain magnitude

    for (int vx=0;vx<par.NVX;vx++) {
        for (int vy=0;vy<par.NVY;vy++) {

	if(par.STRESSFIELD)
	{	


            int hue=240-240*fabs(stress[vx][vy][2])/max_stress;
	    if(stress[vx][vy][2]>max_stress){hue=0;}
            if (hue<0) {
                cerr << "Panic. Hue is: " << hue << endl;
                cerr << "stress at ( " << vx << ", " << vy << ") is " << stress[vx][vy][2] << ", max stress: " << max_stress << endl;
                exit(0);
            }
            c.setHsv(hue,255,255);
	}
        picture->setPen( c );

	
	    for(int pixx = 0; pixx < (par.PIXPERVOX); pixx++){ 	
		for(int pixy = 0; pixy < (par.PIXPERVOX); pixy++){ 	
            		picture->drawPoint((par.PIXPERVOX)*vx+pixx, (par.PIXPERVOX)*vy+pixy);
		}
	    }

	    
        }
    }
EndScene();
}


void PlotCPM::PlotPrincipleStrainField(VOX *pv) {

BeginScene();


    // calculate strains
    //CalculateStrainField(pn);

    // draw strain lines
    for (int vx=0;vx<par.NVX;vx++) { 
        for (int vy=0;vy<par.NVY;vy++) {


	    double max_strain;
	    if(par.COLORBAR)
	    {
		max_strain=par.MAXCOLORBARSTRAIN;
		if(par.MAXCOLORBARSTRAIN==0){max_strain=MaxStrainMagnitude();}

	    }
	    else{max_strain=MaxStrainMagnitude();}

if(max_strain==0){ max_strain=0.0001;}

            double linelength = sqrt(2)*(par.PIXPERVOX)/max_strain; 
            double x1= ((par.PIXPERVOX)*vx+(par.PIXPERVOX)/2-linelength*strain[vx][vy][0]/2);
            double y1= ((par.PIXPERVOX)*vy+(par.PIXPERVOX)/2-linelength*strain[vx][vy][1]/2);
            double x2= ((par.PIXPERVOX)*vx+(par.PIXPERVOX)/2+linelength*strain[vx][vy][0]/2);
            double y2= ((par.PIXPERVOX)*vy+(par.PIXPERVOX)/2+linelength*strain[vx][vy][1]/2);


            if (x1<0) x1=0;
            if (x1>(par.PIXPERVOX)*par.NVX-1) x1=(par.PIXPERVOX)*par.NVX-1;
            if (y1<0) y1=0;
            if (y1>(par.PIXPERVOX)*par.NVY-1) y1=(par.PIXPERVOX)*par.NVY-1;

            if (x2<0) x2=0;
            if (x2>(par.PIXPERVOX)*par.NVX-1) x2=(par.PIXPERVOX)*par.NVX-1;
            if (y2<0) y2=0;
            if (y2>(par.PIXPERVOX)*par.NVY-1) y2=(par.PIXPERVOX)*par.NVY-1;

            // And draw it :-)

            //strain_magnitude[vx][vy] = sqrt(strx*strx + stry*stry);
            int v=vx + vy*par.NVX;

		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
	      		Line(x1+l,y1,x2+l,y2,1);
			Line(x1-l,y1,x2-l,y2,1);
		}


        }
    }


    EndScene();
}
void PlotCPM::PlotPrincipleStressField(VOX *pv) {

BeginScene();


    // calculate strains
    //CalculateStrainField(pn);

    // draw strain lines
    for (int vx=0;vx<par.NVX;vx++) { 
        for (int vy=0;vy<par.NVY;vy++) {


	    double max_stress;
	    if(par.COLORBAR)
	    {
		max_stress=par.MAXCOLORBARSTRESS;
		if(par.MAXCOLORBARSTRESS==0){max_stress=MaxStressMagnitude();}

	    }
	    else{max_stress=MaxStressMagnitude();}

if(max_stress==0){ max_stress=0.0001;}

            double linelength = sqrt(2)*(par.PIXPERVOX)/max_stress; 
            double x1= ((par.PIXPERVOX)*vx+(par.PIXPERVOX)/2-linelength*stress[vx][vy][0]/2);
            double y1= ((par.PIXPERVOX)*vy+(par.PIXPERVOX)/2-linelength*stress[vx][vy][1]/2);
            double x2= ((par.PIXPERVOX)*vx+(par.PIXPERVOX)/2+linelength*stress[vx][vy][0]/2);
            double y2= ((par.PIXPERVOX)*vy+(par.PIXPERVOX)/2+linelength*stress[vx][vy][1]/2);


            if (x1<0) x1=0;
            if (x1>(par.PIXPERVOX)*par.NVX-1) x1=(par.PIXPERVOX)*par.NVX-1;
            if (y1<0) y1=0;
            if (y1>(par.PIXPERVOX)*par.NVY-1) y1=(par.PIXPERVOX)*par.NVY-1;

            if (x2<0) x2=0;
            if (x2>(par.PIXPERVOX)*par.NVX-1) x2=(par.PIXPERVOX)*par.NVX-1;
            if (y2<0) y2=0;
            if (y2>(par.PIXPERVOX)*par.NVY-1) y2=(par.PIXPERVOX)*par.NVY-1;

            // And draw it :-)

            //strain_magnitude[vx][vy] = sqrt(strx*strx + stry*stry);
            int v=vx + vy*par.NVX;

		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
	      		Line(x1+l,y1,x2+l,y2,1);
			Line(x1-l,y1,x2-l,y2,1);
		}


        }
    }


    EndScene();
}

double PlotCPM::MaxStrainMagnitude(void) const {
    double max=0.;
    for (int i=0;i<par.NVX*par.NVY;i++) {
        //for (int vx=0;vx<par.NVX-1;vx++) {
        //for (int vy=0;vy<par.NVY-1;vy++) {


        double magn=strain[0][i][2];
        if (magn>max) {
            max=magn;
        }
    }
    // }
    return max;
}

double PlotCPM::MaxStressMagnitude(void) const {
    double max=0.;
    for (int i=0;i<par.NVX*par.NVY;i++) {
        //for (int vx=0;vx<par.NVX-1;vx++) {
        //for (int vy=0;vy<par.NVY-1;vy++) {


        double magn=fabs(stress[0][i][2]);
        if (magn>max) {
            max=magn;
        }
    }
    // }
    return max;
}


void PlotCPM::CalculateStressField(NOD *pn) const{
    for (int vx=0;vx<par.NVX;vx++) {
        for (int vy=0;vy<par.NVY;vy++) {
            int v=vx + vy*par.NVX;

            double estrains[3],L1,L2,v1[2],v2[2];
            double estress[3];

            get_estrains(pn,v,estrains);
	    get_estress(1,estrains,estress);
		
	double etractionstress[2];

            L1=L2=.0; get_princs(estress,&L1,&L2,v1,v2,0);
		get_etractionstress(pn,v, etractionstress);

            //if(L1>=L2)
            //{

		if(L1>0){
                // calculate strain vector

                stress[vx][vy][0] = L1 * v1[0];
                stress[vx][vy][1] = L1 * v1[1];

		}
            //}
            //if(L2>L1)
            //{

		//if(L2>0){
                //stress[vx][vy][0] = L2 * v2[0];
                //stress[vx][vy][1] = L2 * v2[1];}

            //}
            // magnitude

	    if(!par.TRACTIONSTRESSFIELD)
	    {
		    double strx=stress[vx][vy][0], stry=stress[vx][vy][1];
		    stress[vx][vy][2]=sqrt(strx*strx + stry*stry);
	 	    if(L1<0){stress[vx][vy][2]=-sqrt(strx*strx + stry*stry);}
		    if(par.HYDSTRESS)
		    {
			stress[vx][vy][2]=(estress[0]+estress[1])/2;
			if(stress[vx][vy][2]<0){stress[vx][vy][2]=0;}
		    }

  	    }
	    

	        if(par.TRACTIONSTRESSFIELD){
		double strx=etractionstress[0]; double stry=etractionstress[1];
		stress[vx][vy][2]=sqrt(strx*strx + stry*stry);}



	}

   }
}


void PlotCPM::CalculateForceField(NOD *pn) const{
    for (int nx=1;nx<NNX-1;nx++) {
        for (int ny=1;ny<NNY-1;ny++) {
            int n=nx + ny*NNX;

                // calculate force vector

                force[nx][ny][0] = pn[n].fx;
                force[nx][ny][1] = pn[n].fy;

            // magnitude

            double fx=force[nx][ny][0], fy=force[nx][ny][1];
            force[nx][ny][2]=sqrt(fx*fx + fy*fy)*par.THICKNESS;


	}
   }
}

void PlotCPM::CalculateTensionField(NOD *pn) const{
    for (int nx=1;nx<NNX-1;nx++) {
        for (int ny=1;ny<NNY-1;ny++) {
            int n=nx + ny*NNX;

                // calculate force vector

                tension[nx][ny][0] = pn[n].tx;
                tension[nx][ny][1] = pn[n].ty;

            // magnitude

            double tx=tension[nx][ny][0], ty=tension[nx][ny][1];
            tension[nx][ny][2]=sqrt(tx*tx + ty*ty);


	}
   }
}

double PlotCPM::MaxForceMagnitude(void) const {
    double max=0.;
    for (int i=0;i<NNX*NNY;i++) {
        double magn=force[0][i][2];
        if (magn>max) {
            max=magn;
        }
    }

    return max;
}

double PlotCPM::MaxTensionMagnitude(void) const {
    double max=0.;
    for (int i=0;i<NNX*NNY;i++) {
        double magn=tension[0][i][2];
        if (magn>max) {
            max=magn;
        }
    }
    return max;
}

void PlotCPM::CalculateStrainField(NOD *pn) const{
    for (int vx=0;vx<par.NVX;vx++) {
        for (int vy=0;vy<par.NVY;vy++) {
            int v=vx + vy*par.NVX;

            double estrains[3],L1,L2,v1[2],v2[2];

            get_estrains(pn,v,estrains);

            L1=L2=.0; get_princs(estrains,&L1,&L2,v1,v2,1);

//if(L1<0 & L2 <0){cout << "vx " << vx << endl;cout << "vy " << vy << endl;}

            if(L1>L2)
            {

		if(L1>0)
		{ //plot only compression
                // calculate strain vector

                strain[vx][vy][0] = L1 * v1[0];
                strain[vx][vy][1] = L1 * v1[1];

		}
            }
            else
            {

		if(L2>0)
		{
                strain[vx][vy][0] = L2 * v2[0];
                strain[vx][vy][1] = L2 * v2[1];
		}

            }
            // magnitude

            double strx=strain[vx][vy][0], stry=strain[vx][vy][1];
            strain[vx][vy][2]=sqrt(strx*strx + stry*stry);

	}
   }
}

void PlotCPM::PlotNodalForces(NOD *pn){


BeginScene();

    // draw force lines
    for (int nx=1;nx<NNX-1;nx++) { 
        for (int ny=1;ny<NNY-1;ny++) {



            double linelength = sqrt(2)*(par.PIXPERVOX)/MaxForceMagnitude(); 



	    if(par.MAXFORCE){linelength = sqrt(2)*(par.PIXPERVOX)/par.MAXFORCE;}
            double x1= (par.PIXPERVOX)*nx;
            double y1= (par.PIXPERVOX)*ny;
            double x2= ((par.PIXPERVOX)*nx+linelength*force[nx][ny][0]);
            double y2= ((par.PIXPERVOX)*ny+linelength*force[nx][ny][1]);



            if (x1<(par.PIXPERVOX)) x1=(par.PIXPERVOX);
            if (x1>(par.PIXPERVOX)*par.NVX-1) x1=(par.PIXPERVOX)*par.NVX-1;
            if (y1<(par.PIXPERVOX)) y1=(par.PIXPERVOX);
            if (y1>(par.PIXPERVOX)*par.NVY-1) y1=(par.PIXPERVOX)*par.NVY-1;

            if (x2<(par.PIXPERVOX)) x2=(par.PIXPERVOX);
            if (x2>(par.PIXPERVOX)*par.NVX-1) x2=(par.PIXPERVOX)*par.NVX-1;
            if (y2<(par.PIXPERVOX)) y2=(par.PIXPERVOX);
            if (y2>(par.PIXPERVOX)*par.NVY-1) y2=(par.PIXPERVOX)*par.NVY-1;

            int n=nx + ny*NNX;

if(!(x1==x2) || !(y1==y2)){

		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
	      		Line(x1+l,y1,x2+l,y2,3);
			Line(x1-l,y1,x2-l,y2,3);
		}
		
	int diffx = x2-x1;
	int diffy = y2-y1;
	double length = sqrt(diffx*diffx+diffy*diffy);

	if(length>0)
	{
		double alpha = acos(diffx/length);
		double gamma = (15.0/180)*3.1416;
		double beta = 3.1416-gamma;
		double lengtha = length/5;
		double lengthv = (length-lengtha*cos(gamma))/cos(beta);
		if(lengthv<0){lengthv=-lengthv;}
		if(y2>y1){
		double x3 = x1 - lengthv*cos(alpha-beta);
		double y3 = y1 - lengthv*sin(alpha-beta);
		double x4 = x1 - lengthv*cos(alpha+beta);
		double y4 = y1 - lengthv*sin(alpha+beta);
		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
		Line(x2+l,y2,x3+l,y3,3);
		Line(x2-l,y2,x4-l,y4,3); }}
		if(y2<y1){
		double x3 = x1 - lengthv*cos(alpha-beta);
		double y3 = y1 + lengthv*sin(alpha-beta);
		double x4 = x1 - lengthv*cos(alpha+beta);
		double y4 = y1 + lengthv*sin(alpha+beta);
		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
		Line(x2+l,y2,x3+l,y3,3);
		Line(x2-l,y2,x4-l,y4,3);}}
	}
	
        }
    }
}


    EndScene();


}

void PlotCPM::PlotNodalDeform(NOD *pn){

BeginScene();

    // draw force lines
    for (int nx=1;nx<NNX-1;nx++) { 
        for (int ny=1;ny<NNY-1;ny++) {

	    int n = nx+ny*NNX;
	    double magn = sqrt(pn[n].ux*pn[n].ux+pn[n].uy*pn[n].uy);

            double linelength = sqrt(2)*(par.PIXPERVOX)/magn; 
	    if(par.MAXDEFORM){linelength = sqrt(2)*(par.PIXPERVOX)/par.MAXDEFORM;}
            double x1= (par.PIXPERVOX)*nx;
            double y1= (par.PIXPERVOX)*ny;
            double x2= ((par.PIXPERVOX)*nx+linelength*pn[n].ux);
            double y2= ((par.PIXPERVOX)*ny+linelength*pn[n].uy);

            if (x1<(par.PIXPERVOX)) x1=(par.PIXPERVOX);
            if (x1>(par.PIXPERVOX)*par.NVX-1) x1=(par.PIXPERVOX)*par.NVX-1;
            if (y1<(par.PIXPERVOX)) y1=(par.PIXPERVOX);
            if (y1>(par.PIXPERVOX)*par.NVY-1) y1=(par.PIXPERVOX)*par.NVY-1;

            if (x2<(par.PIXPERVOX)) x2=(par.PIXPERVOX);
            if (x2>(par.PIXPERVOX)*par.NVX-1) x2=(par.PIXPERVOX)*par.NVX-1;
            if (y2<(par.PIXPERVOX)) y2=(par.PIXPERVOX);
            if (y2>(par.PIXPERVOX)*par.NVY-1) y2=(par.PIXPERVOX)*par.NVY-1;

if(!(x1==x2) || !(y1==y2)){

		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
	      		Line(x1+l,y1,x2+l,y2,2);
			Line(x1-l,y1,x2-l,y2,2);
		}
		
	int diffx = x2-x1;
	int diffy = y2-y1;
	double length = sqrt(diffx*diffx+diffy*diffy);

	if(length>0)
	{
		double alpha = acos(diffx/length);
		double gamma = (15.0/180)*3.1416;
		double beta = 3.1416-gamma;
		double lengtha = length/5;
		double lengthv = (length-lengtha*cos(gamma))/cos(beta);
		if(lengthv<0){lengthv=-lengthv;}
		if(y2>y1){
		double x3 = x1 - lengthv*cos(alpha-beta);
		double y3 = y1 - lengthv*sin(alpha-beta);
		double x4 = x1 - lengthv*cos(alpha+beta);
		double y4 = y1 - lengthv*sin(alpha+beta);
		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
		Line(x2+l,y2,x3+l,y3,2);
		Line(x2-l,y2,x4-l,y4,2);}}
		if(y2<y1){
		double x3 = x1 - lengthv*cos(alpha-beta);
		double y3 = y1 + lengthv*sin(alpha-beta);
		double x4 = x1 - lengthv*cos(alpha+beta);
		double y4 = y1 + lengthv*sin(alpha+beta);
		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
		Line(x2+l,y2,x3+l,y3,2);
		Line(x2-l,y2,x4-l,y4,2);}}
	}
	
        }
    }
}


    EndScene();


}

void PlotCPM::PlotNodalFA(NOD *pn,double* FA){


BeginScene();

    // draw FA circles
    for (int nx=0;nx<par.NVX-1;nx++) { 
        for (int ny=0;ny<par.NVY-1;ny++) {

	int radius=0;
//double maxfa=20000;
double maxfa = 50*par.VOXSIZE*par.VOXSIZE/8e-15;
if(par.MAXFA){maxfa=par.MAXFA;}
//maxfa=20000;
double minfa = par.BASEFA;
            int n=nx + ny*par.NVX;
            if(FA[n]>=minfa){radius = 1*(FA[n]-minfa)*par.PIXPERVOX/(maxfa-minfa);
		radius=radius+0.1*par.PIXPERVOX;}
		if(FA[n]>=maxfa){radius = 1*par.PIXPERVOX;
			radius=radius+0.1*par.PIXPERVOX;}

            int x1= (par.PIXPERVOX)*nx+par.PIXPERVOX/2;
            int y1= (par.PIXPERVOX)*ny+par.PIXPERVOX/2;

	    QPoint Qp=QPoint(x1,y1);



	if((nx%par.PATTERNC==1) && (ny%par.PATTERNC==1) && par.PATTERN)
	{
	    double radius2 = 0.2*par.PIXPERVOX;
	    picture->setBrush(Qt::red);
	    this->picture->drawEllipse(Qp, (int)radius2,(int)radius2);

	}


	    if(par.FACOLOUR)
	    {
            int hue=240-240*(FA[n]-5000)/(maxfa-5000);
	    if(FA[n]<=5000){hue=240;}
		QColor c;
            c.setHsv(hue,255,255);
	
            picture->setPen( c );
	    picture->setBrush(c);
	    }
	    if(!par.FACOLOUR){picture->setBrush(Qt::black);}
	    this->picture->drawEllipse(Qp, radius,radius);
/*	if(n==NV/2+par.NVX/2 | n==NV/2+par.NVX/2-11*par.NVX |n==NV/2+par.NVX/2-5*par.NVX | n==NV/2+par.NVX/2+8*par.NVX| n==NV/2+par.NVX/2-8*par.NVX){

		picture->setBrush(Qt::red);
	    this->picture->drawEllipse(Qp, radius,radius);

		}*/
	    this->picture->drawEllipse(Qp, radius,radius);
		/*if(n==18504-1*10| n ==18504+9*200+5*200 | n==18504 |n ==18504+9*200+5*200-10 ){

		picture->setBrush(Qt::red);
	    this->picture->drawEllipse(Qp, radius,radius);

		}*/

//if(stress[nx][ny][2]>0){
//radius=1;	    this->picture->drawEllipse(Qp, radius,radius);
//cout << "stress[nx][ny][2] " << stress[nx][ny][2] << endl;
//}



}
}


    EndScene();


}

void PlotCPM::PlotNodalTension(NOD *pn){


BeginScene();

    // draw force lines
    for (int nx=1;nx<NNX-1;nx++) { 
        for (int ny=1;ny<NNY-1;ny++) {



            double linelength = sqrt(2)*(par.PIXPERVOX)/MaxTensionMagnitude(); 


	    if(par.MAXFORCE){linelength = sqrt(2)*(par.PIXPERVOX)/par.MAXFORCE;}
            double x1= (par.PIXPERVOX)*nx;
            double y1= (par.PIXPERVOX)*ny;
            double x2= ((par.PIXPERVOX)*nx+linelength*tension[nx][ny][0]);
            double y2= ((par.PIXPERVOX)*ny+linelength*tension[nx][ny][1]);


            if (x1<(par.PIXPERVOX)) x1=(par.PIXPERVOX);
            if (x1>(par.PIXPERVOX)*par.NVX-1) x1=(par.PIXPERVOX)*par.NVX-1;
            if (y1<(par.PIXPERVOX)) y1=(par.PIXPERVOX);
            if (y1>(par.PIXPERVOX)*par.NVY-1) y1=(par.PIXPERVOX)*par.NVY-1;

            if (x2<(par.PIXPERVOX)) x2=(par.PIXPERVOX);
            if (x2>(par.PIXPERVOX)*par.NVX-1) x2=(par.PIXPERVOX)*par.NVX-1;
            if (y2<(par.PIXPERVOX)) y2=(par.PIXPERVOX);
            if (y2>(par.PIXPERVOX)*par.NVY-1) y2=(par.PIXPERVOX)*par.NVY-1;

            int n=nx + ny*NNX;

if(!(x1==x2) || !(y1==y2)){

		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
	      		Line(x1+l,y1,x2+l,y2,4);
			Line(x1-l,y1,x2-l,y2,4);
			if(n==NN/2-11*NNX | n==NN/2 | n==NN/2+4 | n == NN/2+4 + 5 + NNX){Line(x1+l,y1,x2+l,y2,5);
			Line(x1-l,y1,x2-l,y2,5);}
		}
		
	int diffx = x2-x1;
	int diffy = y2-y1;
	double length = sqrt(diffx*diffx+diffy*diffy);

	if(length>0)
	{
		double alpha = acos(diffx/length);
		double gamma = (15.0/180)*3.1416;
		double beta = 3.1416-gamma;
		double lengtha = length/5;
		double lengthv = (length-lengtha*cos(gamma))/cos(beta);
		if(lengthv<0){lengthv=-lengthv;}
		if(y2>y1){
		double x3 = x1 - lengthv*cos(alpha-beta);
		double y3 = y1 - lengthv*sin(alpha-beta);
		double x4 = x1 - lengthv*cos(alpha+beta);
		double y4 = y1 - lengthv*sin(alpha+beta);
		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
		Line(x2+l,y2,x3+l,y3,4);
		Line(x2-l,y2,x4-l,y4,4);}}
		if(y2<y1){
		double x3 = x1 - lengthv*cos(alpha-beta);
		double y3 = y1 + lengthv*sin(alpha-beta);
		double x4 = x1 - lengthv*cos(alpha+beta);
		double y4 = y1 + lengthv*sin(alpha+beta);
		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
		Line(x2+l,y2,x3+l,y3,4);
		Line(x2-l,y2,x4-l,y4,4);}}
	}


	
        }
    }

}




    EndScene();


}

void PlotCPM::PlotNodalTension2(NOD *pn){


BeginScene();

    // draw force lines
    for (int nx=1;nx<NNX-1;nx++) { 
        for (int ny=1;ny<NNY-1;ny++) {



            double linelength = sqrt(2)*(par.PIXPERVOX)/MaxTensionMagnitude(); 


	    if(par.MAXFORCE){linelength = sqrt(2)*(par.PIXPERVOX)/par.MAXFORCE;}
            double x1= (par.PIXPERVOX)*nx;
            double y1= (par.PIXPERVOX)*ny;
            double x2= ((par.PIXPERVOX)*nx-linelength*tension[nx][ny][0]);
            double y2= ((par.PIXPERVOX)*ny-linelength*tension[nx][ny][1]);


            if (x1<(par.PIXPERVOX)) x1=(par.PIXPERVOX);
            if (x1>(par.PIXPERVOX)*par.NVX-1) x1=(par.PIXPERVOX)*par.NVX-1;
            if (y1<(par.PIXPERVOX)) y1=(par.PIXPERVOX);
            if (y1>(par.PIXPERVOX)*par.NVY-1) y1=(par.PIXPERVOX)*par.NVY-1;

            if (x2<(par.PIXPERVOX)) x2=(par.PIXPERVOX);
            if (x2>(par.PIXPERVOX)*par.NVX-1) x2=(par.PIXPERVOX)*par.NVX-1;
            if (y2<(par.PIXPERVOX)) y2=(par.PIXPERVOX);
            if (y2>(par.PIXPERVOX)*par.NVY-1) y2=(par.PIXPERVOX)*par.NVY-1;

            int n=nx + ny*NNX;

if(!(x1==x2) || !(y1==y2)){

		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
	      		Line(x1+l,y1,x2+l,y2,4);
			Line(x1-l,y1,x2-l,y2,4);
			if(n==NN/2-11*NNX | n==NN/2 | n==NN/2+4 | n == NN/2+4 + 5 + NNX){Line(x1+l,y1,x2+l,y2,5);
			Line(x1-l,y1,x2-l,y2,5);}
		}
		
	int diffx = x2-x1;
	int diffy = y2-y1;
	double length = sqrt(diffx*diffx+diffy*diffy);

	if(length>0)
	{
		double alpha = acos(diffx/length);
		double gamma = (15.0/180)*3.1416;
		double beta = 3.1416-gamma;
		double lengtha = length/5;
		double lengthv = (length-lengtha*cos(gamma))/cos(beta);
		if(lengthv<0){lengthv=-lengthv;}
		if(y2>y1){
		double x3 = x1 - lengthv*cos(alpha-beta);
		double y3 = y1 - lengthv*sin(alpha-beta);
		double x4 = x1 - lengthv*cos(alpha+beta);
		double y4 = y1 - lengthv*sin(alpha+beta);
		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
		Line(x2+l,y2,x3+l,y3,4);
		Line(x2-l,y2,x4-l,y4,4);}}
		if(y2<y1){
		double x3 = x1 - lengthv*cos(alpha-beta);
		double y3 = y1 + lengthv*sin(alpha-beta);
		double x4 = x1 - lengthv*cos(alpha+beta);
		double y4 = y1 + lengthv*sin(alpha+beta);
		for(int l = 0; l < (par.LINEWIDTH); l++)
		{
		Line(x2+l,y2,x3+l,y3,4);
		Line(x2-l,y2,x4-l,y4,4);}}
	}


	
        }
    }

}




    EndScene();


}

void PlotCPM::PlotNodeConnection(NOD *pn, int *pathx, int* pathy){

    BeginScene();

		//for loop plot the path from node to node calculated in checknodeconnection
		int colour = 5;
		int length = par.NVX;
		for(int i = 0; i < length ; i++)
		{
			int vx = pathx[i]; int vy = pathy[i];
			if(vx>0&&vy>0)
			{
				    for(int pixx = 0; pixx < (par.PIXPERVOX)-1; pixx++){ 	
					for(int pixy = 0; pixy < (par.PIXPERVOX)-1; pixy++){ 	
				    		Point( colour, (par.PIXPERVOX)*vx+pixx, (par.PIXPERVOX)*vy+pixy);
					}
				    }
			}
		}
EndScene();

}

void PlotCPM::StrainColorBar(void){
BeginScene();
QColor c;

double strainhue=0;
double max_strain;

	    if(par.COLORBAR)
	    {
		max_strain=par.MAXCOLORBARSTRAIN;
		if(par.MAXCOLORBARSTRAIN==0){max_strain=MaxStrainMagnitude();}

	    }
	    else{max_strain=MaxStrainMagnitude();}

if(max_strain==0){ max_strain=0.0001;}

double cstep = max_strain/(par.NVY*par.PIXPERVOX);

cout << "max_strain" << max_strain << endl;

for(int i = 0; i < par.NVY*par.PIXPERVOX; i++)
{
	    strainhue = i*cstep;
            int hue=240-240*strainhue/max_strain;
            if (hue<0) {
                cerr << "Panic. Hue is: " << hue << endl;
                exit(0);
            }
            c.setHsv(hue,255,255);
	
        picture->setPen( c );
	
	double x = par.NVX*par.PIXPERVOX + par.WIDTHCOLORBAR;
	double y = par.NVY*par.PIXPERVOX-i;
	for(int pixx = 0; pixx < par.WIDTHCOLORBAR; pixx++)
	{
        	picture->drawPoint(x+pixx,y);
	}
    if(!(i%par.NVY)&&(i>0))
    {
        double val = strainhue;
        if(val>max_strain){val=par.MAXCOLORBARSTRAIN;}
        DrawColorBarLabel(i,val);
    }

}
EndScene();

}

void PlotCPM::ChemColorBar(double *chem){
BeginScene();
QColor c;

double chemhue=0;
	    double max_chem;
	    //if(par.MAXCOLORBAR<=0.01){par.MAXCOLORBAR = MaxStrainMagnitude();}
	    //if(par.MAXCOLORBAR<=0.01){par.MAXCOLORBAR=0.01;}
	    if(par.COLORBAR)
	    {
		max_chem=par.MAXCOLORBARDENS;
		if(par.MAXCOLORBARDENS==0){for (int v=0;v<NV;v++){if(chem[v]>max_chem){max_chem=chem[v];}}}

	    }
	    else{for (int v=0;v<NV;v++){if(chem[v]>max_chem){max_chem=chem[v];}}}

if(max_chem==0){ max_chem=0.0001;}

cout << "max_chem " << max_chem << endl;

double cstep = max_chem/(par.NVY*par.PIXPERVOX);


for(int i = 0; i < par.NVY*par.PIXPERVOX; i++)
{
	    chemhue = i*cstep;
            int hue=240-240*chemhue/max_chem;
            if (hue<0) {
                cerr << "Panic. Hue is: " << hue << endl;
                exit(0);
            }
            c.setHsv(hue,255,255);
	
        picture->setPen( c );
	
	double x = par.NVX*par.PIXPERVOX + par.WIDTHCOLORBAR;
	double y = par.NVY*par.PIXPERVOX-i;
	for(int pixx = 0; pixx < par.WIDTHCOLORBAR; pixx++)
	{
        	picture->drawPoint(x+pixx,y);
	}
    if(!(i%par.NVY)&&(i>0))
    {
        double val = chemhue;
        if(val>max_chem){val=par.MAXCOLORBARDENS;}
        DrawColorBarLabel(i,val);
    }

}
EndScene();

}

void PlotCPM::StressColorBar(void){
BeginScene();
QColor c;

double strainhue=0;
double max_stress;

	    if(par.COLORBAR)
	    {
		max_stress=par.MAXCOLORBARSTRESS;
		if(par.MAXCOLORBARSTRESS==0){max_stress=MaxStressMagnitude();}

	    }
	    else{max_stress=MaxStressMagnitude();}

if(max_stress==0){ max_stress=0.0001;}

double cstep = max_stress/(par.NVY*par.PIXPERVOX);

cout << "max_stress" << max_stress << endl;

for(int i = 0; i < par.NVY*par.PIXPERVOX; i++)
{
	    strainhue = i*cstep;
            int hue=240-240*strainhue/max_stress;
            if (hue<0) {
                cerr << "Panic. Hue is: " << hue << endl;
                exit(0);
            }
            c.setHsv(hue,255,255);
	
        picture->setPen( c );
	
	double x = par.NVX*par.PIXPERVOX + par.WIDTHCOLORBAR;
	double y = par.NVY*par.PIXPERVOX-i;
	for(int pixx = 0; pixx < par.WIDTHCOLORBAR; pixx++)
	{
        	picture->drawPoint(x+pixx,y);
	}
    if(!(i%par.NVY)&&(i>0))
    {
        double val = strainhue;
        if(val>max_stress){val=par.MAXCOLORBARSTRESS;}
        DrawColorBarLabel(i,val);
    }

}
EndScene();

}

void PlotCPM::DrawColorBarLabel(int yval, double val){
string valstr;
ostringstream convert;
convert.precision(4);
convert << fixed;
convert << val;
valstr = convert.str();
const char * c = valstr.c_str();
const QString & s = QString(c);
picture->setPen(QColor::fromRgb(0,0,0));
QFont* font = new QFont("Arial");
font->setPixelSize(100*par.NVX/300);
font->setBold(true);
picture->setFont(*font);
this->picture->drawText(QPointF(par.NVX*par.PIXPERVOX+2*par.WIDTHCOLORBAR,par.NVY*par.PIXPERVOX-yval), s);

}


void PlotCPM::DeleteStrainForce() {
  delete[] strain;
	delete[] force;
	delete[] tension;
	delete[] stress;
}

void PlotCPM::PlotStressTensor(NOD *pn, VOX* pv){


BeginScene();

    // draw FA circles
    for (int nx=0;nx<par.NVX-1;nx++) { 
        for (int ny=0;ny<par.NVY-1;ny++) {

		int n = nx+ny*par.NVX;

            int x1= (par.PIXPERVOX)*nx+par.PIXPERVOX/2;
            int y1= (par.PIXPERVOX)*ny+par.PIXPERVOX/2;


    	double estrains[3],L1,L2,v1[2],v2[2];
	double estress[3]; 

	    QPoint Qp=QPoint(x1,y1);
      		get_estrains(pn,n,estrains);
		get_estress(n, estrains,estress);


        	L1=L2=.0; get_princs(estress,&L1,&L2,v1,v2,0);


		double radius1=5+L1/500;
		double radius2=5+L2/500;

		if(radius1<0){radius1=0;}
		if(radius2<0){radius2=0;}

		double angle =0;
		

		if(L1>=L2){angle = -atan2(v1[1],v1[0])*180/3.14;}
		if(L2>L1){angle = -atan2(v2[1],v2[0])*180/3.14;}

		



	    picture->setBrush(Qt::black);
picture->save();
picture->translate(Qp);
picture->rotate(-angle);
picture->drawEllipse(QPoint(0, 0), (int)radius1,(int)radius2);
picture->restore();

if(pv[n].ctag)
{
double etractionstress[2];
get_etractionstress(pn,n, etractionstress);
double strx=etractionstress[0]; double stry=etractionstress[1];
double tractionstress=sqrt(strx*strx + stry*stry);
double angle2 = -atan2(stry,strx)*180/3.14;
picture->setBrush(Qt::red);
picture->save();
picture->translate(Qp);
picture->rotate(-angle2);
picture->drawEllipse(QPoint(0, 0), int(2+tractionstress/500),2);
picture->restore();

}





}
}


    EndScene();


}







 
