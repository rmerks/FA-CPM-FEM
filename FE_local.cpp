/* file: femlocal.c */
#include "functions.h"

////////////////////////////////////////////////////////////////////////////////
double** set_klocal(void)
{
	double **klocal;
	int i,j,m,n;

	// two-point Gaussian integration
	// local coordinates of the integration points (1/sqrt(3) = 0.57735027)
	double intgrx[4] = {-.57735,  .57735,  .57735, -.57735};
	double intgry[4] = {-.57735, -.57735,  .57735,  .57735};

	double *pD;			// pntr to material matrix
	double D[3][3];		// material matrix
	double *pB;			// pntr to strain displacement matrix
	double B[3][8];		// strain displacement matrix
	double Bt[8][3];	// transpose of strain displacement matrix
	double BD[8][3];	// Bt * D
	double BDB[8][8];	// BD * B
	double dV;			// volume belonging to integration point

	// node positions in local coordinate system
	// double nx[4] = {-1,  1,  1, -1};
	// double ny[4] = {-1, -1,  1,  1};

	// allocate memory for klocal
	klocal = new double*[8];
	for(m=0;m<8;m++)
		klocal[m] = new double[8];

	// set matrix k to zeros
	for(m=0;m<8;m++)
	for(n=0;n<8;n++)
  		klocal[m][n] = 0.0;

	// calculate stifness matrix of the material (linear elastic isotropic)
	pD = &D[0][0];
	material_matrix(pD);
	pB = &B[0][0];

	// Determine local stiffness matrix
	// by integration. Implemented as summation over all integration points
	for(i=0; i<4; i++)      // for all integration points
	{
		// calculate matrix B in intgr pnt i
		set_matrix_B(pB, intgrx[i], intgry[i]);

		// Bt is the transpose of B
		for(m=0; m<8; m++)
		for(n=0; n<3; n++)
         	Bt[m][n] = B[n][m];

		// BD  =  Bt * D
		for(m=0; m<8; m++)
		for(n=0; n<3; n++)
		{
			BD[m][n] = 0;
			for(j=0; j<3; j++)
            	BD[m][n] += Bt[m][j] * D[j][n];
		}

		// BDB  =  BD * B
		for(m=0; m<8; m++)
		for(n=0; n<8; n++)
		{
			BDB[m][n] = 0;
			for(j=0; j<3; j++)
				BDB[m][n] += BD[m][j] * B[j][n];
		}

		// integration over the volume. This leads to adding to local
		// stifness matrix for each integration point
		// dV = dx*dy*dz = det(J) * dr*ds*dt
		// for cubic voxel elements the volume represented by one
		// integration point is equal to dV = (.5*(par.VOXSIZE))^3
		dV = .25 * (par.VOXSIZE) * (par.VOXSIZE);
		for(m=0; m<8; m++)
		for(n=0; n<8; n++)
			klocal[m][n] += BDB[m][n] * dV;
	} // endfor all integration points

cout << "klocal[0][0] " << klocal[0][0] << endl;

	return klocal;


}

////////////////////////////////////////////////////////////////////////////////
void material_matrix(double *pD)
{
	// pD is pntr to stiffness matrix
	int m, n;
	double Es;
	BOOL planestress;

	planestress = TRUE;

	for(m=0; m<3; m++)
	for(n=0; n<3; n++)
		*(pD + m +3*n) = 0;

	if(planestress)
	{
		Es = (par.YOUNGS)/(1-(par.POISSON)*(par.POISSON));
		// fill material matrix
		*(pD+0+3*0) = Es * 1;
		*(pD+1+3*1) = Es * 1;
		*(pD+0+3*1) = Es * (par.POISSON);
		*(pD+1+3*0) = Es * (par.POISSON);
		*(pD+2+3*2) = Es * .5*(1-(par.POISSON));
	}
	else // planestrain
	{
		Es = (par.YOUNGS)/((1+(par.POISSON))*(1-2*(par.POISSON)));
		// fill material matrix
		*(pD+0+3*0) = Es * (1-(par.POISSON));
		*(pD+1+3*1) = Es * (1-(par.POISSON));
		*(pD+0+3*1) = Es * (par.POISSON);
		*(pD+1+3*0) = Es * (par.POISSON);
		*(pD+2+3*2) = Es * .5*(1-2*(par.POISSON));
	}
}

////////////////////////////////////////////////////////////////////////////////
void set_matrix_B(double *pB, double r, double s)
{
	int i;
	// r,s are the local coordinates in the isoparametric element
	// constants in shape functions for node 0 to 7
	double kr[4] = {-1,  1,  1, -1};
	double ks[4] = {-1, -1,  1,  1};
	// shape function and derivatives in point r,s
	double dNdx[4];
	double dNdy[4];

	for(i=0; i<4; i++)   // for all nodes
	{
		// value of the shape function belonging to node i
		// N[i] = .25 * (1 + kx[i] * localx) * (1 + ky[i] * localy);
		// dNdr     dxdr dydr      dNdx      dxdr   0        dNdx
		// dNds  =  dxds dyds   *  dNdy  =     0  dyds    *  dNdy

		// rewriting gives the values of the derivatives of the shape function
		dNdx[i] = (2/(par.VOXSIZE)) * .25 * kr[i] * (1 + ks[i] * s);
		dNdy[i] = (2/(par.VOXSIZE)) * .25 * ks[i] * (1 + kr[i] * r);

		// calculate strain displacement matrix B
		*(pB +    0 +  (2*i)) 	 = dNdx[i];
		*(pB +    0 + ((2*i)+1)) =       0;
		*(pB +    8 +  (2*i)) 	 =       0;
		*(pB +    8 + ((2*i)+1)) = dNdy[i];
		*(pB + 2* 8 +  (2*i)) 	 = dNdy[i];
		*(pB + 2* 8 + ((2*i)+1)) = dNdx[i];
	} // endfor all nodes
}

////////////////////////////////////////////////////////////////////////////////
void set_matrix_U(double *pU, double r, double s)
{
	int i;
	// r,s are the local coordinates in the isoparametric element
	// constants in shape functions for node 0 to 7
	double kr[4] = {-1,  1,  1, -1};
	double ks[4] = {-1, -1,  1,  1};
	// shape function in point r,s
	double N[4];


	for(i=0; i<4; i++)   // for all nodes
	{
	N[i] = double(.25 * (1 + kr[i] * r) * (1 + ks[i] * s));

		// calculate displacement matrix U
		*(pU +    0 +  (2*i)) 	 = N[i];
		*(pU +    0 + ((2*i)+1)) =       0;
		*(pU +    8 +  (2*i)) 	 =       0;
		*(pU +    8 + ((2*i)+1)) = N[i];

	} // endfor all nodes


}


////////////////////////////////////////////////////////////////////////////////
void get_estrains(NOD* pn, int e, double* estrains)
{
	int i, j;
	int vx,vy, n00,n10,n11,n01;
	double *pB;
	double B[3][8];
	double u[8];

	pB = &B[0][0];   set_matrix_B(pB, 0, 0);

	vy = e/(par.NVX); vx = e%(par.NVX);
	// determine corner node numbers of this voxel
	n00 = (vx  ) + (vy  )*NNX;
	n10 = (vx+1) + (vy  )*NNX;
	n11 = (vx+1) + (vy+1)*NNX;
	n01 = (vx  ) + (vy+1)*NNX;

	u[0] = pn[n00].ux;
	u[1] = pn[n00].uy;
	u[2] = pn[n10].ux;
	u[3] = pn[n10].uy;
	u[4] = pn[n11].ux;
	u[5] = pn[n11].uy;
	u[6] = pn[n01].ux;
	u[7] = pn[n01].uy;

	for(i=0;i<3;i++)
		estrains[i]=0;

	// strain displacement relation
	for(i=0;i<3;i++)
	for(j=0;j<8;j++)
		estrains[i] += B[i][j] * u[j];
	
}

void get_deform(NOD* pn, int e, double* deform)
{
	int i, j;
	int vx,vy, n00,n10,n11,n01;
	double *pU;
	double U[2][8];
	double u[8];

	pU = &U[0][0];   set_matrix_U(pU, 0, 0);

	vy = e/(par.NVX); vx = e%(par.NVX);
	// determine corner node numbers of this voxel
	n00 = (vx  ) + (vy  )*NNX;
	n10 = (vx+1) + (vy  )*NNX;
	n11 = (vx+1) + (vy+1)*NNX;
	n01 = (vx  ) + (vy+1)*NNX;

	u[0] = pn[n00].ux;
	u[1] = pn[n00].uy;
	u[2] = pn[n10].ux;
	u[3] = pn[n10].uy;
	u[4] = pn[n11].ux;
	u[5] = pn[n11].uy;
	u[6] = pn[n01].ux;
	u[7] = pn[n01].uy;

	for(i=0;i<2;i++)
		deform[i]=0;



	// strain displacement relation
	for(i=0;i<2;i++)
	for(j=0;j<8;j++)
		deform[i] += U[i][j] * u[j]; 

	//deform[0]=0.25*(u[0]+u[2]+u[4]+u[6]);
	//deform[1]=0.25*(u[1]+u[3]+u[5]+u[7]);


	
}


////////////////////////////////////////////////////////////////////////////////
void get_estress(int e, double* estrains, double* estress)
{
	int i, j;
	double *pD;
	double D[3][3];

	pD = &D[0][0];   material_matrix(pD);
	// if stiffness changes per element I need to modify this here

	for(i=0;i<3;i++)
		estress[i]=0;

	// stress strain relation
	for(i=0;i<3;i++)
	for(j=0;j<3;j++)
		estress[i] += D[i][j] * estrains[j];
}



////////////////////////////////////////////////////////////////////////////////
void get_princs(double* str, double* pL1, double* pL2, double* v1, double* v2, BOOL strain)
{
	double xx,yy,xy;
	double T,D,T2D,sqT2D,Q,R;


	xx=str[0]; yy=str[1]; if(strain) {xy=.5*str[2];} else {xy=str[2];}
	// mat = xx xy
	//       xy yy

	if(xy==0)
	{
    	*pL1 = xx; v1[0]=1; v1[1]=0;
	    *pL2 = yy; v2[0]=0; v2[1]=1;
	}
	else
	{
		T = xx+yy; // trace
		D = xx*yy-xy*xy; // determinant

		T2D = T*T/4-D;

		if(T2D<=0)// can occur if strain very close to isotropic
		{
			*pL1 = xx; v1[0]=1; v1[1]=0;
	    	*pL2 = yy; v2[0]=0; v2[1]=1;
		}
		else
		{
			sqT2D = sqrt(T2D);
			*pL1 = T/2+sqT2D;
    		*pL2 = T/2-sqT2D;

			// eigenvector v must satisfy: mat*v = L*v
    		// xx*v[0]+xy*v[1] = L*v[0] -> if v[0]=1, v[1]=(L-xx)/xy (=Q)
    		// ||v|| = sqrt(1+Q*Q) (=R)   -> v = [1/R;Q/R]
			Q=((*pL1)-xx)/xy; R=sqrt(1+Q*Q); v1[0]=1/R; v1[1]=Q/R;
			Q=((*pL2)-xx)/xy; R=sqrt(1+Q*Q); v2[0]=1/R; v2[1]=Q/R;
		}
	}

}


////////////////////////////////////////////////////////////////////////////////
void get_etractionstress(NOD* pn,int e, double* etractionstress)
{


	int ex=e%par.NVX;			
	int ey=e/par.NVX;
	if(ex<par.NVX-1 & ex>0 & ey<par.NVY-1 & ey >0){		

	//four surrounding nodes of pixel
	int n00 = (ex  ) + (ey  )*NNX;
	int n10 = (ex+1) + (ey  )*NNX;
	int n11 = (ex+1) + (ey+1)*NNX;
	int n01 = (ex  ) + (ey+1)*NNX;
	etractionstress[0]=(pn[n00].fx+pn[n10].fx+pn[n11].fx+pn[n01].fx)/4;
	etractionstress[1]=(pn[n00].fy+pn[n10].fy+pn[n11].fy+pn[n01].fy)/4;
	etractionstress[0]=par.FORCESCALE*par.THICKNESS*etractionstress[0]/(par.VOXSIZE*par.VOXSIZE);
	etractionstress[1]=par.FORCESCALE*par.THICKNESS*etractionstress[1]/(par.VOXSIZE*par.VOXSIZE);
	

	}
		
				

}


/*
////////////////////////////////////////////////////////////////////////////////
double get_sed(int vx, int vy, NOD* pn)
{
	int i, j;
	int n00, n10, n11, n01;
	double *pD;
	double D[3][3];
	double *pB;
	double B[3][8];
	double strain[3];
	double stress[3];
	double u[8];
	double r3,sed;

	pB = &B[0][0];   set_matrix_B(pB, 0, 0);
	pD = &D[0][0];   material_matrix(pD);

	// determine corner node numbers of this voxel
	n00 = (vx  ) + (vy  )*NNX;
	n10 = (vx+1) + (vy  )*NNX;
	n11 = (vx+1) + (vy+1)*NNX;
	n01 = (vx  ) + (vy+1)*NNX;

	u[0] = pn[n00].ux;
	u[1] = pn[n00].uy;
	u[2] = pn[n10].ux;
	u[3] = pn[n10].uy;
	u[4] = pn[n11].ux;
	u[5] = pn[n11].uy;
	u[6] = pn[n01].ux;
	u[7] = pn[n01].uy;

	for(i=0;i<3;i++){ strain[i] = 0; stress[i] = 0;}

	// strain displacement relation
	for(i=0;i<3;i++)
	for(j=0;j<8;j++)
		strain[i] += B[i][j] * u[j];

	// stress strain relation
	for(i=0;i<3;i++)
	for(j=0;j<3;j++)
		stress[i] += D[i][j] * strain[j];

	// determine sed
	sed = 0;
	for(i=0;i<3;i++)
		sed += .5 * stress[i]*strain[i];
	return sed;
}
*/
