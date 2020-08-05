// file: functions.h

#include "myparameters.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "def.h"
#include "structures.h"



extern Parameter par;

// init.c
VOX*		init_voxels(void);
NOD*		init_nodes(void);
int         	init_cells(VOX* pv);
void 		set_forces(NOD* pn,int incr);
void		set_restrictions(NOD* pn);

// cellmoves.c
void 		CPM_moves(VOX* pv, NOD* pn, int* csize, int *csumx, int *csumy, int incr, double *FA, int NRc,int* celltypes,double* sumFA,NOD* pnold);
BOOL 		splitcheckCCR(VOX* pv,  int* csize, int xt, int ttag);
void		CalcPerimeters(VOX* pv, int* cper, int NRc);
void		CalcLengths(VOX* pv,NOD* pn,double* clength,double* ecc,double* cangle,int NRc,int* csize,int* csumx, int* csumy);
int		check_contact(VOX* pv);


// CPM_dH.c
double 		calcdH(VOX* pv, NOD* pn, int* csize, int xt, int xs, int pick, int ttag, int stag, int incr,bool matrix,int NRc,int* celltypes);
double 		calcdHcontact(VOX* pv, int xt, int ttag, int stag,int* celltypes);
double 		contactenergy(int tag1, int tag2,int* celltypes);
double 		calcdHvol(int* csize, int ttag, int stag,int* celltypes);
double          calcdHadh(NOD* pn, int xt, int xs, int pick, int ttag, int stag,int* csize,VOX* pv);

// cellforces.c
BOOL 		CheckCellNodeConnection(VOX* pv, NOD* pn, int n1, int n2, int c);
void 		cell_forces(VOX* pv, NOD* pn, int* csize, int NRc,int* csumx,int* csumy,double* FA, NOD* pnold);


// FE_local.c
double** 	set_klocal(void);
void 		material_matrix(double *pD);
void 		set_matrix_B(double *pB, double x, double y);
void 		set_matrix_U(double *pU, double x, double y);
void 		get_estrains(NOD* pn, int e, double* estrains);
void 		get_deform(NOD* pn, int e , double* deform);
void 		get_estress(int e, double* estrains, double* estress);
void 		get_princs(double* str, double* L1, double* L2, double* v1, double* v2, BOOL strain);
void 		get_etractionstress(NOD* pn,int e, double* etractionstress);


// FE_assembly.c
//int** 		set_topology(void);
void		assembly(int* kcol, double* kval, double** klocal, VOX* pv);
int 		arrange_dofpos(int* dofpos, NOD* pn);
void 		reduce_K(int* kcol, double* kval, int* dofpos);

// FE_nodes2dofs.c
void 		disp_to_nodes(NOD* pn, double* u);
void 		set_disp_of_prev_incr(NOD* pn, double* u);
void 		place_node_forces_in_f(NOD* pn, double* f);
void 		kval_to_nodes(double* Km, double* kval,NOD* pn);

// FE_solver.c
void 		calc_Kdotx(int* kcol, double* kval, double* diag, double* x, double* b, int nrrdof);
void 		solvePCG(int* kcol, double* kval, double* u, double* f, int nrrdof);

// write.c
//RENE
void   		write_increment(int increment);
void 		write_pstrain(VOX* pv, NOD* pn, int increment);
void 		write_pstress(VOX* pv, NOD* pn, int increment);
void 		write_forces(NOD* pn, int increment);
void 		write_disps(NOD* pn, int increment);
//LISANNE
void 		write_cells(VOX* pv,int increment,string sigdir);
void 		write_ratiopa(int* cper,int* csize,int NRc,int increment,FILE* ofp);
void 		write_length(double* length,int NRc,int increment,FILE* ofp);
void		write_area(int* area, int NRC, int increment, FILE* ofp);
void		write_sqdis(double* sqdis, int NRC, int increment, FILE* ofp);
void 		write_eccentricity(double* ecc,int NRc,int increment,FILE* ofp);
void 		write_cangle(double* cangle,int NRc,int increment,FILE* ofp);
void 		write_twocellcontact(int contact,int increment,FILE* ofp);
void 		write_totshape(int* totshape,int increment,string totshapedir);
void		write_fa(double* FA,int increment,string fadir);
// read.c
int   		read_increment(void);
int 		read_cells(VOX* pv, int increment,string sreadcells);

// mylib.c
void 		myitostr(int n, char s[]);
void 		myreverse(char s[]);
unsigned 	mystrlen(const char *s);

// mt.c
void 		mt_init(void);
unsigned long   mt_random();


