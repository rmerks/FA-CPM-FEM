/*
 *
 *  $Id$
 *
 *  This file is part of the Virtual Leaf.
 *
 *  VirtualLeaf is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  VirtualLeaf is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with the Virtual Leaf.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2010 Roeland Merks.
 *
 */

// WARNING: This file is automatically generated by make_parameter_source.pl. Do not edit.
// Do not edit. All edits will be discarded.


#ifndef _PARAMETER_H_
#define _PARAMETER_H_
#include <iostream>
#include <vector>

using namespace std;

 class Parameter {
		
 public: 
   Parameter();
   ~Parameter();
   void CleanUp(void);
   void Read(const char *filename);
   void Write(ostream &os) const;
   //void XMLAdd(xmlNode *root) const;
   //void XMLRead(xmlNode *root);
   void AssignValToPar(const char *namec, const char *valc);
   void AssignValArrayToPar(const char *namec, vector<double> valarray);
  int SEED;
  int NVX;
  int NVY;
  int MCS;
  bool INSERTMEDIUM;
  double VOXSIZE;
  int RELAXTIME;
  int NRINC;
  int MAXNRITER;
  double ACCURACY;
  double YOUNGS;
  double POISSON;
  double THICKNESS;
  double VISC;
  bool GLOBALSTRAIN;
  bool CYCLIC;
  int PERIOD;
  double LOADANGLE;
  double LOAD;
  double MOTILITY;
  bool CLASSICCPM;
  double TARGETVOLUME;
  double INELASTICITY;
  double INELASTICITY2;
  bool LEMMONROMER;
  bool NODECONNECTION;
  double LRTENSION;
  int NBHRAD;
  bool TWOCELLTYPES;
  double NOSTICKJCM;
  double NOSTICKJCC;
  double NOSTICKJCM2;
  double NOSTICKJCC2;
  double NOSTICKJCCb;
  double LAMBDAADHESION;
  double MAXAREA;
  double LAMBDAFA;
  double FAH;
  double LAMBDAPLAQUE;
  double CONFSTRESS;
  bool FORCEFA;
  double LAMBDAFORCEFA;
  bool ACTIN;
  double LAMBDAACTIN;
  double COLLAGEN;
  double COLLAGENSPEED;
  double PDEdt;
  double PDEREPEAT;
  double CAPACITYFA;
  double BASEFA;
  double GROWTHFA;
  double CATCHTENSION;
  double SLIPTENSION;
  double LOGISTICPAR;
  int PIXPERVOX;
  int LINEWIDTH;
  int STRIDE;
  int WSTRIDE;
  bool COLLAGENFIELD;
  bool STRAINFIELD;
  bool STRESSFIELD;
  bool HYDSTRESS;
  bool STRESSTENSOR;
  bool TRACTIONSTRESSFIELD;
  bool FORCEFIELD;
  bool TENSIONFIELD;
  double MAXFORCE;
  bool DEFORMFIELD;
  double MAXDEFORM;
  bool FAFIELD;
  bool FACOLOUR;
  double MAXFA;
  bool PRINCFIELD;
  bool CELLCOLOUR;
  bool ONECELL;
  bool TWOCELL;
  int DISTWOCELLS;
  bool CELLCOL;
  bool READCELLS;
  double CELLDENSITY;
  int CELLDIS;
  int BOUNDARYDIS;
  int FORBIDDENZONE;
  bool PATTERN;
  int PATTERNC;
  bool DUROTAXIS;
  double GRADIENT;
  int NRcf;
  bool COLORBAR;
  double MAXCOLORBARSTRAIN;
  double MAXCOLORBARSTRESS;
  double MAXCOLORBARDENS;
  int WIDTHCOLORBAR;
  bool WRATIOPA;
  bool WLENGTH;
  bool WAREA;
  bool WSQDIS;
  bool WECC;
  bool WANGLE;
  bool WSIGMA;
  bool WTOTSHAPE;
  bool WFA;
  bool WTWOCELLCONTACT;
  bool WTWOCELLANGLECM;
 private:
 };

 ostream &operator<<(ostream &os, Parameter &p);
 const char *sbool(const bool &p);


#endif
