#ifndef PLOTCPM_H
#define PLOTCPM_H
#include <QString>
#include "qtgraph.h"
#include "structures.h"
#include "functions.h"
#include "myparameters.h"

class PlotCPM : public QtGraphics
{
    Q_OBJECT
public:
    explicit PlotCPM(int xfield, int yfield, char *moviefile=0);
    void Plot(VOX* pv, int* celltypes);
    void PlotHueStrainField(bool STRAINFIELD);
    void PlotHueChemField(double *chem);
    void PlotHueStressField(bool STRESSFIELD);
    void PlotHueStressField_CTB(bool STRESSFIELD);
    void PlotPrincipleStrainField(VOX *pv);
    double MaxStrainMagnitude(void) const;
    double MaxForceMagnitude(void) const;
    double MaxTensionMagnitude(void) const;
    void PlotPrincipleStressField(VOX *pv);
    double MaxStressMagnitude(void) const;
    void CalculateStrainField(NOD *pn) const;
    void CalculateForceField(NOD *pn) const;
    void CalculateStressField(NOD *pn) const;
    void PlotNodalForces(NOD *pn,VOX* pv);
    void PlotNodalTension(NOD *pn);
    void PlotNodalTension2(NOD *pn);
    void PlotNodalDeform(NOD *pn);
    void PlotNodalFA(NOD *pn, double *FA);
    void PlotNodalFA2(VOX *pv, double *FA);
    void PlotStressTensor(NOD *pn,VOX* pv);
    void PlotNodeConnection(NOD *pn, int *pathx, int *pathy);
    void StrainColorBar(void);
    void ChemColorBar(double *chem);
    void StressColorBar(void);
    void StressColorBar_CTB(void);
    void DrawColorBarLabel(int yval, double val);
    void DeleteStrainForce();
signals:
    
public slots:
    
private:
    double ***strain;
    double ***force;
    double ***tension;
    double ***stress;
};

#endif // PLOTCPM_H
