
//Runs Mixed and H1-conforming FEMs simulation to numerical verification of Prager-Synge theorem

#include "InputTreatment.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>

int main(int argc, char *argv[]) {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    /*Configures H1-conforming*/
    PreConfig pConfigH1, pConfigMix;
    pConfigH1.k = 1;
    pConfigH1.n = 0;
    pConfigH1.problem = "ESinSin";          // {"ESinSin","EArcTan",ESteklovNonConst","ESteepWave"}
    pConfigH1.approx = "H1";                // {"H1","Hybrid", "Mixed"}
    pConfigH1.topology = "Quadrilateral";   // Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfigH1.refLevel = 2;            // How many refinements
    pConfigH1.postProcess = true;         // Print geometric and computational mesh

    pConfigH1.shouldColor = false;
    pConfigH1.isTBB = false;
    pConfigH1.tData.nThreads = 0;
    
    /*Configures Mixed*/
    pConfigMix.k = 1;
    pConfigMix.n = 0;
    pConfigMix.problem = "ESinSin";          // {"ESinSin","EArcTan",ESteklovNonConst"}
    pConfigMix.approx = "Mixed";                // {"H1","Hybrid", "Mixed"}
    pConfigMix.topology = "Quadrilateral";   // Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfigMix.refLevel = pConfigH1.refLevel;            // How many refinements
    pConfigMix.postProcess = true;
    pConfigMix.shouldColor =false;
    pConfigMix.isTBB = false;
    pConfigMix.tData.nThreads = 0;
    
    pConfigH1.exp *= pow(2,pConfigH1.refLevel);
    pConfigH1.h = 1./pConfigH1.exp;
    pConfigMix.exp *= pow(2,pConfigH1.refLevel);
    pConfigMix.h = 1./pConfigH1.exp;
    
    // Solve H1-conforming and mixed problems, then compute the Prager-Synge equality terms
    VerifyPragerSynge(argc,argv,pConfigH1,pConfigMix);

    return 0;
}
