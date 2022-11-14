
//  ProblemConfig.h
//  ErrorEstimation
//
//  Created by Philippe Devloo on 09/05/18.
//

#ifndef ProblemConfig_h
#define ProblemConfig_h

#include <set>
#include "TPZAnalyticSolution.h"

/// class to guide the error estimator
struct ProblemConfig
{
    /// geometric mesh on which the computational meshes are based
    TPZGeoMesh *gmesh = 0;
    /// polynomial order of the original mesh
    int k = -1;
    int n = -1;
    
    /// option to compute the error based on continuous pressures or not
    bool makepressurecontinuous = true;
    
    /// number of uniform refinements applied to the mesh
    int ndivisions = -1;
    int adaptivityStep = -1;
    int dimension = 0;
    bool prefine = false;
    bool steklovexample = false;
    bool GalvisExample = false;
    bool TensorNonConst = false;
    bool MeshNonConvex = false;
    
    STATE alpha=1;
    STATE Km = 0.;
    STATE coefG = 0.;
    /// directory where the files will be stored
    std::string dir_name = ".";
    /// name identifying the problem
    std::string problemname;
    /// set of materialids in the mesh
    std::set<int> materialids;
    /// set of boundary condition material ids
    std::set<int> bcmaterialids;
    /// exact solution

    TPZAutoPointer<TLaplaceExample1> exact;

    
    ProblemConfig() = default;

    ProblemConfig(const ProblemConfig &cp) = default;

    ProblemConfig &operator=(const ProblemConfig &cp) = default;
};

struct MultithreadData{

    int nThreads;
    int maxThreads, maxRef;
    unsigned long int assembleTime;
    unsigned long int solveTime;

};

struct Statistics{
    int nLoops;
    int iterNum;
    TPZFMatrix<double> timeVec;
    TPZVec<double> avg, spu, cvar, serialTime;
    std::vector<std::ofstream*> csv, txt;
};

struct RSimulation {
    int nThreads, maxThreads, iterNum, maxRef, nRef, nDof;
    double assembleTime, solveTime, totalTime, L2Error,energyError;
};

struct Target{
    bool automated = false;
    bool timeEfficiency = false;
    bool errorMeasurement = false;
};

struct PreConfig{
    std::ofstream Erro, timer;
    std::ofstream *speedUpOfstream, *speeedUpPath;

    TPZVec<REAL> *rate, *Log;
    std::vector<double> errorRate, errorLog;
    int refLevel = -1;

    int k = 1;
    int n = 1;
    int dim = 1;
    int topologyMode = -1;

    std::string problem;
    std::string approx;
    std::string topology;           //Topology' name typed as input
    std::string topologyFileName;   //Simplified name used for naming files/directories

    REAL perm_Q1 = 5;      /// Permeability coefficient of even quadrants (Steklov only)
    REAL perm_Q2 = 1;

    REAL hLog = -1, h = -1000, href0 = 1;
    int numErrors = 4;

    std::string plotfile;
    std::string speedUpFilePath;
    std::string automatedFileName;
    std::string automatedFilePath;
    
    int mode = -1;  // 0 = "H1"; 1 = "Hybrid"; 2 = "Mixed";
    int argc = 1;
    int type= -1;  // 0 = ESinSin; 1 = EArcTan; ...

    bool makeScript = false;
    bool isSimplify = false;
    bool postProcess = true;
    int exp = 2; // Initial exponent of mesh refinement (numElem = 2*2^exp)
    
    bool shouldColor = true;
    bool isTBB = true;
    
    MultithreadData tData;
    Statistics stat;
    RSimulation rSimulation;
    Target target;
    
    int ref2D = -1;
    int ref3D = -1;
    int ref0 = -1;
};

#endif /* ProblemConfig_h */
