#include "InputTreatment.h"
#include "MeshInit.h"
#include "Solver.h"
#include "Output.h"
#include "tpzgeoelrefpattern.h"
#include "DataStructure.h"
#include "Tools.h"
#include <tuple>
#include "pzvisualmatrix.h"
#include <TPZTimer.h>
#include "computStatist.h"
using namespace std;

#ifdef FEMCOMPARISON_TIMER
//#include <chrono>
#include <math.h>//used only for timer statistics
#include <fstream>//used only for timer statistics
#include <vector>//used only for timer statistics

#endif
int main(int argc, char *argv[]) {
#ifdef FEMCOMPARISON_TIMER
    //TPZTimer timer;
    //timer.start();
    
    bool atypical1=false;
    
    bool MKL_contribute;
#ifdef FEMCOMPARISON_USING_MKL
    MKL_contribute=true;
#endif
#endif
    
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    PreConfig pConfig;
    pConfig.k = 1;//
    pConfig.n = 2;
    pConfig.problem = "ESinSin";              //// {"ESinSin","EArcTan",ESteklovNonConst", "ESteepWave"}
    pConfig.approx = "Hybrid";                //// {"H1","Hybrid", "Mixed"}
    pConfig.topology = "Quadrilateral";       //// Triangular, Quadrilateral, Tetrahedral, Hexahedral, Prism
    pConfig.refLevel =5;                     //// How many refinements
    pConfig.postProcess = true;                  //// Print geometric and computational mesh
    pConfig.shouldColor =false;
    pConfig.isTBB = false;
    pConfig.tData.nThreads = 0;
    
    
    if(argc == 2) {
      argc = 1;
      pConfig.tData.nThreads = atoi(argv[1]);
    }
    EvaluateEntry(argc,argv,pConfig);
    InitializeOutstream(pConfig,argv);
    
    
    pConfig.exp *= pow(2,pConfig.refLevel-1);
    pConfig.h = 1./pConfig.exp;
    ProblemConfig config;
    Configure(config,pConfig.refLevel,pConfig,argv);
    Solve(config,pConfig);
    pConfig.hLog = pConfig.h;
    if(pConfig.postProcess){
        std::string command = "cp ErroHybrid.txt " + pConfig.plotfile + "/Erro.txt";
        system(command.c_str());
        FlushTable(pConfig,argv);
    }
    //timer.stop();
    
    cout<<"Number of assembly threads: "<<pConfig.tData.nThreads<<endl;
    cout<<"*********** Statistics for the assembly time *****"<<endl;
    cout<<"Time(seconds): "<<pConfig.tData.assembleTime*1E-9<<endl;
    cout<<"*********** Statistics for the linear system solve time *****"<<endl;
    cout<<"Time(seconds): "<<pConfig.tData.solveTime*1E-9<<endl;
    
    return 0;
}





