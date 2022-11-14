//
// Created by victor on 02/09/2020.
//

#include "Solver.h"
#include <TPZMultiphysicsCompMesh.h>
#include "TPZLinearAnalysis.h"
#include "DataStructure.h"
#include "MeshInit.h"
#include "TPZCompMeshTools.h"
#include "TPZCreateMultiphysicsSpace.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "Tools.h"
#include "Output.h"
#include "pzvisualmatrix.h"
#include "MeshInit.h"
#include "TPZTimer.h"
#include <chrono>
#include "TPZFrontSym.h"
#include "TPZStructMatrixOMPorTBB.h"
#include "InputTreatment.h"
#include <cstdlib>//to perturb uh
#include <ctime>// to perturb uh

#ifdef FEMCOMPARISON_USING_MKL
#include "mkl.h"
#endif

//#include <tbb/parallel_for.h>
//#include <tbb/task_scheduler_init.h>
//#include "omp.h"

#ifdef PZ_LOG
static TPZLogger loggerST("solveTime");
static TPZLogger loggerAT("assembleTime");
#endif


void Solve(ProblemConfig &config, PreConfig &preConfig){

#ifdef FEMCOMPARISON_TIMER
    extern double solveTime;
    extern int nTestsSolve;
#endif
    TPZCompMesh *cmesh = InsertCMeshH1(config,preConfig);
    //std::ofstream outTXT("ZH1_insert_cmeshH1.txt");
    //cmesh->Print(outTXT);
    
    TPZMultiphysicsCompMesh *multiCmesh = new TPZMultiphysicsCompMesh(config.gmesh);
    
    int interfaceMatID = -10;
    int hybridLevel = 1;

    const clock_t start = clock();
    //std::ofstream outTXT2("ZCreatedCondensedElements.txt");

    switch(preConfig.mode){
        case 0: //H1
            //TPZCompMeshTools::CreatedCondensedElements(cmesh, false, false);
            //cmesh->Print(outTXT2);
            SolveH1Problem(cmesh, config, preConfig);
            break;
        case 1: //Hybrid
            CreateHybridH1ComputationalMesh(multiCmesh, interfaceMatID,preConfig, config,hybridLevel);
            SolveHybridH1Problem(multiCmesh, interfaceMatID, config, preConfig);
            break;
        case 2: //Mixed
            CreateMixedComputationalMesh(multiCmesh, preConfig, config);
            {
#ifdef PZ_LOG
                TPZTimer timer;
                if(loggerST.isDebugEnabled())
                timer.start();
#endif
            SolveMixedProblem(multiCmesh, config, preConfig);
            }
            break;
            default:
            DebugStop();
            break;
    }
    FlushTime(preConfig,start);

    if(preConfig.postProcess) DrawCompMesh(config,preConfig,cmesh,multiCmesh);
}

void DrawMesh(ProblemConfig &config, PreConfig &preConfig, TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *multiCmesh) {

    std::stringstream ref;
    ref << "_ref-" << 1/preConfig.h <<" x " << 1/preConfig.h;
    std::string refinement =  ref.str();

    std::ofstream out(preConfig.plotfile + "/gmesh"+ refinement + ".vtk");
    std::ofstream out2(preConfig.plotfile + "/gmesh"+ refinement + "txt");
    std::ofstream out3(preConfig.plotfile + "/cmesh.txt");

    TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
    config.gmesh->Print(out2);

    if (preConfig.mode == 0) cmesh->Print(out3);
    else multiCmesh->Print(out3);
}

void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_Mixed, PreConfig &pConfig, ProblemConfig &config){
    InsertMaterialMixed(cmesh_Mixed, config,pConfig);
    TPZManVector<TPZCompMesh *, 4> meshvector(4);
    CreateMixedAtomicMeshes(meshvector,pConfig, config);

    TPZManVector<int> active(4, 1);

    cmesh_Mixed->BuildMultiphysicsSpace(active,meshvector);
    CreateCondensedMixedElements(cmesh_Mixed);
    cmesh_Mixed->LoadReferences();
    cmesh_Mixed->InitializeBlock();
}
void CreateCondensedMixedElements(TPZMultiphysicsCompMesh *cmesh_Mixed){

    int numActiveSpaces = cmesh_Mixed->MeshVector().size();
    if(numActiveSpaces == 4){
        int64_t nconnects = cmesh_Mixed->NConnects();
        for (int64_t ic = 0; ic<nconnects; ic++) {
            TPZConnect &c = cmesh_Mixed->ConnectVec()[ic];
            if(c.LagrangeMultiplier() == 4) c.IncrementElConnected();
        }
    }
    bool keeponelagrangian = true, keepmatrix = false;
    TPZCompMeshTools::CreatedCondensedElements(cmesh_Mixed, keeponelagrangian, keepmatrix);
}

void CreateHybridH1ComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int &interFaceMatID , PreConfig &pConfig, ProblemConfig &config,int hybridLevel){
    auto spaceType = TPZCreateMultiphysicsSpace::EH1Hybrid;
    if(hybridLevel == 2) {
        spaceType = TPZCreateMultiphysicsSpace::EH1HybridSquared;
    }
    else if(hybridLevel != 1) {
        DebugStop();
    }

    TPZCreateMultiphysicsSpace createspace(config.gmesh, spaceType);
    //TPZCreateMultiphysicsSpace createspace(config.gmesh);
    std::cout << cmesh_H1Hybrid->NEquations();

    createspace.SetMaterialIds({1,2,3}, {-6,-5,-2,-1});
    createspace.fH1Hybrid.fHybridizeBCLevel = 1;//opcao de hibridizar o contorno
    createspace.ComputePeriferalMaterialIds();

    TPZManVector<TPZCompMesh *> meshvec;

    int pOrder = config.n+config.k;
    createspace.CreateAtomicMeshes(meshvec,pOrder,config.k);

    InsertMaterialHybrid(cmesh_H1Hybrid, config,pConfig);
    createspace.InsertPeriferalMaterialObjects(cmesh_H1Hybrid);
    
    {
    std::vector<std::ofstream*> output;
    std::string foldername="PrintMesh";
    std::string command = "mkdir -p " + foldername;
    system(command.c_str());
    for(int ii=0;ii<5;ii++){
        std::stringstream name;
        name << foldername <<"/meshB4" << ii << ".txt";
        std::string oname = name.str();
        std::ofstream* x = new std::ofstream(oname);
        output.push_back(x);
        if(ii==0){
            cmesh_H1Hybrid->Print(*output[ii]);
        }else{
            meshvec[ii-1]->Print(*output[ii]);
            //cmesh_H1Hybrid->MeshVector()[ii-1]->Print(*output[ii]);
        }
    }
    }
    
    cmesh_H1Hybrid->BuildMultiphysicsSpace(meshvec);
    createspace.InsertLagranceMaterialObjects(cmesh_H1Hybrid);

    createspace.AddInterfaceElements(cmesh_H1Hybrid);
    createspace.GroupandCondenseElements(cmesh_H1Hybrid);

    cmesh_H1Hybrid->InitializeBlock();
    cmesh_H1Hybrid->ComputeNodElCon();

    interFaceMatID = createspace.fH1Hybrid.fLagrangeMatid.first;
    {
    std::vector<std::ofstream*> output;
    std::string foldername="PrintMesh";
    std::string command = "mkdir -p " + foldername;
    system(command.c_str());
    for(int ii=0;ii<5;ii++){
        std::stringstream name;
        name << foldername <<"/mesh" << ii << ".txt";
        std::string oname = name.str();
        std::ofstream* x = new std::ofstream(oname);
        output.push_back(x);
        if(ii==0){
            cmesh_H1Hybrid->Print(*output[ii]);
        }else{
            cmesh_H1Hybrid->MeshVector()[ii-1]->Print(*output[ii]);
        }
    }
    }
    
}

void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config, struct PreConfig &pConfig){

#ifndef OPTMIZE_RUN_TIME
    config.exact.operator*().fSignConvention = 1;
#endif
    
    std::cout << "Solving H1 " << std::endl;

    if(1){
        TPZFMatrix<REAL> mat(100,100);
        cmeshH1->ComputeFillIn(100, mat);
        VisualMatrix(mat, "matH1_1.vtk");
    }
    
    TPZLinearAnalysis an(cmeshH1);

    if(1){
        TPZFMatrix<REAL> mat(100,100);
        cmeshH1->ComputeFillIn(100, mat);
        VisualMatrix(mat, "matH1_2.vtk");
    }
#ifdef FEMCOMPARISON_USING_MKL
//    TPZSSpStructMatrix<> strmat(cmeshH1);
//    strmat.SetNumThreads(0);
    
    TPZSSpStructMatrix<STATE,TPZStructMatrixOMPorTBB<STATE>> strmat(cmeshH1);
    strmat.SetNumThreads(pConfig.tData.nThreads);
    auto *strmatPointer = new TPZSSpStructMatrix(strmat);
    
    auto myParInterface = dynamic_cast<TPZStructMatrixOMPorTBB<STATE>*>(strmatPointer);
    if(myParInterface){
        myParInterface->SetShouldColor(pConfig.shouldColor);
        myParInterface->SetTBBorOMP(pConfig.isTBB);
    }
    //        strmat.SetDecomposeType(ELDLt);
#else
    //TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(cmeshH1);
    
    TPZSSpStructMatrix<> strmat(cmeshH1);
    strmat.SetNumThreads(0);
    //        TPZSkylineStructMatrix strmat3(cmesh_HDiv);
    //        strmat3.SetNumThreads(8);
#endif
    
    std::set<int> matids;
    for (auto matid : config.materialids)
        matids.insert(matid);

    for(auto mat : config.bcmaterialids){
        matids.insert(mat);
    }

    strmat.SetMaterialIds(matids);
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE> *direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
    
    an.Assemble();
    an.Solve();
    
    //cmeshH1->Print(std::cout);
    
    //MARK: H1-uh perturbation
//    int64_t ncoefs = cmeshH1->Solution().Rows();
//    TPZFMatrix<STATE> &locmat = cmeshH1->Solution();
//    for (int ind=0; ind < ncoefs; ind++){
//        locmat(ind,0) += 0.;
//
//    }
    
    if(1){
        srand (time(NULL));
        int64_t ncel = cmeshH1->NElements();
        std::set<int> matids = {-1};
        for(int ic=0; ic<ncel; ic++){
            TPZCompEl *cel = cmeshH1->ElementVec()[ic];
            if(!cel){
                continue;
            }
            
            TPZGeoEl *gel = cel->Reference();
            if(!gel){
                continue;
            }
            
            bool isinterior = true;
            
            //TODO: Tour the neighbors
            TPZStack<TPZCompElSide> neighs;
            int nneighs;
            int nsides = gel->NSides();
            for(int side=0; side<nsides; side++){
                TPZGeoElSide geoelside(gel,side);
                if(geoelside.HasNeighbour(matids)){
                    isinterior = false;
                }
                
                if(!isinterior){
                    continue;
                }
                
                TPZConnect &c = cel->Connect(side);
                //int ncon = cel->NConnects();
                TPZFMatrix<STATE> &sol = cmeshH1->Solution();
                TPZBlock &matloc =  cmeshH1->Block();
                int64_t seqnum = c.SequenceNumber();
                int bsize = c.NShape() * c.NState();
                
                for (int i=0; i<bsize; i++){
                    sol(matloc.Index(seqnum, i),0) = 0.5 - (STATE) rand() / RAND_MAX;
                }
                
            }
        }

    }
    
    int64_t nelem = cmeshH1->NElements();
    cmeshH1->LoadSolution(cmeshH1->Solution());
    cmeshH1->ExpandSolution();
    cmeshH1->ElementSolution().Redim(nelem, 10);

    ////Calculo do erro
    std::cout << "Computing Error H1 " << std::endl;

#ifndef OPTMIZE_RUN_TIME
    an.SetExact(config.exact.operator*().ExactSolution());
#endif
    
    StockErrorsH1(an,cmeshH1,pConfig.Erro,pConfig.Log,pConfig);

    ////PostProcess
    if(pConfig.postProcess) {
        TPZStack<std::string> scalnames, vecnames;
        scalnames.Push("Solution");
        vecnames.Push("Derivative");
        scalnames.Push("ExactSolution");

        int dim = cmeshH1->Reference()->Dimension();

        std::string plotname;
        {
            std::stringstream out;
            out << pConfig.plotfile << "/" << config.problemname <<"_" << pConfig.topologyFileName << "_k-" << config.k << "_n-" << config.n;

            if(dim == 2) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";
            if(dim == 3) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";

            plotname = out.str();
        }
        int resolution=5;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution,dim);
    }
    std::cout << "FINISHED!" << std::endl;
}
                      
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,int InterfaceMatId, struct ProblemConfig config,struct PreConfig &pConfig){

#ifndef OPTMIZE_RUN_TIME
    config.exact.operator*().fSignConvention = 1; 
#endif
    NonConformAssemblage(cmesh_H1Hybrid,InterfaceMatId,config,pConfig,1);
}

void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed,struct ProblemConfig config,struct PreConfig &pConfig) {

#ifndef OPTMIZE_RUN_TIME
    config.exact.operator*().fSignConvention = 1;
#endif

    NonConformAssemblage(cmesh_Mixed,-1,config,pConfig,0);

}
void StockErrorsH1(TPZAnalysis &an,TPZCompMesh *cmesh, std::ofstream &Erro, TPZVec<REAL> *Log,PreConfig &pConfig){

    TPZManVector<REAL,6> Errors;
    Errors.resize(pConfig.numErrors);
    bool store_errors = false;

    an.PostProcessError(Errors, store_errors, Erro);

    std::cout << "Errors =(";
    for (int i = 0; i < Errors.size(); i++){
        std::cout << Errors[i];
        if(i<Errors.size()-1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    
//    if ((*Log)[0] != -1) {
//        for (int j = 0; j < 3; j++) {
//            (*pConfig.rate)[j] =
//                    (log10(Errors[j]) - log10((*Log)[j])) /
//                    (log10(pConfig.h) - log10(pConfig.hLog));
//            Erro << "rate " << j << ": " << (*pConfig.rate)[j] << std::endl;
//        }
//    }
//
//    Erro << "h = " << pConfig.h << std::endl;
//    Erro << "DOF = " << cmesh->NEquations() << std::endl;
//    for (int i = 0; i < pConfig.numErrors; i++)
//        (*Log)[i] = Errors[i];
//    Errors.clear();
}

void StockErrors(TPZAnalysis &an,TPZMultiphysicsCompMesh *cmesh, std::ofstream &Erro, TPZVec<REAL> *Log,PreConfig &pConfig){

    TPZManVector<REAL,6> Errors;
    Errors.resize(pConfig.numErrors);
    bool store_errors = false;

    an.PostProcessError(Errors, store_errors, Erro);
    //an.PostProcessError(Errors, store_errors);
    
    std::cout << "Errors =(";
    for (int i = 0; i < Errors.size(); i++){
        std::cout << Errors[i];
        if(i<Errors.size()-1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    
    //std::cout<<"nnnnnnnn"<<std::endl;
    //for(int i=0;i<Errors.size();i++)
        //std::cout<<Errors[i]<<std::endl;
//    if ((*Log)[0] != -1) {
//        for (int j = 0; j < 3; j++) {
//            (*pConfig.rate)[j] =
//                    (log10(Errors[j]) - log10((*Log)[j])) /
//                    (log10(pConfig.h) - log10(pConfig.hLog));
//            Erro << "rate " << j << ": " << (*pConfig.rate)[j] << std::endl;
//        }
//    }

//    Erro << "h = " << pConfig.h << std::endl;
//    Erro << "DOF = " << cmesh->NEquations() << std::endl;
//    for (int i = 0; i < pConfig.numErrors; i++)
//        (*Log)[i] = Errors[i];
//    Errors.clear();
}

void ComputePragerSynge(TPZAnalysis &an,TPZMultiphysicsCompMesh *cmesh, std::ofstream &Erro, TPZVec<REAL> *Log,PreConfig &pConfig){

    TPZManVector<REAL,3> terms;
    int nTerms = 3;
    terms.resize(nTerms);
    bool store_errors = false;

    an.PostProcessError(terms, store_errors, Erro);
    //an.PostProcessError(Errors, store_errors);
    for(int i=0; i<nTerms; i++){
        terms[i] *= terms[i];
    }
    
    std::cout << "Terms =(";
    for (int i = 0; i < terms.size(); i++){
        std::cout << terms[i];
        if(i<nTerms) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    
    std::cout << "term0+term1-term2: " << terms[0]+terms[1]-terms[2] << std::endl;
    
}


void NonConformAssemblage(TPZMultiphysicsCompMesh *multiCmesh,int InterfaceMatId, struct ProblemConfig config,struct PreConfig &pConfig, bool isHybridH1){

//    unsigned long int assembleDuration;
//    unsigned long int solveDuration;
    
    std::cout << "Solving " << pConfig.approx << " " << pConfig.topology << " ref " << pConfig.refLevel << " nThreads " << pConfig.tData.nThreads << " isColoring " << pConfig.shouldColor << " isTBB " << pConfig.isTBB << std::endl;
    {
        TPZFMatrix<REAL> mat(100,100);
        multiCmesh->ComputeFillIn(100, mat);
        VisualMatrix(mat, "arch1.vtk");
    }
    
    TPZLinearAnalysis an(multiCmesh);
    
    {
        TPZFMatrix<REAL> mat(100,100);
        multiCmesh->ComputeFillIn(100, mat);
        VisualMatrix(mat, "arch2.vtk");
    }
#ifdef FEMCOMPARISON_USING_MKL
    TPZSSpStructMatrix<STATE,TPZStructMatrixOMPorTBB<STATE>> strmat(multiCmesh);
    strmat.SetNumThreads(pConfig.tData.nThreads);
    
    auto *strmatPointer = new TPZSSpStructMatrix(strmat);

auto myParInterface = dynamic_cast<TPZStructMatrixOMPorTBB<STATE>*>(strmatPointer);
if(myParInterface){
    myParInterface->SetShouldColor(pConfig.shouldColor);
    myParInterface->SetTBBorOMP(pConfig.isTBB);
}

#else
    //TPZSkylineStructMatrix strmat(cmesh_H1Hybrid);
    
    TPZSSpStructMatrix<> strmat(multiCmesh);
    
    strmat.SetNumThreads(0);
#endif
    
    std::set<int> matIds;
    for (auto matid : config.materialids)
        matIds.insert(matid);
    for (auto matidbc : config.bcmaterialids)
        matIds.insert(matidbc);

    if (isHybridH1){
        matIds.insert(InterfaceMatId);
        strmat.SetMaterialIds(matIds);
    }
    
    an.SetStructuralMatrix(strmat);

    TPZStepSolver<STATE>* direct = new TPZStepSolver<STATE>;
    direct->SetDirect(ELDLt);
    an.SetSolver(*direct);
    delete direct;
    direct = 0;
#ifdef FEMCOMPARISON_TIMER
    auto beginAss = std::chrono::high_resolution_clock::now();

#endif
        
    an.Assemble();
    
    if(0){//MARK: Compute condition number
        double cond = an.MatrixSolver<STATE>().Matrix()->ConditionNumber(2);
        TPZFMatrix<STATE> &frhs = an.Rhs();
        frhs*=-1;
        std::cout << "CondNumber:\t"<< cond << std::endl;
    }
    
#ifdef FEMCOMPARISON_TIMER
    auto endAss = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(endAss - beginAss);
    pConfig.tData.assembleTime = static_cast<unsigned long int>(elapsed.count());
#endif
    
#ifdef FEMCOMPARISON_TIMER
    auto begin = std::chrono::high_resolution_clock::now();
#endif
        
    int effNthreads = pConfig.tData.nThreads;
    if (effNthreads == 0) effNthreads =1;
    
#ifdef FEMCOMPARISON_USING_MKL
    mkl_set_num_threads_local(effNthreads);
#endif
        
    an.Solve();
    
#ifdef FEMCOMPARISON_TIMER
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsedSolve = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        pConfig.tData.solveTime = static_cast<unsigned long int>(elapsedSolve.count());
        std::cout << "nDofs = " << multiCmesh->NEquations() << std::endl;
    
    if (pConfig.target.automated || pConfig.target.timeEfficiency){
        FlushSpeedUpResults(pConfig.tData.assembleTime, pConfig.tData.solveTime, pConfig);
        std::ofstream* dofData  = new std::ofstream;
        dofData->open(pConfig.automatedFilePath + "/dofData.csv",std::ofstream::app);
        *dofData << pConfig.refLevel << "," << multiCmesh->NEquations() << std::endl << std::flush;
    }
#endif
    
    if (pConfig.target.errorMeasurement){
#ifndef OPTMIZE_RUN_TIME
        an.SetExact(config.exact.operator*().ExactSolution());
        TPZManVector<REAL,6> Errors;
        Errors.resize(pConfig.numErrors);
        
        if(pConfig.mode == 1){
            pConfig.numErrors = 4 ;
        } else if (pConfig.mode == 2){
            pConfig.numErrors = 5;
        } else DebugStop();

        Errors.resize(pConfig.numErrors);
        
        bool store_errors = false;

        an.PostProcessError(Errors,store_errors);
        
        double L2error, energyError;
        int nDof = multiCmesh->NEquations();
        if(pConfig.mode == 1){
            L2error = Errors[0];
            energyError = Errors[3];
        } else if (pConfig.mode == 2){
            L2error = Errors[0];
            energyError = Errors[1];
        }
        
        FlushSpeedUpResults(L2error, energyError,nDof, pConfig);
#endif
    }

    if(1){
        std::cout << "Computing Error " << std::endl;
        an.SetExact(config.exact.operator*().ExactSolution());

        ////Calculo do erro
        StockErrors(an,multiCmesh,pConfig.Erro,pConfig.Log,pConfig);
        std::cout << "DOF = " << multiCmesh->NEquations() << std::endl;

    }
    
    if(pConfig.postProcess) {
        std::cout << "Computing Error " << std::endl;
#ifndef OPTMIZE_RUN_TIME
        an.SetExact(config.exact.operator*().ExactSolution());
#endif
        ////Calculo do erro
//        StockErrors(an,multiCmesh,pConfig.Erro,pConfig.Log,pConfig);
//        std::cout << "DOF = " << multiCmesh->NEquations() << std::endl;
    
        ////PostProcess
        TPZStack<std::string> scalnames, vecnames;
        
        if (isHybridH1){
            scalnames.Push("Pressure");
            scalnames.Push("ExactPressure");
            vecnames.Push("Flux");
        }else {
            scalnames.Push("Pressure");
            scalnames.Push("ExactPressure");
            vecnames.Push("Flux");
            vecnames.Push("ExactFlux");
        }
        int dim = pConfig.dim;
        std::string plotname;
        {
            std::stringstream out;
            out << pConfig.plotfile << "/" << config.problemname <<"_" << pConfig.topologyFileName << "_k-" << config.k << "_n-" << config.n;

            if(dim == 2) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";
            if(dim == 3) out  << "_numEl_" << 1/pConfig.h << " x " << 1/pConfig.h << " x " << 1/pConfig.h <<".vtk";

            plotname = out.str();
        }
        int resolution = 1;
        an.DefineGraphMesh(dim, scalnames, vecnames, plotname);
        an.PostProcess(resolution, dim);
    }
}


void VerifyPragerSynge(int argc, char *argv[],PreConfig &pConfigH1, PreConfig &pConfigMix){

    //MARK: H1
    DataInitialization(argc,argv,pConfigH1);
    if (pConfigH1.mode != 0){ //check that the FEM is H1
        DebugStop();
    }

    ProblemConfig configH1;
    Configure(configH1,pConfigH1.refLevel,pConfigH1,argv);

    TPZCompMesh *cmeshH1 = InsertCMeshH1(configH1,pConfigH1);
    SolveH1Problem(cmeshH1, configH1, pConfigH1);

    if(1){
        std::ofstream ofs1("CmeshH1.txt");
        cmeshH1->Print(ofs1);
        cmeshH1->Solution().Print("solU_H1");
        std::string command = "cp Erro.txt FEMH1_Erro.txt";
        system(command.c_str());
    }
    
    //MARK: Mixed
    DataInitialization(argc,argv,pConfigMix);
    if (pConfigMix.mode != 2) DebugStop();

    ProblemConfig configMix;
    Configure(configMix,pConfigMix.refLevel,pConfigMix,argv);
    configMix.gmesh->ResetReference();
    configMix.gmesh = configH1.gmesh;
    
    TPZMultiphysicsCompMesh *cmeshMix = new TPZMultiphysicsCompMesh(configMix.gmesh);

    int interfaceMatID = -10;
    int fluxMatID = -10;

    CreateMixedComputationalMesh(cmeshMix, pConfigMix, configMix);
    SolveMixedProblem(cmeshMix, configMix, pConfigMix);
    
    if(1){
        std::ofstream ofs2("CmeshMix.txt");
        cmeshMix->Print(ofs2);
        //cmeshMix->Solution().Print("solMix");
    }
    
    //MARK: Prager-Synge
    TPZMultiphysicsCompMesh *cmeshH1Mix = new TPZMultiphysicsCompMesh(configH1.gmesh);
    CreateH1MixCompMesh(cmeshH1,cmeshMix,cmeshH1Mix,pConfigH1,configMix);
    TPZLinearAnalysis an(cmeshH1Mix,false);
    
    if(1){
        std::ofstream ofs("CmeshH1Mix.txt");
        cmeshH1Mix->Print(ofs);
    }
    
    PreConfig pConfigH1Mix;
    ComputePragerSynge(an,cmeshH1Mix,pConfigH1Mix.Erro,pConfigH1Mix.Log,pConfigH1Mix);
    //PostProcessH1Mix(cmeshH1Mix,pConfigH1,configMix);
}

void CreateH1MixCompMesh(TPZCompMesh *cmeshH1, TPZMultiphysicsCompMesh *cmeshMix, TPZMultiphysicsCompMesh *cmeshH1Mix,PreConfig &pConfigH1, ProblemConfig &configMix){
    // initialize the post processing mesh
    //fOriginal->CopyMaterials(fPostProcMesh);

    InsertMaterialH1Mix(cmeshH1Mix,pConfigH1,configMix);
    
    TPZGeoMesh* gmesh = cmeshH1Mix->Reference();
    if(cmeshH1->Reference() != gmesh)
        DebugStop();
    if(cmeshMix->Reference() != gmesh)
        DebugStop();
        
    TPZManVector<TPZCompMesh *> mesh_vectors(3, 0);
    TPZManVector<int> active(3, 1);

    mesh_vectors[0] = cmeshH1;      // H1
    mesh_vectors[1] = cmeshMix->MeshVector()[1];      // Mixed potential
    mesh_vectors[2] = cmeshMix->MeshVector()[0];      // HDiv mixed flux

    TPZStack<int64_t> gelIndex;
    for(int64_t i=0; i<cmeshH1->NElements(); i++){
        TPZCompEl* cel = cmeshH1->Element(i);
        if(!cel) continue;
        TPZGeoEl* gel = cel->Reference();
        if(!gel) continue;
        if(gel->MaterialId() != 1) continue;
        gelIndex.Push(gel->Index());
    }
    cmeshH1Mix->BuildMultiphysicsSpace(mesh_vectors,gelIndex);
    std::ofstream ofs("CmeshH1Mix.txt");
    cmeshH1Mix->Print(ofs);
}
