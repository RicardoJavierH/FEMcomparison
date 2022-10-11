//
// Created by victor on 02/09/2020.
//

#ifndef FEMCOMPARISON_ANALYTICS_H
#define FEMCOMPARISON_ANALYTICS_H

#include <TPZMultiphysicsCompMesh.h>
#include "TPZLinearAnalysis.h"
#include "DataStructure.h"

void calcstiffTestTBB(TPZCompMesh *cmesh, int nthreads);
void calcstiffTestOMP(TPZCompMesh *cmesh, int nthreads);


//// Call required methods to build a computational mesh for an Pryymal Hybrid approximation
void CreateHybridH1ComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, int &InterfaceMatId,PreConfig &eData, ProblemConfig &config,int hybridLevel);

//// Call required methods to build a computational mesh for a Mixed approximation
void CreateMixedComputationalMesh(TPZMultiphysicsCompMesh *cmesh_H1Mixed,PreConfig &eData, ProblemConfig &config);
void CreateCondensedMixedElements(TPZMultiphysicsCompMesh *cmesh_Mixed);

//// Solve classical H1 problem
void SolveH1Problem(TPZCompMesh *cmeshH1,struct ProblemConfig &config, struct PreConfig &eData);

//// Solve Primal Hybrid problem
void SolveHybridH1Problem(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, int InterfaceMatId, struct ProblemConfig config, struct PreConfig &eData);

//// Solve Mixed problem
void SolveMixedProblem(TPZMultiphysicsCompMesh *cmesh_Mixed,struct ProblemConfig config,struct PreConfig &eData);

//// Error Management
void StockErrorsH1(TPZAnalysis &an,TPZCompMesh *cmesh,std::ofstream &Erro, TPZVec<REAL> *Log, PreConfig &eData);

//// Error Management
void StockErrors(TPZAnalysis &an,TPZMultiphysicsCompMesh *cmesh,std::ofstream &Erro, TPZVec<REAL> *Log, PreConfig &eData);

void ComputePragerSynge(TPZAnalysis &an,TPZMultiphysicsCompMesh *cmesh, std::ofstream &Erro, TPZVec<REAL> *Log,PreConfig &pConfig);

//// Solve desired problem
void Solve(ProblemConfig &config, PreConfig &preConfig);

//// Draw geometric and computational mesh
void DrawMesh(ProblemConfig &config, PreConfig &preConfig, TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *multiCmesh);

void NonConformAssemblage(TPZMultiphysicsCompMesh *multiCmesh,int InterfaceMatId, struct ProblemConfig config,struct PreConfig &pConfig, bool isHybridH1);

void CreateH1MixCompMesh(TPZCompMesh *cmesh, TPZMultiphysicsCompMesh *multiMix, TPZMultiphysicsCompMesh *multH1Mix,PreConfig &h1Config, ProblemConfig &ConfHyb);


void VerifyPragerSynge(int argc, char *argv[],PreConfig &pConfigH1, PreConfig &pConfigMix);

#endif //FEMCOMPARISON_ANALYTICS_H
