//
// Created by victor on 02/09/2020.
//

#ifndef FEMCOMPARISON_MESHINIT_H
#define FEMCOMPARISON_MESHINIT_H

#include "DataStructure.h"
#include <TPZMultiphysicsCompMesh.h>
void CreateMixedAtomicMeshes(TPZVec<TPZCompMesh *> &meshvec, PreConfig &eData, ProblemConfig &config);

//// Insert volumetric and BC materials on a Primal Hybrid Computational Mesh
void InsertMaterialHybrid(TPZMultiphysicsCompMesh *cmesh, ProblemConfig &config,PreConfig &pConfig);

//// Insert volumetric and BC materials on a Mixed Computational Mesh
void InsertMaterialMixed(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig config,PreConfig &pConfig);

void InsertMaterialH1Mix(TPZMultiphysicsCompMesh *multMesh, PreConfig &pConfig, ProblemConfig &config);

//// Atomic Flux Mesh construction for mixed approximations
void BuildFluxMesh(TPZCompMesh *cmesh_flux, ProblemConfig &config, PreConfig &pConfig);

//// Atomic Flux Mesh construction for mixed approximations
//// Permeability differs on even and odd quadrants
void BuildFluxMesh_MultiK(TPZCompMesh *cmesh_flux, ProblemConfig &config);

//// Atomic Potential Mesh construction for mixed approximations
void BuildPotentialMesh(TPZCompMesh *cmesh_p, ProblemConfig &config, PreConfig &pConfig);

//// Atomic Potential Mesh construction for mixed approximations
//// Permeability differs on even and odd quadrants
void BuildPotentialMesh_MultiK(TPZCompMesh *cmesh_p, ProblemConfig &config);

//// Assign new volumetric material for even and odd quadrants
void SetMultiPermeMaterials(TPZGeoMesh* gmesh);

//// Create different materials for even an odd quadrants
void InsertMaterialHybrid_MultiK(TPZMultiphysicsCompMesh *cmesh_H1Hybrid,ProblemConfig &config, PreConfig &pConfig);

//// Create different materials for even an odd quadrants
void InsertMaterialMixed_MultiK(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config, PreConfig &pConfig);
void InsertNullSpaceMaterialIds(TPZCompMesh *nullspace, ProblemConfig &config);

#ifndef OPTMIZE_RUN_TIME
//// Set exact solution
void SetFExact(TLaplaceExample1 *mat1, TLaplaceExample1 *mat2,PreConfig &pConfig);
#endif

TPZCompMesh* InsertCMeshH1(ProblemConfig &config,PreConfig &pConfig);

#endif //FEMCOMPARISON_MESHINIT_H
