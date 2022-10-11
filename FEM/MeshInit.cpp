//
// Created by victor on 02/09/2020.
//

#include "MeshInit.h"
#include "DataStructure.h"
#include <TPZMultiphysicsCompMesh.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "DarcyFlow/TPZHybridDarcyFlow.h"
#include "TPZNullMaterial.h"
#include "TPZBndCondT.h"
#include "TPZGenGrid2D.h"
#include "TPZCompMeshTools.h"
#include "TPZCompElDisc.h"
#include "LCC_H1MixPSMaterial.h"


void InsertMaterialMixed_MultiK(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig &config, PreConfig &pConfig){
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = config.gmesh->Dimension();
    int dirichlet = 0;
    int neumann = 1;

    int defaultOrder = config.k+config.n;
    cmesh_mixed -> SetDefaultOrder(defaultOrder);
    cmesh_mixed -> SetDimModel(dim);
    cmesh_mixed->SetAllCreateFunctionsMultiphysicElem();
#ifndef OPTMIZE_RUN_TIME
    TLaplaceExample1 *mat1 = new TLaplaceExample1,*mat2 = new TLaplaceExample1;

    SetFExact(mat1,mat2,pConfig);

    TPZFNMatrix<9,REAL> K, invK;
    K.Resize(3,3); invK.Resize(3,3);
    K.Identity();invK.Identity();
    K(0,0) = K(1,1) = pConfig.perm_Q1;
    invK(0,0) = invK(1,1) = 1./pConfig.perm_Q1;
    mat1->setPermeabilyTensor(K,invK);
    K(0,0) = K(1,1) = pConfig.perm_Q2;
    invK(0,0) = invK(1,1) =  1./pConfig.perm_Q2;
    mat2->setPermeabilyTensor(K,invK);
#endif
    
    TPZMixedDarcyFlow *material_Q1 = new TPZMixedDarcyFlow(matID_Q1,dim); //Using standard PermealityTensor = Identity.
    TPZMixedDarcyFlow *material_Q2 = new TPZMixedDarcyFlow(matID_Q2,dim);
#ifndef OPTMIZE_RUN_TIME
    material_Q1->SetForcingFunction(config.exact->ForceFunc(),5);
    material_Q1->SetExactSol(config.exact->ExactSolution(),5);
    material_Q2->SetForcingFunction(config.exact->ForceFunc(),5);
    material_Q2->SetExactSol(config.exact->ExactSolution(),5);
#endif
    
    material_Q1->SetConstantPermeability(pConfig.perm_Q1);
    material_Q2->SetConstantPermeability(pConfig.perm_Q2);

    cmesh_mixed->InsertMaterialObject(material_Q1);
    cmesh_mixed->InsertMaterialObject(material_Q2);

    //Boundary Conditions
    TPZFMatrix<REAL> val1(2,2,0.);
    TPZManVector<REAL> val2(2,0.);

    auto *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
#ifndef OPTMIZE_RUN_TIME
    BCond0_Q1->SetForcingFunctionBC(config.exact->ExactSolution(),5);
#endif
    auto *BCond0_Q2 = material_Q2->CreateBC(material_Q2, -6, dirichlet, val1, val2);
#ifndef OPTMIZE_RUN_TIME
    BCond0_Q2->SetForcingFunctionBC(config.exact->ExactSolution(),5);
#endif
    auto *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    auto *BCond1_Q2 = material_Q1->CreateBC(material_Q2, -9, neumann, val1, val2);

    cmesh_mixed->InsertMaterialObject(BCond0_Q1);
    cmesh_mixed->InsertMaterialObject(BCond0_Q2);
    cmesh_mixed->InsertMaterialObject(BCond1_Q1);
    cmesh_mixed->InsertMaterialObject(BCond1_Q2);
}

void InsertMaterialHybrid_MultiK(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config, PreConfig &pConfig){

    TPZGeoMesh* gmesh = cmesh_H1Hybrid->Reference();
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = gmesh->Dimension();
    int dirichlet = 0;
    int neumann = 1;

#ifndef OPTMIZE_RUN_TIME
    TLaplaceExample1 *mat1 = new TLaplaceExample1,*mat2 = new TLaplaceExample1;
    SetFExact(mat1,mat2,pConfig);

    TPZFNMatrix<9,REAL> K, invK;
    K.Resize(3,3); invK.Resize(3,3);
    K.Identity();invK.Identity();
    K(0,0) = K(1,1) = pConfig.perm_Q1;
    invK(0,0) = invK(1,1) = 1./pConfig.perm_Q1;
    mat1->setPermeabilyTensor(K,invK);
    K(0,0) = K(1,1) = pConfig.perm_Q2;
    invK(0,0) = invK(1,1) =  1./pConfig.perm_Q2;
    mat2->setPermeabilyTensor(K,invK);
#endif
    TPZHybridDarcyFlow *material_Q1 = new TPZHybridDarcyFlow(matID_Q1, dim);
    TPZHybridDarcyFlow *material_Q2 = new TPZHybridDarcyFlow(matID_Q2, dim);

    material_Q1->SetConstantPermeability(pConfig.perm_Q1);
    material_Q2->SetConstantPermeability(pConfig.perm_Q2);

    cmesh_H1Hybrid->InsertMaterialObject(material_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(material_Q2);

#ifndef OPTMIZE_RUN_TIME
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        material_Q1->SetForcingFunction(mat1->ForceFunc(),5);
        material_Q1->SetExactSol(mat1->ExactSolution(),5);

        material_Q2->SetForcingFunction(mat2->ForceFunc(),5);
        material_Q2->SetExactSol(mat2->ExactSolution(),5);
    }
#endif
    // Inserts boundary conditions
    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<REAL> val2(1, 1.);
    auto *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    auto *BCond0_Q2 = material_Q2->CreateBC(material_Q2, -6, dirichlet, val1, val2);
#ifndef OPTMIZE_RUN_TIME
    if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
        BCond0_Q1->SetForcingFunctionBC(mat1->ExactSolution(),5);
        BCond0_Q2->SetForcingFunctionBC(mat2->ExactSolution(),5);
    }
#endif
    val2.Fill(0.);
    auto *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    auto *BCond1_Q2 = material_Q1->CreateBC(material_Q1, -9, neumann, val1, val2);

    cmesh_H1Hybrid->InsertMaterialObject(BCond0_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(BCond0_Q2);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1_Q1);
    cmesh_H1Hybrid->InsertMaterialObject(BCond1_Q2);
}


void InsertMaterialH1Mix(TPZMultiphysicsCompMesh *multMesh, PreConfig &pConfig, ProblemConfig &config){

    int dim = pConfig.dim;
    int matID = 1;
    int dirichlet = -1;
    int neumann = -2;

    if(pConfig.type != 2) {// Dominio nao centrado na origen
        multMesh->SetDefaultOrder(pConfig.k + pConfig.n);
        multMesh->SetDimModel(dim);
        multMesh->SetAllCreateFunctionsMultiphysicElem();

        LCC_H1MixMaterial *material = new LCC_H1MixMaterial(matID, dim);//Using standard PermealityTensor = Identity.
        
#ifndef OPTMIZE_RUN_TIME
            material->SetForcingFunction(config.exact->ForceFunc(),5);
            material->SetExactSol(config.exact->ExactSolution(),5);
#else

        std::function<void (const TPZVec<REAL> &loc, TPZVec<STATE> &result)> sourceFunc = [](const TPZVec<REAL> &loc,TPZVec<STATE> &result)
        {
            for(auto &it:result) it = 1.;
        };

            material->SetForcingFunction(sourceFunc,0);
#endif
        
        multMesh->InsertMaterialObject(material);

        
        
//        //Boundary Conditions
//        TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
//
//        TPZMaterial *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
//        if(pConfig.debugger)
//        {
//            BCond0->SetForcingFunction(config.exact.operator*().Exact());
//        }
//
//        TPZMaterial *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);
//
//        multMesh->InsertMaterialObject(BCond0);
//        multMesh->InsertMaterialObject(BCond1);
    }

    else {
        DebugStop();
    }
}

#ifndef OPTMIZE_RUN_TIME
void SetFExact(TLaplaceExample1 *mat1, TLaplaceExample1 *mat2,PreConfig &pConfig){
    switch(pConfig.type) {
        case 0:
            mat1->fExact = TLaplaceExample1::ESinSin;
            mat2->fExact = TLaplaceExample1::ESinSin;
            break;
        case 1:
            mat1->fExact = TLaplaceExample1::EArcTan;
            mat2->fExact = TLaplaceExample1::EArcTan;
            break;
        case 2:
            mat1->fExact = TLaplaceExample1::ESteklovNonConst;
            mat2->fExact = TLaplaceExample1::ESteklovNonConst;
            break;
        case 3:
//            mat1->fExact = TLaplaceExample1::ESteepWave;
//            mat2->fExact = TLaplaceExample1::ESteepWave;
            DebugStop();
            break;
        default:
            DebugStop();
            break;
    }
    mat1->fSignConvention = 1;
    mat2->fSignConvention = 1;
}
#endif

void SetMultiPermeMaterials(TPZGeoMesh* gmesh){

    TPZAdmChunkVector<TPZGeoEl *> &elvec = gmesh->ElementVec();
    int numEl = elvec.NElements();
    TPZVec<int64_t> nodeInd;
    int numberNodes;
    TPZManVector<REAL,2> elCenterCoord = {0.,0.,0.};
    int meshDim = gmesh ->Dimension();

    if (meshDim == 2) {
        //Getting center coordinates of the element
        for (int ind = 0; ind < numEl; ind++) {
            TPZGeoEl *gel = elvec[ind];
            gel->GetNodeIndices(nodeInd);
            numberNodes = gel->NNodes();
            elCenterCoord[0] = elCenterCoord[1] = 0.;
            for (int nodeIter = 0; nodeIter < numberNodes; nodeIter++) {
                TPZGeoNode nodo = gmesh->NodeVec()[nodeInd[nodeIter]];
                elCenterCoord[0] += nodo.Coord(0);
                elCenterCoord[1] += nodo.Coord(1);
            }
            elCenterCoord[0] /= numberNodes;
            elCenterCoord[1] /= numberNodes;

            // Changing the permeability
            //std::cout <<"el: " << ind  << "x_ave: " << elCenterCoord[0] << "y_ave: " << elCenterCoord[1] << std::endl;
            if ((elCenterCoord[0] >= 0. && elCenterCoord[1] >= 0.) ||
                (elCenterCoord[0] < 0. && elCenterCoord[1] <= 0.)) { // first/third quadrants
                if (gel->Dimension() == gmesh->Dimension()) { /*if(gel->MaterialId() == 1)*/ gel->SetMaterialId(2); }
                else {/*if (gel->MaterialId() < 0 && gel->MaterialId() > -5)*/ gel->SetMaterialId(-5);
                } // Assumes dirichlet boundary condition
            } else {
                if (gel->Dimension() == gmesh->Dimension()) { /*if(gel->MaterialId() == 1)*/ gel->SetMaterialId(3); }
                else {/*if (gel->MaterialId() < 0 && gel->MaterialId() > -5)*/ gel->SetMaterialId(-6);
                }// Assumes dirichlet boundary condition
            }
        }
    }
    if (meshDim == 3) {
        //Getting center coordinates of the element
        for (int ind = 0; ind < numEl; ind++) {
            TPZGeoEl *gel = elvec[ind];
            gel->GetNodeIndices(nodeInd);
            numberNodes = gel->NNodes();
            elCenterCoord[0] = elCenterCoord[1] = 0., elCenterCoord[2] = 0.;
            for (int nodeIter = 0; nodeIter < numberNodes; nodeIter++) {
                TPZGeoNode nodo = gmesh->NodeVec()[nodeInd[nodeIter]];
                elCenterCoord[0] += nodo.Coord(0);
                elCenterCoord[1] += nodo.Coord(1);
                elCenterCoord[2] += nodo.Coord(2);
            }
            elCenterCoord[0] /= numberNodes;
            elCenterCoord[1] /= numberNodes;
            elCenterCoord[2] /= numberNodes;

            // Changing the permeability
            //std::cout <<"el: " << ind  << "x_ave: " << elCenterCoord[0] << "y_ave: " << elCenterCoord[1] << std::endl;
            if ((elCenterCoord[0] > 0. && elCenterCoord[1] > 0. &&  elCenterCoord[2] > 0.)  || //x,y,z > 0
                (elCenterCoord[0] > 0. && elCenterCoord[1] < 0. &&  elCenterCoord[2] < 0.)  || //y,z < 0, x > 0
                (elCenterCoord[0] < 0. && elCenterCoord[1] < 0. &&  elCenterCoord[2] > 0.)  || //x,y < 0, z > 0
                (elCenterCoord[0] < 0. && elCenterCoord[1] > 0. &&  elCenterCoord[2] < 0.))    //x,z < 0, y > 0
                {
                    if (gel->Dimension() == meshDim) gel->SetMaterialId(2);
                    else gel->SetMaterialId(-5);
                }
            else
                {
                    if (gel->Dimension() == meshDim) gel->SetMaterialId(3);
                    else gel->SetMaterialId(-6);
                }
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

void BuildFluxMesh(TPZCompMesh *cmesh_flux, ProblemConfig &config, PreConfig &pConfig){

    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    cmesh_flux->SetDefaultOrder(config.k);
    cmesh_flux->SetDimModel(config.gmesh->Dimension());

    cmesh_flux->SetAllCreateFunctionsHDiv();

    TPZNullMaterial<> *material = new TPZNullMaterial<>(matID);
    material->SetDimension(dim);
    cmesh_flux->InsertMaterialObject(material);

    //(EE)Create Boundary conditions
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<REAL> val2(2,0.);

    //(PROPRIETARY)
    auto *BCond0 = material->CreateBC(material,-1,dirichlet,val1,val2);
    cmesh_flux->InsertMaterialObject(BCond0);

    auto *BCond1 = material->CreateBC(material,-2,neumann,val1,val2);
    cmesh_flux->InsertMaterialObject(BCond1);

    //This part is for multi-K problemas
    if (pConfig.type == 2) BuildFluxMesh_MultiK(cmesh_flux,config);

    cmesh_flux->AutoBuild();
    cmesh_flux->InitializeBlock();
}

void CreateMixedAtomicMeshes(TPZVec<TPZCompMesh *> &meshvec, PreConfig &eData, ProblemConfig &config){
    //Flux mesh creation
    TPZCompMesh *cmesh_flux = new TPZCompMesh(config.gmesh);
    BuildFluxMesh(cmesh_flux,config,eData);
    {
        int64_t nelem = cmesh_flux->NElements();
        int dimgrid = cmesh_flux->Dimension();
        for (int64_t el = 0; el<nelem; el++) {
            TPZCompEl *cel = cmesh_flux->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            int dimgel = gel->Dimension();
            int nconnects = cel->NConnects();
            int nSides = gel->NSides(cmesh_flux->Dimension()-1);
            //if(dimgrid != 2) DebugStop();
            if(dimgel == dimgrid){
                if(nconnects != nSides+1){
                    DebugStop();
                }
                for(int ic = 0 ; ic < nSides; ic ++){
                    cel->Connect(ic).SetLagrangeMultiplier(0);
                }
            }
            if(dimgel == dimgrid-1){
                if(nconnects != nSides){
                    DebugStop();
                }
                cel->Connect(0).SetLagrangeMultiplier(0);
            }
            if (dimgel > dimgrid || dimgel < dimgrid - 1){
                DebugStop();
            }
        }
    }

    //Potential mesh creation
    TPZCompMesh *cmesh_p = new TPZCompMesh(config.gmesh);
    BuildPotentialMesh(cmesh_p,config,eData);
    {
        int64_t nconnects = cmesh_p->NConnects();
        for (int ic=0; ic<nconnects; ic++) {
            cmesh_p->ConnectVec()[ic].SetLagrangeMultiplier(1);
        }
    }


    TPZCompMesh *gspace = new TPZCompMesh(config.gmesh);
    {
        InsertNullSpaceMaterialIds(gspace,config);
        gspace->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        gspace->SetDefaultOrder(0);//sao espacos de pressao media
        gspace->AutoBuild();
        int64_t nconnects = gspace->NConnects();
        for (int ic = 0; ic<nconnects; ic++) {
            gspace->ConnectVec()[ic].SetLagrangeMultiplier(2);
        }
        int64_t nel = gspace->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(gspace->Element(el));
            if(disc) disc->SetFalseUseQsiEta();
        }
    }
    TPZCompMesh *average = new TPZCompMesh(config.gmesh);
    {
        InsertNullSpaceMaterialIds(average,config);
        average->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        average->SetDefaultOrder(0);
        average->AutoBuild();
        int64_t nconnects = average->NConnects();
        for (int ic = 0; ic<nconnects; ic++) {            average->ConnectVec()[ic].SetLagrangeMultiplier(3);
        }
        int64_t nel = average->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(average->Element(el));
            if(disc) disc->SetFalseUseQsiEta();
        }
}

    meshvec[0] = cmesh_flux;
    meshvec[1] = cmesh_p;
    meshvec[2] = gspace;
    meshvec[3] = average;
    TPZCompMeshTools::AdjustFluxPolynomialOrders(meshvec[0], config.n); //Increases internal flux order by "hdivmais"
    TPZCompMeshTools::SetPressureOrders(meshvec[0], meshvec[1]);//Set the pressure order the same as the internal flux
}

void InsertNullSpaceMaterialIds(TPZCompMesh *nullspace, ProblemConfig &config)
{
    for (auto matid:config.materialids) {
        TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(matid);
        nullmat->SetDimension(config.gmesh->Dimension());
        nullmat->SetNStateVariables(1);
        nullspace->InsertMaterialObject(nullmat);
    }
}


void BuildFluxMesh_MultiK(TPZCompMesh *cmesh_flux, ProblemConfig &config) {

    int dim = config.gmesh->Dimension();
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dirichlet = 0;
    int neumann = 1;

    TPZNullMaterial<> *material_Q1 = new TPZNullMaterial<>(matID_Q1);
    TPZNullMaterial<> *material_Q2 = new TPZNullMaterial<>(matID_Q2);
    material_Q1->SetDimension(dim);
    material_Q2->SetDimension(dim);
    cmesh_flux->InsertMaterialObject(material_Q1);
    cmesh_flux->InsertMaterialObject(material_Q2);

    //(EE)Create Boundary conditions
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<REAL> val2(2, 0.);

    //(PROPRIETARY)
    auto *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
    auto *BCond0_Q2 = material_Q1->CreateBC(material_Q2, -6, dirichlet, val1, val2);
    cmesh_flux->InsertMaterialObject(BCond0_Q1);
    cmesh_flux->InsertMaterialObject(BCond0_Q2);

    auto *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
    auto *BCond1_Q2 = material_Q2->CreateBC(material_Q2, -9, neumann, val1, val2);
    cmesh_flux->InsertMaterialObject(BCond1_Q1);
    cmesh_flux->InsertMaterialObject(BCond1_Q2);
}

void BuildPotentialMesh(TPZCompMesh *cmesh_p, ProblemConfig &config, PreConfig &pConfig){
    int matID = 1;
    int dim = config.gmesh->Dimension();

    int potential_order = config.k + config.n;
    cmesh_p->SetDefaultOrder(potential_order);
    cmesh_p->SetDimModel(dim);

    cmesh_p->SetAllCreateFunctionsContinuous(); //H1 functions
    cmesh_p->ApproxSpace().CreateDisconnectedElements(true);

    TPZNullMaterial<> *material = new TPZNullMaterial<>(matID);
    material->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material);

    if (pConfig.type == 2) BuildPotentialMesh_MultiK(cmesh_p, config);

    cmesh_p->AutoBuild();
    cmesh_p->ExpandSolution();

    TPZAdmChunkVector<TPZConnect> &nodeIt = cmesh_p->ConnectVec();
    for(auto &nodo : nodeIt){
        nodo.SetLagrangeMultiplier(1);
    }
}

void BuildPotentialMesh_MultiK(TPZCompMesh *cmesh_p, ProblemConfig &config){
    int matID_Q1 = 2;
    int matID_Q2 = 3;
    int dim = config.gmesh->Dimension();

    TPZNullMaterial<> *material_Q1 = new TPZNullMaterial<>(matID_Q1); material_Q1->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material_Q1);

    TPZNullMaterial<> *material_Q2 = new TPZNullMaterial<>(matID_Q2); material_Q2->SetDimension(dim);
    cmesh_p->InsertMaterialObject(material_Q2);
}

void InsertMaterialMixed(TPZMultiphysicsCompMesh *cmesh_mixed, ProblemConfig config, PreConfig &pConfig){

    int dim = config.gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    if(pConfig.type != 2) {
        cmesh_mixed->SetDefaultOrder(config.k + config.n);
        cmesh_mixed->SetDimModel(dim);
        cmesh_mixed->SetAllCreateFunctionsMultiphysicElem();

        TPZMixedDarcyFlow *material = new TPZMixedDarcyFlow(matID, dim); //Using standard PermealityTensor = Identity.
        
#ifndef OPTMIZE_RUN_TIME
            material->SetForcingFunction(config.exact->ForceFunc(),5);
            material->SetExactSol(config.exact->ExactSolution(),5);
#else

        std::function<void (const TPZVec<REAL> &loc, TPZVec<STATE> &result)> sourceFunc = [](const TPZVec<REAL> &loc,TPZVec<STATE> &result)
        {
            for(auto &it:result) it = 1.;
        };

            material->SetForcingFunction(sourceFunc,0);
#endif
        cmesh_mixed->InsertMaterialObject(material);

        //Boundary Conditions
        TPZFMatrix<STATE> val1(2, 2, 0.);
        TPZManVector<REAL> val2(2, 0.);
        
        auto *BCond0 = material->CreateBC(material, -1, dirichlet, val1, val2);
        
#ifndef OPTMIZE_RUN_TIME
            BCond0->SetForcingFunctionBC(config.exact->ExactSolution(),5);
#endif
        

        auto *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

#ifndef OPTMIZE_RUN_TIME
            BCond1->SetForcingFunctionBC(config.exact->ExactSolution(),5);
#endif
        
        cmesh_mixed->InsertMaterialObject(BCond0);
        cmesh_mixed->InsertMaterialObject(BCond1);
    }

    else {
        InsertMaterialMixed_MultiK(cmesh_mixed, config, pConfig);
    }
}

void InsertMaterialHybrid(TPZMultiphysicsCompMesh *cmesh_H1Hybrid, ProblemConfig &config, PreConfig &pConfig)
{
    TPZGeoMesh *gmesh = cmesh_H1Hybrid->Reference();
    int dim = gmesh->Dimension();
    int matID = 1;
    int dirichlet = 0;
    int neumann = 1;

    // Creates Poisson material
    if(pConfig.type != 2) {
        TPZHybridDarcyFlow *material = new TPZHybridDarcyFlow(matID, dim);
        material->SetConstantPermeability(1.);
        cmesh_H1Hybrid->InsertMaterialObject(material);
        cmesh_H1Hybrid->SetAllCreateFunctionsMultiphysicElem();

#ifndef OPTMIZE_RUN_TIME
            material->SetForcingFunction(config.exact->ForceFunc(),5);
            material->SetExactSol(config.exact->ExactSolution(),5);
#else
            material->SetConstantPermeability(1.);
#endif
        // Inserts boundary conditions
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<REAL> val2(1, 1.);
        auto *BCond0 =
                material->CreateBC(material, -1, dirichlet, val1, val2);
        
#ifndef OPTMIZE_RUN_TIME
            if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
                BCond0->SetForcingFunctionBC(config.exact->ExactSolution(),5);
                
            }
#endif
        
        val2.Fill(0.);
        auto *BCond1 = material->CreateBC(material, -2, neumann, val1, val2);

        cmesh_H1Hybrid->InsertMaterialObject(BCond0);
        cmesh_H1Hybrid->InsertMaterialObject(BCond1);
    }
    else {
        InsertMaterialHybrid_MultiK(cmesh_H1Hybrid, config, pConfig);
    }
}

TPZCompMesh* InsertCMeshH1(ProblemConfig &config, PreConfig &pConfig) {

    TPZCompMesh* cmesh = new TPZCompMesh(config.gmesh);
    TPZMaterialT<STATE>* mat = 0;
    int matID_Q1 = 1;
    int matID_Q2 = 1;
    int dirichlet = -1;
    int neumann = -2;
    int dim = config.gmesh->Dimension();

    if(pConfig.type != 2) {
        for (auto matid : config.materialids) {
            TPZDarcyFlow *mix = new TPZDarcyFlow(matid, cmesh->Dimension());
#ifndef OPTMIZE_RUN_TIME
            mix->SetExactSol(config.exact->ExactSolution(),5);
            mix->SetForcingFunction(config.exact->ForceFunc(),5);
#endif
            if (!mat) mat = mix;
            cmesh->InsertMaterialObject(mix);
        }

        for (auto matid : config.bcmaterialids) {
            TPZFNMatrix<1, REAL> val1(1, 1, 0.);
            TPZManVector<REAL> val2(1, 0.);
            int bctype = 0;
            val2.Fill(0.);
            auto *bc = mat->CreateBC(mat, matid, bctype, val1, val2);
#ifndef OPTMIZE_RUN_TIME
            bc->SetForcingFunctionBC(config.exact->ExactSolution(),5);
#endif

            cmesh->InsertMaterialObject(bc);
        }
    }
    else{
#ifndef OPTMIZE_RUN_TIME
        TLaplaceExample1 *mat1 = new TLaplaceExample1,*mat2 = new TLaplaceExample1;
        SetFExact(mat1,mat2,pConfig);

        TPZFNMatrix<9,REAL> K, invK;
        K.Resize(3,3); invK.Resize(3,3);
        K.Identity();invK.Identity();
        K(0,0) = K(1,1) = pConfig.perm_Q1;
        invK(0,0) = invK(1,1) = 1./pConfig.perm_Q1;
        mat1->setPermeabilyTensor(K,invK);
        K(0,0) = K(1,1) = pConfig.perm_Q2;
        invK(0,0) = invK(1,1) =  1./pConfig.perm_Q2;
        mat2->setPermeabilyTensor(K,invK);
#endif
        
        TPZDarcyFlow *material_Q1 = new TPZDarcyFlow(matID_Q1, dim);
        TPZDarcyFlow *material_Q2 = new TPZDarcyFlow(matID_Q2, dim);

        material_Q1->SetConstantPermeability(pConfig.perm_Q1);
        material_Q2->SetConstantPermeability(pConfig.perm_Q2);

#ifndef OPTMIZE_RUN_TIME
        if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
            material_Q1->SetForcingFunction(mat1->ForceFunc(),5);
            material_Q1->SetExactSol(mat1->ExactSolution(),5);

            material_Q2->SetForcingFunction(mat2->ForceFunc(),5);
            material_Q2->SetExactSol(mat2->ExactSolution(),5);
        }
#endif
        
        cmesh->InsertMaterialObject(material_Q1);
        cmesh->InsertMaterialObject(material_Q2);



        // Inserts boundary conditions
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<REAL> val2(1, 1.);
        auto *BCond0_Q1 = material_Q1->CreateBC(material_Q1, -5, dirichlet, val1, val2);
        auto *BCond0_Q2 = material_Q2->CreateBC(material_Q2, -6, dirichlet, val1, val2);
#ifndef OPTMIZE_RUN_TIME
        if (config.exact.operator*().fExact != TLaplaceExample1::ENone) {
            BCond0_Q1->SetForcingFunctionBC(mat1->ExactSolution(),5);
            BCond0_Q2->SetForcingFunctionBC(mat2->ExactSolution(),5);
        }
#endif
        val2.Fill(0.);
        auto *BCond1_Q1 = material_Q1->CreateBC(material_Q1, -8, neumann, val1, val2);
        auto *BCond1_Q2 = material_Q1->CreateBC(material_Q1, -9, neumann, val1, val2);

        cmesh->InsertMaterialObject(BCond0_Q1);
        cmesh->InsertMaterialObject(BCond0_Q2);
        cmesh->InsertMaterialObject(BCond1_Q1);
        cmesh->InsertMaterialObject(BCond1_Q2);
    }

    cmesh->SetDefaultOrder(config.k);//ordem
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

    cmesh->AutoBuild();


    return cmesh;
}
