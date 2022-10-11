

#ifndef FEMCOMPARISON_PRAGESYNGE_MAT_H
#define FEMCOMPARISON_PRAGESYNGE_MAT_H

#include <Poisson/TPZMatPoisson.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "DataStructure.h"


// 
class LCC_H1MixMaterial: public TPZMixedDarcyFlow  {
    
public:

    LCC_H1MixMaterial(int matid, int dim);

    LCC_H1MixMaterial();

    LCC_H1MixMaterial(TPZMixedDarcyFlow &matdarcy,ProblemConfig &config, PreConfig &pConfig);

    //LCC_H1MixMaterial(const LCC_H1MixMaterial &copy);

    LCC_H1MixMaterial(const TPZMixedDarcyFlow &copy);

    virtual ~LCC_H1MixMaterial();

    LCC_H1MixMaterial &operator=(const LCC_H1MixMaterial &copy);

    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;
    //virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    void ComputePragerSynge(TPZVec<TPZMaterialDataT<STATE>> &datavec,TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<STATE> &terms);
    
    virtual TPZMaterial * NewMaterial() {
        return new LCC_H1MixMaterial(*this);
    }

    int NEvalErrors() const override {return 3;}
    
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
    
    virtual int VariableIndex(const std::string &name);
    virtual int NSolutionVariables(int var);
    //virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
};


#endif //FEMCOMPARISON_PRAGESYNGE_MAT_H
