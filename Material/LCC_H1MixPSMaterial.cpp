//
// Created by victor on 29/03/2021.
//

#include "LCC_H1MixPSMaterial.h"
#include "pzaxestools.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("H1MIX-Material"));
#endif

LCC_H1MixMaterial::LCC_H1MixMaterial(int matid, int dim) : TPZMixedDarcyFlow(matid, dim)
{
}

LCC_H1MixMaterial::LCC_H1MixMaterial() : TPZMixedDarcyFlow()
{

}

LCC_H1MixMaterial::LCC_H1MixMaterial(const TPZMixedDarcyFlow &copy) : TPZMixedDarcyFlow(copy)
{
}

LCC_H1MixMaterial::LCC_H1MixMaterial(TPZMixedDarcyFlow &matdarcy,ProblemConfig &config, PreConfig &pConfig)
{
    this->SetId(matdarcy.Id());
    this->SetDimension(matdarcy.Dimension());
    TLaplaceExample1 *mat1 = new TLaplaceExample1, *mat2 = new TLaplaceExample1;

    TPZFNMatrix<9,STATE> K,invK;
    
#warning fixme
    K.Resize(3,3); invK.Resize(3,3);
    K.Identity();invK.Identity();
    K(0,0) = K(1,1) = pConfig.perm_Q1;
    invK(0,0) = invK(1,1) = 1./pConfig.perm_Q1;
    mat1->setPermeabilyTensor(K,invK);
    K(0,0) = K(1,1) = pConfig.perm_Q2;
    invK(0,0) = invK(1,1) =  1./pConfig.perm_Q2;
    mat2->setPermeabilyTensor(K,invK);
    
//    matdarcy.GetPermeability(K);
//    matdarcy.GetInvPermeability(invK);
//    this->SetPermeabilityTensor(K,invK);

    if (matdarcy.HasForcingFunction()) {
        //this->SetExactSol(matdarcy.GetExactSol());
        this->SetForcingFunction(matdarcy.ForcingFunction(),2);//2?????
    }
}

LCC_H1MixMaterial::~LCC_H1MixMaterial()
{

}

LCC_H1MixMaterial &LCC_H1MixMaterial::operator=(const LCC_H1MixMaterial &copy)
{
    TPZMixedDarcyFlow::operator=(copy);
    return *this;
}

void LCC_H1MixMaterial::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    
    //fem solution for flux and potential
    for(int i =0 ; i< 3; i++){
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
    }
}

/*
void LCC_H1MixMaterial::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec){
//    for(int i =0 ; i< 3; i++){
//        datavec[i].SetAllRequirements(false);
//        datavec[i].fNeedsSol = true;
//        datavec[i].fNeedsNormal = true;
//    }
//}
*/

void LCC_H1MixMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &terms){
    
    //     datavec[0] H1 potential
    //     datavec[1] Mixed potential
    //     datavec[2] Mixed Flux
    
    //    terms[0] = (h1P-u_exact[0])*(h1P-u_exact[0]);
    //    terms[1] = (mixP-u_exact[0])*(mixP-u_exact[0]);
    //    terms[2] = (h1P-mixP)*(h1P-mixP);
    
    //Pressure of H1 and Mixed
    STATE h1P, mixP;
    h1P = data[0].sol[0][0];
    mixP = data[1].sol[0][0];
    
    TPZVec<STATE> u_exact(1,0);
    TPZFMatrix<STATE> du_exact(3,1,0);
    if(this->fExactSol){
        this->fExactSol(data[0].x,u_exact,du_exact);
    }
    
    //Gradient of H1 and flux of mixed
    TPZFNMatrix<3,REAL> gradH1(3,1), mixF(3,1);
    gradH1 = data[0].dsol[0];
    for(int i=0 ; i<fDim; i++){
        mixF(i,0) = data[2].sol[0][i];
    }
    
    //MARK: Compute Prager-Synge terms
    int NTerms = 3;
    terms.Resize(NTerms);
    terms.Fill(0.0);
    
    // |grad(u-uh)|^2
    for(int i=0; i<fDim;i++){
        terms[0] += (du_exact[i]-gradH1[i])*(du_exact[i]-gradH1[i]);
    }
    
    // |grad(u)+sigma_h)|^2,  sigma_h belong to Hdiv
    for(int i=0; i<fDim;i++){
        terms[1] += (du_exact[i]+mixF[i])*(du_exact[i]+mixF[i]);
    }
    
    // |grad(uh)+sigma_h)|^2,  sigma_h belong to Hdiv
    for(int i=0; i<fDim;i++){
        terms[2] += (gradH1[i]+mixF[i])*(gradH1[i]+mixF[i]);
    }
    
//    double dummy = terms[0];
//    terms[0] = terms[1]+terms[0]-terms[2];
    
//    std::ofstream result0("H1MIXED.txt",std::ios_base::app);
//    result0.close();
    std::ofstream result("H1MIXED.txt",std::ofstream::app);

    result << "x: ( " << data[0].x[0] <<", " <<data[0].x[1] << " )\n";
    result << "gradU: ( " << du_exact[0] << ", " << du_exact[1] << " )\n";
    result << "gradH1: ( " << gradH1[0] << ", " << gradH1[1] << " )\n";
    result << "Sigma: ( " << mixF[0] << ", " << mixF[1] << " )\n";
    result << "terms[0]: " << terms[0] << " terms[1]: " << terms[1] << " terms[2]: " << terms[2] << "\n\n";
    
}

/*
void LCC_H1MixMaterial::Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors){

//datavec[0] HybH1 potential
//datavec[1] Mixed potential
//datavec[2] Mixed Flux
//
//error[0] - HybH1 potential exact error
//error[1] - Mixed potential exact error
//error[2] - HybH1 - Mixed error
//error[3] - HybH1 flux exact error (-k\nabla u)
//error[4] - Mixed flux exact error (\sigma)
//error[5] - HybH1 - Mixed flux



errors.Resize(NEvalErrors());
errors.Fill(0.0);

STATE hybP, mixP;

hybP = data[0].sol[0][0];
mixP = data[1].sol[0][0];

TPZVec<STATE> u_exact(1,0);
TPZFMatrix<STATE> du_exact(3,1,0);
if(this->fExactSol){
this->fExactSol->Execute(data[0].x,u_exact,du_exact);
}

errors[0] = (hybP-u_exact[0])*(hybP-u_exact[0]);
errors[1] = (mixP-u_exact[0])*(mixP-u_exact[0]);
errors[2] = (hybP-mixP)*(hybP-mixP);

TPZFNMatrix<3,REAL> gradHyb(3,1), mixF(3,1), hybF(3,1);

gradHyb = data[0].dsol[0];
for(int id=0 ; id<3; id++) {
    mixF(id,0) = data[2].sol[0][id];
}

gradHyb.Resize(3,1);

TPZFNMatrix<9,REAL> PermTensor;
TPZFNMatrix<9,REAL> InvPermTensor;
GetPermeabilities(data[0].x, PermTensor, InvPermTensor);
PermTensor.Resize(3,3);
InvPermTensor.Resize(3,3);

PermTensor.Multiply(gradHyb,hybF);
for(int ip = 0 ; ip < 3 ; ip++){
hybF(ip,0) = -hybF(ip,0);
}

TPZFNMatrix<3,REAL> flux;

{
TPZFNMatrix<9,REAL> minusGradP(3,1);
for (int i=0; i<3; i++) {
minusGradP(i,0) = (-1.)*du_exact[i];
}
PermTensor.Multiply(minusGradP,flux);
}

STATE hybExactF = 0., mixExactF = 0., DiffF = 0.;
for (int i=0; i<3; i++) {
for (int j=0; j<3; j++) {
hybExactF += (hybF[i]-flux(i,0))*InvPermTensor(i,j)*(hybF[j]-flux(j,0));
mixExactF += (mixF[i]-flux(i,0))*InvPermTensor(i,j)*(mixF[j]-flux(j,0));
DiffF += (hybF[i]-mixF[i])*InvPermTensor(i,j)*(hybF[i]-mixF[i]);
}
}

errors[3] = hybExactF;
errors[4] = mixExactF;
errors[5] = DiffF;
}
*/
/*
void LCC_H1MixMaterial::Errors(TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
    
//     datavec[0] HybH1 potential
//     datavec[1] Mixed potential
//     datavec[2] Mixed Flux
//
//      error[0] - HybH1 potential exact error
//      error[1] - Mixed potential exact error
//      error[2] - HybH1 - Mixed error
//      error[3] - HybH1 flux exact error (-k\nabla u)
//      error[4] - Mixed flux exact error (\sigma)
//      error[5] - HybH1 - Mixed flux

     

    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    
    STATE hybP, mixP;
    
    hybP = data[0].sol[0][0];
    mixP = data[1].sol[0][0];

    if(this->fExactSol){
        this->fExactSol->Execute(data[0].x,u_exact,du_exact);
    }

    errors[0] = (hybP-u_exact[0])*(hybP-u_exact[0]);
    errors[1] = (mixP-u_exact[0])*(mixP-u_exact[0]);
    errors[2] = (hybP-mixP)*(hybP-mixP);

    TPZFNMatrix<3,REAL> gradHyb(3,1), mixF(3,1), hybF(3,1);

    gradHyb = data[0].dsol[0];
    for(int id=0 ; id<3; id++) {
        mixF(id,0) = data[2].sol[0][id];
    }

    gradHyb.Resize(3,1);

    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;
    GetPermeabilities(data[0].x, PermTensor, InvPermTensor);
    PermTensor.Resize(3,3);
    InvPermTensor.Resize(3,3);
    
    PermTensor.Multiply(gradHyb,hybF);
    for(int ip = 0 ; ip < 3 ; ip++){
        hybF(ip,0) = -hybF(ip,0);
    }

    TPZFNMatrix<3,REAL> flux;

    {
        TPZFNMatrix<9,REAL> minusGradP(3,1);
        for (int i=0; i<3; i++) {
            minusGradP(i,0) = (-1.)*du_exact[i];
        }
        PermTensor.Multiply(minusGradP,flux);
    }

    STATE hybExactF = 0., mixExactF = 0., DiffF = 0.;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            hybExactF += (hybF[i]-flux(i,0))*InvPermTensor(i,j)*(hybF[j]-flux(j,0));
            mixExactF += (mixF[i]-flux(i,0))*InvPermTensor(i,j)*(mixF[j]-flux(j,0));
            DiffF += (hybF[i]-mixF[i])*InvPermTensor(i,j)*(hybF[i]-mixF[i]);
        }
    }

    errors[3] = hybExactF;
    errors[4] = mixExactF;
    errors[5] = DiffF;
}
*/
int LCC_H1MixMaterial::VariableIndex(const std::string &name)
{
    if(name == "hybP") return 40;
    if(name == "mixP") return 41;
    if(name == "exactP") return 42;
    if(name == "DiffP") return 43;

    if(name == "hybF") return 50;
    if(name == "mixF") return 51;
    if(name == "exactF") return 52;
    if(name == "DiffF") return 53;

    return -1;
}


int LCC_H1MixMaterial::NSolutionVariables(int var)
{
    switch (var) {
        case 40:
        case 41:
        case 42:
        case 43:
            return 3;
            break;
        case 50:
        case 51:
        case 52:
        case 53:
            return 1;
            break;
        default:
            DebugStop();
            break;
    }
    return 0;
}
/*
void LCC_H1MixMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{

    
//     datavec[0] HybH1 potential
//     datavec[1] Mixed potential
//     datavec[2] Mixed Flux
     

    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;

    GetPermeabilities(datavec[0].x, PermTensor, InvPermTensor);

    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> gradu(3,1,0.), fluxinv(3,1);

    if(fExactSol)
    {
        this->fExactSol->Execute(datavec[0].x, pressexact,gradu);
    }

    PermTensor.Multiply(gradu, fluxinv);

    TPZFMatrix<REAL> &dsolaxes = datavec[0].dsol[0];
    TPZFNMatrix<9, REAL> dsol(3, 1);
    TPZFNMatrix<9, REAL> KGradsol(3, 1);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, dsol, datavec[2].axes);

    PermTensor.Multiply(dsol, KGradsol);

    int dim=this->fDim;
    switch (var)
    {
        case 50://hybF
        {
            for (int i = 0; i < 3; i++) Solout[i]  = -KGradsol(i,0);
        }
            break;

        case 51:{// mixF
            for(int id=0 ; id<3; id++) Solout[id] = datavec[2].sol[0][id];
            }
            break;

        case 52://exactF
            for(int i=0; i<dim; i++) Solout[i] = -fluxinv(i);
            break;

        case 53://DiffF
            for(int i=0; i<dim; i++) Solout[i] = -KGradsol(i,0) - datavec[2].sol[0][i];
            break;

        case 40://hybP
            Solout[0] = datavec[0].sol[0][0];
            break;

        case 41://mixP
            Solout[0] = datavec[1].sol[0][0];
            break;

        case 42://exactP
            Solout[0] = pressexact[0];
            break;

        case 43://DiffP
            Solout[0] = datavec[0].sol[0][0] - datavec[1].sol[0][0];
            break;
        default:
            DebugStop();
    }
}
*/
