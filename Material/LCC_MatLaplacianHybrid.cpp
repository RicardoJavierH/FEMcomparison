//
//  LCC_MatLaplacianHybrid.cpp
//  ErrorEstimation
//
//  Created by Philippe Devloo on 14/07/19.
//

#include "LCC_MatLaplacianHybrid.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#ifdef FEMCOMPARISON_USING_MKL
#include "mkl.h"
#endif
#include "TPZTimer.h"
/*
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("MaterialHybrid"));
#endif
#ifdef PZ_LOG
static TPZLogger loggerCTM("contributeTimeVol");
#endif
#ifdef PZ_LOG
static TPZLogger loggerCTB("contributeTimeBoundary");
#endif
*/
#ifdef FEMCOMPARISON_TIMER
    //extern std::vector<unsigned long int> contributeTimeVec;
    extern long long contributeTimeVol;
    extern long long contributeTimeBC;

#endif

LCC_MatLaplacianHybrid::LCC_MatLaplacianHybrid(int matid, int dim)
: TPZRegisterClassId(&LCC_MatLaplacianHybrid::ClassId), TPZHybridDarcyFlow(matid, dim)
{
    
}

LCC_MatLaplacianHybrid::LCC_MatLaplacianHybrid() :
        TPZRegisterClassId(&LCC_MatLaplacianHybrid::ClassId), TPZHybridDarcyFlow()
{
    
}

LCC_MatLaplacianHybrid::LCC_MatLaplacianHybrid(const LCC_MatLaplacianHybrid &copy) :
        TPZRegisterClassId(&LCC_MatLaplacianHybrid::ClassId), TPZHybridDarcyFlow(copy)
{
    
}

LCC_MatLaplacianHybrid::~LCC_MatLaplacianHybrid()
{
    
}

LCC_MatLaplacianHybrid &LCC_MatLaplacianHybrid::operator=(const LCC_MatLaplacianHybrid &copy)
{
    TPZHybridDarcyFlow::operator=(copy);
    return *this;
}

TPZMaterial *LCC_MatLaplacianHybrid::NewMaterial() const
{
    return new LCC_MatLaplacianHybrid(*this);
}

int LCC_MatLaplacianHybrid::ClassId() const
{
    return Hash("LCC_MatLaplacianHybrid") ^ TPZHybridDarcyFlow::ClassId() << 1;
}


void LCC_MatLaplacianHybrid::Write(TPZStream &buf, int withclassid) const
{
    TPZHybridDarcyFlow::Write(buf,withclassid);
}

void LCC_MatLaplacianHybrid::Read(TPZStream &buf, void *context)
{
    TPZHybridDarcyFlow::Read(buf,context);
}

int LCC_MatLaplacianHybrid::VariableIndex(const std::string &name) const
{

    if(name == "Pressure") return 44;
    if(name == "PressureExact") return 45;

    if(name == "Flux") return 10;
    if(name == "ExactFlux") return 13;

    if(name == "ExactFluxShiftedOrigin") return 23;

    return -1;
}

int LCC_MatLaplacianHybrid::NSolutionVariables(int var) const {
    if(var == 44 || var==45) return 1;
    if(var == 10 || var==13 || var == 23) return fDim;
    
    else{
        
        return TPZHybridDarcyFlow::NSolutionVariables(var);
        return 0;
        
    }
}

void LCC_MatLaplacianHybrid::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    /**
     datavec[1] L2 mesh (phi's)
     datavec[0] Hdiv mesh,
     datavec[2] Interface Mesh
     datavec[3] Interface Mesh
     
     Implement the matrix
     |Sk Ck^T |  = |f1|
     |Ck  0   |    |f2|
     Sk = int_K K graduk.gradv dx = int_K K gradphi_i.gradphi_j dx
     CK = int_partialK lambda_k*uk dx = int_K phi_i dx
     f1 = int_K f*v dx = int_K f*phi_j dx
     ck = int_partialK phi_i*mu_j dx
     f2 = int_partialK g*mu_j dx
     **/
#ifdef FEMCOMPARISON_TIMER
      extern bool contributeTest;
//extern double contributeTimeVol;
#ifdef TIMER_CONTRIBUTE
       // auto begin = std::chrono::high_resolution_clock::now();
#endif
//    extern int64_t contributeMaterialCounter;
//    TPZTimer timer;
//    if(contributeTest){
//        timer.start();
//    }
      
#endif
    
    TPZFMatrix<REAL>  &phi = datavec[1].phi;
    TPZFMatrix<REAL> &dphi = datavec[1].dphix;
    TPZVec<REAL>  &x = datavec[1].x;

    int phr = phi.Rows();
    
    STATE fXfLoc = 0;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        fForcingFunction(x, res);
        fXfLoc = res[0];
    }
    
    STATE KPerm = GetPermeability(datavec[0].x);
    
#ifdef FEMCOMPARISON_USING_MKL
    {
        double *A, *B, *C;
        int m, n, k;
        double alpha, beta;
        m = phr, k = fDim, n = phr;
        alpha = weight*KPerm;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDA = dphi.Rows();
        LDB = dphi.Rows();
        LDC = ek.Rows();
        if(LDC != phr+2)
            DebugStop();
        A = &dphi(0,0);
        B = &dphi(0,0);
        C = &ek(0,0);
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                       phr, phr, LDA,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    {
        //saxpy implementation
        int N = phr;
        double alpha = weight*fXfLoc;
        double *X = &phi(0,0);
        int incX = 1;
        double *Y = &ef(0,0);
        int incY = 1;
        cblas_daxpy(N,alpha,X,incX, Y,incY);
    }
    {
        //
        int N = phr;
        double alpha = weight;
        double *X = &phi(0,0);
        int incX = 1;
        double *Y = &ek(0,phr);
        int incY = 1;
        cblas_daxpy(N,alpha,X,incX, Y,incY);
        Y = &ek(phr,0);
        incY = ek.Rows();
        cblas_daxpy(N,alpha,X,incX, Y,incY);
    }
#else
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        
        //matrix Sk
        for( int jn = 0; jn < phr; jn++ ) {
            for(kd=0; kd<fDim; kd++) {
                ek(in,jn) += (STATE)weight*(KPerm*(STATE)(dphi(kd,in)*dphi(kd,jn)));
            }
        }
    }
    for (int in =0; in < phr; in++) {
        ek(phr,in) += weight*phi(in,0);//lambda*phi
        ek(in,phr) += weight*phi(in,0);
    }
#endif
    //equacoes de restricao de pressao media
    ek(phr,phr+1) -= weight;
    ek(phr+1,phr) -= weight;
    
#ifdef FEMCOMPARISON_DEBUG
    if ( !ek.VerifySymmetry(1.e-10) ) {
        std::cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << std::endl;
    }
#endif
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream valuenn;
        ek.Print("ek = ",valuenn,EMathematicaInput);
        ef.Print("ef = ",valuenn,EMathematicaInput);
        LOGPZ_DEBUG(logger,valuenn.str());
    }
#endif
    
#ifdef FEMCOMPARISON_TIMER
#ifdef TIMER_CONTRIBUTE
        //auto end = std::chrono::high_resolution_clock::now();
        //auto contribVolElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        //contributeTimeVol += contribVolElapsed.count();
//        std::cout<< contributeTimeVol<<std::endl;
//        timer.stop();
//        contributeTimeVol += timer.seconds();
//        contributeMaterialCounter++;
#endif
#endif
}

void LCC_MatLaplacianHybrid::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<REAL>  &phi = datavec[1].phi;
    TPZFMatrix<REAL> &dphi = datavec[1].dphix;
    TPZVec<REAL>  &x = datavec[1].x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();
    
    STATE fXfLoc = 0;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        fForcingFunction(x, res);
        fXfLoc = res[0];
    }
    
    STATE KPerm = GetPermeability(datavec[0].x);
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        for(kd=0; kd<fDim; kd++) {
            ef(in,0) -= (STATE)weight*(KPerm*(STATE)(dphi(kd,in)*datavec[1].dsol[0](kd,0)));
        }
    }
    ef(phr,0) += weight*(-datavec[1].sol[0][0]+ datavec[3].sol[0][0]);
    ef(phr+1,0) += weight*datavec[2].sol[0][0];
}

void LCC_MatLaplacianHybrid::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc)
{
#ifdef FEMCOMPARISON_TIMER
    extern bool contributeTest;
#ifdef TIMER_CONTRIBUTE
        auto begin = std::chrono::high_resolution_clock::now();
#endif
//    extern double contributeTimeBoundary;
//    extern int64_t contributeBoundaryCounter;
//    TPZTimer timer;
//    if(contributeTest){
//        timer.start();
//    }
#endif
    
    TPZFMatrix<REAL>  &phi_u = datavec[1].phi;
    TPZFMatrix<REAL>  &phi_flux = datavec[0].phi;
    //    TPZFMatrix<REAL> &axes = data.axes;
    int phr_primal = phi_u.Rows();
    int phr_hybrid = phi_flux.Rows();
    bool primal = true;
    TPZManVector<REAL,3> x(3);
    if(phr_hybrid)
    {
        primal = false;
        x = datavec[0].x;
    }
    else
    {
        x = datavec[1].x;
    }
    short in,jn;
    STATE v2[1];
    v2[0] = bc.Val2()[0];
    
    if(bc.HasForcingFunctionBC()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
        TPZManVector<STATE> res(1);
        TPZFNMatrix<3,STATE> dres(3,1);
        bc.ForcingFunctionBC()(x, res, dres);
        v2[0] = res[0];
    }
    
    if(primal)
    {
        switch (bc.Type()) {
            case 0 :            // Dirichlet condition
                for(in = 0 ; in < phr_primal; in++) {
                    ef(in,0) += (STATE)(fBigNumber* phi_u(in,0) * weight) * v2[0];
                    for (jn = 0 ; jn < phr_primal; jn++) {
                        ek(in,jn) += fBigNumber * phi_u(in,0) * phi_u(jn,0) * weight;
                    }
                }
                break;
            case 1 :            // Neumann condition
                for(in = 0 ; in < phr_primal; in++) {
                    ef(in,0) += v2[0] * (STATE)(phi_u(in,0) * weight);
                }
                break;
            case 2 :        // mixed condition
                for(in = 0 ; in < phr_primal; in++) {
                    ef(in, 0) += v2[0] * (STATE)(phi_u(in, 0) * weight);
                    for (jn = 0 ; jn < phi_u.Rows(); jn++) {
                        ek(in,jn) += bc.Val1()(0,0) * (STATE)(phi_u(in,0) * phi_u(jn,0) * weight);     // peso de contorno => integral de contorno
                    }
                }
                break;
            default:
                DebugStop();
        }
    } else
    {
        switch (bc.Type()) {
            case 0 :            // Dirichlet condition
                for(in = 0 ; in < phr_hybrid; in++) {
                    ef(in,0) += v2[0] * (STATE)(phi_flux(in,0) * weight);
                }
                break;
            case 1 :            // Neumann condition
                for(in = 0 ; in < phr_hybrid; in++) {
                    ef(in,0) += (STATE)(fBigNumber* phi_flux(in,0) * weight) * v2[0];
                    for (jn = 0 ; jn < phr_hybrid; jn++) {
                        ek(in,jn) += fBigNumber * phi_flux(in,0) * phi_flux(jn,0) * weight;
                    }
                }
                break;
            case 2 :        // mixed condition
                DebugStop();
                for(in = 0 ; in < phr_hybrid; in++) {
                    ef(in, 0) += v2[0] * (STATE)(phi_flux(in, 0) * weight);
                    for (jn = 0 ; jn < phi_flux.Rows(); jn++) {
                        ek(in,jn) += 1./bc.Val1()(0,0) * (STATE)(phi_flux(in,0) * phi_flux(jn,0) * weight);     // peso de contorno => integral de contorno
                    }
                }
                break;
        }
    }
#ifdef FEMCOMPARISON_TIMER
#ifdef TIMER_CONTRIBUTE
    auto end = std::chrono::high_resolution_clock::now();
    auto contribBCelapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    contributeTimeBC += contribBCelapsed.count();
#endif
//    if(contributeTest){
//        timer.stop();
//        contributeTimeBoundary+=timer.seconds();
//        contributeBoundaryCounter++;
//    }
#endif
}


void LCC_MatLaplacianHybrid::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout)
{
    /**
     datavec[1] L2 mesh
     datavec[0] Hdiv mesh,
     datavec[2] Interface Mesh
     datavec[3] Interface Mesh
     **/
    if(var == 0)
    {
        TPZHybridDarcyFlow::Solution(datavec,var,Solout);
        return;
    }
    
    STATE KPerm = GetPermeability(datavec[0].x);
    STATE invKPerm = 1./KPerm;
    TPZFNMatrix<9,REAL> PermTensor;
    TPZFNMatrix<9,REAL> InvPermTensor;
    for(int i=0; i< fDim; i++){
        for(int j=0; j<fDim; j++){
            if(i==j){
                PermTensor(i,j)=KPerm;
                InvPermTensor(i,j)= invKPerm;
            } else{
                PermTensor(i,j) = 0;
                InvPermTensor(i,j)= 0;
            }
        }
    }
    
    TPZManVector<STATE,2> pressexact(1,0.);
    TPZFNMatrix<9,STATE> grad(fDim,1,0.), fluxinv(fDim,1),gradu(fDim,1,0);//no TPZAnalytic solution grad é 3x1
    
    if(fExactSol)
    {
        this->fExactSol(datavec[1].x, pressexact,grad);
        
        for(int i = 1; i<fDim ; i++){
            
            gradu(i,0)=grad(i,0);
        }
        
        
    }
    
    PermTensor.Multiply(gradu, fluxinv);
    

    switch (var)
    {
  

        case 44://PressureFem
            Solout[0] = datavec[1].sol[0][0];
            break;
        case 45://Pressure Exact
            Solout[0] = pressexact[0];
            break;
        case 10:
        case 13:
        case 23:
            TPZHybridDarcyFlow::Solution(datavec,var,Solout);
            break;
        default:
            DebugStop();
    }
}



void LCC_MatLaplacianHybrid::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    if(!fExactSol) return;

    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE> u_exact(1);
    TPZFNMatrix<9,STATE> du_exact;


    if(this->fExactSol){

        this->fExactSol(data[1].x,u_exact,du_exact);
    }

    REAL pressure = data[1].sol[0][0];

    // errors[0] norm L2 || u ||_l2

    errors[0] = (pressure-u_exact[0])*(pressure-u_exact[0]);//exact error pressure

    // errors[1] Semi norm H1 || grad u ||_l2

    TPZManVector<STATE,3> sol(1),dsol(3,0.);

    TPZFMatrix<REAL> &dsolaxes = data[1].dsol[0];
    TPZFNMatrix<9,REAL> flux(3,0);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, flux, data[1].axes);

    for(int id=0; id<fDim; id++) {
        REAL diff = fabs(flux(id,0) - du_exact(id,0));
        errors[1]  += diff*diff;
    }

    // error[2] H1 norm

    errors[2] = errors[0] +errors[1];

    // error[3] Energy norm || u ||_e = a(u,u)= int_K K gradu.gradu dx

    STATE KPerm = GetPermeability(data[0].x);
    
    
    TPZFNMatrix<9,REAL> gradpressure(fDim,1),Kgradu(fDim,1);
    for (int i=0; i<fDim; i++) {
        gradpressure(i,0) = du_exact(i,0);
        Kgradu(i,0) = gradpressure(0)*KPerm;
    }
        
    REAL energy = 0.;
    for (int i=0; i<fDim; i++) {
        for (int j=0; j<fDim; j++) {
            double cperm =0.;
            if(i==j)
                cperm = KPerm;
            energy += cperm*fabs(flux(j,0) - du_exact(j,0))*fabs(flux(i,0) - du_exact(i,0));
        }
    }

    errors[3] = energy;
}

void LCC_MatLaplacianHybrid::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++ )
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
        datavec[i].fNeedsHSize = false;
    }
}
