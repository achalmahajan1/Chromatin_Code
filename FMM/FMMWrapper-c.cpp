#include "FMMWrapper-c.h"
#include "../STKFMM.h"

extern "C" {

STKFMM *create_fmm_wrapper(int mult_order,int max_pts,int pbc,int kernel){
  return new STKFMM(mult_order,max_pts,static_cast<PAXIS>(pbc),kernel);
}

//void delete_fmm_wrapper(STKFMM* stkfmm){
//  delete stkfmm;
//}
    
void setBox(STKFMM* myFMM,double xlow,double xhigh,double ylow,double yhigh,double zlow,double zhigh){
  myFMM->setBox(xlow,xhigh,ylow,yhigh,zlow,zhigh);
}

void setupTree(STKFMM* myFMM, int kernel){
  myFMM->setupTree(static_cast<KERNEL>(kernel));
}


void setPoints(STKFMM* myFMM, const int nSL, const double *srcSLCoord, const int nDL, const double *srcDLCoord,
                   const int nTrg, const double *trgCoord){
  myFMM->setPoints(nSL,srcSLCoord,nDL,srcDLCoord,nTrg,trgCoord);
}
//void FMM_UpdateTree(FMM_Wrapper* fmm, const double* trg_coord, const double* src_coord, const int num_trg, const int num_src){
  // Copy arrays to vectors
//  std::vector<double> trg_coord_vec(trg_coord, trg_coord + 3*num_trg);
//  std::vector<double> src_coord_vec(src_coord, src_coord + 3*num_src);

  // Call method to update Tree
//  fmm->FMM_UpdateTree(src_coord_vec, trg_coord_vec);
//}
void clearFMM(STKFMM* myFMM, int kernel){
  myFMM->clearFMM(static_cast<KERNEL>(kernel));
}

void evaluateFMM(STKFMM* myFMM, const int nSL, const double *srcSLValuePtr, const int nDL, const double *srcDLValuePtr,
                         const int nTrg, double *trgValuePtr, const int kernelCh){
  //std::vector<double> trg_value_vec(num_trg * 3);

  // Copy array to vector
  // std::vector<double> src_value_vec(src_value, src_value + 3*num_src);

  // Call method to evaluate FMM
  myFMM->evaluateFMM(nSL,srcSLValuePtr,nDL,srcDLValuePtr,nTrg,trgValuePtr,static_cast<KERNEL>(kernelCh));
  
  // Copy vector to array
  // std::copy(trg_value_vec.begin(), trg_value_vec.end(), trg_value);
}
}
