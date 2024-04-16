
#ifdef __cplusplus
extern "C" {
#endif

typedef struct STKFMM STKFMM;

STKFMM *create_fmm_wrapper(int mult_order,int max_pts,int pbc,int kernel);

//void delete_fmm_wrapper(STKFMM *fmm);

void setBox(STKFMM *myFMM, double xlow, double xhigh, double ylow, double yhigh, double zlow, double zhigh);

void setupTree(STKFMM *myFMM, int kernel);

void setPoints(STKFMM *myFMM, const int nSL, const double *srcSLCoord, const int nDL, const double *srcDLCoord,
                   const int nTrg, const double *trgCoord);
void clearFMM(STKFMM *myFMM, int kernel);

//void FMM_UpdateTree(STKFMM *fmm, const double *trg_coor, const double *src_coord, const int num_trg,
//                    const int num_src);

void evaluateFMM(STKFMM *myFMM,const int nSL, const double *srcSLValuePtr, const int nDL, const double *srcDLValuePtr,
                         const int nTrg, double *trgValuePtr, const int kernelCh);

#ifdef __cplusplus
}
#endif
