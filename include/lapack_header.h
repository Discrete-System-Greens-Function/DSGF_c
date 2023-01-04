// LAPACK
//Description of ZGELS: https://extras.csc.fi/math/nag/mark21/pdf/F08/f08anf.pdf
#define matsize 3*tot_sub_vol // the matrices have 3*tot_sub_vol rows and 3*tot_sub_vol columns
#define m matsize    //m: The number of rows of the matrix A (m≥ 0).
#define lda matsize  //lda: The leading dimension of a; at least max(1, m) for column major layout and at least max(1, n) for row major layout.
#define ldb matsize  //ldb: The leading dimension of b; must be at least max(1, m, n) for column major layout if trans='N' and at least max(1, n) if trans='T' and at least max(1, nrhs) for row major layout regardless of the value of trans.
#define n matsize    //n: The number of columns of the matrix A (n≥ 0).
#define nrhs matsize //nrhs: The number of right-hand sides; the number of columns in B (nrhs≥ 0).
//#define lwork 2*matsize //If LWORK =-1, a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns this value as the first entry of the WORK array, and no error message related to LWORK is issued. 

int info;
//double complex Alapack[lda*n], blapack[ldb*nrhs], work[lwork]; // based on example from https://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=506&p=1692&hilit=zgels#p1692

