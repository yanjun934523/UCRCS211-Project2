#include "../include/for_you_to_do.h"



/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    unsigned i = 0, j = 0, k = 0, t = 0;
    int maxind;
    double MAX;
    double* tempA = (double*)malloc(n * sizeof(double));

    for (i = 0; i < n; i++) {
        MAX = fabs(A[i * n + i]);
        maxind = i;
        for (t = i + 1; t < n; t++) {
            if (fabs(A[t * n + i]) > MAX) { 
                MAX = fabs(A[t * n + i]);
                maxind = t;
            }
        }

        if (MAX == 0) {
            return -1;
        }
        else {
            if (maxind != i) {
                int temp = ipiv[i];
                ipiv[i] = ipiv[maxind];
                ipiv[maxind] = temp;
                memcpy(tempA, A + i * n, n * sizeof(double));
                memcpy(A + i * n, A + maxind * n, n * sizeof(double));
                memcpy(A + maxind * n, tempA, n * sizeof(double));
            }
        }

        for (j = i + 1; j < n; j++) {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            for (k = i + 1; k < n; k++) {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }
        }
    }

    free(tempA);
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    unsigned i = 0, j = 0;
    double* y = (double*)malloc(n * sizeof(double)); 

    if (UPLO == 'L') { 
        y[0] = B[ipiv[0]]; 
        for (i = 1; i < n; i++) {
            double sum = 0.0;
            for (j = 0; j < i; j++) {  
                sum += A[i * n + j] * y[j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
    }
    if (UPLO == 'U') { 
        y[n - 1] = B[n - 1] / A[n * n - 1];
        for (i = n - 2; i >= 0; i--) { 
            double sum = 0.0;
            for (j = i + 1; j < n; j++) {
                sum += A[i * n + j] * y[j];
            }
            y[i] = (B[i] - sum) / A[i * n + i];
        }
    }
    memcpy(B, y, n * sizeof(double)); 
    free(y);
    return;
}


int get_block_size(){
    //return the block size you use in your matrix multiplication code.
    /*add your code here, 128 is an example and can be modified. */
    // With the returned block size the test code will adaptively adjust the input sizes to avoid corner cases.
    return 128;
  
}

//The matrix multiplication function used in blocked GEPP.
// You need to let the mydgemm adapt to non-square inputs.
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* A, B and C are n x n matrices.
    /* This function computes C[:i,:j]+=A[:i,:k]*B[:k,:j] (the first i rows and k columuns of A multiplies the first k rows and j columuns of B added to the the first i rows and j columuns of C)
    /* b is the "block size" used in the dgemm.
    /* In fact this function won't be directly called in the tester code, so you can modify the declaration (parameter list) of mydgemm() if needed. 
    /* you may copy the code from the optimal() function or any of the other functions in your lab1 code (optimized code recommended).*/
    /* add your code here */
   /* int i_1, j_1, k_1;
    for (i_1 = 0; i_1 < n; i_1 += 2) {
        for (j_1 = 0; j_1 < n; j_1 += 2) {
            register int t = i_1 * n + j_1;
            register int tt = t + n;
            register double c00 = C[t];
            register double c01 = C[t + 1];
            register double c10 = C[tt];
            register double c11 = C[tt + 1];
            for (k_1 = 0; k_1 < n; k_1 += 2) {
                register int ta = i_1 * n + k_1;
                register int tta = ta + n;
                register int tb = k_1 * n + j_1;
                register int ttb = tb + n;
                register double a00 = A[ta];
                register double a01 = A[ta + 1];
                register double a10 = A[tta];
                register double a11 = A[tta + 1];
                register double b00 = B[tb];
                register double b01 = B[tb + 1];
                register double b10 = B[ttb];
                register double b11 = B[ttb + 1];
                c00 += a00 * b00 + a01 * b10;
                c01 += a00 * b01 + a01 * b11;
                c10 += a10 * b00 + a11 * b10;
                c11 += a10 * b01 + a11 * b11;
            }
            C[t] = c00;
            C[t + 1] = c01;
            C[tt] = c10;
            C[tt + 1] = c11;
        }
    }*/
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    
    return 0;
}

