#include "mex.h"

/**
 * Computes the classical matrix multiplication between matrix A, which is
 * n-by-p, and matrix B, which is p-by-m, and stores the result in matrix C,
 * which is n-by-m.
 * 
 * @param A An n-by-p real matrix.
 * @param B An p-by-m real matrix.
 * @param C The resulting n-by-m real matrix.
 * @param n The number of rows of A and C.
 * @param m The number of columns of B and C.
 * @param p The number of columns of A and the number of rows of B.
 */
void matmul(const double *A, const double *B, double *C, 
            const size_t n, const size_t m, const size_t p)
{
    register int i, j, k;
    double a, b, c;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            c = 0;
            for (k = 0; k < p; k++)
            {
                a = A[k*n + i];
                b = B[j*p + k];
                c += a*b;
            }
            C[j*n + i] = c;
        }
    }
}

void repmatmul(const double *A, const double *B, double *C,
               const size_t n, const size_t m, const size_t p,
               const size_t nummats)
{
    register int i;
    for (i = 0; i < nummats; i++)
    {
        matmul(A + i * n*p,
               B + i * p*m,
               C + i * n*m,
               n, m, p);
    }
}

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    const mxArray *arg_A = prhs[0];
    const mxArray *arg_B = prhs[1];

    size_t n_dims = mxGetNumberOfDimensions(arg_A);
    if (n_dims != 3)
        mexErrMsgIdAndTxt("CGLib::repmatmul::InvalidDimensions",
            "The number of dimensions of input matrices must be 3.");
    if (mxGetNumberOfDimensions(arg_B) != n_dims)
        mexErrMsgIdAndTxt("CGLib::repmatmul::InvalidDimensions",
                "The number of dimensions of input matrices must be 3.");

    const size_t *dims_A = mxGetDimensions(arg_A);
    const size_t *dims_B = mxGetDimensions(arg_B);
    if (dims_A[1] != dims_B[0])
        mexErrMsgIdAndTxt("CGLib::repmatmul::InvalidDimensions",
                "Each matrix along the third dimensions of the input matrices, must be multipliable.");
    if (dims_A[2] != dims_B[2])
        mexErrMsgIdAndTxt("CGLib::repmatmul::InvalidDimensions",
                "The third dimensions of input matrices must agree.");


    const double *A = mxGetDoubles(arg_A);
    const double *B = mxGetDoubles(arg_B);
    size_t n = dims_A[0];
    size_t m = dims_B[1];
    size_t p = dims_A[1];
    size_t nummats = dims_A[2];

    size_t *resdims = (size_t *)calloc(3, sizeof(size_t));
    resdims[0] = n;
    resdims[1] = m;
    resdims[2] = nummats;
    plhs[0] = mxCreateNumericArray(3, resdims, mxDOUBLE_CLASS, mxREAL);
    double *C = mxGetDoubles(plhs[0]);

    repmatmul(A, B, C, n, m, p, nummats);
}