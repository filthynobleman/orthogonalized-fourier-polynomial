/*=========================================================================
 * compute_gradient.c
 *
 * This function computes the gradient of a function f on a mesh M, given
 * the value of a function at the vertices of each triangle and the 
 * difference vector of each edge.
 *
 * Filippo Maggioli
 */
#include "mex.h"
#include "stdlib.h"

double mydot(mxDouble *v1, mxDouble *v2, size_t n)
{
    double result = 0;
    register size_t i;
    for (i = 0; i < n; i++)
        result += v1[i] * v2[i];
    return result;
}

// Assuming g is 2-by-2
double mydet(mxDouble *g)
{
    return g[0]*g[3] - g[1]*g[2];
}

void myinv(mxDouble *g, mxDouble *ginv)
{
    double detg = mydet(g);
    ginv[0] = g[3] / detg;
    ginv[3] = g[0] / detg;
    ginv[1] = -g[1] / detg;
    ginv[2] = -g[2] / detg;
}

void tri_grad(mxDouble *e12, mxDouble *e13, 
              mxDouble f12, mxDouble f13,
              mxDouble *grad)
{
    // Create the matrix g
    mxArray *garr = mxCreateDoubleMatrix((mwSize)2, (mwSize)2, mxREAL);
    mxDouble *g = mxGetDoubles(garr);
    g[0] = mydot(e12, e12, 3);
    g[1] = mydot(e13, e12, 3);
    g[2] = mydot(e12, e13, 3);
    g[3] = mydot(e13, e13, 3);
    // Create the inverse of g
    mxArray *ginvarr = mxCreateDoubleMatrix((mwSize)2, (mwSize)2, mxREAL);
    mxDouble *ginv = mxGetDoubles(ginvarr);
    myinv(g, ginv);
    // Compute the product E * g^(-1) * F
    // First, E * g^(-1)
    mxArray *egarr = mxCreateDoubleMatrix((mwSize)3, (mwSize)2, mxREAL);
    mxDouble *eg = mxGetDoubles(egarr);
    eg[0] = e12[0]*ginv[0] + e13[0]*ginv[1];
    eg[1] = e12[1]*ginv[0] + e13[1]*ginv[1];
    eg[2] = e12[2]*ginv[0] + e13[2]*ginv[1];
    eg[3] = e12[0]*ginv[2] + e13[0]*ginv[3];
    eg[4] = e12[1]*ginv[2] + e13[1]*ginv[3];
    eg[5] = e12[2]*ginv[2] + e13[2]*ginv[3];
    // Then, (E * g^(-1)) * F
    grad[0] = eg[0]*f12 + eg[3]*f13;
    grad[1] = eg[1]*f12 + eg[4]*f13;
    grad[2] = eg[2]*f12 + eg[5]*f13;
}

void gradient(mxDouble *E12, mxDouble *E13,
              mxDouble *F12, mxDouble *F13,
              mxDouble *grad, size_t m)
{
    register size_t i;
    mxDouble *e12 = (mxDouble *)calloc(3, sizeof(mxDouble));
    mxDouble *e13 = (mxDouble *)calloc(3, sizeof(mxDouble));
    mxDouble *grad_loc = (mxDouble *)calloc(3, sizeof(mxDouble));
    for (i = 0; i < m; i++)
    {
        e12[0] = E12[i];
        e12[1] = E12[m + i];
        e12[2] = E12[m*2 + i];
        e13[0] = E13[i];
        e13[1] = E13[m + i];
        e13[2] = E13[m*2 + i];
        tri_grad(e12, e13, F12[i], F13[i], grad_loc);
        grad[i] = grad_loc[0];
        grad[m + i] = grad_loc[1];
        grad[m*2 + i] = grad_loc[2];
    }
}

/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    if (nlhs != 1)
        mexErrMsgIdAndTxt("CGCOURSE:compute_gradient:invalidNumOutput",
                "compute_gradient has only one output argument.");
    if (nrhs != 5)
        mexErrMsgIdAndTxt("CGCOURSE:compute_gradient:invalidNumInput",
                "compute_gradient has five input arguments.");
    
    // Get the number of triangles
    if (!mxIsUint32(prhs[4]))
        mexErrMsgIdAndTxt("CGCOURSE:compute_gradient:invalidInput",
                "compute_gradient's fifth argument must be int32.");
    size_t m = (size_t)*mxGetUint32s(prhs[4]);
    // Get the inputs
    // E12
    if (!mxIsDouble(prhs[0]))
        mexErrMsgIdAndTxt("CGCOURSE:compute_gradient:invalidInput",
                "compute_gradient's first argument must be double.");
    mxDouble *E12 = mxGetDoubles(prhs[0]);
    // E13
    if (!mxIsDouble(prhs[1]))
        mexErrMsgIdAndTxt("CGCOURSE:compute_gradient:invalidInput",
                "compute_gradient's second argument must be double.");
    mxDouble *E13 = mxGetDoubles(prhs[1]);
    // F12
    if (!mxIsDouble(prhs[2]))
        mexErrMsgIdAndTxt("CGCOURSE:compute_gradient:invalidInput",
                "compute_gradient's third argument must be double.");
    mxDouble *F12 = mxGetDoubles(prhs[2]);
    // F13
    if (!mxIsDouble(prhs[3]))
        mexErrMsgIdAndTxt("CGCOURSE:compute_gradient:invalidInput",
                "compute_gradient's fourth argument must be double.");
    mxDouble *F13 = mxGetDoubles(prhs[3]);
    // Initialize the output
    plhs[0] = mxCreateDoubleMatrix((mwSize)m, (mwSize)3, mxREAL);
    mxDouble *grad = mxGetDoubles(plhs[0]);
    // Compute the gradient
    gradient(E12, E13, F12, F13, mxGetDoubles(plhs[0]), m);
}