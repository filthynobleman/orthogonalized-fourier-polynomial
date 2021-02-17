/**
 * @file        write_ply_verts_rgb.c
 * 
 * @author      Filippo Maggioli (maggioli@di.uniroma1.it)
 *              'La Sapienza' Department of Computer Science
 * 
 * @brief       Writes two matrices VERTS (in double precision floating point) and RGB
 *              (in uint8 precision) to a binary file. The matrices must be both 3-by-n
 *              and will be written to the file alterning the rows.
 */

#include "io64.h"
#include "mex.h"

#define MYREAL mxSingle


void write_matrices(FILE* stream, MYREAL* VERTS, mxUint8* RGB, mxUint32 n)
{
    register mxUint32 i;
    for (i = 0; i < n; i++)
    {
        fwrite(VERTS + (3 * i), sizeof(MYREAL), 3, stream);
        fwrite(RGB + (3 * i), sizeof(mxUint8), 3, stream);
    }
}



void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    MYREAL* VERTS = mxGetSingles(prhs[0]);
    mxUint8* RGB = mxGetUint8s(prhs[1]);
    mxUint32 N = mxGetN(prhs[0]);
    mexPrintf("Number of vertices: %u\n", N);
    mexEvalString("drawnow");

    const char* filename = mxGetChars(prhs[2]);
    FILE *stream = fopen(filename, "ab");
    // FILE* stream;
    // fopen_s(&stream, filename, "ab");
    if (stream == NULL)
    {
        mexPrintf("Error while opening file %s.\n", filename);
        mexEvalString("drawnow");
        return;
    }
    //write_matrices(stream, VERTS, RGB, N);
    register mxUint32 i;
    for (i = 0; i < N; i++)
    {
        fwrite(VERTS + (3 * i), sizeof(MYREAL), 3, stream);
        fwrite(RGB + (3 * i), sizeof(mxUint8), 3, stream);
    }
    fflush(stream);
    fclose(stream);
}