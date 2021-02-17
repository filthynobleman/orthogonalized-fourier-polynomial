/**
 * @file        write_ply_triv.c
 * 
 * @author      Filippo Maggioli (maggioli@di.uniroma1.it)
 *              'La Sapienza' Department of Computer Science
 * 
 * @brief       Writes a int32 matrix 3-by-m to a file and puts a uint8 value equals to 3
 *              at leading of each column.
 */

#include "io64.h"
#include "mex.h"


void write_matrix(FILE* stream, mxInt32* TRIS, mxUint32 m)
{
    register mxUint32 i;
    unsigned char three = 3;
    for (i = 0; i < m; i++)
    {
        fwrite(&three, sizeof(unsigned char), 1, stream);
        fwrite(TRIS + (3 * i), sizeof(mxInt32), 3, stream);
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mxInt32* TRIS = mxGetInt32s(prhs[0]);
    mxUint32 M = mxGetN(prhs[0]);
    mexPrintf("Number of triangles: %u\n", M);
    mexEvalString("drawnow");

    const char* filename = mxGetChars(prhs[1]);
    FILE *stream = fopen(filename, "ab");
    if (stream == NULL)
    {
        mexPrintf("Error while opening file %s.\n", filename);
        mexEvalString("drawnow");
        return;
    }
    write_matrix(stream, TRIS, M);
    fflush(stream);
    fclose(stream);
}