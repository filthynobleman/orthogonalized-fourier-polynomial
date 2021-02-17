/**
 * @file        create_ply_content.c
 * 
 * @author      Filippo Maggioli (maggioli@di.uniroma1.it)
 *              'La Sapienza' Department of Computer Science
 * 
 * @brief       Creates a buffer of mxUint8 containing the bytes to write to the PLY file.
 */

#include "mex.h"
#include "string.h"


void create_buffer(mxSingle* VERTS, mxInt32* TRIS, mxUint8* RGB,
                   mxUint32 N, mxUint32 M, mxArray** buffer_array)
{
    size_t bufsize = 3 * N * sizeof(mxSingle) +
                     3 * M * sizeof(mxInt32) + 
                     (RGB != NULL ? (3 * N * sizeof(mxUint8)) : 0) +
                     M * sizeof(mxUint8);
    
    *buffer_array = mxCreateNumericMatrix(bufsize, 1, mxUINT8_CLASS, mxREAL);
    mxUint8* buffer = mxGetUint8s(*buffer_array);
    mxUint8* bufptr = buffer;

    register int i;
    if (RGB == NULL)
    {
        memcpy(buffer, VERTS, 3 * N * sizeof(mxSingle));
        bufptr += (3 * N * sizeof(mxSingle)) / sizeof(mxUint8);
    }
    else
    {
        for (i = 0; i < N; i++)
        {
            memcpy(bufptr, VERTS + 3 * i, 3 * sizeof(mxSingle));
            bufptr += (3 * sizeof(mxSingle)) / sizeof(mxUint8);
            memcpy(bufptr, RGB + 3 * i, 3 * sizeof(mxUint8));
            bufptr += 3 * sizeof(mxUint8);
        }
    }

    mxUint8 three = 3;
    for (i = 0; i < M; i++)
    {
        memcpy(bufptr, &three, sizeof(mxUint8));
        bufptr += sizeof(mxUint8);
        memcpy(bufptr, TRIS + 3 * i, 3 * sizeof(mxInt32));
        bufptr += (3 * sizeof(mxInt32)) / sizeof(mxUint8);
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mxSingle* VERTS = mxGetSingles(prhs[0]);
    mxUint32 N = mxGetN(prhs[0]);
    mxInt32* TRIS = mxGetInt32s(prhs[1]);
    mxUint32 M = mxGetN(prhs[1]);
    mxUint8* RGB = NULL;
    if (nrhs > 2)
        RGB = mxGetUint8s(prhs[2]);
    
    create_buffer(VERTS, TRIS, RGB, N, M, &(plhs[0]));
}