/**
 * @file        dijkstra.c
 * @author      Filippo Maggioli (maggioli@di.uniroma1.it)
 * @brief       Compute the Dijkstra distance between each pair of vertices.
 */

#include "mex.h"
#include "string.h"

#define QUEUE_STEP_INC      50
#define QUEUE_STEP_RATIO    1.5
#define QUEUE_INC_MUL

#define PARENT(i)   (((i) - 1) / 2)
#define LEFT(i)     (2 * (i) + 1)
#define RIGHT(i)    (2 * (i) + 2)

typedef struct Node
{
    mxUint32        idx;
    mxDouble        val;
} node_t;

typedef struct Queue
{
    mxUint32*   heap;
    mxUint32    length;
    mxUint32    capacity;
} queue_t;


queue_t create_queue(mxUint32 n)
{
    if (n == 0)
        n = QUEUE_STEP_INC;
    queue_t q;
    q.heap = (mxUint32*)calloc(n, sizeof(mxUint32));
    q.length = 0;
    q.capacity = n;
    return q;
}

void heapify(queue_t* q, mxUint32 i, mxDouble* dist)
{
    mxUint32 l = LEFT(i);
    mxUint32 r = RIGHT(i);
    // if (l < q->length && dist[q->heap[l]] < dist[q->heap[i]])
    // {
    //     mxUint32 tmp = q->heap[l];
    //     q->heap[l] = q->heap[i];
    //     q->heap[i] = tmp;
    //     heapify(q, l, dist);
    // }
    // if (r < q->length && dist[q->heap[r]] < dist[q->heap[i]])
    // {
    //     mxUint32 tmp = q->heap[r];
    //     q->heap[r] = q->heap[i];
    //     q->heap[i] = tmp;
    //     heapify(q, r, dist);
    // }
    mxUint32 smallest = i;
    if (l < q->length && dist[q->heap[l]] < dist[q->heap[i]])
        smallest = l;
    if (r < q->length && dist[q->heap[r]] < dist[q->heap[smallest]])
        smallest = r;
    if (smallest != i)
    {
        mxUint32 tmp = q->heap[smallest];
        q->heap[smallest] = q->heap[i];
        q->heap[i] = tmp;
        heapify(q, smallest, dist);
    }
}

void push(queue_t* q, mxUint32 idx, mxDouble* dist)
{
    // Readjust heap capacity
    if (q->length == q->capacity)
    {
        q->heap = (node_t*)realloc(q->heap, (q->capacity + QUEUE_STEP_INC) * sizeof(node_t));
        #ifdef QUEUE_INC_MUL
            q->capacity = (mxUint32)(QUEUE_STEP_RATIO * q->capacity);
        #else
            q->capacity += QUEUE_STEP_INC;
        #endif
    }

    // Insert the new node
    mxUint32 i = q->length++;
    q->heap[i] = idx;

    while (i != 0 && dist[q->heap[PARENT(i)]] > dist[q->heap[i]])
    {
        mxUint32 tmp = q->heap[PARENT(i)];
        q->heap[PARENT(i)] = q->heap[i];
        q->heap[i] = tmp;
        i = PARENT(i);
    }
}

mxUint32 pop(queue_t* q, mxDouble* dist)
{
    if (q->length == 0)
        return -1;
    if (q->length == 1)
    {
        q->length = 0;
        return q->heap[0];
    }

    mxUint32 root = q->heap[0];
    q->heap[0] = q->heap[q->length - 1];
    q->length--;
    heapify(q, 0, dist);

    return root;
}

// mxUint32 find_nearest(queue_t *queue, mxDouble *dist, mxUint32 n)
// {
//     register mxUint32 i;
//     mxUint32 nearest = n;
//     for (i = 0; i < n; i++)
//     {
//         if (!queue[i])
//         {
//             if (nearest == n)
//                 nearest = i;
//             else
//                 nearest = dist[i] < dist[nearest] ? i : nearest;
//         }
//     }

//     queue[nearest] = true;
//     return nearest;
// }

void update_dists(mxUint32 *tris, mxDouble *verts, 
                  mxUint32 n, mxUint32 m,
                  mxDouble *dist, mxUint32 *prev, mxDouble *lengths, mxUint32 U,
                  queue_t* queue, mxLogical* explored)
{
    register mxUint32 i, j;
    for (i = 0; i < m; i++)
    {
        if (tris[i] == U || tris[m + i] == U || tris[2 * m + i] == U)
        {
            for (j = 0; j < 3; j++)
            {
                mxUint32 V = tris[j * m + i];
                mxDouble new_dist = dist[U] + lengths[V * n + U];
                if (new_dist < dist[V])
                {
                    prev[V] = U;
                    dist[V] = new_dist;
                }
            }
        }
    }
    heapify(queue, 0, dist);
}

void dijkstra(mxUint32 *tris, mxDouble *verts, 
              mxUint32 n, mxUint32 m,
              mxDouble *dist, mxUint32 *prev, mxDouble *lengths)
{
    mwSize dims = n;
    mxArray *explored_arr = mxCreateLogicalArray(1, &dims);
    mxLogical *explored = mxGetLogicals(explored_arr);
    mxUint32 U = 0;
    queue_t queue = create_queue(n);
    mxUint32 i;
    for (i = 0; i < n; i++)
        push(&queue, i, dist);
    // push(&queue, U, dist);
    register mxUint32 iter;
    // for (iter = 0; iter < n; iter++)
    //     queue[iter] = false;
    for (iter = 0; iter < n; iter++)
    {
        // U = find_nearest(&queue, dist, n);
        U = pop(&queue, dist);
        explored[U] = true;
        
        update_dists(tris, verts, n, m, dist, prev, lengths, U, &queue, explored);
    }
}

void dijkstra_all(mxUint32 *tris, mxDouble *verts, 
                  mxUint32 n, mxUint32 m,
                  mxDouble *dist, mxUint32 *prev, mxDouble *lengths)
{
    register mxUint32 i;
    for (i = 0; i < 1; i++)
        dijkstra(tris, verts, n, m, dist + (i * n), prev + (i * n), lengths);
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mexPrintf("Starting MEX program...\n");
    mexEvalString("drawnow");
    // Get the vertices and the triangles
    mxDouble *verts = mxGetDoubles(prhs[0]);
    mexPrintf("Retrieved vertices.\n");
    mexEvalString("drawnow");
    mxUint32 *tris = mxGetUint32s(prhs[1]);
    mexPrintf("Retrieved triangles.\n");
    mexEvalString("drawnow");
    // Get the number of vertices and triangles
    mxUint32 n = mxGetM(prhs[0]);
    mexPrintf("Retrieved number of vertices %d.\n", n);
    mexEvalString("drawnow");
    mxUint32 m = mxGetM(prhs[1]);
    mexPrintf("Retrieved number of triangles %d.\n", m);
    mexEvalString("drawnow");
    // Get the data matrices
    mxDouble *dist = mxGetDoubles(prhs[2]);
    mexPrintf("Retrieved distances.\n");
    mexEvalString("drawnow");
    mxUint32 *prev = mxGetUint32s(prhs[3]);
    mexPrintf("Retrieved parent vertices.\n");
    mexEvalString("drawnow");
    mxDouble *lengths = mxGetDoubles(prhs[4]);
    mexPrintf("Retrieved lengths.\n");
    mexEvalString("drawnow");

    dijkstra_all(tris, verts, n, m, dist, prev, lengths);

    // Return the distances
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    mxDouble *dijk = mxGetDoubles(plhs[0]);
    memcpy(dijk, dist, n * n * sizeof(mxDouble));
}