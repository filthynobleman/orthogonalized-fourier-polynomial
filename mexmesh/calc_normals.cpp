#include "mex.h"
#include <cmath>
#include <vector>
#include <queue>
#include <iostream>
#include "mesh.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	if (nrhs != 2 || nlhs != 1)
      mexErrMsgTxt("Usage: normals = calc_normals(vertices, triangles).");

	const double* const pts = mxGetPr(prhs[0]);
	const int np = int( mxGetN(prhs[0]) );
	const double* const tri = mxGetPr(prhs[1]);
	const int nt = int( mxGetN(prhs[1]) );
	
	if (np == 3)
		mexErrMsgTxt("It seems like you only have 3 vertices. Please try to transpose the input matrix.");
		
	if (nt == 3)
		mexErrMsgTxt("It seems like you only have 3 triangles. Please try to transpose the input matrix.");
	
	//std::cout << np << " vertices, " << nt << " triangles" << std::endl;
	
	// Load the mesh

	mesh_t mesh;
	
	std::vector< vec3d<double> > vertices(np);
	
	for (int i=0; i<np; ++i)
	{
		vec3d<double>& pt = vertices[i];
		pt.x = pts[i*3];
		pt.y = pts[i*3+1];
		pt.z = pts[i*3+2];
		//std::cout << pt << std::endl;
	}
	
	mesh.put_vertices(vertices);
	
	for (int i=0; i<nt; ++i)
	{
		int a, b, c;
		a = tri[i*3] - 1; // 1-based to 0-based
		b = tri[i*3+1] - 1;
		c = tri[i*3+2] - 1;
		mesh.add_triangle(a,b,c);
		//std::cout << a << " " << b << " " << c << std::endl;
	}
	
	// Return the normals per vertex

	plhs[0] = mxCreateDoubleMatrix(3, np, mxREAL);
	double* normals = mxGetPr(plhs[0]);
	
	for (size_t k = 0; k < np; ++k)
	{
		mesh_t::vertex_data& vert = mesh.vertices[k];
		vec3d<double> n(0,0,0);

		std::list<int>::const_iterator it = vert.tri_indices.begin(), it_end = vert.tri_indices.end();
		for (; it!=it_end; ++it)
		{
			 n += mesh.triangles[ *it ].n;
		}

		double mag = magnitude(n);
		n /= (mag != 0. ? mag : 1.); // some vertices might be outliers

		normals[k*3] = n.x;
		normals[k*3+1] = n.y;
		normals[k*3+2] = n.z;

	} // next vertex
}
