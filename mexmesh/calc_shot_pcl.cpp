/*
 * Emanuele Rodola
 * USI Lugano, Apr 2016
 */
 
// mex -v calc_shot_pcl.cpp shot_descriptor_pcl.cpp -DUSE_FLANN -I"C:\Program Files\flann\include" -I../Eigen

#include "mex.h"
#include <vector>
#include "shot_descriptor_pcl.h"

//DEBUG
#include <fstream>
#include <iostream>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	if (nrhs != 6 || nlhs != 1)
      mexErrMsgTxt("Usage: descr = calc_shot_pcl(vertices, normals, indices, n_bins, radius, min_neighs).");

	// Parse input

	const double* const vertices = mxGetPr(prhs[0]);
	const int nv = int( mxGetN(prhs[0]) );
	
	if (int( mxGetM(prhs[0]) ) != 3)
		mexErrMsgTxt("Vertices should be given as a 3xN matrix.");
		
	const double* const normals = mxGetPr(prhs[1]);
	
	if (int( mxGetM(prhs[1]) ) != 3)
		mexErrMsgTxt("Normals should be given as a 3xN matrix.");
	
	if ( int( mxGetN(prhs[1]) ) != nv )
		mexErrMsgTxt("One normal vector per vertex is required.");
	
	// SHOT descriptors will be computed only for these points (1-based indices)
	const double* const idx = mxGetPr(prhs[2]);
	const int np = std::max( int( mxGetM(prhs[2]) ) , int( mxGetN(prhs[2]) ) );

	// Create point cloud structure
	
	pcl_t pcl;
	{
		pcl.points.resize(nv);
		pcl.normals.resize(nv);
		
		for (int i=0; i<nv; ++i)
		{
			vec3d<double>& p = pcl.points[i];
			p.x = vertices[i*3];
			p.y = vertices[i*3+1];
			p.z = vertices[i*3+2];
			vec3d<double>& n = pcl.normals[i];
			n.x = normals[i*3];
			n.y = normals[i*3+1];
			n.z = normals[i*3+2];
			//mexPrintf("%.2f %.2f %.2f - %.2f %.2f %.2f\n", p.x, p.y, p.z, n.x, n.y, n.z);
		}
		
	}
	
	
	// Compute SHOT descriptors at the desired point indices
	
	unibo::SHOTParams params;
	params.radius = *mxGetPr(prhs[4]);
	params.localRFradius = params.radius;
	params.minNeighbors = *mxGetPr(prhs[5]);
	params.bins = *mxGetPr(prhs[3]);
	
	unibo::SHOTDescriptor sd(params);
	const size_t sz = sd.getDescriptorLength();
	
	plhs[0] = mxCreateDoubleMatrix(sz, np, mxREAL);
	double* D = mxGetPr(plhs[0]);
	
	std::cout << "Computing SHOTs on " << np << " points... " << std::flush;

	for (size_t i=0; i<np; ++i)
	{
		unibo::shot s;
		sd.describe(pcl, static_cast<int>(idx[i]-1), s);
		for (size_t j=0; j<sz; ++j)
			D[i*sz+j] = s(j);
	}
	
	std::cout << "done." << std::endl;
}
