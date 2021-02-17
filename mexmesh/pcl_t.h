/*
 * Emanuele Rodola
 * USI Lugano, Apr 2016
 */
 
#ifndef PCL_T_H
#define PCL_T_H

#include <vector>
#include <list>

#ifdef USE_FLANN
#include <flann/flann.hpp>
#endif

template <typename T>
struct vec3d
{
  T x,y,z;

  explicit vec3d() : x(0), y(0), z(0) {}
  vec3d(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

  template <typename S>
  vec3d(const vec3d<S>& s) : x(T(s.x)), y(T(s.y)), z(T(s.z)) {}
  
  vec3d<T> operator+(const vec3d<T>& p) const
  {
	 return vec3d<T>(x + p.x, y + p.y, z + p.z);
  }

  vec3d<T> operator-(const vec3d<T>& p) const
  {
	 return vec3d<T>(x - p.x, y - p.y, z - p.z);
  }

  vec3d<T> operator-() const
  {
	 return vec3d<T>(-x, -y, -z);
  }
  
  vec3d<T>& operator+=(const vec3d<T>& p)
  {
	 x += p.x;
	 y += p.y;
	 z += p.z;
	 return *this;
  }
  
  vec3d<T>& operator-=(const vec3d<T>& p)
  {
	 x -= p.x;
	 y -= p.y;
	 z -= p.z;
	 return *this;
  }
  
  template <typename Scalar>
  vec3d<T>& operator*=(const Scalar& s) {
	 x *= s;
	 y *= s;
	 z *= s;
	 return *this;
  }
  
  template <typename Scalar>
  vec3d<T> operator/(const Scalar& v) const
  {
	 return vec3d<T>(x/v, y/v, z/v);
  }
  
  template <typename Scalar>
  friend vec3d<T> operator*(const vec3d<T>& p, const Scalar& s) {
	 return vec3d<T>(s * p.x, s * p.y, s * p.z);
  }

  template <typename Scalar>
  friend vec3d<T> operator*(const Scalar& s, const vec3d<T>& p) {
	 return p*s;
  }
  
  template <typename Scalar>
  vec3d<T>& operator/=(const Scalar& s) {
	 const T i = T(1) / T(s);
	 x *= i;
	 y *= i;
	 z *= i;
	 return *this;
  }
};

template <typename T>
inline T dot_product(const vec3d<T>& v1, const vec3d<T>& v2) 
{
	return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
}

template <typename T>
inline vec3d<T> cross_product(const vec3d<T>& v1, const vec3d<T>& v2)
{
  return vec3d<T>(
		(v1.y*v2.z) - (v1.z*v2.y),
		(v1.z*v2.x) - (v1.x*v2.z),
		(v1.x*v2.y) - (v1.y*v2.x)
		);
}

template <typename T>
inline double magnitude(const vec3d<T>& v1)
{
	return std::sqrt((v1.x*v1.x) + (v1.y*v1.y) + (v1.z*v1.z));
}

template <typename T>
inline double normalize(vec3d<T>& p)
{
  const double n = magnitude(p);
  if (n==0.) {
	 p.x = 0;
	 p.y = 0;
	 p.z = 0;
  }
  else {
	 p.x /= n;
	 p.y /= n;
	 p.z /= n;
  }
  return n;
}

class oriented_point : public vec3d<double> {
   public:
      vec3d<double> n;

      explicit oriented_point() : vec3d<double>() {}
      oriented_point(const vec3d<double>& p) : vec3d<double>(p) {}
      oriented_point(double xx, double yy, double zz) : vec3d<double>(xx,yy,zz) {}
      oriented_point(const oriented_point& p) : vec3d<double>(p), n(p.n) {}
      oriented_point(const vec3d<double>& p, const vec3d<double>& nn) : vec3d<double>(p), n(nn) {}
   };

class pcl_t
{
private:

#ifdef USE_FLANN
   flann::Index< flann::L2<double> >* flann_index;
   double* flann_data;
   flann::Matrix<double>* flann_data_mat;
#endif

public:
	std::vector< vec3d<double> > points;
	std::vector< vec3d<double> > normals;
	
#ifdef USE_FLANN
	pcl_t() : flann_index(0), flann_data(0), flann_data_mat(0) {}
	
	~pcl_t() throw()
	{
	   if (flann_index)
	   {
			delete[] flann_data_mat->ptr();
			delete flann_data_mat;
			delete flann_index;
	   }
	}
#endif

#ifdef USE_FLANN
   void update_kd_tree()
   {
		if (flann_index)
	   {
			delete[] flann_data_mat->ptr();
			delete flann_data_mat;
			delete flann_index;
	   }
		
		flann_data = new double[3*points.size()];
	
		for (size_t i=0; i<points.size(); ++i)
		{
			const vec3d<double>& p = points[i];
			flann_data[3*i] = p.x;
			flann_data[3*i+1] = p.y;
			flann_data[3*i+2] = p.z;
		}
		
		flann_data_mat = new flann::Matrix<double>(flann_data, points.size(), 3);
		
		// NOTE: L2 actually gives the *squared* euclidean distance
		flann_index = new flann::Index< flann::L2<double> >(*flann_data_mat, flann::KDTreeSingleIndexParams());
		flann_index->buildIndex();
   }
#endif
	
	const vec3d<double>& get_vertex(int p) const { return points.at(p); }
	
	const vec3d<double>& get_normal(int p) const { return normals.at(p); }
	
#ifdef USE_FLANN
	/**
	 * NOTE: * The query point is NOT included in the neighbors list.
	 *       * The function returns SQUARED distances for efficiency reasons.
	 */
	void nearest_neighbors_with_dist(int p, double radius, std::vector<int>& neighs, std::vector<double>& dists)
	{
	   if (!flann_index)
			update_kd_tree();

		const vec3d<double>& pt = points.at(p);
		
		double q[3] = {pt.x, pt.y, pt.z};
		flann::Matrix<double> query((double*)&q, 1, 3);
		
		std::vector< std::vector<int> > out_idx;
		std::vector< std::vector<double> > out_dist;
		
		flann_index->radiusSearch(query, out_idx, out_dist, radius*radius, flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
		
		if (out_idx.size() != 1 || out_dist.size() != 1)
			std::cout << "[ERROR] mesh_t::nearest_neighbors_with_dist()" << std::endl;
			
		const std::vector<int>& ns = out_idx.front();
		const std::vector<double>& ds = out_dist.front();
		
		if (ns.front() != p)
		{
			std::cout << "[WARNING] mesh_t::nearest_neighbors_with_dist(): The first neighbor should be the query itself (point " << p << ")." << std::endl;
			neighs = ns;
			dists = ds;
		}
		else
		{
			neighs.resize(ns.size() - 1);
			std::copy(ns.begin()+1, ns.end(), neighs.begin());
			
			dists.resize(ns.size() - 1);
			std::copy(ds.begin()+1, ds.end(), dists.begin());
		}
	}
#endif
};

#endif
