#ifndef __SPHERICALPARTICLE_HPP_
#define __SPHERICALPARTICLE_HPP_

template <class T>
struct sPoint3D {
	T x;
	T y;
	T z;
};

template <class T>
struct sParticleparams3D {
	T t;
	sPoint3D<T> r;
	sPoint3D<T> v;
	sPoint3D<T> a;
};

template <class T>
struct sParticleparams1D{
	T t;
	T r;
	T v;
	T a;
};

#endif