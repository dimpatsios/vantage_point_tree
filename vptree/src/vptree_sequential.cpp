/* Project no.1 on Parallel and Destributed Systems @AUTh - ECE Dept.
*  PATSIOS DIMITRIOS
*  DATE: 15/11/2019
* Version no.1 : sequential implementation
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../inc/quick_select.h"
#include "../inc/vptree.h"

#if 1
namespace
{
	// Calculates distances between two d-points.
	template < class T >
	T EuclideanDistance(T *X, T *Y, int d)
	{
		T sum(0.);
		for (int i(0); i < d; ++i) {
			sum += (X[i] - Y[i])*(X[i] - Y[i]);
		}
		return std::sqrt(sum);
	}

	double *DistanceMatrix(vptree *T, int *xID, double *X, int n, int d)
	{
		double *distances = new double[n - 1];
		if ( !distances ) return NULL;
		for ( int i(0); i < n - 1; ++i ) {
			double *tmp(&X[xID[i]*d]);
			distances[i] = EuclideanDistance<double>(T->vp, tmp, d);
		}
		return distances;
	}

	vptree *InitVptree(double *X, int n, int d, int *xID)
	{
		vptree *T = new vptree();
		int i(n - 1); // last one
//		int i(rand() % n);
		int index(xID[i]);

		if ( n == 1 ) {
			T->vp = X;
			T->idx = index;
			T->md = 0.;
			T->inner = NULL;
			T->outer = NULL;
			return T;
		}

		T->vp  = &X[index * d];
		T->idx = index;

		return T;
	}

	int *InitXid(int n)
	{
		int *Xid = new int[n];
		for( int i(0); i < n; ++i ) { Xid[i] = i; }
		return Xid;
	}

	vptree *buildVpRecursively(double *X, int n, int d, int *xID)
	{
		if ( n < 1) return NULL;

		//	printRandomArray(n, d, X);

		vptree *T = InitVptree(X, n, d, xID);

		if ( n == 1 ) return T;

		double *distances(DistanceMatrix(T, xID, X, n, d));
		if ( !distances ) return T;
		int k(static_cast<int>(std::floor(n * 0.5)));
		qcksel::QuickSelect<double> quickselect(distances, xID, 0, n - 2, k);
		T->md = quickselect.Result();
		delete [] distances;

		// call buildvpRecursively() recursively - Pre-order tree traversal
		T->inner = buildVpRecursively(X, k, d, &xID[0]);
		T->outer = buildVpRecursively(X, n - k - 1, d, &xID[k]);

		return T;
	}

}

vptree *buildvp(double *X, int n, int d)
{
	int *Xid(InitXid(n));
	vptree *T(buildVpRecursively(X, n, d, Xid));
	delete [] Xid;
	return T;
}

vptree * getInner(vptree * T)
{
	return T->inner;
}

vptree * getOuter(vptree * T)
{
	return T->outer;
}

double getMD(vptree * T)
{
	return T->md;
}

double *getVP(vptree *T)
{
	return T->vp;
}

int getIDX(vptree *T)
{
	return T->idx;
}

void destroy(vptree *T)
{
	if (T == NULL ) return;
	delete(T->inner);
	delete(T->outer);
	delete(T);
}

void printRandomArray(int const &reserve_size, int const reserve_dim, double *array)
{
	for ( int i(0); i < reserve_size; ++i) {
		for ( int j(0); j < reserve_dim; ++j ) {
			std::cout << array[i*reserve_dim + j] << "\t";
		}
		std::cout << std::endl;
	}
}
#endif
