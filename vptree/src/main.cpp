#if 0
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <time.h>
#include <cstdlib>
#include <sys/time.h>
#include "../inc/quick_select.h"
#include "../inc/vptree.h"

//#include <chrono>

namespace
{
	void randomArray(int const reserve_size, int const reserve_dim, double *array)
	{
		PRINT();
		srand(static_cast<unsigned int>(time(NULL)));

		for ( int i(0); i < reserve_size; ++i) {
			for ( int j(0); j < reserve_dim; ++j ) {
				double number = rand() / (double(RAND_MAX)*1e-04);
				array[i*reserve_dim + j] = number;
			}
		}
		PRINT();
	}

}

int main(void)
{
	int D = 3, N = 10;

	printf("Enter Dimension: \n");
	scanf("%du", &D);

	printf("Size: \n");
	scanf("%du", &N);

	if ( D < 1 || N < 1 ) {
		printf("Exit !!! Invalid input data.");
		return 0;
	}

	double *array = new double[N*D+1];
	if ( !array ) return 0;
	randomArray(N, D, array);
//	printRandomArray(N, D, array);

#if 0 // Problems with cilk and -std=gnu++11, thus changed to C03 the whole compilation
	auto start = std::chrono::system_clock::now();

	vptree *vptr = buildvp(array, N, D);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds(end - start);
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::cout << "finished computation at " << std::ctime(&end_time)
	          << "TOTAL time vptree:: "     << elapsed_seconds.count() << "s\n";
#else
	struct timeval swtime, ewtime;
	gettimeofday(&swtime, NULL);
	vptree *vptr = buildvp(array, N, D);
	gettimeofday(&ewtime, NULL);
	double time = (double)((ewtime.tv_usec - swtime.tv_usec) / 1.0e6 + ewtime.tv_sec - swtime.tv_sec);
	printf("TOTAL time vptree:: %f sec\n", time);
#endif

	destroy(vptr);
	delete [] array;

	return 0;
}
#else

/*!
  \file   tester.c
  \brief  Validate Vantage Point implementation.

  \author Dimitris Floros
  \date   2019-10-19
*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "../inc/vptree.h"
#include "../inc/quick_select.h"

//#define VERBOSE 1

static char * STR_CORRECT_WRONG[] = {"WRONG", "CORRECT"};


// =========================================================
// === STACK IMPLEMENTATION TO TEST RECUSRION INVARIANCE ===
// =========================================================

// Declare linked list node

typedef struct node {
  vptree *T;
  int isInner;
  double md;
  struct node* link;
} node;


// Utility function to add an element data in the stack
// insert at the beginning
void push(node **top, vptree *T, double md, int isInner)
{
  // create new node temp and allocate memory
  node *temp = (node *) malloc( sizeof(node) );

#ifdef VERBOSE
  fprintf( stdout, "PUSH ENTERED\n");
#endif

  // check if stack (heap) is full. Then inserting an element would
  // lead to stack overflow
  if (!temp) {
	fprintf( stderr, "\npush: Heap Overflow\n" );
	exit(1);
  }

  // initialize data into temp data field
  temp->T       = T;
  temp->isInner = isInner;
  temp->md      = md;

  // put top pointer reference into temp link
  temp->link = *top;

  // make temp as top of Stack
  *top = temp;

}

// Remove top element from stack
void pop(node **top)
{

#ifdef VERBOSE
	fprintf( stdout, "POP ENTERED\n");
#endif

	// check for stack underflow
	if (*top == NULL) {
		fprintf( stderr, "\npop: Stack Underflow\n" );
	  exit(1);
	}

	// top assign into temp
	struct node* temp = *top;

	// return second as top
	*top = (*top)->link;

	// destroy connection between first and second
	temp->link = NULL;

	// release memory of top node
	free(temp);

}




// =================
// === UTILITIES ===
// =================

double dist(double *X, double *Y, int d){
  double dist2 = 0;
  for (int i = 0; i < d; i++){
	dist2 += (X[i] - Y[i])*(X[i] - Y[i]);
  }
  return sqrt(dist2);
}

// Function to print all the
// elements of the stack
int verifyLeafPlace(node **top, double *X, int d)
{


  // check for stack underflow
  if (*top == NULL) {
	fprintf( stderr, "\nverifyLeafPlace: Stack Underflow\n" );
	exit(1);
  }

  struct node* temp = *top;

  // iterate the ancestors in stack
  while (temp != NULL) {

#ifdef VERBOSE
	fprintf( stdout, "%f | %f | %d == %d\n", dist(X, getVP(temp->T), d), temp->md, temp->isInner, dist(X, getVP(temp->T), d) <= temp->md );
#endif

	// check whether point should be inside or outside
	int isInner =  dist(X, getVP(temp->T), d) <= temp->md;

	// if the direction is invalid, break and return false
	if ( isInner != temp->isInner ) return 0;

	// assign temp link to temp
	temp = temp->link;
  }

  return 1;

}




// ==================
// === VALIDATION ===
// ==================

int *foundInTree;

int verifyTree(vptree *T, double *vp, node **stack, double md, int isInner,
               int n, int d){

  int isValid = 1;

  // if empty, return
  if (T == NULL) return isValid;

  int isValidAncestor = 1;
  int isLeaf = (getInner(T) == NULL && getOuter(T) == NULL);

#ifdef VERBOSE
  fprintf( stdout, "%x %x\n", stack, *stack);
#endif

  // if leaf check ancestor path
  if (isLeaf) isValidAncestor = verifyLeafPlace(stack, getVP(T), d);


  // if inner, radius must be smaller than parent's diameter
  if (isInner && getMD(T) > 2*md) return 0;

  // update list of indices
  int idx = getIDX(T);
  if (idx < n)
	foundInTree[ idx ] = 1;

  // validate distance to parent
  if (isInner)
	isValid = dist(vp, getVP(T), d) <= md;
  else
	isValid = dist(vp, getVP(T), d) > md;

  // recurse if not leaf
  if (!isLeaf) {

	// add to stack as inner and recurse, then pop
	push( stack, T, getMD(T), 1 );
	int isValidInn = verifyTree( getInner( T ), getVP(T), stack, getMD(T), 1, n, d );
	pop( stack );

	// add to stack as outer and recurse, then pop
	push( stack, T, getMD(T), 0 );
	int isValidOut = verifyTree( getOuter( T ), getVP(T), stack, getMD(T), 0, n, d );
	pop( stack );

	// all conditions must be true
	isValid = isValid && isValidInn && isValidOut;

  } else {

	// make sure ancestory path is correct
	isValid = isValid && isValidAncestor;

  }

  // return
  return isValid;

}

int main()
{

  int n=100;//data
  int d=4;//dimensions

  double  * dataArr = (double * ) malloc( n*d * sizeof(double) );
  double  * zeros   = (double * ) calloc( d   , sizeof(double) );

  foundInTree = (int *) calloc( n, sizeof(int) );

  for (int i=0;i<n*d;i++)
	dataArr[i]=rand()%100;

  vptree *root=buildvp(dataArr,n,d);

  node *stack = NULL;

  int isValid = verifyTree(root, zeros, &stack, 1.0/0.0, 1, n, d );
  int foundAll = 1;

  assert( stack == NULL );

  for (int i = 0; i<n; i++)
	if (!foundInTree[i]){
		foundAll = 0;
	  break;
	}

  printf("Tester validation: %s PROPERTIES | %s INDICES\n",
           STR_CORRECT_WRONG[isValid], STR_CORRECT_WRONG[foundAll]);

  return 0;

}

#endif
