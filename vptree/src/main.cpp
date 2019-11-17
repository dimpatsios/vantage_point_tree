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

