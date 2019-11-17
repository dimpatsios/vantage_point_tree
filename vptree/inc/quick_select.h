#ifndef QUICK_SELECT_H
#define QUICK_SELECT_H

#include <iostream>
#include <vector>

#define PRINT() \
do { \
	std::cout << __FUNCTION__ << "\t" << __LINE__ << std::endl;\
} while(0)

#define NUM_POINTS 500000

namespace qcksel
{
	template <class T>
	class QuickSelect
	{
		public:

			explicit QuickSelect(T *A, int *xID, int low, int high, int k)
			{
				Result() = kthSmallest(xID, A, low, high, k);
			}

			explicit QuickSelect(std::vector<T> &arr, size_t const &k)
			{
				Result() = kth_smallest(arr.data(), 0, arr.size() - 1, arr.size() - k + 1);
			}

			~QuickSelect() { }

			// Accessors & mutators
			T &Result() { return m_result; }
			const T &Result() const { return m_result; }

		private:

			int partition(int *xID, T *arr, int l, int r)
			{
				double x = arr[r];
				int i = l;
				for (int j = l; j <= r - 1; j++) {
					if (arr[j] <= x) {
						std::swap(arr[i], arr[j]);
						std::swap(xID[i], xID[j]);
						i++;
					}
				}
				std::swap(arr[i], arr[r]);
				std::swap(xID[i], xID[r]);
				return i;
			}

			T kthSmallest(int *xID, T *arr, int l, int r, int k)
			{
				if (k > 0 && k <= r - l + 1) {

					int index = partition(xID,arr, l, r);

					if (index - l == k - 1) {
						return arr[index]; //it is equal return median value
					}

					if (index - l > k - 1) {
						return kthSmallest(xID,arr, l, index - 1, k); //if index greater than k set the boundaries to the left side of the array
					}

					return kthSmallest(xID,arr, index + 1, r,(k - index + l - 1)); //set the boundaries to the right side and set the k correct
				}

				return -1.;
			}

			T kth_smallest(T *A, size_t const &low, size_t const &high, size_t const &k)
			{
				if ( low == high ) return A[low];

				size_t pivotIndex = (low + high) / 2;
				pivotIndex = partition(A, low, high, pivotIndex);
				if ( k == pivotIndex ) return A[k];

				if ( k < pivotIndex ) {
					return kth_smallest(A, low, pivotIndex - 1, k);
				} else {
					return kth_smallest(A, pivotIndex + 1, high, k);
				}
			}


		private:

			T m_result;
	};

}
/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------*/
#endif // QUICK_SELECT_H
