/**
 * @file utility.cuh
 * @brief Provides utility functions and classes for CUDA operations.
 *
 * This file includes common CUDA helper functions, error checking macros,
 * memory management wrappers, a timer class, and mathematical utilities
 * for both host and device code.
 */
#pragma once
#include <bit> // bit_width
#include <cstddef> // size_t
#include <cstdint> // uintmax_t
#include <cstdio> // fprintf
#include <cuda_runtime.h> // CUDA runtime
#include <limits> // std::numeric_limits

/**
 * @def BLOCKS
 * @brief Default number of blocks to launch in CUDA kernels.
 */
#define BLOCKS 64
/**
 * @def THREADS
 * @brief Default number of threads per block in CUDA kernels.
 */
#define THREADS 1024

/**
 * @def gpuErrchk(ans)
 * @brief Macro to check CUDA API calls for errors.
 *
 * Calls gpuAssert to check the result of a CUDA API call.
 * @param ans The result of the CUDA API call (cudaError_t).
 */
#define gpuErrchk(ans) { gpuAssert(ans, __FILE__, __LINE__); }

/**
 * @brief Asserts if a CUDA error code indicates failure.
 *
 * Checks the provided CUDA error code. If it's not `cudaSuccess`, it prints
 * an error message including the error string, file name, and line number.
 * If `abort` is true, it terminates the program.
 *
 * @param code The CUDA error code to check.
 * @param file The source file where the check occurs (__FILE__).
 * @param line The line number where the check occurs (__LINE__).
 * @param abort If true (default), exits the program upon error.
 */
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);

		if (abort)
		{
			exit(code);
		}
	}
}

/**
 * @brief Compile-time flag to enable/disable CUDA error checking via gpuErrchk/gpuAssert.
 * @note Set to `true` to enable runtime checks (potentially slower), `false` to disable them.
 */
static constexpr bool CHECK_CUDA_ERRORS = false;

/**
 * @brief Computes the floor of log base 2 for a host integer using std::bit_width.
 * @param n The input value. Must be greater than 0.
 * @return size_t The largest integer `k` such that 2^k <= n.
 * @note Host code only.
 */
__host__ inline size_t host_fast_log2(size_t n)
{
	// std::bit_width(n) returns number of bits needed to represent n.
	// If n is power of 2, say 8 (1000), bit_width is 4. log2(8) = 3. So bit_width - 1.
	// If n is not power of 2, say 7 (0111), bit_width is 3. log2(7) is approx 2.8. floor(log2(7)) = 2. So bit_width - 1.
	return std::bit_width(n) - 1;
}

/**
 * @brief Computes the ceiling of log base 2 for a host integer using std::bit_width.
 * @param n The input value. Must be greater than 0.
 * @return size_t The smallest integer `k` such that 2^k >= n.
 * @note Host code only.
 */
__host__ inline size_t host_fast_log2_ceil(size_t n)
{
	// std::bit_width(n-1) returns number of bits needed for n-1.
	// If n is power of 2, say 8 (1000), n-1 is 7 (0111). bit_width(7) = 3. ceil(log2(8)) = 3.
	// If n is not power of 2, say 7 (0111), n-1 is 6 (0110). bit_width(6) = 3. ceil(log2(7)) = 3.
	// Handles n=1 case correctly: n-1=0, bit_width(0)=0. ceil(log2(1))=0.
	return std::bit_width(n - 1);
}

/**
 * @brief Computes the floor of log base 2 for a device integer using count leading zeros.
 * @param n The input value. Must be greater than 0.
 * @return size_t The largest integer `k` such that 2^k <= n.
 * @note Device code only. Uses `__clzll` (count leading zeros for long long).
 */
__device__ __forceinline__ size_t device_fast_log2(size_t n)
{
	// For size_t (unsigned long long on many 64-bit systems), digits is 64.
	// __clzll(n) counts leading zeros.
	// If n=8 (0...01000), clzll=60. digits-1-clzll = 64-1-60 = 3. log2(8)=3.
	// If n=7 (0...00111), clzll=61. digits-1-clzll = 64-1-61 = 2. floor(log2(7))=2.
	return std::numeric_limits<size_t>::digits - 1 - __clzll(n);
}

/**
 * @brief Computes the ceiling of log base 2 for a device integer using count leading zeros.
 * @param n The input value. Must be greater than 0.
 * @return size_t The smallest integer `k` such that 2^k >= n.
 * @note Device code only. Uses `__clzll` (count leading zeros for long long).
 */
__device__ __forceinline__ size_t device_fast_log2_ceil(size_t n)
{
	// If n=8 (0...01000), n-1=7 (0...0111). clzll(7)=61. digits-clzll(n-1) = 64-61 = 3. ceil(log2(8))=3.
	// If n=7 (0...00111), n-1=6 (0...0110). clzll(6)=61. digits-clzll(n-1) = 64-61 = 3. ceil(log2(7))=3.
	// Handles n=1: n-1=0. clzll(0) is undefined but often returns digits (64). digits-digits = 0. ceil(log2(1))=0.
	return std::numeric_limits<size_t>::digits - __clzll(n - 1);
}

/**
 * @brief A simple CUDA timer class using CUDA events.
 *
 * Measures the execution time of CUDA operations between calls to
 * `start_timer()` and `stop_timer()`.
 */
class CudaTimer
{
public:
	/**
	 * @brief Constructs the CudaTimer and creates CUDA events.
	 */
	CudaTimer()
	{
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
	}

	/**
	 * @brief Destroys the CUDA events.
	 */
	~CudaTimer()
	{
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
	}

	/**
	 * @brief Records the start event in the default CUDA stream.
	 */
	void start_timer()
	{
		cudaEventRecord(start);
	}

	/**
	 * @brief Records the stop event, synchronizes, and calculates elapsed time.
	 * @return float The elapsed time in milliseconds between the last calls to
	 * `start_timer()` and `stop_timer()`.
	 */
	float stop_timer()
	{
		cudaEventRecord(stop);
		cudaEventSynchronize(stop); // Wait for the stop event to complete
		cudaEventElapsedTime(&time, start, stop); // Calculate time difference
		return time;
	}

private:
	cudaEvent_t start, stop; ///< CUDA events to mark start and stop points.
	float time;              ///< Stores the calculated elapsed time in milliseconds.
};

/**
 * @brief Allocates memory on the CUDA device.
 *
 * Wrapper around `cudaMalloc` with optional error checking.
 *
 * @tparam T The data type of the array elements.
 * @param size The number of elements of type T to allocate.
 * @return T* Pointer to the allocated device memory.
 * @note Uses `gpuErrchk` if `CHECK_CUDA_ERRORS` is true.
 */
template <typename T>
inline T* alloc_to_device(size_t size)
{
	T* d_arr;

	if constexpr (CHECK_CUDA_ERRORS)
	{
		gpuErrchk(cudaMalloc(&d_arr, size * sizeof(T)));
	}
	else
	{
		cudaMalloc(&d_arr, size * sizeof(T));
	}

	return d_arr;
}

/**
 * @brief Copies data from host memory to existing device memory.
 *
 * Wrapper around `cudaMemcpy` (HostToDevice) with optional error checking.
 *
 * @tparam T The data type of the array elements.
 * @param d_arr Pointer to the destination device memory (must be pre-allocated).
 * @param arr Pointer to the source host memory.
 * @param size The number of elements of type T to copy.
 * @return T* The pointer to the destination device memory (`d_arr`).
 * @note Uses `gpuErrchk` if `CHECK_CUDA_ERRORS` is true.
 */
template <typename T>
inline T* copy_to_device(T* d_arr, const T* arr, size_t size)
{
	if constexpr (CHECK_CUDA_ERRORS)
	{
		gpuErrchk(cudaMemcpy(d_arr, arr, size * sizeof(T), cudaMemcpyHostToDevice));
	}
	else
	{
		cudaMemcpy(d_arr, arr, size * sizeof(T), cudaMemcpyHostToDevice);
	}

	return d_arr;
}

/**
 * @brief Allocates device memory and copies data from host memory to it.
 *
 * Combines allocation (`alloc_to_device`) and copy (`copy_to_device`).
 *
 * @tparam T The data type of the array elements.
 * @param arr Pointer to the source host memory.
 * @param size The number of elements of type T to allocate and copy.
 * @return T* Pointer to the newly allocated and populated device memory.
 * @note Uses `gpuErrchk` indirectly if `CHECK_CUDA_ERRORS` is true.
 */
template <typename T>
inline T* copy_to_device(const T* arr, size_t size)
{
	T* d_arr = alloc_to_device<T>(size);
	copy_to_device(d_arr, arr, size);
	return d_arr;
}

/**
 * @brief Sets a region of device memory to a specific byte value.
 *
 * Wrapper around `cudaMemset` with optional error checking.
 *
 * @tparam T The data type of the array elements (used for size calculation).
 * @param d_arr Pointer to the device memory region to set.
 * @param size The number of elements of type T in the region.
 * @param value The integer byte value to set (e.g., 0).
 * @return T* The pointer to the device memory region (`d_arr`).
 * @note Uses `gpuErrchk` if `CHECK_CUDA_ERRORS` is true.
 */
template <typename T>
inline T* memset_within_device(T* d_arr, size_t size, int value)
{
	if constexpr (CHECK_CUDA_ERRORS)
	{
		gpuErrchk(cudaMemset(d_arr, value, size * sizeof(T)));
	}
	else
	{
		cudaMemset(d_arr, value, size * sizeof(T));
	}

	return d_arr;
}

/**
 * @brief Allocates device memory and sets it to a specific byte value.
 *
 * Combines allocation (`alloc_to_device`) and memory set (`memset_within_device`).
 *
 * @tparam T The data type of the array elements.
 * @param size The number of elements of type T to allocate and set.
 * @param value The integer byte value to set (e.g., 0).
 * @return T* Pointer to the newly allocated and initialized device memory.
 * @note Uses `gpuErrchk` indirectly if `CHECK_CUDA_ERRORS` is true.
 */
template <typename T>
inline T* memset_to_device(size_t size, int value)
{
	T* d_arr = alloc_to_device<T>(size);
	memset_within_device(d_arr, size, value);
	return d_arr;
}

/**
 * @brief Copies data from one device memory location to another.
 *
 * Wrapper around `cudaMemcpy` (DeviceToDevice) with optional error checking.
 *
 * @tparam T The data type of the array elements.
 * @param d_dest Pointer to the destination device memory.
 * @param d_src Pointer to the source device memory.
 * @param size The number of elements of type T to copy.
 * @return T* The pointer to the destination device memory (`d_dest`).
 * @note Uses `gpuErrchk` if `CHECK_CUDA_ERRORS` is true.
 */
template <typename T>
inline T* copy_within_device(T* d_dest, const T* d_src, size_t size)
{
	if constexpr (CHECK_CUDA_ERRORS)
	{
		gpuErrchk(cudaMemcpy(d_dest, d_src, size * sizeof(T), cudaMemcpyDeviceToDevice));
	}
	else
	{
		cudaMemcpy(d_dest, d_src, size * sizeof(T), cudaMemcpyDeviceToDevice);
	}

	return d_dest;
}

/**
 * @brief Synchronizes the host thread with the default CUDA device stream.
 *
 * Waits until all preceding commands in the default stream have completed.
 * Wrapper around `cudaDeviceSynchronize` with optional error checking.
 *
 * @note Uses `gpuErrchk` if `CHECK_CUDA_ERRORS` is true.
 */
inline void synchronize_device()
{
	if constexpr (CHECK_CUDA_ERRORS)
	{
		gpuErrchk(cudaDeviceSynchronize());
	}
	else
	{
		cudaDeviceSynchronize();
	}
}

/**
 * @brief Allocates host memory and copies data from device memory to it.
 *
 * Wrapper around `cudaMemcpy` (DeviceToHost) with optional error checking.
 * Allocates host memory using `new[]`.
 *
 * @tparam T The data type of the array elements.
 * @param d_arr Pointer to the source device memory.
 * @param size The number of elements of type T to copy.
 * @return T* Pointer to the newly allocated host memory containing the copied data.
 * @note The caller is responsible for freeing the returned host memory using `delete[]`.
 * @note Uses `gpuErrchk` if `CHECK_CUDA_ERRORS` is true.
 */
template <typename T>
inline T* copy_from_device(const T* d_arr, size_t size)
{
	T* arr = new T[size]; // Allocate host memory

	if constexpr (CHECK_CUDA_ERRORS)
	{
		gpuErrchk(cudaMemcpy(arr, d_arr, size * sizeof(T), cudaMemcpyDeviceToHost));
	}
	else
	{
		cudaMemcpy(arr, d_arr, size * sizeof(T), cudaMemcpyDeviceToHost);
	}

	return arr;
}

/**
 * @brief Frees memory on the CUDA device.
 *
 * Wrapper around `cudaFree` with optional error checking.
 *
 * @tparam T The data type of the array elements (used for pointer type).
 * @param d_arr Pointer to the device memory to free.
 * @note Uses `gpuErrchk` if `CHECK_CUDA_ERRORS` is true.
 */
template <typename T>
inline void device_free(T* d_arr)
{
	if constexpr (CHECK_CUDA_ERRORS)
	{
		gpuErrchk(cudaFree(d_arr));
	}
	else
	{
		cudaFree(d_arr);
	}
}

/**
 * @brief Gets the global thread ID within a CUDA kernel launch.
 *
 * Calculates the unique ID for the current thread across all blocks in the grid.
 * Assumes a 1D grid and 1D block configuration.
 *
 * @return size_t The global thread identifier.
 * @note Device code only.
 */
__device__ __forceinline__ size_t get_tid()
{
	return threadIdx.x + blockIdx.x * blockDim.x;
}

// Very marginally slower (about 0.5% than get_num_threads() below)
// __device__ __forceinline__ size_t get_num_threads() {
//   return blockDim.x * gridDim.x;
// }


/**
 * @brief Gets the total number of threads in the kernel launch (compile-time).
 * @details Calculates the total number of threads based on the predefined
 * `BLOCKS` and `THREADS` macros.
 * @return constexpr size_t The total number of threads.
 * @note This is a `constexpr` function, evaluated at compile time.
 */
constexpr size_t get_num_threads()
{
	return BLOCKS * THREADS;
}
