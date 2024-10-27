#pragma once
#include <bit> // bit_width
#include <cstddef> // size_t
#include <cstdint> // uintmax_t
#include <cstdio> // fprintf
#include <cuda_runtime.h> // CUDA runtime

#define BLOCKS 64
#define THREADS 1024

#define gpuErrchk(ans)                                                         \
	{ gpuAssert(ans, __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
						line);
		if (abort)
			exit(code);
	}
}

static constexpr bool CHECK_CUDA_ERRORS = false;

__host__ inline size_t host_fast_log2(size_t n)
{
	return std::bit_width(n) - 1;
}

__host__ inline size_t host_fast_log2_ceil(size_t n)
{
	return std::bit_width(n - 1);
}

__device__ __forceinline__ size_t device_fast_log2(size_t n)
{
	return std::numeric_limits<size_t>::digits - 1 - __clzll(n);
}

__device__ __forceinline__ size_t device_fast_log2_ceil(size_t n)
{
	return std::numeric_limits<size_t>::digits - __clzll(n - 1);
}

class CudaTimer
{
	public:
		CudaTimer()
		{
			cudaEventCreate(&start);
			cudaEventCreate(&stop);
		}

		~CudaTimer()
		{
			cudaEventDestroy(start);
			cudaEventDestroy(stop);
		}

		void start_timer()
		{
			cudaEventRecord(start);
		}

		float stop_timer()
		{
			cudaEventRecord(stop);
			cudaEventSynchronize(stop);
			cudaEventElapsedTime(&time, start, stop);
			return time;
		}

	private:
		cudaEvent_t start, stop;
		float time;
};

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

template <typename T>
inline T* copy_to_device(const T* arr, size_t size)
{
	T* d_arr = alloc_to_device<T>(size);

	copy_to_device(d_arr, arr, size);

	return d_arr;
}

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

template <typename T>
inline T* memset_to_device(size_t size, int value)
{
	T* d_arr = alloc_to_device<T>(size);

	memset_within_device(d_arr, size, value);

	return d_arr;
}

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

template <typename T>
inline T* copy_from_device(const T* d_arr, size_t size)
{
	T* arr = new T[size];

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

__device__ __forceinline__ size_t get_tid()
{
	return threadIdx.x + blockIdx.x * blockDim.x;
}

// __device__ __forceinline__ size_t get_num_threads() {
//   return blockDim.x * gridDim.x;
// }


// this gives a very small speedup, about 0.5%
constexpr size_t get_num_threads()
{
	return BLOCKS * THREADS;
}
