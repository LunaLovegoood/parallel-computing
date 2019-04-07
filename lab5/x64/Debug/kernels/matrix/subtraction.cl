#if __x__64__
	#define SIZE_T ulong
#else
	#define SIZE_T uint
#endif

#define MATRIX_SUBTRACT(type_name)           \
__kernel void matrix_subtract_ ## type_name( \
	__global const type_name *A,             \
	__global const type_name *B,             \
	__global type_name *C                    \
) {                                          \
	SIZE_T i = get_global_id(0);             \
	C[i] = A[i] - B[i];                      \
}


MATRIX_SUBTRACT(float)
MATRIX_SUBTRACT(double)
