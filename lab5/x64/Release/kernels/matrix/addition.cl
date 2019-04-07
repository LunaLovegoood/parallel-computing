#if __x__64__
	#define SIZE_T ulong
#else
	#define SIZE_T uint
#endif

#define MATRIX_ADD(type_name)           \
__kernel void matrix_add_ ## type_name( \
	__global const type_name *A,        \
	__global const type_name *B,        \
	__global type_name *C               \
) {                                     \
	SIZE_T i = get_global_id(0);        \
	C[i] = A[i] + B[i];                 \
}

#define MATRIX_SCALAR_ADD(type_name)           \
__kernel void matrix_scalar_add_ ## type_name( \
	__global const type_name *A,               \
	__global type_name *C,                     \
	const type_name scalar                     \
) {                                            \
	SIZE_T i = get_global_id(0);               \
	C[i] = A[i] + scalar;                      \
}


MATRIX_ADD(float)
MATRIX_ADD(double)

MATRIX_SCALAR_ADD(float)
MATRIX_SCALAR_ADD(double)
