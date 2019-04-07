#if __x__64__
	#define SIZE_T ulong
#else
	#define SIZE_T uint
#endif

#define MATRIX_MULTIPLY(type_name)                  \
__kernel void matrix_multiply_ ## type_name(        \
	__global const type_name *A,                    \
	__global const type_name *B,                    \
	__global type_name *C,                          \
	const SIZE_T colsA,                              \
	const SIZE_T colsB                               \
) {                                                 \
	const SIZE_T i = get_global_id(0);               \
	const SIZE_T j = get_global_id(1);               \
                                                    \
	type_name acc = 0.0;                            \
	for (SIZE_T k = 0; k < colsA; ++k) {             \
		acc += A[i * colsA + k] * B[k * colsB + j]; \
	}                                               \
                                                    \
	C[i * colsB + j] = acc;                         \
}

#define MATRIX_SCALAR_MULTIPLY(type_name)           \
__kernel void matrix_scalar_multiply_ ## type_name( \
	__global const type_name *A,                    \
	__global type_name *C,                          \
	const type_name scalar                          \
) {                                                 \
	SIZE_T i = get_global_id(0);                    \
	C[i] = A[i] * scalar;                           \
}


MATRIX_MULTIPLY(float)
MATRIX_MULTIPLY(double)

MATRIX_SCALAR_MULTIPLY(float)
MATRIX_SCALAR_MULTIPLY(double)
