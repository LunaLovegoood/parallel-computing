#if __x__64__
	#define SIZE_T ulong
#else
	#define SIZE_T uint
#endif

#define NEGATIVE(type_name)           \
__kernel void negative_ ## type_name( \
	__global const type_name *A,      \
	__global type_name *negativeA     \
) {                                   \
	SIZE_T i = get_global_id(0);      \
	negativeA[i] = -A[i];             \
}


NEGATIVE(float)
NEGATIVE(double)
