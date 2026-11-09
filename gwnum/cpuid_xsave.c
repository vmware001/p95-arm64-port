#include <x86intrin.h>
#include <immintrin.h>
#include <cpuid.h>
#include "cpuid.h"

void exgetbv (struct cpuid_data *t) {
	uint64_t a = _xgetbv(t->ECX);
	t->EAX = (uint32_t)a;
	t->EDX = a >> 32;
}
