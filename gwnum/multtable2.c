#include "multtable2.h"

#define STAT(A,B,F)
#define ENT1(A,P1,B) \
void P1(void *);
#define ENT2(A,P1,P2,B) \
void P1(void *); \
void P2(void *);
#define D8(A,B,C,D,E,F,G,H)
#define DZERO

LIST_xjmptable
LIST_xjmptablep
LIST_yjmptable
LIST_yjmptablep
LIST_zjmptable
LIST_zjmptablep



#undef STAT
#undef ENT1
#undef ENT2
#undef D8
#undef DZERO

typedef void (*func1)(void *);

#define VAR2(A) a##A
#define VAR(A) VAR2(A)

#define STAT(A,B,F) \
	int VAR(__COUNTER__); \
	int VAR(__COUNTER__); \
	float VAR(__COUNTER__);
#define ENT1(A,P1,B) \
	unsigned int VAR(__COUNTER__); \
	func1 VAR(__COUNTER__); \
	int VAR(__COUNTER__);

#define ENT2(A,P1,P2,B) \
	unsigned int VAR(__COUNTER__); \
	func1 VAR(__COUNTER__); \
	func1 VAR(__COUNTER__); \
	int VAR(__COUNTER__);

#define D8(A,B,C,D,E,F,G,H) int VAR(__COUNTER__) [8];
#define DZERO int VAR(__COUNTER__);

// Use packed struct to prevent padding.
#pragma pack(push, 1)

struct S_xjmptable {
	LIST_xjmptable
};
struct S_xjmptablep {
	LIST_xjmptablep
};
struct S_yjmptable {
	LIST_yjmptable
};
struct S_yjmptablep {
	LIST_yjmptablep
};
struct S_zjmptable {
	LIST_zjmptable
};
struct S_zjmptablep {
	LIST_zjmptablep
};

#pragma pack(pop)

#undef STAT
#undef ENT1
#undef ENT2
#undef D8
#undef DZERO

#define STAT(A,B,F) A,B,F,
#define ENT1(A,P1,B) A,P1,B,
#define ENT2(A,P1,P2,B) A,P1,P2,B,
#define D8(A,B,C,D,E,F,G,H) {A,B,C,D,E,F,G,H},
#define DZERO 0,

const struct S_xjmptable xjmptable = {
	LIST_xjmptable
};
const struct S_xjmptablep xjmptablep = {
	LIST_xjmptablep
};
const struct S_yjmptable yjmptable = {
	LIST_yjmptable
};
const struct S_yjmptablep yjmptablep = {
	LIST_yjmptablep
};
const struct S_zjmptable zjmptable = {
	LIST_zjmptable
};
const struct S_zjmptablep zjmptablep = {
	LIST_zjmptablep
};

