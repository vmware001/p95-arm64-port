#include <stddef.h>
#include "gwnum.h"

static const unsigned int VERSION_NUMBER1 =
	GWNUM_MAJOR_VERSION * 100 + GWNUM_MINOR_VERSION;

#define VIS_INTERN __attribute__((visibility("internal")))
VIS_INTERN extern int xjmptable[];
VIS_INTERN extern int xjmptablep[];
VIS_INTERN extern int yjmptable[];
VIS_INTERN extern int yjmptablep[];
VIS_INTERN extern int zjmptable[];
VIS_INTERN extern int zjmptablep[];

struct gwinfo1s {
	void *d0,*d1,*d2,*d3,*d4,*d5,*d6,*d7;
	unsigned int v;
};

void gwinfo1(struct gwinfo1s *c) {
	c->d0 = xjmptable;
	c->d1 = xjmptablep;
	c->d2 = 0;
	c->d3 = 0;
	c->d4 = yjmptable;
	c->d5 = yjmptablep;
	c->d6 = zjmptable;
	c->d7 = zjmptablep;
	c->v = VERSION_NUMBER1;
}
