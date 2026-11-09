#include <stdlib.h>
#include <ctype.h>

#include "../gwnum.h"
#include "../gwtables.h"
#include "../gwnum_func.h"
#include "../gwcommon.h"
#include "../unravel.h"
#include "../xarch.h"
#include "../xbasics.h"
#include "xnormal.h"

#include "xmult1a.h"


void non2zpdn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3) {
	//puts("non2zpdn");
	uintptr_t rax=0, rbp, rdi, rsi;
	rsi = (uintptr_t)g->DESTARG;
	rdi = (uintptr_t)g->norm_biglit_array;
	rbp = (uintptr_t)g->norm_col_mults;
	xnorm012_1d_zpad(noexec, noexec);
	//cmnend;
}
void non2zpcdn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3) {
	//puts("non2zpcdn");
	uintptr_t rax=0, rbp, rdi, rsi;
	rsi = (uintptr_t)g->DESTARG;
	rdi = (uintptr_t)g->norm_biglit_array;
	rbp = (uintptr_t)g->norm_col_mults;
	xnorm012_1d_zpad(exec, noexec);
	//cmnend;
}
void zpcdn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3) {
	//puts("zpcdn");
	uintptr_t rax=0, rbp, rdi, rsi;
	rsi = (uintptr_t)g->DESTARG;
	rdi = (uintptr_t)g->norm_biglit_array;
	rbp = (uintptr_t)g->norm_col_mults;
	xnorm012_1d_zpad(exec, exec);
	//cmnend;
}
void zpdn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3) {
	//puts("zpdn");
	uintptr_t rax=0, rbp, rdi, rsi;
	rsi = (uintptr_t)g->DESTARG;
	rdi = (uintptr_t)g->norm_biglit_array;
	rbp = (uintptr_t)g->norm_col_mults;
	xnorm012_1d_zpad(noexec, exec);
	//cmnend;
}
void non2dn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3) {
	//puts("non2dn");
	uintptr_t rax, rcx, rbp, rdi, rsi;
	rsi = (uintptr_t)g->DESTARG;
	xnorm_top_carry_1d;				// Adjust top carry when k > 1
	rdi = (uintptr_t)g->norm_biglit_array;
	rbp = (uintptr_t)g->norm_col_mults;
	xnorm012_1d(noexec);
	//cmnend;
}
void idn(struct gwasm_data *__restrict g, vec2f64 xmm2, vec2f64 xmm3) {
	//printf("idn xmm2 %f %f xmm3 %f %f MAXERR %f\n", xmm2[0], xmm2[1], xmm3[0], xmm3[1], g->MAXERR);
	uintptr_t rax, rcx, rbp, rdi, rsi;
	rsi = (uintptr_t)g->DESTARG;	// Address of squared number
	xnorm_top_carry_1d;				// Adjust top carry when k > 1
	rdi = (uintptr_t)g->norm_biglit_array; // Address of the big/little flags array
	rbp = (uintptr_t)g->norm_col_mults; // Restart the column multipliers
	xnorm012_1d(exec);		// Add in carries
	//cmnend;
}


inorm(xr1, noexec, noexec, noexec, exec, noexec)
inorm(xr1e, noexec, exec, noexec, exec, noexec)
inorm(xr1c, noexec, noexec, exec, exec, noexec)
inorm(xr1ec, noexec, exec, exec, exec, noexec)
inorm(xi1, exec, noexec, noexec, exec, noexec)
inorm(xi1e, exec, exec, noexec, exec, noexec)
inorm(xi1c, exec, noexec, exec, exec, noexec)
inorm(xi1ec, exec, exec, exec, exec, noexec)
zpnorm(xr1zp, noexec, noexec, noexec, exec, noexec, exec, noexec, noexec)
zpnorm(xr1zpc1, noexec, noexec, noexec, exec, noexec, exec, exec, noexec)
zpnorm(xr1zpcm1, noexec, noexec, noexec, exec, noexec, exec, noexec, exec)
zpnorm(xr1zpe, noexec, exec, noexec, exec, noexec, exec, noexec, noexec)
zpnorm(xr1zpec1, noexec, exec, noexec, exec, noexec, exec, exec, noexec)
zpnorm(xr1zpecm1, noexec, exec, noexec, exec, noexec, exec, noexec, exec)
zpnorm(xr1zpc, noexec, noexec, exec, exec, noexec, exec, noexec, noexec)
zpnorm(xr1zpec, noexec, exec, exec, exec, noexec, exec, noexec, noexec)
zpnorm(xi1zp, exec, noexec, noexec, exec, noexec, exec, noexec, noexec)
zpnorm(xi1zpc1, exec, noexec, noexec, exec, noexec, exec, exec, noexec)
zpnorm(xi1zpcm1, exec, noexec, noexec, exec, noexec, exec, noexec, exec)
zpnorm(xi1zpe, exec, exec, noexec, exec, noexec, exec, noexec, noexec)
zpnorm(xi1zpec1, exec, exec, noexec, exec, noexec, exec, exec, noexec)
zpnorm(xi1zpecm1, exec, exec, noexec, exec, noexec, exec, noexec, exec)
zpnorm(xi1zpc, exec, noexec, exec, exec, noexec, exec, noexec, noexec)
zpnorm(xi1zpec, exec, exec, exec, exec, noexec, exec, noexec, noexec)
zpnorm(xr1zpk, noexec, noexec, noexec, exec, noexec, noexec, noexec, noexec)
zpnorm(xr1zpkc1, noexec, noexec, noexec, exec, noexec, noexec, exec, noexec)
zpnorm(xr1zpkcm1, noexec, noexec, noexec, exec, noexec, noexec, noexec, exec)
zpnorm(xr1zpek, noexec, exec, noexec, exec, noexec, noexec, noexec, noexec)
zpnorm(xr1zpekc1, noexec, exec, noexec, exec, noexec, noexec, exec, noexec)
zpnorm(xr1zpekcm1, noexec, exec, noexec, exec, noexec, noexec, noexec, exec)
zpnorm(xr1zpck, noexec, noexec, exec, exec, noexec, noexec, noexec, noexec)
zpnorm(xr1zpeck, noexec, exec, exec, exec, noexec, noexec, noexec, noexec)
zpnorm(xi1zpk, exec, noexec, noexec, exec, noexec, noexec, noexec, noexec)
zpnorm(xi1zpkc1, exec, noexec, noexec, exec, noexec, noexec, exec, noexec)
zpnorm(xi1zpkcm1, exec, noexec, noexec, exec, noexec, noexec, noexec, exec)
zpnorm(xi1zpek, exec, exec, noexec, exec, noexec, noexec, noexec, noexec)
zpnorm(xi1zpekc1, exec, exec, noexec, exec, noexec, noexec, exec, noexec)
zpnorm(xi1zpekcm1, exec, exec, noexec, exec, noexec, noexec, noexec, exec)
zpnorm(xi1zpck, exec, noexec, exec, exec, noexec, noexec, noexec, noexec)
zpnorm(xi1zpeck, exec, exec, exec, exec, noexec, noexec, noexec, noexec)

inorm(xr1b, noexec, noexec, noexec, noexec, noexec)
inorm(xr1eb, noexec, exec, noexec, noexec, noexec)
inorm(xr1cb, noexec, noexec, exec, noexec, noexec)
inorm(xr1ecb, noexec, exec, exec, noexec, noexec)
inorm(xi1b, exec, noexec, noexec, noexec, noexec)
inorm(xi1eb, exec, exec, noexec, noexec, noexec)
inorm(xi1cb, exec, noexec, exec, noexec, noexec)
inorm(xi1ecb, exec, exec, exec, noexec, noexec)
zpnorm(xr1zpb, noexec, noexec, noexec, noexec, noexec, exec, noexec, noexec)
zpnorm(xr1zpbc1, noexec, noexec, noexec, noexec, noexec, exec, exec, noexec)
zpnorm(xr1zpbcm1, noexec, noexec, noexec, noexec, noexec, exec, noexec, exec)
zpnorm(xr1zpeb, noexec, exec, noexec, noexec, noexec, exec, noexec, noexec)
zpnorm(xr1zpebc1, noexec, exec, noexec, noexec, noexec, exec, exec, noexec)
zpnorm(xr1zpebcm1, noexec, exec, noexec, noexec, noexec, exec, noexec, exec)
zpnorm(xr1zpcb, noexec, noexec, exec, noexec, noexec, exec, noexec, noexec)
zpnorm(xr1zpecb, noexec, exec, exec, noexec, noexec, exec, noexec, noexec)
zpnorm(xi1zpb, exec, noexec, noexec, noexec, noexec, exec, noexec, noexec)
zpnorm(xi1zpbc1, exec, noexec, noexec, noexec, noexec, exec, exec, noexec)
zpnorm(xi1zpbcm1, exec, noexec, noexec, noexec, noexec, exec, noexec, exec)
zpnorm(xi1zpeb, exec, exec, noexec, noexec, noexec, exec, noexec, noexec)
zpnorm(xi1zpebc1, exec, exec, noexec, noexec, noexec, exec, exec, noexec)
zpnorm(xi1zpebcm1, exec, exec, noexec, noexec, noexec, exec, noexec, exec)
zpnorm(xi1zpcb, exec, noexec, exec, noexec, noexec, exec, noexec, noexec)
zpnorm(xi1zpecb, exec, exec, exec, noexec, noexec, exec, noexec, noexec)
zpnorm(xr1zpbk, noexec, noexec, noexec, noexec, noexec, noexec, noexec, noexec)
zpnorm(xr1zpbkc1, noexec, noexec, noexec, noexec, noexec, noexec, exec, noexec)
zpnorm(xr1zpbkcm1, noexec, noexec, noexec, noexec, noexec, noexec, noexec, exec)
zpnorm(xr1zpebk, noexec, exec, noexec, noexec, noexec, noexec, noexec, noexec)
zpnorm(xr1zpebkc1, noexec, exec, noexec, noexec, noexec, noexec, exec, noexec)
zpnorm(xr1zpebkcm1, noexec, exec, noexec, noexec, noexec, noexec, noexec, exec)
zpnorm(xr1zpcbk, noexec, noexec, exec, noexec, noexec, noexec, noexec, noexec)
zpnorm(xr1zpecbk, noexec, exec, exec, noexec, noexec, noexec, noexec, noexec)
zpnorm(xi1zpbk, exec, noexec, noexec, noexec, noexec, noexec, noexec, noexec)
zpnorm(xi1zpbkc1, exec, noexec, noexec, noexec, noexec, noexec, exec, noexec)
zpnorm(xi1zpbkcm1, exec, noexec, noexec, noexec, noexec, noexec, noexec, exec)
zpnorm(xi1zpebk, exec, exec, noexec, noexec, noexec, noexec, noexec, noexec)
zpnorm(xi1zpebkc1, exec, exec, noexec, noexec, noexec, noexec, exec, noexec)
zpnorm(xi1zpebkcm1, exec, exec, noexec, noexec, noexec, noexec, noexec, exec)
zpnorm(xi1zpcbk, exec, noexec, exec, noexec, noexec, noexec, noexec, noexec)
zpnorm(xi1zpecbk, exec, exec, exec, noexec, noexec, noexec, noexec, noexec)
