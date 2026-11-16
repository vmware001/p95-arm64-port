/* Hand written AArch64 NEON assembly for the hottest inner loops of
 * gwnum's radix-4 pass 1 FFT (the "r4_pass1sc64" kernels, clm=4).
 *
 * The converted-to-C original performs, per call,
 *   for z in 0..z_cnt-1:
 *     for y in 0..mid_cnt-1:
 *       copy 512 bytes src -> dst      (2*clm 64-byte cache lines)
 *       src += src_mid; dst += dst_mid        (net advances)
 *     src += src_z; dst += dst_z
 * which is a pure cache-line swizzle between the pass 1 scratch area and
 * the FFT data area (bit exact by construction: no arithmetic at all).
 * clang -O3 issues ldp/stp pairs with an immediate load->store dependency,
 * so on the in-order Cortex-A53 every store waits for its load.  Here the
 * loads of the next 128 bytes are issued before the stores of the current
 * 128 bytes, keeping 8 loads in flight (v0-v15).
 *
 * void r4_pass1_copy_asm(void *src, void *dst,
 *                        long src_mid, long dst_mid,
 *                        long src_z, long dst_z,
 *                        long mid_cnt, long z_cnt);
 *   x0=src x1=dst x2=src_mid x3=dst_mid x4=src_z x5=dst_z x6=mid_cnt x7=z_cnt
 */
	.text
	.globl	r4_pass1_copy_asm
	.type	r4_pass1_copy_asm,%function
	.align	4
r4_pass1_copy_asm:
	cbz	x7, .Lcopy_done		/* z_cnt == 0 -> nothing */
	cbz	x6, .Lcopy_done		/* mid_cnt == 0 -> nothing */
	mov	x9, x0			/* src cursor */
	mov	x10, x1			/* dst cursor */
	sub	x2, x2, #512		/* extra advance after the 512 copied bytes */
	sub	x3, x3, #512
.Lcopy_z:
	mov	x11, x6			/* mid counter */
.Lcopy_mid:
	/* ---- 0 .. 255 : issue all 8 loads, then all 8 stores ---- */
	ldp	q0, q1, [x9]
	ldp	q2, q3, [x9, #0x20]
	ldp	q4, q5, [x9, #0x40]
	ldp	q6, q7, [x9, #0x60]
	ldp	q8, q9, [x9, #0x80]
	ldp	q10, q11, [x9, #0xa0]
	ldp	q12, q13, [x9, #0xc0]
	ldp	q14, q15, [x9, #0xe0]
	add	x9, x9, #0x100
	stp	q0, q1, [x10]
	stp	q2, q3, [x10, #0x20]
	stp	q4, q5, [x10, #0x40]
	stp	q6, q7, [x10, #0x60]
	stp	q8, q9, [x10, #0x80]
	stp	q10, q11, [x10, #0xa0]
	stp	q12, q13, [x10, #0xc0]
	stp	q14, q15, [x10, #0xe0]
	add	x10, x10, #0x100
	/* ---- 256 .. 511 ---- */
	ldp	q0, q1, [x9]
	ldp	q2, q3, [x9, #0x20]
	ldp	q4, q5, [x9, #0x40]
	ldp	q6, q7, [x9, #0x60]
	ldp	q8, q9, [x9, #0x80]
	ldp	q10, q11, [x9, #0xa0]
	ldp	q12, q13, [x9, #0xc0]
	ldp	q14, q15, [x9, #0xe0]
	add	x9, x9, #0x100
	stp	q0, q1, [x10]
	stp	q2, q3, [x10, #0x20]
	stp	q4, q5, [x10, #0x40]
	stp	q6, q7, [x10, #0x60]
	stp	q8, q9, [x10, #0x80]
	stp	q10, q11, [x10, #0xa0]
	stp	q12, q13, [x10, #0xc0]
	stp	q14, q15, [x10, #0xe0]
	add	x10, x10, #0x100
	/* advance by the mid strides (the 512 copied bytes are already counted) */
	add	x9, x9, x2
	add	x10, x10, x3
	subs	x11, x11, #1
	b.ne	.Lcopy_mid
	/* end of one z iteration */
	add	x9, x9, x4
	add	x10, x10, x5
	subs	x7, x7, #1
	b.ne	.Lcopy_z
.Lcopy_done:
	ret
	.size	r4_pass1_copy_asm, .-r4_pass1_copy_asm
