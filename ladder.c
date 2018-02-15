#include <stdio.h>
#include <stdlib.h>

#include "gmp.h"

#define size 32

//mpn_cnd_add_n ...

void xDBL(mpz_t N, mpz_t A, mpz_t X, mpz_t Z, mpz_t Xm, mpz_t Zm)
// Q * R * S * (R + S * (A + 2) / 4)
// Q = (X + Z)²
// R = (X - Z)²
// S = Q - R
{
	mpz_t Q, R, S;
	mp_limb_t* tmp;
	mp_limb_t* scratch;
	mpz_inits(Q, R, S, NULL);
	
	mp_limb_t* limbs_X = mpz_limbs_read(X);
	mp_limb_t* limbs_Z = mpz_limbs_read(Z);
	mp_limb_t* limbs_A = mpz_limbs_read(A);
	mp_limb_t* limbs_N = mpz_limbs_read(N);
	mp_limb_t* limbs_Q = mpz_limbs_write(Q, size * 2);
	mp_limb_t* limbs_R = mpz_limbs_write(R, size * 2);
	mp_limb_t* limbs_S = mpz_limbs_write(S, size * 2);
	
	
//	printf(".");
	gmp_printf("Q = %Zd\n", Q);
	
	mpn_add_n(limbs_Q, limbs_X, limbs_Z, size + 1);
	mpz_limbs_finish(Q, size + 1);
	
//	printf(".");
	gmp_printf("Q = %Zd\n", Q);
	
	scratch = malloc(mpn_sec_div_r_itch(size + 1, size) * sizeof(mp_limb_t));
	mpn_sec_div_r(limbs_Q, size + 1, limbs_N, size, scratch);
	mpz_limbs_finish(Q, size * 2);
	free(scratch);
	
//	printf(".");
	gmp_printf("Q = %Zd\n", Q);
	
	limbs_Q = mpz_limbs_modify(Q, size * 2);
	tmp = mpz_limbs_read(Q);
	scratch = malloc(mpn_sec_sqr_itch(size * 2) * sizeof(mp_limb_t));
	mpn_sec_sqr(limbs_Q, tmp, size, scratch);
	free(scratch);
	
//	printf(".");
	gmp_printf("Q = %Zd\n", Q);
	gmp_printf("N = %Zd\n", N);
	
	scratch = malloc(mpn_sec_div_r_itch(size * 2, size) * sizeof(mp_limb_t));
	mpn_sec_div_r(limbs_Q, size * 2, limbs_N, size, scratch);
	mpz_limbs_finish(Q, size * 2);
	free(scratch);
	
//	printf(".");
	gmp_printf("Q = %Zd\n", Q);
	
	mpz_set(Xm, Q);
	
//	printf(".");
	
	mpz_clears(Q, R, S, NULL);
}

void ladder(mpz_t N, mpz_t A, mpz_t m, mpz_t X, mpz_t Z, mpz_t Xm, mpz_t Zm)
{
	
}

int main()
{
	mpz_t n, a, x, z, xm, zm;
	mpz_inits(n, a, x, z, xm, zm, NULL);
	
//	printf(".");
	
	mpz_set_ui(n, 256);
	mpz_set_ui(a, 0);
	mpz_set_ui(x, 320);
	mpz_set_ui(z, 10);
	
//	printf(".");
	
	mpz_realloc(n, size);
	mpz_realloc(a, size);
	mpz_realloc(x, size);
	mpz_realloc(z, size);
	
//	printf(".");
	
	xDBL(n, a, x, z, xm, zm);
	
//	printf(".");
	
	gmp_printf("dbl = %Zd\n", xm);
	
	mpz_clears(n, a, x, z, xm, zm, NULL);
	
	return 0;
}
