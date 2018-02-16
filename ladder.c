#include <stdio.h>
#include <stdlib.h>

#include "gmp.h"

int size = 32 / sizeof(const long unsigned int);

//mpn_cnd_add_n ...

void sec_div_r(mpz_t a, mpz_t b, int size_b)
{
	const mp_limb_t* limbs_b;
	mp_limb_t* limbs_a;
	mp_limb_t* scratch;
	int size_a = mpz_size(a);
	limbs_b = mpz_limbs_read(b);
	limbs_a = mpz_limbs_modify(a, size_a);
	scratch = malloc(mpn_sec_div_r_itch(size_a, size_b) * sizeof(mp_limb_t));
	mpn_sec_div_r(limbs_a, size_a, limbs_b, size_b, scratch);
	mpz_limbs_finish(a, size_b);
	free(scratch); 
}

void sec_div_q(mpz_t a, mpz_t b, int size_b)
{
	mpz_t r;
	mpz_init(r);
	const mp_limb_t* limbs_b;
	mp_limb_t* limbs_a;
	mp_limb_t* scratch;
	mp_limb_t* limbs_r;
	mp_limb_t tmp;
	int size_a = mpz_size(a);
	limbs_b = mpz_limbs_read(b);
	limbs_a = mpz_limbs_modify(a, size_a);
	limbs_r = mpz_limbs_write(r, size_a);
	scratch = malloc(mpn_sec_div_qr_itch(size_a, size_b) * sizeof(mp_limb_t));
	tmp = mpn_sec_div_qr(limbs_a, limbs_r, size_a, limbs_b, size_b, scratch);
	mpz_limbs_finish(a, size_a - size_b);
	free(scratch);
	mpz_clear(r);
}

void sec_sqr(mpz_t a)
{
	mp_limb_t* limbs_a;
	mp_limb_t* scratch;
	const mp_limb_t* tmp;
	int size_a;
	size_a = mpz_size(a);
	limbs_a = mpz_limbs_modify(a, size_a * 2);
	tmp = mpz_limbs_read(a);
	scratch = malloc(mpn_sec_sqr_itch(size_a) * sizeof(mp_limb_t));
	mpn_sec_sqr(limbs_a, tmp, size_a, scratch);
	mpz_limbs_finish(a, size_a * 2);
	free(scratch);
}

void sec_mul(mpz_t a, mpz_t b, mpz_t res)
{
	const mp_limb_t *limbs_a, *limbs_b;
	mp_limb_t *limbs_res, *scratch;
	int size_a;
	size_a = mpz_size(a);
	limbs_a = mpz_limbs_read(a);
	limbs_b = mpz_limbs_read(b);
	limbs_res = mpz_limbs_write(res, size_a * 2);
	scratch = malloc(mpn_sec_mul_itch(size_a, size_a));
	mpn_sec_mul(limbs_res, limbs_a, size_a, limbs_b, size_a, scratch);
	mpz_limbs_finish(res, size_a * 2);
}

void sec_sub_mod(mpz_t a, mpz_t b, mpz_t res, mpz_t n)
{
	mpz_t tmp1, tmp2;
	const mp_limb_t *limbs_a, *limbs_b, *limbs_n;
	mp_limb_t *limbs_res, *limbs_tmp1, *limbs_tmp2, *scratch;
	int sizea;
	gmp_printf("a = %Zd\nb = %Zd\n", a, b);
	mpz_inits(tmp1, tmp2, NULL);
	sizea = mpz_size(a);
	printf("sizea = %d\n", sizea);
	limbs_a = mpz_limbs_read(a);
	limbs_b = mpz_limbs_read(b);
	limbs_n = mpz_limbs_read(n);
	limbs_res = mpz_limbs_write(res, sizea);
	mpz_limbs_finish(res, size);
	limbs_tmp1 = mpz_limbs_write(tmp1, sizea);
	limbs_tmp2 = mpz_limbs_write(tmp2, sizea);
	mpn_sub_n(limbs_res, limbs_a, limbs_b, sizea);
	mpz_limbs_finish(res, sizea);
/*	mpn_add_n(limbs_tmp1, limbs_a, limbs_n, sizea);
	mpz_limbs_finish(tmp1, sizea);
	gmp_printf("tmp1 = %Zd\n", tmp1);
	mpn_sub_n(limbs_tmp2, limbs_tmp1, limbs_b, sizea);
	mpz_limbs_finish(tmp2, sizea);
	gmp_printf("tmp2 = %Zd\n", tmp2);
	sec_div_r(tmp2, n, size);
	gmp_printf("res = %Zd\ntmp2 = %Zd\n", res, tmp2);*/
	mpn_sub_n(limbs_tmp1, limbs_a, limbs_b, sizea);
	mpn_add_n(limbs_tmp2, limbs_tmp1, limbs_n, sizea);
	mpn_cnd_swap(mpz_cmp(a, b) < 0, limbs_res, limbs_tmp2, size);
	mpz_limbs_finish(res, size);
	mpz_clears(tmp1, tmp2, NULL);
}

void xDBL(mpz_t N, mpz_t A, mpz_t X, mpz_t Z, mpz_t Xm, mpz_t Zm)
{
	mpz_t Q, R, S, x, y;
	const mp_limb_t* tmp;
	mp_limb_t* scratch, *limbs_x, *limbs_y;
	mpz_inits(Q, R, S, x, y, NULL);
	
	const mp_limb_t* limbs_X = mpz_limbs_read(X);
	const mp_limb_t* limbs_Z = mpz_limbs_read(Z);
	const mp_limb_t* limbs_A = mpz_limbs_read(A);
	const mp_limb_t* limbs_N = mpz_limbs_read(N);
	mp_limb_t* limbs_Q = mpz_limbs_write(Q, size * 2);
	mp_limb_t* limbs_R = mpz_limbs_write(R, size * 2);
	mp_limb_t* limbs_S = mpz_limbs_write(S, size * 2);
	mp_limb_t* limbs_Xm = mpz_limbs_write(Xm, size * 2);
	mp_limb_t* limbs_Zm = mpz_limbs_write(Zm, size * 2);
	
	
	
	mpn_add_n(limbs_Q, limbs_X, limbs_Z, size + 1);
	mpz_limbs_finish(Q, size + 1);
	sec_div_r(Q, N, size);
	sec_sqr(Q);
	sec_div_r(Q, N, size);
	
	
	sec_sub_mod(X, Z, R, N);
	/*
	limbs_x = mpz_limbs_write(x, size);
	mpn_sub_n(limbs_R, limbs_X , limbs_Z, size);
	mpn_sub_n(limbs_x, limbs_Z, limbs_X, size);
	mpn_cnd_swap(mpz_cmp(X, Z) < 0, limbs_R, limbs_x, size);
	mpz_limbs_finish(R, size);
	sec_div_r(R, N, size);
	sec_sqr(R);
	sec_div_r(R, N, size);*/
	gmp_printf("R = %Zd\n", R);
	
	
	
	mpn_sub_n(limbs_S, limbs_Q, limbs_R, size);
	mpz_limbs_finish(S, size + 1);
	scratch = malloc(mpn_sec_div_r_itch(size + 1, size) * sizeof(mp_limb_t));
	mpn_sec_div_r(limbs_S, size + 1, limbs_N, size, scratch);
	mpz_limbs_finish(S, size);
	free(scratch);
	
	
	
	sec_mul(Q, R, Xm);
	sec_div_r(Xm, N, size);
	
	

	sec_mul(A, S, x);
	limbs_x = mpz_limbs_read(x);
	limbs_y = mpz_limbs_write(y, size + 2);
	mpn_rshift(limbs_y, limbs_x, size + 2, 2);
	mpz_limbs_finish(y, size);
	sec_div_r(y, N, size);
	limbs_x = mpz_limbs_write(x, size);
	limbs_y = mpz_limbs_read(y);
	limbs_R = mpz_limbs_read(R);
	mpn_add_n(limbs_x, limbs_y, limbs_R, size);
	mpz_limbs_finish(x, size);
	sec_mul(x, S, Zm);
	sec_div_r(Zm, N, size);
	
	
	
	mpz_clears(Q, R, S, x, y, NULL);
}

void xADD(mpz_t N, mpz_t A, mpz_t Xp, mpz_t Zp, mpz_t Xq, mpz_t Zq, mpz_t Xminus, mpz_t Zminus, mpz_t Xm, mpz_t Zm)
{
	mpz_t U, V, tmp1, tmp2;
	const mp_limb_t *limbs_Xp, *limbs_Xminus, *limbs_Zminus, *limbs_Zp, *limbs_Xq, *limbs_Zq, *limbs_A, *limbs_N;
	mp_limb_t *limbs_U, *limbs_V, *limbs_Xm, *limbs_Zm, *limbs_tmp1, *limbs_tmp2;
	
	mpz_inits(U, V, tmp1, tmp2, NULL);
	limbs_Xp = mpz_limbs_read(Xp);
	limbs_Zp = mpz_limbs_read(Zp);
	limbs_Xq = mpz_limbs_read(Xq);
	limbs_Zq = mpz_limbs_read(Zq);
	limbs_A = mpz_limbs_read(A);
	limbs_N = mpz_limbs_read(N);
	limbs_U = mpz_limbs_write(U, size * 2);
	limbs_V = mpz_limbs_write(V, size * 2);
	limbs_Xminus = mpz_limbs_write(Xminus, size * 2);
	limbs_Zminus = mpz_limbs_write(Zminus, size * 2);
	limbs_Xm = mpz_limbs_write(Xm, size * 2);
	limbs_Zm = mpz_limbs_write(Zm, size * 2);
	
	
	
	sec_sub_mod(Xp, Zp, tmp1, N);
/*	limbs_tmp1 = mpz_limbs_write(tmp1, size + 1);
	mpn_sub_n(limbs_tmp1, limbs_Xp, limbs_Zp, size + 1);
	mpn_sub_n(limbs_tmp2, limbs_Zp, limbs_Xp, size + 1);
	mpn_cnd_swap(mpz_cmp(Xp, Zp) < 0, limbs_tmp1, limbs_tmp2, size + 1);
	mpz_limbs_finish(tmp1, size + 1);*/
	mpz_realloc(tmp1, size + 1);
	
	
	
	limbs_tmp2 = mpz_limbs_write(tmp2, size + 1);
	mpn_add_n(limbs_tmp2, limbs_Xq, limbs_Zq, size + 1);
	mpz_limbs_finish(tmp2, size + 1);
	
	
	
	sec_mul(tmp1, tmp2, U);
	sec_div_r(U, N, size);
	
	
	
	
	
	limbs_tmp1 = mpz_limbs_write(tmp1, size + 1);
	mpn_add_n(limbs_tmp1, limbs_Xp, limbs_Zp, size + 1);
	mpz_limbs_finish(tmp1, size + 1);
	
	limbs_tmp2 = mpz_limbs_write(tmp2, size + 1);
	mpn_sub_n(limbs_tmp2, limbs_Xq, limbs_Zq, size + 1);
	mpz_limbs_finish(tmp2, size + 1);
	
	sec_mul(tmp1, tmp2, V);
	sec_div_r(V, N, size);
	
	gmp_printf("V = %Zd\n", V);
	
	
	
	limbs_tmp1 = mpz_limbs_write(tmp1, size + 1);
	limbs_U = mpz_limbs_read(U);
	limbs_V = mpz_limbs_read(V);
	mpn_add_n(limbs_tmp1, limbs_U, limbs_V, size + 1);
	mpz_limbs_finish(tmp1, size + 1);
	sec_sqr(tmp1);
	sec_div_r(tmp1, N, size);
	sec_mul(tmp1, Zminus, Xm);
	sec_div_r(Xm, N, size);

//	gmp_printf("sizeof U = %d\nsizeof V = %d\n", mpz_size(U), mpz_size(V));
	
	sec_sub_mod(U, V, tmp1, N);
	
	gmp_printf("U - V = %Zd - %Zd = %Zd\n", U, V, tmp1);
	
	sec_sqr(tmp1);
	sec_div_r(tmp1, N, size);
	sec_mul(tmp1, Xminus, Zm);
	sec_div_r(Zm, N, size);
	
	
	
	mpz_clears(U, V, tmp1, tmp2, NULL);
}

void ladder(mpz_t N, mpz_t A, mpz_t m, mpz_t X, mpz_t Z, mpz_t Xm, mpz_t Zm)
{
	
}

int main()
{
	mpz_t n, a, x1, z1, x2, z2, xm, zm, xu, zu;
	mpz_inits(n, a, x1, z1, x2, z2, xm, zm, xu, zu, NULL);
	
	mpz_set_ui(n, 243);
	size = mpz_size(n);
	mpz_realloc(a, size);
	mpz_realloc(x1, size);
	mpz_realloc(z1, size);
	mpz_realloc(x2, size);
	mpz_realloc(z2, size);
	mpz_realloc(xu, size);
	mpz_realloc(zu, size);
	mpz_set_ui(a, 6);
	mpz_set_ui(x1, 7);
	mpz_set_ui(z1, 15);
	mpz_set_ui(x2, 3);
	mpz_set_ui(z2, 4);
	mpz_set_ui(xu, 22);
	mpz_set_ui(zu, 1);
	
	mpz_add_ui(a, a, 2);
	
	xADD(n, a, x1, z1, x2, z2, xu, zu, xm, zm);
	
	gmp_printf("U = %Zd\nV = %Zd\n", xm, zm);
	
	mpz_clears(n, a, x1, z1, x2, z2, xm, zm, xu, zu, NULL);
	
	return 0;
}
