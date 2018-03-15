#include <stdio.h>
#include <stdlib.h>

#include "gmp.h"

int size = 32 / sizeof(const long unsigned int);

//mpn_cnd_add_n ...

void sec_div_r(mpz_t a, mpz_t b, int size_b)
{
	mpz_t a2;
	const mp_limb_t* limbs_b, *limbs_a2;
	mp_limb_t* limbs_a;
	mp_limb_t* scratch;
	int size_a = mpz_size(a);
	mpz_init(a2);
	mpz_set(a2, b);
	mpz_realloc(a2, size_a + 1);
	limbs_b = mpz_limbs_read(b);
	limbs_a = mpz_limbs_modify(a, size_a + 1);
	limbs_a2 = mpz_limbs_read(a2);
	mpn_add_n(limbs_a, limbs_a2, limbs_b, size_a + 1);
	mpz_limbs_finish(a, size_a + 1);
	scratch = malloc(mpn_sec_div_r_itch(size_a + 1, size_b) * sizeof(mp_limb_t));
	mpn_sec_div_r(limbs_a, size_a + 1, limbs_b, size_b, scratch);
	mpz_limbs_finish(a, size_b);
	free(scratch); 
	mpz_clear(a2);
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

void sec_mul(mpz_t a, mpz_t b, mpz_t res, int size_ops)
{
	const mp_limb_t *limbs_a, *limbs_b;
	mp_limb_t *limbs_res, *scratch;
	int size_a, size_b;
	mpz_realloc(a, size_ops);
	mpz_realloc(b, size_ops);
//	size_a = mpz_size(a);
//	size_b = mpz_size(b);
	limbs_a = mpz_limbs_read(a);
	limbs_b = mpz_limbs_read(b);
	limbs_res = mpz_limbs_write(res, size_ops * 2);
	scratch = malloc(mpn_sec_mul_itch(size_ops, size_ops) * sizeof(mp_limb_t));
	mpn_sec_mul(limbs_res, limbs_a, size_ops, limbs_b, size_ops, scratch);
	mpz_limbs_finish(res, size_ops * 2);
}

void sec_sub_mod(mpz_t a, mpz_t b, mpz_t res, mpz_t n)
{
	mpz_t tmp1, tmp2;
	const mp_limb_t *limbs_a, *limbs_b, *limbs_n;
	mp_limb_t *limbs_res, *limbs_tmp1, *limbs_tmp2, *scratch;
	int sizea;
	mpz_inits(tmp1, tmp2, NULL);
	sizea = mpz_size(a);
	limbs_a = mpz_limbs_read(a);
	limbs_b = mpz_limbs_read(b);
	limbs_n = mpz_limbs_read(n);
	limbs_res = mpz_limbs_write(res, sizea);
	mpz_limbs_finish(res, size);
	limbs_tmp1 = mpz_limbs_write(tmp1, sizea);
	limbs_tmp2 = mpz_limbs_write(tmp2, sizea);
	mpn_sub_n(limbs_res, limbs_a, limbs_b, sizea);
	mpz_limbs_finish(res, sizea);
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
	mp_limb_t* scratch, *limbs_x, *limbs_y, *op1, *op2;
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
	
	
	
	mpn_sub_n(limbs_S, limbs_Q, limbs_R, size);
	mpz_limbs_finish(S, size + 1);
	sec_div_r(S, N, size);
	
	
	
	sec_mul(Q, R, Xm, size);
	sec_div_r(Xm, N, size);
	
	
	
	sec_mul(S, A, x, size);
	sec_div_r(x, N, size);
	limbs_y = mpz_limbs_write(y, size);
	limbs_x = mpz_limbs_read(x);
	limbs_R = mpz_limbs_read(R);
	mpn_add_n(limbs_y, limbs_x, limbs_R, size);
	mpz_limbs_finish(y, size);
	sec_mul(y, S, Zm, size);
	sec_div_r(Zm, N, size);
	
	
	mpz_clears(Q, R, S, x, y, NULL);
}

void xADD(mpz_t N, mpz_t A, mpz_t Xp, mpz_t Zp, mpz_t Xq, mpz_t Zq, mpz_t Xminus, mpz_t Zminus, mpz_t Xm, mpz_t Zm)
{
	mpz_t U, V, tmp1, tmp2;
	const mp_limb_t *limbs_Xp, *limbs_Xminus, *limbs_Zminus, *limbs_Zp, *limbs_Xq, *limbs_Zq, *limbs_A, *limbs_N;
	mp_limb_t *limbs_U, *limbs_V, *limbs_Xm, *limbs_Zm, *limbs_tmp1, *limbs_tmp2;
	
	printf(".");
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
	
	printf(".");
	
	sec_sub_mod(Xp, Zp, tmp1, N);
//	mpz_realloc(tmp1, size + 1);
	printf(".");
	
	
	
	limbs_tmp2 = mpz_limbs_write(tmp2, size + 1);
	mpn_add_n(limbs_tmp2, limbs_Xq, limbs_Zq, size + 1);
	mpz_limbs_finish(tmp2, size + 1);
	printf(".");
	
	
	
	sec_mul(tmp1, tmp2, U, size + 1);
	sec_div_r(U, N, size);
	printf(".");
	
	
	
	
	
	limbs_tmp1 = mpz_limbs_write(tmp1, size + 1);
	mpn_add_n(limbs_tmp1, limbs_Xp, limbs_Zp, size + 1);
	mpz_limbs_finish(tmp1, size + 1);
	printf(".");
	
	sec_sub_mod(Xq, Zq, tmp2, N);
	sec_mul(tmp1, tmp2, V, size);
	sec_div_r(V, N, size);
	printf(".");
	


	sec_sub_mod(U, V, tmp1, N);
	printf(".");
	
	sec_sqr(tmp1);
	sec_div_r(tmp1, N, size);
	printf(".");
	
	
	
	sec_mul(Zminus, tmp1, Xm, size);
	sec_div_r(Xm, N, size);
	printf(".");


	
	limbs_tmp1 = mpz_limbs_write(tmp1, size + 1);
	limbs_U = mpz_limbs_read(U);
	limbs_V = mpz_limbs_read(V);
	mpn_add_n(limbs_tmp1, limbs_U, limbs_V, size + 1);
	mpz_limbs_finish(tmp1, size + 1);
	printf(".");
	
	
	
	sec_sqr(tmp1);
	sec_div_r(tmp1, N, size);
	sec_mul(tmp1, Xminus, Zm, size);
	sec_div_r(Zm, N, size);
	printf(".");
	
	
	
	mpz_clears(U, V, tmp1, tmp2, NULL);
}

int mod2(mpz_t m)
{
	mp_limb_t out, one;
	mp_limb_t *limbs_read, *limbs_write;
	
	one = 1 << sizeof(int) * 8 - 1;
	limbs_write = mpz_limbs_modify(m, size);
	
	mpn_and_n(&out, &one, limbs_write, 1);
	
	mpn_lshift(limbs_write, limbs_write, size, 1);
	mpz_limbs_finish(m, size);
	
	return out;
}

void ladder(mpz_t N, mpz_t A, mpz_t m, mpz_t X, mpz_t Z, mpz_t Xm, mpz_t Zm)
{
	int i, n;
	mpz_t X0, Z0, X1, Z1, Xu, Zu, xtmp1, xtmp2, ztmp1, ztmp2, xtmp3, ztmp3;
	mp_limb_t *limbs_tmp1, *limbs_tmp2, *limbs_tmp3;
	mpz_inits(X0, Z0, X1, Z1, Xu, Zu, xtmp1, xtmp2, ztmp1, ztmp2, xtmp3, ztmp3, NULL);
	mpz_set(Xu, X);
	mpz_set_ui(Zu, 1);
	mpz_set_ui(X0, 1);
	mpz_set_ui(Z0, 0);
	mpz_set(X1, Xu);
	mpz_set(Z1, Zu);
	
	mpz_realloc(xtmp1, size);
	mpz_realloc(ztmp1, size);
	mpz_realloc(xtmp2, size);
	mpz_realloc(ztmp2, size);
	mpz_realloc(xtmp3, size);
	mpz_realloc(ztmp3, size);
	mpz_realloc(X0, size);
	mpz_realloc(Z0, size);
	mpz_realloc(X1, size);
	mpz_realloc(Z1, size);
	mpz_realloc(Xu, size);
	mpz_realloc(Zu, size);
	
	for(i = mpz_size(m) * sizeof(int) - 1; i > 0; i--)
	{
		printf("0\n");
		xDBL(N, A, X0, Z0, xtmp1, ztmp1);
		printf("1\n");
		xADD(N, A, X0, Z0, X1, Z1, Xu, Zu, xtmp2, ztmp2);
		printf("2\n");
		limbs_tmp1 = mpz_limbs_modify(xtmp1, size);
		limbs_tmp2 = mpz_limbs_modify(xtmp2, size);
		n = mod2(m);
		mpn_cnd_swap(n, limbs_tmp1, limbs_tmp2, size);
		limbs_tmp1 = mpz_limbs_modify(ztmp1, size);
		limbs_tmp2 = mpz_limbs_modify(ztmp2, size);
		mpn_cnd_swap(n, limbs_tmp1, limbs_tmp2, size);
		
		xADD(N, A, X0, Z0, X1, Z1, Xu, Zu, xtmp2, ztmp2);
		xDBL(N, A, X1, Z1, xtmp3, ztmp3);
		limbs_tmp2 = mpz_limbs_modify(xtmp2, size);
		limbs_tmp3 = mpz_limbs_modify(xtmp3, size);
		mpn_cnd_swap(n, limbs_tmp2, limbs_tmp3, size);
		limbs_tmp2 = mpz_limbs_modify(ztmp2, size);
		limbs_tmp3 = mpz_limbs_modify(ztmp3, size);
		mpn_cnd_swap(n, limbs_tmp2, limbs_tmp3, size);
		printf("3\n");
		
		mpz_limbs_finish(xtmp1, size);
		mpz_limbs_finish(ztmp1, size);
		mpz_limbs_finish(xtmp2, size);
		mpz_limbs_finish(ztmp2, size);
		
		mpz_set(X0, xtmp1);
		mpz_set(Z0, ztmp1);
		mpz_set(X1, xtmp2);
		mpz_set(Z1, ztmp2);
	}
	mpz_set(Xm, X0);
	mpz_set(Zm, Z0);
	
	mpz_clears(X0, Z0, X1, Z1, Xu, Zu, xtmp1, xtmp2, ztmp1, ztmp2, xtmp3, ztmp3, NULL);
}

void calc_A(mpz_t a, mpz_t n)
{
	mp_limb_t *limbs_inv, *limbs_tmp, *scratch;
	const mp_limb_t *limbs_n, *limbs_four, *limbs_a;
	mpz_t four, inv, tmp;
	mpz_inits(four, inv, tmp, NULL);
	
	mpz_set_ui(four, 4);
	limbs_four = mpz_limbs_modify(four, size);
	limbs_inv = mpz_limbs_write(inv, size);
	limbs_n = mpz_limbs_read(n);
	limbs_a = mpz_limbs_read(a);
	limbs_tmp = mpz_limbs_write(tmp, size);
	scratch = malloc(mpn_sec_invert_itch(size) * sizeof(mp_limb_t));
	mpn_sec_invert(limbs_inv, limbs_four, limbs_n, size, 2 * size * GMP_NUMB_BITS, scratch);
	mpz_limbs_finish(inv, size);
	
	mpz_sub_ui(four, four, 2);
	limbs_four = mpz_limbs_modify(four, size);
	
	mpn_add_n(limbs_tmp, limbs_a, limbs_four, size);
	mpz_limbs_finish(tmp, size);
	sec_mul(tmp, inv, a, size);
	sec_div_r(a, n, size);
	mpz_limbs_finish(a, size);
	
	free(scratch);
	mpz_clears(four, inv, tmp, NULL);
}

int main()
{
	mpz_t n, a, x1, z1, x2, z2, xm, zm, xu, zu, m;
	mpz_inits(n, a, x1, z1, x2, z2, xm, zm, xu, zu, m, NULL);
	
	mpz_set_ui(n, 101);
	size = mpz_size(n);
	mpz_realloc(a, size);
	mpz_realloc(x1, size);
	mpz_realloc(z1, size);
	mpz_realloc(x2, size);
	mpz_realloc(z2, size);
	mpz_realloc(xu, size);
	mpz_realloc(zu, size);
	mpz_realloc(m, size);
	mpz_set_ui(a, 49);
	mpz_set_ui(x1, 2);
	mpz_set_ui(z1, 1);
	mpz_set_ui(m, 2);
	mpz_set_ui(x2, 0);
//	mpz_set_ui(z2, 4);
//	mpz_set_ui(xu, 22);
//	mpz_set_ui(zu, 1);
	
//	gmp_printf("%Zd + 2 / 4 ^ -1 mod %Zd = ", a, n);
	
	calc_A(a, n);
	
//	gmp_printf("x2 = %Zd\n", x2);
	
//	sec_mul(a, x1, z2);
	
//	gmp_printf("%Zd * %Zd = %Zd\n", a, x1, z2);
	
	ladder(n, a, m, x1, z1, xm, zm);
	
//	gmp_printf("X = %Zd\nZ = %Zd\n", xm, zm);
	
	mpz_clears(n, a, x1, z1, x2, z2, xm, zm, xu, zu, m, NULL);
	
	return 0;
}
