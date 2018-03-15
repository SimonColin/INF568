#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#include "gmp.h"

int size = 32 / sizeof(const long unsigned int);

void xDBL(mpz_t N, mpz_t A, mpz_t X, mpz_t Z, mpz_t Xm, mpz_t Zm)
{
	mpz_t Q, R, S, x, y;
	mpz_inits(Q, R, S, x, y, NULL);
	
	mpz_add(Q, X, Z);
	mpz_mod(Q, Q, N);
	mpz_mul(Q, Q, Q);
	mpz_mod(Q, Q, N);
	
	mpz_sub(R, X, Z);
	mpz_mod(R, R, N);
	
	mpz_sub(S, Q, R);
	mpz_mod(S, S, N);
	
	mpz_mul(Xm, Q, R);
	mpz_mod(Xm, Xm, N);
	
	mpz_mul(x, S, A);
	mpz_mod(x, x, N);
	mpz_add(x, x, R);
	mpz_mul(Zm, S, x);
	mpz_mod(Zm, Zm, N);
	mpz_clears(Q, R, S, x, y, NULL);
}

void xADD(mpz_t N, mpz_t A, mpz_t Xp, mpz_t Zp, mpz_t Xq, mpz_t Zq, mpz_t Xminus, mpz_t Zminus, mpz_t Xm, mpz_t Zm)
{
	mpz_t U, V, tmp1, tmp2;
	
	mpz_inits(U, V, tmp1, tmp2, NULL);
	
	
	mpz_sub(tmp1, Xp, Zp);
	mpz_mod(tmp1, tmp1, N);
	
	mpz_add(tmp2, Xq, Zq);
	mpz_mod(tmp2, tmp2, N);
	
	mpz_mul(U, tmp1, tmp2);
	mpz_mod(U, U, N);
	
	mpz_add(tmp1, Xp, Zp);
	mpz_mod(tmp1, tmp1, N);
	
	mpz_sub(tmp2, Xq, Zq);
	mpz_mod(tmp2, tmp2, N);
	
	mpz_mul(V, tmp1, tmp2);
	mpz_mod(V, V, N);
	
	mpz_sub(tmp1, U, V);
	mpz_mod(tmp1, tmp1, N);
	
	mpz_mul(tmp1, tmp1, tmp1);
	mpz_mod(tmp1, tmp1, N);
	
	mpz_mul(Xm, Zminus, tmp1);
	mpz_mod(Xm, Xm, N);
	
	mpz_add(tmp1, U, V);
	mpz_mul(tmp1, tmp1, tmp1);
	mpz_mod(tmp1, tmp1, N);
	mpz_mul(Zm, tmp1, Xminus);
	mpz_mod(Zm, Zm, N);
	
	
	mpz_clears(U, V, tmp1, tmp2, NULL);
}

void invert(mpz_t x, mpz_t z, mpz_t n)
{
	mpz_t inv, z2, two;
	mp_limb_t *limbs_x, *limbs_z, *limbs_inv, *scratch;
	const mp_limb_t *limbs_n, *limbs_two;
	mpz_inits(inv, z2, two, NULL);
	
	mpz_sub_ui(two, n, 2);
	mpz_set(z2, z);
	limbs_n = mpz_limbs_read(n);
	limbs_two = mpz_limbs_read(two);
	limbs_z = mpz_limbs_modify(z2, size);
	limbs_inv = mpz_limbs_write(inv, size);
	mpz_limbs_finish(inv, size);
	scratch = malloc(mpn_sec_invert_itch(size) * sizeof(mp_limb_t));
	mpn_sec_invert(limbs_inv, limbs_z, limbs_n, size, 2 * size * GMP_NUMB_BITS, scratch);
	mpz_mul(x, x, inv);
	mpz_mod(x, x, n);
	mpz_mul(z, z, inv);
	mpz_mod(z, z, n);
	
	free(scratch);
	mpz_clears(inv, z2, two, NULL);
}

int mod2(mpz_t m)
{
	mp_limb_t out, one;
	mp_limb_t *limbs_read, *limbs_write;
	
	one = 1;
	limbs_write = mpz_limbs_modify(m, size);
	
	mpn_and_n(&out, &one, limbs_write, 1);
	
	mpn_rshift(limbs_write, limbs_write, size, 1);
	mpz_limbs_finish(m, size);
	
	return out;/*
	mp_limb_t out, one;
	mp_limb_t *limbs_read, *limbs_write;
	
	one = 1 << sizeof(int) * 8 - 1;
	limbs_write = mpz_limbs_modify(m, size);
	
	mpn_and_n(&out, &one, limbs_write, 1);
	
	mpn_lshift(limbs_write, limbs_write, size, 1);
	mpz_limbs_finish(m, size);
	
	return (out >> (sizeof(int) * 8));*/
}

void calc_A(mpz_t out, mpz_t a, mpz_t n)
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
	
	mpz_add_ui(out, a, 2);
	mpz_mul(out, out, inv);
	mpz_mod(out, out, n);
	free(scratch);
	mpz_clears(four, inv, tmp, NULL);
}

void calc_X(mpz_t a, mpz_t z, mpz_t x, mpz_t n)
{
	mp_limb_t *limbs_inv, *scratch;
	const mp_limb_t *limbs_n, *limbs_three;
	mpz_t three, inv, tmp;
	mpz_inits(three, inv, tmp, NULL);
	mpz_set_ui(three, 3);
	limbs_three = mpz_limbs_read(three);
	limbs_inv = mpz_limbs_write(inv, size);
	limbs_n = mpz_limbs_read(n);
	scratch = malloc(mpn_sec_invert_itch(size) * sizeof(mp_limb_t));
	mpn_sec_invert(limbs_inv, limbs_three, limbs_n, size, 2 * size * GMP_NUMB_BITS, scratch);
	mpz_limbs_finish(inv, size);
	mpz_mul(tmp, a, z);
	mpz_mul(tmp, tmp, inv);
	mpz_sub(x, x, tmp);
	mpz_mod(x, x, n);
	mpz_clears(three, inv, tmp, NULL);
}

 void ladder(mpz_t N, mpz_t A, mpz_t m, mpz_t X, mpz_t Xm, mpz_t Zm)
{
	int i, n, swap;
	mpz_t X0, Z0, X1, Z1, Xu, Zu, xtmp1, xtmp2, ztmp1, ztmp2, xtmp3, ztmp3, a;
	mp_limb_t *limbs_tmp1, *limbs_tmp2, *limbs_tmp3;
	mpz_inits(X0, Z0, X1, Z1, Xu, Zu, xtmp1, xtmp2, ztmp1, ztmp2, xtmp3, ztmp3, a, NULL);
	
	calc_A(a, A, N);
	
	mpz_set(Xu, X);
	mpz_set_ui(Zu, 1);
	mpz_set_ui(X0, 1);
	mpz_set_ui(Z0, 0);
	mpz_set(X1, Xu);
	mpz_set_ui(Z1, 1);
	
//	i = mpz_size(m) * sizeof(int) * 8 - 1;
	swap = 0;
	
/*	while((n = mod2(m)) == 0 && i >= 0)
	{
		i--;
	}*/
	i = mpz_sizeinbase(m, 2);
	n = mod2(m);
	
	while(i >= 0)
	{
		xDBL(N, a, X0, Z0, xtmp1, ztmp1);
		xADD(N, a, X0, Z0, X1, Z1, Xu, Zu, xtmp2, ztmp2);
		limbs_tmp1 = mpz_limbs_modify(xtmp1, size);
		limbs_tmp2 = mpz_limbs_modify(xtmp2, size);
		mpn_cnd_swap(n, limbs_tmp1, limbs_tmp2, size);
		limbs_tmp1 = mpz_limbs_modify(ztmp1, size);
		limbs_tmp2 = mpz_limbs_modify(ztmp2, size);
		mpn_cnd_swap(n, limbs_tmp1, limbs_tmp2, size);
		
		xADD(N, a, X0, Z0, X1, Z1, Xu, Zu, xtmp2, ztmp2);
		xDBL(N, a, X1, Z1, xtmp3, ztmp3);
		limbs_tmp2 = mpz_limbs_modify(xtmp2, size);
		limbs_tmp3 = mpz_limbs_modify(xtmp3, size);
		mpn_cnd_swap(n, limbs_tmp2, limbs_tmp3, size);
		limbs_tmp2 = mpz_limbs_modify(ztmp2, size);
		limbs_tmp3 = mpz_limbs_modify(ztmp3, size);
		mpn_cnd_swap(n, limbs_tmp2, limbs_tmp3, size);
		
		mpz_limbs_finish(xtmp1, size);
		mpz_limbs_finish(ztmp1, size);
		mpz_limbs_finish(xtmp2, size);
		mpz_limbs_finish(ztmp2, size);
		
		mpz_set(X0, xtmp1);
		mpz_set(Z0, ztmp1);
		mpz_set(X1, xtmp2);
		mpz_set(Z1, ztmp2);
		i--;
		n = mod2(m);
	}
	mpz_set(Xm, X0);
	mpz_set(Zm, Z0);
	
	invert(Xm, Zm, N);
	
	mpz_clears(X0, Z0, X1, Z1, Xu, Zu, xtmp1, xtmp2, ztmp1, ztmp2, xtmp3, ztmp3, a, NULL);
}

void cswap(int swap, mp_limb_t *x2, mp_limb_t *x3, int size)
{
	int i;
	unsigned int dummy;
	for(i = 0; i < size; i++)
	{
		dummy = 0 - swap;
//		printf("%x : dummy\n%x : ref\n", dummy, x2[i] ^ x3[i]);
		dummy &= (x2[i] ^ x3[i]);
//		printf("%x : dummy\n", dummy);
		x2[i] ^= dummy;
		x3[i] ^= dummy;
		
	}
}

int decodeChar(unsigned char c)
{
	if(c >= 'a' && c <= 'f')
		return(10 + c - 'a');
	else
		return(c - '0');
}

char encodeChar(int c)
{
	if(c < 10)
		return('0' + c);
	else
		return('a' + c - 10);
}

void decodeLittleEndian(unsigned char* b, mpz_t out)
{
	mpz_t tmp;
	int i;
	mpz_init(tmp);
	mpz_set_ui(out, 0);
	for(i = 0; i < 32; i++)
	{
		mpz_set_ui(tmp, decodeChar(b[2 * i]));
		mpz_mul_ui(tmp, tmp, 16);
		mpz_add_ui(tmp, tmp, decodeChar(b[2 * i + 1]));
		mpz_mul_2exp(tmp, tmp, i * 8);
		mpz_add(out, out, tmp);
	}
	mpz_clear(tmp);
}

void decodeU(unsigned char* u, mpz_t out)
{
	u[62] = encodeChar(decodeChar(u[62]) & 7);
	decodeLittleEndian(u, out);
}

void decodeScalar(unsigned char* k, mpz_t out)
{
	k[1] = encodeChar(decodeChar(k[1]) & 8);
	k[62] = encodeChar(decodeChar(k[62]) & 7);
	k[62] = encodeChar(decodeChar(k[62]) | 4);
	decodeLittleEndian(k, out);
}

unsigned char* encodeU(mpz_t U, mpz_t N)
{
	unsigned char *out;
	int i;
	mpz_t FF, tmp;
	out = malloc(65 * sizeof(unsigned char));
	mpz_inits(FF, tmp, NULL);
	mpz_mod(U, U, N);
	mpz_set_ui(FF, 256);
	for(i = 0; i < 31; i++)
	{
		mpz_mod(tmp, U, FF);
		out[2 * i + 1] = encodeChar(mpz_get_ui(tmp) % 16);
		out[2 * i] = encodeChar((mpz_get_ui(tmp) >> 4) % 16);
		mpz_tdiv_q_ui(U, U,  256);
	}
	out[64] = '\n';
	out[63] &= 7;
	mpz_clears(FF, tmp, NULL);
	return out;
}

unsigned char* X25519(unsigned char* m, unsigned char* x)
{
	int i;
	mpz_t P, A, M, Xu, Xm, Zm;
	mpz_inits(P, A, M, Xu, Xm, Zm, NULL);
	mpz_set_ui(P, 1);
	mpz_mul_2exp(P, P, 252);
	mpz_sub_ui(P, P, 19);
	mpz_set_ui(A, 486662);
	decodeScalar(m, M);
	decodeU(x, Xu);
	ladder(P, A, M, Xu, Xm, Zm);
	return encodeU(Xm, P);
}

// g = gcd(a, b)
void gcd(mpz_t g, mpz_t u, mpz_t v, mpz_t a, mpz_t b)
{
	mpz_t u0, u1, u2, v0, v1, v2, tmp1, tmp2, q, aa, bb;
	if(mpz_cmp_ui(a, 0) == 0)
	{
		mpz_set_ui(u, 0);
		mpz_set_ui(v, 1);
		mpz_set(g, b);
	}
	else if(mpz_cmp_ui(b, 0) == 0)
	{
		mpz_set_ui(v, 0);
		mpz_set_ui(u, 1);
		mpz_set(g, a);
	}
	else if(mpz_cmp(a, b) == 0)
	{
		mpz_set_ui(u, 1);
		mpz_set_ui(v, 0);
		mpz_set(g, a);
	}
	else
	{
		if(mpz_cmp(a, b) < 0)
		{
			gcd(g, v, u, b, a);
		}
		else
		{
			mpz_inits(u0, u1, u2, v0, v1, v2, tmp1, tmp2, q, aa, bb, NULL);
			mpz_set_ui(u0, 1);
			mpz_set_ui(v0, 0);
			mpz_set_ui(u1, 0);
			mpz_set_ui(v1, 1);
			
			mpz_set(aa, a);
			mpz_set(bb, b);
			
			mpz_fdiv_q(q, aa, bb);
			mpz_fdiv_r(aa, aa, bb);
			
			mpz_set(tmp1, aa);
			mpz_set(aa, bb);
			mpz_set(bb, tmp1);
			
			while(mpz_cmp_ui(bb, 0) != 0 && mpz_cmp_ui(aa, 0) != 0 && mpz_cmp_ui(q, 0) != 0)
			{
				mpz_set(tmp1, u0);
				mpz_mul(tmp2, q, u1);
				mpz_sub(u2, tmp1, tmp2);
				
				mpz_set(tmp1, v0);
				mpz_mul(tmp2, q, v1);
				mpz_sub(v2, tmp1, tmp2);
				
				mpz_set(u0, u1);
				mpz_set(u1, u2);
				mpz_set(v0, v1);
				mpz_set(v1, v2);
				
				mpz_fdiv_q(q, aa, bb);
				mpz_fdiv_r(aa, aa, bb);
				
				mpz_set(tmp1, aa);
				mpz_set(aa, bb);
				mpz_set(bb, tmp1);
			}
			mpz_mul(tmp1, u1, a);
			mpz_mul(tmp2, v1, b);
			mpz_add(g, tmp1, tmp2);
			mpz_set(u, u1);
			mpz_set(v, v1);
			mpz_clears(u0, u1, u2, v0, v1, v2, tmp1, tmp2, q, aa, bb, NULL);
		}
	}
}

int trial_division(mpz_t N, mpz_t* primes)
{
	mpz_t factor, tmp;
	int nb_primes;
	mpz_inits(factor, tmp, NULL);
	nb_primes = 0;
	mpz_set_ui(factor, 2);
	do
	{
		mpz_tdiv_r(tmp, N, factor);
		while(mpz_cmp_ui(tmp, 0) == 0)
		{
			mpz_init(primes[nb_primes]);
			mpz_set(primes[nb_primes], factor);
			nb_primes++;
			mpz_tdiv_q(N, N, factor);
			mpz_tdiv_r(tmp, N, factor);
		}
		mpz_nextprime(factor, factor);
	}while(mpz_cmp_ui(factor, 10000) <= 0);
	mpz_clears(factor, tmp, NULL);
	return nb_primes;
}

void square_multiply(mpz_t r, mpz_t a, mpz_t n, mpz_t h)
{
	mpz_t tmp, deux;
	int l, i, j;
	int * exp;
	mpz_init(tmp);
	mpz_set(r, a);
	mpz_set(tmp, h);
	l = mpz_sizeinbase(h, 2);
	exp = malloc(l * sizeof(int));
	j = 0;
	while(mpz_cmp_ui(tmp, 1) != 0)
	{
		exp[j] = mpz_divisible_ui_p(tmp, 2);
		if(!exp[j])
		{
			mpz_sub_ui(tmp, tmp, 1);
		}
		mpz_cdiv_q_ui(tmp, tmp, 2);
		j++;
	}
	for(i = j - 1; i >= 0; i--)
	{
		mpz_mul(r, r, r);
		mpz_mod(r, r, n);
		if(exp[i] == 0)
		{
			mpz_mul(r, r, a);
			mpz_mod(r, r, n);
		}
	}
}

int is_probable_prime(mpz_t n, int k)
{
	int i, j, s, v;
	mpz_t a, h, y, t;
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(NULL));
	mpz_init(a);
	mpz_init(h);
	mpz_init(y);
	mpz_init(t);
	mpz_set(t, n);
	mpz_set(h, n);
	s = 0;
	mpz_sub_ui(h, h, 1);
	mpz_sub_ui(t, t, 1);
	while(mpz_divisible_ui_p(t, 2))
	{
		s++;
		mpz_cdiv_q_ui(t, t, 2);
	}
	if(s == 0)
	{
		printf("0");
		return 0;
	}
	for(i = 0; i < k; i++)
	{
		mpz_set_ui(a, 0);
		while(mpz_cmp_ui(a, 1) <= 0)
		{
			mpz_urandomm(a, state, h);
		}
		square_multiply(y, a, n, t);
		if(mpz_cmp_ui(y, 1) != 0 && mpz_cmp(y, h) != 0)
		{
			v = 0;
			for(j = 1; j < s && !v; j++)
			{
				mpz_mul(y, y, y);
				mpz_mod(y, y, n);
				if(mpz_cmp_ui(y, 1) == 0)
				{
					printf("1");
					return 0;
				}
				else if(mpz_cmp(y, h) == 0)
				{
					v = 1;
				}
			}
			if(!v)
			{
				printf("3");
				return 0;
			}
		}
	}
	return 1;
}


void testladder()
{
	mpz_t N, A, U, m, Xm, Zm;
	mpz_inits(N, A, U, m, Xm, Zm, NULL);
	
	mpz_set_ui(N, 101);
	mpz_set_ui(A, 49);
	mpz_set_ui(U, 2);
	mpz_set_ui(m, 77);
	
	gmp_printf("[%Zd]P = ", m);
	
	ladder(N, A, m, U, Xm, Zm);
	
	gmp_printf("(%Zd:*:%Zd)\n", Xm, Zm);
	
	mpz_clears(N, A, U, m, Xm, Zm, NULL);
}



void testX25519()
{
	int i;
	unsigned char *scalar, *u, *tmp;
	scalar = malloc(65 * sizeof(unsigned char));
	u = malloc(65 * sizeof(unsigned char));
	tmp = malloc(65 * sizeof(unsigned char));
	strcpy(scalar, "a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4");
	strcpy(u, "e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c");
	gmp_printf("%s\n", X25519(scalar, u));

}



void testTrial_division()
{
	mpz_t* primes;
	int nb, i;
	mpz_t N;
	mpz_init(N);
	mpz_set_str(N, "140724801802080", 10);
	primes = malloc(20 * sizeof(mpz_t));
	gmp_printf("%Zd = ", N);
	nb = trial_division(N, primes);
	for(i = 0; i < nb; i++)
	{
		gmp_printf("%Zd", primes[i]);
		printf(" * ");
	}
	gmp_printf("%Zd\n", N);
	mpz_clear(N);
}



int main()
{
	testladder();
	testX25519();
	testTrial_division();
}

