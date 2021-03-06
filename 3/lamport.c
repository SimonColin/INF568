#include "stdlib.h"
#include "stdio.h"
#include "time.h"
#include "string.h"
#include "keccak.h"

void lamportKeygen(u8 * r0, u8 * r1, u8 * hr0, u8 * hr1, int msgl, int hashl)
{
	int i, j;
	u8 *tmp;
	srand(ctime(NULL));
	tmp = malloc(hashl * sizeof(u8));
	for(i = 0; i < hashl * msgl * 8; i++)
	{
		r0[i] = rand() % 256;
		r1[i] = rand() % 256;
	}
	for(i = 0; i < msgl * 8; i ++)
	{
		memcpy(tmp, r0 + i * hashl, hashl);
		FIPS202_SHAKE128(tmp, hashl, hr0 + i * hashl, hashl);
		memcpy(tmp, r1 + i * hashl, hashl);
		FIPS202_SHAKE128(tmp, hashl, hr1 + i * hashl, hashl);
	}
}

void lamportSign(u8 * msg, u8 * sign, u8 * r0, u8 * r1, int msgl, int hashl)
{
	int i, j;
	u8 * use;
	for(i = 0; i < msgl * 8; i++)
	{
		if((msg[i] >> (i % 8)) % 2 == 1)
			use = r1;
		else
			use = r0;
		memcpy(sign + i * hashl, use + i * hashl, hashl);
	}
}

int lamportVerify(u8 * msg, u8 * sign, u8 * hr0, u8 * hr1, int msgl, int hashl)
{
	int i, j;
	u8 * use, * hash;
	hash = malloc(hashl * sizeof(u8));
	for(i = 0; i < msgl * 8; i ++)
	{
		if((msg[i] >> (i % 8)) % 2 == 1)
			use = hr1;
		else
			use = hr0;
		FIPS202_SHAKE128(sign + i * hashl, hashl, hash, hashl);
		if(memcmp(hash, use + i * hashl, hashl) != 0)
			return 0;
	}
	return 1;
}

int main()
{
	u8 * msg, * s0, * s1, * v0, * v1, *sign;
	int length, i, j;
	length = 32;
	s0 = malloc(length * length * 8 * sizeof(u8));
	s1 = malloc(length * length * 8 * sizeof(u8));
	v0 = malloc(length * length * 8 * sizeof(u8));
	v1 = malloc(length * length * 8 * sizeof(u8));
	msg = malloc(length * sizeof(u8));
	sign = malloc(length * length * 8 * sizeof(u8));
	memcpy(msg, "abcdabcdabcdabcdabcdabcdabcdabcd", length);
	lamportKeygen(s0, s1, v0, v1, length, length);
	lamportSign(msg, sign, s0, s1, length, length);
//	sign[2] ^= 255;
	if(lamportVerify(msg, sign, v0, v1, length, length))
		printf("success\n");
	else
		printf("failure\n");
}
	
