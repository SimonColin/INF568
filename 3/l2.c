#include "stdlib.h"
#include "stdio.h"
#include "time.h"
#include "string.h"
#include "keccak.h"

void lamportKeygen(u8 ** s0, u8 ** s1, u8 ** v0, u8 ** v1, int inputLength, int outputLength)
{
	int i, j;
	srand(ctime(NULL));
	s0 = malloc(inputLength * sizeof(u8*));
	s1 = malloc(inputLength * sizeof(u8*));
	v0 = malloc(inputLength * sizeof(u8*));
	v1 = malloc(inputLength * sizeof(u8*));
	printf("gen\n");
	for(i = 0; i < inputLength; i++)
	{
		s0[i] = malloc(sizeof(u8) * outputLength / 8);
		s1[i] = malloc(sizeof(u8) * outputLength / 8);
		v0[i] = malloc(sizeof(u8) * outputLength / 8);
		v1[i] = malloc(sizeof(u8) * outputLength / 8);
//		printf("%d = ", i);
		for(j = 0; j < outputLength / 8; j++)
		{
			s0[i][j] = rand() % 256;
			s1[i][j] = rand() % 256;
//			printf("%d : ", s0[i][j]);
		}
//		printf("%d\n", i);
		printf("%d : ", s1[i][0]);
		printf("\n");
		FIPS202_SHAKE128(s0[i], outputLength / 8, v0[i], outputLength / 8);
		FIPS202_SHAKE128(s1[i], outputLength / 8, v1[i], outputLength / 8);
	}
}

void lamportSign(u8 * msg, u8 ** sign, u8 ** s0, u8 ** s1, int inputLength, int outputLength)
{
	int i, j;
	u8 ** use;
	sign = malloc(outputLength * sizeof(u8*) / 8);
	printf("sign\n");
	printf("length = %d\n", outputLength / 8);
	for(i = 0; i < inputLength / 8; i++)
	{
		for(j = 0; j < 8; j++)
		{
			sign[j + 8 * i] = malloc(outputLength * sizeof(u8) / 8);
//			memcpy(sign[j + (8 * i)], s1[j + (8 * i)], outputLength / 8);
			printf("%d - %d\n", j + 8 * i, s1[j + (8 * i)][0]);
		/*	if((msg[i] >> j) % 2)
				
			else
				memcpy(sign[j + 8 * i], s0[j + 8 * i], outputLength / 8);*/
		}
	}
}

int lamportVerify(u8 * msg, u8 ** sign, u8 ** v0, u8 ** v1, int inputLength, int outputLength)
{
	int i, j;
	u8 ** use, * hash;
//	hash = malloc(outputLength * sizeof(u8*));
	for(i = 0; i < inputLength; i++)
	{
		for(j = 0; j < 8; j++)
		{
/*			if((msg[i] >> j) % 2)
				use = v1;
			else
				use = v0;
			FIPS202_SHAKE128(sign[j + 8 * i], outputLength, hash, outputLength);
			if(!memcmp(hash, use[j + 8 * i], inputLength))
				return 0;
*/		}
	}
	return 1;
}

int main()
{
	u8 * msg, ** s0, ** s1, ** v0, ** v1, **sign;
	int length;
	length = 32;
	msg = malloc(length * sizeof(u8) / 8);
	memcpy(msg, "abcdabcdabcdabcdabcdabcdabcdabcd", 32);
	lamportKeygen(s0, s1, v0, v1, length, length);
	lamportSign(msg, sign, s0, s1, length, length);
/*	if(lamportVerify(msg, sign, v0, v1, length, length))
		printf("success\n");
	else
		printf("failure\n");*/
}
	
