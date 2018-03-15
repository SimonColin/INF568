#include "stdlib.h"
#include "stdio.h"
#include "time.h"
#include "string.h"
#include "keccak.h"

void chop(u8 * msgin, u8 * msgout, int msglength)
{
	int i;
	for(i = 0; i < msglength; i++)
	{
		msgout[2 * i] = msgin[i] / 16;
		msgout[2 * i + 1] = msgin[i] % 16;
		printf("%d : %d\n%d : %d\n", 2 * i, msgout[2 * i], 2 * i + 1, msgout[2 * i + 1]);
	}
}

void genRandVals(u8 * values, int lengthvals)
{
	int i;
	srand(ctime(NULL));
	for(i = 0; i < lengthvals; i++)
	{
		values[i] = rand() % 256;
	}
}

void computePublic(u8 * valuesin, u8 * valuesout, int valuesize, int valueamount)
{
	int i, j;
	u8 * hash, * tmp;
	hash = malloc(valuesize * sizeof(u8));
	tmp = malloc(valuesize * sizeof(u8));
	for( i = 0; i < valueamount; i++)
	{
		memcpy(tmp, valuesin + i * valuesize, valuesize);
		for(j = 0; j < 16; j++)
		{
			FIPS202_SHAKE128(tmp, valuesize, hash, valuesize);
			memcpy(tmp, hash, valuesize);
		}
		memcpy(valuesout + i * valuesize, hash, valuesize);
	}
}

void computeSignature(u8 * valuesin, u8 * valuesout, u8 * msgin, int valuesize, int valueamount)
{
	int i, j;
	u8 * hash, * tmp;
	hash = malloc(valuesize * sizeof(u8));
	tmp = malloc(valuesize * sizeof(u8));
	for( i = 0; i < valueamount; i++)
	{
		memcpy(tmp, valuesin + i * valuesize, valuesize);
		for(j = 0; j < msgin[i]; j++)
		{
			FIPS202_SHAKE128(tmp, valuesize, hash, valuesize);
			memcpy(tmp, hash, valuesize);
		}
		memcpy(valuesout + i * valuesize, hash, valuesize);
	}
}

int check(u8 * msgchopped, u8 * public, u8 * signature, int hashlength, int msglength)
{
	int i, j;
	u8 * hash, * tmp;
	hash = malloc(hashlength * sizeof(u8));
	tmp = malloc(hashlength * sizeof(u8));
	for( i = 0; i < msglength; i++)
	{
		memcpy(tmp, signature + i * hashlength, hashlength);
		for(j = 0; j < 16 - msgchopped[i]; j++)
		{
			FIPS202_SHAKE128(tmp, hashlength, hash, hashlength);
			memcpy(tmp, hash, hashlength);
		}
		if(memcmp(public + i * hashlength, hash, hashlength) != 0)
		{
			return 0;
		}
	}
	return 1;
}

int main()
{
	int i;
	u8 * msg, * msgchopped, * randvals, * public, * sign;
	int msglength, hashlength, size;
	msglength = 4;//32;
	hashlength = 4;//32;
	size = msglength * 2 * hashlength;
	msg = malloc(msglength * sizeof(u8));
	memcpy(msg, "abcdabcdabcdabcdabcdabcdabcdabcd", msglength);
	msgchopped = malloc(msglength * 2 * sizeof(u8));
	randvals = malloc(size * sizeof(u8));
	public = malloc(size * sizeof(u8));
	sign = malloc(size * sizeof(u8));
	chop(msg, msgchopped, msglength);
	genRandVals(randvals, size);
	computePublic(randvals, public, hashlength, msglength);
	computeSignature(randvals, sign, msgchopped, hashlength, msglength);
	printf("before : %x\n", sign[3]);
	sign[3] ^= 255;
	printf("after : %x\n", sign[3]);
		for(i = 0; i < msglength * 2; i++)
	{
		printf("%d : %d\n", i, msgchopped[i]);
	}
	if(check(msgchopped, public, sign, hashlength, msglength) == 1)
	{
		printf("success\n");
	}
	else
	{
		printf("failure\n");
	}
}

