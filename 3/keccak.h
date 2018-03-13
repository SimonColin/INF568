#define FOR(i,n) for(i=0; i<n; ++i)
typedef unsigned char u8;
typedef unsigned long long int u64;
typedef unsigned int ui;

void FIPS202_SHAKE128(const u8 *in, u64 inLen, u8 *out, u64 outLen);
