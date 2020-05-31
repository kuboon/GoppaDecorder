/* generate GF(2^n) using irreducible polynomial */
//ゼフ対数表を作るためのプログラム。正規基底を生成します。
#include <stdio.h>
#include <stdlib.h>
#ifndef M
  #define M 8192
#endif

unsigned int gf[M], fg[M];

unsigned int intlog2(unsigned int x)
{
  unsigned int ans = 0;
  while (x >>= 1)
    ans++;
  return ans;
}

void make_gf_fg(int n)
{
  /* generate Galois Field over GF(2^?) */
  const char *normal[14] = {
      "1011",
      //"11001", /* GF(16) */
      "10011",
      "110111",
      "1100001",
      "11000001",
      "110101001",
      //"100011011",
      "1100110001",
      "11000010011",
      "110000001101",
      "1100101000001",  /* 4096 */
      "11011000000001", /* 8192 */
      "110000100010001",
      "1100000000000001",
      "11010000000010001"};
  unsigned int MM = strtoul(normal[intlog2(n) - 3], (char **)NULL, 2) ^ n;

  gf[0] = 0;
  int bit = 1;
  for (int i = 1; i < n; i++)
  {
    gf[i] = bit;
    bit = (bit << 1);
    if (bit > n - 1)
      bit = (bit - n) ^ MM;
  }

  for (int i = 0; i < n; i++)
    fg[gf[i]] = i;
}
