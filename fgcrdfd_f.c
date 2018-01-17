#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <imageio.h>
#include <nd_malloc.h>
#include <dft.h>
#include "fgcrdfd.h"
#include "bpf.h"
#include "pre_filter.h"

/*---------------------------------------------------------------------------*/
void mzbpf(double **x, int xs, int ys, int fr, double **y)
{
  int i, j, ii, jj;
  double tmp;
  
  for(j = -fr ; j < ys + fr ; j++) {
    for(i = -fr ; i < xs + fr ; i++) {
      tmp = 0.0;
      for(jj = -BPF_HS ; jj <= BPF_HS ; jj++) {
	for(ii = -BPF_HS ; ii <= BPF_HS ; ii++) {
	  tmp += x[i-ii][j-jj] * bpf[BPF_N][ii + BPF_HS][jj + BPF_HS];
	}
      }
      y[i][j] = tmp;
    }
  }
}

/*---------------------------------------------------------------------------*/
void pref(double **x, int xs, int ys, int fr, int prfn, double **y)
{
  int i, j, ii, jj;
  double tmp;
  
  for(j = -fr ; j < ys + fr ; j++) {
    for(i = -fr ; i < xs + fr ; i++) {
      tmp = 0.0;
      for(jj = -PRF_HS ; jj <= PRF_HS ; jj++) {
	for(ii = -PRF_HS ; ii <= PRF_HS ; ii++) {
	  tmp += x[i-ii][j-jj] * pre[prfn][ii + PRF_HS][jj + PRF_HS];
	}
      }
      y[i][j] = tmp;
    }
  }
}

/*---------------------------------------------------------------------------*/
void mzdfd(double ***bpfinp, int nz, int xbc, int ybc, int bs, 
	   double *dini, double *kini, int *maxz)
{
  int hbs, z, i, j;
  double *m, *logm, sum2, max;
  double Sx4, Sx3, Sx2, Sx1, Sx0, Syx2, Syx1, Syx0;
  double aa[3][3], bb[3], ans[3], a, b, c; 
  int sol3x3(double a[][3], double *, double *);

  /***** 初期化 *****/
  hbs = bs / 2;

  /***** 領域確保 *****/
  m    = malloc_double_1d(nz, 0);  
  logm = malloc_double_1d(nz, 0);  
  
  /***** 合焦点評価値 *****/
  max = -1.0;
  *maxz = 0;
  for(z = 0 ; z < nz ; z++) {
    sum2 = 0;
    for(j = ybc - hbs ; j <= ybc + hbs ; j++) {
      for(i = xbc - hbs ; i <= xbc + hbs ; i++) {
	sum2 += bpfinp[z][i][j] * bpfinp[z][i][j];
      }
    }
    m[z] = sum2;
    logm[z] = log(m[z]);
    if (max < m[z]) {
      max = m[z];
      *maxz = z;
    }
  }
  if (*maxz == 0) (*maxz)++;
  if (*maxz == nz - 1) (*maxz)--;
 
  /***** 最小2乗近似 *****/
  Sx4 = Sx3 = Sx2 = Sx1 = Sx0 = Syx2 = Syx1 = Syx0 = 0.0;
  for(z = *maxz - 1 ; z <= *maxz + 1 ; z++) {
    Sx4 += (double)(z * z * z * z);
    Sx3 += (double)(z * z * z);
    Sx2 += (double)(z * z);
    Sx1 += (double) z ;
    Sx0 += 1.0;
    Syx2 += logm[z] * (double)(z * z);
    Syx1 += logm[z] * (double) z;
    Syx0 += logm[z];
  }
  aa[0][0] = Sx4;
  aa[0][1] = aa[1][0] = Sx3;
  aa[0][2] = aa[1][1] = aa[2][0] = Sx2;
  aa[1][2] = aa[2][1] = Sx1;
  aa[2][2] = Sx0;
  bb[0] = Syx2;
  bb[1] = Syx1;
  bb[2] = Syx0;

  if (sol3x3(aa, bb, ans) < 0) {
    fprintf(stderr, "\nError : apprx\n\n"); exit(1);
  } else {
    a = ans[0];
    b = ans[1];
    c = ans[2];
    *dini = -b / 2.0 / a;
    *kini = sqrt(-a) / (2.0 * M_PI * fc[BPF_N]);  /* sgm^2 * omg^2 = -a */
#if DEBUG == 1
    printf("    (%3d, %3d) --> dini = %f, kini = %f (maxz = %d)\n", xbc, ybc, *dini, *kini, *maxz);
#endif
  }

  free_double_1d(m, nz, 0);  
  free_double_1d(logm, nz, 0);  
}

/*---------------------------------------------------------------------------*/
int sol3x3(double a[][3], double b[], double x[])
{ 
  double d[3][3], bunbo, bunshi;
  int i, ii, jj;
  double det3x3(double d[][3]);
    
  bunbo = det3x3(a);
  if (fabs(bunbo) < EPS)  {
    fprintf(stderr, "   ***** Singular matrix *****\n");
    fprintf(stderr,      "%f  %f  %f\n", d[0][0], d[1][0], d[2][0]);
    fprintf(stderr,      "%f  %f  %f\n", d[0][1], d[1][1], d[2][1]);
    fprintf(stderr,      "%f  %f  %f\n", d[0][2], d[1][2], d[2][2]);
    return -1;
  }
  for(i = 0 ; i < 3 ; i++) {
    for(ii = 0 ; ii < 3 ; ii++) {
      for(jj = 0 ; jj < 3 ; jj++) {
        if (jj == i) {
          d[ii][jj] = b[ii];
        } else {
          d[ii][jj] = a[ii][jj];
        }
      }
    }
    bunshi = det3x3(d);
    x[i] = bunshi / bunbo;
  }
    
  return 0;
}
   
/*---------------------------------------------------------------------------*/
double det3x3(double a[][3])
{
  return
    a[0][2]*a[2][1]*a[1][0] +
      a[0][1]*a[1][2]*a[2][0] +
        a[0][0]*a[1][1]*a[2][2] -
          a[0][0]*a[2][1]*a[1][2] -
            a[0][1]*a[1][0]*a[2][2] -
              a[0][2]*a[1][1]*a[2][0];
}

/*---------------------------------------------------------------------------*/
void fgcrdfd(double ***inp, int maxz, int xbc, int ybc, int bs, 
	     double dini, double kini, int wgt_flg, double *dest, double *kest)
{
  static int first = 1;
  static COMPLEX **tf0C, **tf1C, **tf2C;
  static COMPLEX **tF0C, **tF1C, **tF2C;
  static double  **tA0, **tA1, **tA2;
  static double  **tP0, **tP1, **tP2;
  static double  **tdP0dd, **tdP1dd, **tdP2dd, **tdP0dk, **tdP1dk, **tdP2dk; 
  static double  **tA0P1, **tA1P0, **tA1P2, **tA2P1, **tA0P2, **tA2P0;
  static double  **tA0dP1dd, **tA1dP0dd, **tA1dP2dd, **tA2dP1dd, **tA0dP2dd, **tA2dP0dd;
  static double  **tA0dP1dk, **tA1dP0dk, **tA1dP2dk, **tA2dP1dk, **tA0dP2dk, **tA2dP0dk;
  static double  **tD01, **tDd01, **tDk01, **tD12, **tDd12, **tDk12, **tD20, **tDd20, **tDk20;
  static int N, hbs;
  int z0, z1, z2, converge, iter;
  int i, j, ii, jj;
  double a[2][2], b[2], ans[2];
  double d, k, E, Eprev;
  void trwin_2D(COMPLEX **tf0, int N);
  void PSF(double z, double d, double k, int N, double **P);
  void PSFddk(double z, double d, double k, int N, double **P, double **dPdd, double **dPdk);
  void CABS(COMPLEX **a, int N, double **b);
  void mul_double_double(double **a, double **b, int N, double **c);
  void dif_double_double(double **a, double **b, int N, double **c);
  double S_double_double(double **a, double **b, int N);
  int sol2x2(double a[][2], double b[], double x[]);


  /***** 準備 *****/
  if (first == 1) {
    N = bs;   /* DFT サイズ */
    hbs = bs / 2;
    tf0C     = malloc_COMPLEX_2d(bs, 0, bs, 0);
    tf1C     = malloc_COMPLEX_2d(bs, 0, bs, 0);
    tf2C     = malloc_COMPLEX_2d(bs, 0, bs, 0);
    tF0C     = malloc_COMPLEX_2d(bs, 0, bs, 0);
    tF1C     = malloc_COMPLEX_2d(bs, 0, bs, 0);
    tF2C     = malloc_COMPLEX_2d(bs, 0, bs, 0);
    tA0      = malloc_double_2d(bs, 0, bs, 0);
    tA1      = malloc_double_2d(bs, 0, bs, 0);
    tA2      = malloc_double_2d(bs, 0, bs, 0);
    tP0      = malloc_double_2d(bs, 0, bs, 0);
    tP1      = malloc_double_2d(bs, 0, bs, 0);
    tP2      = malloc_double_2d(bs, 0, bs, 0);
    tdP0dd   = malloc_double_2d(bs, 0, bs, 0);
    tdP1dd   = malloc_double_2d(bs, 0, bs, 0);
    tdP2dd   = malloc_double_2d(bs, 0, bs, 0);
    tdP0dk   = malloc_double_2d(bs, 0, bs, 0);
    tdP1dk   = malloc_double_2d(bs, 0, bs, 0);
    tdP2dk   = malloc_double_2d(bs, 0, bs, 0);
    tA0P1    = malloc_double_2d(bs, 0, bs, 0);
    tA1P0    = malloc_double_2d(bs, 0, bs, 0);
    tA1P2    = malloc_double_2d(bs, 0, bs, 0);
    tA2P1    = malloc_double_2d(bs, 0, bs, 0);
    tA0P2    = malloc_double_2d(bs, 0, bs, 0);
    tA2P0    = malloc_double_2d(bs, 0, bs, 0);
    tA0dP1dd = malloc_double_2d(bs, 0, bs, 0);
    tA1dP0dd = malloc_double_2d(bs, 0, bs, 0);
    tA1dP2dd = malloc_double_2d(bs, 0, bs, 0);
    tA2dP1dd = malloc_double_2d(bs, 0, bs, 0);
    tA0dP2dd = malloc_double_2d(bs, 0, bs, 0);
    tA2dP0dd = malloc_double_2d(bs, 0, bs, 0);
    tA0dP1dk = malloc_double_2d(bs, 0, bs, 0);
    tA1dP0dk = malloc_double_2d(bs, 0, bs, 0);
    tA1dP2dk = malloc_double_2d(bs, 0, bs, 0);
    tA2dP1dk = malloc_double_2d(bs, 0, bs, 0);
    tA0dP2dk = malloc_double_2d(bs, 0, bs, 0);
    tA2dP0dk = malloc_double_2d(bs, 0, bs, 0);
    tD01     = malloc_double_2d(bs, 0, bs, 0);
    tDd01    = malloc_double_2d(bs, 0, bs, 0);
    tDk01    = malloc_double_2d(bs, 0, bs, 0);
    tD12     = malloc_double_2d(bs, 0, bs, 0);
    tDd12    = malloc_double_2d(bs, 0, bs, 0);
    tDk12    = malloc_double_2d(bs, 0, bs, 0);
    tD20     = malloc_double_2d(bs, 0, bs, 0);
    tDd20    = malloc_double_2d(bs, 0, bs, 0);
    tDk20    = malloc_double_2d(bs, 0, bs, 0);
    
    first = 0;
  }

  /***** 初期化 *****/
  z0 = maxz - 1;
  z1 = maxz;
  z2 = maxz + 1;

  /***** ブロックにコピー *****/
  for(jj = 0, j = ybc - hbs ; j <= ybc + hbs ; j++, jj++) {
    for(ii = 0, i = xbc - hbs ; i <= xbc + hbs ; i++, ii++) {
      tf0C[ii][jj].r = inp[z0][i][j];
      tf1C[ii][jj].r = inp[z1][i][j];
      tf2C[ii][jj].r = inp[z2][i][j];
      tf0C[ii][jj].i = 0.0;
      tf1C[ii][jj].i = 0.0;
      tf2C[ii][jj].i = 0.0;
    }
  }
  
  /***** 2D DFT *****/
  trwin_2D(tf0C, N);
  trwin_2D(tf1C, N);
  trwin_2D(tf2C, N);
  DFT_2D(tf0C, tF0C, N);
  DFT_2D(tf1C, tF1C, N);
  DFT_2D(tf2C, tF2C, N);
  CABS(tF0C, N, tA0);
  CABS(tF1C, N, tA1);
  CABS(tF2C, N, tA2);
  
  /***** Gauss-Newton 法ループ *****/
  d = dini;
  k = kini;
  converge = 0;
  Eprev = 0.0;
  for(iter = 0 ; iter < GN_MAX ; iter++) {
    /*** PSF ***/
    PSFddk((double)z0, d, k, N, tP0, tdP0dd, tdP0dk);
    PSFddk((double)z1, d, k, N, tP1, tdP1dd, tdP1dk);
    PSFddk((double)z2, d, k, N, tP2, tdP2dd, tdP2dk);
    /*** tA0P1, tA1P0, ... ***/
    mul_double_double(tA0, tP1, N, tA0P1);
    mul_double_double(tA1, tP0, N, tA1P0);
    dif_double_double(tA0P1, tA1P0, N, tD01);
    mul_double_double(tA1, tP2, N, tA1P2);
    mul_double_double(tA2, tP1, N, tA2P1);
    dif_double_double(tA1P2, tA2P1, N, tD12);
    mul_double_double(tA0, tP2, N, tA0P2);
    mul_double_double(tA2, tP0, N, tA2P0);
    dif_double_double(tA0P2, tA2P0, N, tD20);
    /*** E 計算 ***/
    E = 0.0;
    E += S_double_double(tD01, tD01, N);
    E += S_double_double(tD12, tD12, N);
    E += S_double_double(tD20, tD20, N);
    /*** 終了判定 ***/    
    if (fabs(Eprev - E) / E < GN_EPS) {
      converge = 1;
      break;
    }
    /*** 次のループ準備 ***/    
    Eprev = E;
    /*** tF0dP1dd, tF1dP0dd ... ***/
    mul_double_double(tA0, tdP1dd, N, tA0dP1dd);
    mul_double_double(tA1, tdP0dd, N, tA1dP0dd);
    dif_double_double(tA0dP1dd, tA1dP0dd, N, tDd01);
    mul_double_double(tA1, tdP2dd, N, tA1dP2dd);
    mul_double_double(tA2, tdP1dd, N, tA2dP1dd);
    dif_double_double(tA1dP2dd, tA2dP1dd, N, tDd12);
    mul_double_double(tA0, tdP2dd, N, tA0dP2dd);
    mul_double_double(tA2, tdP0dd, N, tA2dP0dd);
    dif_double_double(tA0dP2dd, tA2dP0dd, N, tDd20);
    /*** tF0dP1dk, tF1dP0dk ... ***/
    mul_double_double(tA0, tdP1dk, N, tA0dP1dk);
    mul_double_double(tA1, tdP0dk, N, tA1dP0dk);
    dif_double_double(tA0dP1dk, tA1dP0dk, N, tDk01);
    mul_double_double(tA1, tdP2dk, N, tA1dP2dk);
    mul_double_double(tA2, tdP1dk, N, tA2dP1dk);
    dif_double_double(tA1dP2dk, tA2dP1dk, N, tDk12);
    mul_double_double(tA0, tdP2dk, N, tA0dP2dk);
    mul_double_double(tA2, tdP0dk, N, tA2dP0dk);
    dif_double_double(tA0dP2dk, tA2dP0dk, N, tDk20);
    /*** 連立方程式 ***/
    a[0][0]  = S_double_double(tDd01, tDd01, N);
    a[0][0] += S_double_double(tDd12, tDd12, N);
    a[0][0] += S_double_double(tDd20, tDd20, N);
    a[0][1]  = S_double_double(tDd01, tDk01, N);
    a[0][1] += S_double_double(tDd12, tDk12, N);
    a[0][1] += S_double_double(tDd20, tDk20, N);
    a[1][1]  = S_double_double(tDk01, tDk01, N);
    a[1][1] += S_double_double(tDk12, tDk12, N);
    a[1][1] += S_double_double(tDk20, tDk20, N);
    b[0]     = S_double_double(tD01,  tDd01, N);
    b[0]    += S_double_double(tD12,  tDd12, N);
    b[0]    += S_double_double(tD20,  tDd20, N);
    b[1]     = S_double_double(tD01,  tDk01, N);
    b[1]    += S_double_double(tD12,  tDk12, N);
    b[1]    += S_double_double(tD20,  tDk20, N);
    a[1][0]  = a[0][1];
    /*** 方程式解く ***/    
    if (sol2x2(a, b, ans) < 0) {
      converge = -1;
      break;
    } else {
    /*** 更新 ***/    
      d -= ans[0];
      k -= ans[1];
    }
    /*** 表示 ***/    
#if DEBUG == 1
    fprintf(stderr, "   #%03d d = %6.4f, k = %6.4f, E = %12.8f\n", iter, d, k, E);
    fprintf(stderr, "      |%12.8f  %12.8f| | d | = |%12.8f|\n", a[0][0], a[0][1], b[0]);
    fprintf(stderr, "      |%12.8f  %12.8f| | k | = |%12.8f|\n", a[1][0], a[1][1], b[1]);
#endif
  }
  /***** Gauss-Newton 法ループ終了処理 *****/
  if (converge == 1) {
    *dest = d;
    *kest = k;
  } else {    
    if (converge == 0) {
      fprintf(stderr, "(%3d, %3d) : no convergence --> dini = %5.3f returned\n", xbc, ybc, dini);
    } else if (converge == -1) {
      fprintf(stderr, "(%3d, %3d) : sol2x2 singular --> dini = %5.3f returned\n", xbc, ybc, dini);
    }
    *dest = dini;
    *kest = kini;
  }

#if DEBUG == 1
  printf("    (%3d, %3d) --> dest = %f, kest = %f (Emin = %f)\n", xbc, ybc, *dest, *kest, E);
#endif
    
}

/*---------------------------------------------------------------------------*/
void  trwin_2D(COMPLEX **f, int N)
{
  static int first = 1;
  int i, j, ii, jj, nc = (N - 1) / 2, absii, absjj;

  for(j = 0, jj = -nc ; j < N ; j++, jj++) {
    absjj = abs(jj);
    for(i = 0, ii = -nc ; i < N ; i++, ii++) {
      absii = abs(ii);
      if (absii == nc || absjj == nc) {
	f[i][j].r = 0.4 * f[i][j].r;
	f[i][j].i = 0.4 * f[i][j].i;
      } else if (absii == nc - 1 || absjj == nc - 1 ) {
	f[i][j].r = 0.5 * f[i][j].r;
	f[i][j].i = 0.5 * f[i][j].i;
      } else if (absii == nc - 2 || absjj == nc - 2 ) {
	f[i][j].r = 0.6 * f[i][j].r;
	f[i][j].i = 0.6 * f[i][j].i;
      } else if (absii == nc - 3 || absjj == nc - 3 ) {
	f[i][j].r = 0.7 * f[i][j].r;
	f[i][j].i = 0.7 * f[i][j].i;
      } else if (absii == nc - 4 || absjj == nc - 4 ) {
	f[i][j].r = 0.8 * f[i][j].r;
	f[i][j].i = 0.8 * f[i][j].i;
      } else if (absii == nc - 5 || absjj == nc - 5 ) {
	f[i][j].r = 0.9 * f[i][j].r;
	f[i][j].i = 0.9 * f[i][j].i;
      }
    }
  }

  if (first) {
    FILE *fp;

    fp = fopen("w.dat", "w");
    for(j = 0, jj = -nc ; j < N ; j++, jj++) {
      for(i = 0, ii = -nc ; i < N ; i++, ii++) {
	fprintf(fp, "%2d %2d,  %f\n", i, j, f[i][j].r);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
    first = 0;
  }    
}

/*---------------------------------------------------------------------------*/
void PSF(double z, double d, double k, int N, double **P)
{
  int i, j;
  double kzd, sgm2s2, pi2N, wx, wy, w2; 
  
  kzd = k * (z - d);
  sgm2s2 = kzd * kzd / 2.0;
  pi2N = 2.0 * M_PI / (double) N;

  for(j = 0 ; j < N ; j++) {
    wy = pi2N * (double) j;
    if (wy > M_PI) wy = 2.0 * M_PI - wy;
    for(i = 0 ; i < N ; i++) {
      wx = pi2N * (double) i;
      if (wx > M_PI) wx = 2.0 * M_PI - wx;
      w2 = wx * wx + wy * wy;
      P[i][j] = exp(-sgm2s2 * w2);
    }
  }
}  

/*---------------------------------------------------------------------------*/
void PSFddk(double z, double d, double k, int N, double **P, double **dPdd, double **dPdk)
{
  int i, j;
  double kzd, sgm2s2, pi2N, wx, wy, w2;
  
  kzd = k * (z - d);
  sgm2s2 = kzd * kzd / 2.0;
  pi2N = 2.0 * M_PI / (double) N;

  for(j = 0 ; j < N ; j++) {
    wy = pi2N * (double) j;
    if (wy > M_PI) wy = 2.0 * M_PI - wy;
    for(i = 0 ; i < N ; i++) {
      wx = pi2N * (double) i;
      if (wx > M_PI) wx = 2.0 * M_PI - wx;
      w2 = wx * wx + wy * wy;
      P[i][j] = exp(-sgm2s2 * w2);
      dPdd[i][j] =  P[i][j] * k * kzd * w2;
      dPdk[i][j] = -P[i][j] * kzd * (z - d) * w2;
    }
  }
}

/*---------------------------------------------------------------------------*/
void CABS(COMPLEX **a, int N, double **b)
{
  int i, j;

  for(j = 0 ; j < N ; j++) {
    for(i = 0 ; i < N ; i++) {
      b[i][j] = sqrt(a[i][j].r * a[i][j].r + a[i][j].i * a[i][j].i);
    }
  }
}

/*---------------------------------------------------------------------------*/
void mul_double_double(double **a, double **b, int N, double **c)
{
  int i, j;

  for(j = 0 ; j < N ; j++) {
    for(i = 0 ; i < N ; i++) {
      c[i][j] = a[i][j] * b[i][j];
    }
  }
}

/*---------------------------------------------------------------------------*/
void dif_double_double(double **a, double **b, int N, double **c)
{
  int i, j;

  for(j = 0 ; j < N ; j++) {
    for(i = 0 ; i < N ; i++) {
      c[i][j] = a[i][j] - b[i][j];
    }
  }
}

/*---------------------------------------------------------------------------*/
double S_double_double(double **a, double **b, int N)
{
  static int first = 1;
  static double **w;
  int i, j, ii, jj;
  double E, c, r;
  
  /***** 定数 *****/
  c = 0.7;

  /***** 重み w[][] 設定 *****/
  if (first == 1) {
    w = malloc_double_2d(N, 0, N, 0);
    for(j = 0  ; j < N ; j++) {
      if (j < N / 2) {
	jj = j;
      } else {
	jj = N - j - 1;
      }
      for(i = 0 ; i < N ; i++) {
	if (i < N / 2) {
	  ii = i;
	} else {
	  ii = N - i - 1;
	}
	r = sqrt((double)(ii * ii + jj * jj));
	w[i][j] = tanh(c * r);
      }
    }
    
    /*{    
      FILE *fp;
      
      fp = fopen("w.dat", "w");
      for(j = 0 ; j < N ; j++) {
	for(i = 0 ; i < N ; i++) {
	  fprintf(fp, "%d  %d  %f\n", i, j, w[i][j]);
	}
	fprintf(fp, "\n");
      }
      fclose(fp);
      }*/
   
      
    first = 0;
  }
  
  E = 0.0;
  for(j = 0 ; j < N ; j++) {
    for(i = 0 ; i < N ; i++) {
      E += w[i][j] * a[i][j] * b[i][j];
    }
  }

  return E;
}

/*---------------------------------------------------------------------------*/
int sol2x2(double a[][2], double b[], double x[])
{ 
  double bunbo, bunshi0, bunshi1;
    
  bunbo   = a[0][0] * a[1][1] - a[0][1] * a[1][0];
  bunshi0 = b[0] * a[1][1] - b[1] * a[0][1];
  bunshi1 = b[1] * a[0][0] - b[0] * a[1][0];

  if (fabs(bunbo) < EPS) {
    return -1;
  } else {
    x[0] = bunshi0 / bunbo;
    x[1] = bunshi1 / bunbo;
    return 0;
  }
}

/*---------------------------------------------------------------------------*/
void calc_err(double **depth, double **model, 
	      int xs, int ys, int bp, int bs, int ex_blk,
              double *mae, double *rmse)
{
  int ib, jb, xbc, ybc, hbs, nblk;
  double tmp;

  hbs = bs / 2;
  nblk = 0;
  *rmse = 0.0;
  *mae  = 0.0;
  for(jb = ex_blk, ybc = hbs + ex_blk * bp ; ybc < ys - hbs - ex_blk * bp ; ybc += bp, jb++) {
    for(ib = ex_blk, xbc = hbs + ex_blk * bp ; xbc < xs - hbs - ex_blk * bp ; xbc += bp, ib++) {
      tmp = model[ib][jb] - depth[ib][jb];
      if (fabs(tmp) < ERR_MAX) {
	*rmse += tmp * tmp;
	*mae  += fabs(tmp);
	nblk++;
      }
    }
  }
  *rmse = sqrt(*rmse / (double) nblk);
  *mae = *mae / (double) nblk;

#if DEBUG == 1
  printf("\n");
  printf("mae = %7.5f, rmse = %7.5f  (n = %d)", *mae, *rmse, nblk);
#endif
}

/*---------------------------------------------------------------------------*/
void quant(double **f, int xs, int ys, int frm)
{
  int i, j;

  for(j = -frm ; j < ys + frm ; j++) {
    for(i = -frm ; i < xs + frm ; i++) {	
      if (f[i][j] < 0) {
	f[i][j] = 0.0;
      } else if (f[i][j] > 255.0) {
	f[i][j] = 255.0;
      } else {
	f[i][j] = floorf(f[i][j] + 0.5);
      }
    }
  }
}  

/*---------------------------------------------------------------------------*/
void add_noise(double **f, int xs, int ys, int frm, double nsgm)
{
  int i, j;
  double gauss();

  for(j = -frm ; j < ys + frm ; j++) {
    for(i = -frm ; i < xs + frm ; i++) {	
      f[i][j] += (double)nsgm * gauss();
      if (f[i][j] < 0) {
	f[i][j] = 0.0;
      } else if (f[i][j] > 255.0) {
	f[i][j] = 255.0;
      } else {
	f[i][j] = floorf(f[i][j] + 0.5);
      }
    }
  }
}  

/*--------------------------------------------------------------------------*/  
double gauss(void)
{
  double rnd(void);

  return rnd() + rnd() + rnd() + rnd() + rnd() + rnd() +
    rnd() + rnd() + rnd() + rnd() + rnd() + rnd() - 6.0;
}

/*--------------------------------------------------------------------------*/  
double rnd(void)
{
  return (float)random() / (float) RAND_MAX;
}
