#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <imageio.h>
#include <nd_malloc.h>
#include "fgcrdfd.h"

/*---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  FILE *fp;
  char com[256], fullbase[256], outbase[256], base[256], flnm[256], c, buf[256], *cp;
  int bp, bs, hbs, frm;
  int xs, ys, nz, nbx, nby, xbc, ybc;
  int quant_flg, noise_flg, wgt_flg, out_flg, prfn, exblk;
  int ib, jb, z;
  double noise_sgm, za, zb, lens_k, dz, tdini, tkini, tdest, tkest, mae, rmse;
  double ***inp, ***preinp, ***bpfinp, **pmodel, **model, **dini, **kini, **dest, **kest;
  int **maxz, tmaxz;

  /***** 引数の処理 *****/
  strcpy(com, argv[0]);
  bp  = BP_DEFAULT;
  bs = BS_DEFAULT;
  quant_flg = 0;
  noise_flg = 0;
  wgt_flg = 0;
  out_flg = 0;
  noise_sgm = NOISE_SGM_DEFAULT;
  exblk = EXBLK_DEFAULT;
  prfn = PRFN_DEFAULT;
  while( (c = getopt(argc, argv, "f:b:p:s:qn:R:wovhH") ) != -1) {
    switch(c) {
    case 'f': /* prefilter 番号 */
      prfn = atoi(optarg);
      break;
    case 'b': /* 推定ブロクサイズ(odd) */
      bs = atof(optarg); 
      break;
    case 'p': /* 推定ブロックピッチ */
      bp = atof(optarg); 
      break;
    case 'q': /* 量子化 */
      quant_flg = 1;
      break;
    case 'n': /* 観測雑音と sgm */
      noise_flg = 1;
      noise_sgm = atof(optarg);
      break;
    case 'R': /* 誤差測定に除外する外側推定点の個数 */
      exblk = atoi(optarg);
      break;
    case 'w': /* weight */
      wgt_flg = 1;
      break;
    case 'o': /* データ出力 */
      out_flg = 1;
      break;
    case 'v':  /* version 表示 */
      fprintf(stderr, "\n%s  version: %s\n\n", com, VER);
      exit(1);
      break;
    case 'h':  /* ヘルプ */
    case 'H':
    default:
      usage(com);
    }
  }
  argc -= optind - 1;
  argv += optind - 1;
  if (argc < 3) usage(com);
  if (bs % 2 == 0) usage(com);

  /***** 初期化 *****/
  strcpy(fullbase, argv[1]);
  strcpy(outbase, argv[2]);
  hbs = bs / 2;
  frm = BPF_HS;
  for(cp = argv[1] ; *cp != '\0' ; cp++);  /* path 情報除去 */
  for( ; *cp != '/' ; cp--);
  strcpy(base, ++cp);

  /***** データファイル読み込み *****/
  sprintf(flnm, "%s.dat", fullbase);
  if ( (fp = fopen(flnm, "r")) == NULL) {
    fprintf(stderr, "\nError : %s cannot be found\n\n", flnm); exit(1);
  }
  fgets(buf, 256, fp);  /* XS */
  strtok(buf, " \t\n");
  xs = atoi(strtok(NULL, " \t\n"));
  fgets(buf, 256, fp);  /* YS */
  strtok(buf, " \t\n");
  ys = atoi(strtok(NULL, " \t\n"));
  fgets(buf, 256, fp);  /* ZA */
  strtok(buf, " \t\n");
  za = atof(strtok(NULL, " \t\n"));
  fgets(buf, 256, fp);  /* ZB */
  strtok(buf, " \t\n");
  zb = atof(strtok(NULL, " \t\n"));
  fgets(buf, 256, fp);  /* NZ */
  strtok(buf, " \t\n");
  nz = atoi(strtok(NULL, " \t\n"));
  fgets(buf, 256, fp);  /* LENS_K */
  strtok(buf, " \t\n");
  lens_k = atof(strtok(NULL, " \t\n"));
  fclose(fp);

  /***** 諸定数 *****/
  dz = (zb - za) / (double)(nz - 1);
  for(nbx = 0, xbc = hbs ; xbc < xs - hbs ; xbc += bp, nbx++);
  for(nby = 0, ybc = hbs ; ybc < ys - hbs ; ybc += bp, nby++);
  
#if DEBUG == 1  
  printf("\n");
  printf("(xs, ys)  = (%d, %d)\n", xs, ys);
  printf("bp        = %d\n", bp);
  printf("bs(hbs)   = %d(%d)\n", bs,hbs);
  printf("frm       = %d\n", frm);
  printf("flgs      = (quant, noise, block_weight, out)=(%d, %d, %d, %d)\n", 
	 quant_flg, noise_flg, wgt_flg, out_flg);
  printf("noise sgm = %5.3f\n", noise_sgm);
  printf("lens_k    = %5.3f\n", lens_k);
  printf("za -- zb  = %5.3f -- %5.3f\n", za, zb);
  printf("nz        = %d\n", nz);
  printf("dz        = %5.3f\n", dz);
  printf("exblk     = %d\n", exblk);
  printf("\n");
#endif

  /***** 領域確保 *****/
  inp    = malloc_double_3d(nz, 0, 2*frm + xs + 2*frm, 2*frm, 2*frm + ys + 2*frm, 2*frm);
  preinp = malloc_double_3d(nz, 0, frm + xs + frm, frm, frm + ys + frm, frm);
  bpfinp = malloc_double_3d(nz, 0, xs, 0, ys, 0);
  pmodel = malloc_double_2d(xs, 0, ys, 0);
  model  = malloc_double_2d(nbx, 0, nby, 0);  
  dini   = malloc_double_2d(nbx, 0, nby, 0);
  kini   = malloc_double_2d(nbx, 0, nby, 0);
  dest   = malloc_double_2d(nbx, 0, nby, 0);
  kest   = malloc_double_2d(nbx, 0, nby, 0);
  maxz   = malloc_int_2d(nbx, 0, nby, 0);
  
  /***** 形状モデル生成 *****/  
  get_model_bin(fullbase, model, xs, ys, bp, bs);
  save_model(model, xs, ys, bp, bs, outbase);

  /***** 多焦点画像群読み込み *****/  
  for(z = 0 ; z < nz ; z++) {
    sprintf(flnm, "%s_ra_%02d.bin", fullbase, z);
    get_bin(flnm, inp[z], xs, ys, frm);
  } 

  /***** 量子化 ＆ 観測雑音 *****/  
  if (noise_flg == 1) {
    if (quant_flg == 1) {
      for(z = 0 ; z < nz ; z++) {
	add_noise(inp[z], xs, ys, frm, noise_sgm);
	quant(inp[z], xs, ys, frm);
      } 
    } else {
      for(z = 0 ; z < nz ; z++) {
	add_noise(inp[z], xs, ys, frm, noise_sgm);
      }
    } 
  } else {
    if (quant_flg == 1) {
      for(z = 0 ; z < nz ; z++) {
	quant(inp[z], xs, ys, frm);
      } 
    }
  }

  /***** pre filter 適用 *****/
  for(z = 0 ; z < nz ; z++) {
    pref(inp[z], xs, ys, frm, prfn, preinp[z]);
  }

  /***** mzdfd による初期 depth 推定用の BPF 適用 *****/
  for(z = 0 ; z < nz ; z++) {
    mzbpf(inp[z], xs, ys, 0, bpfinp[z]);
  }

  /***** ブロック単位に depth 推定  *****/
  for(jb = 0, ybc = hbs ; ybc < ys - hbs ; ybc += bp, jb++) {
    for(ib = 0, xbc = hbs ; xbc < xs - hbs ; xbc += bp, ib++) {
      mzdfd(bpfinp, nz, xbc, ybc, bs, &tdini, &tkini, &tmaxz);
      /*** 初期 depth 推定 ***/
      dini[ib][jb] = tdini;
      kini[ib][jb] = tkini;
      maxz[ib][jb]  = tmaxz;
      /*** fgcrdfd で depth 推定 ***/
      fgcrdfd(preinp, tmaxz, xbc, ybc, bs, tdini, tkini, wgt_flg, &tdest, &tkest);
//      fgcrdfd(inp, tmaxz, xbc, ybc, bs, tdini, tkini, wgt_flg, &tdest, &tkest);
      dest[ib][jb] = tdest;
      kest[ib][jb] = tkest;
    }
  }
  fprintf(stderr, "\n");
  save_depth(dest, model, dini, kest, xs, ys, bp, bs, exblk, outbase);

  /***** 誤差測定 *****/
  calc_err(model, dest, xs, ys, bp, bs, exblk, &mae, &rmse);  
  
  /***** 高さ表示の gp ファイル保存 *****/
  save_gp(outbase);

  /***** データファイル保存 *****/
  save_data(com, fullbase, lens_k, bs, bp, 
	    quant_flg, noise_flg, wgt_flg, noise_sgm, exblk, mae, rmse, outbase);

  return 0;
}

/*---------------------------------------------------------------------------*/
void usage(char *com)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "  Usage : %s  [-fbpqnRwovhH]  <full_path_md_base_name>  <outbase>\n\n", com);
  fprintf(stderr, "     -f <#filter>        : prefilter number (%d)\n", PRFN_DEFAULT);
  fprintf(stderr, "     -b <blk size (ODD)> : block size (odd) (%d)\n", BS_DEFAULT);
  fprintf(stderr, "     -p <blk pitch>      : block pitch (%d)\n", BP_DEFAULT);
  fprintf(stderr, "     -q                  : 8 bit quantize\n");
  fprintf(stderr, "     -n <noise sgm>      : noise added and sgm \n");
  fprintf(stderr, "     -R <Nblk>           : outer blks excluded in RMS calc.\n");
  fprintf(stderr, "     -w                  : block weight on\n");
  fprintf(stderr, "     -o                  : Data & Image output\n");
  fprintf(stderr, "     -v                  : shows ver.\n");
  fprintf(stderr, "     -h, -H              : shows this\n");
  fprintf(stderr, "\n");

  exit(1);
}

