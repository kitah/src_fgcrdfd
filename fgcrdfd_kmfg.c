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

  /***** $B0z?t$N=hM}(B *****/
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
    case 'f': /* prefilter $BHV9f(B */
      prfn = atoi(optarg);
      break;
    case 'b': /* $B?dDj%V%m%/%5%$%:(B(odd) */
      bs = atof(optarg); 
      break;
    case 'p': /* $B?dDj%V%m%C%/%T%C%A(B */
      bp = atof(optarg); 
      break;
    case 'q': /* $BNL;R2=(B */
      quant_flg = 1;
      break;
    case 'n': /* $B4QB,;(2;$H(B sgm */
      noise_flg = 1;
      noise_sgm = atof(optarg);
      break;
    case 'R': /* $B8m:9B,Dj$K=|30$9$k30B&?dDjE@$N8D?t(B */
      exblk = atoi(optarg);
      break;
    case 'w': /* weight */
      wgt_flg = 1;
      break;
    case 'o': /* $B%G!<%?=PNO(B */
      out_flg = 1;
      break;
    case 'v':  /* version $BI=<((B */
      fprintf(stderr, "\n%s  version: %s\n\n", com, VER);
      exit(1);
      break;
    case 'h':  /* $B%X%k%W(B */
    case 'H':
    default:
      usage(com);
    }
  }
  argc -= optind - 1;
  argv += optind - 1;
  if (argc < 3) usage(com);
  if (bs % 2 == 0) usage(com);

  /***** $B=i4|2=(B *****/
  strcpy(fullbase, argv[1]);
  strcpy(outbase, argv[2]);
  hbs = bs / 2;
  frm = BPF_HS;
  for(cp = argv[1] ; *cp != '\0' ; cp++);  /* path $B>pJs=|5n(B */
  for( ; *cp != '/' ; cp--);
  strcpy(base, ++cp);

  /***** $B%G!<%?%U%!%$%kFI$_9~$_(B *****/
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

  /***** $B=tDj?t(B *****/
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

  /***** $BNN0h3NJ](B *****/
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
  
  /***** $B7A>u%b%G%k@8@.(B *****/  
  get_model_bin(fullbase, model, xs, ys, bp, bs);
  save_model(model, xs, ys, bp, bs, outbase);

  /***** $BB?>GE@2hA|72FI$_9~$_(B *****/  
  for(z = 0 ; z < nz ; z++) {
    sprintf(flnm, "%s_ra_%02d.bin", fullbase, z);
    get_bin(flnm, inp[z], xs, ys, frm);
  } 

  /***** $BNL;R2=(B $B!u(B $B4QB,;(2;(B *****/  
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

  /***** pre filter $BE,MQ(B *****/
  for(z = 0 ; z < nz ; z++) {
    pref(inp[z], xs, ys, frm, prfn, preinp[z]);
  }

  /***** mzdfd $B$K$h$k=i4|(B depth $B?dDjMQ$N(B BPF $BE,MQ(B *****/
  for(z = 0 ; z < nz ; z++) {
    mzbpf(inp[z], xs, ys, 0, bpfinp[z]);
  }

  /***** $B%V%m%C%/C10L$K(B depth $B?dDj(B  *****/
  for(jb = 0, ybc = hbs ; ybc < ys - hbs ; ybc += bp, jb++) {
    for(ib = 0, xbc = hbs ; xbc < xs - hbs ; xbc += bp, ib++) {
      mzdfd(bpfinp, nz, xbc, ybc, bs, &tdini, &tkini, &tmaxz);
      /*** $B=i4|(B depth $B?dDj(B ***/
      dini[ib][jb] = tdini;
      kini[ib][jb] = tkini;
      maxz[ib][jb]  = tmaxz;
      /*** fgcrdfd $B$G(B depth $B?dDj(B ***/
      fgcrdfd(preinp, tmaxz, xbc, ybc, bs, tdini, tkini, wgt_flg, &tdest, &tkest);
//      fgcrdfd(inp, tmaxz, xbc, ybc, bs, tdini, tkini, wgt_flg, &tdest, &tkest);
      dest[ib][jb] = tdest;
      kest[ib][jb] = tkest;
    }
  }
  fprintf(stderr, "\n");
  save_depth(dest, model, dini, kest, xs, ys, bp, bs, exblk, outbase);

  /***** $B8m:9B,Dj(B *****/
  calc_err(model, dest, xs, ys, bp, bs, exblk, &mae, &rmse);  
  
  /***** $B9b$5I=<($N(B gp $B%U%!%$%kJ]B8(B *****/
  save_gp(outbase);

  /***** $B%G!<%?%U%!%$%kJ]B8(B *****/
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

