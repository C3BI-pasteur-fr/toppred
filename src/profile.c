/*   File:                /home/edeveaud/Work/toppred/src/profile.c
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#include <time.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <ctype.h>

#include "profile.h"
#include "error.h"
#include "output.h"

/* as sequence verification was performed before. Thus I don't care
   about checking, it is assumed all char (aa) submited will be
   correct and is comprise between correc t bounds. */
#define aa2H(c,d) (d)[(c)-'A']

/* check aa syntax from sequence  based on IUPAC alphabet */
static int is_aa(int aa);


/* retreive Hphobes datas from file */
void read_Hphobes_datas (char *file, double *hphobes_datas)
{

  int i, indx;
  char *BUFF, *p, *q, *val;
  int aa;
  double Hval;

  FILE *IN;

  /* we initialise the storage to 0 */
  for(i = 0; i < 26; i++)
    hphobes_datas[i] = 0;

  if ((IN = fopen(file, "r")) == NULL)
    error_fatal(file, NULL);

  if((BUFF = (char *)malloc(sizeof(char)*BUFFLEN)) == NULL)
    error_fatal("memory", NULL);

  if((val = (char *)malloc(sizeof(char)*(MAXSIZE + 1))) == NULL)
    error_fatal("memory", NULL);

  while (fgets(BUFF, BUFFLEN, IN) != NULL) {

    /* we skip the comment lines */
    if (*BUFF == '\0')
      continue;

    if (*BUFF == '#')
      continue;

    indx = -1;
    p = BUFF;

    /* get aa definition */
    aa = (int)*p;
    if (!is_aa(aa)){
      error_warn(file, "Unknown aminoacid.");
      continue;
    }

    p++;

    while (*p != '\0'  && isalpha((int)*p)) p++;
    while (*p != '\0' && isspace((int)*p)) p++;

    if(*p == '\0')
      error_fatal(file, "Incorrect format.");

    /* get hydrophobic value */
    q = val;
    *q = '\0';
    while (*p != '\0' && !isspace((int)*p)) {
      *q++ = *p++; }
    *q = '\0';

    /*store the datas */
    Hval = atof(val);
    indx = (aa -'A');
    hphobes_datas[indx] = Hval;
  }

  if(fclose(IN) == EOF){
    error_fatal(file, NULL);
  }

  free(BUFF);
  free(val);

  return ;
}




/* retrieve the hydropbobic calcul for each positions */
/* ___________________________________________________________________________
 *
 *  based on  G. von Heijne J. Mol. Biol. 1992 225,487-494
 *
 *       -->  l=2q+1  <--
 *          +--------+
 *     +----|        |----+
 *     +----+--------+----+
 *  -->       l=2n+1       <--
 *
 * for each given window position on sequence, i,
 * the hydrophobicity value of the window is calculated this way
 * heach hi for a given aa on the window is multiplied by Wi
 *  Wi = | 1/S              1<= i <= (n - q + 1)
 *       | (n - q + 1)/S    (n - q + 1) < i < (n + q + 1)
 *       | (2n + 2 - i)/S   (n + q + 1) <= i <= (2n + 1)
 *        with S = (1 + n)^2 - q^2
 * ___________________________________________________________________________
 *
 * rewritten in order to use n as length of the center window and q
 * as length of the wedge window
 *
 *  ie
 *       -->  l=n     <--
 *          +--------+
 *     +----|        |----+
 *     +----+--------+----+
 *  --> l=q  <--  --> l=q  <--
 *
 *  -->       l=2q+n       <--
 *
 * thus ponderation calculus become
 *
 * Wi = | 1/S                 1<= i <= q
 *      | (q + 1)/S           (q + 1) < i < (q + n)
 *      | (2q + n + 1 - i)/S  (q + n +1) <= i <= (q + 2n)
 *       with S = (q + 1)* (n + q)
 */

void calc_profile(seq_t *seq_holder, param_t params, double *profile) {

  char *seq, *p;
  int i, j, n, q;
  int R1, R2;
  int L, size;
  double S;
  double Hval;
  double *Hdatas;
  double d, wi;

  n = params.n;			/* 11  */
  q = params.q;			/*  5  */
  Hdatas = params.Hdatas;

  L = (2*q)+n;
  size = seq_holder->size - L ;
  seq = seq_holder->seq;
  p = seq;

  S = (double) ((q+1)*(n+q));

  R1 = q;
  R2 = n + q + 1;
  for(i = 0; i<= size; i++ ) { /* loop on sequence */
    Hval = 0.0;
    d = 0.0;
    p = seq+i;

    for(j = 1; j<= R1; j++) { /*loop on window 1*/
      wi = ((double) j) / S;
      d = aa2H(p[j-1], Hdatas) * wi;
      Hval = Hval + d;
    }

    wi = ((double) (q + 1)) / S;
    for(j=R1+1 ; j < R2; j++){ /*loop on window 2*/
      d = aa2H(p[j-1], Hdatas) * wi;
      Hval = Hval + d;
    }


    for(j = R2; j<= L; j++) { /*loop on window 3*/
      wi = ((double) (2 * q + n + 1 - j)) /  S;
      d = aa2H(p[j-1], Hdatas) * wi;
      Hval = Hval + d;
    }
    profile[i] = Hval;
  }
}



/* print hydrophobic values, in a useable manner */
void plot_values(double *plot, seq_t *seq, param_t params) {

  char *id, *filename, *dir;
  FILE *OUT;
  time_t timer;
  int size;
  char  *scale;
  int core_win, full_win;
  int len;
  int i;

  id = seq->id;
  size = seq->size;
  dir = params.out_dir;
  scale = params.data_file;

  core_win = params.n;
  full_win = core_win + (2 * params.q);

  /* setting the creation time */
  if (time(&timer) == (time_t) -1) error_warn("time", NULL);

  /* creating the output file-name and associated FILE */
  len =  strlen(dir) + 1 + strlen(id) + 6 + 1;
  if ((filename = (char *)malloc((size_t)len)) == NULL){
     error_fatal("memory", NULL);
   }
   (void)sprintf(filename, "%s/%s.hydro",dir, id);
   if((OUT=fopen(filename, "w")) == NULL){
     error_fatal(filename, NULL);
   }

   /* generating header */
   (void)fprintf(OUT, "# sequence:\t%s (%d res.)\n", id, size);
   (void)fprintf(OUT, "# hydrophobicity values generated: %s", ctime(&timer));
   (void)fprintf(OUT, "# core window: %d\n", core_win);
   (void)fprintf(OUT, "# full window: %d\n", full_win);
   (void)fprintf(OUT, "# using values from file: %s\n", scale);
   (void)fprintf(OUT, "# Position Hydrophobicity\n");

   for(i = 0 ; i <= size -  full_win; i++) {
     (void)fprintf(OUT, "%d\t%6.2f\n", i+1, plot[i]);
   }
   if(fclose(OUT) == EOF){
     error_fatal(filename, NULL);
   }

   free(filename);
}


/* return TRUE if c is one of the 20 amino acid known on IUPAC */
int is_aa(int aa) {
    char c;
    c = (char)aa;

    if ((c >= 'A') && (c <= 'Z'))
      if (strchr("BJOUXZ", c) == NULL)
	return TRUE;

    return FALSE;
}


#ifdef HAVE_GNUPLOT
void gplot (double *plot, segment_t **segments,  seq_t *seq, int nbr,
	     param_t params) {

  char *id, *plot_format, *plot_outfile, *dir;
  int size, x_max;
  FILE  *OUT;
  double p_cut, c_cut;
  double y_max, y_min, val;
  double *p;
  segment_t *seg;
  char tempfilename[]=TEMPFILENAME;

  char *gnuplot_command;
  int j, i, win_len;
  dir = params.out_dir;

  p_cut = params.p_cut;
  c_cut = params.c_cut;

  id = seq->id;
  size = seq->size;

  i = nbr;
  /* win_len = params.n + (2 * params.q) + 1; */
  win_len = params.n + (2 * params.q);
  seg = *segments;

  /* generate tempory file name */
  if((j = mkstemp(tempfilename)) == -1){
    error_fatal("generating temporary file", NULL);
  }

  if ((OUT = fdopen(j, "w")) == NULL) {
    error_fatal("writing temporary file", NULL);
  }

  plot_format = params.plot_format;
  plot_outfile = params.plot_outfile;

  /* getting axis range */
  y_min = y_max = val = 0;
  x_max = size - win_len + 1;
  p = plot;
  for(i = 0; i < size - win_len - 1; i++, p++) {
    val = *p;
    if (y_max < val) y_max = val;
    if (y_min > val) y_min = val;
  }
  y_min = y_min - 0.25;
  y_max = y_max + 0.5;

  /* writting gnuplot file definition */
  (void)fprintf(OUT, "set title '%s'\n", id);
  (void)fprintf(OUT, "set time\n");
  (void)fprintf(OUT, "set xlabel 'start position of window in sequence'\n");
  (void)fprintf(OUT, "set ylabel 'hydrophobicity value'\n");
  (void)fprintf(OUT, "set key right bottom\n");
  (void)fprintf(OUT, "%s\n", plot_format);
  if(plot_outfile){
    (void)fprintf(OUT, "set output \"%s/%s.%s\"\n",dir, id, plot_outfile);
  }
  /* setting plot axis */
  (void)fprintf(OUT, "set xrange [%d:%d]\n", 0, x_max);
  (void)fprintf(OUT, "set yrange [%.1f:%.1f]\n", y_min, y_max);
  (void)fprintf(OUT, "set y2range [%.1f:%.1f]\n", y_min - y_max, 0.2);

  /* tagging the tm segments */
  for(i = 0; i < nbr; i++, seg++) {
    int x;
    x = seg->pos;
    (void)fprintf(OUT, "set arrow from second %d,0 to second %d,-0.2 nohead lw 2\n", x, x);
  }

  (void)fprintf(OUT, "plot %f title 'Upper cutoff', %f title 'Lower cutoff', '%s/%s.hydro' title 'hydrophobicity' with lines\n", c_cut, p_cut, dir, id);
  (void)fprintf(OUT, "%s\n", params.plot_pause);

  if(fclose(OUT) == EOF) {
    error_fatal("gpl file 3", NULL);
  }

  /* gnuplot command */
  if((gnuplot_command = (char*)malloc(sizeof(char) *(TEMPFILENAMELEN + 8 +1))) == NULL)
    error_fatal ("memory", NULL);
  (void)sprintf(gnuplot_command, "gnuplot %s", tempfilename);

  if(system(gnuplot_command) == -1){
  	error_fatal("plotting" , NULL);
  }

    free(gnuplot_command);

    if (unlink(tempfilename) == -1) {
      error_fatal("gpl file", "deleting");
    }

}
#endif /* HAVE_GNUPLOT */
