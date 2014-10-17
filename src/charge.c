/* -----------------------------------------------------------------
 file      : /home/schuerer/toppred_com/src/top_loop.c

 author    : Schuerer
 description :


-------------------------------------------------------------------- */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <ctype.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_LIBM
#include <math.h>
#endif

#include "charge.h"
#include "error.h"


#define MAXSIZE 5

/* internal function prototypes */
static int is_aa(int c);

/* calc the Arg+Lys content of a loop

   to calculate over :
      _____           _____
   ...     |\       /|      ...
      _____|_\     /_|_____
           |<------->|                 */
int countkr (char *seq, int soff, int eoff) {

  int c = 0;
  int i;

  if (soff == -1) return 0; /* undefined sequence element */

  if (soff < 0)
    error_fatal (seq, "start position smaller than zero");
  if (eoff > strlen (seq))
    error_fatal (seq, "end position greater than sequence lenght");

  /* if (eoff - soff > 70) return 0; peut etre un niveau plus haut */

  for(i=soff; i<eoff; i++)
    if (seq[i] == 'K' || seq[i] == 'R') c++;

  return c;
}

/* calc the number of negative aa (Asp+Glu)
   in a subsequence (soff, eoff) of seq

   to calculate over :
      _____           _____
   ...     |\       /|      ...
      _____|_\     /_|_____
           |<------->|                 */
int countneg (char *seq, int soff, int eoff) {

  int c = 0;
  int i;

  if (soff == -1) return 0; /* undefined sequence element */

  if (soff < 0)
    error_fatal (seq, "start position smaller than zero");
  if (eoff > strlen (seq)){
       error_fatal (seq, "end position greater than sequence lenght");
  }
  for(i=soff; i<eoff; i++)
    if (seq[i] == 'D' || seq[i] == 'E') c++;

  return c;
}

/* charge of N-terminus =
   1 if Met start is cleaved after translation
   0 other one */
int ncharge (char *seq, int kingdom) {

  if (kingdom == EUKARYOTE) return 1;
  if (strlen(seq) > 2 && strchr("KRLFI", seq[1]) == NULL) return 1;
  return 0;
}

/* calc the difference between the cyt ext compositions
   a value of smaller than zero indicates a cytoplasmatic location
   should only be used for loops longer than 50 aa
   FEBS Lett. 303:141-46 (1992)

   to calculate over :
      _____           _____
   ...     |\       /|      ...
      _____|_\     /_|_____
              |<->|                 */
double distance (char *seq, int soff, int eoff, scale_t *scale) {

  double *cyt = scale->cyt;
  double *ext = scale->ext;
  double *sd  = scale->sd;

  int i;
  int naa[26];
  double freq, sumcyt, sumext, temp;

  if (soff == eoff) return 0.0;

  for (i=0; i<26; i++) naa[i] = 0;
  for (i=soff; i<eoff; i++) naa[seq[i] - 'A']++;

  sumext = sumcyt = 0.0;
  for (i=0; i<26; i++) {
    if (is_aa (i + 'A')) {
      freq = naa[i] * 100.0 / (eoff - soff);
      temp = (cyt[i] - freq) / sd[i];
      sumcyt += temp * temp;
      temp = (ext[i] - freq) / sd[i];
      sumext += temp * temp;
    }
  }

  return sqrt(sumcyt) - sqrt(sumext);
}


/* calc the net charge difference of the N-terminal segment including 15 aa
   up- and 15 aa downward
   Hartmann et al., PNAS

   to calculate over :
                      __________
                    /|    S0    |\
          | 15 aa |/_|__________|_\| 15 aa |
          |<-------->|          |<-------->|    */
int nterminus (char *seq, int soff, int eoff, int delta) {

  int start, end, charge;
  int i, len;

  charge = 0;
  /* calc charge avant the segment */
  start = soff - 15;
  start = (start < 0) ? 0 : start;
  end = soff + delta;

  for (i=start; i<end; i++)
    if (strchr("KR", seq[i]) != NULL) charge++;
    else if (strchr("DE", seq[i]) != NULL) charge--;

  /* calc charge after the segment */
  start = eoff - 1 - delta;
  end = eoff - 1 + 15;
  len = strlen(seq);
  end = (end > len) ? len : end;

  for (i=start; i<end; i++)
    if (strchr("KR", seq[i]) != NULL) charge--;
    else if (strchr("DE", seq[i]) != NULL) charge++;

  return charge;
}

static int is_aa(int aa) {

    char c;
    c = (char)aa;

    if ((c >= 'A') && (c <= 'Z'))
      if (strchr("BJOUXZ", c) == NULL)
        return TRUE;

    return FALSE;

}

/* retreive cyt-ext scale datas from file */
int read_cytext_datas (char *file, scale_t *scales) {

  double *cyt = scales->cyt;
  double *ext = scales->ext;
  double *sd = scales->sd;

  int i;
  char *BUFF, *p, *q, *val;
  int aa;

  FILE *IN;

  /* we initialise the storage to 0 */
  for(i = 0; i < 26; i++) {
    cyt[i] = ext[i] = sd[i] = 0;
  }

  if ((IN = fopen(file, "r")) == NULL)
    error_fatal(file, NULL);

  if((BUFF = (char *)malloc(BUFFLEN)) == NULL)
    error_fatal("memory", NULL);

  if((val = (char *)malloc(sizeof(char) *  (MAXSIZE + 1))) == NULL)
    error_fatal("memory", NULL);

  while (fgets(BUFF, BUFFLEN, IN) != NULL) {

    /* we skip the comment lines */
    if (*BUFF == '\0')
      continue;

    if (*BUFF == '#')
      continue;

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


    /* get cyt value */
    q = val;
    *q = '\0';
    while (*p != '\0' && !isspace((int)*p)) {
      *q++ = *p++; }
    *q = '\0';
     cyt[aa - 'A'] = atof(val);

     while (*p != '\0' && isspace((int)*p)) p++;
     /* get ext value */
     q = val;
     *q = '\0';
     while (*p != '\0' && !isspace((int)*p)) {
       *q++ = *p++; }
     *q = '\0';
     ext[aa - 'A'] = atof(val);

     while (*p != '\0' && isspace((int)*p)) p++;
     /* get ext value */
     q = val;
     *q = '\0';
     while (*p != '\0' && !isspace((int)*p)) {
       *q++ = *p++; }
     *q = '\0';
     sd[aa - 'A'] = atof(val);

  }

  if(fclose(IN) == EOF)
    error_fatal(file, NULL);

  free(BUFF);
  free(val);

  /* print_Hphobes_datas(hphobes_datas); */

  return 0;
}

/* determine the N-terminus orientation */
char *orientation(double val) {
  if (val > 0.0) return "N-in";
  else if (val < 0.0) return "N-out";
  else return "undecided";
}

/* determine the N-terminus orientation */
char *new_orientation(double val) {
  if (val > 0.0) return "N-in";
  else if (val < 0.0) return "N-out";
  else return "?";
}
