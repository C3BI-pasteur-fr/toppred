/*   File:                /home/edeveaud/Work/toppred/src/charge.h
 *   Author:              Katja Schuerer
 */


#ifndef __CHARGE_H_
#define __CHARGE_H_


#include "params.h"

#define EUKARYOTE 1
#define PROKARYOTE 2

/* retreive cyt-ext scale datas from file */
int read_cytext_datas (char *file, scale_t *scales);

/* calc the Arg+Lys content over a subsequence from soff to eoff */
int countkr (char *seq, int soff, int eoff);

/* calc the Asp+Glu content over a subsequence from soff to eoff */
int countneg (char *seq, int soff, int eoff);

/* calc the charge of the N-terminus, depend on cleavage of Met start */
int ncharge (char *seq, int kingdom);

/* calc the difference of cyt-ext values of a subsequence */
double distance (char *seq, int soff, int eoff, scale_t *scale);

/* calc the charge difference over adjacent loops of the first segment
   including 15 aa at each side */
int nterminus (char *seq, int soff, int eoff, int delta);

/* determine the N-terminus orientation */
char *orientation(double val);

/* determine the N-terminus orientation
   if undecided return "?" */
char *new_orientation(double val);

#endif /* __CHARGE_H_ */
