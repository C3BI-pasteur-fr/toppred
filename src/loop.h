/*   File:                /home/edeveaud/Work/toppred/src/loop.h
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 */


#ifndef __LOOP_H_
#define __LOOP_H_



#include "params.h"
#include "charge.h"
#include "seq-reader.h"



typedef struct segment_S {
  double max;
  int pos;
  int stop;
  int keep;
}segment_t;


typedef struct  seg_S {
  int start;
  int stop;
  int kind;
  int delta;
  double H;
  double probTM;
}  seg_t;


typedef struct loop_S {
  int start;
  int stop;
  int delta;
} loop_t;

#define MAXSEGMENT 20
#define TMSPACER 2
#define PUTATIVE 0
#define CERTAIN 1
#define LOOP -1

#define MAXSIZE 5

/* retrieve the position of ech peak on the curve, with is associated
   H value  */
int get_segments(double *Hplot, segment_t **res,int nb, param_t params);

int calc_loop(seq_t *sequence, segment_t **segment, loop_t **loop,
	       seg_t **seg, param_t params, int nbr);
#endif /* __LOOP_H_ */
