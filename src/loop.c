/*   File:                /home/edeveaud/Work/toppred/src/loop.c
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "loop.h"
#include "error.h"



/* eliminate segments overlapping */
static int purge_segments(segment_t **seg, param_t params, int nbr);

/* sort segments by Hvalues, used by purge_segments*/
static  int seg_H_sort(const void *p1, const void *p2);
/* sort segments by positions, used by purge_segments*/
static int seg_pos_sort(const void *p1, const void *p2);

/* sort segments regardin their Hydrophobicity  */
static int seg_H_sort(const void *p1, const void *p2) {
  const segment_t *seg1, *seg2 ;

  seg1 = (const segment_t *)p1;
  seg2 = (const segment_t *)p2;

  if(seg1->max > seg2->max)
    return -1 ;

  if(seg1->max < seg2->max)
    return 1;

  return 0 ;
}


/* get segments from the hydrophobicity profile
 * extract all maximum whose value is greater than the putative cut-of value
 * on the curve, and check for their compliance to the
 * parameters, ie eliminate overlaping segments
 * allocte the correspondig structure, and returns the number of "correct"
 * segments found
 */
int get_segments(double *Hplot, segment_t **res,int nb, param_t params)  {

  int n;
  int i, n_tot;
  size_t len;
  double p_cut;
  double prev, curr, next;

  segment_t curr_seg;
  segment_t *pos, *p;

  p_cut = params.p_cut;
  len = sizeof(segment_t);
  prev = curr = next = 0.0;
  pos = p = NULL;

  if ((pos = (segment_t *)malloc(MAXSEGMENT*len)) == NULL)
    error_fatal("Memory", NULL);

  n_tot = MAXSEGMENT;

  p = pos;
  i = 0;
  n = 0;
  while(i < nb ) {
    /* take segment at beginning of the sequence */
    if (i == 0) prev = params.p_cut - 1.0;
    else prev = Hplot[i-1];
    curr = Hplot[i];
    /* take segment at end of the sequence */
    if (i == nb-1) next = curr - 1.0;
    else next = Hplot[i+1];
    /* <= and >= allow to get segments in position 0 */
    /* and to get segments at the begining of a "plateau" */
    if((curr > p_cut) && ((prev <= curr) && (curr >= next))) {
      curr_seg.max = curr;
      curr_seg.pos = i;
      curr_seg.keep = 1;
      /* resize the segment holder if needed */
      if (n >= n_tot) {
        n_tot =  n_tot + MAXSEGMENT;
        if ((pos=(segment_t *)realloc(pos, n_tot*len)) == NULL)
          error_fatal("Memory", "Reallocating segment holder");
        p = pos + n;
      }
      *p=curr_seg;
      n++;
      p++;
    } /* end if max */
    i++;
  } /* end while */

  if (n != 0) {
    qsort(pos, (size_t)n, len, seg_H_sort) ;
    n = purge_segments(&pos, params, n);

    /* allocating the real segment holder */
    if ((pos = (segment_t *)realloc(pos, n*len)) == NULL)
      error_fatal("Memory", "resizing results");
  }

  *res = pos;
  return n;
}


/* clean all the "incorretc" segments,
 * overlaping ones,
 * contiguous ones */
static int purge_segments(segment_t **seg, param_t params, int nbr) {

  int i, j;
  int n, len;
  int pos1, pos2;
  int tmspacer;
  segment_t *p, *q;

  n = 0;
  tmspacer = params.tmspacer;

  /* don't allow transmembrane segments to be contiguous*/
  len = (2* params.q) + params.n + tmspacer;
  i = j = 0;
  p = q = * seg;
  for (i = 0; i <nbr; i++) {
    pos1 = p->pos;
    q = *seg;
    for (j=0; j< nbr; j++) {
      pos2 = q->pos;

      if(p->keep == FALSE){
	q++;
	continue;
      }

      if((pos2 < pos1) && (pos1 - len < pos2)){
	q->keep = FALSE;
	q++;
	continue;
      }

      if((pos2 > pos1) && (pos1 + len > pos2)) {
	q->keep = FALSE;
	q++;
	continue;
      }
      q++;
    }/*end for j */

    p++;
  }/*end for i */

  p = *seg;
  n = 0;
  for(i = 0; i < nbr; i++, n++, p++){
    if(p->keep == 0){
      p->pos = 9999999;
      n--;
    }
  }
  q = *seg;
  qsort(q, (size_t)nbr, sizeof(segment_t), seg_pos_sort) ;

  return n;
}


static int seg_pos_sort(const void *p1, const void *p2) {
  const segment_t *seg1, *seg2;

  seg1 = (const segment_t *)p1;
  seg2 = (const segment_t *)p2;

  if(seg1->pos > seg2->pos)
    return 1 ;

  if(seg1->pos < seg2->pos)
    return -1;

  return 0 ;

}

/* allocate seg_t and loop_t structures from the segment_t structure
 * and returns the number of loops */
int calc_loop(seq_t *sequence, segment_t **segment, loop_t **loopKS,
	       seg_t **segKS, param_t params, int nbr)
{
  loop_t *l_res, *l;
  seg_t  *s_res, *s;
  segment_t *seg;

  int i, n, q_win;
  int start, stop, x, y;
  int state;
  int win_len;
  int seq_len;
  char *seq;

  start = stop = 0;

  win_len = (2*params.q) + params.n;  /* full window length  */
  q_win = params.q;		      /* flanking region length  */

  seq = sequence->seq;
  seq_len = sequence->size;
  seg = *segment;
  state = -1;

  n = nbr+1;

  if ((l_res = (loop_t *)malloc(n*sizeof(loop_t))) == NULL){
    error_fatal("memory", "allocating loops");
  }
  if ((s_res = (seg_t *)malloc(nbr*sizeof(seg_t))) == NULL) {
    error_fatal("memory", "allocating segments");
  }
  l = l_res;
  s = s_res;

  if(seg[0].pos == 0) { /* sequence starts with a segment */
    l[0].start = -1;
    l++;
    state = 1;
  }

  for(i = 0; i < nbr; i++, seg++) {

    /* we are in a loop */
    if(state == -1){
      x = start = stop;
      y = stop = seg->pos;

      l->start = start;
      l->stop  = stop;
      if ((start - q_win) >= 0){
	x = start - q_win;
      }
      if ((stop + q_win) <= seq_len) {
	y = stop + q_win;
      }
      l->delta = countkr(seq, x, y);
      state *= -1;
      l++;
    }
    /* we are now in a menbrane segment */
    if(state == 1);{
      start = stop;
      stop = start + win_len;
      s->start = start;
      s->stop = stop;
      s->H = seg->max;
      if (seg->max >= params.c_cut) {
	s->kind = CERTAIN; /* 1 */
	s->probTM = 1.0;
      }
      else {
	s->kind = PUTATIVE; /* 0 */
	s->probTM = (seg->max - params.p_cut) / (params.c_cut - params.p_cut);
      }
      x = start + q_win;
      y =  stop - q_win;
      s->delta = countkr(seq, x, y);
      s++;
      state *= -1;
    }
  }

  /* dealing with the last segment if it exists */
  if (stop < seq_len ) {
    start = stop;
    stop = sequence->size;
    l->start = start;
    l->stop = stop;
    l->delta = countkr(seq, start - q_win, stop);
  }
  else { /* sequence ends with a segment */
    l->start = -1;
  }

  *segKS = s_res;
  *loopKS = l_res;

  return n++;
}
