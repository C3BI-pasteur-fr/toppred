/* -----------------------------------------------------------------
 file      : /home/schuerer/toppred_com/src/top_topology.c
 author    : Schuerer
 description :


-------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#ifdef HAVE_LIBM
#include <math.h>
#endif

#include "error.h"

#include "topology.h"

#include "charge.h"


/* calculation of the total kr bias or each possible topology
   stockage of the para->tmax best topologies */
int tp_calc (topo_t *topos, elem_t *elems, seq_t *seq_h, param_t *para) {

  int tmax = para->n_topos;
  int cllen = para->seg_len;

  int nsegments = elems->nsegs;
  seg_t *segments = elems->segs;
  loop_t *loops = elems->loops;
  int ntopos = (int) pow(2.0, (double) elems->nputatives);

  int i, k, n, saved;
  int kr, kr_l, start_l, side, first_tm;
  double probtm;
  int putatives, integrate, kingdom;
  char * seq;

  seq = seq_h->seq;

  kingdom = para->kingdom;

  saved = 0; k = 0;
  for (n=0; n < ntopos; n++) {
    putatives = n; side = 1;
    /* init kr with 1 if N-terminal met is cleaved elsewhere 0 */
    kr = 0;
    kr_l = 0; start_l = 0;
    probtm = 1.0;
    first_tm = 1;
    for (i=0; i < nsegments; i++) {

      /* segmentment integrated ? */
      if (segments[i].kind == PUTATIVE) {
	integrate = putatives % 2;
	putatives /= 2;
      }
      else integrate = 1;

      /* calculation of global transmembrane probability */
      if (integrate) probtm *= segments[i].probTM;

      /* calculation of total kr bias */
      if (loops[i].start != -1) kr_l += loops[i].delta;
      if(integrate) {
	if (loops[i].start != -1 && loops[i].stop - start_l <= cllen) {
	  kr += kr_l * side;
	  if (first_tm) kr += ncharge(seq, kingdom);
	}

	side *= (-1);
	kr_l = 0;
	start_l = segments[i].stop;
	first_tm = 0;
      }
      else kr_l += segments[i].delta;

    } /* end for to calculate total kr bias */
    if (loops[nsegments].start == -1 &&
	segments[nsegments-1].stop - start_l <= cllen )
      kr += kr_l * side;
    else if (loops[nsegments].stop - start_l <= cllen)
      kr += (kr_l + loops[nsegments].delta) * side;

    /* if there is at least one segment save the topology */
    if (!first_tm) {
      /* find topology with minimum kr bias to replace */
      if (saved >= tmax) {
	k = 0;
	for (i=1; i<tmax; i++)
	  if (abs(topos[i].kr) < abs(topos[k].kr)) k = i;
      }
      else { k = saved; saved++;}
      /* init topology object */
      topos[k].putatives = n;
      topos[k].kr    = kr;
      topos[k].prob  = probtm;
    }
  }

  return saved;
}

/* sort topologies in increasing order */
int tp_compare (const void *first, const void *second) {
  const topo_t *topof = (const topo_t *) first;
  const topo_t *topos = (const topo_t *) second;

  int res = abs(topos->kr) - abs(topof->kr);

  if (res) { return (abs(topos->kr) - abs(topof->kr)); }
  if (topof->prob > topos->prob) return -1;
  return 1;
}


/* extract loops and segments as 2 distinct structures */
/*
void KS_wrap (KSsegment_t *segs, KSloop_t *loops,
 	      boucle_t *boucle, int n) {
   int i, s, l;
   boucle_t *p;

   KSseg_t seg;
   KSloop_t loop;

   i = s = l = 0;
   seg = NULL;
   loop = NULL;

   P = boucle;
    getting the sizes for both struct to create
   for(i = 0; i< n; i++) {
     if(p[i].kind == PUTATIVE || p[i].kind == CERTAIN) {
       s++;
     }
     else {
       l++;
     }
   }

    allocate both structures
   if ((seg = (KSseg_t)malloc(sizeof(KSseg_t) *s)) == NULL) {
     error_fatal("memory", NULL);
   }
   if ((loop = (KSloop_t)malloc(sizeof(KSloop_t) *s)) == NULL) {
     error_fatal("memory", NULL);
   }
   s = l = 0;
   for(i = 0; i < n; i++) {
     if(p[i].kind == LOOP) {
       loop[l].start = p[i].start;
       loop[l].stop = p[i].stop;
       loop[l].charge = p[i].delta;
       l++;
     }  end if LOOP

     else {  it's a segment...
       seg[s].start = p[i].start;
       seg[s].stop = p[i].stop;
       seg[s].kind = p[i].kind;
       seg[s].hpval = p[i].H;
       seg[s].charge = p[i].delta;
     }

   }  end for

   segs = seg;
   loops =loop;

}
*/
