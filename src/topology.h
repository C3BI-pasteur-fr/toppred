/*   File:                /home/edeveaud/Work/toppred/src/topology.h
 *   Author:              Katja Schuerer
 */

#ifndef __TOP_TOPOLOGY_
#define __TOP_TOPOLOGY_


#include "params.h"
#include "loop.h"


typedef struct topo {
  int putatives; /* je sais pas  */
  int kr; /* total Arg + Lys bias */
  double prob; /* global transmembrane probability */
  int cytext; /* total cyt - ext difference */
} topo_t;

typedef struct elem {
  int nputatives; /* number of putative segments */
  int nsegs;      /* number of all segements */

  seg_t *segs; /* segment structures */
  loop_t *loops;  /* loop structures */
} elem_t;

#define CYTOPLASMIC 1
#define PERIPLASMIC -1
#define INDETERMINED 0
#define MAXPUTATIVES_CALC 10

int tp_calc (topo_t *topos, elem_t *elems, seq_t *seq, param_t *para);
int tp_compare (const void *first, const void *second);

#endif





