/*   File:                /home/edeveaud/Work/toppred/src/topoprint.h
 *   Author:              Katja Schuerer
 */

#ifndef __TOP_TOPOPRINT_
#define __TOP_TOPOPRINT_

#include "params.h"

#include "topology.h"

/* types for full length topology  printing */

typedef struct tmprint
{
  int nr;
  int start;
  int stop;
  int len;
  double prob;
  double H;
  char *type;
} tmprint_t;

typedef struct loprint
{
  int start;
  int stop;
  int len;
  int kr;
  double cytext;
  char *type;
} loprint_t;

typedef union elprint
{
  loprint_t lo;
  tmprint_t tm;
} elprint_t;

typedef struct topoprint
{
  double prob;
  double cytext;
  int pos;
  int neg;
  double negpos;
  int nr;
  int putatives;
  int kr;
  int nterm;
  elprint_t *elps;
  char *image;
  int ncharge;
} topoprint_t;

/* transforme topologie as topo_t object to topology_t with
   informatigons of all elements  */
int tp_decode (topoprint_t *topo, elem_t *elems, seq_t *sq, param_t *para);

void tp_toppred_fprintf (topoprint_t *topo, int nel, param_t *para);

/* print topology information in new output format */
void tp_new_print (topoprint_t *topo, int nel, param_t *para);

#endif


