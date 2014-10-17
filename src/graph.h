/*   File:                /home/edeveaud/Work/toppred/src/top-graph.h
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 */


#ifndef __GRAPH_H__
#define __GRAPH_H__


#include "topoprint.h"

				
#define IMAGE_LENGTH 800           /* image size is fixed */
#define IMAGE_HEIGTH 400
#define Y_CENTER IMAGE_HEIGTH/2
#define X_CENTER IMAGE_LENGTH/2
#define SEG_L 20                     /* segment representation is fixed */
#define SEG_H 60
#define LOOP_H 40
#define MARGIN 20
#define MAX_SEGMENTS 37		/* maximum segments we can represent */

#define CYT -1
#define MBR 0
#define EXT +1

#define FIRST -1
#define LAST   1
#define INNER 0

void topo_graph_print (topoprint_t *topo,  int n, param_t *params, seq_t *seq);

#endif /* __GRAPH_H__ */
