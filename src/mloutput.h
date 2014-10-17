/*   File:                /home/schuerer/toppred/src/mloutput.h
 *   Author:              Schuerer
 */

#ifndef __MLOUTPUT_H_
#define __MLOUTPUT_H_

#include <stdlib.h>
#include "loop.h"
#include "params.h"
#include "topoprint.h"
#include "seq-reader.h"

/* html output functions */
#define start_phrase() (void) fprintf(OUT, "<PRE><P>\n")
#define end_phrase()   (void) fprintf(OUT, "</P></PRE>\n")

FILE *init_html(char *name, char *dir);
void html_header (FILE *OUT, char *name);
void html_parameters (FILE *OUT, param_t *p);

void html_sequence (FILE *OUT, seq_t *seq);
void html_plot (FILE *OUT, char *seqname, param_t *p);
void html_topology (topoprint_t *topo, int nel, param_t *para);

#endif /* _MLOUTPUT_H_ */
