/*   File:                /home/schuerer/toppred/src/output.h
 *   Author:              Schuerer
 */

#ifndef __OUTPUT_H_
#define __OUTPUT_H_

#include <stdlib.h>
#include "loop.h"

void new_print_parameters (param_t *p);
void old_print_firstparameters (param_t *p);
void old_print_secondparameters (param_t *p);

void print_sequence (FILE *OUT, char *seq);
void print_tmsummary (FILE *OUT, int nseg, seg_t *KSseg);
void old_print_tmsummary (FILE *OUT, int nseg, seg_t *KSseg,char *seq);

#endif /* __OUTPUT_H_ */
