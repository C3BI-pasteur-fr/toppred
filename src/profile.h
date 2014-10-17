/*   File:                /home/edeveaud/Work/toppred/src/profile.h
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 */


#ifndef __PROFILE_H_
#define __PROFILE_H_

#include "params.h"
#include "seq-reader.h"
#include "loop.h"



/* WARNING if changing TEMPFILENAME value, adjust TEMPFILENAMELEN */
#define TEMPFILENAME "/tmp/top-XXXXXX"
#define TEMPFILENAMELEN 15



/* retreive Hphobes datas from file */
void read_Hphobes_datas (char *file, double *hphobes_datas);


/* print hydrophobic values, in a useable manner */
void plot_values(double *plot, seq_t *seq, param_t params);

/* retrieve the hydropbobic calcul for each positions */
void calc_profile(seq_t *seq, param_t params, double *res);
void calc_profile_bug(seq_t *seq_holder, param_t params, double *profile);

#ifdef HAVE_GNUPLOT
/* display or produce the hydrophobic profile using gnuplot */
void gplot (double *plot, segment_t **segments, seq_t *seq, int i,
	    param_t params);
#endif /* HAVE_GNUPLOT */

#endif
