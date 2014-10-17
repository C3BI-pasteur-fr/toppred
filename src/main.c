/*   File:                /home/edeveaud/Work/toppred/src/main.c
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <libgen.h>
#include <math.h>

#include <sys/stat.h>
#include <sys/types.h>

#ifdef HAVE_LIBGD
#include "graph.h"
#endif

/* #include "params.h" */
#include "main.h"
#include "error.h"
#include "usage.h"
#include "seq-reader.h"
#include "profile.h"
#include "loop.h"
#include "topology.h"
#include "topoprint.h"
#include "output.h"
#include "mloutput.h"

static char *prog;

static void process_seq(FILE *IN, param_t params);

int main(int argc, char **argv) {

  /* variables and initialisation */
  param_t params;
  int i, uplot;
  char  *out_file, *p;
  FILE *IN;
  char *data_dir;
  char *buf;
  size_t len;

  /* default values */
  out_file = NULL;
  params.web = 0;
  params.n = 11;
  params.q = 5;
  params.p_cut = 0.6; /* 0.5; */
  params.c_cut = 1.0;
  params.n_topos = 16;
  params.seg_len = 60; /* 70; */
  params.OUT = stdout;
  params.gplot = TRUE;
  params.hydro_file = TRUE;
  params.plot_format = "x11";
  params.plot_pause = "";
  params.output = NEW;
  params.topo_format = "png";
  params.tmspacer = 2;
  params.kingdom = PROKARYOTE;
  uplot = FALSE;

  /* check to see if a TOPPREDDATA environment variable is available */
  /* yes overwritte the DATADIR default, and use it*/
  if ((data_dir = getenv("TOPPREDDATA")) == NULL){
    data_dir = DATADIR;
  }

  params.data_file = "GES-scale";
  params.ce_file = "CYTEXT-scale";

  /* get progname */
  prog = argv[0];
  if((p = strrchr(prog, '/')) != NULL)
    prog = ++p;


  /* check syntax option on command line */
  i = 0;
  while((i = getopt(argc, argv, "c:d:eg:hH:n:N:o:O:p:q:s:t:vy")) != -1) {
    switch(i) {

    case 'c':
      params.c_cut = atof(optarg);
      break;

    case 'd':
      params.tmspacer = atoi(optarg);
      break;

    case 'e':
      params.kingdom = EUKARYOTE;
      break;

    case 'g':
#ifdef HAVE_GNUPLOT
      params.gplot = TRUE;
      params.plot_format = optarg;
      uplot = TRUE;
#else
      if (strcmp(optarg, "none") != 0) {
	(void)fprintf(stderr, "Missing gnuplot support: -g option unavailable");
	return EXIT_FAILURE;
      }
      params.gplot = FALSE;
      params.plot_format = "none";
#endif /* HAVE_GNUPLOT */
      break;

    case 'h':
      usage(prog);
      return EXIT_SUCCESS;

    case 'H':
      params.data_file = optarg;
      break;

    case 'n':
      params.n = atoi(optarg);
      break;

    case 'N':
      params.n_topos = atoi(optarg);
      break;

    case 'o':
      out_file = optarg;
      break;

    case 'O':
      params.output = check_output_format(prog, optarg);
      break;

    case 'p':
      params.p_cut = atof(optarg);
      break;

    case 'q':
      params.q = atoi(optarg);
      break;

    case 's':
      params.seg_len = atoi(optarg);
      break;

    case 't':
#ifdef HAVE_LIBGD
      params.topo_format = optarg;
#else
      if (strcmp(optarg, "none") != 0) {
	(void)fprintf(stderr, "Missing libgd support: -t option unavailable");
	return EXIT_FAILURE;
      }
      params.topo_format = "none";
#endif
      break;

    case 'v':
      (void)fprintf(stdout, "%s (%s %s)\n", prog, PACKAGE, VERSION);
      return EXIT_SUCCESS;

    case 'y':
      params.hydro_file = FALSE;
      break;

    default :
      usage(prog);
      return EXIT_FAILURE;
    }
  }

  /* change default values for different output formats */
  if (params.output == HTML && !uplot) {
    params.plot_format = PNG;
  }

  /* checking validity of arguments */
  if (params.n < 0)
    error_fatal(prog, "n should be a positive number.");
  if (params.q < 0)
    error_fatal(prog, "q should be a positive number.");
  if(params.n + params.q == 0)
    error_fatal(prog, "Humm, n or q can be null, but not simultaneously.");
  if(params.c_cut < params.p_cut)
    error_fatal(prog, "certain cutoff should be greater than putative cutoff.");
  if(params.seg_len <= 0)
    error_fatal(prog, "critical segment length s should be a positive integer.");
  if(params.tmspacer <=0)
    error_fatal(prog, "critical tm distance should be a positive integer.");
  if(params.n_topos <= 0)
    error_fatal(prog, "N should be a positive number.");

  /* output file */
  if(out_file != NULL && (params.OUT = fopen(out_file, "w")) == NULL)
    error_fatal(out_file, NULL);

  /* get the directory holder for files*/
  params.out_dir = strdup(dirname(out_file));

  /* check the plot format */
  check_plot_format(prog, &params);

  /* check the topo format */
  check_topo_format(prog, &params);

  /* check if there's some files to deal with */
  if (optind >= argc) {
    usage(prog);
    return EXIT_FAILURE;
  }

  /* load the hphobes values once */
  len = strlen(data_dir) + 1 + strlen(params.data_file) + 1;
  if ((buf = (char *)malloc(len)) == NULL)
    error_fatal("memory", NULL);
  (void)snprintf(buf, len, "%s/%s", data_dir, params.data_file);
  read_Hphobes_datas(buf, params.Hdatas);
  free(buf);

  /* load cyt-ext values from cefile */
  len = strlen(data_dir) + 1 + strlen(params.ce_file) + 1;
  if ((buf = (char *)malloc(len)) == NULL)
    error_fatal("memory", NULL);
  (void)snprintf(buf, len, "%s/%s", data_dir, params.ce_file);
  (void) read_cytext_datas (buf, &(params.scales));
  free(buf);

  /* print parameter values */

  switch (params.output) {
  case HTML:
    break;
  case OLD:
    old_print_firstparameters(&params);
    break;
  default: /* NEW */
    new_print_parameters(&params);
  }

  IN = stdin;
  /* process all the input files */
  for (i = optind; i <argc; i++) {
    if (*argv[i] != '-' && (IN = fopen(argv[i], "r")) == NULL)
      error_fatal(argv[i], NULL);

    if (params.output == OLD) {
      (void) fprintf (params.OUT, "Using sequence file: %s\n\n", argv[i]);
    }

    process_seq(IN, params);

    if(fclose(IN) !=0)
      error_fatal(argv[i], NULL);

  }

  if(out_file != NULL && (fclose(params.OUT) == EOF)){
    error_fatal(out_file, NULL); }

  free(params.out_dir);

  return 0;
}


static void process_seq(FILE *IN, param_t params) {

  int i, s;
  int nb_plot;
  int ntopos, maxtopos;
  size_t len;
  double *Hprofile;
  seq_t seq;
  segment_t *segments;
  loop_t *KSloop;
  seg_t *KSseg;
  elem_t KSelem;
  topo_t *KStopo;
  int nel2prn;
  topoprint_t topo2prn;
  FILE *OUT;
  int skip;
  char *alphabet;
  char *err_msg;

  OUT = params.OUT;
  Hprofile = NULL;
  segments = NULL;


  alphabet = alphabet_maker(DEFAULT_ALPHABET);
  while(read_seq(IN, &seq, alphabet) !=0) {

/**** verify sequence ****/
    skip = 0;

    if(seq.size > MAXSEQLEN) {
      err_msg = "sequence too long -- skipped";
      skip = 1;
    }
    if(seq.size < ((2 * params.q) + params.n)) {
      err_msg = "sequence too short -- skipped";
      skip = 1;
    }
    if(seq.size == 0) {
      err_msg = "empty sequence -- skipped";
      skip = 1;
    }

    if(skip) {
      error_warn(seq.id, err_msg);
      free_seq(&seq);
      continue;
    }

/**** print sequence ****/
    switch (params.output) {
    case HTML:
      OUT = params.OUT = init_html (seq.id, params.out_dir);
      html_header (OUT, seq.id);
      html_parameters (OUT, &params);
      /* sequence is printed below the plot */
      break;
    case OLD:
    default:
      (void) fprintf(OUT, "\nSequence : %s  (%d res)\n", seq.id, seq.size);
      print_sequence(OUT, seq.seq);
      (void) fprintf(OUT, "\n");
    }

    if (params.output == OLD) { old_print_secondparameters(&params); }

/**** generation profile ****/

    /* allocate memory of the hydrophobic profile */
    nb_plot = seq.size - ((2 * params.q) + params.n) + 1;
    len = nb_plot * sizeof(double);
    if ((Hprofile = (double *)malloc(len)) == NULL) {
      error_fatal("memory", NULL);
    }

    /* calc Hval for all window positions */
    calc_profile(&seq, params, Hprofile);

    /*****  find transmembran segments  *****/
    s = get_segments(Hprofile, &segments, nb_plot , params);

/*****  produce profile plot with transmembran segment indication  *****/

    /* we produce the hydrophobic profile datas */
    if (params.hydro_file == TRUE)
      plot_values(Hprofile, &seq, params);

    switch (params.output) {
    case HTML:
      html_plot(OUT, seq.id, &params);
      html_sequence (OUT, &seq);
      break;
    /* case OLD: */
    /* default: */ /* NEW */
    }

#ifdef HAVE_GNUPLOT
    if(!params.plot_outfile && strlen(params.plot_pause) != 0) {
      (void)fprintf(stdout, "hit return to continue\n");
    }
    /* produce de gnuplot image */
    if(params.gplot){
      gplot(Hprofile, &segments, &seq, s, params);
    }
#endif

/**** transform segments structure to segment-loop structure   ****/

    if (s != 0) {
      (void)calc_loop(&seq, &segments, &KSloop,  &KSseg, params, s);
    }

    /* as profile and segments-storage struct is no longer needed purge it */
    free(Hprofile);
    free(segments);

/**** print transmembran segment summary ****/

    switch (params.output) {
    case OLD:
      break;
    case HTML:
      (void)fprintf (OUT, "<H4><CENTER>Transmembran segments</CENTER></H4>\n");
      start_phrase();
    default:
      (void) fprintf(OUT, "Found: %d segments\n\n", s);
      if (params.output == HTML) { end_phrase(); start_phrase(); }
      if (s) { print_tmsummary(OUT, s, KSseg); }
    }

    if (params.output == HTML) { end_phrase(); }

/**** construct topologies ****/

    if (s != 0) {

      KSelem.nsegs = s;
      KSelem.nputatives = 0;
      for (i=0; i< s; i++)
	if(KSseg[i].kind == 0) KSelem.nputatives++;
      KSelem.segs = KSseg;
      KSelem.loops = KSloop;

      if (KSelem.nputatives > MAXPUTATIVES_CALC) {
	error_warn(seq.id,
		   "too many putative segments to calculate best topologies");
      }
      else {

	maxtopos = ntopos = (int) pow(2.0, (double) KSelem.nputatives);

/**** print topology summary ****/

	if (params.output == HTML) {
	  (void) fprintf (OUT, "<H4><CENTER>Topologies</CENTER></H4>\n");
	  start_phrase();
	}


	(void)fprintf(OUT, "\nTotal of %d structures are to be tested\n\n",
		      maxtopos);
	if (params.output == OLD) {
	  (void) fprintf (OUT, "\n");
	  old_print_tmsummary(OUT, s, KSseg, seq.seq);
	}

	if (params.output == HTML) { end_phrase(); start_phrase(); }

	if (ntopos > params.n_topos) {
	  error_warn(seq.id, "more topologies than printed");
	  ntopos = params.n_topos;
	}

/**** calculate topologies and stock the ntopos best ones ****/

	if ((KStopo = (topo_t *) malloc (ntopos*sizeof(topo_t))) == NULL)
	  error_fatal("memory", NULL);

	ntopos = tp_calc(KStopo, &KSelem, &seq, &params);
	/* sort topologies by highest Arg+Lys bias */
	qsort((void *)KStopo, (size_t)ntopos, sizeof(topo_t), tp_compare);

/**** print the best topologies ****/

	/* allocate memory for max number of structure elements to print */
	nel2prn =  2*KSelem.nsegs + 1;
	if ((topo2prn.elps = (elprint_t *) malloc (nel2prn*sizeof(elprint_t))) == NULL)
	  error_fatal("memory", NULL);
	topo2prn.image = NULL;

	for (i=0; i<ntopos; i++) {
	  /* construct all topology information */
	  topo2prn.nr = i+1;
	  topo2prn.putatives = KStopo[i].putatives;
	  topo2prn.kr = KStopo[i].kr;
	  nel2prn = tp_decode (&topo2prn, &KSelem, &seq, &params);

/**** print topology image ****/
#ifdef HAVE_LIBGD
	  if(strcmp(params.topo_format, "none") != 0) {
	    int image_len = strlen(seq.id) + strlen(params.topo_format) + 6 + 3;
	    if ((topo2prn.image = (char *) realloc (topo2prn.image, image_len*sizeof(char))) == NULL)
	      error_fatal("memory", NULL);
	    (void)sprintf(topo2prn.image, "%s-%d.%s",
			  seq.id, topo2prn.nr, params.topo_format);

	    /* graphic print go here */
	    (void)topo_graph_print(&topo2prn, nel2prn, &params, &seq);
	  }
#endif

/**** printing nongraphic output ****/

	  switch (params.output) {
	  case HTML:
	    html_topology (&topo2prn, nel2prn, &params);
	    break;
	  case OLD:
	    (void)fprintf(OUT, "\n-----------------------------------------------------------------------\n");
	    tp_toppred_fprintf (&topo2prn, nel2prn, &params);
	    break;
	  default: /* NEW */
	    tp_new_print (&topo2prn, nel2prn, &params);
	  }

	} /* end for */

/**** nettoyage ****/

	/* free element structure for printing */
	if (topo2prn.image != NULL) {
	  free(topo2prn.image);
	}
	free(topo2prn.elps);
	free(KStopo);

      } /* end else of if (KSelem.nputatives > MAXPUTATIVES_CALC) */

      /* free KS structs even if no topologie calculated */
      free(KSloop);
      free(KSseg);

    } /* end if (s != 0) */

    /* freeing the seq-storage struct */
    free_seq(&seq);

    /* print sequence end informations and close html output file */

    switch(params.output) {
    case OLD:
      (void)fprintf(OUT, "\n-----------------------------------------------------------------------\n");
      break;
    case HTML:
      (void)fprintf(OUT, "\n");
      (void)fprintf(OUT, "</BODY>\n");
      (void)fprintf(OUT, "</HTML>\n");
      break;
      /* default : */
    }

    if (params.output == HTML) {
      if (fclose(params.OUT) == EOF) {
	error_fatal(prog, NULL);
      }
      OUT=stdout;
    }

  }/* end while j=read_seq */

  free(alphabet);

  return;
}


