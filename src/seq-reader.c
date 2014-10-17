/*   File:                /home/edeveaud/Work/dnatool/src/readseq.c
 *   Author:              Eric Deveaud <edeveaud@pasteur.fr>
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#include <string.h>
#endif

#include <ctype.h>

#include "seq-reader.h"
#include "error.h"



static void process_header(seq_t *seq, char *header) ;
static int clean_buff(char **buffer, char *alphabet) ;


int read_seq (FILE *IN, seq_t *seq_holder, char *alphabet) {

  char *BUFF, *stock, *p, *q;
  int i, state, l;
  size_t buffsize, bufflen, seq_len;

  /* check if there is something to read */
  if ((i = fgetc(IN)) == EOF && feof(IN) != 0) {
    return NO_SEQ; }

  ungetc(i, IN);

  state = SEQ_NONE;
  seq_len = bufflen = 0;
  buffsize = BUFFSIZE;
  seq_holder->seq = NULL;

  if ((BUFF = (char *) malloc(sizeof(char) * buffsize+1 )) == NULL)
    error_fatal ("memory1" , NULL);
  if ((stock = (char *) malloc(sizeof(char) * buffsize+1 )) == NULL)
    error_fatal ("memory2" , NULL);

  *stock = '\0';
  q = stock;

  while ((i = fgetc(IN)) != EOF) {
    /* skip empty lines */
    if (isspace(i)) {continue;}
    
    if ( ungetc(i, IN) == EOF ) error_fatal ("ungetc", NULL);

    /* end entry or start entry */
    if ( i == '>' ) {
      if (state == SEQ_NONE )  state = HEADER;
      if (state == SEQ) break ;
    }

    if (fgets(BUFF, BUFFSIZE + 1, IN) == NULL)  break;

    /* header processing */
    if ( state == HEADER ) {

      /* memcopy with '/0' inclusion */
      memcpy(q, BUFF, BUFFSIZE+1);
      /* is header line completly read */
      if (strrchr(BUFF, '\n') == NULL) {
	if((stock = realloc(stock, buffsize + BUFFSIZE + 1)) == NULL)
	  error_fatal("realloc", NULL);
	q = stock + buffsize ;
	buffsize += BUFFSIZE ;
	continue;
      }
      process_header(seq_holder, stock);
      state = SEQ;
      buffsize = BUFFSIZE;
      p = stock;
      *p ='\0';
      continue;
    }

    if ( state == SEQ || state == DUMP ){
      /* we clean the buffer before any further use */
      if ((l = clean_buff(&BUFF,  alphabet)) == -1)
	error_fatal(seq_holder->id, "sequence contain spurious characters");
    }

    if (state == SEQ ) {
      bufflen = l;
      /* is stock buffer long enough */
      if ( seq_len + bufflen >= buffsize ) {
	buffsize += BUFFSIZE;
	if ((stock = realloc(stock, sizeof(char) * buffsize+1)) == NULL)
	  error_fatal("Memory3", "Reallocating seq");
      }
      q = stock + seq_len;
      strncpy (q, BUFF, bufflen);
      seq_len += bufflen;
    }

    if ( state == SEQ_NONE ) {
      error_fatal("Sequence", "is NOT fasta formated");
    }
  }


  free(BUFF);
  if ( state == SEQ ) {
    *(stock+seq_len) = '\0';
    seq_holder->seq = stock;
  }
  seq_holder->size = seq_len;

  return state;

}

char *alphabet_maker (char *alphabet) {

  char  *res, *p;
  int j;

  if ((res = calloc(255, sizeof(char))) == NULL) {
    error_fatal("memory", NULL);
  }

  p = alphabet;
  while(*p) {
    j = (int)*p;
    if ( islower(j) ){
      res[j] = toupper((int)*p);
      res[j-32] = toupper((int)*p);
    }
    else{
      res[j] = *p;
      res[j+32] = *p;
    }
    p++;
  }

  return res;
}


static void process_header(seq_t *seq, char *header) {
  char *id, *comm, *p, *q;
  int i, n, l;

  n = 0;
  p = header;

  /* skip > */
  p++;

  while(*p && isascii((int)*p) && !isspace((int)*p)) {
    p++; n++;
  }

  /* check if sequence is named or not */
  if (n == 1 && strchr( ".,;", (int)header[1])) n = 0;
  l = ( n == 0) ? 9 : n;

  if((id = malloc(sizeof(char)*(l+1))) == NULL){
    error_fatal("memory", NULL);
  }
  p = header;
  p++;
  
  /* check if sequence is named or not */
  if ( n == 0 ){
    error_warn("anonymous sequence",
	       "name will be forced to \"anonymous\"");
    snprintf(id,10, "anonymous");
  }
  else {
    q = id;
    for (i=0; i<n; i++) {
      /* replace all non alnum character by '_' */
      if (isalnum((int)*p))  *q++ = *p++;
      else { *q++ = '_'; p++; }
    }
    *q = '\0';
  }


  /* skip spaces */
  while(*p && (isspace((int)*p))) { p++; n++; }

  /* get comment */

  if((comm = malloc(sizeof(char) * (strlen(header) - n))) == NULL){
    error_fatal("memory", "COMLEN");
  }
  q = comm;
  while(*p && *p != '\n') { *q++ = *p++; }
  *q = '\0';

  seq->id = id;
  seq->comment = comm;
}



static signed int clean_buff(char **buffer, char *alphabet) {
  char *p, *c;
  int bufflen;

  p = *buffer;
  c = *buffer;
  bufflen = 0;
  while (*p && *p != '\n') {
    if (isspace ((int)*p)) { p++; continue; }
    if( alphabet[(int)*p] != '\0') {
      *c++ = alphabet[(int)*p];
      p++;
      bufflen++;
    }
    else {
      if (isalpha((int)*p) || *p == '*') {
	char err_aa[2];
	err_aa[0] = *p;
	err_aa[1] = '\0';
	error_warn(err_aa, "not a valid aa, skipped");
	p++;
      }
      else {
	return -1;
      }
    }
  }
  *c = '\0';

  return bufflen;
}

void free_seq(seq_t *seq) {
  free(seq->id);
  free(seq->comment);
  if ( seq->seq != NULL ) free(seq->seq);
}
