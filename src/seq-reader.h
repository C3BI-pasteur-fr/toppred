/*   File:                /home/edeveaud/Work/dnatool/src/readseq.h
 *   Author:              Eric Deveaud edeveaud@pasteur.fr
 */


#ifndef __SEQ_READER_H_
#define __SEQ_READER_H_

/* sequence size switch between memory and file processing */
#define MAX_MEM_SIZE 1000000

/* default allowed characters in sequence */
#define DEFAULT_ALPHABET "ARNDCQEGHILKMFPSTWYV"

#define BUFFSIZE 100
#define SEQ_NONE -1
#define HEADER 0
#define SEQ 1
#define DUMP 2

#define NO_SEQ 0
#define SEQ_OK 1
#define SEQ_OVERMEM 2



typedef struct seq_S {
  char *id;
  char *comment;
  char *seq;
  size_t size;
} seq_t;



/* retreive sequence in fasta format from file descriptor */
int read_seq(FILE *IN, seq_t  *seq_holder, char *alphabet) ;

/* generate alphabet table used to check the sequences */
char *alphabet_maker (char *alphabet) ;

/* free the sequence holder */
void free_seq(seq_t *seq) ;

#endif /* __SEQ_READER_ */
