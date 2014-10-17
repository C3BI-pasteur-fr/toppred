TESTS = toppred.test hydro.test \
        detect_segments.test seqlen.test construct_topos.test \
	naming.test math.test
#	float.test naming.test

SEQS = seq-test \
       first_seg.fasta last_seg.fasta no_seg.fasta \
       min_seqlen.fasta \
       too_many_put.fasta more_calc.fasta only_put.fasta \
       seq_anonymous.fasta seq_zero_div.fasta
#       seq_float.fasta seq_float2.fasta seq_float3.fasta


OUT =  first_seg.out last_seg.out no_seg.out \
       min_seqlen.out \
       too_many_put.out more_calc.out only_put.out

ERR =  first_seg.err last_seg.err no_seg.err \
       min_seqlen.err seq_anonymous.err \
       too_many_put.err more_calc.err only_put.err \
       no_error.err
       

EXTRA_DIST = $(TESTS) $(SEQS) $(OUT) $(ERR)

CLEANFILES = *.hydro *_tmp.out *_tmp.err *.png
