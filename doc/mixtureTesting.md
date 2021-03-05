Mixture testing
================

The server can test whether particular sequences are likely to contain mixed base calls,
by computing the number of Ns, Ms (a term used to include all IUPAC mixed base designations), or both
in the variant bases in a collection of similar sequences,
relative to various null hypotheses.  Four hypotheses are tested by default; it is up to the user
to choose which is an appropriate test for their data set.

All the *p* values reported are derived from exact, one-sided binomial tests as implemented in python's scipy.stats.binom_test().

For *Mycobacterium tuberculosis* as mapped by our pipeline, testing using test1, 2 or 3 gives very similar answers.   Test 1 is recommended (i.e. using p_value1 as the mixture criterion).  
For *Salmonella enterica* as mapped by our pipeline, testing using test4 may be appropriate when applied to clusters,
because the proportion of Ns in the multisequence alignment of different SNV-based clusters differs markedly, presumably because
some sequences match the reference better than others.

The clustering algorithms which are SNP based can use estimates of whether the sequence is mixed.  
Which test is used has to be prespecified when the clustering algorithm is defined.  If you wish to compare the results of (say) 50 SNV clustering with
different p-value criteria, you can define multiple clustering algorithms using different criteria and compare them in parallel on the same trial dataset.  

The cutoff value also needs careful choice.  Bear in mind that a test is applied to each sequence in a cluster.
Thus, if there are 100k samples, up to 100k tests will be performed.  At present, no FDR adjustment is performed, although we could implement this.  
At present adjustment by Bonferroni is therefore recommended.  Thus, if you anticipate there will be 100k samples in the database, then
up to 100k tests could be performed.  A p-value of 1x10<sup>-7</sup> might be an appropriate value for the cutoff.


## Terminology 
*Variant bases* means the positions, containing A,C,T,G (i.e. confidently called bases) which vary between a set of similar sequences relative to  
The *alignment* means the collection of variant bases.  
The *whole population of sequences* means all sequences present in the server.  

## TEST 1:
This tests the hypothesis that the number of Ns in the *alignment*
is GREATER than those expected from the expected_N in the population of whole sequences.

Does so by comparing the observed number of Ns in the alignment (alignN),
given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
The expected number of Ns is estimated by
i) randomly sampling sample_size guids from those stored in the server and
observing the number of Ns per base across the genome.  The estimate_expected_N() function performs this.
ii) randomly sampling sample_size guids from those stored in the server and
observing the number of Ns per base across the relevant  genome.  The estimate_expected_N() function performs this.
  
This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
is a relatively unbiased estimate of the median number of Ns in the good quality samples.

If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
reported as None.

## TEST 2:
This tests the hypothesis that the number of Ns in the *alignment*
is GREATER than those expected from the expected_N in the population of whole sequences
*at the bases examined in the alignment*.
This might be relevant if these particular bases are generally hard to call.

Does so by comparing the observed number of Ns in the alignment (alignN),
given the alignment length  and an expectation of the proportion of bases which will be N.
The expected number of Ns is estimated by randomly sampling *sample_size* (a configurable parameter) samples from those stored in the server and
observing the number of Ns per base at the relevant sites.  The estimate_expected_N_sites() function performs this.

This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
is a relatively unbiased estimate of the median number of Ns in the good quality samples.

If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
reported as None.
          
## TEST 3:  
tests whether the proportion of Ns in the alignment is greater
than in the bases not in the alignment, for this sequence.

## TEST 4:  
tests whether the proportion of Ns in the alignment  for this sequence
is greater than the median proportion of Ns in the alignment for all other sequences.
This is a sensible option if the test is performed per-cluster, and the proportion of Ns in the
aligned sequences differs markedly by cluster.

This test is computed if there are two or more samples in the cluster.

## Note
This testing is performed by the *seqComparer._msa()* module.