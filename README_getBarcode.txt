README for code NGSMapping_getBarcode_main.cpp

this is the program to get unique of the barcodes/indexes in the read.
so far we intend to do on undetermined reads to characterize the distribution
of the barcode/indexes. We only implement to do this on read data, to get
barcodes from the sequence names, but not from index reads.
we will not demux, but simply do stats. To do demux, we need to run
NGSMapping_barcode_main.cpp to get partial match etc.

to call
ngsmapping_getBarcode -n -s file.fastq(.gz) -d

  #-n to do read from sequence name, -s is the data file
  #-d is the dualindex mode.
there are 3 outputs, unique barcode in fasta format (r1 and r2),
      stats/counts of each barcods
      table format summarizing together including above two.
