# Genome:

I acquired the human genome gene annotation through [CoGe](https://genomevolution.org/coge/GenomeInfo.pl?gid=25747). And I acquired a TE annotation through the [UCSC Genome Browser](https://genome-euro.ucsc.edu/cgi-bin/hgTables) using the options groups: repeats, track: RepeatMasker, and output format: all fields. I had to do some reformatting for this TE annotation to get it to be compliant with TE Density and those scripts are present in: xyz (**CITE**))

Reformatting of human gene annotation file, *gencode.v37.annotation.gff3*, use commands *python examples/Human/import_scripts/import_human_gene_anno.py* on the gene file.
