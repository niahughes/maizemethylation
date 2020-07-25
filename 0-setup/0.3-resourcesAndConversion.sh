###################################################################
### RESOURCE DOWNLOADS AND CONVERSIONS 
###################################################################

# downloading conversion utilities; saved in /bin
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed

# download gene annotation file
wget ftp://ftp.ensemblgenomes.org/pub/release-44/plants/gff3/zea_mays/Zea_mays.B73_RefGen_v4.44.gff3.gz
gunzip Zea_mays.B73_RefGen_v4.44.gff3.gz

# download TE annotations
wget https://github.com/mcstitzer/maize_TEs/raw/master/B73.structuralTEv2.fulllength.2018-09-19.gff3.gz # most recent
gunzip B73.structuralTEv2.fulllength.2018-09-19.gff3.gz

wget ftp://ftp.gramene.org/pub/gramene/CURRENT_RELEASE/gff3/zea_mays/repeat_annotation/B73v4.TE.filtered.gff3.gz # legacy file with MITE annotations
gunzip B73v4.TE.filtered.gff3.gz
mv B73v4.TE.filtered.gff3 legacy.gff3

# convert gene gff3 to BED
gff3ToGenePred Zea_mays.B73_RefGen_v4.44.gff3 intermediateGene.genePred
genePredToBed intermediateGene.genePred geneAnnotationb73.bed

# download genomic FASTA files
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-36/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.chromosome.Pt.fa.gz # chloroplast
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz # whole genome