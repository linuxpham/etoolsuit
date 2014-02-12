Usage
=====
NAME
        Etool Suit Application for Biologist

VERSION
        1.0.10

SYNOPSIS
        etool   <command>
                + pre <reference fasta> <left fastQ genome> <right fastQ genome> <output directory> <thread number>: Prepare processing
                + bwa <command> [options]: BWA tool
                + sam <command> [options]: Sam tool
                + bed <command> [options]: Bed tool
                + bmv : Bam viewer
                + igv : IGV browser
                + sib <options>: Sibelia
                + csb <options>: C-Sibelia

EXAMPLES
        etool pre reference.fasta left.fastq.gz right.fastq.gz /tmp/test 4 (Prepare processing)
        etool bwa index reference.fasta (Index reference genome)
        etool bwa mem -a -t 10 reference.fasta left.fastq.gz right.fastq.gz > out.sam (FastQ alignment)
        etool sam view -bS out.sam > out.bam (Create BAM)
        etool sam sort out.bam out.sorted (Sort BAM)
        etool sam index out.sorted.bam (Index sorted BAM)
        etool bed bamtobed -i out.sorted.bam > out.bed (BAM to BED)
        etool bmv (BAM viewer)
        etool igv (IGV browser)
        etool sib (Sibelia) 
        etool csb (C-Sibelia)

DIAGNOSTICS
        The etool utility exits 0 on success, and >0 if an error occurs.
