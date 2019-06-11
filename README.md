# Generate BED file describing probe footprints +/- 5kb with bedtools [1]

```bash
$ bedtools slop \
    -b 5000 \
    -g hs37d5.fasta.fai \
    -i IDT_PacBio_Duke_SNCA-capture.bed | \
  bedtools merge > target.bed
```

# Short variant detection and phasing workflow

1. SMRT Analysis 6.0.0 [2] Demultiplex Barcodes analysis
  ```
  PD-1: gcagtcgaacatgtagctgactcaggtcacCACATATCAGAGTGCG
  PD-2: gcagtcgaacatgtagctgactcaggtcacACACACAGACTGTGAG
  PD-3: gcagtcgaacatgtagctgactcaggtcacACACATCTCGTGAGAG
  PD-4: gcagtcgaacatgtagctgactcaggtcacCACGCACACACGCGCG
  N-1: gcagtcgaacatgtagctgactcaggtcacCACTCGACTCTCGCGT
  N-2: gcagtcgaacatgtagctgactcaggtcacCATATATATCAGCTGT
  N-3: gcagtcgaacatgtagctgactcaggtcacTCTGTATCTCTATGTG
  N-4: gcagtcgaacatgtagctgactcaggtcacACAGTCGAGCGCTGCG
  DLB-1: gcagtcgaacatgtagctgactcaggtcacACACACGCGAGACAGA
  DLB-2: gcagtcgaacatgtagctgactcaggtcacACGCGCTATCTCAGAG
  DLB-3: gcagtcgaacatgtagctgactcaggtcacCTATACGTATATCTAT
  DLB-4: gcagtcgaacatgtagctgactcaggtcacACACTAGATCGCGTGT
  ```
2. SMRT Analysis 6.0.0 Circular Consensus Sequence analysis
3. Merge consensus read sets from multiple runs.  Make symlinks with convenient names. [3]
  ```bash
  $ parallel 'dataset merge {}.merged.consensusreadset.xml {}.*.consensusreadset.xml' ::: {PD,N,DLB}-{1..4}
  ```
4.  Align merged consensus read sets to hs37d5. [4,5]
  ```bash
  $ export REF=hs37d5.fasta
  $ SAMPLE=PD-1 CCSSET=PD-1.merged.consensusreadset.xml \
    OUTSET=PD-1.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=PD-2 CCSSET=PD-2.merged.consensusreadset.xml \
    OUTSET=PD-2.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=PD-3 CCSSET=PD-3.merged.consensusreadset.xml \
    OUTSET=PD-3.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=PD-4 CCSSET=PD-4.merged.consensusreadset.xml \
    OUTSET=PD-4.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=N-1 CCSSET=N-1.merged.consensusreadset.xml \
    OUTSET=N-1.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=N-2 CCSSET=N-2.merged.consensusreadset.xml \
    OUTSET=N-2.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=N-3 CCSSET=N-3.merged.consensusreadset.xml \
    OUTSET=N-3.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=N-4 CCSSET=N-4.merged.consensusreadset.xml \
    OUTSET=N-4.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=DLB-1 CCSSET=DLB-1.merged.consensusreadset.xml \
    OUTSET=DLB-1.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=DLB-2 CCSSET=DLB-2.merged.consensusreadset.xml \
    OUTSET=DLB-2.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=DLB-3 CCSSET=DLB-3.merged.consensusreadset.xml \
    OUTSET=DLB-3.consensusalignmentset.xml qsub pbmm2.sh
  $ SAMPLE=DLB-4 CCSSET=DLB-4.merged.consensusreadset.xml \
    OUTSET=DLB-4.consensusalignmentset.xml qsub pbmm2.sh
    
  # after alignments are complete
  $ parallel 'samtools index {}' ::: {PD,N,DLB}-{1..4}.bam
  ```
5. Deduplicate aligned BAMs.
  ```bash
  $ parallel 'python pcr_duplicates.py \
                       --outBAM {.}.dedup.bam \
                       --wiggle 4 {} > {.}.dup && \
              samtools index {.}.dedup.bam' ::: *.bam
  ```
6. Produce a gVCF using GATK4 HaplotypeCaller. [6]
  ```
  $ export REF=hs37d5.fasta
  $ parallel 'BED=target.bed BAMLIST={}.dedup.bam SAMPLE={} PLOIDY=2 PCRINDELMODEL=HOSTILE MAPQ=50 \
      qsub gatk_hc_gvcf_bed.sh' ::: {PD,N,DLB}-{1..4}
  ```
7. Combine gVCFs.
  ```bash
  $ export REF=hs37d5.fasta
  $ source Broad.env
  $ $GATK CombineGVCFs \
      --reference $REF \
      --output cohort.g.vcf.gz \
      --variant DLB-1.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant DLB-2.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant DLB-3.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant DLB-4.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant N-1.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant N-2.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant N-3.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant N-4.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant PD-1.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant PD-2.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant PD-3.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz \
      --variant PD-4.2ploid.HOSTILE.HaplotypeCaller.g.vcf.gz
  ```
8. Produce a VCF from combined gVCFs.
  ```bash
  $ export REF=hs37d5.fasta
  $ BED=target.bed GVCF=cohort.g.vcf.gz DBSNP=00-common_all.vcf.gz qsub gatk_genotypeGVCFs_bed.sh
  ```
9. Filter the VCF.
  ```bash
  $ VCF=cohort.vcf.gz qsub filter_SNPs.sh
  
  # generate BED with all homopolymer stretches >= 7
  $ python repeat_bed.py --threads 12 --minhrun 7 hs37d5.fasta
  $ mv hs37d5.fasta.dinuc_repeats_ge0.hrun_ge7.bed hs37d5.hrun_ge7.bed
  $ export BED=hs37d5.hrun_ge7.bed
  
  $ VCF=cohort.vcf.gz qsub filter_INDELs.sh
  $ VCF=cohort.vcf.gz qsub filter_MIXED.sh
  ```
10. Merge filtered VCFs.
  ```bash
  $ $GATK MergeVcfs \
      --INPUT cohort.vcf.filtered_snps.vcf \
      --INPUT cohort.vcf.filtered_indels.vcf \
      --INPUT cohort.vcf.filtered_mixed.vcf \
      --OUTPUT cohort.filtered.vcf
  ```
11. Remove filtered sites from VCF.
  ```bash
  $ $GATK SelectVariants \
      --output cohort.filtered.selected.vcf \
      --variant cohort.filtered.vcf \
      --exclude-filtered \
      --exclude-non-variants \
      --remove-unused-alternates \
      --set-filtered-gt-to-nocall
  ```
12. Manually curate indels, mixed variants, and variants without annotated rsid.  Add rsid for variants equivalent to those already in dbSNP.  Remove any indels (or indel components of mixed variants) that do not phase with nearby SNPs in at least one sample.  See `manual_curation.xslx`.  Final product: `cohort.curated.vcf.gz`
13. Phase with Whatshap.  Generate statistics and haplotag BAMs. [7,8]
  ```bash
  $ export REF=hs37d5.fasta
  $ samtools merge -@ 11 cohort.bam {PD,N,DLB}-{1..4}.dedup.bam && samtools index -@ 11 cohort.bam
  $ whatshap phase \
      --indels \
      --output cohort.phased.vcf.gz \
      --reference $REF \
      cohort.curated.vcf.gz \
      cohort.bam
  $ for i in {PD,N,DLB}-{1..4}; do
      whatshap stats --sample ${i} \
        --gtf ${i}.phase_blocks.gtf \
        --tsv ${i}.phased.tsv \
        --block-list ${i}.phased.blocklist \
        cohort.phased.vcf.gz
    done
  $ for i in {PD,N,DLB}-{1..4}; do
      whatshap haplotag --output ${i}.haplotagged.bam \
        --reference $REF cohort.phased.vcf.gz \
        ${i}.dedup.bam
    done
  ```
  
# Haplotype determination for CT-rich region chr4:90742331-90742559 (hs37d5)

1. Extract subsequence from aligned, deduplicated reads.
   ```bash
   $ parallel 'extract_aligned_fasta.py {} {.}.fasta 4 90742331 90742559' ::: *.bam
   ```
2. Examine distribution of extracted subsequence sizes per sample.
   ```python
    import matplotlib.pyplot as plt
    from Bio import SeqIO
    
    samples = ['-'.join([x,str(y)]) for x in ['PD', 'N', 'DLB'] for y in range(1,5)]
    fig, axarr = plt.subplots(12, sharex=True)
    fig.set_size_inches(12, 24)
    for i, sample in enumerate(samples):
        data = []
        for record in SeqIO.parse(sample + '.rmdup.fasta', 'fasta'):
            data.append(len(record))
        axarr[i].hist(data, bins=range(220,360,2), label=sample)
        axarr[i].legend()
    axarr[-1].set_xlabel('size (bp)');
   ```
3. Based on separation of subsequence size distribution, divide haplotypes by size.
   ```python
    import matplotlib.pyplot as plt
    from Bio import SeqIO
    
    samples = ['-'.join([x,str(y)]) for x in ['PD', 'N', 'DLB'] for y in range(1,5)]
    breakpoints = {
        'PD-1': 300,
        'PD-2': 300, 
        'PD-3': 240, 
        'PD-4': 300, 
        'N-1': 300,
        'N-2': 300, 
        'N-3': 300, 
        'N-4': 300, 
        'DLB-1': 300, 
        'DLB-2': 241, 
        'DLB-3': 240, 
        'DLB-4': 300}
    for i, sample in enumerate(samples):
        lowrecords = []
        highrecords = []
        for record in SeqIO.parse(sample + '.rmdup.fasta', 'fasta'):
            if len(record) < breakpoints[sample]:
                lowrecords.append(record)
            else:
                highrecords.append(record)
            SeqIO.write(lowrecords, sample + '_lt' + str(breakpoints[sample]) + 'bp.fasta', 'fasta')
            SeqIO.write(highrecords, sample + '_ge' + str(breakpoints[sample]) + 'bp.fasta', 'fasta')
   ```
4. For each size-clustered fasta, cluster and generate consensus with MUSCLE[9].

[1] Bioinformatics, Volume 26, Issue 6, 15 March 2010, Pages 841–842, https://doi.org/10.1093/bioinformatics/btq033

[2] https://www.pacb.com/support/software-downloads/

[3] O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47.

[4] https://github.com/PacificBiosciences/pbmm2

[5] https://github.com/lh3/minimap2

[6] From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M, 2013 CURRENT PROTOCOLS IN BIOINFORMATICS 43:11.10.1-11.10.33

[7] Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]

[8] Marcel Martin, Murray Patterson, Shilpa Garg, Sarah O. Fischer, Nadia Pisanti, Gunnar W. Klau, Alexander Schoenhuth, Tobias Marschall. WhatsHap: fast and accurate read-based phasing, bioRxiv 085050, doi: 10.1101/085050

[9] Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput.Nucleic Acids Res. 32(5):1792-1797. doi:10.1093/nar/gkh340
