---
title: Sanger-tol/Variantcalling Workflow
---

%%{ init: {
'gitGraph': {'mainBranchName': 'BAM/CRAM'
},
'themeVariables': {
'commitLabelFontSize': '18px'
}
}
}%%
gitGraph TB:
commit id: "START"
branch Fasta order: 4
commit id: "SAMTOOLS_FAIDX"
checkout BAM/CRAM
commit id: "SAMPLESHEET_CHECK"
branch AlignedReads order: 3
branch UnAlignedReads order: 2
commit id: "SAMTOOLS_COLLATE"
commit id: "SAMTOOLS_FASTA"
commit id: "BLAST_BLASTN"
commit id: "PACBIO_FILTER"
commit id: "PACBIO_SAMTOOLS_FILTER"
commit id: "SAMTOOLS_FASTQ"
commit id: "MINIMAP2_ALIGN"
commit id: "SAMTOOLS_MERGE_BY_SAMPLE 1" type: HIGHLIGHT
commit id: "SAMTOOLS_STATS" type: HIGHLIGHT
commit id: "SAMTOOLS_FLAGSTAT" type: HIGHLIGHT
commit id: "SAMTOOLS_IDXSTATS" type: HIGHLIGHT
checkout BAM/CRAM
merge UnAlignedReads
checkout Fasta
branch SplitFasta order: 5
commit id: "FASTA_SPLIT"
branch DeepVariant order: 6
checkout AlignedReads
merge Fasta
commit id: "SAMTOOLS_SORT"
commit id: "SAMTOOLS_MERGE_BY_SAMPLE 2"
checkout BAM/CRAM
merge AlignedReads
commit id: "SAMTOOLS_FILTER"
checkout DeepVariant
merge BAM/CRAM
commit id: "DEEPVARIANT"
commit id: "BCFTOOLS_CONCAT_VCF" type: HIGHLIGHT
branch "VCFtools" order: 7
commit id: "VCFTOOLS_SITE_PI" type: HIGHLIGHT
commit id: "VCFTOOLS_HET" type: HIGHLIGHT
checkout DeepVariant
commit id: "BCFTOOLS_CONCAT_GVCF" type: HIGHLIGHT
checkout BAM/CRAM
merge VCFtools
merge DeepVariant
commit id: "DUMPSOFTWAREVERSIONS" type: HIGHLIGHT
