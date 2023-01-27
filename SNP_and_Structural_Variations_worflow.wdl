version 1.0

import "wdls/samtools_faidx.wdl" as samtools
import "wdls/pbmm2.wdl" as pbmm2
import "wdls/pbsv.wdl" as pbsv
import "wdls/deepvariant.wdl" as deepvariant
import "wdls/bgzip.wdl" as bgzip
import "wdls/tabix.wdl" as tabix
import "wdls/bcftools_merge.wdl" as bcftools_merge
import "wdls/unzip.wdl" as unzip
import "wdls/survivor.wdl" as survivor

workflow snp_and_snv_analysis {
	input {
		File reference_fasta # samtools pbmm2 pbsv deepvariant
		Array[File] fastqs # pbmm2 pbsv deepvariant
		String reference_name # pbmm2
		String sample_name = "some_name" #pbmm2 pbsv
		Array[String] regions_pbsv
		String regions_deepvariant
		File tr_bed #pbsv
		String model_type = "PACBIO" # deepvariant
		## String output_vcf_path # deepvariant
	}

	call samtools.call_samtools_faidx as samtools_faidx { 
		input:
		reference_fasta=reference_fasta
	}

	File reference_index = samtools_faidx.fai

	call pbmm2.run_pbmm2 as align {
		input:
		sample_name = sample_name,
		reference_name = reference_name,
		reference_fasta = reference_fasta,
		reference_index = reference_index,
		movies = fastqs
	}

	Array[File] bams = align.bams
	Array[File] bais = align.bais

	call pbsv.run_pbsv as call_svc {
		input:
		sample_name = sample_name,
		bams = bams,
		bais = bais,
		reference_name = reference_name,
		reference_fasta = reference_fasta,
		reference_index = reference_index,
		regions = regions_pbsv,
		tr_bed = tr_bed
	}

	Array[Int] indexes = range(length(bams))

	scatter (i in indexes) {
		call deepvariant.call_deepvariant as call_snps {
			input:
			model_type = model_type,
			reference_fasta = reference_fasta,
			reference_index = reference_index,
			reads_bam = bams[i],
			reads_index = bais[i],
			output_vcf_path =basename(sub(bams[i], ".bam", ".vcf.gz")),
			regions = regions_deepvariant,
			num_shards = 2 ## redo deepvariant.wdl
		}
	}
	
	Array[File] all_vcf_gzs = call_snps.snp_vcf

	Array[Int] indexes2 = range(length(all_vcf_gzs))

	scatter (i in indexes2) {
		call tabix.tabix as call_tabix {
			input:
			vcf_gz = all_vcf_gzs[i]
		}
	}

	Array[File] all_vcf_gz_tbis = call_tabix.vcf_gz_tbi

	call bcftools_merge.merge_vcfs as call_merge_vcf {
		input:
		vcfs = all_vcf_gzs,
		tbis = all_vcf_gz_tbis
	}

	call unzip.unzip as unzip_ {
		input:
		vcf_gz = call_svc.vcf
	}

	call survivor.survivor as survivor_ {
		input:
		vcf = unzip_.vcf
	}

	output {
		File vcf_svcs = unzip_.vcf
		File vcf_snps = call_merge_vcf.vcf
		Array[File] aligned_bams = align.bams
		Array[File] html_reports = call_snps.html_report
		File svc_report_txt = survivor_.report_txt
	}
}
