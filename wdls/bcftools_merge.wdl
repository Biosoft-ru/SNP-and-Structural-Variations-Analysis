version 1.0

task merge {
	input {
		Array[File] vcfs
		Array[File] tbis
		String? collapse="snps|indels"
		String? regions
	}
	
	Int vcfs_length = length(vcfs)

	command <<<
		if [[~{vcfs_length} -gt 1]]
		then
		bcftools merge --force-samples ~{sep=" " vcfs} > deepvariant_snp.vcf
		else
		cp ~{sep="" vcfs} deepvariant_snp.vcf.gz
		bgzip -d deepvariant_snp.vcf.gz
		fi
	>>>

	output {
		File merged_vcf = "deepvariant_snp.vcf"
	}

	runtime {
		docker: "developmentontheedge/bcftools:1.14"
	}
}

workflow merge_vcfs {
	input {
		Array[File] vcfs
		Array[File] tbis
		String? collapse
		String? regions
	}

	call merge {
		input:
		vcfs=vcfs,
		tbis=tbis,
		collapse=collapse,
		regions=regions
	}
	
	output {
		File vcf = merge.merged_vcf
	}
}