version 1.0

task tabix_ {
	input {
		File vcf_gz
		String? preset = "vcf"
	}

	command <<<
		tabix -p ~{preset} ~{vcf_gz}
		mv ~{vcf_gz}.tbi $(basename ~{vcf_gz}.tbi)
	>>>

	output {
		File vcf_gz_tbi = "${vcf_gz}.tbi"
	}

	runtime {
		docker: "developmentontheedge/vcftools:v0.0.2"
	}
}

workflow tabix {
	input {
		File vcf_gz
		String? preset
	}

	call tabix_ {
		input:
		vcf_gz=vcf_gz,
		preset=preset
	}

	output {
		File vcf_gz_tbi = tabix_.vcf_gz_tbi
	}
}
