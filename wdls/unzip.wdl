version 1.0

task unzip_ {
	input {
		File vcf_gz
	}

	command <<<
		bgzip -dc ~{vcf_gz} > pbsv_svc.vcf
	>>>
	
	output {
		File vcf = "pbsv_svc.vcf"
	}

	runtime {
		docker: "developmentontheedge/vcftools:v0.0.2"
	}
}

workflow unzip{
	input {
		File vcf_gz
	}

	call unzip_ {
		input:
		vcf_gz = vcf_gz
	}

	output {
		File vcf = unzip_.vcf
	}
}
