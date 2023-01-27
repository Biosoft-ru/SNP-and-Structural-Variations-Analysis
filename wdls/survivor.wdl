version 1.0

task stats {
	input {
		File vcf
	}

	command <<<
		/SURVIVOR/Debug/SURVIVOR stats ~{vcf} -1 -1 -1 "svc_stats.txt"
	>>>

	output {
		File report_txt = "svc_stats.txt"
	}

	runtime {
		docker: "developmentontheedge/survivor:latest"
	}
}

workflow survivor {
	input {
		File vcf
	}

	call stats {
		input:
		vcf = vcf
	}

	output {
		File report_txt = stats.report_txt
	}
}
