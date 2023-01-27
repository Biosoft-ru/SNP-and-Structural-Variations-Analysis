version 1.0
task samtools_faidx {
        input {
                File reference_fasta
        }

        command <<<
                samtools faidx ~{reference_fasta}
                mv ~{reference_fasta}.fai $(basename ~{reference_fasta}.fai)
        >>>

        output {
                File reference_index = "${reference_fasta}.fai"
        }

        runtime {
                docker: "developmentontheedge/samtools:latest"
        }
}

workflow call_samtools_faidx {
        input {
                File reference_fasta
        }

        call samtools_faidx {
                input:
                reference_fasta=reference_fasta,
        }

        output {
                File fai=samtools_faidx.reference_index
        }

        meta {
                description: "##making index of reference fasta"
        }
}
