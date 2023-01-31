rule liftover_germline:
	input: "results/{run}/germline/vcf/{sample}.vcf.gz"
	output: 
		intermediate=temp("results/{run}/germline/vcf/{sample}.37.tmp.vcf.gz"),
		vcf="results/{run}/germline/vcf/{sample}.37.sorted.vcf.gz"
	params: 
		crossmap=config['tools']['crossmap'],
		chain=config['liftover']['chain'],
		reference=config['liftover']['reference_fasta'],
		bcftools=config['tools']['bcftools']
	shell: """
			{params.crossmap} vcf \
			{params.chain} {input} \
			{params.reference} {output.intermediate} \
			--compress --chromid s && \
			{params.bcftools} sort {output.intermediate} -Oz -o {output.vcf}
			"""

use rule liftover_germline as liftover_germline_joint with:
	input: 
		"results/{run}/germline/vcf/cohort.vcf.gz"
	output:
		intermediate="results/{run}/germline/vcf/cohort.37.tmp.vcf.gz",
		vcf="results/{run}/germline/vcf/cohort.37.sorted.vcf.gz"

use rule liftover_germline as liftover_somatic with:
	input: 
		'results/{run}/somatic/{patient}/mutect2.final.vcf.gz'
	output: 
		intermediate='results/{run}/somatic/{patient}/mutect2.37.tmp.vcf.gz',
		vcf='results/{run}/somatic/{patient}/mutect2.37.sorted.vcf.gz'

# rule liftover_sv:
# 	input: 'results/{run}/germline/sv/results/variants/diploidSV.inv_converted.vcf.gz'
# 	output: 'results/{run}/germline/sv/results/variants/sv.liftovered_37.unsorted.vcf.gz'
# 	params: 
# 		chain=config['chain'],
# 		reference=config['reference_37']
# 	shell: """
# 			CrossMap.py vcf \
# 			{params.chain} {input} \
# 			{params.reference} {output} \
# 			--compress"""

# rule preannotation_prep_vcf_sv:
# 	input: 'results/{run}/germline/sv/results/variants/sv.liftovered_37.unsorted.vcf.gz'
# 	output: 'results/{run}/germline/sv/results/variants/sv.liftovered_37.vcf.gz'
# 	shell: """
# 			for i in {{1..22}} X Y MT; \
# 				do \
# 					echo chr${{i}} $i | \
# 					sed 's/chrMT/chrM/'; \
# 				done | \
# 			bcftools annotate --rename-chrs - {input} -Ou | \
# 			bcftools sort -Oz -o {output}
# 			"""
