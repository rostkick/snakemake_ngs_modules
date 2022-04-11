rule liftover_germline:
	input: "results/{run}/germline/vcf/cohort.filtered.vcf.gz"
	output: "results/{run}/germline/vcf/cohort.filtered.liftovered_37.unsorted.vcf.gz"
	params: 
		chain=config['chain'],
		reference=config['reference_37']
	shell: """
			CrossMap.py vcf \
			{params.chain} {input} \
			{params.reference} {output} \
			--compress"""

rule preannotation_prep_vcf_germline:
	input: "results/{run}/germline/vcf/cohort.filtered.liftovered_37.unsorted.vcf.gz"
	output: "results/{run}/germline/vcf/cohort.filtered.liftovered_37.vcf.gz"
	shell: """
			for i in {{1..22}} X Y MT; \
				do \
					echo chr${{i}} $i | \
					sed 's/chrMT/chrM/'; \
				done | \
			bcftools annotate --rename-chrs - {input} -Ou | \
			bcftools sort -Oz -o {output}
			"""

rule liftover_somatic:
	input: 'results/{run}/somatic/{patient}/mutect2.filtered.pass.vcf.gz'
	output: 'results/{run}/somatic/{patient}/mutect2.liftovered_37.unsorted.vcf.gz'
	params: 
		chain=config['chain'],
		reference=config['reference_37']
	shell: """
			CrossMap.py vcf \
			{params.chain} {input} \
			{params.reference} {output} \
			--compress"""

rule preannotation_prep_vcf_somatic:
	input: 'results/{run}/somatic/{patient}/mutect2.liftovered_37.unsorted.vcf.gz'
	output: 'results/{run}/somatic/{patient}/mutect2.liftovered_37.vcf.gz'
	shell: """
			for i in {{1..22}} X Y MT; \
				do \
					echo chr${{i}} $i | \
					sed 's/chrMT/chrM/'; \
				done | \
			bcftools annotate --rename-chrs - {input} -Ou | \
			bcftools sort -Oz -o {output}
			"""


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
