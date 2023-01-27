rule liftover_germline:
	input: "results/{run}/germline/vcf/{sample}.vcf.gz"
	output: "results/{run}/germline/vcf/{sample}.37.unsorted.vcf.gz"
	params: 
		chain=config['liftover']['chain'],
		reference=config['liftover']['reference_fasta']
	shell: """
			CrossMap.py vcf \
			{params.chain} {input} \
			{params.reference} {output} \
			--compress"""

rule preannotation_prep_vcf_germline:
	input: "results/{run}/germline/vcf/{sample}.37.unsorted.vcf.gz"
	output: "results/{run}/germline/vcf/{sample}.37.sorted.vcf.gz"
	shell: """
			for i in {{1..22}} X Y MT; \
				do \
					echo chr${{i}} $i | \
					sed 's/chrMT/chrM/'; \
				done | \
			bcftools annotate --rename-chrs - {input} -Ou | \
			bcftools sort -Oz -o {output}
			"""

use rule liftover_germline as liftover_germline_joint with:
	input: 
		"results/{run}/germline/vcf/cohort.vcf.gz"
	output: 
		"results/{run}/germline/vcf/cohort.37.unsorted.vcf.gz"
	# params: 
	# 	chain=config['liftover']['chain'],
	# 	reference=config['liftover']['reference_fasta']
	# shell: """
	# 		CrossMap.py vcf \
	# 		{params.chain} {input} \
	# 		{params.reference} {output} \
	# 		--compress"""

use rule preannotation_prep_vcf_germline as preannotation_prep_vcf_germline_join with:
	input: 
		"results/{run}/germline/vcf/cohort.37.unsorted.vcf.gz"
	output: 
		"results/{run}/germline/vcf/cohort.37.sorted.vcf.gz"
	# shell: """
	# 		for i in {{1..22}} X Y MT; \
	# 			do \
	# 				echo chr${{i}} $i | \
	# 				sed 's/chrMT/chrM/'; \
	# 			done | \
	# 		bcftools annotate --rename-chrs - {input} -Ou | \
	# 		bcftools sort -Oz -o {output}
	# 		"""

use rule liftover_germline as liftover_somatic with:
	input: 
		'results/{run}/somatic/{patient}/mutect2.final.vcf.gz'
	output: 
		'results/{run}/somatic/{patient}/mutect2.37.unsorted.vcf.gz'
	# params: 
	# 	chain=config['liftover']['chain'],
	# 	reference=config['liftover']['reference_fasta']
	# shell: """
	# 		CrossMap.py vcf \
	# 		{params.chain} {input} \
	# 		{params.reference} {output} \
	# 		--compress"""

use rule preannotation_prep_vcf_germline as preannotation_prep_vcf_somatic with:
	input: 
		'results/{run}/somatic/{patient}/mutect2.37.unsorted.vcf.gz'
	output: 
		'results/{run}/somatic/{patient}/mutect2.37.sorted.vcf.gz'
	# shell: """
	# 		for i in {{1..22}} X Y MT; \
	# 			do \
	# 				echo chr${{i}} $i | \
	# 				sed 's/chrMT/chrM/'; \
	# 			done | \
	# 		bcftools annotate --rename-chrs - {input} -Ou | \
	# 		bcftools sort -Oz -o {output}
	# 		"""


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
