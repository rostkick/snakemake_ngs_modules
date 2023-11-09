log <- file(snakemake@log[[1]], open="wt")
sink(log)

defaultW <- getOption("warn")
options(warn = -1)

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(vcfR))


vcf <- read.vcfR(snakemake@input[[1]])

fields = c('CONTQ', 'DP', 'GERMQ',
			'NALOD', 'POPAF', 'ROQ',
			'SEQQ', 'STRANDQ', 'STRQ',
			'TLOD') # PON was excluded

df <- suppressMessages(vcfR2tidy(vcf, 
								single_frame = TRUE,
								info_fields = fields))

na_fields <- df$dat %>% dplyr::select(all_of(fields)) %>% 
	dplyr::select(which(colMeans(is.na(.)) == 1)) %>% colnames()

plot_data <- df$dat %>% dplyr::select(all_of(fields)) %>% 
	dplyr::select(!which(colMeans(is.na(.)) == 1)) %>% 
	mutate(across(everything(.), ~as.numeric(.))) %>% 
	tidyr::gather()

pdf(snakemake@output[[1]], width = 7, height = 5)
ggplot(plot_data, aes(x=value, group=key, fill=key))+
						geom_histogram(bins=30)+
						facet_wrap(~key, scales='free',
									strip.position='bottom')+
						theme_bw()+
						scale_y_continuous('Count')+
						scale_x_continuous('')+
						scale_fill_npg(name='')+
						labs(title = paste0('Patient: ', snakemake@wildcards$patient))+
						theme(strip.background = element_blank(),
							strip.placement = "outside")
dev.off()

options(warn = defaultW)
