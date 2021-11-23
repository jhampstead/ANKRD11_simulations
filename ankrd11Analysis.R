library(data.table)
library(tidyverse)
library(ggpubr)

#Per nucleotide mutation rates taken from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1461236/pdf/10978293.pdf
mutation_rates <- fread('/Users/jh/Projects/penetranceANKRD11/data/per_nucleotide_mutation_rates.txt')
mutation_rates <- mutation_rates[,Scaled_Mutation_Rate:=Mutation_Rate*1000000]
infile <- fread('/Users/jh/Projects/penetranceANKRD11/data/ANKRD11_ENST00000301030.10_4')
total_variants_rd2 = 17
total_variants_other = 8
observed_arginines_rd2 = 12
observed_arginines_other = 0

## Parse gnomAD input
amino_acids <- data.table(symbol=c('R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W','*'),
                          name=c('Arg','His','Lys','Asp','Glu','Ser','Thr','Asn','Gln','Cys','Gly','Pro','Ala','Val','Ile','Leu','Met','Phe','Tyr','Trp','Ter'))

parseAminoAcid <- function( protein_consequence ) {
  
  aa_position <- unlist(strsplit(protein_consequence,"\\D+"))[2]
  ref_name <- gsub('p.','',unlist(strsplit(protein_consequence,"\\d+")))[1]
  alt_name <- gsub('p.','',unlist(strsplit(protein_consequence,"\\d+")))[2]
  
  ref_aa <- amino_acids[which(name==ref_name),symbol]
  alt_aa <- amino_acids[which(name==alt_name),symbol]
  
  return(list(aa_position,
              ref_aa,
              alt_aa))
  
}

## Parse simulated input

infile_parsed <- infile %>%
  filter(position_aa!='-') %>%
  mutate(lead_ref = dplyr::lead(ref,n = 3,order_by=as.numeric(position_aa)),
         is_cpg = ifelse(ref=='C' & lead_ref=='G',TRUE,FALSE)) %>%
  dplyr::select(-lead_ref) %>%
  filter(vep_csq=='missense_variant') %>%
  mutate(data_type='Simulated mutations') %>%
  data.table()

getScaledMutationRate <- function( ref, alt, is_cpg ) {
  
  cpg_status <- ifelse(is_cpg==TRUE,'CpG','Non-CpG')
  scaled_mutation_rate <- mutation_rates[Ref==ref & Alt==alt & CpG_Status==cpg_status,Scaled_Mutation_Rate]
  return(scaled_mutation_rate)

}

sampleMutations <- function( data, iterations, num_mutations ) {
  
  data_sampled <- vector('list',length=iterations)
  
  for( i in 1:iterations ) {
    
    data_sampled[[i]] <- data %>%
      sample_n(num_mutations,weight=scaled_mr,replace=TRUE) %>%
      mutate(iteration = i) %>%
      data.table()
    
  }
  
  data_out <- do.call(rbind,data_sampled)
  return(data_out)
  
}

infile_parsed <- infile_parsed[,scaled_mr:=getScaledMutationRate(ref,alt,is_cpg),by=1:nrow(infile_parsed)]

infile_rd2 <- infile_parsed %>%
  filter(position_aa >= 2369 & position_aa <= 2663) %>%
  data.table()

infile_other <- infile_parsed %>%
  filter(position_aa < 2369) %>%
  data.table()

## Differences in arginine mutations between region of interest (2369-2663) and other parts of the gene

set.seed(1234)

rd2_sampled <- sampleMutations(infile_rd2,10000,total_variants_rd2) %>%
  mutate(class='Expected_RD2') %>%
  data.table()

other_sampled <- sampleMutations(infile_other,10000,total_variants_other) %>%
  mutate(class='Expected_Other') %>%
  data.table()

out <- rbind(rd2_sampled,other_sampled) %>%
  group_by(iteration,class) %>%
  mutate(num_arg = sum(ref_aa=='R')) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(mean = mean(num_arg),
         sd = sd(num_arg),
         min = mean - sd,
         max = mean + sd) %>%
  data.table()

getPermutationP <- function( data, observed_class, observed_value ) {
  
  observed_data <- data.table(class=observed_class,
                              num_arg=observed_value,
                              iteration=0)
  
  temp <- rbind(data[,c('class','num_arg','iteration')],observed_data) %>%
    unique() %>%
    arrange(desc(num_arg),iteration) %>%
    filter(class==observed_class) %>%
    data.table()
  
  permutation_p <- which(temp[,iteration==0])/(max(temp[,iteration]) + 1)
  return(permutation_p)
  
}
p_rd2 <- getPermutationP(out,'Expected_RD2',observed_arginines_rd2)
p_other <- getPermutationP(out,'Expected_Other',observed_arginines_other)
label_dt <- data.table(class=c('RD2','Other'),
                       label=c(paste0('p = ',round(p_rd2,5)),paste0('p = ',round(p_other,5))),
                       xmin = c(1.3170,0.4514))
  
plot2 <- out %>%
  select(class,mean,min,max) %>%
  unique() %>%
  data.table()

plot2 <- plot2[,c('class','loc'):=as.list(unlist(strsplit(class,'_'))),by=1:nrow(plot2)]

observed_arg <- data.table(class=c('Observed','Observed'),
                           mean=c(12,0),
                           min=c(0,0),
                           max=c(0,0),
                           loc=c('RD2','Other'))

plot_out <- rbind(plot2,observed_arg)

#Regarding supplementary figure lay-out (see example Figure S5 attached):
#Figure title: Arial bold, 12
#A/B: Arial bold, 12
#Figure description: Arial 8 regular (above figure)
#Colours: in attached excel we listed colours we used in other figures, so best to choose any from these
#Your figure title: Figure S3: ANKRD11 missense variants affecting arginine residues in RD2 are overrepresented in the cohort and predicted to be more damaging.
#The arginine frequency cohort vs permutated is section A; the CADD comparison is section B.

plot2_labels <- c('RD2' = 'Repressor domain 2\n(amino acids 2369-2663)',
                  'Other' = 'All other amino acids')

p2 <- ggplot(data=plot_out,aes(x=mean,y=loc,colour=class)) +
  geom_pointrange(data=plot_out[class=='Expected'],aes(xmin=min,xmax=max),size=2) +
  geom_point(data=plot_out[class=='Observed'],size=9,shape=18) +
  geom_text(data =label_dt,aes(x=xmin,y=class,label=label),position = position_nudge(y=0.15),size=4,colour='black') +
  theme_pubr() +
  xlab('Mean number of arginine residues mutated') +
  ylab('') +
  theme(legend.title = element_blank(),
        axis.text = element_text(size=12),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        plot.title = element_text(size=12,face='bold'),
        legend.position = 'bottom') +
  scale_colour_manual(values=c('darkred','black'),labels=plot2_labels,
                      guide=guide_legend(override.aes = list(linetype=c(0,0),
                                                             shape=c(19,18)))) +
  scale_y_discrete(labels=plot2_labels)
  #ggtitle('Figure S3: ANKRD11 missense variants affecting arginine residues in RD2 are overrepresented in the cohort')

ggsave('ANKRD11_simulations.png',plot=last_plot(),device='png',dpi=300,width=10,height=7)
ggsave('ANKRD11_simulations.svg',plot=last_plot(),device='svg',dpi=300,width=10,height=7)