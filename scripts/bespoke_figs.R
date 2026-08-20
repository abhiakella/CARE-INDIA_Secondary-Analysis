# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Build bespoke figures (antibiogram heatmap, therapeutic tier map, carbapenemase gene prevalence) from long-format antibiogram and gene-trend data.
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressPackageStartupMessages({library(tidyverse)})
th <- theme_minimal(base_size=11)+theme(panel.grid=element_blank(),plot.title=element_text(face="bold",size=12))
dl <- read_csv("data/data_full_long.csv",show_col_types=FALSE)
drug_ord <- c("Piperacillin-tazobactam","Cefotaxime","Ceftazidime","Cefepime","Ertapenem","Imipenem",
              "Meropenem","Amikacin","Gentamicin","Tobramycin","Ciprofloxacin","Levofloxacin","Colistin","Minocycline")
org_ord <- c("E. coli","K. pneumoniae","A. baumannii","P. aeruginosa","Enterobacter spp.")

## FIG 2: 2024 antibiogram heatmap
d24 <- dl %>% filter(year==2024) %>% mutate(drug=factor(drug,levels=rev(drug_ord)),organism=factor(organism,levels=org_ord))
f2 <- ggplot(d24,aes(organism,drug,fill=susc_pct))+geom_tile(colour="white",linewidth=.4)+
  geom_text(aes(label=round(susc_pct)),size=3)+
  scale_fill_gradientn(colours=c("#B2182B","#EF8A62","#FDDBC7","#D1E5F0","#2166AC"),limits=c(0,100),name="Susc %")+
  labs(title="2024 cumulative antibiogram (% susceptibility)",x=NULL,y=NULL)+th+
  theme(axis.text.x=element_text(angle=20,hjust=1))
ggsave("output/fig2_antibiogram_heatmap.png",f2,width=7.5,height=6,dpi=300)

## FIG 3: therapeutic tier map (2024)
tier <- d24 %>% mutate(tier=cut(susc_pct,c(-1,10,30,60,101),labels=c("Inactive","Not recommended","Targeted","Viable")))
f3 <- ggplot(tier,aes(organism,drug,fill=tier))+geom_tile(colour="white",linewidth=.4)+
  scale_fill_manual(values=c("Inactive"="#B2182B","Not recommended"="#EF8A62","Targeted"="#FEE08B","Viable"="#1A9850"),name="Tier",drop=FALSE)+
  labs(title="Therapeutic tier classification, 2024",x=NULL,y=NULL)+th+
  theme(axis.text.x=element_text(angle=20,hjust=1))
ggsave("output/fig3_tier_map.png",f3,width=7.5,height=6,dpi=300)

## FIG 4: India carbapenemase gene prevalence (latest reported year)
gt <- read_csv("data/gene_trends.csv",show_col_types=FALSE)
carb <- c("NDM","OXA-48","KPC","VIM","IMP")
gl <- gt %>% filter(gene %in% carb) %>% group_by(organism) %>% filter(year==max(year)) %>% ungroup() %>%
  mutate(organism=factor(organism,levels=org_ord),gene=factor(gene,levels=carb))
f4 <- ggplot(gl,aes(gene,prevalence_pct,fill=gene))+geom_col()+facet_wrap(~organism,nrow=1)+
  scale_fill_brewer(palette="Set2",guide="none")+
  labs(title="India carbapenemase gene prevalence, AMRSN (latest reported year)",x=NULL,y="Prevalence (%)")+th+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8))
ggsave("output/fig4_india_carbapenemase.png",f4,width=9,height=3.6,dpi=300)
cat("bespoke figures written: fig2_antibiogram_heatmap, fig3_tier_map, fig4_india_carbapenemase\n")
