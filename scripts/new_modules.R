# CARE-INDIA (Carbapenem resistance among Gram-negative ESKAPEE pathogens in India, 2017-2024)
# Supplementary figure/table modules: specimen-stratified antibiogram, isolation share, resistance-gene trends, and K. pneumoniae genotype-phenotype
# Authors: Abhishek Akella, Anand Srinivasan  |  Dept of Pharmacology, AIIMS Bhubaneswar, India
# License: MIT (see LICENSE)

suppressPackageStartupMessages({library(tidyverse)})
root <- getwd()
for (d in c("output/specimen","output/burden","output/genes"))
  dir.create(d, showWarnings=FALSE, recursive=TRUE)
th <- theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold",size=13), legend.position="top")
wslope <- function(d){ if(nrow(d)<3||length(unique(d$year))<3) return(NA_real_)
  f<-lm(res_pct~year,data=d,weights=tested); unname(coef(f)[2]) }

sp <- read_csv("data/specimen_stratified.csv",show_col_types=FALSE)
mero <- sp %>% filter(drug=="Meropenem")
slopes <- mero %>% group_by(organism,specimen) %>% summarise(slope=wslope(cur_data()),
            r2024=res_pct[year==2024], .groups="drop")
write_csv(slopes,"output/specimen/table_specimen_meropenem_slopes.csv")
cat("\n== SPECIMEN meropenem resistance slopes & 2024 ==\n"); print(slopes)
p1 <- ggplot(mero,aes(year,res_pct,colour=specimen))+geom_line(linewidth=1)+geom_point(size=2)+
  facet_wrap(~organism)+scale_colour_manual(values=c(Systemic="#C1272D",Urine="#0071BC"))+
  labs(title="Meropenem resistance: urine vs systemic isolates",x=NULL,y="Resistance (%)",colour="Specimen")+th
ggsave("output/specimen/fig_specimen_meropenem.png",p1,width=8,height=4.2,dpi=300)

oral <- c("Nitrofurantoin","Fosfomycin","Trimethoprim-sulfamethoxazole","Ciprofloxacin","Amikacin","Meropenem")
uti <- sp %>% filter(year==2024,specimen=="Urine",drug %in% oral) %>%
  mutate(drug=factor(drug,levels=oral))
write_csv(uti %>% select(organism,drug,tested,susc_pct),"output/specimen/table_oral_uti_agents_2024.csv")
p2 <- ggplot(uti,aes(drug,susc_pct,fill=organism))+geom_col(position="dodge")+
  geom_hline(yintercept=80,linetype="dashed",colour="grey40")+
  scale_fill_manual(values=c("E. coli"="#2E86AB","K. pneumoniae"="#A23B72"))+
  labs(title="2024 urinary-isolate susceptibility (oral & key agents)",x=NULL,y="Susceptibility (%)",fill=NULL)+
  th+theme(axis.text.x=element_text(angle=25,hjust=1))
ggsave("output/specimen/fig_specimen_oral_uti_2024.png",p2,width=8,height=4.5,dpi=300)
keyd<-c("Meropenem","Imipenem","Ertapenem","Amikacin","Ciprofloxacin","Piperacillin-tazobactam","Nitrofurantoin","Fosfomycin","Trimethoprim-sulfamethoxazole")
tab24 <- sp %>% filter(year==2024,drug %in% keyd) %>%
  select(organism,specimen,drug,susc_pct) %>%
  pivot_wider(names_from=specimen,values_from=susc_pct)
write_csv(tab24,"output/specimen/table_specimen_2024_antibiogram.csv")
cat("\n== 2024 urine vs systemic susceptibility (key drugs) ==\n"); print(tab24,n=Inf)

iso <- read_csv("data/isolation_trends.csv",show_col_types=FALSE)
burden <- iso %>% group_by(organism) %>%
  summarise(y2017=isolation_pct[year==2017], peak=max(isolation_pct),
            peak_yr=year[which.max(isolation_pct)], y2023=isolation_pct[year==2023],.groups="drop")
write_csv(burden,"output/burden/table_isolation_share.csv")
cat("\n== ISOLATION SHARE (organism % of culture-positives) ==\n"); print(burden)
ab <- iso %>% filter(organism=="A. baumannii")
ab_slope <- unname(coef(lm(isolation_pct~year,ab))[2])
cat(sprintf("A. baumannii isolation-share slope: %+.2f pp/yr\n",ab_slope))
p3 <- ggplot(iso,aes(year,isolation_pct,colour=organism))+geom_line(linewidth=1)+geom_point(size=2)+
  labs(title="Organism share of Gram-negative culture-positives, 2017-2023",
       x=NULL,y="Isolation share (%)",colour=NULL)+th
ggsave("output/burden/fig_isolation_share.png",p3,width=8,height=4.5,dpi=300)

gt <- read_csv("data/gene_trends.csv",show_col_types=FALSE)
carb <- c("NDM","OXA-48","KPC","VIM","IMP")
kpg <- gt %>% filter(organism=="K. pneumoniae",gene %in% carb)
p4 <- ggplot(kpg,aes(year,prevalence_pct,colour=gene))+geom_line(linewidth=1)+geom_point(size=2)+
  labs(title="K. pneumoniae carbapenemase gene prevalence (AMRSN)",x=NULL,y="Prevalence (%)",colour="Gene")+th
ggsave("output/genes/fig_kp_carbapenemase_trends.png",p4,width=8,height=4.5,dpi=300)
dom <- gt %>% filter(gene %in% carb) %>% group_by(organism) %>%
  filter(year==max(year)) %>% arrange(organism,desc(prevalence_pct))
write_csv(dom,"output/genes/table_dominant_carbapenemase_latest.csv")
cat("\n== Dominant carbapenemase genes (latest yr) ==\n"); print(dom,n=Inf)
kp_ndm <- gt %>% filter(organism=="K. pneumoniae",gene=="NDM") %>% select(year,ndm=prevalence_pct)
kp_mero <- read_csv("data/data_full_long.csv",show_col_types=FALSE) %>%
  filter(organism=="K. pneumoniae",drug=="Meropenem") %>% select(year,mero_res=res_pct)
gp <- inner_join(kp_ndm,kp_mero,by="year")
write_csv(gp,"output/genes/table_genotype_phenotype_kp.csv")
cat("\n== K.p NDM prevalence vs meropenem resistance ==\n"); print(gp)
if(nrow(gp)>=3){ cc<-cor(gp$ndm,gp$mero_res); cat(sprintf("Pearson r (NDM vs mero-res): %.2f\n",cc)) }
p5 <- ggplot(gp,aes(year))+
  geom_line(aes(y=mero_res,colour="Meropenem resistance"),linewidth=1)+geom_point(aes(y=mero_res,colour="Meropenem resistance"))+
  geom_line(aes(y=ndm,colour="NDM prevalence"),linewidth=1,linetype="dashed")+geom_point(aes(y=ndm,colour="NDM prevalence"))+
  scale_colour_manual(values=c("Meropenem resistance"="#C1272D","NDM prevalence"="#F15A24"))+
  labs(title="K. pneumoniae: NDM prevalence vs meropenem resistance",x=NULL,y="Percent",colour=NULL)+th
ggsave("output/genes/fig_genotype_phenotype_kp.png",p5,width=8,height=4.5,dpi=300)

cat("\n=== NEW MODULES COMPLETE ===\n")
