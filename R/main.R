#加载R包
library(TwoSampleMR)
library(ggplot2)
library(MRPRESSO)
library(gwasvcf)
library(VariantAnnotation)
library(dplyr)

print("============================== 开始运行 ==============================")


####Day6-工具变量信息提取(方法一：openGWAS网站在线提取)####
#####6.1提取和暴露因素强相关的SNP作为工具变量#####
exposure_schizophrenia <- extract_instruments(
  outcomes='ieu-a-22', 
  clump=T,
  r2=0.001,
  kb=10000
) 

#####6.2提取SNP在结局中的信息#####
outcome_schizophrenia <- extract_outcome_data(
  snps=exposure_schizophrenia$SNP, 
  outcomes='ukb-b-19255',
  proxies = TRUE
)

####Day8-MR计算#####
#####8.1.harmonise#####
dat <- harmonise_data(
  exposure_dat = exposure_schizophrenia, 
  outcome_dat = outcome_schizophrenia
)

#####8.2.计算F值并去除弱工具变量#####
dat$eaf.exposure <- dat$eaf.outcome
PVEfx <- function(BETA, MAF, SE, N){
  pve <- (2*(BETA^2)*MAF*(1 - MAF))/((2*(BETA^2)*MAF*(1 - MAF)) + ((SE^2)*2*N*MAF*(1 - MAF)))
  return(pve) 
}
dat$EAF2 <- (1 - dat$eaf.exposure)
dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)
dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.exposure)
dat$FSTAT <- ((dat$samplesize.exposure - 1 - 1)/1)*(dat$PVE/(1 - dat$PVE)) 

write.csv(dat,"results/instrumental_dat.csv")
dat <- dat[dat$FSTAT > 10, ] 
saveRDS(dat,"results/dat.rds")


print("======================== 从这里开始会非常耗时 ========================")


#####8.3MR计算#####
res <- mr(dat)
res
write.csv(res,"results/results.csv")
OR <-generate_odds_ratios(res)
OR
write.csv(OR,"results/results_OR.csv")

######8.4MR结果可视化#####

mr_scatter_plot(res,dat)
ggsave("results/scatter_plot.pdf",width = 6,height = 6)
res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)
ggsave("results/forest_plot.pdf",width = 8,height = 8)


######Day9敏感性检验######
#异质性-Q检验
mr_heterogeneity(dat)

#漏斗图
res_single <- mr_singlesnp(dat)
p1 <- mr_funnel_plot(res_single)
p1[[1]]
ggsave("results/funnel_plot.pdf",width = 6,height = 6)

#leave-one-out plot
res_loo <- mr_leaveoneout(dat)
p2 <- mr_leaveoneout_plot(res_loo)
p2[[1]]
ggsave("results/leave_one_out_plot.pdf",width = 8,height = 8)



#多效性检验（Egger回归）
mr_pleiotropy_test(dat)

# Load a simulated toy dataset
data(SummaryStats)



head(dat)

# Run MR-PRESSO global method
mr_presso(BetaOutcome = "beta.outcome", 
          BetaExposure = "beta.exposure", 
          SdOutcome = "se.outcome", 
          SdExposure = "se.exposure", 
          OUTLIERtest = TRUE, 
          DISTORTIONtest = TRUE, 
          data = dat, 
          NbDistribution = 1000,  
          SignifThreshold = 0.05)


print("========================= 开始提取 ieu-a-22 =========================")
####Day7-工具变量信息提取(方法二：本地GWAS数据集提取)####
#####7.1本地提取暴露强相关SNP作为工具变量#####
vcf_file <- "local/ieu-a-22.vcf.gz"
vcf_data <- readVcf(vcf_file)
schizo_vcf <- vcf_to_tibble(vcf_data, id = NULL)
head(schizo_vcf)
schizo_vcf$P <- 10^(-schizo_vcf$LP)
write.csv(schizo_vcf,"results/schizo_vcf.csv")
saveRDS(schizo_vcf,"results/schizo_vcf.rds")
rm(vcf_data)
schizo_vcf <- readRDS("results/schizo_vcf.rds")

schizo_vcf_filtered <- schizo_vcf %>%  filter(P < 5e-8)#筛选出p< 5e-8的SNPs
#整理为TwoSampleMR包暴露数据的格式
exposure_dat <- format_data(
  schizo_vcf_filtered,
  type='exposure',
  snp_col = "rsid",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P"
)
exposure_schizophrenia <-clump_data(exposure_dat,clump_r2=0.001,clump_kb=10000)



print("======================== 开始提取 ubk-b-19255 ========================")
#####7.2本地提取结局数据#####
vcf_file <- "local/ukb-b-19255.vcf.gz"
vcf_data <- readVcf(vcf_file)
outcome_vcf <- vcf_to_tibble(vcf_data, id = NULL)
head(outcome_vcf)

outcome_vcf$P <- 10^(-outcome_vcf$LP)
write.csv(outcome_vcf,"results/gwas_summary.csv")
rm(vcf_data)
# outcome_vcf <- readRDS("outcome_vcf.rds")

outcome_dat <- read_outcome_data(
  snps = exposure_schizophrenia$SNP,
  filename = "results/gwas_summary.csv",
  sep = ",",
  snp_col = "rsid",
  beta_col = "ES",
  se_col = "SE",
  eaf_col = "AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "P")
str(outcome_dat)
head(outcome_dat)

print("============================== 运行完毕 ==============================")