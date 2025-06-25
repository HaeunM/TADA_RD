
# Original code : Haeun Moon (2025-06-25)

require(dplyr); require(openxlsx); require(ebmc)

### Change to working directory where the folder was cloned

source('functions.R')
genetable_info=read.xlsx("data/genetable_info.xlsx")

########################################################################
### Loading gene count tables 
### Family_based data : ASC published data (see Supplementary material of Fu et al.(2022))
### case, control data : simulated data (output from realistic data generation/simulation-of-denovo-and-inherited-ptv-in-probands-and-siblings_2025-06-09)
########################################################################

### denovo + inherited data, ASC published (Fu et al.)

asc_family_ptv <- read.xlsx("data/fu_suppl_table5.xlsx")%>%select(gene_id, dn.ptv, dn.ptv.sib, in.t.ptv, in.u.ptv)
n_prob_asc <- 7287+288; n_sib_asc <- 2348+11
genetable_d=left_join(genetable_info, asc_family_ptv, by=c("Gene_ID"="gene_id"))

### control data, simulated, for TADA_CC 

var_sibling= read.table("data/SIBLINGS-simulated-ptv-variants.txt", sep="\t", head=T)
control_var=var_sibling%>%
  group_by(gene)%>%
  summarize(control=n())
genetable_c=left_join(genetable_d, control_var, by=c("Gene_ID"="gene"))
genetable_c$control[which(is.na(genetable_c$control))]=0

### case data, simulated 

var_proband= read.table("data/PROBANDS-simulated-ptv-variants.txt", sep="\t", head=T)

case_var=var_proband%>%
  select("gene", "ccr_pct", "gnomAD.MAF")%>%
  rename(gene_id=gene, CCR=ccr_pct, Freq=gnomAD.MAF)%>%
  filter(gene_id %in% genetable_c$Gene_ID)
case_var$CCR[which(is.na(case_var$CCR))]=0

### Attach gene-level covariates 

case_var_com=left_join(case_var, data.frame(genetable_c$Gene_ID, genetable_c$LOEUF, genetable_c$exp_lof, genetable_c$obs_lof, genetable_c$FDR_TADA_DD), by=c("gene_id"="genetable_c.Gene_ID"))%>%
  rename(LOEUF=genetable_c.LOEUF, exp_lof=genetable_c.exp_lof, obs_lof=genetable_c.obs_lof,FDR_TADA_DD=genetable_c.FDR_TADA_DD )

### Infer de novo with ClassDn with a threshold 0.7

load("data/ClassDn_underbagging.Rdata")
PTV_ub=my_model
case_var_com$pred_rb=predict(PTV_ub, case_var_com, type = "prob")

c=0.7
case_gene=case_var_com%>%
          group_by(gene_id)%>%
          summarize(case=n(), ldn=sum(pred_rb>c), lin=case-ldn)

### Combine the counts into the table

genetable=left_join(genetable_c, case_gene, by=c("Gene_ID"="gene_id"))
genetable$case[which(is.na(genetable$case))]=0
genetable$ldn[which(is.na(genetable$ldn))]=0
genetable$lin[which(is.na(genetable$lin))]=0

########################################################################
### TADA_CC outcome calculation
### Source : https://github.com/talkowski-lab/TADA_2022/blob/main/run_TADA.R
########################################################################

### Setting sample sizes

n_asc_trio_case=6430; n_asc_trio_control=2179;  
n_case= 8000; n_control=2500

###Bayes factors calculation

beta.dn <- 0.2

BF_dn_ptv_asc <- BF_DN_SNV(count_case=genetable$dn.ptv, count_con=genetable$dn.ptv.sib, n_case=n_asc_trio_case, n_con=n_asc_trio_control, mut=genetable$mut.ptv, gamma.dn=genetable$prior.dn.ptv)
BF_in_ptv_asc <- BF_CC_SNV(count_case=genetable$in.t.ptv, count_con=genetable$in.u.ptv, n_case=n_asc_trio_case, n_con=n_asc_trio_control, mut=genetable$mut.ptv, gamma.cc=genetable$prior.in.ptv)
BF_cc_ptv <- BF_CC_SNV(count_case=genetable$case, count_con=genetable$control, n_case, n_control, mut=genetable$mut.ptv, gamma.cc=genetable$prior.cc.ptv)

genetable$BF_TADA_family=pmax(BF_dn_ptv_asc,1)*pmax(BF_in_ptv_asc,1)
genetable$BF_TADA_CC=pmax(BF_dn_ptv_asc,1)*pmax(BF_in_ptv_asc,1)*pmax(BF_cc_ptv,1)

###FDR calculation

genetable$qval_TADA_CC=Bayesian.FDR(genetable$BF_TADA_CC, pi0 =1-0.06)
genetable$qval_TADA_family=Bayesian.FDR(genetable$BF_TADA_family, pi0 =1-0.06)

##################################################
### TADA_CC outcome calculation
##################################################

### Load sensitivity and specificity
load("data/accuracy_underbagging")
w1= 0.335; w2=0.990

genetable$BF_TADA_RD=genetable$BF_TADA_family*BF_RD_CC(genetable$ldn, genetable$lin, genetable$alpha_risk, genetable$alpha_nonrisk, genetable$beta_risk, genetable$beta_nonrisk, w1, w2)
genetable$qval_TADA_RD=Bayesian.FDR(genetable$BF_TADA_RD, pi0 =1-0.06)


##################################################
### Result
### qval_TADA_family : FDR value from family-based data only
### qval_TADA_CC : FDR value computed using the TADA_CC method (Xin et al (2013))
### qval_TADA_RD : FDR value computed using the TADA_RD method
### qval_Fu : FDR value with multiple types of variant using the TADA_CC method (Fu.et al (2022))
##################################################

genetable_result=genetable%>%
  select(Gene_ID,gene, dn.ptv, dn.ptv.sib, in.t.ptv, in.u.ptv,  case, ldn, lin, BF_TADA_family,BF_TADA_CC, BF_TADA_RD, qval_TADA_family, qval_TADA_CC, qval_TADA_RD, qval_Fu)%>%
  rename(proband_dn=dn.ptv, sibling_dn=dn.ptv.sib, proband_in=in.t.ptv, sibling_in=in.u.ptv)

genetable_result%>%
  filter(qval_TADA_RD<0.05&qval_TADA_CC>0.05)%>%
  select(gene, proband_dn, proband_in, case, ldn, qval_TADA_CC, qval_TADA_RD)%>%
  mutate(qval_TADA_CC=round(qval_TADA_CC,3), qval_TADA_RD=round(qval_TADA_RD))

