
# Original code : Haeun Moon (2025-06-25)

################################################################################################
### Wrapper function to obtain Bayes Factors for case data, over a list of genes for a given variant type
### Makes calls to log.bayes.factor.dn for each gene
## Input: 
# count_ldn: a vector of likely de novo variant counts over each gene in cases
# count_lin: a vector of likely inherited variant counts over each gene in cases
# Prior distribution of de novo ratio: 
# If risk gene, ~ beta(alpha_risk, beta_risk)
# If nonrisk gene, ~ beta(alpha_nonrisk, beta_nonrisk)
# w1 : sensitivity of the applied classifier to infer de novo status
# w2 : specificity of the applied classifier to infer de novo status
## Output:
# BF: a vector of Bayes factors of risk gene evidence for each gene in a case data 
################################################################################################

BF_RD_CC <- function(count_ldn, count_lin, alpha_risk, alpha_nonrisk, beta_risk, beta_nonrisk, w1, w2){
  
  ibeta0=seq(0.01,0.5,0.005); ibeta1=seq(0.3,0.9,0.02); 
  
  int_value_num=matrix(0, length(ibeta1), length(count_ldn))
  int_density_num=matrix(0, length(ibeta1), length(count_ldn))
  int_value_deno=matrix(0, length(ibeta0), length(count_ldn))
  int_density_deno=matrix(0, length(ibeta0), length(count_ldn))
  
  for (j1 in 1:length(ibeta1)){
    p1=ibeta1[j1]; PD=w1*p1+(1-w2)*(1-p1); 
    int_value_num[j1,]=PD^count_ldn*(1-PD)^count_lin
    int_density_num[j1,]=dbeta(p1,genetable$alpha_risk,genetable$beta_risk)
  }
  
  for (j0 in 1:length(ibeta0)){
    p0=ibeta0[j0]; PDnot=w1*p0+(1-w2)*(1-p0); 
    int_value_deno[j0,]=PDnot^count_ldn*(1-PDnot)^count_lin
    int_density_deno[j0,]=dbeta(p0,genetable$alpha_nonrisk,genetable$beta_nonrisk)
  }
  
  BF=(colSums(int_value_num*int_density_num)/colSums(int_density_num))/(colSums(int_value_deno*int_density_deno)/colSums(int_density_deno))
  BF[BF<1] <- 1
  BF[is.na(BF)] <- 1
  
  return(BF)
}





## Functions used to apply the TADA framework

################################################################################################
### Wrapper function to obtain de novo Bayes Factors for SNV/indels, over a list of genes for a given variant type
### Makes calls to log.bayes.factor.dn for each gene
## Input: 
# count_case: a vector of de novo variant counts over each gene in cases
# count_con: a vector of de novo variant counts over each gene in controls
# n_case: the number of cases
# n_con: the number of controls
# mut: a vector of mutation rates for each gene
# Prior distribution of RR: gamma ~ Gamma(gamma.dn*beta.dn, beta.dn)
## Output:
# BF: a vector of Bayes factors of de novo evidence from SNV/indel for each gene for a given variant type
################################################################################################
BF_DN_SNV <- function(count_case, count_con, n_case, n_con, mut, gamma.dn, beta.dn=.2){
  ### If there are siblings, evaluate BF evidence for siblings to offset evidence for probands
  if(n_con>0){
    BF_con <- exp(sapply(1:length(mut), function(i) log.bayes.factor.dn(x.dn=count_con[i], n.dn=n_con, mu=mut[i], beta.dn=beta.dn, gamma.dn=gamma.dn[i])))
  }else{
    BF_con <- rep(1, length(mut))
  }
  BF_con[BF_con<1] <- 1
  
  ### BF evidence calculation for probands
  BF_case <- exp(sapply(1:length(mut), function(i) log.bayes.factor.dn(x.dn=count_case[i], n.dn=n_case, mu=mut[i], beta.dn=beta.dn, gamma.dn=gamma.dn[i])))
  BF_case[BF_case<1] <- 1
  
  ### Combine proband and sibling BF evidence
  ### Enforce BF back to 1 if no observations in probands or siblings
  BF <- pmax(BF_case/BF_con, 1) 
  BF[which(count_case==0 & count_con==0)] <- 1
  BF[is.na(BF)] <- 1
  return(BF)
}

################################################################################################
### functions downloaded from https://github.com/talkowski-lab/TADA_2022/blob/main/functions.R
################################################################################################

################################################################################################
### log Bayes factor calculation for de novo variant contribution to a gene
## Input:
# x.dn: the de novo variant count 
# n.dn: the sample size (number of families)
# mu: the mutation rate (of this type of mutational events in this gene)
# Prior distribution of risk: gamma ~ Gamma(gamma.dn*beta.dn, beta.dn)
## Output:
# lBF: log Bayes factor evidence for de novo variants of this type in this gene
################################################################################################
log.bayes.factor.dn <- function(x.dn, n.dn, mu,gamma.dn, beta.dn){
  marg.lik0 <- dpois(x.dn, 2*n.dn*mu, log=TRUE)
  marg.lik1 <- dnbinom(x.dn, gamma.dn*beta.dn, beta.dn/(beta.dn+2*n.dn*mu), log=TRUE)
  lBF <- marg.lik1-marg.lik0
  return (lBF=lBF)
}



################################################################################################
### Bayes factor calculation for case-control/inherited contribution to a gene
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# gamma.cc is a vector of the prior mean on risk
# Prior distribution of q|H1: Gamma(rho1, nu1)
# Prior distribution of q|H0: Gamma(rho0, nu0)
## Output:
# BF: Bayes factor evidence for case-control/inherited contribution to a gene
################################################################################################
bayes.factor.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1, rho0, nu0){
  marglik0.cc <- evidence.null.cc(x.cc, n.cc, rho0, nu0)
  marglik1.cc <- evidence.alt.cc(x.cc, n.cc, gamma.cc, rho1, nu1)
  BF.cn <- marglik1.cc$cn / marglik0.cc$cn
  ## Assuming null and alt of the control model is identical
  if(is.na(BF.cn)){BF.cn<-1} 
  BF.ca <- marglik1.cc$ca / marglik0.cc$ca
  BF <- BF.cn * BF.ca
  return(BF=BF)
}



################################################################################################
### Helper function to bayes.factor.cc
### Calculates the evidence of case-control/inherited counts under null model
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# Prior distribution of q|H0: Gamma(rho0, nu0)
## Output:
# cn: marginal likelihood of control/untransmitted variant counts under null model
# ca: marginal likelihood of case/transmitted variant counts under null model
################################################################################################
evidence.null.cc <- function(x.cc, n.cc, rho0, nu0) {
  marglik0.ctrl.log <- log(dnbinom(x.cc$cn, rho0, nu0/(nu0+n.cc$cn)))
  marglik0.case.log <- log(dnbinom(x.cc$ca, rho0+x.cc$cn, (nu0+n.cc$cn)/(nu0+n.cc$cn+n.cc$ca)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  
  return (list(cn=exp(marglik0.ctrl.log), ca=exp(marglik0.case.log)))#, total=exp(marglik0.log)))
}

################################################################################################
### Helper function to bayes.factor.cc
### Calculates the evidence of case-control/inherited counts under alternate model
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# gamma.cc is a vector of the prior mean on risk
## Output:
# cn: marginal likelihood of control/untransmitted variant counts under alternate model
# ca: marginal likelihood of case/transmitted variant counts under alternate model
################################################################################################
evidence.alt.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1){
  marglik1.ctrl <- dnbinom(x.cc$cn, rho1, nu1/(nu1+n.cc$cn))
  marglik1.case <- dnbinom(x.cc$ca, rho1+x.cc$cn, (nu1+n.cc$cn)/(nu1+n.cc$ca*gamma.cc+n.cc$cn))
  marglik1 <- marglik1.ctrl * marglik1.case
  return (list(cn=marglik1.ctrl, ca=marglik1.case))#, total=marglik1))
}


################################################################################################
### Wrapper function to obtain case-control/inherited Bayes Factors for SNV/indels, over a list of genes for a given variant type
### Makes calls to bayes.factor.cc for each gene
## Input: 
# count_case: a vector of case/transmitted variant counts over each gene in cases
# count_con: a vector of control/untransmitted variant counts over each gene in controls
# n_case: the number of cases
# n_con: the number of controls (case-control analysis) or number of cases (inherited analysis)
# mut: a vector of mutation rates for each gene
# gamma.cc: a vector of the prior mean on risk
# nu: a hyperparameter that controls nu0 and nu1 in determining prior distributions, passed to bayes.factor.cc
## Output:
# BF: a vector of Bayes factors of case-control/inherited evidence from SNV/indel for each gene for a given variant type
################################################################################################
BF_CC_SNV <- function(count_case, count_con, n_case, n_con, mut, gamma.cc, nu=5000){
  rho.in = nu * sum(count_con, na.rm=T) / (2*length(mut)*n_con)
  rho.in = as.numeric(rho.in)*mut/mean(mut, na.rm=TRUE)
  BF = sapply((1:length(mut)), function(i){
    bayes.factor.cc(x.cc = data.frame(ca=count_case[i], cn=count_con[i]), n.cc=data.frame(ca=n_case, cn=n_con), gamma.cc=gamma.cc[i], rho1=rho.in[i], nu1=nu, rho0=rho.in[i], nu0=nu)
  })
  
  ### Enforce BF back to 1 if no mutation rate or no observed variant counts
  BF[mut==0] <- 1
  BF[which((count_case+count_con)==0)] <- 1
  return(BF) 
}


################################################################################################
### Function to transform total BF evidence to FDR values
### Bayesian FDR control (PMID:19822692, Section2.3)
## Input: 
# BF: total Bayes factor matrix (each row is a gene, multiple columns allowed)
# pi0: estimated proportion of non-risk genes
## Output:
# FDR: a vector of FDR values corresponding to the rows of the input matrix
################################################################################################
Bayesian.FDR <- function(BF, pi0) {
  # order the BF in decreasing order, need to retain order to get results back in proper order 
  i.order=order(BF, decreasing = T)
  BF=BF[i.order]
  # convert BFs to PPA (posterior probability of alternative model)
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model
  
  # the FDR at each PPA cutoff
  FDR=cumsum(q0)/(1:length(BF))
  
  # reorder to the original order
  FDR[i.order]=FDR
  
  return (FDR=FDR)
}

