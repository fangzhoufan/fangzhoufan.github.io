
# DM2 and OP mendelian randomization
**Author:** 范方舟 Fangzhou Fan  
**Date:** 2023-07-14


# 1.包的载入

```{r echo=TRUE}
options(warn = -1)
library(MRInstruments)
library(remotes)
library(TwoSampleMR)
library(ggplot2)
library(ggsci)
library(MendelianRandomization)
library(MRPRESSO)
```

# 2.数据收集

```{r}
ao <- available_outcomes()
##exposure
ao1<-ao%>%
  subset(grepl('diabetes',trait,ignore.case = TRUE))
##outcomes
ao2<-ao%>%
  subset(grepl('osteoporosis',trait,ignore.case = TRUE))
```

```{r echo=TRUE}
ao1$id[ao1$trait=='Type 2 diabetes, strict (exclude DM1)']
```

```{r results = 'asis'}
ao2$id[ao2$trait=='Non-cancer illness code, self-reported: osteoporosis']
```

```{r echo=TRUE}
##exposure
DM2_exp_dat_all <- extract_instruments(outcomes='finn-b-E4_DM2_STRICT',clump=FALSE, access_token = NULL)

DM2_exp_dat <- clump_data(DM2_exp_dat_all,clump_r2=0.0001,clump_kb=10000,clump_p1 = 5*10^-8)
##>默认有统计显著性的阈值设置为“P < 5 × 10 − 8;LD r2<0.0001，kb = 10000“来识别与DM2相关的SNPs
dim(DM2_exp_dat_all)
dim(DM2_exp_dat)
```

```{r echo=TRUE}
##Outcomes
OP_out_dat <- extract_outcome_data(snps = DM2_exp_dat$SNP, outcomes = 'ukb-b-12141')
dat <- harmonise_data(DM2_exp_dat, OP_out_dat)
```

```{r}
dim(OP_out_dat)
dim(dat)
```

```{r echo=TRUE}
head(dat)
```

# 3.数据检验，MR模型分析

```{r echo=TRUE}
###MR评估模型——常见五种
knitr::kable(mr_method_list())
```

```{r echo=TRUE}
##异质性分析
mr_heterogeneity(dat)
mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw"))
```

```{r echo=TRUE}
##>如果发现存在p小于0.05，存在异质性
##>可以使用随机效应模型
mr(dat,method_list=c('mr_ivw_mre')) #使用随机效应模型
```

```{r echo=TRUE}
##水平多效性分析
c<-mr_pleiotropy_test(dat)
knitr::kable(c)
run_mr_presso(dat)
```

# 4.MR可视化

```{r echo=TRUE}
res <- mr(dat)
```

```{r}
knitr::kable(res)
```

```{r}
res$pval<-formatC(res$pval,digits = 2)
```

```{r}
res$pval<-formatC(res$pval,digits = 2)
res$lo<-round(exp(res$b - 1.96*res$se), 4)
res$uo<-round(exp(res$b + 1.96*res$se), 4)
res$OR<-round(exp(res$b), 4)

res$method<-factor(res$method,levels = res$method)

p <- ggplot(res, aes(OR, method)) 

annotation<-data.frame(x=c(1.10, 1.10, 1.10, 1.10, 1.10,1.10),
                       y=c(5.49,5, 4, 3, 2, 1),
                       label=c('P Value',res$pval[1],res$pval[2],res$pval[3],res$pval[4],res$pval[5]))
annotation$z<-c(rep(1.07,6))

for (i in 1:5){
  res$ORlabel[i]<-paste(res$OR[i],'(',res$lo[i],"-",res$uo[i],")",sep = "")
}
annotation$OR[1]<-'OR(95%)'
annotation$OR[2:6]<-res$ORlabel



p<-p + geom_point(size=3.6) +
  geom_errorbarh(aes(xmax =uo, xmin = lo), height = 0.4) +
  scale_x_continuous(limits= c(0.99, 1.10), breaks= seq(0.9, 1.1, 0.01)) +
  geom_vline(aes(xintercept = 0.99)) +
  xlab('Odds Ratio ') + ylab(' ')+theme_minimal()+
  geom_text(data=annotation,aes(x=x,y=y,label=label))+
  geom_text(data=annotation,aes(x=z,y=y,label=OR))
p
```

```{r}
p1 <- mr_scatter_plot(res, dat)
p1[[1]]+scale_fill_lancet()+scale_color_lancet()
```

```{r echo=TRUE}
res_single <- mr_singlesnp(dat, single_method = "mr_meta_fixed")##没有异质性的时候
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]+scale_fill_lancet()+scale_color_lancet()
```

```{r results='markup'}
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]+scale_fill_lancet()+scale_color_lancet()
```

```{r}
res_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(res_single)
p4[[1]]+scale_fill_lancet()+scale_color_lancet()
```
