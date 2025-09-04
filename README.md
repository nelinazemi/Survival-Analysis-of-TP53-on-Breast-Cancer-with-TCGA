````markdown

The intuition behind Survival Analysis is how long individuals are likely to survive before a specific event occurs. The following questions are often the main focus of survival analysis experiments:  
1- What is the probablity that a certain patient will live until a specific timepoint?  
2- What factors influence the lifespan of patients?  
3- Is there a difference in survival between patients who recieved treatment A vs. those who recieved treatment B?  
4- What factors contribute to the differences in survival in this specific disease?  

---

## Censoring

One important factor in survival analysis is 'censoring'. Censoring is a type of missing data problem which is caused when the observation of the survival of specific patients goes unknown due to the following reasons:  
- Some patients may simply drop out of the study.  
- Some patients may lose follow-up. That means they survive until the end of the experiment but we don't know how long they lived after. So we don't know the exact survival time.  
- Some other event makes follow-up impossible.  

In all the mentioned cases, we need to censor the data since there's no exact survival time available.  

---

## Experiment Goal

In this experiment we want to see if high levels of TP53 can help us in predicting breast cancer.  

---

## Loading Libraries

We start by loading the necessary libraries:  

```r
library(BiocManager)
library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(tidyr)
library(dplyr)
library(GenomicDataCommons)
library(DESeq2)
````

---

## Downloading Data

We download the data from GDC portal:

```r
clinical_brca <- GDCquery_clinic("TCGA-BRCA")
```

or

```r
remotes::install_github("https://github.com/BioinformaticsFMRP/TCGAbiolinks",ref = "devel")
```

---

## Clinical Factors

We need three main factors to deal with in this experiment:
1- vital status: whether the patient was alive or dead at the time of their last follow-up.
2- days to last follow up: is the time from diagnosis until the most recent check.
3- days to death: is the time from diagnosis until the patient’s death if they had passed away.

So we extract the columns associated with these info:

```r
any(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_brca[,c(10,47,55)]
```

---

## Applying Censoring

If the status of a patient is 'Alive' at the time of the last follow-up, we extract the time of their last follow-up. If their status is 'Dead', we extract their time of death.

```r
clinical_brca$overall_survival <- ifelse(clinical_brca$vital_status == "Alive", clinical_brca$days_to_last_follow_up,
clinical_brca$days_to_death)
```

This way we apply censoring by only preserving data which captures the overall survival info of each patient. For now, we have a list of patients and their exact survival time.

---

## RNA-Expression Metadata

For the next step we need to get the RNA-expression metadata for our patients. We obtain it by using the following query:

```r
tumoric_rna_query = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")
```

```r
tumoric_rna_brca <- getResults(tumoric_rna_query)
```

This basically gives us all primary tumor RNA-Seq metadata for breast cancer (BRCA). In order to make the analysis faster, we only take the first 100 patients (cases):

```r
tumoric_cases <- tumoric_rna_brca$cases[1:100]
```

---

## Matching Normal Samples

In order to make comparisons between normal vs. tumoric tissues, we need matching normal tissue samples as well. Here, the barcode gives us the expression counts for matching normal tissues of the same 100 patients.

```r
query_brca <- GDCquery(
normal_rna_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumoric_cases)
```

```r
GDCdownload(normal_rna_query)
```

---

## Summarized Data

We can get the summerized info as well by running the following command:

```r
tcga_brca_data <- GDCprepare(normal_rna_query, summarizedExperiment = TRUE)
tcga_brca_matrix <- assay(tcga_brca_data, "unstranded")
tcga_brca_matrix[1:5,1:5]
```

The tgca\_brca\_data consists of genes, samples and their raw expression counts. We can also get the metadata of each gene and patient by running the following command:

```r
gene_metadata <- as.data.frame(rowData(tcga_brca_data))
cases <- as.data.frame(colData(tcga_brca_data))
```

---

## Normalization

In order to normalize our counts, we create a DESeq2 object, filter out weakly expressed genes, and transform the data so it’s comparable across samples.

```r
dds <- DESeqDataSetFromMatrix(countData = brca_matrix,
                              colData = cases,
                              design = ~ 1)

dds <- dds[rowSums(counts(dds)) >= 10]

vsd <- vst(dds, blind=FALSE)
brca_matrix_vst <- assay(vsd)
brca_matrix_vst[1:10,1:10]
```

design = \~ 1 means: no groups, no conditions to compare. This object is only being used to store and normalize counts (e.g., variance-stabilizing transform) so samples can be compared fairly in downstream steps.

---

## TP53 Expression

We create another dataframe out of the normalized counts and use genes as a seperate column named 'gene\_id'. We also use the sample names as a different column called 'case\_id'. The expression values are also stored in 'counts' column. We merge the gene metadata and brca\_matrix\_vst matrix by gene\_ids.

```r
brca_tp53 <- brca_matrix_vst %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = 'gene_id') %>%
  tidyr::pivot_longer(
    cols = -gene_id,
    names_to = "case_id",
    values_to = "counts"
  ) %>%
  dplyr::left_join(gene_metadata, by = "gene_id") %>%
  dplyr::filter(gene_name == "TP53")
```

---

## Stratifying Patients

We can get the strata of our patients using the median() function. If the counts are above the median, they're considered highly expressed whereas the others are low expressed.

```r
median_value <- median(brca_tp53$counts)
brca_tp53$strata <- ifelse(brca_tp53$counts >= median_value, "HIGH", "LOW")
```

Our brca\_tp53 dataframe consists of case\_ids, genes and their expression level. We already stored our patients with survival info in clinical dataset. Now we have to merge them based on patients.

The format of patients in brca\_tp53 (case\_id): TCGA-GM-A2DL-01A-11R-A18M-07
The format of patients in clinical\_brca (submitter\_id): TCGA-A7-A0DC

```r
brca_tp53$case_id <- gsub('-01.*', '', brca_tp53$case_id)
brca_tp53 <- merge(brca_tp53, clinical_brca, by.x='case_id', by.y='submitter_id')
```

---

## Kaplan–Meier Survival Curves

At the end of our experiment, we use survfit() function which builds Kaplan–Meier survival curves for HIGH vs LOW TP53 expression. ggsurvplot() plots them with p-values and risk tables. survdiff() also performs a statistical test (log-rank test) to see if survival differs significantly between the two groups.

```r
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = brca_tp53)
ggsurvplot(fit,
           data = brca_tp53,
           pval = T,
           risk.table = T)
fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = brca_tp53)
```

Here's what we can interpret out of the plot:
Red (HIGH TP53 expression): Patients with higher TP53 expression. Their curve drops faster, meaning a higher proportion died sooner.
Blue (LOW TP53 expression): Patients with lower TP53 expression. Their survival curve stays higher, meaning they lived longer overall.
Y-axis = probability of being alive.
At the beginning (time = 0), both groups start at 1.0 (100%).
Over time, the HIGH group falls to \~0.5, while the LOW group stays above \~0.75.
X-axis = follow-up time in days.
You see more separation after \~1000 days, with HIGH dropping more steeply.

---

## Cox Proportional Hazards Regression

Cox Proportional Hazards Regression is another statistical method that we can use in our analysis. It models the hazard rate (risk of death/event at a given time). Unlike Kaplan–Meier (which only compares groups like HIGH vs LOW), Cox regression lets you test how strongly a variable (like TP53 expression, age, stage, etc.) influences survival.

```r
cox_fit <- coxph(Surv(overall_survival, deceased) ~ strata, data = brca_tp53)
summary(cox_fit)

ggforest(cox_fit, data = brca_tp53)
```

OUTPUT:

```
Call:
coxph(formula = Surv(overall_survival, deceased) ~ strata, data = brca_tp53)

  n= 100, number of events= 12 

             coef exp(coef) se(coef)      z Pr(>|z|)  
strataLOW -1.6424    0.1935   0.7833 -2.097    0.036 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

          exp(coef) exp(-coef) lower .95 upper .95
strataLOW    0.1935      5.167   0.04169    0.8984

Concordance= 0.697  (se = 0.057 )
Likelihood ratio test= 5.8  on 1 df,   p=0.02
Wald test            = 4.4  on 1 df,   p=0.04
Score (logrank) test = 5.46  on 1 df,   p=0.02
```

According to the results, patients in the LOW TP53 group have about 19% of the risk of dying compared to the HIGH group. Put differently, their hazard of death is \~80% lower.

```