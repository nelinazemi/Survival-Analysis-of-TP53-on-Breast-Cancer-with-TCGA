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

clinical_brca <- GDCquery_clinic("TCGA-BRCA")

any(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_brca[,c(10,47,55)]

clinical_brca$deceased <- ifelse(clinical_brca$vital_status == "Alive", FALSE, TRUE)

table(clinical_brca$deceased)

clinical_brca$overall_survival <- ifelse(clinical_brca$vital_status == "Alive",
                                         clinical_brca$days_to_last_follow_up,
                                         clinical_brca$days_to_death)

table(clinical_brca$overall_survival)

tumoric_rna_query = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")

tumoric_rna_brca <- getResults(tumoric_rna_query)

tumoric_cases <- tumoric_rna_brca$cases[1:100]

normal_rna_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumoric_cases)

GDCdownload(normal_rna_query)
tcga_brca_data <- GDCprepare(normal_rna_query, summarizedExperiment = TRUE)
tcga_brca_matrix <- assay(tcga_brca_data, "unstranded")
tcga_brca_matrix[1:5,1:5]

gene_metadata <- as.data.frame(rowData(tcga_brca_data))
cases <- as.data.frame(colData(tcga_brca_data))

dds <- DESeqDataSetFromMatrix(countData = tcga_brca_matrix,
                              colData = cases,
                              design = ~ 1)

dds <- dds[rowSums(counts(dds)) >= 10]

vsd <- vst(dds, blind=FALSE)
brca_matrix_vst <- assay(vsd)
brca_matrix_vst[1:10,1:10]

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

median_value <- median(brca_tp53$counts)
brca_tp53$strata <- ifelse(brca_tp53$counts >= median_value, "HIGH", "LOW")

brca_tp53$case_id <- gsub('-01.*', '', brca_tp53$case_id)
brca_tp53 <- merge(brca_tp53, clinical_brca, by.x = 'case_id', by.y = 'submitter_id')

fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = brca_tp53)

ggsurvplot(fit,
           data = brca_tp53,
           pval = T,
           risk.table = T)

fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = brca_tp53)

cox_fit <- coxph(Surv(overall_survival, deceased) ~ strata, data = brca_tp53)
summary(cox_fit)

ggforest(cox_fit, data = brca_tp53)

cox_fit2 <- coxph(Surv(overall_survival, deceased) ~ strata + age_at_index, 
                  data = brca_tp53)
summary(cox_fit2)