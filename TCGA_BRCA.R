# source("https://bioconductor.org/biocLite.R")
# biocLite("TCGAbiolinks")
rm(list=ls())
library(TCGAbiolinks)
clin.query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical")
json<-tryCatch(GDCdownload(clin.query), 
                  error = function(e) GDCdownload(clin.query, method = "client"))
clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
clinical.patient.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up")
clinical.drug <- GDCprepare_clinic(clin.query, clinical.info = "drug")
clinical.admin <- GDCprepare_clinic(clin.query, clinical.info = "admin")
clinical.radiation <- GDCprepare_clinic(clin.query, clinical.info = "radiation")
clinical.stage_event <- GDCprepare_clinic(clin.query, clinical.info = "stage_event")
clinical.new_tumor_event <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
clinical.all <- GDCprepare_clinic(clin.query)


GBM_path_subtypes <- TCGAquery_subtype(tumor = "brca")






clin.query <- GDCquery(project = "TCGA-BLCA", data.category = "Clinical", barcode = "TCGA-FD-A5C0")
json  <- tryCatch(GDCdownload(clin.query), error = function(e) GDCdownload(clin.query, method = "client"))

z=TCGAbiolinks:::getGDCprojects()$project_id
query <- GDCquery(project = "TCGA-STAD", 
                  data.category = "Clinical")
GDCdownload(query)
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")

clinical <- GDCprepare_clinic(query, clinical.info = "patient")

LGG_path_subtypes <- TCGAquery_subtype(tumor = "brca")