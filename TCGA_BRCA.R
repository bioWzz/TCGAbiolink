#Processing the data from the TCGA using the TCGAbiolinks R package
#Install teh TCGAbiolinks package
source("http://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")

#load the TCGAbiolinks package
library(TCGAbiolinks)
clin.query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical")
json<-tryCatch(GDCdownload(clin.query),error = function(e) GDCdownload(clin.query, method = "client"))

clinical.patient <- GDCprepare_clinic(clin.query, clinical.info = "patient")
write.table(clinical.patient, file = "/Users/wzz/TCGA/brca_patient.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
clinical.patient.followup <- GDCprepare_clinic(clin.query, clinical.info = "follow_up")
clinical.drug <- GDCprepare_clinic(clin.query, clinical.info = "drug")
clinical.admin <- GDCprepare_clinic(clin.query, clinical.info = "admin")
clinical.radiation <- GDCprepare_clinic(clin.query, clinical.info = "radiation")
clinical.stage_event <- GDCprepare_clinic(clin.query, clinical.info = "stage_event")
clinical.new_tumor_event <- GDCprepare_clinic(clin.query, clinical.info = "new_tumor_event")
clinical.all <- GDCprepare_clinic(clin.query)


# Check with subtypes from TCGAprepare and update examples
BRCA_path_subtypes <- TCGAquery_subtype(tumor = "brca")
write.table(BRCA_path_subtypes, file = "/Users/wzz/TCGA/brca.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
#miRNA data query
miRNA.query <- GDCquery(project = "TCGA-BRCA", data.category ="Gene expression",data.type="miRNA gene quantification",legacy=TRUE)

# Mutation data processing
acc.muse.maf <- GDCquery_Maf("BRCA", pipelines = "muse")
acc.varscan2.maf <- GDCquery_Maf("BRCA", pipelines = "varscan2")
acc.somaticsniper.maf <- GDCquery_Maf("BRCA", pipelines = "somaticsniper")
acc.mutect.maf <- GDCquery_Maf("BRCA", pipelines = "mutect")

write.table(acc.muse.maf, file = "/Users/wzz/TCGA/brca_Mutation.txt", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
