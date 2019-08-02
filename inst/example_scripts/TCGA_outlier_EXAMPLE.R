# library(TCGAbiolinks)
# 
# GDCquery()
# 
# temp1 <- getGDCprojects()
# temp2 <- TCGAbiolinks:::getProjectSummary("TCGA-BRCA")
# 
# query <- GDCquery(project = "TCGA-BRCA",
#                   data.category = "Gene expression",
#                   data.type = "Gene expression quantification",
#                   platform = "Illumina HiSeq",
#                   file.type  = "normalized_results",
#                   experimental.strategy = "RNA-Seq",
#                   legacy = TRUE)
# 
# GDCdownload(query)


library(GenomicDataCommons)

## Checks the status of GDC - needs to be OK
#status()

qfiles = files()
files(filter( type == "gene_expression") )
files() %>% filter( type == "gene_expression")
qfiles = files() %>% filter( ~type == 'gene_expression')

grep('pro',available_fields('files'),value=TRUE) %>% 
    head()

qfiles = files() %>%
    filter( ~cases.project.project_id == 'TCGA-OV' & type == 'gene_expression')


files(filter( cases.project.project_id == 'TCGA-OV' & type == 'gene_expression'))



pquery = projects()
str(pquery)
presults = results_all(pquery)


qcases = cases()
qcases1 = cases(available_fields("cases"))
qcases2 = cases() %>% GenomicDataCommons::select(available_fields('cases'))


qfiles = files()
files(type=='gene_expression')
qfiles = files() %>% filter( type == 'gene_expression')


qfiles = files() %>%
    filter( cases.project.project_id == 'TCGA-OV' & type == 'gene_expression')


resp = cases() %>% filter(~ project.project_id=='TCGA-BRCA' &
                              samples.sample_type=='Solid Tissue Normal') %>%
    GenomicDataCommons::select(c(default_fields(cases()),'samples.sample_type')) %>%
    response_all()
count(resp)

# the number of all the projects queried
pcount = count(pquery)
# JSON file that is presented as a nested list that has all the info on our query
presults = results(pquery)


default_fields('files')
default_fields(pquery)

