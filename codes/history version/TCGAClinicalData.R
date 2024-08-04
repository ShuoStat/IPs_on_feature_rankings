
#-------------------------------------------------------------------------------
#- TCGA
#-------------------------------------------------------------------------------

library(Biobase)      #- version, 2.56.0
library(TCGAbiolinks) #- version, 2.24.3
library(SummarizedExperiment) #- version, 1.26.0
library(dplyr)        #- version, 1.0.9

#- functions
#
# trans_ordinal <- function(x, ord = list(), name = NULL){
#
#     trans <- function(x){
#         x <- gsub("-", "_", x)
#         gsub(" ", "", x)
#     }
#
#     x <- trans(x)
#     ord <- lapply(ord, trans)
#
#     n <- length(ord)
#     re <- mapply(function(i){
#         ifelse(x %in% unlist(ord[i:n]), 1, 0)
#     }, i = 2:n, SIMPLIFY = T)
#
#     na <- setdiff(unique(x), unlist(ord))
#     re[x %in% na,] <- NA
#
#     nam <- lapply(ord, function(x) gsub("\\s+", "_", paste0(x, collapse = "_")))
#     if (!is.null(name))
#         colnames(re) <- paste0("clin_", name, "_", unlist(nam)[-1])
#
#     if (is.null(name))
#         colnames(re) <- paste0("clin_", unlist(nam)[-1])
#
#     re
# }
#
# #- generate dummy variables
# trans_dummy <- function(x, grps = list(), name = NULL){
#
#     #- replace "-" to "_"
#
#     trans <- function(x){
#         x <- gsub("-", "_", x)
#         gsub(" ", "", x)
#     }
#
#     x <- trans(x)
#     grps <- lapply(grps, trans)
#
#     n <- length(grps)
#     re <- mapply(function(i){
#         ifelse(x %in% unlist(grps[i]), 1, 0)
#     }, i = 1:n, SIMPLIFY = T)
#
#     na <- setdiff(unique(x), unlist(grps))
#     re[x %in% na,] <- NA
#
#     nam <- lapply(grps, function(x) gsub("\\s+", "_", paste0(x, collapse = "_")))
#     if (!is.null(name))
#         colnames(re) <- paste0("clin_", name, "_", unlist(nam))
#
#     if (is.null(name))
#         colnames(re) <- paste0("clin_", unlist(nam))
#
#     re
# }
#
#
# getos <- function(...){
#
#     input <- list(...)
#     re <- mapply(function(x){
#         x %>%
#             dplyr::select(bcr_patient_barcode,
#                           vital_status,
#                           last_contact_days_to,
#                           death_days_to) %>%
#             mutate(ID = bcr_patient_barcode,
#                    status = recode(vital_status, "Alive" = 0, "Dead" = 1),
#                    time = as.numeric(ifelse(status == 0,
#                                             last_contact_days_to,
#                                             death_days_to))) %>%
#             dplyr::select(ID, status, time)
#     }, x = input, SIMPLIFY = F)
#
#     Reduce(rbind, re) %>%na.omit %>%
#         arrange(desc(time)) %>%
#         distinct(ID, .keep_all = T) %>%
#         filter(time > 0)
# }
#
# filterMrna <- function(rnaseq){
#
#     X <- assays(rnaseq)[["tpm_unstrand"]]
#     cNames <- colnames(X)
#     sel    <- substr(cNames, 14, 15) == "01"
#     colnames(X) <- substr(cNames, 1, 12)
#     X <- X[,sel]
#
#     #- filter genes
#
#     Z <- assays(rnaseq)[["unstranded"]]
#     Z <- Z[,sel]
#     genes.in <- edgeR::filterByExpr(Z)
#
#     #- filter in X
#     X[genes.in, ]
# }

#-- Clinical -------------------------------------------------------------------
#- BLCA
#-------------------------------------------------------------------------------
#
# s <- "BLCA"
# path <- paste0("D:/backup/TCGA/", s, "/clin/")
# dir.create(path, recursive = T)
#
# query <- GDCquery(project = paste0("TCGA-", s),
#                   data.category = "Clinical",
#                   data.type = "Clinical Supplement",
#                   data.format = "BCR Biotab")
#
# GDCdownload(query, directory = path)
# dat <- GDCprepare(query, directory = path)
#
# #- clinical information
# #- must included: age gender, cancer stage
# #- alternative height and weights
#
# tmp <- as.data.frame(dat$clinical_patient_blca[-c(1, 2),])
#
# apply(tmp, 2, function(x) {
#     if(length(table(x)) < 20)
#         table(x)
# })
#
# #- clinical
# clin <- dat$clinical_patient_blca %>%
#     as.data.frame() %>%
#     dplyr::slice(-c(1:2)) %>%
#     dplyr::select(bcr_patient_barcode,
#                   age_at_diagnosis, gender,
#                   ajcc_pathologic_tumor_stage,
#                   # new_tumor_event_dx_indicator,
#                   ajcc_metastasis_pathologic_pm,
#                   height_cm_at_diagnosis,
#                   weight_kg_at_diagnosis,
#                   history_other_malignancy,
#                   noninvasive_bladder_history,
#                   # tumor_status,
#                   tobacco_smoking_history_indicator,
#                   histologic_subtype,
#                   incidental_prostate_cancer_indicator) %>%
#
#     mutate(ID = bcr_patient_barcode,
#            clin_age = as.numeric(age_at_diagnosis),
#            clin_gender = recode(gender, "MALE" = 0, "FEMALE" = 1),
#            height = as.numeric(height_cm_at_diagnosis),
#            weight = as.numeric(weight_kg_at_diagnosis),
#            # clin_tumor_status = recode(tumor_status, "TUMOR FREE" = 0, "WITH TUMOR" = 1),
#            clin_BMI = ifelse(is.na(weight) | is.na(height), NA, weight / height^2 * 10000),
#            clin_history_other_malignancy = recode(history_other_malignancy, "No" = 0, "Yes" = 1),
#            clin_noninvasive_bladder_history =
#                recode(noninvasive_bladder_history, "NO" = 0, "YES" = 1),
#            clin_histologic_subtype = recode(histologic_subtype, "Non-Papillary" = 0,
#                                             "Papillary" = 1),
#            clin_incidental_prostate_cancer_indicator = recode(incidental_prostate_cancer_indicator,
#                                                               "NO" = 0, "YES" = 1)
#            # clin_new_tumor_event_dx_indicator = recode(new_tumor_event_dx_indicator, "NO" = 0,
#            #                                       "YES" = 1)
#     ) %>%
#     bind_cols(trans_ordinal(dat$clinical_patient_blca$ajcc_pathologic_tumor_stage[-c(1:2)],
#                             list(c("Stage I", "Stage II"),
#                                  c("Stage III"),
#                                  c("Stage IV"))),
#               trans_dummy(dat$clinical_patient_blca$tobacco_smoking_history_indicator[-c(1:2)],
#                           list(c("1"), c("2"), c("3"), c("4")),
#                           name = "tobacco_smoking_history_indicator")) %>%
#
#     dplyr::select(ID, starts_with("clin"))
#
# #- follow-up
# os <- getos(dat$clinical_patient_blca,
#             dat$clinical_follow_up_v2.0_blca,
#             dat$clinical_follow_up_v4.0_blca)
# clin <- inner_join(os, clin, by = "ID")
#
# #- missing data imputation
#
# library(mice)
# mis <- clin %>%
#     dplyr::select(starts_with("clin")) %>%
#     mice(m = 5, seed = 1) %>%
#     complete(1)
#
# clin[,colnames(mis)] <- mis
#
# #- RNA ---
#
# path <- paste0("D:/backup/TCGA/", s, "/mrna/")
# dir.create(path, recursive = T)
# query <- GDCquery(project = paste0("TCGA-", s),
#                   data.category = "Transcriptome Profiling",
#                   data.type = "Gene Expression Quantification",
#                   workflow.type = "STAR - Counts",
#                   legacy = F)
# GDCdownload(query, directory = path)
# rnaseq <- GDCprepare(query, directory = path)
# # save(list = "rnaseq", file = paste0(path, "mRNA.RData"))
#
# load(paste0(path, "mRNA.RData"))
# X <- filterMrna(rnaseq)
#
# #- match clinical and moleculars
# ind <- intersect(clin$ID, colnames(X))
# X   <- X[,ind]
# clin <- clin[match(ind, clin$ID),]
# clin$ID == colnames(X)
#
# X <- cbind(clin, t(X))
# setwd(path)
# save(list = "X", file = paste0(paste0("../", s, ".RData")))

#-------------------------------------------------------------------------------
#- COAD - Not included
#-------------------------------------------------------------------------------

s <- "COAD"
path <- paste0("D:/backup/TCGA/", s, "/clin/")
dir.create(path, recursive = T)

query <- GDCquery(project = paste0("TCGA-", s),
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query, directory = path)
dat <- GDCprepare(query, directory = path)
dat <- dat$clinical_patient_coad[-c(1, 2),]

# load used data
load(paste0("../data/TCGA-COAD.RData"))
usedCases <- colnames(X)
otherID <- substr(usedCases, 1, 12)[-73]
ipsID   <- substr(usedCases, 1, 12)[73]

library(dplyr)

dat <- dat %>%
    filter(bcr_patient_barcode == ipsID) %>%
    bind_rows(filter(dat, bcr_patient_barcode %in% otherID))

write.csv(dat, "../output/clin_coad.csv")

#-------------------------------------------------------------------------------

s <- "LUSC"
path <- paste0("D:/backup/TCGA/", s, "/clin/")
dir.create(path, recursive = T)

query <- GDCquery(project = paste0("TCGA-", s),
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query, directory = path)
dat <- GDCprepare(query, directory = path)
dat <- dat$clinical_patient_lusc[-c(1, 2),]

# load used data
load(paste0("../data/TCGA-LUSC.RData"))
usedCases <- colnames(X)
otherID <- substr(usedCases, 1, 12)[-51]
ipsID   <- substr(usedCases, 1, 12)[51]

library(dplyr)

dat <- dat %>%
    filter(bcr_patient_barcode == ipsID) %>%
    bind_rows(filter(dat, bcr_patient_barcode %in% otherID))

write.csv(dat, "../output/clin_lusc.csv")
#-------------------------------------------------------------------------------

s <- "HNSC"
path <- paste0("D:/backup/TCGA/", s, "/clin/")
dir.create(path, recursive = T)

query <- GDCquery(project = paste0("TCGA-", s),
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query, directory = path)
dat <- GDCprepare(query, directory = path)
dat <- dat$clinical_patient_hnsc[-c(1, 2),]

# load used data
load(paste0("../data/TCGA-HNSC.RData"))
usedCases <- colnames(X)
otherID <- substr(usedCases, 1, 12)[-93]
ipsID   <- substr(usedCases, 1, 12)[93]

library(dplyr)

dat <- dat %>%
    filter(bcr_patient_barcode == ipsID) %>%
    bind_rows(filter(dat, bcr_patient_barcode %in% otherID))

write.csv(dat, "../output/clin_hnsc.csv")

#-------------------------------------------------------------------------------

s <- "KIRC"
path <- paste0("D:/backup/TCGA/", s, "/clin/")
dir.create(path, recursive = T)

query <- GDCquery(project = paste0("TCGA-", s),
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query, directory = path)
dat <- GDCprepare(query, directory = path)
dat <- dat$clinical_patient_kirc[-c(1, 2),]

# load used data
load(paste0("../data/TCGA-KIRC.RData"))
usedCases <- colnames(X)
otherID <- substr(usedCases, 1, 12)[-133]
ipsID   <- substr(usedCases, 1, 12)[133]

library(dplyr)

dat <- dat %>%
    filter(bcr_patient_barcode == ipsID) %>%
    bind_rows(filter(dat, bcr_patient_barcode %in% otherID))

write.csv(dat, "../output/clin_kirc.csv")

#-------------------------------------------------------------------------------

s <- "PRAD"
path <- paste0("D:/backup/TCGA/", s, "/clin/")
dir.create(path, recursive = T)

query <- GDCquery(project = paste0("TCGA-", s),
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query, directory = path)
dat <- GDCprepare(query, directory = path)
dat <- dat$clinical_patient_prad[-c(1, 2),]

# load used data
load(paste0("../data/TCGA-PRAD.RData"))
usedCases <- colnames(X)
otherID <- substr(usedCases, 1, 12)[-107]
ipsID   <- substr(usedCases, 1, 12)[107]

library(dplyr)

dat <- dat %>%
    filter(bcr_patient_barcode == ipsID) %>%
    bind_rows(filter(dat, bcr_patient_barcode %in% otherID))

write.csv(dat, "../output/clin_prad.csv")


#-------------------------------------------------------------------------------

s <- "THCA"
path <- paste0("D:/backup/TCGA/", s, "/clin/")
dir.create(path, recursive = T)

query <- GDCquery(project = paste0("TCGA-", s),
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query, directory = path)
dat <- GDCprepare(query, directory = path)
dat <- dat$clinical_patient_thca[-c(1, 2),]

# load used data
load(paste0("../data/TCGA-THCA.RData"))
usedCases <- colnames(X)
otherID <- substr(usedCases, 1, 12)[-134]
ipsID   <- substr(usedCases, 1, 12)[134]

library(dplyr)

dat <- dat %>%
    filter(bcr_patient_barcode == ipsID) %>%
    bind_rows(filter(dat, bcr_patient_barcode %in% otherID))

write.csv(dat, "../output/clin_thca.csv")

#-------------------------------------------------------------------------------

s <- "STAD"
path <- paste0("D:/backup/TCGA/", s, "/clin/")
dir.create(path, recursive = T)

query <- GDCquery(project = paste0("TCGA-", s),
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query, directory = path)
dat <- GDCprepare(query, directory = path)
dat <- dat$clinical_patient_stad[-c(1, 2),]

# load used data
load(paste0("../data/TCGA-STAD.RData"))
usedCases <- colnames(X)
otherID <- substr(usedCases, 1, 12)[-22]
ipsID   <- substr(usedCases, 1, 12)[22]

library(dplyr)

dat <- dat %>%
    filter(bcr_patient_barcode == ipsID) %>%
    bind_rows(filter(dat, bcr_patient_barcode %in% otherID))

write.csv(dat, "../output/clin_stad.csv")

#-------------------------------------------------------------------------------

s <- "BRCA"
path <- paste0("D:/backup/TCGA/", s, "/clin/")
dir.create(path, recursive = T)

query <- GDCquery(project = paste0("TCGA-", s),
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query, directory = path)
dat <- GDCprepare(query, directory = path)
dat <- dat$clinical_patient_brca[-c(1, 2),]

# load used data
load(paste0("../data/TCGA-BRCA.RData"))
usedCases <- colnames(X)
otherID <- substr(usedCases, 1, 12)[-338]
ipsID   <- substr(usedCases, 1, 12)[338]

library(dplyr)

dat <- dat %>%
    filter(bcr_patient_barcode == ipsID) %>%
    bind_rows(filter(dat, bcr_patient_barcode %in% otherID))

write.csv(dat, "../output/clin_brca.csv")

#-------------------------------------------------------------------------------

load(paste0("../data/TCGA-LUSC.RData"))
icf <- EPIC::EPIC(X) # icf, immune cell fraction
cf  <- icf$cellFractions # cf, cell fraction
mp  <- icf$mRNAProportions

# percentage stacked barplot
# objD <- mp %>%
#     as.data.frame() %>%
#     mutate(grp = y,
#            grp = ifelse(row_number() == 51, "00", grp)) %>%
#     tidyr::pivot_longer(cols = 1:8) %>%
#     group_by(grp, name) %>%
#     summarise(value = round(mean(value), 3))
#
# ggplot(objD, aes(fill = name, y = value, x = grp)) +
#     geom_bar(position="fill", stat="identity") +
#     theme_bw()


histPurity <- function(datName, obs) {

    require(dplyr)
    require(ggplot2)

    load(paste0("../data/TCGA-", datName, ".RData"))
    # icf <- EPIC::EPIC(X) # icf, immune cell fraction
    # cf  <- icf$cellFractions # cf, cell fraction
    # # mp  <- icf$mRNAProportions

    useX <- X %>%
        as.data.frame() %>%
        tibble::rownames_to_column()
    re <- tidyestimate::estimate_score(useX , is_affymetrix = FALSE)

    objD <- re %>%
        as.data.frame() %>%
        mutate(grp = y) %>%
        filter(row_number() != obs) %>%
        dplyr::select(purity, grp) %>%
        na.omit()

    ibsValue <- re %>%
        as.data.frame() %>%
        mutate(grp = y) %>%
        filter(row_number() == obs) %>%
        dplyr::select(purity) %>%
        pull()


    # objD <- cf %>%
    #     as.data.frame() %>%
    #     mutate(grp = y) %>%
    #     filter(row_number() != obs) %>%
    #     select(otherCells, grp)
    #
    # ibsValue <- cf %>%
    #     as.data.frame() %>%
    #     filter(row_number() == obs) %>%
    #     select(otherCells) %>%
    #     pull()


    ggplot(objD, aes(purity , fill = grp)) +
        geom_histogram(alpha = 0.4, aes(y = ..density..), position = 'identity') +
        scale_fill_manual(values = c("11" = "darkcyan", "01" = "darkgoldenrod1")) +
        theme_bw() +
        theme(legend.position = "none",
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()) +
        xlab("Estimated Purity") +
        geom_vline(xintercept = ibsValue,
                   linetype = "dashed",
                   color = ifelse(y[obs] == "01", "darkgoldenrod1", "darkcyan"),
                   size = 1) +
        ggtitle(datName)
}

p1 <- histPurity("LUSC", 51)
p2 <- histPurity("STAD", 22)

p3 <- histPurity("BRCA", 338)
p4 <- histPurity("COAD", 73)

p5 <- histPurity("HNSC", 93)
p6 <- histPurity("KIRC", 133)

p7 <- histPurity("PRAD", 107)
p8 <- histPurity("THCA", 134)

p <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, p7, p8,
               ncol = 2, nrow = 4)

ggsave("../output/FigureHist.tiff", plot = p, width = 5, height = 7,
       compression = "lzw",
       units = "in",
       dpi = 600)

# library(IOBR)
# load(paste0("../data/TCGA-BRCA.RData"))
# D <- CIBERSORT(sig_matrix = lm22,
#                mixture_file =as.data.frame(X),
#                perm = 200,
#                QN = TRUE,
#                absolute = FALSE,
#                abs_method = "sig.score")
#


purityPlot <- function(datName, obs) {

    load(paste0("../data/TCGA-", datName, ".RData"))
    sampleID <- substr(colnames(X), 1, 16)
    ipsID <- sampleID[obs]

    re <- Tumor.purity[na.omit(match(sampleID, Tumor.purity$Sample.ID)), ]

    hisPlot <- function(method) {

        dat <- re[, c("Sample.ID", method)] %>%
            rename_with(~ c("ID", "value")) %>%
            mutate(value = as.numeric(gsub(",", ".", as.character(value))))

        ggplot(dat, aes(value)) +
            geom_histogram(alpha = 0.5, aes(y = ..density..),
                           position = 'identity',
                           fill = "darkcyan",
                           color = "grey") +
            theme_bw() +
            theme(legend.position = "none",
                  axis.title.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_blank()) +
            xlab(method) +
            geom_vline(xintercept = dat[dat$ID == ipsID, "value"],
                       linetype = "dashed",
                       color ="darkgoldenrod1",
                       size = 1) +
            scale_x_continuous(limits =  c(0, 1)) +
            ggtitle(datName)
    }

    p1 <- hisPlot("ESTIMATE")
    p2 <- hisPlot("ABSOLUTE") + ggtitle("")
    p3 <- hisPlot("LUMP") + ggtitle("")
    p4 <- hisPlot("IHC") + ggtitle("")

    ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1)

}

p1 <- purityPlot("LUSC", obs = 51)
p2 <- purityPlot("COAD", obs = 73)
p3 <- purityPlot("KIRC", obs = 133)

p <- ggpubr::ggarrange(p1, p2, p3, ncol = 1, nrow = 3)

ggsave("../output/FigurePurity.tiff", plot = p, width = 8, height = 6,
       compression = "lzw",
       units = "in",
       dpi = 600)

