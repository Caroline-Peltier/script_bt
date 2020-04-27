# Code for Responder  in IL2 treated patient // Bosco Taddei

# General Set Up --------------------------------------------------------------------

# Set Work directory
setwd("./..")
getwd()

# Load functions
source("script_bt/iMAP_analysis.R")

#Load Packages
load_pkg()

#General Parameters
crit  <- c("Responder","ARM")
cols  <- c("#FF3333", "#00CC00")

# Load Clinical data -----------------------------------------------------------------

df_clini <- get_clini()

# Extract paper subjects
paper_subjects <- df_clini[df_clini[,"Paper"] == T,"Subject"]



# Load Cytokine data -----------------------------------------------------------------

# Not Normalized file
#file_cytokine <- "Lupil2_cytokines_191105/lupil2_conc_log2imp_20191022.csv"

# Normalized file
file_cytokine <- "Lupil2_cytokines_191105/lupil2_conc_log2impnorm_20191022.csv"

df <- get_cytokine() %>%
      manage_HS()

varnum <- get_varnum(df)

# Long cytokine data frame
lg_df <- make_long(df,c("Responder","ARM"))



# Load Cytometrie data ---------------------------------------------------------------

#General Set up
all_panels <- c("01","02","03","04","05","08","09","10","12")

df <- get_cytometrie()

varnum <- get_varnum(df)

lg_df <- make_long(df,c("Responder","ARM")) %>%
         manage_panel()


# Per time point analysis -----------------------------------------------------------

visit    <- "V29"

# Get appropriate data
dat <- get_data(df, visit, crit)


dat_IL2 <- dat[dat$ARM == "ILT-101",]
dat_IL2$ARM <- NULL
dat_pcb <- dat[dat$ARM == "Placebo",]
dat_pcb$ARM <- NULL          


#t.test and boxplot if pval significant
plot_w.test(dat_IL2, bonferroni = FALSE, cols)
plot_w.test(dat_pcb, bonferroni = FALSE, cols)


#PCA
pca_IL2 <- plot_pca(dat_IL2, cols)
pca_pcb <- plot_pca(dat_pcb, cols)


# PLS-DA
plot_pls(dat_IL2,cols)
plot_pls(dat_pcb,cols)



# Evolution from baseline ----------------------------------------------------------

dat0 <- get_data(df, "V01", crit)
visit <- "V05"

dat  <- get_data(df, visit, crit) %>%
        dat_diff(dat0)

dat_IL2 <- dat[dat$ARM == "ILT-101",]
dat_IL2$ARM <- NULL
dat_pcb <- dat[dat$ARM == "Placebo",]
dat_pcb$ARM <- NULL 

#t.test and boxplot if pval < .05 / 62
plot_w.test(dat_IL2, bonferroni = FALSE, cols)
plot_w.test(dat_pcb, bonferroni = FALSE, cols)


#t.test and boxplot if pval < .05 / 62
plot_w.test(dat_IL2, pval.max, cols)

#PCA
pca_IL2 <- plot_pca(dat_IL2, cols)
pca_pcb <- plot_pca(dat_pcb, cols)

#PLS
plot_pls(dat_IL2, type = "cor", 2)
plot_pls(dat_pcb, type = "cor", 2)
