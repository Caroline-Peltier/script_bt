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

# Long cytokine data frame
lg_df <- make_long(df,c("Responder","ARM"))



# Load Cytometrie data ---------------------------------------------------------------

all_panels <- c("01","02","03","04","05","08","09","10","12")

df <- get_cytometrie()

# Long cytometrie data frame
lg_df <- make_long(df,c("Responder","ARM")) %>%
         manage_panel()


# Per time point analysis   -----------------------------------------------------------

# Set visit
visit    <- "V01"

# Get appropriate data
dat <- get_data(df, visit, crit)

# Split data into IL2 and pcb treated
dat_IL2 <- dat[dat$ARM == "ILT-101",]
dat_IL2$ARM <- NULL
dat_pcb <- dat[dat$ARM == "Placebo",]
dat_pcb$ARM <- NULL          


# Wilcoxon test and boxplot if pval significant
plot_w.test(dat_IL2, bonferroni = FALSE, cols)
plot_w.test(dat_pcb, bonferroni = FALSE, cols)


# PCA
pca_IL2 <- plot_pca(dat_IL2, cols)
pca_pcb <- plot_pca(dat_pcb, cols)


# PLS-DA
plot_pls(dat_IL2,cols)
plot_pls(dat_pcb,cols)



# Evolution from baseline ----------------------------------------------------------

# Get baseline data
dat0  <- get_data(df, "V01", crit)

# Set visit
visit <- "V05"

# Get appropriate data
dat  <- get_data(df, visit, crit) %>%
        dat_diff(dat0)

# Split data into IL2 and pcb treated
dat_IL2 <- dat[dat$ARM == "ILT-101",]
dat_IL2$ARM <- NULL
dat_pcb <- dat[dat$ARM == "Placebo",]
dat_pcb$ARM <- NULL 



# Wilcoxon test and boxplot if pval significant
plot_w.test(dat_IL2, bonferroni = FALSE, cols)
plot_w.test(dat_pcb, bonferroni = FALSE, cols)

# PCA
pca_IL2 <- plot_pca(dat_IL2, cols)
pca_pcb <- plot_pca(dat_pcb, cols)

# PLS-DA
plot_pls(dat_IL2,cols)
plot_pls(dat_pcb,cols)