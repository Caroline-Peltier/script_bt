# iMAP Analysis // Bosco TADDEI

# This script is analyzing cytokine and cytometric data
# collected from patient treated for lupus.


# General Set-up Functions -----------------------------------------------------------

# Packages acquisition 
load_pkg <- function(){
  library("plyr")
  library("dplyr")
  library("stringr")
  library("ggplot2")
  library("factoextra")
  library("gridExtra")
  library("ggrepel")
  library("reshape2")
  library("RGCCA")
}


# Creating the different data frames -------------------------------------------------

# Get clinical data
get_clini <- function(){
 
   df <- read.table("lupil2_populations.txt", 
                    header=TRUE, sep="\t")    %>%
         dplyr::rename(Subject = SUBJID)      %>%
         arrange(Subject)
  
  # Add quality of life columns
  df["Subj_id"]    <- factor(row.names(df))
  df["Country_id"] <- substr(df$Subject,1,3) %>%
                      as.factor()
  
  return(df)
}

# Get Cytometrie Panels data 
get_cytometrie <- function(subjects = paper_subjects, panels = all_panels){
  
  # Initialize the df
  df <- data.frame(paper_subjects) %>%
        dplyr::rename(Subject = paper_subjects)
  
  # Collect data from the panels
  for (id in panel_ids){
    file_panel <- paste0("lupilCytometrie/Panel_", id, "_norm.txt")
    panel <- read.table(file_panel, header=TRUE, sep="\t") %>%
             dplyr::rename(Subject = Sample.name)
    
    df <- join(df, panel)
  }
  
  #Reorder the df regarding Subjects
  df <- arrange(df, Subject, Visit)
  
  return(df)
}

# Get Cytokine data
get_cytokine <- function(file = file_cytokine){
  
  df <- read.table(file, header=TRUE, sep=",")
  
  # Extracting Subject id and visit id
  id <- as.vector(df[,"X"])
  id <- strsplit(id,"V")
  subject <- sapply(id, "[[", 1)
  subject <- str_pad(subject, 9, "left", "0") 
  subject <- paste0(substr(subject,1,3),"-",substr(subject,4,6),"-",substr(subject,7,9))
  visit   <- sapply(id, "[[", 2)
  visit   <- str_pad(visit, 2, "left", "0")
  visit   <- paste0("V",visit)
  df["Subject"] <- as.factor(subject)
  df["Visit"]   <- as.factor(visit)
  
  df <- df[-1] %>%
        subset(Subject %in% paper_subjects) %>%
        arrange(Subject, Visit) %>%
        dplyr::select("Subject","Visit",everything())
  
  return(df)
}


# Other information extracting functions ---------------------------------------------

# Names of the cytometrie numerical variables 
get_varnum <- function(df){
  
  varnum <- df %>%
            select_if(is.numeric) %>%
            colnames()
  
  return(varnum)
}

get_crit <- function(dat){
  
  crit <- dat %>%
          select_if(is.factor) %>%
          colnames()
  
  return(crit)
}

# Long data frame      
make_long <- function(df,crits){
  
  factors <- df %>%
             select_if(is.factor) %>%
             colnames()
  
  lg_df <- df %>%
           melt(id.vars = factors) %>%
           na.omit() %>%
           join(df_clini[,c("Subject",crits)]) %>%
           arrange(Subject,Visit)
}

# Add panel info
manage_panel <- function(lg_df){
  
  lg_df["Panel"]    <- lg_df$variable %>%
                        substr(1,8) %>%
                        as.factor()
  lg_df["variable"] <- lg_df$variable %>%
                       substring(10,1000) %>%
                       as.factor()
  
  #Re-order cols
  lg_df <- lg_df %>%
           dplyr::select("Subject","Visit","Panel",everything())
  
  return(lg_df)
}

# Remove cols that have a HS equivalent
manage_HS <- function(df){
  
  varnum  <- get_varnum(df)
  badcols <- varnum[grepl('HS',varnum)] %>%
              substring(3)
  df[,c(badcols)] <- list(NULL)
  
  return(df)
}

#Remove cols that have var == 0
drop_nul_var <- function(df){
  
  badcols <- df[,which(apply(df,2,var) == 0)] %>%
             colnames()
  print(paste("Zero variance cols removed :",badcols))
  
  df[,badcols] <- list(NULL)
  
  return(df)
}


# Analysis Functions -----------------------------------------------------------------

# Prepare data 
get_data <- function(df, visit, crit){
  
  varnum <- get_varnum(df)
  
  dat <- df[df$Visit == visit,] %>%
         join(df_clini[,c("Subject", crit, "Subj_id")]) %>%
         dplyr::select(c(varnum, crit, "Subj_id"))
  row.names(dat) <- dat$Subj_id
  dat$Subj_id <- NULL
  
  return(dat)
}

#Compute the difference of 2 dat
dat_diff <- function(dat1, dat2){
  
  rows   <- intersect(rownames(dat1),rownames(dat2))
  varnum <- get_varnum(dat1)
  
  dat <- dat1[rows,]
  dat[,varnum] <- dat[,varnum] - dat2[rows,varnum]
  
  return(dat)
  
}

#PCA
plot_pca <- function(dat, cols, scale = TRUE){
  
  varnum <- get_varnum(dat)
  crit   <- get_crit(dat)
  
  pca <- dat[,varnum] %>%
         na.omit() %>%
         drop_nul_var() %>%
         prcomp(scale = scale)
  
  group <- dat %>%
           na.omit() %>%
           dplyr::select(crit)
  
  eig <- fviz_eig(pca)
  
  ind <- fviz_pca_ind(pca,
                      col.ind = group[,1]) +
         scale_color_manual(values = cols)
  
  var <- fviz_pca_var(pca,
                      select.var = list(contrib = 5),
                      repel = TRUE)
  print(eig)
  print(var)
  print(ind)
  
  
  return(pca)
}

# t.test
plot_t.test <- function(dat, pval.max, cols){
  
  varnum <- dat %>%
            drop_nul_var() %>%
            get_varnum()
  crit   <- get_crit(dat)
  y <- dat[,crit]
      
  i <- 0
  vars <- c()
  for (var in varnum){
    x = dat[,var]
    
    t <- t.test(x ~ y)
    pval <- t$p.value
    
    if (pval < pval.max){
      i <- i+1
      vars[[i]] <- var
      
      bp <- ggplot(data = dat, aes(x = y, y = x, fill = y)) + 
            geom_boxplot() +
            geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.5, color = "white", alpha = 0.4) +
            labs(title = var, y = "", 
                 x = paste("Student test pval =", pval)) +
            theme_classic() +
            scale_fill_manual(values = cols)
      print(bp)
    }
    
  }
  
  return(vars)
}

#Plot Wilcoxon test
plot_w.test <- function(dat, pval.max, cols){
  
  varnum <- dat %>%
            drop_nul_var() %>%
            get_varnum()
  crit   <- get_crit(dat)
  y <- dat[,crit]
  
  i <- 0
  vars <- c()
  for (var in varnum){
    x = dat[,var]
    
    t <- wilcox.test(x ~ y)
    pval <- t$p.value
    
    if (pval < pval.max){
      i <- i+1
      vars[[i]] <- var
      
      bp <- ggplot(data = dat, aes(x = y, y = x, fill = y)) + 
        geom_boxplot() +
        geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.5, color = "white", alpha = 0.4) +
        labs(title = var, y = "", 
             x = paste("Wilcoxon test pval =", pval)) +
        theme_classic() +
        scale_fill_manual(values = cols)
      print(bp)
    }
    
  }
  
  return(vars)
}

# PLS DA
plot_pls <- function(dat, cols = NULL, type = "both", compx = 1, n_mark = 10){
  
  varnum <- get_varnum(dat)
  crit   <- get_crit(dat)
  
  group <- dat[,crit]; names(group) <- rownames(dat)
  y <- data.frame(model.matrix( ~  group-1, data = group))
  blocks <- list(x = dat[,varnum], y = y)
  pls <-rgcca(block = blocks, ncomp = c(3,1), type = "pls")
  
  plot(pls, i_block=1, type=type, resp = group, compx = compx, compy = 2,n_mark = n_mark,colors = cols)
  
  return(pls)
}



