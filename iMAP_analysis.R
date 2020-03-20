# iMAP Analysis // Bosco TADDEI

# This script is analyzing cytokine and cytometric data
# collected from patient treated for lupus.



# Packages acquisition --------------------------------------------------------------
library("plyr")
library("dplyr")
library("stringr")
library("factoextra")
library("gridExtra")
library("ggrepel")
library("reshape2")


# Loading data ----------------------------------------------------------------------
setwd("./..")
getwd()

# Get population --------------------------------------------------------------------
df_pop <- read.table("lupil2_populations.txt", 
                     header=TRUE, sep="\t") %>%
          dplyr::rename(Subject = SUBJID) %>%
          arrange(Subject)
df_pop["Subj_id"] <- rownames(df_pop)
df_pop["Establishment"] <- substr(df_pop$Subject,1,3) %>%
                           as.factor()
summary(df_pop)

# Extract paper subjects
paper_subjects <- df_pop[df_pop[,"Paper"] == T,"Subject"]


# Get Cytometrie Panels data --------------------------------------------------------
# Initialize the df
df_cytometrie <- data.frame(paper_subjects) %>%
                 dplyr::rename(Subject = paper_subjects)



# Collect data from the different panels
panel_ids <- c("02","03","08","09","10","12","13")
for (id in panel_ids){
  file_panel = paste0("lupilCytometrie/Panel_", id, ".txt")
  panel <- read.table(file_panel,
                      header=TRUE, sep="\t") %>%
           dplyr::rename(Subject = Sample.name)
  df_cytometrie  <- join(df_cytometrie, panel)
}

df_cytometrie <- arrange(df_cytometrie, Subject, Visit)
summary(df_cytometrie)

# Names of the cytometrie numerical variables 
varnum_cytometrie <- df_cytometrie %>%
                     select(-c("Subject","Visit")) %>%
                     colnames()
                
# Long cytometrie data frame      
lg_cytometrie <- melt(data = df_cytometrie, id.vars = c("Subject", "Visit")) %>%
                 join(df_pop[,c("Subject","ARM","Responder")]) %>%
                 arrange(Subject,Visit)
lg_cytometrie["Panel"] <- lg_cytometrie$variable %>%
                          as.vector() %>%
                          strsplit(split = "_P") %>%
                          sapply("[[",1)

# Get Cytokine data ------------------------------------------------------------------
file_cytokine <- "Lupil2_cytokines_191105/lupil2_conc_log2impnorm_20191022.csv"
df_cytokine <- read.table(file_cytokine,
                          header=TRUE, sep=",")

# Extracting Subject id and visit id
id <- as.vector(df_cytokine[,1])
id <- strsplit(id,"V")
subject <- sapply(id, "[[", 1)
subject <- str_pad(subject, 9, "left", "0") 
subject <- paste0(substr(subject,1,3),"-",substr(subject,4,6),"-",substr(subject,7,9))
visit   <- sapply(id, "[[", 2)
visit   <- str_pad(visit, 2, "left", "0")
visit   <- paste0("V",visit)
df_cytokine["Subject"] <- as.factor(subject)
df_cytokine["Visit"] <- as.factor(visit)

df_cytokine <- df_cytokine[-1] %>%
               subset(Subject %in% paper_subjects) %>%
               arrange(Subject, Visit)

# Name of the cytokine numerical variables
varnum_cytokine <- df_cytokine %>%
                   select(-c("Subject","Visit")) %>%
                   colnames()


# Remove the Outlier
df_cytokine <- df_cytokine[,"Subj_id" != 90]

summary(df_cytokine)





# Long cytokine data frame
lg_cytokine <- melt(data = df_cytokine, id.vars = c("Subject", "Visit")) %>%
               join(df_pop[,c("Subject","ARM","Responder")]) %>%
               arrange(Subject,Visit)



# Analyzing data ---------------------------------------------------------------------

# PCA --------------------------------------------------------------------------------

# Function to plot infos related to pca of the previously constructed df
pca_iMAP <- function(df,category = "Responder"){
  # Compute pca
  pca<- select(df, - c("Subject","Visit")) %>%
        na.omit() %>%
        prcomp(scale. = TRUE)
  # Info on whether a subject is responder or not
  groups <- na.omit(df) %>%
            join(df_pop) %>%
            select(category)
  # Create the different plots
  eig <- fviz_eig(pca)
  ind <- fviz_pca_ind(pca,
                      col.ind = groups[,1])
  var <- fviz_pca_var(pca)
                      # select.var = list(c(contrib = 1,"IL.2")))
  # Organize the different plots
  plots <- list(eig,ind,var)
  lay   <- rbind(c(1,1),
                 c(2,3))
  grid.arrange(grobs = plots,
               layout_matrix = lay)
  
  return(pca)
}

# More PCA info
main_contrib <- function(df, pca, n, add_il2 = FALSE){
  contribs <- get_pca_var(pca_cytokine)$contrib[,1] %>%
    sort(decreasing = TRUE) %>%
    names()
  contribs <- contribs[1:n]
  if (add_il2){
    contribs <- append(contrib_dim1,"IL.2")
  }
  data = df[,contribs]
  if (length(contribs) == 1){
    hist(data,
         xlab = contribs)
  } else {
    plot(data)
  }
}


# Cytometrie PCA
pca_iMAP(df_cytometrie)

# Cytokine PCA
pca_cytokine <- pca_iMAP(df_cytokine,"Establishment")

main_contrib(df_cytokine,pca_cytokine,0,TRUE)

contrib_dim1 <- get_pca_var(pca_cytokine)$contrib[,1] %>%
                sort(decreasing = TRUE) %>%
                names()
contrib_dim1 <- contrib_dim1[1:3]

plot(df_cytokine[,c(contrib_dim1[1:3],"IL.2")])
var_1 = contrib_dim1[1]
t.test(df_cytokine[,var_1])   
hist(df_cytokine[,var_1],
     xlab = var_1)

# Per visit pca
# visits <- levels(df_cytokine$Visit)
visits <- "V05"
for (visit in visits){
  filter <- df_cytokine$Visit == visit
  data <- df_cytokine[filter,]
  row.names(data) <- join(data,df_pop)[,"Subj_id"]
  pca_iMAP(data)
  contrib_dim1 <- get_pca_var(pca_cytokine)$contrib[,1] %>%
                  sort(decreasing = TRUE)
}


# Analyse monovariée
# Plot Cytokine data distribution
vars <- select(df_cytokine, -c("Subject","Visit")) %>%
        colnames() %>%
        reorder(contrib_dim1)
for (var in vars){
  data <- df_cytokine[, var]
  png(paste0("./plots/hist_",var,"_all.png"))
  hist(data,
       main = paste("Histogram of", var, "(all)"),
       xlab = var,
       col = "cadetblue3")
  dev.off()
}

# Plot Means
mean <- select(df_cytokine,vars) %>%
        colMeans(na.rm = TRUE)
sem  <- select(df_cytokine,vars) %>%
        apply(2,sd,na.rm = TRUE)
p <- join(df_cytokine,df_pop)[,c(vars,"Responder")] %>%
     ggplot() +
     geom_bar(aes(x=vars, y = mean), 
              stat = "identity",
              fill = "cadetblue3") +
     geom_errorbar((aes(x=vars, ymin=mean - sem, ymax = mean + sem))) +
     facet_grid( ~ Responder) +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     ggtitle("Cytokine data (Mean and SEM)")
p
