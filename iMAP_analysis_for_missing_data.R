# Loading libraries
#---------------------
library("plyr")
library("dplyr")
library("stringr")
library("factoextra")
library("gridExtra")
library("ggrepel")
library("reshape2")

# Sourcing external functions (you should be in script_bt)
#-------------------------------------------------------
source("functionForNAAnalysis.r")


# Getting dataframes 
#--------------------------------

df_pop <- read.table("../lupil2_populations.txt", 
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
#panel_ids <- c("02","03","08","09","10","12","13")

#---> Here I replaced by the new panels, I also changed the path in file_panel
panel_ids<- c("01_norm","02_norm","03_norm","04_norm","05_norm","08_norm","09_norm","10_norm","12_norm")
#----
for (id in panel_ids){
  file_panel = paste0("../lupilCytometrie/Panel_", id, ".txt")
  panel <- read.table(file_panel,
                      header=TRUE, sep="\t") %>%
           dplyr::rename(Subject = Sample.name)
  df_cytometrie  <- join(df_cytometrie, panel)
}

df_cytometrie <- arrange(df_cytometrie, Subject, Visit)
summary(df_cytometrie)

library(ggplot2)
p1 <- ggplot(df, aes(factor(year), Y, fill = browser)) + 
  geom_bar(stat = "identity", width = 1, size = 1, color = "white") +
  coord_polar("y") + 
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors)
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
file_cytokine <- "../Lupil2_cytokines_191105/lupil2_conc_log2impnorm_20191022.csv"
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
df_cytokine_wo90 <- df_cytokine[,"Subj_id" != 90]

summary(df_cytokine)

# Searching where the subject 010-001-033 is strange ?
for(i in 1:75) #loop on the attributes
{ 
  df_cytokine=as.data.frame(df_cytokine)
  #Calculation of the mean of the others for this attribute
  means=mean(df_cytokine[df_cytokine[,"Subject"]!="010-001-033" &df_cytokine[,"Visit"]=="V01" ,i],na.rm=T)
  # Calculation of the standard deviations of the others
  sds=sd(df_cytokine[df_cytokine[,"Subject"]!="010-001-033" &df_cytokine[,"Visit"]=="V01" ,i],na.rm=T)
  # Score of the subject
  score_suj=df_cytokine[df_cytokine[,"Subject"]=="010-001-033" &df_cytokine[,"Visit"]=="V01" ,][,i]
  if(abs(score_suj-means)>3*sds)
  {
    print(colnames(df_cytokine)[i])
  }
}  

# Long cytokine data frame
lg_cytokine <- melt(data = df_cytokine, id.vars = c("Subject", "Visit")) %>%
               join(df_pop[,c("Subject","ARM","Responder")]) %>%
               arrange(Subject,Visit)

# Analyzing data ---------------------------------------------------------------------
library(RGCCA)
visits=union(unique(df_cytometrie[,"Visit"]),unique(df_cytokine[,"Visit"]))
pop_2=df_pop[df_pop[,"Subject"]%in%paper_subjects,]
# Deux repetitions pour le patient 003-006-003 :lignes 8 et 9 /!\ plus vrai pour nouvelles données
dfList=list(cytometrie=df_cytometrie, cytokine=df_cytokine)
# Producing html for missing data analysis
selection=list(cytometrie=varnum_cytometrie,cytokine=varnum_cytokine)
produceHtmlMissingByBlock(dfList,selection)
produceHtmlMissingByVisit()



#----- Code from the net for multilayer pieplots
# using your data
df <- structure(list(Animal = structure(c(2L, 2L, 2L, 2L, 2L, 2L, 2L,
                                          2L, 2L, 2L, 2L, 1L, 1L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L), .Label = c("Buffalo",
# add special attribute "dist" using "->" as sep, this will increase the size of the terminal nodes to make space for the cheese names
df$Name <- paste(df$Name, "dist:3", sep="->")

df1 <- df %>% 
  mutate(Colour = ifelse(.$Animal == "Goat", "#CD9B1D", ifelse(.$Animal == "Sheep", "#EEC900", "#FFD700"))) %>% 
  mutate(Index = 1) %>% 
  group_by(Animal)
First <- ggplot(df1) + geom_bar(aes(x=1, y=Animal, fill=Animal, 
                                    label = Animal), position='stack', stat='identity', size=0.15) + theme(panel.grid = element_blank(), axis.title=element_blank(), 
        legend.position="none", axis.ticks=element_blank(), 
        axis.text = element_blank())

Second <- First + geom_bar(data=df1, aes(x=2, y=Animal, fill=Texture, group=Animal),
           position='stack', stat='identity', size=0.15, colour = "black") + scale_color_brewer(palette = "YlOrBr")+ scale_fill_brewer(palette = "YlOrBr") + theme(axis.title=element_blank(), legend.position="none",    
        axis.ticks=element_blank(), axis.text = element_blank())

Third <- Second + geom_bar(data=df1, aes(x=3, y=Animal, fill=Name), 
                           position='stack', stat='identity', size=0.15, colour = "black") + scale_fill_manual(values = c("#EEC900", "#FFD700", "#CD9B1D",
                               "#FFD700", "#DAA520", "#EEB422", "#FFC125", "#8B6914", "#EEC591", 
                               "#FFF8DC", "#EEDFCC", "#FFFAF0", "#EEC900", "#FFD700", "#CDAD00", 
                               "#FFF68F", "#FFEC8B", "#FAFAD2", "#FFFFE0", "#CD853F", "#EED8AE", 
                               "#F5DEB3", "#FFFFFF", "#FFFACD", "#D9D9D9", "#EE7600", "#FF7F00",
                               "#FFB90F", "#FFFFFF")) + theme(axis.title=element_blank(), legend.position="none",
        axis.ticks=element_blank(), axis.text.y = element_blank(), 
        panel.background = element_rect(fill = "black"))
Fourth <- Third + geom_bar(data=df1,aex(x=4,y))

Third + coord_polar('y')
#---- End  of the code

# Test for panel 2
file_panel = paste0("lupilCytometrie/Panel_", panel_ids[2], ".txt")
panel <- read.table(file_panel,
                    header=TRUE, sep="\t") %>%
  dplyr::rename(Subject = Sample.name)
newcols=sapply(colnames(panel)[3:length(colnames(panel))],function(x){return(substr(x,10,nchar(x)))})
panel2=panel
colnames(panel2)=c("Subject","Visit",newcols)
panelP=cbind(panel2[,c("Subject","Visit")],panel2[,substr(colnames(panel2),1,1)=="P"])
newcolsP=sapply(colnames(panel2)[3:length(colnames(panelP))],function(x){return(substr(x,3,nchar(x)))})
colnames(panelP)=c("Subject","Visit",newcolsP)
panelN=cbind(panel2[,c("Subject","Visit")],panel2[,substr(colnames(panel2),1,1)=="N"])
newcolsN=sapply(colnames(panelN)[3:length(colnames(panelN))],function(x){return(substr(x,3,nchar(x)))})
colnames(panelN)=c("Subject","Visit",newcolsN)

#
panelN[1,]
# 58/1374
#[1] 0.04221252

#7.379493/58.05049 
panelP[1,]

# PCA on panel 2

resPCA=PCA(cbind(merge(panelN,df_pop)[,"Responder"],panelN[,-c(1:2)]),quali.sup=1)
plot(resPCA,habillage=1)

resPCA=PCA(cbind(merge(panelN,df_pop)[,"ARM"],t(scale(t(panelN[,-c(1:2)]),center=FALSE))),quali.sup=1)
plot(resPCA,habillage=1,label="none")

summary(apply(panelN[,-c(1,2)],1,sum))





pheatmap(as.matrix(scale(t(panelN[,-c(1,2,3)]),center=FALSE)))

# Panel 2: problems ! ! qu'est ce qui devrait se sommer à 100% ?  
for(i in 1:length(newcols))
{
  
}
       