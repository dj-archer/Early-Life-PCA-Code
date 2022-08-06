#Generate outputs of PCA analysis (plots and PERMANOVA)

install.packages("factoextra")
install.packages('vegan')

library(factoextra)
library(vegan)

setwd("~/Desktop/PCAFiles/Data Tables")
#Replace with your working directory

PCA_output <- function(filename){
  
  setwd("~/Desktop/PCAFiles/Data Tables")
  #Replace with your working directory
  
  datafile <- paste(filename)
  
  data <- read.csv(datafile, na = 'NA', strip.white = TRUE)
  #Replace with your excel csv file. I have the variables defined in the columns, and the mouse ids as the rows. 
  
  variables <- data[1:nrow(data), 6:ncol(data)]
  #Only include columns (4:18) with numerical data, or your response variables.
  #Col #1 = 1; Row #1 = 0

  groups = as.factor(data$Treatment[1:nrow(data)])
  #This defines the condition/variable you're interested in running the PCA on (eg: treatment condition).

  PCA = prcomp(variables,
               center = TRUE, 
               scale. = T)
  #Center and scale your variables, ensure these arguments are set to TRUE.
  
  treatment <- data$Treatment[1:nrow(data)]
  exp <- data$Experiment[1:nrow(data)]
  sex <- data$Sex[1:nrow(data)]
  
  #Extract R^2 and P-value from PERMANOVA output to put on PCA plot
  PERMANOVAtreat <- adonis(variables ~ treatment, data = data, method="manhattan")
  PERMANOVAtreatR2 <-round(PERMANOVAtreat$aov.tab$R2[1], digits = 4)
  PERMANOVAtreatP_val <-toString(round(PERMANOVAtreat$aov.tab$"Pr(>F)"[1], digits = 4))
  
  
  #fviz_eig(PCA)
  #Generate a scree plot to look at variance of each PC
  
  #dev.new()
  
  #To make publication quality figures,
  
  #Generating axix(/PC) labels
  eig.val <- get_eigenvalue(PCA)
  eig.val
  
  #Get PC1 (x axis) variance 
  PCnum1 = round(eig.val[1,2], digits = 2)
  #Make label for x axis fpr PC1
  PC1 = paste("PC1 ", "(", PCnum1, "%", ")", sep="")
  print(PC1)
  
  #Get PC2 (y axis) variance
  PCnum2 = round(eig.val[2,2], digits = 2)
  #Make label for y axis for PC2
  PC2 = paste("PC2 ", "(", PCnum2, "%", ")", sep="")
  print(PC2)
  
  #Get PC3 (y axis) variance
  PCnum3 = round(eig.val[3,2], digits = 2)
  #Make label for y axis for PC3
  PC3 = paste("PC3 ", "(", PCnum3, "%", ")", sep="")
  print(PC3)
  
  #The ind.p code below is the main code that will generate the PCA plot
  ind.p = fviz_pca_biplot(
    PCA,
    geom="point", #Makes the sample names invisble
    #If you have specific colours you want to attribute to specific groups use the habillage and pallete functions
    #Otherwise get rid of these functions for random colours
    habillage = data$Treatment[1:nrow(data)],
    palette = c("Early Colonization" = "#4B8AF3",#blue
                "Delayed Colonization" = "#FF2D38",#red
                "L. reuteri R2lc + Delayed Colonization" = "#00EAE1",#turquoise
                "L. reuteri PB-W1 + Delayed Colonization" = "#570094",#purple
                "Conventional" = "#FFAD08",#orange
                "Germ-free" = "#000000"),#black
    axes = c(1,2),
    col.ind = groups,
    addEllipses = TRUE,
    ellipse.level = 0.95, #95 percent confidence interval for ellipses
    ellipse.type = "confidence",
    legend.title = "",
    select.var = list(contrib = 50),
    col.var = "Steel Blue", 
    invisible = "none",
    ggtheme = theme_gray(),
    pointshape = 16, 
    labelsize = 10,
    pointsize = 10,
    mean.point = "false",
    repel = TRUE
  )
  #ind.p
  
  #Use this prettyplot code to annotate the plot with text (eg: R^2 values and/or P-values)
  prettyplot = ggpubr::ggpar(ind.p,
                title = "",
                subtitle = "",
                font.x=24, font.y = 24, 
                font.caption = 14, font.legend = 14, font.main = 14,
                font.tickslab = 20,
                xlab = PC1, ylab = PC2, #Change this based on what the axis titles are in the original plot (eg: PC1 vs PC2 or PC2 vs PC3...)
                ggtheme = theme_classic())
  
  xrange <- range(ind.p$data$x)
  yrange <- range(ind.p$data$y)
  
  prettyplot + annotate("text", 
                        x=((xrange[2])), 
                        y=((yrange[1])+1):(yrange[1]),
                        label = c((paste("italic(R) ^ 2 ==",PERMANOVAtreatR2)),(paste("italic(P) ==",PERMANOVAtreatP_val))), 
                        size = 6,
                        parse = TRUE) + coord_cartesian(clip = 'off')
  
  
  #dev.set(dev.prev())
  #dev.set(dev.next())
  date <- format(Sys.Date(), format="%b,%e,%Y")
  plotfile <- paste(filename,"(",date,")",".pdf",sep="")
  
  #Use this to save the PCA plots to the PCA Plots folder (as defined in the directory below)
  #The file name is defined as "plotfile" which will use the file name of the data table as the plot file name
  ggsave(file=file.path("~/Desktop/PCAFiles/PCA Plots", plotfile), width=400, height=300, units=c("mm"), dpi=300)
  
  graphics.off()
  #Use this to clear all graphical windows
  
  #Perform PERMANOVA to determine whether the centroid and distributions of the ellipses are significantly different.
  
  #Assign a value to each variable in your data set (eg: treatment, sex of mice, etc.)
  treatment <- data$Treatment[1:nrow(data)]
  exp <- data$Experiment[1:nrow(data)]
  sex <- data$Sex[1:nrow(data)]
  cage <- data$Cage[1:nrow(data)]
  cage <- data$Cage[1:nrow(data)]
  
  #Run the PERMANOVA for each variable
  #eg: "PERMANOVAtreatout" will output the R^2 and P-value of the PERMANOVA with respect to the treatment condition
  PERMANOVAtreat <- adonis(variables ~ treatment, data = data, method="manhattan")
  PERMANOVAtreatout <-list(filename,PERMANOVAtreat$aov.tab)
  
  PERMANOVAexp <- adonis(variables ~ exp, data = data, method="manhattan")
  PERMANOVAexpout <-PERMANOVAexp$aov.tab
  
  PERMANOVAsex <- adonis(variables ~ sex, data = data, method="manhattan")
  PERMANOVAsexout <-PERMANOVAsex$aov.tab
  
  PERMANOVAcage <- adonis(variables ~ cage, data = data, method="manhattan")
  PERMANOVAcageout <-list(filename,PERMANOVAcage$aov.tab)
  
  # x*y = x + y + x:y (x:y = x interacts with y (: = interacts with))
  #Use these lines if you are interested in the interaction between two variables (eg: sex and treatment)
  PERMANOVAtreatcage <- adonis(variables ~ treatment*cage, data = data, method="manhattan")
  PERMANOVAtreatcageout <-list(filename,PERMANOVAtreatcage$aov.tab)
  
  PERMANOVAtreatsex <- adonis(variables ~ treatment*sex, data = data, method="manhattan")
  PERMANOVAtreatsexout <-list(filename,PERMANOVAtreatsex$aov.tab)
  
  PERMANOVAtreatexp <- adonis(variables ~ treatment*exp, data = data, method="manhattan")
  PERMANOVAtreatexpout <-list(filename,PERMANOVAtreatexp$aov.tab)
  
  #Use this next line (get rid of the #) if you want multiple PERMANOVA tests exported from the function 
  #instead of one test results exported and change return(PERMANOVAtreatexpout) to return(PERMANOVAall)
  #PERMANOVAall <- list(filename," ","Treatment",PERMANOVAtreatout," ","Experiment",PERMANOVAexpout," ","Sex",PERMANOVAsexout," ","Cage",PERMANOVAcageout)
  
  return(PERMANOVAtreatexpout)

}

#The code below will run the above function for all files in a folder.
date <- format(Sys.Date(), format="%b,%e,%Y")
#Change the directory in the next line to match where your file(s) are located.
listoffiles <- list.files(path = "~/Desktop/PCAFiles/Data Tables")
filenum <- 1
ALLPERM <- c()

while (filenum <= length(listoffiles)) {
    file <- listoffiles[filenum]
    PCA_output(filename = file)
    ALLPERM <- c(ALLPERM, PCA_output(filename = file))
    filenum <- filenum + 1 
   if (filenum == length(listoffiles)) { 
     file <- listoffiles[filenum]
     PCA_output(filename = file)
     filenum <- filenum + 1 
     ALLPERM <- c(ALLPERM, PCA_output(filename = file))
     print("done")
     break
  } 
} 

#This will generate a csv file with all of the PERMANOVA stats that were run in the above function 
#and save it to the PCA Plots folder (as defined in the directory below)
write.csv(ALLPERM, file=file.path("~/Desktop/PCAFiles/PERMANOVA Files", paste("ALLPermanovaTests","(",date,")",".csv")))
