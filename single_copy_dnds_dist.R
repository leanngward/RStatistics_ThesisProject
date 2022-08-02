
## Set Working Directory ####
  setwd("C:/Users/leann/OneDrive/grad research_newer/butterfly_aging/run5_asproteins_orthofinder/longevitygenesonly_results")

## Load Libraries ####
  library("readxl")
  library(ggplot2)
  library(reshape2)
  library(plyr)
  library(qvalue)
  #library(svglite)
  library(tidyr)

## Set in file names ####

  infile <- "branch_output_heliconius_clade.tsv"
  mapfile <- "unique_go_term_mapping_sc_hmel.txt"
  list <- "genage_go_terms.txt"

## Find the subset of orthogroups in the longevity list ####

  listdf <- read.delim(list,header=TRUE,sep="\t")
  longevitylist <- listdf$V5
  head(longevitylist)
  
  mapdf <- read.delim(mapfile,header=FALSE,sep="\t")

  inlist <- rep("NA",1)
  for (row in 1:nrow(mapdf)){
    go <- unlist(strsplit(mapdf[row,2],","))
    for (i in go) {
      if (i %in% longevitylist) {
        inlist <- append(inlist,mapdf[row,1])
      }
    }
  
  }
  
  head(inlist)
  inlist <- na.omit(inlist)
  inlist <- unique(inlist)
  print(inlist)
  
## Read in branch model data ####
  ld1 = read.delim(infile, header = TRUE, sep = "\t")
  head(ld1)
  
  ld <- ld[ld$Group %in% inlist, ]
  nrow(ld)
  
#Calculate log ratio 
  ld$LR <- -2 * (as.numeric(ld$Null.lnL) - as.numeric(ld$Alt.lnL))
    head(ld)
#Change column names
    colnames(ld) <- c("Group", "nulllnl", "nullomega", "altlnl", "bg", "fg", "LR")
#Calculate significance of the log ratio
  ld$pval <- pchisq(ld$LR, df = 1, lower.tail = FALSE)
  ld$padj <- p.adjust(ld$pval, method = 'fdr')
  
    print("# of genes with omega higher in Heliconius")
    head(ld)
  nrow(ld)

    #number of omegas where the background is higher than the foreground
  nrow(subset(ld,bg  > fg))


#Sort Results by LRT
  sorted_ld = ld[order(-ld$LR),]
    head(sorted_ld)

#Convert dataframe into long format
  wide_ld = subset(sorted_ld, select = c(Group,bg,fg,LR,pval,padj))
    head(wide_ld)
  long_ld = gather(wide_ld, level, omega, bg, fg)
    head(long_ld)

#remove "omega" form bg and fg
  long_ld$level <- replace(long_ld$level, long_ld$level=='bg', 'Background')
  long_ld$level <- replace(long_ld$level, long_ld$level=='fg', 'Heliconius')
head(long_ld)  

#bin the omegas
  c = 3
  long_ld$omega_bin <- long_ld$omega
  long_ld$omega_bin[long_ld$omega_bin > c] <- c
  head(long_ld)
  nrow(long_ld)
  nrow(long_ld[long_ld$omega > .25, ])

#Graph the results.
  mu <- ddply(subset(long_ld), "level", summarise, grp.mean=mean(omega_bin))
  mu

avg_bg = mu[1,2]
favg_bg = format(round(avg_bg,2), nsmall=2)

avg_fg = mu[2,2]
favg_fg = format(round(avg_fg,2), nsmall=2)

myplot <- ggplot(subset(long_ld), aes(x=omega_bin, fill = level, color=level)) +
  ggtitle("Density Distribution of dN/dS Estimates for Longevity-Associated genes in Heliconius clade and Background Species") +
  geom_density(alpha = 0.6) +
  xlab('dN/dS') + 
  ylab('Density') +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=level),
             linetype="dashed") +
  #geom_text(aes(x = .1, label = paste0("omega = ",favg_bg), y=10), colour = "red") +
  #geom_text(aes(x = .18, label = paste0("omega = ",favg_fg), y=8), colour = "blue") +
  scale_color_manual(values=c(Background="red", Heliconius="darkblue")) +
  scale_fill_manual(values=c(Background="red", Heliconius="darkblue")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line =   element_line(colour = "black"))


myplot

#Perform stats - Kolgomorov-Smirnoff test (D=0.10349, P = 4.515e-06).
  print(favg_fg)
  print(avg_bg)
  bg_omega = subset(long_ld, level == "Background")$omega_bin
  fg_omega = subset(long_ld, level == "Heliconius")$omega_bin
  ks.test(fg_omega, bg_omega)
