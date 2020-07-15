#####################################
##### Summarize GISAID metadata ######
#####################################

### Author: Elizabeth Tseng (etseng@pacb.com)
### Date: March 2020

#********************** Taking arguments

args <- commandArgs(trailingOnly = TRUE)
csv.file <- args[1];

submitter.min_count <- 10;
if (length(args)>=2) {
  if (args[2]=='--mincount') {
    submitter.min_count <- args[3];
  }
}

report.file <- "gisaid_metadata_report.by_submitter.pdf"

#********************** Packages (install if not found)

list_of_packages <- c("ggplot2", "scales", "reshape", "gridExtra", "grid", "dplyr")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

library(ggplot2)
library(scales)
library(reshape)
library(gridExtra)
library(grid)
library(dplyr)

# -----------------------------------------
# humans only
# -----------------------------------------
x.all <- read.table(csv.file,sep=',',header=T)
x <- x.all[x.all$Host=='Human',]
x$Msg1 <- as.character(x$Msg1)
x[is.na(x$Msg1),"Msg1"] <- "PASS"

# -----------------------------------------
# group by submitter, sort by count
# -----------------------------------------
t <- x %>% group_by(Submitter) %>% summarise(count=n()) %>% arrange(desc(count));
t.top <- filter(t, count >= submitter.min_count);

###** Output plots

pdf(file=report.file, width = 6.5, height = 6.5)

# cover
grid.newpage()
cover <- textGrob("GISAID submitter report",
    gp=gpar(fontface="italic", fontsize=40, col="orangered"))
grid.draw(cover)

for (i in 1:dim(t.top)[1]) {

  submitter.name <- t.top$Submitter[i];
  p <- x[x$Submitter==submitter.name,];

  tmp.lab <- strwrap(as.character(p$Submitting.lab), width = 20, simplify = FALSE) # modify 30 to your needs
  p$Submitting.lab <- sapply(tmp.lab, paste, collapse = "\n")

  tmp.authors <- as.character(p$Authors);
  n <- length(tmp.authors);
  if (n > 100) {
    tmp.authors <- paste(substr(tmp.authors, 1, 50), substr(tmp.authors, n-50, n), sep=' ... ');
  }
  tmp.authors <- strwrap(tmp.authors, width = 40, simplify = FALSE) # modify 30 to your needs
  p$Authors <- sapply(tmp.authors, paste, collapse = "\n")
  tmp.addr <- strwrap(as.character(p$Address.2), width = 40, simplify = FALSE) # modify 30 to your needs
  p$Address.2 <- sapply(tmp.addr, paste, collapse = "\n")

  info.table <- p %>% group_by(Submitting.lab, Authors, Address.2) %>% summarise(count=n()) %>% arrange(desc(count));
  if (dim(info.table)[1]>=3) { # just show the top 3
    info.table <- info.table[1:3,]
  }

  table1 <- tableGrob(info.table, rows=NULL, cols=c("Lab", "Authors", "Address", "Count"), theme=ttheme_default(base_size = 8));
  text1 <- textGrob(paste("Submitter:", submitter.name, sep=' '), gp=gpar(fontface="bold", fontsize=12), vjust = 0);

  p.tech <- ggplot(p, aes(Sequencing.technology, fill=Msg1)) + geom_bar() + xlab("") + ylab("Count")

  grid.arrange(text1, table1, p.tech, ncol=1, heights=c(1,15,5))

}


dev.off()


# make CSV
t <- x %>% group_by(Submitter, Submitting.lab, Authors, Address.2) %>% summarise(count=n()) %>% arrange(desc(count))
write.table(t, "gisaid.submitting_lab_authors_address.csv", sep=',',quote=T, row.names=F)