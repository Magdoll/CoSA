#####################################
##### Summarize GISAID metadata ######
#####################################

### Author: Elizabeth Tseng (etseng@pacb.com)
### Date: March 2020

#********************** Taking arguments

args <- commandArgs(trailingOnly = TRUE)
csv.file <- args[1]

report.file <- "gisaid_metadata_report.pdf"
report.figdir <- "gisaid_figures"

if (!file.exists(report.figdir)){
    dir.create(report.figdir)
}

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
# plot by specimen source
# -----------------------------------------
x.all <- read.table(csv.file,sep=',',header=T)
x <- subset(x.all, Host=='Human')
p.source <- ggplot(x, aes(x=Specimen.source, fill=Specimen.source)) + geom_bar() + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none") + labs(title="Specimen Source") + xlab("") + ylab("Count")
ggsave(file.path(report.figdir,'Rplot.specimen_source.png'), dpi=200, width=6, height=4)
# -----------------------------------------
# plot geography
# -----------------------------------------
p.continent <- ggplot(x, aes(Continent, fill=Continent)) + geom_bar() + coord_flip() +
  xlab("") + ylab("") + labs(title="Sequences by Continent") + theme(legend.position="none")
ggsave(file.path(report.figdir,'Rplot.geography_byContinent.png'), dpi=200, width=6, height=4)

p.asia <- ggplot(x[x$Continent=='Asia',], aes(Country, fill=Country)) + geom_bar() + coord_flip() +
  xlab("") + ylab("") + labs(title="Sequences by Country, Asia") + theme(legend.position="none")
ggsave(file.path(report.figdir,'Rplot.geography_asia.png'), dpi=200, width=6, height=4)

p.europe <- ggplot(x[x$Continent=='Europe',], aes(Country, fill=Country)) + geom_bar() + coord_flip() +
  xlab("") + ylab("") + labs(title="Sequences by Country, Europe") + theme(legend.position="none")
ggsave(file.path(report.figdir,'Rplot.geography_europe.png'), dpi=200, width=6, height=4)

p.na <- ggplot(x[x$Continent=='North America',], aes(Location, fill=Location)) + geom_bar() + coord_flip() +
  xlab("") + ylab("") + labs(title="Sequences in NA") + theme(legend.position="none")
ggsave(file.path(report.figdir,'Rplot.geography_na.png'), dpi=200, width=6, height=4)

# -----------------------------------------
# plot sequence qualities (only if already run through FAIL/PASS
# -----------------------------------------
FAIL_PASS_FIELDS <- c("Pass1", "Msg1", "GapCount", "NCount")

if (all(FAIL_PASS_FIELDS %in% colnames(x))) {

  # fail/pass and reasons for fail
  x$Msg1 <- as.character(x$Msg1)
  x[is.na(x$Msg1),"Msg1"] <- "PASS"
  x$Status=x$Msg1
  x[x$Status=='Too many ambiguous bases', "Status"] <- "10+ ambiguous bases"
  x[x$Status=='Too many gaps', "Status"] <- "2+ stretches of 'N's"
  x[x$Status=='Too short', "Status"] <- "Less than 28kb"
  p.passfail <- ggplot(x, aes(x=Status, fill=Status)) + geom_bar() +
    xlab("") + ylab("Count") + labs(title="Sequence quality")
  ggsave(file.path(report.figdir,'Rplot.pass_fail.png'), dpi=200, width=6.5, height=4)

  # number of 'N's
  x$NCountCat <- "0"
  x[x$NCount>=1,"NCountCat"] <- "1-10"
  x[x$NCount>=11,"NCountCat"] <- "11-100"
  x[x$NCount>=101,"NCountCat"] <- "101-1000"
  x[x$NCount>=1001,"NCountCat"] <- ">1000"
  x$NCountCat <- factor(x$NCountCat)
  x$NCountCat <- ordered(x$NCountCat, levels=c("0", "1-10", "11-100", "101-1000", ">1000"))
  p.Ns <- ggplot(x, aes(NCountCat, fill=NCountCat)) + geom_bar() + xlab("") + ylab("Count") +
    labs(title="Number of 'N's in sequences")
  ggsave(file.path(report.figdir,'Rplot.count_Ns.png'), dpi=200, width=6.5, height=4)

  # number of segments
  x$Segments <- "1"
  x[x$GapCount>=1, "Segments"] <- "2"
  x[x$GapCount>=2, "Segments"] <- "3-5"
  x[x$GapCount>=5, "Segments"] <- "6+"
  x$Segments <- as.factor(x$Segments)
  x$Segments <- ordered(x$Segments, levels=c("1", "2", "3-5", "6+"))
  p.segments <- ggplot(x, aes(Sequencing.technology, fill=Segments)) + geom_bar() + xlab("") + ylab("Count") +
    labs(title="Number of continuous (non-N) sequence segments")
  ggsave(file.path(report.figdir,'Rplot.count_segments.png'), dpi=200, width=6.5, height=4)

}

###** Output plots

pdf(file=report.file, width = 6.5, height = 6.5)

# cover
grid.newpage()
cover <- textGrob("GISAID metadata report",
    gp=gpar(fontface="italic", fontsize=40, col="orangered"))
grid.draw(cover)

grid.arrange(p.species, p.continent, ncol=1)
print(p.asia)
print(p.europe)
print(p.na)

grid.arrange(p.passfail, p.Ns, p.segments, ncol=1)
dev.off()
