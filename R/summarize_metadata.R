#####################################
##### Summarize GISAID metadata ######
#####################################

### Author: Elizabeth Tseng (etseng@pacb.com)
### Date: March 2020

#********************** Taking arguments

args <- commandArgs(trailingOnly = TRUE)
csv.file <- args[1]
flag.by_tech <- FALSE;
if (length(args)>=2) {
  if (args[2]=='--bytech') {
    flag.by_tech <- TRUE;
  }
}

report.file <- "gisaid_metadata_report.pdf"
#report.figdir <- "gisaid_figures"
#
#if (!file.exists(report.figdir)){
#    dir.create(report.figdir)
#}

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
x <- x.all[x.all$Host=='Human',]
p.source <- ggplot(x, aes(x=Specimen.source, fill=Specimen.source)) + geom_bar() + theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none") + labs(title="Specimen Source") + xlab("") + ylab("Count")
#ggsave(file.path(report.figdir,'Rplot.specimen_source.png'), dpi=200, width=6, height=4)
# -----------------------------------------
# plot geography
# -----------------------------------------
p.continent <- ggplot(x, aes(Continent, fill=Continent)) + geom_bar() + coord_flip() +
  xlab("") + ylab("") + labs(title="Sequences by Continent") + theme(legend.position="none")
#ggsave(file.path(report.figdir,'Rplot.geography_byContinent.png'), dpi=200, width=6, height=4)

p.asia <- ggplot(x[x$Continent=='Asia',], aes(Country, fill=Country)) + geom_bar() + coord_flip() +
  xlab("") + ylab("") + labs(title="Sequences by Country, Asia") + theme(legend.position="none")
#ggsave(file.path(report.figdir,'Rplot.geography_asia.png'), dpi=200, width=6, height=4)

p.europe <- ggplot(x[x$Continent=='Europe',], aes(Country, fill=Country)) + geom_bar() + coord_flip() +
  xlab("") + ylab("") + labs(title="Sequences by Country, Europe") + theme(legend.position="none")
#ggsave(file.path(report.figdir,'Rplot.geography_europe.png'), dpi=200, width=6, height=4)

p.na <- ggplot(x[x$Continent=='North America',], aes(Location, fill=Location)) + geom_bar() + coord_flip() +
  xlab("") + ylab("") + labs(title="Sequences in NA") + theme(legend.position="none")
#ggsave(file.path(report.figdir,'Rplot.geography_na.png'), dpi=200, width=6, height=4)

# -----------------------------------------
# plot sequence qualities (only if already run through FAIL/PASS
# -----------------------------------------
FAIL_PASS_FIELDS <- c("Pass1", "Msg1", "GapCount", "NCount")

if (all(FAIL_PASS_FIELDS %in% colnames(x))) {

  # fail/pass and reasons for fail
  x$Msg1 <- as.character(x$Msg1)
  x[is.na(x$Msg1),"Msg1"] <- "PASS"
  x$Status=x$Msg1
  p.passfail <- ggplot(x, aes(x=Status, fill=Msg1)) + geom_bar() +
    xlab("") + ylab("Count") + labs(title="Sequence quality")
  #ggsave(file.path(report.figdir,'Rplot.pass_fail.png'), dpi=100, width=4.5, height=2.5)

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
  #ggsave(file.path(report.figdir,'Rplot.count_Ns.png'), dpi=200, width=6, height=4)

  # number of segments
  x$Segments <- "1"
  x[x$GapCount>=1, "Segments"] <- "2"
  x[x$GapCount>=2, "Segments"] <- "3-5"
  x[x$GapCount>=5, "Segments"] <- "6+"
  x$Segments <- as.factor(x$Segments)
  x$Segments <- ordered(x$Segments, levels=c("1", "2", "3-5", "6+"))
  p.segments <- ggplot(x, aes(Segments, fill=Segments)) + geom_bar() + xlab("") + ylab("Count") +
    labs(title="Number of continuous (non-N) sequence segments")
  #ggsave(file.path(report.figdir,'Rplot.count_segments.png'), dpi=100, width=4.5, height=2.5)

}


if (flag.by_tech) {
  # ---------------
  # summarize sequencing technology
  # ---------------
  ggplot(x, aes(Sequencing.technology, fill=Sequencing.technology)) + geom_bar() + xlab("") + ylab("") + labs(title="Sequencing Technology") + theme(legend.position='none')
  ggsave('Rplot.seq_tech.png', dpi=200, width=6, height=4)
  # ---------------
  # summarize high sites by tech
  # ---------------
  t <- x %>% group_by(Country, Sequencing.technology) %>% summarise(count=n())
  ggplot(subset(t,count>=10), aes(x=Country, y=count, fill=Sequencing.technology)) + geom_bar(stat='identity') + coord_flip() + xlab("") + ylab("") + labs(title="Top sequencing countries")
  ggsave('Rplot.seq_tech_by_top_country.png', width=6, height=4, dpi=200)

  t2 <- x %>% group_by(Country) %>% summarise(count=n())
  ggplot(subset(t,count>=10), aes(x=Country, y=count)) + geom_bar(stat='identity') + coord_flip() + xlab("") + ylab("") + labs(title="Top sequencing countries")
  ggsave('Rplot.top_country.png', width=6, height=4, dpi=200)

  if (all(FAIL_PASS_FIELDS %in% colnames(x))) {
    p.passfail.by_tech <- ggplot(x, aes(x=Sequencing.technology, fill=Msg1)) + geom_bar() +
      xlab("") + ylab("Count") + labs(title="Sequence quality") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="top")
    ggsave('Rplot.pass_fail.by_tech.png', dpi=200, width=6, height=4)

    p.Ns.by_tech <- ggplot(x, aes(Sequencing.technology, fill=NCountCat)) + geom_bar() + xlab("") + ylab("Count") +
      labs(title="Number of 'N's in sequences")
    ggsave('Rplot.count_Ns.by_tech.png', dpi=200, width=6, height=4)

    p.segments.by_tech <- ggplot(x, aes(Sequencing.technology, fill=Segments)) + geom_bar() + xlab("") + ylab("Count") +
      labs(title="Number of continuous (non-N) sequence segments")
    ggsave('Rplot.count_segments.by_tech.png', dpi=200, width=6, height=4)
  }
}

###** Output plots
pdf(file=report.file, width = 6.5, height = 6.5)

# cover
grid.newpage()
cover <- textGrob("GISAID metadata report",
    gp=gpar(fontface="italic", fontsize=40, col="orangered"))
grid.draw(cover)

grid.arrange(p.source, p.continent, ncol=1)
print(p.asia)
print(p.europe)
print(p.na)

grid.arrange(p.passfail, p.Ns, p.segments, ncol=1)
dev.off()
