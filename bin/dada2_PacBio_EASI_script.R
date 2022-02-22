#.libPaths("/work_ifs/sukmb276/Microbiome/clean_data_from_dada2/R_libraries_new/1.14")
library("dada2"); packageVersion("dada2")
## [1] '1.12.1'
library("Biostrings"); packageVersion("Biostrings")
## [1] '2.46.0'
library("ShortRead"); packageVersion("ShortRead")
## [1] '1.36.1'
library("ggplot2"); packageVersion("ggplot2")
## [1] '3.1.0'
library("reshape2"); packageVersion("reshape2")
## [1] '1.4.3'
#library(gridExtra); packageVersion("gridExtra")
## [1] '2.3'
#library(phyloseq); packageVersion("phyloseq")
## [1] '1.25.2'
#library(dplyr); packageVersion("dplyr")
## [1] '1.4.3'

###################################
####                          #####
####        FUNCTIONS         #####
####                          #####
###################################

get_stats <- function (fl, n = 5e+05, aggregate = FALSE) 
{
  require(ShortRead)
  statdf <- data.frame(Cycle = integer(0), Mean = numeric(0), 
                       Q25 = numeric(0), Q50 = numeric(0), Q75 = numeric(0), 
                       Cum = numeric(0), file = character(0))
  anndf <- data.frame(minScore = numeric(0), label = character(0), 
                      rclabel = character(0), rc = numeric(0), file = character(0))
  FIRST <- TRUE
  for (f in fl[!is.na(fl)]) {
    srqa <- qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    rc <- sum(srqa[["readCounts"]]$read)
    if (rc >= n) {
      rclabel <- paste("Reads >= ", n)
    }
    else {
      rclabel <- paste("Reads: ", rc)
    }
    means <- rowsum(df$Score * df$Count, df$Cycle)/rowsum(df$Count, 
                                                          df$Cycle)
    get_quant <- function(xx, yy, q) {
      xx[which(cumsum(yy)/sum(yy) >= q)][[1]]
    }
    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.25), simplify = TRUE)
    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.5), simplify = TRUE)
    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, 
                                                     foo$Count, 0.75), simplify = TRUE)
    cums <- by(df, df$Cycle, function(foo) sum(foo$Count), 
               simplify = TRUE)
    if (!all(sapply(list(names(q25s), names(q50s), names(q75s), 
                         names(cums)), identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    if (FIRST) {
      plotdf <- cbind(df, file = basename(f))
      FIRST <- FALSE
    }
    else {
      plotdf <- rbind(plotdf, cbind(df, file = basename(f)))
    }
    statdf <- rbind(statdf, data.frame(Cycle = as.integer(rownames(means)), 
                                       Mean = means, Q25 = as.vector(q25s), Q50 = as.vector(q50s), 
                                       Q75 = as.vector(q75s), Cum = 10 * as.vector(cums)/rc, 
                                       file = basename(f)))
    anndf <- rbind(anndf, data.frame(minScore = min(df$Score), 
                                     label = basename(f), rclabel = rclabel, rc = rc, 
                                     file = basename(f)))
  }
  anndf$minScore <- min(anndf$minScore)
  statdf
}

make_qual_plot <- function(x){
  p <- ggplot(data = x, aes(x = Cycle, y = Score)) + geom_line(data = x, aes(y = Mean), 
                                                               color = "#66C2A5") + geom_line(data = x, aes(y = Q25), 
                                                                                              color = "#FC8D62", size = 0.25, linetype = "dashed") + 
    geom_line(data = x, aes(y = Q50), color = "#FC8D62", 
              size = 0.25) + geom_line(data = statdf, aes(y = Q75), 
                                       color = "#FC8D62", size = 0.25, linetype = "dashed") + 
    ylab("Quality Score") + xlab("Cycle") + theme_bw() + 
    theme(panel.grid = element_blank()) + guides(fill = FALSE) + 
    ylim(c(0, NA))
  if (length(unique(x$Cum)) > 1) {
    p <- p + geom_line(data = x, aes(y = Cum), color = "red", 
                       size = 0.25, linetype = "solid") + scale_y_continuous(sec.axis = sec_axis(~. * 
                                                                                                   10, breaks = c(0, 100), labels = c("0%", "100%"))) + 
      theme(axis.text.y.right = element_text(color = "red"), 
            axis.title.y.right = element_text(color = "red"))
  }
  p
}





#################
###           ###
###   CODE    ###
###           ###
#################

# ARGUMENTS: output directory, threads, keep intermediate files,
args = commandArgs(trailingOnly=TRUE)
#args=c("m54349Ue",24)

libraryid = args[1]
threads=as.numeric(args[2])
numbases=1e+09

biosamp=read.csv("biosample.csv", stringsAsFactors=F, head=T)
biosamp$SeqID = paste0(libraryid,".",biosamp$Barcode)



path2 <- "raw" # CHANGE ME to location of the First Replicate fastq files
path.out <- "Figures/"
path.rds <- "RDS/"
dir.create(path.out)
dir.create(path.rds)


#statdf = get_stats(filts2[1])
#make_qual_plot(statdf)


fns2 <- list.files(path2, pattern="fastq.gz", full.names=TRUE)
#fns2 = sample(fns2, 10)
#fns2<-fns2[grep("ccs",fns2)]
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
rc <- dada2:::rc
theme_set(theme_bw())

sample.names = gsub("[.]fastq[.]gz$", "", basename(fns2))

fileframe=data.frame(sample.names=sample.names, fns2=fns2, nops2 =file.path("noprimers", paste0(sample.names, ".noprimer.fastq.gz")), stringsAsFactors=F)

library(parallel)
mclapply(seq_along(fileframe$sample.names), function(x) removePrimers(fileframe$fns2[x], fileframe$nops2[x], primer.fwd=F27, primer.rev=dada2:::rc(R1492), orient=TRUE), mc.cores=threads)

#lens.fn <- mclapply(fileframe$nops2, function(fn) nchar(getSequences(fn)), mc.cores=threads)
#lens <- do.call(c, lens.fn)
#summary(lens)

fileframe$filts2 <- file.path("filtered", paste0(sample.names, ".noprimer.filtered.fastq.gz"))
track2 <- filterAndTrim(fileframe$nops2[file.exists(fileframe$nops2)], fileframe$filts2[file.exists(fileframe$nops2)], minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2, multithread = threads, qualityType="FastqQuality")

drp2 <- derepFastq(fileframe$filts2[file.exists(fileframe$filts2)], verbose=TRUE, qualityType="FastqQuality")
names(drp2) <- fileframe$sample.names[file.exists(fileframe$filts2)]

err2 <- learnErrors(drp2, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=threads, qualityType="FastqQuality", nbases=numbases, randomize=T)

saveRDS(err2, file.path(path.rds, "err.rds"))

#plotErrors(err2)
#dd2 <- dada(drp2, err=err2, BAND_SIZE=32, multithread=threads)
dd2 <- dada(drp2, err=err2, BAND_SIZE=32, multithread=threads, pool="pseudo")

saveRDS(dd2, file.path(path.rds, "dd2.rds"))


st2 <- makeSequenceTable(dd2); dim(st2)

st.final <- removeBimeraDenovo(st2, method="consensus", minFoldParentOverAbundance=3.5, multithread=threads)

#seqcounts=mclapply(unlist(fileframe[,-1]), function(x) if(file.exists(x)){ShortRead::countLines(x)/4})

trackreads = data.frame(ccs=unlist(mclapply(fileframe$fns2, function(x) return(ifelse(file.exists(x),ShortRead::countLines(x)/4,0)))),
                primers=unlist(mclapply(fileframe$nops2, function(x) return(ifelse(file.exists(x),ShortRead::countLines(x)/4,0)))), 
                filtered=unlist(mclapply(fileframe$filts2, function(x) return(ifelse(file.exists(x),ShortRead::countLines(x)/4,0)))),
                denoised=sapply(dd2, function(x) sum(x$denoised))[fileframe$sample.names], tabled=rowSums(st2)[fileframe$sample.names], nonchim=rowSums(st.final)[fileframe$sample.names])

trackreads[is.na(trackreads)]<-0

rownames(trackreads) <- fileframe$sample.names

trackreads.rel <- t(apply(trackreads, 1, function(x) x/x[1]))

#trackreads %>% reshape2::melt(.) %>% ggplot(aes(x=Var2, y=value, group=Var1)) + geom_point() + geom_line() 
#trackreads.rel %>% reshape2::melt(.) %>% ggplot(aes(x=Var2, y=value, group=Var1)) + geom_point() + geom_line() + geom_line(data= . %>% group_by(Var2) %>% summarise(value=mean(value)) %>% mutate(Var1 = "Mean"), col="green")

#### Tax assignments
asv.names = data.frame(name = paste0("ASV_",sprintf("%06d", seq_along(colnames(st.final)))), seq = colnames(st.final), stringsAsFactors = F)

st.final.renamed <- st.final
colnames(st.final.renamed) <- asv.names$name

tax_rdp <- assignTaxonomy(st.final, "/rdp_train_set_18.fa.gz", multithread=threads) # Slowest part
tax_rdp.sp <- addSpecies(tax_rdp, "/rdp_species_assignment_18.fa.gz")

#tax_silva <- assignTaxonomy(st.final, "/work_ifs/sukmb276/Microbiome/clean_data_from_dada2/reference_data//silva_nr99_v138.1_train_set.fa.gz", multithread=threads) # Slowest part
#tax_silva.sp <- addSpecies(tax_silva, "/work_ifs/sukmb276/Microbiome/clean_data_from_dada2/reference_data//silva_species_assignment_v138.1.fa.gz")

tax_gtdb <- assignTaxonomy(st.final, "/GTDB_bac120_arc122_ssu_r202_Genus.fa.gz", multithread=threads) # Slowest part
tax_gtdb.sp <- addSpecies(tax_gtdb, "/GTDB_bac120_arc122_ssu_r202_Species.fa.gz")

#rownames(tax_rdp.sp) <- rownames(tax_silva.sp) <- rownames(tax_gtdb.sp) <- asv.names$name
rownames(tax_rdp.sp) <- rownames(tax_gtdb.sp) <- asv.names$name

tax_rdp.sp.final <- t(apply(tax_rdp.sp, 1, function(x){x=gsub(" ","_", x); x[7]=ifelse(!is.na(x[7]), paste0(x[6],"_",x[7]), NA); x=ifelse(!is.na(x), paste0(substr(colnames(tax_rdp.sp), 1, 1), "_", x), NA); return(x)}))
#tax_silva.sp.final <- t(apply(tax_silva.sp, 1, function(x){x=gsub(" ","_", x); x[7]=ifelse(!is.na(x[7]), paste0(x[6],"_",x[7]), NA); x=ifelse(!is.na(x), paste0(substr(colnames(tax_silva.sp), 1, 1), "_", x), NA); return(x)}))
tax_gtdb.sp.final <- t(apply(tax_gtdb.sp, 1, function(x){x=gsub(" ","_", x); x[7]=ifelse(!is.na(x[7]), paste0(x[6],"_",x[7]), NA); x=ifelse(!is.na(x), paste0(substr(colnames(tax_gtdb.sp), 1, 1), "_", x), NA); return(x)}))

#####

prefilt = sapply(fileframe$nops2[file.exists(fileframe$nops2)], get_stats, simplify = F)
postfilt = sapply(fileframe$filts2[file.exists(fileframe$filts2)], get_stats, simplify = F)
names(prefilt) <- fileframe$sample.names[file.exists(fileframe$nops2)]
names(postfilt) <- fileframe$sample.names[file.exists(fileframe$filts2)]

out = list(
  library_id = libraryid,
  sample_ids = biosamp,
  asv.names = asv.names, 
  asv = st.final.renamed, 
#  tax = list("SILVAv138" = tax_silva.sp.final, "GTDBr202" = tax_gtdb.sp.final, "RDP18" = tax_rdp.sp.final),
  tax = list("GTDBr202" = tax_gtdb.sp.final, "RDP18" = tax_rdp.sp.final),
  stats = list(absolute = trackreads, relative = trackreads.rel),
  error = err2,
  qual = list(prefilt=prefilt, postfilt=postfilt)
  )

saveRDS(out, "allout.rds")

dir.create("tables",recursive=F)

st.final.renamed2=st.final.renamed
rownames(st.final.renamed2) = biosamp[match(rownames(st.final.renamed), biosamp$SeqID),"Bio.Sample"]
        write.table(st.final.renamed2, paste0("tables/",libraryid,"_","ASV","_",names(out$tax)[1],".tsv"), sep="\t")

taxtab = data.frame(t(apply(out$tax[[1]], 1, function(x){x[is.na(x)] <- paste0(rev(x[!is.na(x)])[1], "_unclassified"); return(x)})), stringsAsFactors=F)
for(tl in colnames(taxtab)){
        taxlab = taxtab[, tl]
        abu_lev = t(data.frame(aggregate(. ~ taxlab, data.frame(t(st.final.renamed)), sum), row.names = 1))
        rownames(abu_lev) = biosamp[match(rownames(abu_lev), gsub("-",".",biosamp$SeqID)),"Bio.Sample"]
        write.table(abu_lev, paste0("tables/",libraryid,"_",tl,"_",names(out$tax)[1],".tsv"), sep="\t")
}

write.table(cbind(asv.names, taxtab[match(asv.names$name, rownames(taxtab)),]), paste0("tables/",libraryid,"_","ASV","_",names(out$tax)[1],"_taxonomy.tsv"), sep="\t", quote=F)




