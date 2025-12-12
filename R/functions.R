# palettes ####
pal <- list(
  pal.earthtones = c("#4E6172","#D57500","#8F3B1B","#404F24","#613318","#668D3C"),
  pal.okabe = c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
)

# %ni% ####
# "not in" logical expression
'%ni%' <- Negate('%in%')


# ra() ####
# relative abundance transformation
ra <- function(x){x/sum(x)}

# EE() ####
# calculate expected errors in a sequence from a vector of quality scores
EE <- function(qual.scores){
  sum(10^(-qual.scores/10))
}

# plot_bar2() ####
# phyloseq bar plot without lines for each ASV
plot_bar2 <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
          title = NULL, facet_grid = NULL, width = 0.9) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack", width = width)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

# remove_primers() ####
# Function to remove primers from raw amplicon files

remove_primers <- function(metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
                           amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
                           amplicon = "ITS", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
                           sampleid.colname = "library_id", # column name in metadata containing unique sample identifier
                           fwd.fp.colname = "fwd_filepath", # name of column in metadata indicating fwd filepath to raw data
                           rev.fp.colname = "rev_filepath",
                           fwd_pattern="_R1_",
                           rev_pattern="_R2_",
                           fwd_primer="CTTGGTCATTTAGAGGAAGTAA",
                           rev_primer="GCTGCGTTCTTCATCGATGC",
                           multithread=parallel::detectCores()-1){
  
  library(tidyverse); packageVersion("tidyverse")
  library(dada2); packageVersion("dada2")
  library(purrr); packageVersion("purrr")
  library(Biostrings); packageVersion("Biostrings")
  library(ShortRead); packageVersion("ShortRead")
  library(parallel); packageVersion("parallel")
  
  # tests
  if(all(is.na(metadata[[fwd.fp.colname]]))){
    stop("Fwd filepath column has no filenames in it.")
  }
  
  # File parsing
  
  # subset metadata to just that amplicon and get fwd and rev file paths
  x <- metadata[metadata[[amplicon.colname]] == amplicon & !is.na(metadata[[fwd.fp.colname]]),]
  fnFs <- x[[fwd.fp.colname]]
  fnRs <- x[[rev.fp.colname]]

  FWD <- fwd_primer # Sequence of FWD primer
  REV <- rev_primer  # Sequence of REV primer
  
  allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  
  # Prefilter to remove reads with ambiguous (N) bases ####
  
  # build new directory names
  filtN_paths <- dirname(fnFs)
  filtN_paths <- filtN_paths[!is.na(filtN_paths)]
  filtN_paths <- file.path(filtN_paths,"filtN")
  
  # get sample_names
  sample_names <- x[[sampleid.colname]]
  # make output file paths
  fnFs.filtN <- file.path(filtN_paths, paste0(sample_names,"_filtN_fwd.fastq.gz")) # Put N-filterd files in filtN/ subdirectory
  fnRs.filtN <- file.path(filtN_paths, paste0(sample_names,"_filtN_rev.fastq.gz"))
  
  # create directories, if needed
  # create new directories as needed
  for(i in unique(filtN_paths)){
    if(!file_test("-d", i)){dir.create(i)}
  }
  
  dir.name <- dirname(fnFs) %>% unique
  
  # Check for missing files and remove them from fnFs, fnFs.filtN, fnRs, and fnRs.filtN
  extant_files <- file.exists(fnFs) & file.exists(fnRs)
  fnFs <- fnFs[extant_files]
  fnRs <- fnRs[extant_files]
  fnFs.filtN <- fnFs.filtN[extant_files]
  fnRs.filtN <- fnRs.filtN[extant_files]
  # c(fnFs,fnRs,fnFs.filtN,fnRs.filtN)[which(duplicated(c(fnFs,fnRs,fnFs.filtN,fnRs.filtN)))]
  # Remove reads with any Ns
  
  # only run this if the output files don't exist already from a previous pass through
  fnFs <- fnFs[!file.exists(fnFs.filtN)]; fnFs.filtN <- fnFs.filtN[!file.exists(fnFs.filtN)]
  fnRs <- fnRs[!file.exists(fnRs.filtN)]; fnRs.filtN <- fnRs.filtN[!file.exists(fnRs.filtN)]
  
  if(all(length(fnFs) > 0 & length(fnRs) > 0 & length(fnFs.filtN) > 0 & length(fnRs.filtN) > 0)){
    filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, 
                  maxN = 0, 
                  multithread = ifelse(multithread>1,TRUE,FALSE)) # on Windows, set multithread = FALSE
  }
  
  # rebuild fnFs.filtN and fnRs.filtN
  # make output file paths
  fnFs.filtN <- list.files(unique(filtN_paths),pattern="_filtN_fwd.fastq.gz",full.names = TRUE)
  fnRs.filtN <- list.files(unique(filtN_paths),pattern="_filtN_rev.fastq.gz",full.names = TRUE)
  
  # build cutadapt file structure
  path.cut <- file.path(dir.name,"cutadapt")
  # create directories as needed
  for(i in unique(path.cut)){
    if(!dir.exists(i)) dir.create(i)
  }
  
  # build file names for cutadapt output
  fnFs.cut <- file.path(path.cut, paste0(sample_names,"_cutadapt_fwd.fastq.gz"))
  fnRs.cut <- file.path(path.cut, paste0(sample_names,"_cutadapt_rev.fastq.gz"))
  file.exists(fnFs.cut)
  file.info(fnFs.cut)
  FWD.RC <- dada2:::rc(FWD)
  REV.RC <- dada2:::rc(REV)
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  R1.flags <- paste("-g", FWD, "-a", REV.RC) 
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  R2.flags <- paste("-G", REV, "-A", FWD.RC) 
  
  # Run Cutadapt
  
  # clean up cutadapt command log file from previous attempts, if present
  if(file.exists(paste("./R/cutadapt_commands_",amplicon,".sh",sep=""))){
    file.remove(paste("./R/cutadapt_commands_",amplicon,".sh",sep=""))
  }
  
  for(i in seq_along(fnFs.cut)) {
    if(!file.exists(fnFs.cut[i])){
      
      args = c(R1.flags, 
               R2.flags, 
               "-n", 2, 
               "--minimum-length 100",
               "--cores 0",
               "--nextseq-trim=20",
               "-o", fnFs.cut[i],
               "-p", fnRs.cut[i], 
               fnFs.filtN[i], 
               fnRs.filtN[i])
      
      y <- paste("cutadapt",paste(args,collapse = " "),collapse = " ")
      sink(paste("./R/cutadapt_commands_",amplicon,".sh",sep=""),append = TRUE)
      cat(y,"\n")
      sink(NULL)
      
      system2("cutadapt", args = args)
    } else {next}
    
  }
  
}

# run_itsxpress() ####
# Isolate ITS region 

# Valid taxa_group arguments: 
# Alveolata,Bryophyta,Bacillariophyta,Amoebozoa,Euglenozoa,Fungi,Chlorophyta,Rhodophyta,Phaeophyceae,Marchantiophyta,Metazoa,Oomycota,Haptophyceae,Raphidophyceae, Rhizaria,Synurophyceae,Tracheophyta,Eustigmatophyceae,All

run_itsxpress <- function(directory="./data/raw/cutadapt", # where cutadapted reads live
                          itsregion="ITS1", # must be "ITS1" or "ITS2"
                          taxa_group="All",
                          nthreads=(parallel::detectCores()-1),
                          fwd_pattern="cutadapt_fwd.fastq.gz",
                          rev_pattern="cutadapt_rev.fastq.gz",
                          itsxpress.path="/home/gzahn/.local/bin/itsxpress", #path to executable
                          fwd.only=TRUE){
  
  # my paths, for easy reference:
  # itsxpress.path="/home/gzahn/.local/bin/itsxpress"
  # itsxpress.path="/uufs/chpc.utah.edu/common/home/u6033249/.local/bin/itsxpress"
  
  # find the "cutadapted" files
  fwds <- list.files(directory,pattern = fwd_pattern,full.names = TRUE)
  revs <- list.files(directory,pattern = rev_pattern,full.names = TRUE)

  
run7 <-
  run7 %>%
  filter(amplicon == "ITS")
file.size(fnFs.cut)
  
    fwds <- fwds[fwds %in% fnFs.cut]
    revs <- revs[revs %in% fnRs.cut]

    
    
  # build names for outfiles
  outs_fwd <- paste0(tools::file_path_sans_ext(fwds) %>% 
                   tools::file_path_sans_ext(), 
                 "_ITSxpress.fastq.gz")
  outs_rev <- paste0(tools::file_path_sans_ext(revs) %>% 
                     tools::file_path_sans_ext(), 
                   "_ITSxpress.fastq.gz")
  
  its_dir <- directory %>% str_replace('cutadapt','ITSx')

  if(!dir.exists(its_dir)){dir.create(its_dir)}
  
  outs_fwd <- file.path(its_dir,basename(outs_fwd))
  outs_rev <- file.path(its_dir,basename(outs_rev))
  file.exists(outs_fwd)
  # build the ITSxpress command and run it on each file in turn

  if(fwd.only){
    for(i in 1:length(fwds)){
      itsxpress <- paste0(itsxpress.path,
                          " --fastq ",fwds[i],
                          " --outfile ",outs_fwd[i],
                          " --region ",itsregion,
                          " --taxa ",taxa_group,
                          " --threads ",nthreads,
                          " --log ",outs_fwd[i],".log",
                          " --single_end")
      # write commands to file
      sink("./R/itsxpress_commands_run6.sh",append = TRUE)
      cat(itsxpress,"\n")
      sink(NULL)

      cat(outs_fwd[i])
      # if(!file.exists(outs_fwd[i])){
      #   system(command = itsxpress)
      # }
    }
  }

  if(!fwd.only){
    for(i in 1:length(fwds)){
      itsxpress <- paste0(itsxpress.path,
                          " --fastq ",fwds[i],
                          " --fastq2 ",revs[i],
                          " --outfile ",outs_fwd[i],
                          " --outfile2 ",outs_rev[i],
                          " --region ",itsregion,
                          " --taxa ",taxa_group,
                          " --threads ",nthreads,
                          " --log ",outs_fwd[i],".log"
      )
      system(command = itsxpress)
    }
  }
    
    
}



# build_asv_table() ####
# function to run dada2 pipeline on trimmed amplicon reads
# error profiling and correction should be done by sequencing run
# this function will run dada2 on files from a single sequencing run, given a metadata sheet that has samples from multiple runs

build_asv_table <- function(metadata, # metadata object for multi-seq-run samples; must contain "run" column and fwd/rev filepath columns
                            run.id.colname = "run_id", # name of column in metadata indicating which sequencing run a sample comes from
                            run.id, # the run ID to perform function on, from the run.id.colname column. Enter as a character, e.g., "1"
                            amplicon.colname = "amplicon", # column name that contains the amplicon info for each sample
                            amplicon = "SSU", # which amplicon from the run are you processing (ITS, SSU, LSU, etc)?
                            sampleid.colname = "library_id", # column name in metadata containing unique sample identifier
                            fwd.fp.colname = "fwd_filepath", # name of column in metadata indicating fwd filepath to trimmed data (e.g., cutadapt)
                            rev.fp.colname = "rev_filepath", # name of column in metadata indicating rev filepath to trimmed data (e.g., cutadapt)
                            fwd.pattern = "_R1_", # pattern in filepath column indicating fwd reads (to be safe)
                            rev.pattern = "_R2_", # pattern in filepath column indicating rev reads (to be safe),
                            maxEE = c(2,2), # max expected errors for filtration step of dada2 (for single-end cases like ITS, will default to maxEE=2)
                            trim.right = NA, # amount to trim off of 3' end after quality truncation
                            truncQ = 2, # special value denoting "end of good quality sequence"
                            rm.phix = TRUE, # remove phiX sequences?
                            compress = TRUE, # gzip compression of output?
                            multithread = (parallel::detectCores() -1), # how many cores to use? Set to FALSE on windows
                            single.end = TRUE, # use only forward reads and skip rev reads and merging (e.g., for ITS data)?
                            filtered.dir = "filtered", # name of output directory for all QC filtered reads. will be created if not extant. subdirectory of trimmed filepath
                            asv.table.dir = "./ASV_Tables", # path to directory where final ASV table will be saved
                            control.col = "treatment", # col containing negative control designation
                            control.indicator = "negative_ctl", # string indicating neg control
                            random.seed = 666){
  
  # tests
  stopifnot("data.frame" %in% class(metadata))
  if(is.null(metadata[[run.id.colname]])){stop("run.id.colname is empty or not found in metadata.")}
  if(class(run.id) != 'character'){stop("run.id must be a character vector of length 1.")}
  if(class(fwd.fp.colname) != 'character'){stop("fwd.fp.colname must be a character vector of length 1, matching the column name in metadata.")}
  if(single.end & class(rev.fp.colname) != 'character'){stop("rev.fp.colname must be a character vector of length 1, matching the column name in metadata.")}
  
  # parse options
  if(single.end){
    paired <- FALSE
    maxEE <- maxEE[1]
  } else {
    paired <- TRUE
    maxEE <- maxEE
  }
  set.seed(random.seed)

  # subset to sequencing run and only rows with file names
  metadata <- 
  metadata[as.character(metadata[[run.id.colname]]) == as.character(run.id) & 
             !is.na(metadata[[fwd.fp.colname]]) & metadata[[amplicon.colname]] == amplicon,]
  
  # parse filepaths and sample names
  fns <- metadata[[fwd.fp.colname]]
  if(paired){rns <- metadata[[rev.fp.colname]]}
  
  sample_names <- metadata[[sampleid.colname]]
  names(fns) <- sample_names
  if(paired){names(rns) <- sample_names}
  
  filtered_paths <- file.path(dirname(fns),filtered.dir)
  
  filts_f <- file.path(filtered_paths,paste0(sample_names,"_filtered_fwd.fastq.gz"))
  names(filts_f) <- sample_names
  if(paired){
    filts_r <- file.path(filtered_paths,paste0(sample_names,"_filtered_rev.fastq.gz"))
    names(filts_r) <- sample_names
  }
  
  
  # create new directories as needed
  for(i in unique(filtered_paths)){
    if(!file_test("-d", i)){dir.create(i)}
  }
  
  
  fns %>% length; filts_f %>% length;rns %>% length; filts_r %>% length
  
  # # filter and trim
  # if(paired){
  #   out <- filterAndTrim(fns, filts_f, rns, filts_r, 
  #                      maxN=0, 
  #                      maxEE=maxEE, 
  #                      truncQ=truncQ,
  #                      trimRight = ifelse(any(is.na(trim.right)),0,trim.right),
  #                      rm.phix=rm.phix, 
  #                      compress=compress,
  #                      multithread=multithread)
  # } else {
  #   out <- filterAndTrim(fns, filts_f,
  #                        maxN=0,
  #                        maxEE=maxEE[1],
  #                        truncQ=truncQ,
  #                        trimRight = ifelse(any(is.na(trim.right)),0,trim.right[1]),
  #                        rm.phix=rm.phix,
  #                        compress=compress,
  #                        multithread=multithread)
  #   }
  #   saveRDS(out,"./output/QC_Stats.RDS")
    out <- readRDS("./output/QC_Stats.RDS")
  # In case some samples may have had zero reads pass QC, reassign filts
  # for each unique filtered path, go in, find all files, stick them together in one big vector,
  # then sort them to match the metadata sheet by finding any missing and removing those from the filts_f filts_r and samplenames

  if(paired){
    new_filts_f <- c()
    new_filts_r <- c()
    for(i in unique(filtered_paths)){
      x <- list.files(i,full.names = TRUE,pattern = "_filtered_fwd.fastq.gz")  
      y <- list.files(i,full.names = TRUE,pattern = "_filtered_rev.fastq.gz")  
      new_filts_f <- c(new_filts_f,x)
      new_filts_r <- c(new_filts_r,y)
    }
  } else {
    new_filts_f <- c()
    # new_filts_r <- c()
    for(i in unique(filtered_paths)){
      x <- list.files(i,full.names = TRUE,pattern = "_filtered_fwd.fastq.gz")  
      # y <- list.files(i,full.names = TRUE,pattern = "_filtered_rev.fastq.gz")  
      new_filts_f <- c(new_filts_f,x)
      # new_filts_r <- c(new_filts_r,x)
    }
  }
  
  
  # find any missing and pull them out
  filts_f <- filts_f[which((filts_f %in% new_filts_f))]
  if(paired){filts_r <- filts_r[which((filts_r %in% new_filts_r))]}
  
  # make metadata match files that made it through all previous steps!
  metadata <- meta[which((filts_f %in% new_filts_f)),]
  
  
  # learn errors
  errF <- learnErrors(filts_f, multithread=ifelse(multithread>1,TRUE,FALSE),
                      MAX_CONSIST = 10,verbose = 1,
                      randomize = TRUE) # set multithread = FALSE on Windows
  errF_out <- paste0(amplicon,"_err_Fwd.RDS")
  if(!dir.exists(asv.table.dir)){dir.create(asv.table.dir)}
  saveRDS(errF,file.path(asv.table.dir,errF_out))
  errF <- readRDS("./data/ASV_Tables/16S_V3V4_err_Fwd.RDS")
  
  if(paired){
  errR <- learnErrors(filts_r, multithread=ifelse(multithread>1,TRUE,FALSE),
                      MAX_CONSIST = 10,verbose = 1,
                      randomize = TRUE) # set multithread = FALSE on Windows
  errR_out <- paste0(amplicon,"_err_Rev.RDS")
  saveRDS(errR,file.path(asv.table.dir,errR_out))
  }
  errR <- readRDS("./data/ASV_Tables/16S_V3V4_err_Rev.RDS")  
  
  # dereplication
  
  derepF <- derepFastq(filts_f, verbose=TRUE)
  saveRDS(derepF,file.path(asv.table.dir,paste0(amplicon,"_derep_Fwd.RDS")))
  if(paired){
    derepR <- derepFastq(filts_r, verbose=TRUE)
    saveRDS(derepR,file.path(asv.table.dir,paste0(amplicon,"_derep_Rev.RDS")))
    }
  
  # SAMPLE INFERRENCE ####
  dadaFs <- dada(derepF, err=errF, multithread=ifelse(multithread>1,TRUE,FALSE), 
                 selfConsist = FALSE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows
  saveRDS(dadaFs,file.path(asv.table.dir,paste0(amplicon,"_dada_Fwd.RDS")))
  
  if(paired){
    dadaRs <- dada(derepR, err=errR, multithread=ifelse(multithread>1,TRUE,FALSE), 
                   selfConsist = FALSE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows
    saveRDS(dadaRs,file.path(asv.table.dir,paste0(amplicon,"_dada_Rev.RDS")))
  }

    
  # MERGE FWD and REV READS ####
  if(paired){
    mergers <- mergePairs(dadaFs, filts_f, dadaRs, filts_r, verbose=FALSE)
  } 
    
  # MAKE SEQUENCE TABLE ####
  if(paired){
    seqtab <- makeSequenceTable(mergers)
  } else {
    seqtab <- makeSequenceTable(dadaFs)
    }
  
  
  # REMOVE CHIMERAS ####
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=ifelse(multithread>1,TRUE,FALSE), verbose=TRUE)
  saveRDS(seqtab.nochim,"./data/ASV_Tables/16S_V3V4_seqtab_nochim.RDS")
  # TRACK READS ####
  getN <- function(x) sum(getUniques(x))
  if(nrow(out) > 1){
    if(paired){
      track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
      colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
      rownames(track) <- sample_names
    } else {
      track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
      colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
      rownames(track) <- sample_names
    }
    
    # make name for read tracking output
    track.out <- paste0(asv.table.dir,"/",amplicon,"_trackreads.RDS")
    # export tracking results
    saveRDS(track,track.out)
    
  } else {
    warning("Warning: Only one sample is present in the seq table. Sequence stat tracking disabled. Something is very wrong.")
  }
    
  
  
  
  # REMOVE CONTAMINANTS ####
  
  # find negative control samples, if any
  metadata[["control"]] <- metadata[[control.col]] == control.indicator
  
  
  # only run if there are negative control(s) that have at least some reads
  if(any(metadata[["control"]]) & sum(seqtab.nochim[metadata[["control"]],]) > 0){
    # Find and remove contaminants
    contams = decontam::isContaminant(seqtab.nochim, neg = metadata$control, normalize = TRUE)

    # remove contaminant sequences and control samples from both tables, respectively ####
    seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
    seqtab.nochim = seqtab.nochim[!metadata$control,]
    print(paste0(sum(contams$contaminant)," likely contaminants were removed from the ASV table of Run ",as.character(run.id),"."))
  }
  
  # make output name for ASV table
  asv_out <- paste0(asv.table.dir,"/",amplicon,"_ASV_Table.RDS")
  
  saveRDS(seqtab.nochim,asv_out)

}



# assign_taxonomy_to_asv_table()
# function to assign taxonomy to an asv table from the dada2 pipeline
# inputs: asv table | path-to-database | 

assign_taxonomy_to_asv_table <- function(asv.table, # asv table object name
                                         tax.database, # path to taxonomic database (fasta)
                                         multithread=(parallel::detectCores()-1), # set to FALSE on Windows
                                         random.seed=666,
                                         try.rc = TRUE, # attempt revComplement assignments as well? (doubles time)
                                         min.boot=50 # bootstrap of 50% recommended for seqs shorter than 250nt
                                         ){
  x <- assignTaxonomy(seqs = asv.table,
                      refFasta = tax.database,
                      minBoot = min.boot,
                      tryRC = try.rc,
                      outputBootstraps = FALSE,
                      multithread = multithread,
                      verbose = FALSE)
  
  return(x)
  
}


# find_gps_dists() ####
# given two data.frames (x,y) of lat/lon points, finds the distance from each point 
# in x to each point in y. By default, returns the distance to only the closest point in y
# good for finding nearest points
find_gps_dists <- 
  function(points1,points2,min.only=TRUE){
    
    # tests
    stopifnot(any(class(points1) == "data.frame"))
    stopifnot(any(class(points2) == "data.frame"))
    
    if(ncol(points1) != 2 | ncol(points2) != 2){
      stop("data frames must have 2 columns only")
    }
    
    if(!apply(points1,2,class) %>% unique() %in% c("numeric","integer")){
      stop("columns must be numeric; col1=longitude,col2=latitude")
    }
    
    if(!apply(points2,2,class) %>% unique() %in% c("numeric","integer")){
      stop("columns must be numeric; col1=longitude,col2=latitude")
    }
    
    # actual function
    mylist <- list()
    for(i in 1:nrow(points1)){
      mylist[[i]] <- points1[i,] %>% unlist
    }
    
    distances <- list()
    for(i in 1:nrow(points2)){
      
      mydistfunction <- function(x){geosphere::distHaversine(x,points2[i,])}
      colname <- paste0("dist_to_",i)
      distances[[colname]] <- map_dbl(mylist, mydistfunction)
    }
    
    x <- as.data.frame(distances)
    mins <- apply(x,1,min)
    
    if(min.only){
      return(mins)
    } else {
      return(x)
    }
    
  }


# auto_permanova() ####
auto_permanova <- function(physeq,
                            data,
                            pred.cols,
                            strata.col,
                            mod.type = "additive"){
  
  ra_comm <- 
    physeq %>% 
    transform_sample_counts(ra) %>% 
    otu_table() %>% 
    as.matrix()
  
  cc <- data %>% 
    dplyr::select(all_of(c(pred.cols,strata.col))) %>% 
    complete.cases()
  
  df <- data[cc,]

  if(mod.type == "additive"){
    mod.formula <- as.formula(paste0("ra_comm[cc,]","~",paste(pred.cols,collapse=" + ")))
  }
  
  if(mod.type == "interactive"){
    mod.formula <- as.formula(paste0("ra_comm[cc,]","~",paste(pred.cols,collapse=" * ")))
  }
  
  
  # run simple permanova
  if(!is.na(strata.col)){
    mod <- 
      adonis2(data = df,
              formula = mod.formula,
              strata = df[[strata.col]])
  } else {
    mod <- 
      adonis2(data = df,
              formula = mod.formula)
  }
  
  
 return(mod) 
}


# clean_ps_taxonomy() ####
# get rid of the annoying "k__" stuff at the beginning of taxonomy assignments for each tax level
# these are usually only an issue with fungal assignments (e.g., UNITE format)


clean_ps_taxonomy <- function(physeq,
                              n.ranks=7 # number of taxonomic ranks in physeq object
                              ){
  
  # get rank names
  ranks <- rank_names(physeq)
  prefix <- str_sub(ranks,end=1) %>% str_to_lower() %>% paste0("__")
  
  x <- tax_table(physeq)@.Data
  
  for(i in seq_along(ranks)){
    x[,i] <- x[,i] %>% str_remove(prefix[i]) %>% str_remove("\\(([^)]*)\\)")
  }
  
  out <- phyloseq(otu_table(physeq,taxa_are_rows = FALSE),
                  tax_table(x),
                  sample_data(physeq))
  
  return(out)
}


# simplify_fungal_guilds() ####
# extract major guild groupings
# needs a data.frame with a "Guild" column that contains the results from FunGuild assignment
# will return the data.frame with a new column "major_guild"
simplify_fungal_guilds <- 
  function(x){
    x %>% 
      mutate(major_guild = case_when(grepl("Orchid Mycorrhizal",Guild,ignore.case = TRUE) ~ "Orchid Mycorrhizal",
                                     grepl("Ectomycorrhizal",Guild,ignore.case = TRUE) ~ "Ectomycorrhizal",
                                     grepl("ericoid",Guild,ignore.case = TRUE) ~ "Ericoid mycorrhizal",
                                     grepl("arbuscular",Guild,ignore.case=TRUE) ~ "Arbuscular mycorrhizal",
                                     grepl("Plant Pathogen",Guild,ignore.case=TRUE) &
                                       !grepl("mycorrhizal",Guild,ignore.case=TRUE) ~ "Plant pathogen",
                                     grepl("Animal Pathogen",Guild,ignore.case=TRUE) &
                                       !grepl("mycorrhizal",Guild,ignore.case=TRUE) ~ "Animal pathogen",
                                     grepl("Saprotroph",Guild,ignore.case=TRUE) &
                                       !grepl("mycorrhizal",Guild,ignore.case=TRUE) &
                                       !grepl("pathogen",Guild,ignore.case=TRUE) ~ "Saprotroph",
                                     grepl("lichenized",Guild,ignore.case=TRUE) &
                                       !grepl("mycorrhizal",Guild,ignore.case=TRUE) ~ "Lichenized",
                                     grepl("Animal Parasite",Guild,ignore.case = TRUE) &
                                       !grepl("mycorrhizal",Guild,ignore.case=TRUE) ~ "Animal Parasite",
                                     grepl("Algal Parasite|Plant Parasite",Guild,ignore.case = TRUE) &
                                       !grepl("mycorrhizal",Guild,ignore.case=TRUE) ~ "Plant Parasite",
                                     grepl("Animal Symbiotroph",Guild,ignore.case = TRUE) &
                                       !grepl("mycorrhizal",Guild,ignore.case=TRUE) ~ "Animal Symbiotroph"
      ))
  }

# googlemap_json_to_string()
# convert json google map styling to api string
googlemap_json_to_string <- 
  function (style_list) 
  {
    style_string <- ""
    for (i in 1:length(style_list)) {
      if ("featureType" %in% names(style_list[[i]])) {
        style_string <- paste0(style_string, "feature:", 
                               style_list[[i]]$featureType, "|")
      }
      if ("elementType" %in% names(style_list[[i]])) {
        style_string <- paste0(style_string, "element:", 
                               style_list[[i]]$elementType, "|")
      }
      elements <- style_list[[i]]$stylers
      a <- lapply(elements, function(x) paste0(names(x), ":", 
                                               x)) %>% unlist() %>% paste0(collapse = "|")
      style_string <- paste0(style_string, a)
      if (i < length(style_list)) {
        style_string <- paste0(style_string, "&style=")
      }
    }
    style_string <- gsub("#", "0x", style_string)
    style_string <- gsub("[|]", "%7C", style_string)
    return(style_string)
  }

# convert SINTAX taxonomic output to table form
parse_sintax_vsearch <- function(path, asv_order,
                                 ranks = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                 cutoffs = c(Kingdom=.5, Phylum=.6, Class=.7, Order=.7, Family=.8, Genus=.85, Species=.9)) {
  # Map any rank code (lower/upper) to our rank names
  rank_map <- c(
    d = "Kingdom", k = "Kingdom",
    p = "Phylum",
    c = "Class",
    o = "Order",
    f = "Family",
    g = "Genus",
    s = "Species",
    D = "Kingdom", K = "Kingdom",
    P = "Phylum",
    C = "Class",
    O = "Order",
    F = "Family",
    G = "Genus",
    S = "Species"
  )
  
  # Read TSV with unknown number of columns, no names
  raw <- read_tsv(path, col_names = FALSE, show_col_types = FALSE, comment = "#")
  stopifnot(ncol(raw) >= 2)
  
  # Heuristically pick the column that contains rank tokens (k:/d:/p:/...)
  has_rank_tokens <- function(x) any(str_detect(x, "(^|[,\\s])(tax=)?[dkpcfgsoDKPCFGSO]:"))
  cand_cols <- which(vapply(raw, has_rank_tokens, logical(1)))
  stopifnot(length(cand_cols) >= 1)
  sintax_col <- cand_cols[1]             # first column with tokens
  qseqid_col <- 1                        # first column is your header (the ASV sequence)
  
  df <- tibble(
    qseqid = raw[[qseqid_col]],
    sintax = raw[[sintax_col]]
  )
  
  # Normalize to a common token list:
  # - remove optional 'tax=' prefix
  # - split on commas or whitespace between tokens
  # - tokens look like: 'd:Name(0.97)' OR 'k:Name,0.97'
  tokens <- df %>%
    mutate(s = str_remove(sintax, "^tax=")) %>%
    separate_rows(s, sep = "[,\\s]+") %>%              # split on comma or whitespace
    filter(str_detect(s, "^[dkpcfgsoDKPCFGSO]:")) %>% # keep rank tokens only
    mutate(code = str_sub(s, 1, 1),
           # taxon: after "code:" up to first '(' or ',' (whichever comes first)
           taxon = str_remove(str_extract(s, ":[^\\(,]+"), "^:"),
           # confidence: either inside (...) or after a trailing comma
           conf  = coalesce(
             suppressWarnings(as.numeric(str_match(s, "\\(([0-9.]+)\\)")[,2])),
             suppressWarnings(as.numeric(str_match(s, ",([0-9.]+)$")[,2]))
           ),
           rank = unname(rank_map[code])) %>%
    filter(!is.na(rank)) %>%
    transmute(qseqid = df$qseqid[match(qseqid, df$qseqid)], rank, taxon, conf)
  
  # Choose top per-rank (usually one) and apply per-rank cutoffs
  tax_wide <- tokens %>%
    group_by(qseqid, rank) %>%
    slice_max(order_by = conf, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(taxon = ifelse(!is.na(conf) & conf >= cutoffs[rank], taxon, NA_character_)) %>%
    select(qseqid, rank, taxon) %>%
    distinct(qseqid, rank, .keep_all = TRUE) %>%
    pivot_wider(names_from = rank, values_from = taxon)
  
  # Ensure all expected rank columns exist
  for (r in ranks) if (!r %in% names(tax_wide)) tax_wide[[r]] <- NA_character_
  
  # Left-join to full ASV list to preserve order and include missing classifications
  full <- tibble(qseqid = asv_order) %>% left_join(tax_wide, by = "qseqid")
  
  M <- as.matrix(full[, ranks])
  rownames(M) <- full$qseqid
  M
}

