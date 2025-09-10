read_4DNdata <- function(file){
  # read files downloaded from 4DN datasets
  lines <- readLines(file)
  batch_size <- 50
  step <- batch_size
  total_lines <- length(lines)
  last_hash_row <- 0
  start <- 1
  while (step <= total_lines) { 
    matches <- which(grepl("^\\s*[^A-Za-z0-9]", lines[start:step]))
    if (length(matches) > 0) {
      last_hash_row <- max(matches)  # Get the max row number within the current step
      if (step == total_lines) break
      start <- step + 1
      step <- min(step + batch_size, total_lines)
    }else{break}
  }
  column_names <- trimws(strsplit(sub("^##", "", lines[last_hash_row]), ",")[[1]])
  column_names[1] <- sub("^columns=\\(", "", column_names[1])
  column_names[length(column_names)] <- sub("\\)$", "", column_names[length(column_names)])
  filtered_lines <- lines[(last_hash_row + 1):length(lines)]
  temp_file <- tempfile()
  writeLines(filtered_lines, temp_file)
  data <- data.table::fread(temp_file, header = FALSE)
  colnames(data) <- column_names
  #colnames(data) <- c("Spot_ID", "Trace_ID", "X", "Y", "Z", "Chrom", "Chrom_Start", "Chrom_End")
  outdata <- data.frame(Trace_ID = as.character(data$Trace_ID),
                        X = as.numeric(data$X), Y = as.numeric((data$Y)), Z = as.numeric(data$Z),
                        Chrom = as.character(data$Chrom), Chrom_Start = as.numeric(data$Chrom_Start), Chrom_End = as.numeric(data$Chrom_End))
  file.remove(temp_file)
  return(outdata)
}


read_csvdataframe <- function(file){
  data <- read.csv(file, header = TRUE)
  if(all(c("Trace_ID", "X", "Y", "Z", "Chrom", "Chrom_Start", "Chrom_End") %in% colnames(data))){
    outdata <- data.frame(Trace_ID = as.character(data$Trace_ID),
                          X = as.numeric(data$X), Y = as.numeric((data$Y)), Z = as.numeric(data$Z),
                          Chrom = as.character(data$Chrom), Chrom_Start = as.numeric(data$Chrom_Start), Chrom_End = as.numeric(data$Chrom_End))
  }else if(all(c("Trace_ID", "X", "Y", "Z", "Region_ID") %in% colnames(data))){
    outdata <- data.frame(Trace_ID = as.character(data$Trace_ID), Region_ID = as.numeric(data$Region_ID),
                          X = as.numeric(data$X), Y = as.numeric((data$Y)), Z = as.numeric(data$Z))
  }else{
    stop("The input dataset is missing required columns or column names. Ensure that the dataset contains one of the followging sets of columns (with column names):
         1. 'Trace_ID', 'X', 'Y', 'Z', 'Chrom', 'Chrom_Start', 'Chrom_End', or
         2. 'Trace_ID', 'X', 'Y', 'Z', 'Region_ID'.
         Please check your dataset (column names) and try again.")
  }
  return(outdata)
}

call_region <- function(dataset){
  if (!all(c("Chrom", "Chrom_Start", "Chrom_End") %in% colnames(dataset))) {
    stop("The input dataset must contain 'Chrom', 'Chrom_Start', and 'Chrom_End' columns.")
  }
  region_in <- data.frame(Chrom=as.character(dataset$Chrom), Chrom_Start=as.integer(dataset$Chrom_Start), Chrom_End=as.integer(dataset$Chrom_End))
  region_1 <- dplyr::distinct(region_in, Chrom, Chrom_Start, Chrom_End)
  region_2 <- region_1[order(region_1[,1], as.numeric(region_1[,2]), as.numeric(region_1[,3])),]
  region_out <- data.frame(region_2, region = paste(region_2[,1], region_2[,2], region_2[,3], sep = "_"), Region_ID=seq(1,dim(region_2)[1]))
  return(region_out)
}


calculate_euclidean_distance <- function(dataset, region_num){
  dataset <- dataset[which(dataset$X != 0 & dataset$Y != 0 & dataset$Z != 0), ]
  coords <- dataset
  traces <- sort(unique(coords$Trace_ID))
  dismat <- matrix(ncol=length(traces), nrow=(region_num*(region_num-1))/2) # region by trace matrix
  for(n in 1:length(traces)){
    tccoord <- coords[which(coords$Trace_ID==traces[n]),]
    i <- 1
    m <- 1
    while(i < region_num){
      j <- i+1
      while(j <= region_num){
        if(i %in% tccoord$Region_ID & j %in% tccoord$Region_ID){
          dismat[m,n] <- sqrt((tccoord$X[which(tccoord$Region_ID==i)]-tccoord$X[which(tccoord$Region_ID==j)])^2
                              +(tccoord$Y[which(tccoord$Region_ID==i)]-tccoord$Y[which(tccoord$Region_ID==j)])^2
                              +(tccoord$Z[which(tccoord$Region_ID==i)]-tccoord$Z[which(tccoord$Region_ID==j)])^2)
        }
        j <- j+1
        m <- m+1
      }
      i <- i+1
    }
  }
  return(dismat)
}

WilcoxonP <- function(k, dismat1, dismat2){
  if(length(which(!is.na(as.numeric(dismat1[k,]))))>0 & 
     length(which(!is.na(as.numeric(dismat2[k,]))))>0){
    pgreat <- wilcox.test(as.numeric(dismat1[k,]), as.numeric(dismat2[k,]), paired=FALSE, alternative=c("greater"))$p.value
    pless <- wilcox.test(as.numeric(dismat1[k,]), as.numeric(dismat2[k,]), paired=FALSE, alternative=c("less"))$p.value
  }else{
    pgreat <- NA
    pless <- NA
  }
  wp <- data.frame(pgreat = pgreat, pless = pless)
  return(wp)
}

WilcoxonPMat <- function(wpboth){
  region_num <- (1 + sqrt(1 + 8 * dim(wpboth)[1])) / 2
  mat <- matrix(ncol = region_num, nrow=region_num)
  i <- 1
  m <- 1
  while(i < region_num){
    j <- i + 1
    while(j <= region_num){
      if(!is.na(wpboth$pgreat[m]) & !is.na(wpboth$pless[m])){
        a <- c(-log10(wpboth$pgreat[m]), log10(wpboth$pless[m]))
        mat[i,j] <- a[which.max(abs(a))]
        mat[j,i] <- a[which.max(abs(a))]
      }
      j <- j + 1
      m <- m + 1
    }
    i <- i + 1
  }
  return(mat)
} 

PairCellTypeP <- function(n, opt, dismat1, dismat2, objsmat){
  ct <- as.matrix(objsmat[n,])
  pr_num <- dim(dismat1)[1]
  wpboth <- furrr::future_map_dfr(1:pr_num, function(k) WilcoxonP(k, dismat1, dismat2))
  wpmat <- WilcoxonPMat(wpboth)
  write.table(wpmat, file=paste0(opt$res2path, '/DiffScoreMatrix_', ct[1], "VS", ct[2],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  # plot
  colnames(wpmat) <- paste0("region", seq(1:dim(wpmat)[1]))
  rownames(wpmat) <- paste0("region", seq(1:dim(wpmat)[1]))
  data <- reshape2::melt(dplyr::mutate(as.data.frame(wpmat), index=row.names(wpmat)), id="index")
  colnames(data) <- c("T1", "T2", "value")
  data$T1 <- factor(data$T1, levels=paste0("region", seq(1,dim(wpmat)[1],1)))
  data$T2 <- factor(data$T2, levels=paste0("region", seq(dim(wpmat)[1],1,-1)))
  pdf(paste0(opt$res2path,'/DiffScoreHeatmap_', ct[1], "VS", ct[2],'.pdf'))
  p <- ggplot(data,aes(x=T1,y=T2,fill=value))+ 
    scale_fill_gradient2(low="#d60b0e", mid="white", high="#0f70bf", midpoint = 0, na.value = "grey", limits=c(-max(abs(data$value), na.rm = TRUE), max(abs(data$value), na.rm = TRUE)))+
    geom_raster()+
    labs(fill="DiffScore", title = paste0(ct[1], " VS ", ct[2]),
        x="Region ID", y="Region ID")+
    scale_x_discrete(breaks = paste0("region", seq(5, pr_num, by = 5)), labels=seq(5, pr_num, by=5))+
    scale_y_discrete(breaks = paste0("region", seq(5, pr_num, by = 5)), labels=seq(5, pr_num, by=5))
  print(p)
  dev.off()
  cat("complete the ", n, "th pairs of cell types. \n", file = logfile, append = TRUE)
}


process_EachChr <- function(chrind,chrset, foredataall, backdataall, opt, logfile, outpath, objs){
  opt$res0path <- paste0(outpath, "/S0_DataInfo/chr_", chrset[chrind])
  opt$res1path <- paste0(outpath, "/S1_Distance/chr_", chrset[chrind])
  opt$res2path <- paste0(outpath, "/S2_DiffScore/chr_", chrset[chrind])
  if(!dir.exists(opt$res0path)){
    dir.create(opt$res0path, recursive = TRUE)
  }else{
    cat("Step0 directory already exists! \n", file = stdout())
  }
  if(!dir.exists(opt$res1path)){
    dir.create(opt$res1path, recursive = TRUE)
  }else{
    cat("Step1 directory already exists! \n", file = stdout())
  }
  if(!dir.exists(opt$res2path)){
    dir.create(opt$res2path, recursive = TRUE)
  }else{
    cat("Step2 directory already exists! \n", file = stdout())
  }
  foredata <- foredataall[which(foredataall$Chrom == chrset[chrind]), ]
  backdata <- backdataall[which(backdataall$Chrom == chrset[chrind]), ]
  ## step0 output: genomic region and region ID table; trace count in each state fig.
  fbdata <- c("foredata", "backdata")
  tracesallct <- c()
  for(i in 1:length(fbdata)){
    indata0 <- get(fbdata[i])
    # trace
    tracesallct <- rbind(tracesallct, length(sort(unique(indata0$Trace_ID))))
    # region
    if (!all(c("Region_ID") %in% colnames(indata0))) {
      region_out <- call_region(indata0)
      assign(fbdata[i], dplyr::left_join(indata0, region_out, by=c("Chrom","Chrom_Start", "Chrom_End")))
    }else{
      region_out <- data.frame(Region_ID = sort(unique(indata0$Region_ID)))
    }
    assign(paste0(fbdata[i], 'rg'), region_out)
  }
  rownames(foredatarg) <- NULL
  rownames(backdatarg) <- NULL
  if(identical(foredatarg, backdatarg)){
    region_out <- foredatarg
    region_num <- dim(region_out)[1]
    write.table(region_out, file=paste0(opt$res0path, '/RegionIn', chrind, '.tsv'), sep="\t", quote=FALSE, row.names = FALSE)
  }else{
    stop("The genomic regions in these two celltypes/states/clusters are not the same")
  }
  dataplot <-data.frame(celltype = objs, tracecount = tracesallct)
  pdf(paste0(opt$res0path,'/TraceCount', chrind, '.pdf'))
  pct <- ggplot(dataplot, aes(x = celltype, y = tracecount)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_bw() +
    labs(title = chrind, x = "", y = "Trace count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(pct)
  dev.off()
  ## step1 output: region paris' distance
  dismat1 <- calculate_euclidean_distance(foredata, region_num)
  dismat2 <- calculate_euclidean_distance(backdata, region_num) 
  write.table(dismat1, file=paste0(opt$res1path, '/RegionPairsByTrace_', chrind, objs[1],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(dismat2, file=paste0(opt$res1path, '/RegionPairsByTrace_', chrind, objs[2],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  cat(chrset[chrind], "complete the step1 'Distance Calculation'. \n", file = logfile, append = TRUE)
  ## step2 output: DiffScore matrix
  objsmat <- matrix(c(objs[1], objs[2]), nrow=1, ncol=2)
  pair_celltype_p <- PairCellTypeP(1, opt, dismat1, dismat2, objsmat)
  cat(chrset[chrind], "complete the step2 'DiffScore Generation'. \n", file = logfile, append = TRUE)
}




time <- Sys.time()
future::plan("multisession")
options(future.globals.maxSize = 1024 * 1024 * 1024)  # 1 GB

### read command line 
option_list <- list(
  optparse::make_option(c("-a", "--inputPath1"), type = "character", help = "Required. The foreground input file"),
  optparse::make_option(c("-b", "--inputPath2"), type = "character", help = "Required. The background input file"),
  optparse::make_option(c("-o", "--outputPath"), type = "character", help = "Optional. Path to output directory"),
  optparse::make_option(c("-m", "--fileformat"), type = "character", help = "Required. Format of input data ['4DN'] or ['csv']"),
  optparse::make_option(c("-c", "--chrnum"), type = "character", help = "Optional. Options: 'one' or 'more'. 'one' means all regions are located on one chromosome; 'more' means regions span multiple chromosomes, typically across the whole genome.")

)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
### creat output directory for step0: data information inlcuding genomic regions and trace count in each state; step1: region paris' distance; and step2: DiffScore matrix.
outpath <- opt$outputPath
if(substr(outpath, nchar(outpath), nchar(outpath)) == "/"){
  outpath <- substr(outpath, 1, nchar(outpath)-1)
}
## read input data
library(ggplot2)
logfile <- paste0("Dory_logfile_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")

listsdata <- data.frame(strsplit(opt$inputPath1, "[._/]"), strsplit(opt$inputPath2, "[._/]")) 
objs <- as.character(t(listsdata[dim(listsdata)[1]-1,])) ## the object name is the file name without '.csv'

if(opt$fileformat == '4DN'){
  foredataall <- read_4DNdata(opt$inputPath1)
  backdataall <- read_4DNdata(opt$inputPath2)
}else if(opt$fileformat == 'csv'){
  foredataall <- read_csvdataframe(opt$inputPath1)
  backdataall <- read_csvdataframe(opt$inputPath2)
}else{
  stop("Please indicate the data format of input files: '-m 4DN' or '-m csv' ")
}
chrset <- unique(sort(c(foredataall$Chrom, backdataall$Chrom)))

if(is.null(opt$chrnum)|| opt$chrnum == ""){
  if(length(chrset) == 1){
    opt$chrnum <- 'one'
  }else if(length(chrset) > 1){
    opt$chrnum <- 'more'
  }else{
    stop("No chromosome information detected. Please verify the 'Chrom' column in the coordinate file.
    If the input data contains only region IDs and not chromosome names, please indicate whether these regions belong to one chromosome using '-c one', or span more chromosomes by '-c more'.")
  }
}

if(opt$chrnum == 'one'){
  opt$res0path <- paste0(outpath, "/S0_DataInfo/")
  opt$res1path <- paste0(outpath, "/S1_Distance")
  opt$res2path <- paste0(outpath, "/S2_DiffScore")
  if(!dir.exists(opt$res0path)){
    dir.create(opt$res0path, recursive = TRUE)
  }else{
    cat("Step0 directory already exists! \n", file = stdout())
  }
  if(!dir.exists(opt$res1path)){
    dir.create(opt$res1path)
  }else{
    cat("Step1 directory already exists! \n", file = stdout())
  }
  if(!dir.exists(opt$res2path)){
    dir.create(opt$res2path)
  }else{
    cat("Step2 directory already exists! \n", file = stdout())
  }
  foredata <- foredataall
  backdata <- backdataall
  ## step0 output: genomic region and region ID table; trace count in each state fig.
  fbdata <- c("foredata", "backdata")
  tracesallct <- c()
  for(i in 1:length(fbdata)){
    indata0 <- get(fbdata[i])
    # trace
    tracesallct <- rbind(tracesallct, length(sort(unique(indata0$Trace_ID))))
    # region
    if (!all(c("Region_ID") %in% colnames(indata0))) {
      region_out <- call_region(indata0)
      assign(fbdata[i], dplyr::left_join(indata0, region_out, by=c("Chrom","Chrom_Start", "Chrom_End")))
    }else{
      region_out <- data.frame(Region_ID = sort(unique(indata0$Region_ID)))
    }
    assign(paste0(fbdata[i], 'rg'), region_out)
  }
  rownames(foredatarg) <- NULL
  rownames(backdatarg) <- NULL
  if(identical(foredatarg, backdatarg)){
    region_out <- foredatarg
    region_num <- dim(region_out)[1]
    write.table(region_out, file=paste0(opt$res0path, '/Region.tsv'), sep="\t", quote=FALSE, row.names = FALSE)
  }else{
    stop("The genomic regions in these two celltypes/states/clusters are not the same")
  }
  dataplot <-data.frame(celltype = objs, tracecount = tracesallct)
  pdf(paste0(opt$res0path,'/TraceCount.pdf'))
  pct <- ggplot(dataplot, aes(x = celltype, y = tracecount)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_bw() +
    labs(x = "", y = "Trace count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(pct)
  dev.off()
  ## step1 output: region paris' distance
  dismat1 <- calculate_euclidean_distance(foredata, region_num)
  dismat2 <- calculate_euclidean_distance(backdata, region_num)
  write.table(dismat1, file=paste0(opt$res1path, '/RegionPairsByTrace_', objs[1],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(dismat2, file=paste0(opt$res1path, '/RegionPairsByTrace_', objs[2],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  cat("complete the step1 'Distance Calculation'. \n", file = logfile, append = TRUE)
  cat("time: ", as.numeric(Sys.time() - time), attr(Sys.time() - time, "units"), "\n", file = logfile, append = TRUE)
  ## step2 output: DiffScore matrix
  objsmat <- matrix(c(objs[1], objs[2]), nrow=1, ncol=2)
  pair_celltype_p <- PairCellTypeP(1, opt, dismat1, dismat2,  objsmat)
  cat("complete the step2 'DiffScore Generation'. \n", file = logfile, append = TRUE)
}else if(opt$chrnum == 'more'){
  ### for regions spanning more chromosomes, typically across the whole genome
  allchrsout <- furrr::future_map_dfr(1:length(chrset), function(n) process_EachChr(n, chrset, foredataall, backdataall, opt, logfile, outpath, objs), .options = furrr::furrr_options(seed = TRUE))
}else{
    stop("Please indicate the chromosome number: '-c one' or '-c more' ")
}

cat("time: ", as.numeric(Sys.time() - time), attr(Sys.time() - time, "units"), "\n", file = logfile, append = TRUE)

