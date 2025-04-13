read_4DNgroup <- function(coordfile, cellfile, ctcolname){
  if(is.null(opt$coordfile)|| opt$coordfile == ""){
    stop("No coordinate file detected. Please input it with '-i inputCoordFile'.")
  }
  if(is.null(opt$cellfile)|| opt$cellfile == ""){
    stop("No cell label file detected. Please input it with '-l labelFile'.")
  }
  # read files downloaded from 4DN datasets
  lines <- readLines(coordfile)
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
  column_names <- strsplit(sub("^##", "", lines[last_hash_row]), ",")[[1]]
  column_names[1] <- sub("^columns=\\(", "", column_names[1])
  column_names[length(column_names)] <- sub("\\)$", "", column_names[length(column_names)])
  filtered_lines <- lines[(last_hash_row + 1):length(lines)]
  temp_file <- tempfile()
  writeLines(filtered_lines, temp_file)
  data <- data.table::fread(temp_file, header = FALSE)
  colnames(data) <- column_names
  outdata <- data.frame(Trace_ID = as.character(data$Trace_ID),
                        X = as.numeric(data$X), Y = as.numeric((data$Y)), Z = as.numeric(data$Z),
                        Chrom = as.character(data$Chrom), Chrom_Start = as.numeric(data$Chrom_Start), Chrom_End = as.numeric(data$Chrom_End),
                        Cell_ID = as.numeric(data$Cell_ID))
  ct_lines <- readLines(cellfile)
  ct_hash_lines <- grep("^##", ct_lines)
  if (length(ct_hash_lines) > 0){
    ct_last_hash_row <- max(grep("^\\s*[^A-Za-z0-9]", ct_lines))
    ct_column_names <- strsplit(sub("^##", "", ct_lines[ct_last_hash_row]), ",")[[1]]
    ct_column_names[1] <- sub("^columns=\\(", "", ct_column_names[1])
    ct_column_names[length(ct_column_names)] <- sub("\\)$", "", ct_column_names[length(ct_column_names)])
    ct_filtered_lines <- ct_lines[(ct_last_hash_row + 1):length(ct_lines)]
    data_ct <- read.csv(text = paste(ct_filtered_lines, collapse = "\n"), header = FALSE, col.names = ct_column_names, sep = ",")
  }else{
    data_ct <- read.csv(cellfile, header = TRUE)
  }
  if (ctcolname %in% colnames(data_ct)) {
    if(!("Cell_ID" %in% colnames(data_ct))){
      stop("'Cell_ID' not found in 'cellfile'. Please indicate the 'Cell_ID' column.")
    }
    cell_label <- data_ct[, c("Cell_ID", ctcolname)]
  } else {
    stop("No cell label column found in 'cellfile'. Please specify it using '-c ColNameOfCellLabel'.
    Although '-c' is optional, its default value is 'Cell_Type'. If your cells are labeled under a different column name, please provide that column name explicitly in the command line.")
  }
  colnames(cell_label) <- c("Cell_ID", "Cell_Type")
  outcelltype <- data.frame(Cell_ID = as.numeric(cell_label$Cell_ID), Cell_Type = as.character(cell_label$Cell_Type))
  outdata_celltype <- dplyr::left_join(outdata, cell_label, by=c("Cell_ID"))
  file.remove(temp_file)
  return(outdata_celltype)
}


read_csvgroup <- function(coordfile, cellfile, ctcolname){
  if(is.null(opt$coordfile)|| opt$coordfile == ""){
    stop("No coordinate file detected. Please input it with '-i inputCoordFile'.")
  }
  if(is.null(opt$cellfile)|| opt$cellfile == ""){
    stop("No cell label file detected. Please input it with '-l labelFile'.")
  }
  data <- read.csv(coordfile, header = TRUE)
  if(all(c("Trace_ID", "X", "Y", "Z", "Chrom", "Chrom_Start", "Chrom_End", "Cell_ID") %in% colnames(data))){
    outdata <- data.frame(Trace_ID = as.numeric(data$Trace_ID),
                          X = as.numeric(data$X), Y = as.numeric((data$Y)), Z = as.numeric(data$Z),
                          Chrom = as.character(data$Chrom), Chrom_Start = as.numeric(data$Chrom_Start), Chrom_End = as.numeric(data$Chrom_End),
                          Cell_ID = as.numeric(data$Cell_ID))
  }else if(all(c("Trace_ID", "X", "Y", "Z", "Region_ID", "Cell_ID") %in% colnames(data))){
    outdata <- data.frame(Trace_ID = as.numeric(data$Trace_ID), Region_ID = as.numeric(data$Region_ID),
                          X = as.numeric(data$X), Y = as.numeric((data$Y)), Z = as.numeric(data$Z),
                          Cell_ID = as.numeric(data$Cell_ID))
  }else{
    stop("The input dataset is missing required columns or column names. Ensure that the dataset contains one of the followging sets of columns (with column names):
         1. 'Trace_ID', 'X', 'Y', 'Z', 'Chrom', 'Chrom_Start', 'Chrom_End', 'Cell_ID' or
         2. 'Trace_ID', 'X', 'Y', 'Z', 'Region_ID', 'Cell_ID'.
         Please check your dataset (column names) and try again.")
  }
  cell_label <- read.csv(cellfile, header = TRUE)
  if (ctcolname %in% colnames(cell_label)) {
      outcelltype <- data.frame(Cell_ID = as.numeric(cell_label$Cell_ID), Cell_Type = as.character(cell_label[ctcolname]))
  } else {
    stop("No cell label column found in 'cellfile'. Please specify it using '-c ColNameOfCellLabel'.
    Although '-c' is optional, its default value is 'Cell_Type'. If your cells are labeled under a different column name, please provide that column name explicitly in the command line.")
  }
  outdata_celltype <- dplyr::left_join(outdata, cell_label, by=c("Cell_ID"))
  return(outdata_celltype)
}

read_csvgroup_withCT <- function(coordfile, ctcolname){
  data <- read.csv(coordfile, header = TRUE)
  if(all(c("Trace_ID", "X", "Y", "Z", "Chrom", "Chrom_Start", "Chrom_End") %in% colnames(data))){
    outdata <- data.frame(Trace_ID = as.character(data$Trace_ID),
                          X = as.numeric(data$X), Y = as.numeric((data$Y)), Z = as.numeric(data$Z),
                          Chrom = as.character(data$Chrom), Chrom_Start = as.numeric(data$Chrom_Start), Chrom_End = as.numeric(data$Chrom_End))
  }else if(all(c("Trace_ID", "X", "Y", "Z", "Region_ID") %in% colnames(data))){
    outdata <- data.frame(Trace_ID = as.character(data$Trace_ID), Region_ID = as.character(data$Region_ID),
                          X = as.numeric(data$X), Y = as.numeric((data$Y)), Z = as.numeric(data$Z))
  }else{
    stop("If the input xxxxxx. The input dataset is missing required columns or column names. Ensure that the dataset contains one of the followging sets of columns (with column names):
         1. 'Trace_ID', 'X', 'Y', 'Z', 'Chrom', 'Chrom_Start', 'Chrom_End', 'Cell_Type' or
         2. 'Trace_ID', 'X', 'Y', 'Z', 'Region_ID', 'Cell_Type'.
         Please check your dataset (column names) and try again.")
  }
  if (ctcolname %in% colnames(data)) {
      outdata$Cell_Type <- as.character(data[ctcolname])
  } else {
    stop("No cell label column found in the input file. If the cell label is not included in the coordinate file, please input cell label file using '-l labelFile'. 
    If the cell label is included in the coordinate file, please specify it using '-c ColNameOfCellLabel'.
    Although '-c' is optional, its default value is 'Cell_Type'. If your cells are labeled under a different column name, please provide that column name explicitly in the command line.")
  }
  outdata_celltype <- outdata
  return(outdata_celltype)
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

WilcoxonP <- function(k, ct){
  if(length(which(!is.na(as.numeric(get(ct[1])[k,]))))>0 & 
     length(which(!is.na(as.numeric(get(ct[2])[k,]))))>0){
    pgreat <- wilcox.test(as.numeric(get(ct[1])[k,]), as.numeric(get(ct[2])[k,]), paired=FALSE, alternative=c("greater"))$p.value
    pless <- wilcox.test(as.numeric(get(ct[1])[k,]), as.numeric(get(ct[2])[k,]), paired=FALSE, alternative=c("less"))$p.value
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

PairCellType <- function(objs){
  num <- length(objs)
  objsmat <- matrix(ncol = 2, nrow = num * (num - 1)/2)
  i <- 1
  k <- 1
  while(i < num){
    j <- i + 1
    while(j <= num){
      objsmat[k,] <- c(objs[i], objs[j])
      j <- j + 1
      k <- k + 1
    }
    i <- i + 1
  }
  return(objsmat)
}

PairCellTypeP <- function(n, opt, objs, objsmat){
  library(ggplot2)
  ct <- objsmat[n,]
  ct <- as.matrix(ct)
  pr_num <- dim(get(ct[1]))[1]
  wpboth <- furrr::future_map_dfr(1:pr_num, function(k) WilcoxonP(k, ct))
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
        x="Region ID",y="Region ID")+
    #theme(axis.ticks = element_blank())+
    scale_x_discrete(breaks = paste0("region", seq(5, dim(wpmat)[1], by = 5)), labels=seq(5, dim(wpmat)[1], by=5))+
    scale_y_discrete(breaks = paste0("region", seq(5, dim(wpmat)[1], by = 5)), labels=seq(5, dim(wpmat)[1], by=5))
  print(p)
  dev.off()  
  cat("complete the ", n, "th pairs of cell types. \n", file = logfile, append = TRUE)
}


process_EachChr <- function(chrind, allchrs, opt, outpath){
  cat(allchrs[chrind], "time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = logfile, append = TRUE)
  library(ggplot2)
  opt$res0path <- paste0(outpath, "/S0_DataInfo/", allchrs[chrind])
  opt$res1path <- paste0(outpath, "/S1_Distance/", allchrs[chrind])
  opt$res2path <- paste0(outpath, "/S2_DiffScore/", allchrs[chrind])
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
  indata0 <- get(allchrs[chrind])
  if (!all(c("Region_ID") %in% colnames(indata0))) {
    region_out <- call_region(indata0)
    indata <- dplyr::left_join(indata0, region_out, by=c("Chrom","Chrom_Start", "Chrom_End"))
  }else{
    region_out <- data.frame(Region_ID = sort(unique(indata0$Region_ID)))
    indata <- indata0
  }
  region_num <- dim(region_out)[1]
  write.table(region_out, file=paste0(opt$res0path, '/RegionIn', allchrs[chrind], '.tsv'), sep="\t", quote=FALSE, row.names = FALSE)
  celltype <- unique(sort(indata$Cell_Type))
  objs <- paste0("CT_", celltype)
  tracesallct <- c()
  for(i in 1:length(objs)){
    assign(objs[i],  calculate_euclidean_distance(indata[which(indata$Cell_Type == celltype[i]), ], region_num), envir = .GlobalEnv)
    write.table(get(objs[i]), file=paste0(opt$res1path, '/RegionPairsByTrace_', objs[i],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
    tracesallct <- rbind(tracesallct, length(sort(unique(indata$Trace_ID[which(indata$Cell_Type == celltype[i])]))))
  }
  dataplot <-data.frame(celltype = celltype, tracecount = tracesallct)
  pdf(paste0(opt$res0path,'/TraceCount', allchrs[chrind],'.pdf'))
  pct <- ggplot(dataplot, aes(x = reorder(celltype, -tracecount), y = tracecount)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_bw() +
    labs(title = allchrs[chrind], x = "CellType/State/Cluster", y = "Trace count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(pct)
  dev.off()
  cat(allchrs[chrind], "complete the step1 'Distance Calculation'. \n", file = logfile, append = TRUE)
  ## step2 output: DiffScore matrix
  objsmat <- PairCellType(objs)
  objsmat <- as.data.frame(objsmat)
  pair_celltype_p <- furrr::future_map_dfr(1:dim(objsmat)[1], function(n) PairCellTypeP(n, opt, objs, objsmat), .options = furrr::furrr_options(seed = TRUE, globals = c(PairCellTypeP, WilcoxonP, WilcoxonPMat, objsmat, objs, opt, logfile)))
  cat(allchrs[chrind], "complete the step2 'DiffScore Generation'. \n", file = logfile, append = TRUE)
}





time <- Sys.time()
logfile <- paste0("Dory_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
##### parallel
future::plan("multisession")

cat("parallel::detectCores: ", parallel::detectCores(), "\n", file = logfile, append = TRUE)
cat("future::availableCores: ", future::availableCores(), "\n", file = logfile, append = TRUE)
cat("Workers (used cores): ", future::nbrOfWorkers(), "\n", file = logfile, append = TRUE)
cat("Current plan:", capture.output(future::plan()), "\n", file = logfile, append = TRUE)
cores_info <- future::availableCores(methods = "all")
cat(paste(names(cores_info), cores_info, sep = ": ", collapse = "\n"), "\n", file = logfile, append = TRUE)

options(future.globals.maxSize = 1024 * 1024 * 1024 * 3)  # 3 GB

### read command line 
option_list <- list(
  optparse::make_option(c("-i", "--inputCoordFile"), type = "character", help = "Required. Input file containing region coordinates"),
  optparse::make_option(c("-l", "--labelFile"), type = "character", help = "Optional if CoordFile include CellType/State label. Input file containing cell type or state labels for cells"),
  optparse::make_option(c("-n", "--labelcolname"), type = "character", help = "Optional. Please indicate the column name of cell-type/state label. The default is 'Cell_Type'."),
  optparse::make_option(c("-o", "--outputPath"), type = "character", help = "Optional. Path to output files"),
  optparse::make_option(c("-m", "--fileformat"), type = "character", help = "Required. Format of input data ['4DN' or 'csv']"),
  optparse::make_option(c("-c", "--chrnum"), type = "character", help = "Optional. Options: 'one' or 'more'. 'one' means all regions are located on one chromosome; 'more' means regions span multiple chromosomes, typically across the whole genome. ")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
outpath <- opt$outputPath
if(substr(outpath, nchar(outpath), nchar(outpath)) == "/"){
  outpath <- substr(outpath, 1, nchar(outpath)-1)
}

library(ggplot2)
### read input files
if(opt$fileformat == '4DN'){
  indataall <- read_4DNgroup(opt$inputCoordFile, opt$labelFile, opt$labelcolname)
}else if(opt$fileformat == 'csv'){
  if (!is.null(opt$labelFile)|| opt$labelFile == "") {
    indataall <- read_csvgroup_withCT(opt$inputCoordFile, opt$labelcolname)
  }else {
    indataall <- read_csvgroup(opt$inputCoordFile, opt$labelFile, opt$labelcolname)
  }
}else{
  stop("Please indicate the data format of input files: '-m 4DN' or '-m csv' ")
}

chrset <- unique(sort(indataall$Chrom))
if(is.null(opt$chrnum)|| opt$chrnum == ""){
  if(length(chrset) == 1){
    opt$chrnum <- 'one'
  }else if(length(chrset) > 1){
    opt$chrnum <- 'more'
  }else{
    stop("No chromosome detected. Please check the 'Chrom' column in the coordinate file.")
  }
}

if(opt$chrnum == 'one'){
  ### for regions in one chromosome
  opt$res0path <- paste0(outpath, "/S0_DataInfo")
  opt$res1path <- paste0(outpath, "/S1_Distance")
  opt$res2path <- paste0(outpath, "/S2_DiffScore")
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
  indata0 <- indataall
  # region
  if (!all(c("Region_ID") %in% colnames(indata0))) {
    region_out <- call_region(indata0)
    indata <- dplyr::left_join(indata0, region_out, by=c("Chrom","Chrom_Start", "Chrom_End"))
  }else{
    region_out <- data.frame(Region_ID = sort(unique(indata0$Region_ID)))
    indata <- indata0
  }
  region_num <- dim(region_out)[1]
  write.table(region_out, file=paste0(opt$res0path, '/Region.tsv'), sep="\t", quote=FALSE, row.names = FALSE)
  # trace
  celltype <- unique(sort(indata$Cell_Type))
  objs <- paste0("CT_", celltype)
  tracesallct <- c()
  for(i in 1:length(objs)){
    assign(objs[i],  calculate_euclidean_distance(indata[which(indata$Cell_Type == celltype[i]), ], region_num), envir = .GlobalEnv)
    write.table(get(objs[i]), file=paste0(opt$res1path, '/RegionPairsByTrace_', objs[i],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
    tracesallct <- rbind(tracesallct, length(sort(unique(indata$Trace_ID[which(indata$Cell_Type == celltype[i])]))))
  }
  assign("objs", objs, envir = .GlobalEnv)
  for (obj in objs) {
    if (!exists(obj, envir = .GlobalEnv)) {
      stop(paste("Object", obj, "not found in global environment!"))
    }
  }
  dataplot <-data.frame(celltype = celltype, tracecount = tracesallct)
  pdf(paste0(opt$res0path,'/TraceCount.pdf'))
  pct <- ggplot(dataplot, aes(x = reorder(celltype, -tracecount), y = tracecount)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_bw() +
    labs(x = "CellType/State/Cluster", y = "Trace count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(pct)
  dev.off()
  cat("complete the step1 'Distance Calculation'. \n", file = logfile, append = TRUE)
  ## step2 output: DiffScore matrix
  objsmat <- PairCellType(objs)
  pair_celltype_p <- furrr::future_walk(1:dim(objsmat)[1], function(n) PairCellTypeP(n, opt, objs, objsmat), .options = furrr::furrr_options(seed = TRUE, globals = c(PairCellTypeP, WilcoxonP, WilcoxonPMat, objsmat, objs, opt, logfile)))
  cat("complete the step2 'DiffScore Generation'. \n", file = logfile, append = TRUE)

}else if(opt$chrnum == 'more'){
  allchrs <- paste0('chr_', chrset)
  for(chrind in 1:length(chrset)){
    assign(allchrs[chrind], indataall[which(indataall$Chrom == chrset[chrind]), ] , envir = .GlobalEnv)
  }
  allchrsout <- furrr::future_map_dfr(1:length(chrset), function(n) process_EachChr(n, allchrs, opt, outpath), .options = furrr::furrr_options(seed = TRUE, globals = c(process_EachChr, PairCellType, PairCellTypeP, call_region, calculate_euclidean_distance, WilcoxonP, WilcoxonPMat, objsmat, opt, allchrs, logfile, outpath)))

}else{
    stop("Please indicate the chromosome number: '-c one' or '-c more' ")
}

cat("time: ", as.numeric(Sys.time() - time), attr(Sys.time() - time, "units"), "\n", file = logfile, append = TRUE)
