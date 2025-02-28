read_4DNgroup <- function(coordfile, cellfile){
  # read files downloaded from 4DN datasets
  lines <- readLines(coordfile)
  filtered_lines <- lines[!grepl("^\\s*[^A-Za-z0-9]", lines)]
  temp_file <- tempfile()
  writeLines(filtered_lines, temp_file)
  data <- read.csv(temp_file, header = FALSE)
  colnames(data) <- c("Spot_ID", "Trace_ID", "X", "Y", "Z", "Chrom", "Chrom_Start", "Chrom_End", "Cell_ID")
  outdata <- data.frame(Trace_ID = as.numeric(data$Trace_ID),
                        X = as.numeric(data$X), Y = as.numeric((data$Y)), Z = as.numeric(data$Z),
                        Chrom = as.character(data$Chrom), Chrom_Start = as.numeric(data$Chrom_Start), Chrom_End = as.numeric(data$Chrom_End),
                        Cell_ID = as.numeric(data$Cell_ID))
  cell_label <- read.csv(cellfile, header = TRUE)
  colnames(cell_label) <- c("Cell_ID", "Cell_Type")
  outcelltype <- data.frame(Cell_ID = as.numeric(cell_label$Cell_ID), Cell_Type = as.character(cell_label$Cell_Type))
  outdata_celltype <- dplyr::left_join(outdata, cell_label, by=c("Cell_ID"))
  file.remove(temp_file)
  return(outdata_celltype)
}


read_csvgroup <- function(coordfile, cellfile){
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
  colnames(cell_label) <- c("Cell_ID", "Cell_Type")
  outcelltype <- data.frame(Cell_ID = as.numeric(cell_label$Cell_ID), Cell_Type = as.character(cell_label$Cell_Type))
  outdata_celltype <- dplyr::left_join(outdata, cell_label, by=c("Cell_ID"))
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


calculate_euclidean_distance <- function(dataset){
  dataset <- dataset[which(dataset$X != 0 & dataset$Y != 0 & dataset$Z != 0), ]
  if (!all(c("Region_ID") %in% colnames(dataset))) {
    region_out <- call_region(dataset)
    coords <- dplyr::left_join(dataset, region_out, by=c("Chrom","Chrom_Start", "Chrom_End"))
  }else{
    region_out <- data.frame(Region_ID = sort(unique(dataset$Region_ID)))
    coords <- dataset
  }
  region_num <- dim(region_out)[1]
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
      #print(cbind(i,j,m))
    }
    i <- i + 1
  }
  return(mat)
} 

PairCellTypeP <- function(n){
  ct <- objsmat[n,]
  pr_num <- dim(get(ct[1]))[1]
  wpboth <- furrr::future_map_dfr(1:pr_num, function(k) WilcoxonP(k, ct))
  #write.table(wpboth, file=paste0(opt$res2path, '/WpCauchyPGreatLess_', ct[1], "VS", ct[2],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
  wpmat <- WilcoxonPMat(wpboth)
  write.table(wpmat, file=paste0(opt$res2path, '/DiffScoreMatrix_', ct[1], "VS", ct[2],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  # plot
  colnames(wpmat) <- paste0("region", seq(1:dim(wpmat)[1]))
  rownames(wpmat) <- paste0("region", seq(1:dim(wpmat)[1]))
  data <- reshape2::melt(dplyr::mutate(as.data.frame(wpmat), index=row.names(wpmat)), id="index")
  colnames(data) <- c("T1", "T2", "value")
  data$T1 <- factor(data$T1, levels=paste0("region", seq(1,dim(wpmat)[1],1)))
  data$T2 <- factor(data$T2, levels=paste0("region", seq(dim(wpmat)[1],1,-1)))
  pdf(paste0(opt$res2path,'/DiffScoreHeatmap_', ct[1], "VS", ct[2],'.pdf'), height = 2.6, width = 3.5)
  p <- ggplot(data,aes(x=T1,y=T2,fill=value))+ 
    scale_fill_gradient2(low="#d60b0e", mid="white", high="#0f70bf", midpoint = 0, na.value = "grey", limits=c(-max(abs(data$value), na.rm = TRUE), max(abs(data$value), na.rm = TRUE)))+
    geom_raster()+
    labs(fill="DiffScore",
        x="Region ID",
        y="Region ID")+
    #theme(axis.ticks = element_blank())+
    scale_x_discrete(breaks = paste0("region", seq(5, dim(wpmat)[1], by = 5)), labels=seq(5, dim(wpmat)[1], by=5))+
    scale_y_discrete(breaks = paste0("region", seq(5, dim(wpmat)[1], by = 5)), labels=seq(5, dim(wpmat)[1], by=5))
  print(p)
  dev.off()  
  cat("complete the ", n, "th pairs of cell types. \n")
  cat("time: ", as.numeric(Sys.time() - time), attr(Sys.time() - time, "units"),  "\n")
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


time <- Sys.time()
future::plan("multicore", workers = min(8, max(1, (parallel::detectCores() - 1))))
### read command line 
option_list <- list(
  optparse::make_option(c("-i", "--inputCoordFile"), type = "character", help = "Required. Input file containing region coordinates"),
  optparse::make_option(c("-l", "--labelFile"), type = "character", help = "Required. Input file containing cell type or state labels for cells"),
  optparse::make_option(c("-o", "--outputPath"), type = "character", help = "Optional. Path to output files"),
  optparse::make_option(c("-m", "--format"), type = "character", help = "Required. Format of input data ['4DN' or 'csv']")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
### creat output directory for step1: region paris' distance; and step2: DiffScore matrix.
outpath <- opt$outputPath
if(substr(outpath, nchar(outpath), nchar(outpath)) == "/"){
  outpath <- substr(outpath, 1, nchar(outpath)-1)
}
opt$res1path <- paste0(outpath, "/S1_Distance")
opt$res2path <- paste0(outpath, "/S2_DiffScore")
if(!dir.exists(paste0(outpath, "/S1_Distance"))){
  dir.create(paste0(outpath, "/S1_Distance"))
}else{
  cat("Step1 directory already exists! \n")
}
if(!dir.exists(paste0(outpath, "/S2_DiffScore"))){
  dir.create(paste0(outpath, "/S2_DiffScore"))
}else{
  cat("Step2 directory already exists! \n")
}


library(ggplot2)

if(opt$format == '4DN'){
  indata <- read_4DNgroup(opt$inputCoordFile, opt$labelFile)
  celltype <- unique(sort(indata$Cell_Type))
  objs <- paste0("CT_", celltype)
  for(i in 1:length(objs)){
    assign(objs[i],  calculate_euclidean_distance(indata[which(indata$Cell_Type == celltype[i]), ]))
    write.table(get(objs[i]), file=paste0(opt$res1path, '/RegionPairsByTrace_', objs[i],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  }
}else if(opt$format == 'csv'){
  indata <- read_csvgroup(opt$inputCoordFile, opt$labelFile)
  celltype <- unique(sort(indata$Cell_Type))
  objs <- paste0("CT_", celltype)
  for(i in 1:length(objs)){
    assign(objs[i],  calculate_euclidean_distance(indata[which(indata$Cell_Type == celltype[i]), ]))
    write.table(get(objs[i]), file=paste0(opt$res1path, '/RegionPairsByTrace_', objs[i],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  }
}else{
  stop("Please indicate the data format of input files: '-m 4DN' or '-m csv' ")
}
cat("complete the step1 'Distance Calculation'. \n")
cat("time: ", as.numeric(Sys.time() - time), attr(Sys.time() - time, "units"), "\n")
## step2 output: DiffScore matrix
objsmat <- PairCellType(objs)
pair_celltype_p <- furrr::future_map_dfr(1:dim(objsmat)[1], function(n) PairCellTypeP(n), .options = furrr::furrr_options(seed = TRUE))
cat("complete the step2 'DiffScore Generation'. \n")
cat("time: ", as.numeric(Sys.time() - time), attr(Sys.time() - time, "units"), "\n")

