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
    pgreat <- wilcox.test(as.numeric(dismat1[k,]), as.numeric(dismat2[k,]), paired=FALSE, alternative=c("greater"), exact=FALSE)$p.value
    pless <- wilcox.test(as.numeric(dismat1[k,]), as.numeric(dismat2[k,]), paired=FALSE, alternative=c("less"), exact=FALSE)$p.value
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


PairCellTypeP <- function(n, opt, dismat1, dismat2, objsmat, chrname){
  ct <- as.matrix(objsmat[n,])
  pr_num <- dim(dismat1)[1]
  wpboth <- furrr::future_map_dfr(1:pr_num, function(k) WilcoxonP(k, dismat1, dismat2))
  wpmat <- WilcoxonPMat(wpboth)
  #write.table(wpmat, file=paste0(opt$res2path, '/DiffScoreMatrix_', chrname, "_", ct[1], "VS", ct[2],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  inds <- which(lower.tri(wpmat), arr.ind = TRUE)
  lower_tri_df <- data.frame(regionID1 = inds[, 1], regionID2 = inds[, 2],DiffScore = wpmat[inds])
  DiffScore <- lower_tri_df[order(lower_tri_df$DiffScore, decreasing = FALSE), ]
  write.table(DiffScore, file=paste0(opt$res2path, '/DiffScore_', chrname, "_", ct[1], "VS", ct[2],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
  DiffScore_less <- data.frame(lower_tri_df[order(lower_tri_df$DiffScore, decreasing = FALSE), ], rank = seq(1, dim(lower_tri_df)[1], 1))
  DiffScore_greater <- data.frame(lower_tri_df[order(lower_tri_df$DiffScore, decreasing = TRUE), ], rank = seq(1, dim(lower_tri_df)[1],  1))
  thrd=0.05
  DiffScore_l_thrd <- DiffScore_less[which(DiffScore_less$DiffScore < -abs(-log10(thrd))), ]
  DiffScore_g_thrd <- DiffScore_greater[which(DiffScore_greater$DiffScore > abs(-log10(thrd))), ]
  write.table(DiffScore_l_thrd, file=paste0(opt$res2path, '/DRP_less_', chrname, "_", ct[1], "VS", ct[2],'.tsv'), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  write.table(DiffScore_g_thrd, file=paste0(opt$res2path, '/DRP_greater_', chrname, "_", ct[1], "VS", ct[2],'.tsv'), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  # plot
  colnames(wpmat) <- paste0("region", seq(1:dim(wpmat)[1]))
  rownames(wpmat) <- paste0("region", seq(1:dim(wpmat)[1]))
  data <- reshape2::melt(dplyr::mutate(as.data.frame(wpmat), index=row.names(wpmat)), id="index")
  colnames(data) <- c("T1", "T2", "value")
  data$T1 <- factor(data$T1, levels=paste0("region", seq(1,dim(wpmat)[1],1)))
  data$T2 <- factor(data$T2, levels=paste0("region", seq(dim(wpmat)[1],1,-1)))
  # pdf
  pdf(paste0(opt$res2path,'/DiffScoreHeatmap_', chrname, "_", ct[1], "VS", ct[2],'.pdf'))
  p <- ggplot(data,aes(x=T1,y=T2,fill=value))+ 
    scale_fill_gradient2(low="#d60b0e", mid="white", high="#0f70bf", midpoint = 0, na.value = "grey", limits=c(-max(abs(data$value), na.rm = TRUE), max(abs(data$value), na.rm = TRUE)))+
    geom_raster()+
    labs(fill="DiffScore", title = paste0(chrname,": ", ct[1], " VS ", ct[2]),
        x="Region ID", y="Region ID")+
    scale_x_discrete(breaks = paste0("region", seq(5, pr_num, by = 5)), labels=seq(5, pr_num, by=5))+
    scale_y_discrete(breaks = paste0("region", seq(5, pr_num, by = 5)), labels=seq(5, pr_num, by=5))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
  print(p)
  dev.off()
  # png
  png(paste0(opt$res2path,'/DiffScoreHeatmap_', chrname, "_", ct[1], "VS", ct[2],'.png'), width = 6, height = 5, units = "in", res=300)
  p <- ggplot(data,aes(x=T1,y=T2,fill=value))+ 
    scale_fill_gradient2(low="#d60b0e", mid="white", high="#0f70bf", midpoint = 0, na.value = "grey", limits=c(-max(abs(data$value), na.rm = TRUE), max(abs(data$value), na.rm = TRUE)))+
    geom_raster()+
    labs(fill="DiffScore", title = paste0(chrname,": ", ct[1], " VS ", ct[2]),
        x="Region ID", y="Region ID")+
    scale_x_discrete(breaks = paste0("region", seq(5, pr_num, by = 5)), labels=seq(5, pr_num, by=5))+
    scale_y_discrete(breaks = paste0("region", seq(5, pr_num, by = 5)), labels=seq(5, pr_num, by=5))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
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
      if (!all(c("Chrom", "Chrom_Start", "Chrom_End") %in% colnames(indata0))) {
        stop("Input must contain either 'Region_ID' or all of 'Chrom', 'Chrom_Start', and 'Chrom_End'.")
      }
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
    write.table(region_out, file=paste0(opt$res0path, '/RegionIn_', chrset[chrind], '.tsv'), sep="\t", quote=FALSE, row.names = FALSE)
  }else{
    stop("The genomic regions in these two celltypes/states/clusters are not the same")
  }
  dataplot <-data.frame(celltype = objs, tracecount = tracesallct)
  pdf(paste0(opt$res0path,'/TraceCount', chrset[chrind], '.pdf'))
  pct <- ggplot(dataplot, aes(x = celltype, y = tracecount)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_bw() +
    labs(title = chrset[chrind], x = "", y = "Trace count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(pct)
  dev.off()
  ## step1 output: region paris' distance
  dismat1 <- calculate_euclidean_distance(foredata, region_num)
  dismat2 <- calculate_euclidean_distance(backdata, region_num) 
  write.table(dismat1, file=paste0(opt$res1path, '/RegionPairsByTrace_', chrset[chrind],'_', objs[1],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  write.table(dismat2, file=paste0(opt$res1path, '/RegionPairsByTrace_', chrset[chrind], '_', objs[2],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  cat(chrset[chrind], "complete the step1 'Distance Calculation'. \n", file = logfile, append = TRUE)
  ## step2 output: DiffScore matrix
  objsmat <- matrix(c(objs[1], objs[2]), nrow=1, ncol=2)
  pair_celltype_p <- PairCellTypeP(1, opt, dismat1, dismat2, objsmat, chrset[chrind])
  cat(chrset[chrind], "complete the step2 'DiffScore Generation'. \n", file = logfile, append = TRUE)
}


########## bootstrapping functions

PairCellTypePBoot <- function(n, opt, dismat1, dismat2, objsmat, chrname, bt_outpath, bboot){
  ct <- as.matrix(objsmat[n,])
  pr_num <- dim(dismat1)[1]
  wpboth <- furrr::future_map_dfr(1:pr_num, function(k) WilcoxonP(k, dismat1, dismat2))
  wpmat <- WilcoxonPMat(wpboth)
  inds <- which(lower.tri(wpmat), arr.ind = TRUE)
  lower_tri_df <- data.frame(regionID1 = inds[, 1], regionID2 = inds[, 2],DiffScore = wpmat[inds])
  DiffScore <- lower_tri_df[order(lower_tri_df$DiffScore, decreasing = FALSE), ]
  write.table(DiffScore, file=paste0(bt_outpath, '/DiffScore_', chrname, "_BT", bboot, "_", ct[1], "VS", ct[2],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
  cat("complete the ", n, "th pairs of cell types. \n", file = logfile, append = TRUE)
}


process_EachChrBoot <- function(chrind,chrset, foredataall, backdataall, opt, logfile, objs, bt_outpath, bboot){
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
      if (!all(c("Chrom", "Chrom_Start", "Chrom_End") %in% colnames(indata0))) {
        stop("Input must contain either 'Region_ID' or all of 'Chrom', 'Chrom_Start', and 'Chrom_End'.")
      }
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
    #write.table(region_out, file=paste0(opt$res0path, '/RegionIn_', chrset[chrind], '.tsv'), sep="\t", quote=FALSE, row.names = FALSE)
  }else{
    stop("The genomic regions in these two celltypes/states/clusters are not the same")
  }
  ## step1 output: region paris' distance
  dismat1 <- calculate_euclidean_distance(foredata, region_num)
  dismat2 <- calculate_euclidean_distance(backdata, region_num) 
  ## step2 output: DiffScore matrix
  objsmat <- matrix(c(objs[1], objs[2]), nrow=1, ncol=2)
  pair_celltype_p <- PairCellTypePBoot(1, opt, dismat1, dismat2, objsmat, chrset[chrind], bt_outpath, bboot)
  cat(chrset[chrind], "complete the step2 'DiffScore Generation'. \n", file = logfile, append = TRUE)
}


Run_Dory_Once <- function(bboot, foredataall, backdataall, opt, bt_outpath, logfile, objs, inner_parallel = TRUE){
  chrset <- unique(sort(c(foredataall$Chrom, backdataall$Chrom)))
  if(is.null(opt$chrnum)|| opt$chrnum == ""){
    if(length(chrset) == 1){
      opt$chrnum <- 'one'
    }else if(length(chrset) > 1){
      opt$chrnum <- 'more'
    }else{
      cat("Warning: No chromosome information was detected. The input will be processed as one chromosome by default, and ChrNum will be set to '-c one'. 
      This may occur if the chromosome column is missing or incorrectly named (for example, 'chrom' instead of 'Chrom'), or if the input uses Region_ID only.\n")
      opt$chrnum <- "one"
    }
  }

  if(opt$chrnum == 'one'){
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
        if (!all(c("Chrom", "Chrom_Start", "Chrom_End") %in% colnames(indata0))) {
          stop("Input must contain either 'Region_ID' or all of 'Chrom', 'Chrom_Start', and 'Chrom_End'.")
        }
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
      #write.table(region_out, file=paste0(opt$res0path, '/Region.tsv'), sep="\t", quote=FALSE, row.names = FALSE)
    }else{
      stop("The genomic regions in these two celltypes/states/clusters are not the same")
    }
    ## step1 output: region paris' distance
    dismat1 <- calculate_euclidean_distance(foredata, region_num)
    dismat2 <- calculate_euclidean_distance(backdata, region_num)
    cat("time: ", as.numeric(Sys.time() - time), attr(Sys.time() - time, "units"), "\n", file = logfile, append = TRUE)
    ## step2 output: DiffScore matrix
    objsmat <- matrix(c(objs[1], objs[2]), nrow=1, ncol=2)
    chrname <- 'onechr'
    pair_celltype_p <- PairCellTypePBoot(1, opt, dismat1, dismat2,  objsmat, chrname, bt_outpath, bboot)
    cat("complete the step2 'DiffScore Generation'. \n", file = logfile, append = TRUE)
  }else if(opt$chrnum == 'more'){
    ### for regions spanning more chromosomes, typically across the whole genome
    allchrsout <- furrr::future_map_dfr(1:length(chrset), function(n) process_EachChrBoot(n, chrset, foredataall, backdataall, opt, logfile, objs, bt_outpath, bboot), .options = furrr::furrr_options(seed = TRUE))
  }else{
      stop("Please indicate the chromosome number: '-c one' or '-c more' ")
  }

}



prep_trace_bootstrap <- function(df, trace_col = "Trace_ID") {
  trace_vals <- df[[trace_col]]
  idx_by_trace <- split(seq_len(nrow(df)), trace_vals)
  list(df = df, trace_col = trace_col, idx_by_trace = idx_by_trace, trace_ids = names(idx_by_trace))
}

bootstrap_from_prep <- function(prep) {
  sampled_ids <- sample(prep$trace_ids, length(prep$trace_ids), replace = TRUE)
  idx_list <- prep$idx_by_trace[sampled_ids]
  out_idx <- unlist(idx_list, use.names = FALSE)
  res <- prep$df[out_idx, , drop = FALSE]
  res[[prep$trace_col]] <- rep(seq_along(idx_list), lengths(idx_list))
  rownames(res) <- NULL
  res
}


BootstrapCI <- function(n, opt, objs, logfile, bt_outpath, outpath, chrname){
  objsmat <- matrix(c(objs[1], objs[2]), nrow=1, ncol=2)
  #chrname <- 'onechr'
  ct <- as.matrix(objsmat[n,])
  res_list <- lapply(seq_len(as.numeric(opt$nboot)), function(id) {
    data <- read.table(paste0(bt_outpath, '/DiffScore_', chrname, "_BT", id, "_", ct[1], "VS", ct[2],'.tsv'), header = TRUE, sep = "\t")
    data <- data[, c("regionID1", "regionID2", "DiffScore")]
    colnames(data)[3] <- paste0("DiffScore", id)
    data
  })
  res_all <- Reduce(function(x, y) merge(x, y, by = c("regionID1", "regionID2"), all = TRUE), res_list)
  score_mat <- as.matrix(res_all[, c(-1, -2)])
  ci_level <- as.numeric(opt$CIpct)
  lower_p <- (1 - ci_level) / 2
  upper_p <- (1 + ci_level) / 2

  ci_mat <- t(apply(score_mat, 1, function(x) {
    quantile(x, probs = c(lower_p, upper_p), na.rm = TRUE)
  }))
  colnames(ci_mat) <- c("CI_lower", "CI_upper")
  rgci <- cbind(res_all[, c(1, 2)], ci_mat)

  bootpath <- paste0(outpath, "/S3_BootStrap")
  if(!dir.exists(bootpath)){
    dir.create(bootpath, recursive = TRUE)
  }else{
    cat("Step3 directory already exists! \n", file = stdout())
  }
  outrgci <- rgci[order(as.numeric(rgci$regionID2), as.numeric(rgci$regionID1), decreasing = FALSE), ]
  write.table(outrgci, file=paste0(bootpath, '/DiffScore', ci_level,'CI_', chrname, "_", ct[1], "VS", ct[2],'.tsv'), sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
  #return(outrgci)
}

#args_full <- commandArgs(trailingOnly = FALSE)
#file_arg <- grep("^--file=", args_full, value = TRUE)
#if (length(file_arg) == 0) {
#  stop("Cannot determine script path. Please run with Rscript.")
#}
#script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
#script_dir  <- dirname(script_path)
#source(file.path(script_dir, "bootstrap_single_core.R"))

time <- Sys.time()
future::plan("multisession")
options(future.globals.maxSize = 1024 * 1024 * 1024)  # 1 GB
### read command line 
option_list <- list(
  optparse::make_option(c("-a", "--inputPath1"), type = "character", help = "Required. The foreground input file"),
  optparse::make_option(c("-b", "--inputPath2"), type = "character", help = "Required. The background input file"),
  optparse::make_option(c("-o", "--outputPath"), type = "character", help = "Optional. Path to output directory"),
  optparse::make_option(c("-m", "--fileformat"), type = "character", help = "Required. Format of input data ['4DN'] or ['csv']"),
  optparse::make_option(c("-c", "--chrnum"), type = "character", help = "Optional. Options: 'one' or 'more'. 'one' means all regions are located on one chromosome; 'more' means regions span multiple chromosomes, typically across the whole genome."),
  optparse::make_option(c("--nboot"), type = "numeric", help = "Optional. Resampling times for bootstrapping. Default is 1000"),
  optparse::make_option(c("--CIpct"), type = "numeric", help = "Optional. Confidence interval for bootstrapping [numeric number from 0 to 1]. Default of 0.9"),
  optparse::make_option(c("--seed"), type = "numeric", help = "Optional. Seed for bootstrapping. Default is 1")
 
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

listsdata <- data.frame(strsplit(opt$inputPath1, "[./]"), strsplit(opt$inputPath2, "[./]")) 
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
    cat("Warning: No chromosome information was detected. The input will be processed as one chromosome by default, and ChrNum will be set to '-c one'. 
    This may occur if the chromosome column is missing or incorrectly named (for example, 'chrom' instead of 'Chrom'), or if the input uses Region_ID only.\n")
    opt$chrnum <- "one"
  }
}

#### step0-2

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
      if (!all(c("Chrom", "Chrom_Start", "Chrom_End") %in% colnames(indata0))) {
        stop("Input must contain either 'Region_ID' or all of 'Chrom', 'Chrom_Start', and 'Chrom_End'.")
      }
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
  chrname <- 'onechr'
  pair_celltype_p <- PairCellTypeP(1, opt, dismat1, dismat2,  objsmat, chrname)
  cat("complete the step2 'DiffScore Generation'. \n", file = logfile, append = TRUE)
}else if(opt$chrnum == 'more'){
  ### for regions spanning more chromosomes, typically across the whole genome
  allchrsout <- furrr::future_map_dfr(1:length(chrset), function(n) process_EachChr(n, chrset, foredataall, backdataall, opt, logfile, outpath, objs), .options = furrr::furrr_options(seed = TRUE))
}else{
    stop("Please indicate the chromosome number: '-c one' or '-c more' ")
}


############# Step3 Bootstrapping 
set.seed(as.numeric(opt$seed))
fore_prep <- prep_trace_bootstrap(foredataall, "Trace_ID")
back_prep <- prep_trace_bootstrap(backdataall, "Trace_ID")

bt_outpath <- file.path(outpath, "BootDiffScore")
if(!dir.exists(bt_outpath)){
  dir.create(bt_outpath, recursive = TRUE)
}else{
  cat("BootDiffScore directory already exists! \n", file = stdout())
}
for (b in seq_len(as.numeric(opt$nboot))) {
  fore_bt <- bootstrap_from_prep(fore_prep)
  back_bt <- bootstrap_from_prep(back_prep)
  Run_Dory_Once(bboot = b, foredataall = fore_bt, backdataall = back_bt, opt = opt, bt_outpath = bt_outpath, objs = objs, logfile = logfile,inner_parallel = TRUE)
}

#### bootstrap confidence interval
if(opt$chrnum == 'one'){
  chrname <- 'onechr'
  BootstrapCI(1, opt, objs, logfile, bt_outpath, outpath, chrname)
  #cat("complete the step2 'DiffScore Generation'. \n", file = logfile, append = TRUE)
}else if(opt$chrnum == 'more'){
  ### for regions spanning more chromosomes, typically across the whole genome
  allchrsout <- furrr::future_map_dfr(1:length(chrset), function(n) BootstrapCI(1, opt, objs, logfile, bt_outpath, outpath, chrset[n]), .options = furrr::furrr_options(seed = TRUE))
}else{
    stop("Please indicate the chromosome number: '-c one' or '-c more' ")
}



unlink(bt_outpath, recursive = TRUE, force = TRUE)

cat("time: ", as.numeric(Sys.time() - time), attr(Sys.time() - time, "units"), "\n", file = logfile, append = TRUE)



