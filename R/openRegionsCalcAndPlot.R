#' The function calcOpenSignal takes the input BED file, overlaps it with all defined
#' open chromatin regions across cell types and returns a matrix, where each row is 
#' the input genomic region (if overlap was found), each column is a cell type, and the 
#' value is a normalized ATAC-seq signal.
#
#' @param  bedInput Genomic regions to be analyzed in form of data.table, data.frame, or 
#'            GRanges object.
#' @param cellMatrix Matrix with open chromatin signal values, rows are genomic regions, 
#'            columns are cell types. First column contains information about the genomic 
#'            region in following form: chr_start_end. Can be either data.frame or data.table.
#' @export
#' @examples
#' \dontrun{
#' signalMatrix = calcOpenSignal(bedInput, cellMatrix)
#' }

calcOpenSignal = function(bedInput,
                          cellMatrix){
  
  # if the cellMatrix is in data.frame format, convert it to data.table
  if(class(cellMatrix)[1] == "data.frame"){
    cellMatrix = as.data.table(cellMatrix)
  }
  
  # get the genomic coordinates from the open chromatin signal matrix - the first column
  # conver the chr_start_end into a three column data.table
  openRegions = cellMatrix[,1]
  colnames(openRegions) = "V1"
  openRegions[, c("chr", "start", "end") := tstrsplit(V1, "_", fixed=TRUE)]
  numericColumns = c("start", "end")
  openRegions[, (numericColumns ) := lapply(.SD, as.numeric), .SDcols = numericColumns]
  openRegions = openRegions[,.(chr, start, end)]
  
  # if the class of query BED file differs from data.table, convert it to data.table
  if (class(bedInput)[1] == "data.table") {
    bedInput = bedInput
  } else if (class(bedInput) == "data.frame") {
    bedInput = as.data.table(bedInput)
  } else if (class(bedInput) == "GRanges"){
    bedInput = as.data.table(data.frame(chr = seqnames(bedInput),
                                        start = start(bedInput),
                                        end = end(bedInput)))
  } else {
    stop("Genomic coordinates must be passed as data.frame or data.table or GRanges object.")
  }
  
  # select just the fist three columns and give them name chr, start, end
  # create a 4th column with peak name in following format: chr_start_end
  bedInput = bedInput[, 1:3]
  colnames(bedInput) = c("chr", "start", "end")
  bedInput[, peakName := paste(bedInput[,chr], bedInput[,start], bedInput[,end], sep = "_")]
  
  # find which regions of the genomic regions of interest overlap with the open
  # chromatin signal matrix regions
  setkey(bedInput, chr, start, end)
  overlaps = foverlaps(openRegions,bedInput, which = T)
  overlaps = na.omit(overlaps)
  
  # extract the regions which overlap with query, assign the query peaks to them 
  # and calculate the sum of the signal within these regions
  signalMatrix = cellMatrix[overlaps[,xid]] 
  signalMatrix[,V1:=NULL]
  signalMatrix[, queryPeak := bedInput[overlaps[,yid], peakName]]
  signalMatrix = signalMatrix[, lapply(.SD, sum), by = .(queryPeak)]
  
  return(signalMatrix)
}

#' The function plotOpenSignal visualizes the signalMatrix obtained from calcOpenSignal.
#'
#' @param signalMatrix - output data.frame from calcOpenSignal function
#' @param  plotType - what plot type should be used to visualize the results, options are:
#'           jitter (default) - jitter plot with box plot on top / boxPlot - box plot
#'           without individual points and outliers / barPlot - bar height represents the
#'           median signal value for a given cell type
#' @param cellGroup - this option allows to selcet a group of cells to be plotted, if NA (default)
#'           all available cell groups are ploted, available options: {"blood", "bone", "CNS", 
#'           "embryonic", "eye", "foreskin", "gastrointestinal", "heart", "liver", "lymphatic", 
#'           "mammaryGland", "mouth", "respiratorySystem", "skeletalMuscle", "skin", 
#'            "urinarySystem", "vasculature"}, can be passed as a singe character string or vector
#'            of strings
#' @export
#' @examples
#' \dontrun{
#' signalMatrix = calcOpenSignal(bedInput)
#' plotOpenSignal(signalMatrix)
#' plotOpenSignal(signalMatrix, plotType = "barPlot", cellGroup = "blood")
#' }
plotOpenSignal = function(signalMatrix, 
                          plotType = "jitter", 
                          cellGroup = NA){
  
  # hg 19 has 67 cell typer, hg38 has 74 cell types - 
  # based on this fact distinguish between geneome and upload metadata
  # for color scheme
  if(ncol(signalMatrix) == 68){
    # upload metadata for coloring
    metadata = read.csv("data/cellTissue_metadata_hg19.csv", as.is = T)
  } else if (ncol(signalMatrix) == 75){
    metadata = read.csv("data/cellTissue_metadata_hg38.csv", as.is = T)
  } else {
    stop("Signal matrix dimensions do not match the number of available cell types for either hg19 or hg38.")
  }
  
  
  # reshape the signal matrix into ggplot usable form 
  # attach the metadata for coloring
  # sort table alphabetically by tissue-cellType
  plotSignalMatrix = signalMatrix %>% 
    gather(key = "cellType", value = "signal", -one_of("queryPeak")) %>% 
    left_join(metadata, by = "cellType") %>% 
    arrange(tissue, cellType) %>% 
    unite("mixedVar", tissue, cellType, sep = "_", remove = F)
  
  # if user defines cell group, filter the data
  if (length(cellGroup) == 1){
    if (is.na(cellGroup)){
      plotSignalMatrix = plotSignalMatrix
    } else if (cellGroup %in% levels(factor(metadata$tissue))){
      plotSignalMatrix = plotSignalMatrix %>% 
        filter(tissue == cellGroup)
    } else {
      stop("The input cell group is not in predefined list of options.")
    }
  } else if (all(cellGroup %in% levels(factor(metadata$tissue)))) {
    plotSignalMatrix = plotSignalMatrix %>% 
      filter(tissue %in% cellGroup)
  } else {
    stop("At least one of the input cell groups is not in predefined list of options.")
  }
  
  # arrange labels in a way, that corresponding groups are plotted together
  myLabels = plotSignalMatrix  %>% 
    ungroup() %>% 
    select(mixedVar, cellType) %>% 
    distinct() %>% 
    arrange(mixedVar) %>% 
    mutate(spaceLabel = str_replace_all(cellType, "_", " ")) %>% 
    mutate(spaceLabel = str_replace(spaceLabel, "xx", ",")) %>% 
    mutate(spaceLabel = str_replace(spaceLabel, " ,", ","))
  
  # set up a color palette for plotting
  myColors = c("#E31A1C","#666666","#B3DE69","#A65628","#33A02C","#E6AB02","#F0027F",
               "#FDC086","#FFFF99","#B3E2CD",
               "#B3CDE3","#66A61E","#F4CAE4","#80B1D3","#FFED6F","#B15928","#999999")
 
   # get box plot staistics to set plot limits - outliers throw off the limits a lot
  boxStats = plotSignalMatrix %>% 
    group_by(mixedVar) %>% 
    summarise(boxStats = list(boxplot.stats(signal)$stats)) %>% 
    unnest(cols = boxStats) 
  
  # get mediand of signal values to make a bar plot
  barPlotStats = plotSignalMatrix %>% 
    group_by(mixedVar) %>% 
    summarise(medianBar = median(signal)) %>% 
    left_join(plotSignalMatrix %>% 
                select(-c("signal", "queryPeak")) %>% 
                distinct(), by = "mixedVar")
  
  # do the plotting
  p = ggplot(plotSignalMatrix,
             aes(x = mixedVar, y = signal))
  
  if (plotType == "jitter"){
    jitterPlot = p + geom_jitter(alpha = 0.5, height = 0, width = 0.35, aes(color = tissue)) + 
      geom_boxplot(outlier.colour = NA, fill = NA) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size=15)) +
      xlab("") +
      ylab("normalized signal") + 
      scale_x_discrete(labels = myLabels$spaceLabel) +
      scale_fill_manual(values=myColors) + 
      scale_color_manual(values=myColors)
    print(jitterPlot)
  } else if (plotType == "boxPlot") {
    boxPlot = p + geom_boxplot(outlier.colour = NA, aes(fill = tissue), alpha = 0.9) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size=15)) +
      xlab("") +
      ylab("normalized signal") + 
      scale_x_discrete(labels = myLabels$spaceLabel) +
      scale_fill_manual(values=myColors) + 
      scale_color_manual(values=myColors) +
      ylim(min(boxStats$boxStats), max(boxStats$boxStats))
    suppressWarnings(print(boxPlot))
  } else if (plotType == "barPlot") {
    barPlot = ggplot(barPlotStats, aes(x = mixedVar, y = medianBar, fill = tissue)) +
      geom_col(alpha = 0.9)+
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size=15)) +
      xlab("") +
      ylab("median normalized signal") + 
      scale_x_discrete(labels = myLabels$spaceLabel) +
      scale_fill_manual(values=myColors)
    print(barPlot)
  } else {
    stop("Plot type does not match any of the available options. Available options: jitter, boxPlot, barPlot. ")
  }
  
  
}

