#' The function calcOpenSignal takes the input BED file, overlaps it with all defined
#' open chromatin regions across cell types and returns a matrix, where each row is 
#' the input genomic region (if overlap was found), each column is a cell type, and the 
#' value is a normalized ATAC-seq signal.
#
#' @param  query Genomic regions to be analyzed in form of data.table, data.frame, or 
#'            GRanges object.
#' @param cellMatrix Matrix with open chromatin signal values, rows are genomic regions, 
#'            columns are cell types. First column contains information about the genomic 
#'            region in following form: chr_start_end. Can be either data.frame or data.table.
#' @return 
#' A data.table with cell specific open chromatin signal values for query regions.
#' 
#' @export

calcOpenSignal = function(query,
                          cellMatrix){
  
  # if the cellMatrix is in data.frame format, convert it to data.table
  if(class(cellMatrix)[1] == "data.frame"){
    cellMatrix = as.data.table(cellMatrix)
  }
  
  # get the genomic coordinates from the open chromatin signal matrix - the first column
  # convert the chr_start_end into a three column data.table
  openRegions = cellMatrix[,1]
  colnames(openRegions) = "V1"
  openRegions[, c("chr", "start", "end") := tstrsplit(V1, "_", fixed=TRUE)]
  numericColumns = c("start", "end")
  openRegions[, (numericColumns ) := lapply(.SD, as.numeric), .SDcols = numericColumns]
  openRegions = openRegions[,.(chr, start, end)]
  
  # if the class of query BED file differs from data.table, convert it to data.table
  if (class(query)[1] == "data.table") {
    query = query
  } else if (class(query) == "data.frame") {
    query = as.data.table(query)
  } else if (class(query) == "GRanges"){
    query = as.data.table(data.frame(chr = seqnames(query),
                                        start = start(query),
                                        end = end(query)))
  } else {
    stop("Genomic coordinates must be passed as data.frame or data.table or GRanges object.")
  }
  
  # select just the fist three columns and give them name chr, start, end
  # create a 4th column with peak name in following format: chr_start_end
  query = query[, 1:3]
  colnames(query) = c("chr", "start", "end")
  query[, peakName := paste(query[,chr], query[,start], query[,end], sep = "_")]
  
  # find which regions of the genomic regions of interest overlap with the open
  # chromatin signal matrix regions
  setkey(query, chr, start, end)
  overlaps = foverlaps(openRegions,query, which = T)
  overlaps = na.omit(overlaps)
  
  # extract the regions which overlap with query, assign the query peaks to them 
  # and calculate the sum of the signal within these regions
  signalMatrix = cellMatrix[overlaps[,xid]] 
  signalMatrix[,-1]
  signalMatrix[, queryPeak := query[overlaps[,yid], peakName]]
  signalMatrix = signalMatrix[, lapply(.SD, sum), by = .(queryPeak)]
  
  return(signalMatrix)
}

#' The function plotOpenSignal visualizes the signalMatrix obtained from calcOpenSignal.
#'
#' @param signalMatrix Output data.table from \code{calcOpenSignal} function.
#' @param  plotType Options are: jitter (default) - jitter plot with box plot on top / 
#'           boxPlot - box plot without individual points and outliers / 
#'           barPlot - bar height represents the median signal value for a given cell type.
#' @param cellGroup - This option allows to selcet a group of cells to be plotted, if NA (default)
#'           all available cell groups are ploted, available options: {"blood", "bone", "CNS", 
#'           "embryonic", "eye", "foreskin", "gastrointestinal", "heart", "liver", "lymphatic", 
#'           "mammaryGland", "mouth", "respiratorySystem", "skeletalMuscle", "skin", 
#'            "urinarySystem", "vasculature"}, can be passed as a singe character string or vector
#'            of strings.
#' @param cellTypeMetadata Metadata for cell type - tissue association. This option is for users, 
#'          who provide their own open region signal matrix. The cellTypeMetadata matrix mast contain
#'          two columns called cellType and tissue. cellType column containes the cell type names 
#'          in provided signalMatrix column names. The tissue columns provides an information, which 
#'          tissue the cell type comes from.
#' @param colorScheme Provide color values for each tissue if you want to change the default colors.
#' @export
#' @examples
#' \dontrun{
#' signalMatrix = calcOpenSignal(query, cellMatrix)
#' plotOpenSignal(signalMatrix)
#' plotOpenSignal(signalMatrix, plotType = "barPlot", cellGroup = "blood")
#' }
plotOpenSignal = function(signalMatrix, 
                          plotType = "jitter", 
                          cellGroup = NA,
                          cellTypeMetadata = NA, 
                          colorScheme = c("#E31A1C","#666666","#B3DE69",
                                          "#A65628","#33A02C","#E6AB02",
                                          "#F0027F", "#FDC086","#FFFF99",
                                          "#B3E2CD", "#B3CDE3","#66A61E",
                                          "#F4CAE4","#80B1D3","#FFED6F",
                                          "#B15928","#999999")){
  
  if(is.na(cellTypeMetadata)){
    # upload metadata for coloring
    cellTypeMetadata = data("cellTypeMetadata")
  } 
  
  # reshape the signal matrix into ggplot usable form 
  # attach the metadata for coloring
  # sort table alphabetically by tissue-cellType
  plotSignalMatrix = melt(signalMatrix, id.vars = "queryPeak", variable.name = "cellType", value.name = "signal")
  setkey(plotSignalMatrix, cellType)
  setkey(cellTypeMetadata, cellType)
  plotSignalMatrix = merge(plotSignalMatrix, cellTypeMetadata, all = F)
  plotSignalMatrix[, lowerCaseTissue := tolower(tissue)]
  setorder(plotSignalMatrix, lowerCaseTissue, cellType)
  plotSignalMatrix[, mixedVar := paste(plotSignalMatrix[,tissue], plotSignalMatrix[,cellType], sep = "_")]
  
  # if user defines cell group, filter the data
  if (length(cellGroup) == 1){
    if (is.na(cellGroup)){
      plotSignalMatrix = plotSignalMatrix
    } else if (cellGroup %in% levels(factor(cellTypeMetadata$tissue))){
      plotSignalMatrix = plotSignalMatrix[tissue == cellGroup]
    } else {
      stop("The input cell group is not in predefined list of options.")
    }
  } else if (all(cellGroup %in% levels(factor(cellTypeMetadata$tissue)))) {
    plotSignalMatrix = plotSignalMatrix[tissue %in% cellGroup]
  } else {
    stop("At least one of the input cell groups is not in predefined list of options.")
  }
  
  # arrange labels in a way, that corresponding groups are plotted together
  myLabels = unique(plotSignalMatrix[, .(mixedVar, cellType)])
  myLabels[, spaceLabel := gsub("_", " ", myLabels[, cellType])]
  
  # get box plot staistics to set plot limits - outliers throw off the limits a lot
  boxStats = plotSignalMatrix[, .(boxStats = list(boxplot.stats(signal)$stats)), by = mixedVar]
  minBoxLimit = min(unlist(boxStats$boxStats))
  maxBoxLimit = max(unlist(boxStats$boxStats))
  
  # get mediand of signal values to make a bar plot
  barPlotStats = plotSignalMatrix[, .(medianBar = median(signal)), by = mixedVar]
  tableToMerge = unique(plotSignalMatrix[, .(mixedVar, cellType, tissue, group)])
  setkey(barPlotStats, mixedVar)
  setkey(tableToMerge, mixedVar)
  barPlotStats = merge(barPlotStats, tableToMerge, all = F)
  
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
      scale_fill_manual(values=colorScheme) + 
      scale_color_manual(values=colorScheme)
    print(jitterPlot)
  } else if (plotType == "boxPlot") {
    boxPlot = p + geom_boxplot(outlier.colour = NA, aes(fill = tissue), alpha = 0.9) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size=15)) +
      xlab("") +
      ylab("normalized signal") + 
      scale_x_discrete(labels = myLabels$spaceLabel) +
      scale_fill_manual(values=colorScheme) + 
      scale_color_manual(values=colorScheme) +
      ylim(minBoxLimit, maxBoxLimit)
    suppressWarnings(print(boxPlot))
  } else if (plotType == "barPlot") {
    barPlot = ggplot(barPlotStats, aes(x = mixedVar, y = medianBar, fill = tissue)) +
      geom_col(alpha = 0.9)+
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size=15)) +
      xlab("") +
      ylab("med (normalized signal)") + 
      scale_x_discrete(labels = myLabels$spaceLabel) +
      scale_fill_manual(values=colorScheme)
    print(barPlot)
  } else {
    stop("Plot type does not match any of the available options. Available options: jitter, boxPlot, barPlot. ")
  }
}