#' The function calcOpenSignal takes the input BED file(s) in form of GRanges 
#' or GRangesList object, overlaps it with all defined open chromatin regions 
#' across cell types and returns a matrix, where each row is the input genomic
#' region (if overlap was found), each column is a cell type, and the value 
#' is a normalized ATAC-seq signal.
#
#' @param query Genomic regions to be analyzed. Can be GRanges or GRangesList 
#'     object.
#' @param cellMatrix Matrix with open chromatin signal values, rows are genomic
#'     regions, columns are cell types. First column contains 
#'     information about the genomic region in following form: 
#'     chr_start_end. Can be either data.frame or data.table object.
#' @return  A data.table with cell specific open chromatin signal values for 
#'     query regions.
#' 
#' @export
#' @examples
#' signalMatrix = calcOpenSignal(vistaEnhancers, exampleOpenSignalMatrix_hg19)
calcOpenSignal = function(query, cellMatrix){
  .validateInputs(list(query=c("GRanges","GRangesList")))
  if (is(query, "GRangesList")) {
    # Recurse over each GRanges object
    x = lapply(query, calcOpenSignal, cellMatrix)
    nameList = names(query)
    if(is.null(nameList)) {
      nameList = seq_along(query) # Fallback to sequential numbers
    }
    # Append names
    xb = rbindlist(x)
    xb$name = rep(nameList, vapply(x, nrow, integer(1)))
    return(xb)
  }
  
  # if the cellMatrix is in data.frame format, convert it to data.table
  if(!is(cellMatrix, "data.table") && is(cellMatrix, "data.frame")){
    cellMatrix = as.data.table(cellMatrix)
  } else if (!is(cellMatrix, "data.table")){
    stop("The cellMatrix object is in incorrect format - must be data.table or
         data.frame.")
  }
  
  # get the genomic coordinates from the open chromatin signal matrix - 
  # the first column convert the chr_start_end into a three column data.table
  openRegions = cellMatrix[,1]
  colnames(openRegions) = "V1"
  openRegions[, c("chr", "start", "end") := tstrsplit(V1, "_", fixed=TRUE)]
  numericColumns = c("start", "end")
  openRegions[, (numericColumns ) := lapply(.SD, as.numeric), 
              .SDcols = numericColumns]
  openRegions = openRegions[,.(chr, start, end)]
  
  # check if query BED file is GRanges object and convert it to data.table
  # otherwise give error
  query = as.data.table(data.frame(chr = seqnames(query),
                                   start = start(query),
                                   end = end(query)))
  
  # select just the fist three columns and give them name chr, start, end
  # create a 4th column with peak name in following format: chr_start_end
  query = query[, 1:3]
  colnames(query) = c("chr", "start", "end")
  query[, peakName := paste(query[,chr], query[,start], query[,end], 
                            sep = "_")]
  
  # find which regions of the genomic regions of interest overlap with the open
  # chromatin signal matrix regions
  setkey(query, chr, start, end)
  overlaps = foverlaps(openRegions,query, which = TRUE)
  overlaps = na.omit(overlaps)
  
  # extract the regions which overlap with query, assign the query peaks to 
  # them and calculate the sum of the signal within these regions
  signalMatrix = cellMatrix[overlaps[,xid]] 
  signalMatrix = signalMatrix[,-1]
  signalMatrix[, queryPeak := query[overlaps[,yid], peakName]]
  signalMatrix = signalMatrix[, lapply(.SD, sum), by = .(queryPeak)]
  
  return(signalMatrix)
}

#' The function plotOpenSignal visualizes the signalMatrix obtained from
#' calcOpenSignal.
#'
#' @param signalMatrix Output data.table from \code{calcOpenSignal} function.
#' @param  plotType Options are: jitter - jitter plot with box plot on top
#'     boxPlot - box plot without individual points and outliers
#'     barPlot (default) - bar height represents the median signal value
#'     for a given cell type.
#' @param cellGroup - This option allows to selcet a group of cells to be 
#'     plotted, if NA (default) all available cell groups are ploted, 
#'     available options: {"blood", "bone", "CNS", "embryonic", "eye", 
#'     "foreskin", "gastrointestinal", "heart", "liver", "lymphatic", 
#'     "mammaryGland", "mouth", "respiratorySystem", "skeletalMuscle",
#'     "skin", "urinarySystem", "vasculature"}, can be passed as a 
#'     character string or vector of strings.
#' @param cellTypeMetadata Metadata for cell type - tissue association. This
#'     option is for users, who provide their own open region signal 
#'     matrix. The cellTypeMetadata matrix must contain two columns called
#'     cellType and tissue. cellType column containes the cell type names 
#'     in the provided signalMatrix column names. The tissue columns 
#'     provides an information, which tissue the cell type comes from.
#' @param colorScheme Provide color values for each tissue if you want to 
#'     change the default colors.
#' @return A ggplot object.
#' 
#' @export
#' @examples
#' \dontrun{
#' signalMatrix = calcOpenSignal(vistaEnhancers, exampleOpenSignalMatrix_hg19)
#' plotSignal = plotOpenSignal(signalMatrix)
#' plotSignal = plotOpenSignal(signalMatrix, plotType = "jitter", cellGroup = "blood")
#' }
plotOpenSignal = function(signalMatrix, 
                          plotType = "barPlot", 
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
    cellTypeMetadata = getReferenceData("","cellTypeMetadata")
  } 
  
  # reshape the signal matrix into ggplot usable form 
  # attach the metadata for coloring
  # sort table alphabetically by tissue-cellType
  if ("name" %in% names(signalMatrix)){
    plotSignalMatrix = reshape2::melt(signalMatrix, id.vars = c("queryPeak", "name"), 
                            variable.name = "cellType", value.name = "signal")
  } else {
    plotSignalMatrix = reshape2::melt(signalMatrix, id.vars = "queryPeak", 
                            variable.name = "cellType", value.name = "signal")
  }
  data.table::setDT(plotSignalMatrix)
  data.table::setkey(plotSignalMatrix, cellType)
  data.table::setkey(cellTypeMetadata, cellType)
  plotSignalMatrix = merge(plotSignalMatrix, cellTypeMetadata, all = FALSE)
  plotSignalMatrix[, lowerCaseTissue := tolower(tissue)]
  data.table::setorder(plotSignalMatrix, lowerCaseTissue, cellType)
  plotSignalMatrix[, mixedVar := paste(plotSignalMatrix[,tissue], 
                                       plotSignalMatrix[,cellType], sep = "_")]
  
  
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
    stop("At least one of the input cell groups is not in predefined list of 
         options.")
  }
  
  # arrange labels in a way, that corresponding groups are plotted together
  myLabels = unique(plotSignalMatrix[, .(mixedVar, cellType)])
  myLabels[, spaceLabel := gsub("_", " ", myLabels[, cellType])]
  
  # get box plot staistics to set plot limits
  boxStats = plotSignalMatrix[, .(boxStats = 
                                    list(boxplot.stats(signal)$stats)),
                              by = mixedVar]
  minBoxLimit = min(unlist(boxStats$boxStats))
  maxBoxLimit = max(unlist(boxStats$boxStats))
  
  # get mediand of signal values to make a bar plot
  if ("name" %in% names(plotSignalMatrix)){
    barPlotStats = plotSignalMatrix[, .(medianBar = median(signal)), 
                                    by = c("mixedVar", "name")]
  } else {
    barPlotStats = plotSignalMatrix[, .(medianBar = median(signal)), 
                                    by = mixedVar]
  }
  tableToMerge = unique(plotSignalMatrix[, .(mixedVar, cellType, 
                                             tissue, group)])
  setkey(barPlotStats, mixedVar)
  setkey(tableToMerge, mixedVar)
  barPlotStats = merge(barPlotStats, tableToMerge, all = FALSE)
  
  # do the plotting
  p = ggplot(plotSignalMatrix,
             aes(x = mixedVar, y = signal))
  
  if (plotType == "jitter"){
    
    if ("name" %in% names(plotSignalMatrix)){
      p = p + facet_grid(name ~ .)
    }
    
    jitterPlot = p + 
      geom_jitter(alpha = 0.5, height = 0, width = 0.35, aes(color = tissue)) +
      geom_boxplot(outlier.colour = NA, fill = NA) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size=10)) +
      xlab("") +
      ylab("normalized signal") + 
      scale_x_discrete(labels = myLabels$spaceLabel) +
      scale_fill_manual(values=colorScheme) + 
      scale_color_manual(values=colorScheme)
    return(jitterPlot)
  } else if (plotType == "boxPlot") {
    
    if ("name" %in% names(plotSignalMatrix)){
      p = p + facet_grid(name ~ .)
    }
    boxPlot = p + geom_boxplot(outlier.colour = NA, aes(fill = tissue), alpha = 0.9) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size=10)) +
      xlab("") +
      ylab("normalized signal") + 
      scale_x_discrete(labels = myLabels$spaceLabel) +
      scale_fill_manual(values=colorScheme) + 
      scale_color_manual(values=colorScheme) +
      ylim(minBoxLimit, maxBoxLimit)
    return(boxPlot)
  } else if (plotType == "barPlot") {
    barPlot = ggplot(barPlotStats, aes(x = mixedVar, y = medianBar, fill = tissue))
    
    if ("name" %in% names(barPlotStats)){
      barPlot = barPlot + facet_grid(name ~ .)
    }
    barPlot = barPlot +
      geom_col(alpha = 0.9)+
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size=10)) +
      xlab("") +
      ylab("med (normalized signal)") + 
      scale_x_discrete(labels = myLabels$spaceLabel) +
      scale_fill_manual(values=colorScheme)
    return(barPlot)
  } else {
    stop("Plot type does not match any of the available options. 
         Available options: jitter, boxPlot, barPlot. ")
  }
}
