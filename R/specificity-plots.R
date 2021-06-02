#' The function calcOpenSignal takes the input BED file(s) 
#' in form of GRanges or GRangesList object, overlaps 
#' it with all defined open chromatin regions across 
#' cell types and returns a matrix, where each row is 
#' the input genomic region (if overlap was found), 
#' each column is a cell type, and the value 
#' is a normalized ATAC-seq signal.
#
#' @param query Genomic regions to be analyzed. Can be GRanges or GRangesList 
#'     object.
#' @param cellMatrix Matrix with open chromatin signal values, rows are genomic
#'     regions, columns are cell types. First column contains 
#'     information about the genomic region in following form: 
#'     chr_start_end. Can be either data.frame or data.table object.
#' @return  A list with named components:
#'            signalMatrix - data.table with cell specific open chromatin signal
#'                           values for query regions
#'            matrixStats - data.frame containing boxplot stats for individual 
#'                           cell type
#' 
#' @export
#' @examples
#' openRegionSummary = calcOpenSignal(vistaEnhancers, exampleOpenSignalMatrix_hg19)
calcOpenSignal = function(query, cellMatrix){
  .validateInputs(list(query=c("GRanges","GRangesList")))
  if (is(query, "GRangesList")) {
    # Recurse over each GRanges object
    regionSummaryList = lapply(query, calcOpenSignal, cellMatrix)
    nameList = names(query)
    if(is.null(nameList)) {
      nameList = seq_along(query) # Fallback to sequential numbers
    }
    # Extract signal matrices and boxplot matrices
    # rbind them and add a name column
    signalList = lapply(regionSummaryList, `[[`, 1)
    statsList = lapply(regionSummaryList, `[[`, 2)
    signalMatrix = rbindlist(signalList)
    matrixStats = rbindlist(statsList)
    signalMatrix$name = rep(nameList, vapply(signalList, nrow, integer(1)))
    matrixStats$name = rep(nameList, vapply(statsList, nrow, integer(1)))
    statNames = c("lowerWhisker", "lowerHinge", "median", 
                  "upperHinge",  "upperWhisker")
    matrixStats$boxStats = rep(statNames, length(statsList))
    
    regionSummaries = list(signalMatrix = signalMatrix, matrixStats = matrixStats)
    return(regionSummaries)
  }
  
  # if the cellMatrix is in data.frame format, convert it to data.table
  if(!is(cellMatrix, "data.table") && is(cellMatrix, "data.frame")){
    cellMatrix = as.data.table(cellMatrix)
  } else if (!is(cellMatrix, "data.table")){
    stop("The cellMatrix object is in incorrect format - must be data.table or
         data.frame.")
  }
  # get the genomic coordinates from the open chromatin signal matrix 
  openRegions = signalMatrixToOpenRegions(cellMatrix)
  
  # convert GRanges query to data.table
  queryTable = queryToDataTable(query)
  
  # find overlaps between query and cellMatrix regions
  signalMatrix = getSignalMatrix(queryTable, openRegions, cellMatrix)
  
  # get box plot statistics about the signal matrix
  matrixStatsTable = getMatrixStats(signalMatrix)
  
  # create list containing matrix with signal values, and matrix with boxplot stats
  openRegionSummary = list(signalMatrix = signalMatrix, matrixStats = matrixStatsTable)
  return(openRegionSummary)
  }

#' The function plotOpenSignal visualizes the signalMatrix obtained from
#' calcOpenSignal.
#'
#' @param openRegionSummary Output list from \code{calcOpenSignal} function.
#' @param  plotType Options are: "jitter" - jitter plot with box plot on top,
#'     "boxPlot" - box plot without individual points and outliers,
#'     "barPlot" (default) - bar height represents the median signal value
#'     for a given cell type, 
#'     "violinPlot" - violin plot with medians.
#' @param cellGroup - This option allows to selcet a tissue type to be 
#'     plotted, if NA (default) all available tissue types are ploted, 
#'     available options: {"blood", "bone", "CNS", "embryonic", "eye", 
#'     "foreskin", "gastrointestinal", "heart", "liver", "lymphatic", 
#'     "mammaryGland", "mouth", "respiratorySystem", "skeletalMuscle",
#'     "skin", "urinarySystem", "vasculature"}, can be passed as a 
#'     character string or vector of strings.
#' @param cellTypeMetadata Metadata for cell type - tissue association. This
#'     option is for users, who provide their own open region signal 
#'     matrix. The cellTypeMetadata matrix must contain two columns called
#'     cellType and tissueType. cellType column containes the cell type names 
#'     in the provided signalMatrix column names. The tissueType columns 
#'     provides an information, which tissue the cell type comes from.
#' @param colorScheme Provide color values for each tissueType if you want to 
#'     change the default colors.
#' @return A ggplot object.
#' 
#' @export
#' @examples
#' openRegionSummary = calcOpenSignal(vistaEnhancers, exampleOpenSignalMatrix_hg19)
#' plotSignal = plotOpenSignal(openRegionSummary)
#' plotSignal = plotOpenSignal(openRegionSummary, plotType = "jitter", 
#' cellGroup = "blood")
plotOpenSignal = function(openRegionSummary, 
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
  
  # reshape signal matrix and stats into matrix usable by ggplot
  # and attach cell type metadata
  plotTables = reshapeDataToPlot(openRegionSummary, cellTypeMetadata)
  
  # if user defines cell group, filter the data
  filteredPlotTables = filterGroups(cellGroup, plotTables)
  plotSignalMatrix = filteredPlotTables$plotSignalMatrix
  plotBoxStats = filteredPlotTables$plotBoxStats
  
  # arrange labels in a way, that corresponding groups are plotted together
  myLabels = unique(plotSignalMatrix[, .(mixedVar, cellType)])
  myLabels[, spaceLabel := gsub("_", " ", myLabels[, cellType])]
  
  # do the plotting
  if (plotType == "jitter"){
    
    jitterPlot = OpenSignalJitterPlot(plotSignalMatrix, myLabels, colorScheme)
    return(jitterPlot)
    
  } else if (plotType == "boxPlot") {
    
    boxPlot = OpenSignalBoxPlot(plotSignalMatrix, plotBoxStats, myLabels, colorScheme)
    return(boxPlot)
    
  } else if (plotType == "barPlot") {
    
    barPlot = OpenSignalBarPlot(plotBoxStats, myLabels, colorScheme)
    return(barPlot)
    
  } else if (plotType == "violinPlot"){
    
    violinPlot = OpenSignalViolinPlot(plotSignalMatrix,plotBoxStats, myLabels, colorScheme)
    return(violinPlot)
    
  } else {
    stop("Plot type does not match any of the available options. 
         Available options: jitter, boxPlot, barPlot, violinPlot. ")
  }
  }


# Internal helper function to extract genomic regions
# from matrix with open chromatin signal values
#
# @param cellMatrix Matrix with open chromatin signal values, rows are genomic
#     regions, columns are cell types. First column contains 
#     information about the genomic region in following form: 
#     chr_start_end. Can be either data.frame or data.table object.
# @return A data table with genomic coordinates of regions defined 
# in 'cellMatrix'
signalMatrixToOpenRegions = function(cellMatrix){
  # get the genomic coordinates from the open chromatin signal matrix - 
  # the first column convert the chr_start_end into a three column data.table
  openRegions = cellMatrix[,1]
  colnames(openRegions) = "V1"
  openRegions[, c("chr", "start", "end") := tstrsplit(V1, "_", fixed=TRUE)]
  numericColumns = c("start", "end")
  openRegions[, (numericColumns ) := lapply(.SD, as.numeric), 
              .SDcols = numericColumns]
  openRegions = openRegions[,.(chr, start, end)]
  return(openRegions)
}

# Internal helper function to convert GRanges query into 
# data.table and adds 4th column in form chr_start_end
#
# @param query Genomic regions to be analyzed in form of GRanges object.
# @return A data table with genomic coordinates from 'query' and 
# 4th column containing peak identifiers
queryToDataTable = function(query){
  # convert GRanges query to data.table
  queryTable = as.data.table(data.frame(chr = seqnames(query),
                                        start = start(query),
                                        end = end(query)))
  
  # select just the fist three columns and give them name chr, start, end
  # create a 4th column with peak name in following format: chr_start_end
  queryTable = queryTable[, c(1, 2, 3)]
  colnames(queryTable) = c("chr", "start", "end")
  queryTable[, peakName := paste(queryTable[,chr], 
                                 queryTable[,start], 
                                 queryTable[,end], 
                                 sep = "_")]
  return(queryTable)
}

# Internal helper function to overlap two data.table objects and
# and summarize signal values from cellMatrix which fall into 
# queryTable regions
#
# @param queryTable data.table created by \code{queryToDataTable}
#  function
# @param openRegions data.table created by \code{signalMatrixToOpenRegions}
#  function
# @param cellMatrix Matrix with open chromatin signal values, rows are genomic
#     regions, columns are cell types. First column contains 
#     information about the genomic region in following form: 
#     chr_start_end. Can be either data.frame or data.table object.
# @return A data.table with signal values for 
getSignalMatrix = function(queryTable, openRegions, cellMatrix){
  # find overlaps between query and cellMatrix regions
  setkey(queryTable, chr, start, end)
  overlaps = foverlaps(openRegions, queryTable, which=TRUE)
  overlaps = na.omit(overlaps)
  
  # extract the regions which overlap with query, assign the query peaks to 
  # them and calculate the sum of the signal within these regions
  signalMatrix = cellMatrix[overlaps[,xid]] 
  signalMatrix = signalMatrix[,-1]
  signalMatrix[, queryPeak := queryTable[overlaps[,yid], peakName]]
  signalMatrix = signalMatrix[, lapply(.SD, max), by = .(queryPeak)]
  
  return(signalMatrix)
}

# Internal helper function to obtain summary statistics about each 
# column in the matrix
# 
# @param signalMatrix Data.table object conatianing signal values 
# for genomic regions (rows) across cell types (columns). The first
# column carries information about the genomic regions in form
# chr_start_end
# @return A data frame with boxplot statistics
getMatrixStats = function(signalMatrix){
  # get box plot statistics about the signal matrix
  matrixStatsList = apply(signalMatrix[,-1], 2, boxplot.stats)
  matrixStats = unlist(lapply(matrixStatsList, `[[`, 1))
  matrixStatsTable = data.frame(matrix(matrixStats, nrow = 5))
  rownames(matrixStatsTable) = c("lowerWhisker", "lowerHinge", "median", "upperHinge", "upperWhisker")
  colnames(matrixStatsTable) = names(matrixStatsList)
  
  return(matrixStatsTable)
}

# Internal helper function to reshape signal matrix and stats matrix
# into a form, which can be plotted by ggplot
# 
# @param openRegionSummary Results from \code{calcOpenSignal}.
# @param cellTypeMetadata Metadata for cell type - tissue association. Reqired
#        input in \code{plotOpenSignal} function
# 
# @return A list with two matrices reshaped for ggplot functions
reshapeDataToPlot = function(openRegionSummary, cellTypeMetadata){
  # reshape the signal matrix ans boxplotStats matrices into ggplot usable form
  # attach the metadata for coloring sort table alphabetically by 
  # tissueType-cellType
  signalMatrix = openRegionSummary[["signalMatrix"]]
  boxStats = openRegionSummary[["matrixStats"]]
  
  if ("name" %in% names(signalMatrix)){
    plotSignalMatrix = reshape2::melt(signalMatrix, 
                                      id.vars=c("queryPeak", "name"), 
                                      variable.name="cellType", value.name="signal")
    plotBoxStats = reshape2::melt(boxStats, 
                                  id.vars=c("boxStats", "name"), 
                                  variable.name="cellType", value.name="value")
  } else {
    plotSignalMatrix = reshape2::melt(signalMatrix, id.vars="queryPeak", 
                                      variable.name="cellType", value.name="signal")
    boxStats$boxStats = rownames(boxStats)
    plotBoxStats = reshape2::melt(boxStats, id.vars="boxStats",
                                  variable.name="cellType", value.name="value")
  }
  
  data.table::setkey(cellTypeMetadata, cellType)
  data.table::setDT(plotSignalMatrix)
  data.table::setkey(plotSignalMatrix, cellType)

  plotSignalMatrix = plotSignalMatrix[cellTypeMetadata, on = "cellType", nomatch=0]
  plotSignalMatrix[, lowerCaseTissue := tolower(tissueType)]
  data.table::setorder(plotSignalMatrix, lowerCaseTissue, cellType)
  plotSignalMatrix[, mixedVar := paste(plotSignalMatrix[,tissueType], 
                                       plotSignalMatrix[,cellType], sep="_")]
  
  data.table::setDT(plotBoxStats)
  data.table::setkey(plotBoxStats, cellType)
  plotBoxStats = plotBoxStats[cellTypeMetadata, on = "cellType", nomatch=0]
  plotBoxStats[, lowerCaseTissue := tolower(tissueType)]
  data.table::setorder(plotBoxStats, lowerCaseTissue, cellType)
  plotBoxStats[, mixedVar := paste(plotBoxStats[,tissueType], 
                                   plotBoxStats[,cellType], sep ="_")]
  
  plotTables = list(plotSignalMatrix = plotSignalMatrix, plotBoxStats = plotBoxStats)
  return(plotTables)
}

# Internal helper function to filter out tissue type predefined
# by user
# 
# @param cellGroup Input in \code{plotOpenSignal} function specifying
#        tissue type to be plotted
# @param plotTables Result from \code{reshapeDataToPlot} function
# 
# @return A list with two matrices filtered for user defined tissue type
filterGroups = function(cellGroup, plotTables){
  plotSignalMatrix = plotTables$plotSignalMatrix
  plotBoxStats = plotTables$plotBoxStats
  
  if (length(cellGroup) == 1){
    if (is.na(cellGroup)){
      plotSignalMatrix = plotSignalMatrix
      plotBoxStats = plotBoxStats
    } else if (cellGroup %in% levels(factor(cellTypeMetadata$tissueType))){
      plotSignalMatrix = plotSignalMatrix[tissueType == cellGroup]
      plotBoxStats = plotBoxStats[tissueType == cellGroup]
    } else {
      stop("The input cell group is not in predefined list of options.")
    }
  } else if (all(cellGroup %in% levels(factor(cellTypeMetadata$tissueType)))) {
    plotSignalMatrix = plotSignalMatrix[tissueType %in% cellGroup]
    plotBoxStats = plotBoxStats[tissueType %in% cellGroup]
  } else {
    stop("At least one of the input cell groups is not in predefined list of 
         options.")
  }
  
  plotTables = list(plotSignalMatrix = plotSignalMatrix, plotBoxStats = plotBoxStats)
  return(plotTables)
  }

# Internal helper function to plot jitter plot
# 
# @param plotSignalMatrix Result from \code{reshapeDataToPlot} function followed 
#        by \code{filterGroups} function
# @param myLabels Labels created to group corresponding tissue type together
# @param colorScheme Colors for tissue types - input to \code{plotOpenSignal}
#        function
#
# @return ggplot object - jitter plot
OpenSignalJitterPlot = function(plotSignalMatrix, myLabels, colorScheme){
  
  jitterPlot = ggplot(plotSignalMatrix, aes(x = mixedVar, y = signal)) +
    theme_blank_facet_label() +
    theme(strip.text.y.right = element_text(angle = 0))
  
  if ("name" %in% names(plotSignalMatrix)){
    jitterPlot = jitterPlot + facet_grid(name ~ .)
  }
  
  jitterPlot = jitterPlot + 
    geom_jitter(alpha=0.5, height=0, width=0.35, aes(color=tissueType)) +
    geom_boxplot(outlier.colour=NA, fill=NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          text = element_text(size=10)) +
    xlab("") +
    ylab("normalized signal") + 
    scale_x_discrete(labels=myLabels$spaceLabel) +
    scale_fill_manual(values=colorScheme) + 
    scale_color_manual(values=colorScheme)
  return(jitterPlot)
}

# Internal helper function to plot box plot
# 
# @param plotSignalMatrix Result from \code{reshapeDataToPlot} function followed 
#        by \code{filterGroups} function
# @param plotBoxStats Result from \code{reshapeDataToPlot} function followed 
#        by \code{filterGroups} function
# @param myLabels Labels created to group corresponding tissue type together
# @param colorScheme Colors for tissue types - input to \code{plotOpenSignal}
#        function
#
# @return ggplot object - box plot
OpenSignalBoxPlot = function(plotSignalMatrix, plotBoxStats, myLabels, colorScheme){
  
  boxPlot = ggplot(plotSignalMatrix, aes(x = mixedVar, y = signal)) +
    theme_blank_facet_label() +
    theme(strip.text.y.right = element_text(angle = 0))
  
  if ("name" %in% names(plotSignalMatrix)){
    boxPlot = boxPlot + facet_grid(name ~ .)
  }
  boxPlot = boxPlot + geom_boxplot(outlier.colour=NA, 
                                   aes(fill=tissueType), alpha=0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          text = element_text(size=10)) +
    xlab("") +
    ylab("normalized signal") + 
    scale_x_discrete(labels=myLabels$spaceLabel) +
    scale_fill_manual(values=colorScheme) + 
    scale_color_manual(values=colorScheme) +
    ylim(min(plotBoxStats$value), max(plotBoxStats$value))
  return(boxPlot)
}

# Internal helper function to plot bar plot
# 
# @param plotBoxStats Result from \code{reshapeDataToPlot} function followed 
#        by \code{filterGroups} function
# @param myLabels Labels created to group corresponding tissue type together
# @param colorScheme Colors for tissue types - input to \code{plotOpenSignal}
#        function
#
# @return ggplot object - bar plot

OpenSignalBarPlot = function(plotBoxStats, myLabels, colorScheme){
  barPlot = ggplot(plotBoxStats[boxStats == "median"], 
                   aes(x=mixedVar, 
                       y=value, 
                       fill=tissueType))
  
  if ("name" %in% names(plotBoxStats)){
    barPlot = barPlot + facet_grid(name ~ .)
  }
  barPlot = barPlot +
    geom_col(alpha=0.9)+
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          text = element_text(size=10)) +
    xlab("") +
    ylab("med (normalized signal)") + 
    scale_x_discrete(labels=myLabels$spaceLabel) +
    scale_fill_manual(values=colorScheme)  +
    theme_blank_facet_label() +
    theme(strip.text.y.right = element_text(angle=0))
  return(barPlot)
}


# Internal helper function to plot bar plot
# 
# @param plotSignalMatrix Result from \code{reshapeDataToPlot} function followed 
#        by \code{filterGroups} function
# @param plotBoxStats Result from \code{reshapeDataToPlot} function followed 
#        by \code{filterGroups} function
# @param myLabels Labels created to group corresponding tissue type together
# @param colorScheme Colors for tissue types - input to \code{plotOpenSignal}
#        function
#
# @return ggplot object - violin plot

OpenSignalViolinPlot = function(plotSignalMatrix,plotBoxStats, myLabels, colorScheme){
  violinPlot = ggplot(plotSignalMatrix, aes(x = mixedVar, y = signal)) +
    theme_blank_facet_label() +
    theme(strip.text.y.right = element_text(angle = 0))
  
  if ("name" %in% names(plotSignalMatrix)){
    violinPlot = violinPlot + facet_grid(name ~ .)
  }
  violinPlot = violinPlot + 
    geom_violin(aes(fill=tissueType), alpha=0.8, scale = "width", trim = TRUE) +
    geom_point(data = plotBoxStats[plotBoxStats$boxStats == "median",], 
               aes(x = mixedVar, y = value), color = "black", size = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          text = element_text(size=10)) +
    xlab("") +
    ylab("normalized signal") + 
    scale_x_discrete(labels=myLabels$spaceLabel) +
    scale_fill_manual(values=colorScheme) + 
    scale_color_manual(values=colorScheme)
  return(violinPlot)
}
