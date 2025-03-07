#' The function calcSummarySignal takes the input BED file(s) 
#' in form of GRanges or GRangesList object, overlaps 
#' it with all defined open chromatin regions across 
#' conditions (e.g. cell types) and returns a matrix, 
#' where each row is the input genomic region 
#' (if overlap was found), each column is a condition, 
#' and the value is a meam signal from regions where
#' overlap was found.
#
#' @param query Genomic regions to be analyzed. Can be GRanges or GRangesList 
#'     object.
#' @param signalMatrix Matrix with signal values in predfined regions, where
#'     rows are predefined genomic regions, columns are conditions 
#'     (e.g. cell types in which the signal was measured). 
#'     First column contains information about the genomic region in 
#'     following form: chr_start_end. 
#'     Can be either data.frame or data.table object.
#' @param openRegions Extracted genomic regions by signalMatrixToOpenRegions.
#'      It is only calculated for the first bed file.
#' @return  A list with named components:
#'            signalSummaryMatrix - data.table with cell specific open chromatin signal
#'                           values for query regions
#'            matrixStats - data.frame containing boxplot stats for individual 
#'                           cell type
#' 
#' @export
#' @examples
#' signalSummaryList = calcSummarySignal(vistaEnhancers, exampleOpenSignalMatrix_hg19)
calcSummarySignal = function(query, signalMatrix, openRegions = NULL) {
  .validateInputs(list(query = c("GRanges", "GRangesList")))

  if (is(query, "GRangesList") && is.null(openRegions)) {
    openRegions <- signalMatrixToOpenRegions(signalMatrix)
  }

  if (is(query, "GRangesList")) {
    regionSummaryList = lapply(query, function(q) calcSummarySignal(q, signalMatrix, openRegions))  
    nameList = names(query)
    
    if (is.null(nameList)) nameList = seq_along(query)
    
    signalSummaryMatrix = rbindlist(lapply(regionSummaryList, `[[`, 1), use.names = TRUE, fill = TRUE)
    matrixStats = rbindlist(lapply(regionSummaryList, `[[`, 2), use.names = TRUE, fill = TRUE)
    
    signalSummaryMatrix[, name := rep(nameList, vapply(regionSummaryList, function(x) nrow(x[[1]]), integer(1)))]
    matrixStats[, name := rep(nameList, vapply(regionSummaryList, function(x) nrow(x[[2]]), integer(1)))]

    statNames = c("lowerWhisker", "lowerHinge", "median", "upperHinge", "upperWhisker")
    matrixStats[, boxStats := rep(statNames, length(regionSummaryList))]
    
    return(list(signalSummaryMatrix = signalSummaryMatrix, matrixStats = matrixStats))
  }

  if (!is(signalMatrix, "data.table")) {
    if (is(signalMatrix, "data.frame")) {
      setDT(signalMatrix)
    } else {
      stop("The signalMatrix object is in incorrect format - must be data.table or data.frame.")
    }
  }

  if (nrow(openRegions) != nrow(signalMatrix)) {
    openRegions[, peakName := sprintf("%s_%s_%s", chr, start, end)]
    signalMatrix = signalMatrix[V1 %in% openRegions$peakName]
  }
  
  queryTable = queryToDataTable(query)
  
  signalSummaryMatrix = getSignalMatrix(queryTable, openRegions, signalMatrix)
  matrixStatsTable = getMatrixStats(signalSummaryMatrix)
  
  return(list(signalSummaryMatrix = signalSummaryMatrix, matrixStats = matrixStatsTable))
}

#' The function plotSummarySignal visualizes the signalSummaryMatrix obtained from
#' \code{calcSummarySignal}.
#'
#' @param signalSummaryList Output list from \code{calcSummarySignal} function.
#' @param  plotType Options are: "jitter" - jitter plot with box plot on top,
#'     "boxPlot" - box plot without individual points and outliers,
#'     "barPlot" (default) - bar height represents the median signal value
#'     for a given cell type, 
#'     "violinPlot" - violin plot with medians.
#' @param metadata (optional) data.table used for grouping columns from 
#'    'signalMatrix' into categories, that are then plotted with different colors. 
#'    Must contain variable 'colName' that contains all the condition column names 
#'    from 'signaMatrix'.
#' @param colorColumn (optional only if metadata provided) columns name from 
#'    'metadata' table that will be used as grouping variable for coloring.
#' @param filterGroupColumn (optional only if metadata provided and 
#'    'filterGroup' specified) allows user to plot specified subgroups only. 
#'    String specifying the column name in 'metadata' from which groups will 
#'    be filtered (groups are specified in as 'filterGroups)
#' @param filterGroup (optional only if 'metadata' and 'filterGroupColumn' 
#'    provided) - string (or vector of strings) of groups from 
#'    'filterGroupColumn' to be plottted.
#' @return A ggplot object.
#' 
#' @export
#' @examples
#' signalSummaryList = calcSummarySignal(vistaEnhancers, exampleOpenSignalMatrix_hg19)
#' metadata = cellTypeMetadata
#' plotSignal = plotSummarySignal(signalSummaryList)
#' 
#' plotSignalTissueColor = plotSummarySignal(signalSummaryList = signalSummaryList, 
#' plotType = "jitter", metadata = metadata, colorColumn = "tissueType")
#' 
#' plotSignalFiltered = plotSummarySignal(signalSummaryList = signalSummaryList,
#' plotType = "violinPlot", metadata = metadata, colorColumn = "tissueType", 
#' filterGroupColumn = "tissueType", filterGroup = c("skin", "blood"))
plotSummarySignal = function(signalSummaryList, 
                             plotType="barPlot", 
                             metadata=NULL, 
                             colorColumn=NULL,
                             filterGroupColumn=NULL,
                             filterGroup=NULL){
  
  if(is.null(metadata) & (any(c(!is.null(colorColumn), 
                                !is.null(filterGroupColumn), 
                                !is.null(filterGroup))))){
    stop("In order to color or filter groups, you must provide metadata.")
  }
  
  # reshape signal matrix and stats into matrix usable by ggplot
  # and attach cell type metadata
  plotTables = reshapeDataToPlot(signalSummaryList, metadata, colorColumn)
  
  # if user defines cell group, filter the data
  filteredPlotTables = filterGroups(plotTables, metadata, 
                                    filterGroupColumn, filterGroup)
  plotSignalMatrix = filteredPlotTables$plotSignalMatrix
  plotBoxStats = filteredPlotTables$plotBoxStats
  
  # do the plotting
  if (plotType == "jitter"){
    
    jitterPlot = SummarySignalJitterPlot(plotSignalMatrix, colorColumn)
    return(jitterPlot)
    
  } else if (plotType == "boxPlot") {
    
    boxPlot = SummarySignalBoxPlot(plotSignalMatrix, plotBoxStats, colorColumn)
    return(boxPlot)
    
  } else if (plotType == "barPlot") {
    
    barPlot = SummarySignalBarPlot(plotBoxStats, colorColumn)
    return(barPlot)
    
  } else if (plotType == "violinPlot"){
    
    violinPlot = SummarySignalViolinPlot(plotSignalMatrix,plotBoxStats, colorColumn)
    return(violinPlot)
    
  } else {
    stop("Plot type does not match any of the available options. 
         Available options: jitter, boxPlot, barPlot, violinPlot. ")
  }
  }

# Internal helper function to extract genomic regions
# from matrix with open chromatin signal values
#
# @param signalMatrix Matrix with signal values in predfined regions, where
#     rows are predefined genomic regions, columns are conditions 
#     (e.g. cell types in which the signal was measured). 
#     First column contains information about the genomic region in 
#     following form: chr_start_end. 
#     Can be either data.frame or data.table object.
# @return A data table with genomic coordinates of regions defined 
# in 'signalMatrix'
signalMatrixToOpenRegions = function(signalMatrix) {

  openRegions = setDT(tstrsplit(signalMatrix[[1]], '[:_]', type.convert = TRUE, keep = 1:3))
  
  setnames(openRegions, c("chr", "start", "end"))
  
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
# and summarize signal values from signalMatrix which fall into 
# queryTable regions
#
# @param queryTable data.table created by \code{queryToDataTable}
#  function
# @param openRegions data.table created by \code{signalMatrixToOpenRegions}
#  function
# @param signalMatrix Matrix with signal values in predfined regions, where
#     rows are predefined genomic regions, columns are conditions 
#     (e.g. cell types in which the signal was measured). 
#     First column contains information about the genomic region in 
#     following form: chr_start_end. 
#     Can be either data.frame or data.table object.
# @return A data.table with signal values for 
getSignalMatrix = function(queryTable, openRegions, signalMatrix) {

  setkey(openRegions, chr, start, end) 
  setkey(queryTable, chr, start, end)

  overlaps = foverlaps(openRegions, queryTable, which=TRUE)
  overlaps = na.omit(overlaps)

  overlaps[, queryPeak := queryTable$peakName[yid]]

  signalSummaryMatrix = signalMatrix[overlaps$xid]
  signalSummaryMatrix[, queryPeak := overlaps$queryPeak]
  signalSummaryMatrix = signalSummaryMatrix[, lapply(.SD, max), by = queryPeak, .SDcols = setdiff(names(signalSummaryMatrix), c("V1", "queryPeak"))]

  return(signalSummaryMatrix)
}

# Internal helper function to obtain summary statistics about each 
# column in the matrix
# 
# @param signalSummaryMatrix Data.table object conatianing signal values 
# for genomic regions (rows) across cell types (columns). The first
# column carries information about the genomic regions in form
# chr_start_end
# @return A data frame with boxplot statistics
getMatrixStats = function(signalSummaryMatrix){
  # get box plot statistics about the signal matrix
  matrixStatsList = apply(signalSummaryMatrix[,-1], 2, boxplot.stats)
  matrixStats = unlist(lapply(matrixStatsList, `[[`, 1))
  matrixStatsTable = data.frame(matrix(matrixStats, nrow = 5))
  rownames(matrixStatsTable) = c("lowerWhisker", "lowerHinge", "median", "upperHinge", "upperWhisker")
  colnames(matrixStatsTable) = names(matrixStatsList)
  
  return(matrixStatsTable)
}


# Internal helper function to reshape signal matrix and stats matrix
# into a form, which can be plotted by ggplot
# 
# @param signalSummaryList Results from \code{calcSummarySignal}.
# @param metadata Metadata as provided to \code{plotSummarySignal}
# @param colorColumn Column name specifying varibale for coloring 
#   as provided to \{plotSummarySignal}
# 
# @return A list with two matrices reshaped for ggplot functions
reshapeDataToPlot = function(signalSummaryList, metadata, colorColumn){
  # reshape the signal matrix ans boxplotStats matrices into ggplot usable form
  # attach the metadata for coloring sort table alphabetically by 
  # tissueType-cellType
  signalSummaryMatrix = signalSummaryList[["signalSummaryMatrix"]]
  boxStats = signalSummaryList[["matrixStats"]]
  
  if ("name" %in% names(signalSummaryMatrix)){
    plotSignalMatrix = reshape2::melt(signalSummaryMatrix, 
                                      id.vars=c("queryPeak", "name"), 
                                      variable.name="colName", value.name="signal")
    plotBoxStats = reshape2::melt(boxStats, 
                                  id.vars=c("boxStats", "name"), 
                                  variable.name="colName", value.name="value")
  } else {
    plotSignalMatrix = reshape2::melt(signalSummaryMatrix, id.vars="queryPeak", 
                                      variable.name="colName", value.name="signal")
    boxStats$boxStats = rownames(boxStats)
    plotBoxStats = reshape2::melt(boxStats, id.vars="boxStats",
                                  variable.name="colName", value.name="value")
  }
  
  if(!is.null(metadata)){
    data.table::setkey(metadata, colName)
    
    data.table::setDT(plotSignalMatrix)
    data.table::setkey(plotSignalMatrix, colName)
    # attach metadata table to the plotSignalMatrix
    plotSignalMatrix = plotSignalMatrix[metadata, on = "colName", nomatch=0]
    
    data.table::setDT(plotBoxStats)
    data.table::setkey(plotBoxStats, colName)
    plotBoxStats = plotBoxStats[metadata, on = "colName", nomatch=0]
  }
  
  if(!is.null(colorColumn)){
    plotSignalMatrix[, lowerColorColumn := tolower(get(colorColumn))]
    data.table::setorder(plotSignalMatrix, lowerColorColumn, colName)
    plotSignalMatrix[, mixedVar := paste(plotSignalMatrix[,get(colorColumn)], 
                                         plotSignalMatrix[,colName], sep="_")]
    
    plotBoxStats[, lowerColorColumn := tolower(get(colorColumn))]
    data.table::setorder(plotBoxStats, lowerColorColumn, colName)
    plotBoxStats[, mixedVar := paste(plotBoxStats[,get(colorColumn)], 
                                     plotBoxStats[,colName], sep ="_")]
  }
  
  plotTables = list(plotSignalMatrix = plotSignalMatrix, plotBoxStats = plotBoxStats)
  return(plotTables)
}

# Internal helper function to filter out tissue type predefined
# by user
# @param plotTables Result from \code{reshapeDataToPlot} function
# @param metadata Metadata as provided to \code{plotSummarySignal}
# @param filterGroupColumn Input to \code{plotSummarySignal} specifying
#    from which 'metadata' column will be the filtering done.
# @param filterGroup string (or vector of strings) of groups from 
#    'filterGroupColumn' to be plottted
# 
# @return A list with two matrices filtered for user defined tissue type
filterGroups = function(plotTables, 
                        metadata, 
                        filterGroupColumn, 
                        filterGroup){
  plotSignalMatrix = plotTables$plotSignalMatrix
  plotBoxStats = plotTables$plotBoxStats
  
  if(length(filterGroupColumn) >1){
    stop("Only one column from metadata can be specified 
         for filtering.")
  } else if (is.null(filterGroupColumn) & !is.null(filterGroup)){
    stop("You must provide the column name from 
         which groups are selected.")
  } else if (!is.null(filterGroupColumn) & is.null(filterGroup)){
    stop(paste0("Please specify which groups from ", filterGroupColumn, 
                "do you want plotted."))
  } else if (is.null(filterGroupColumn) & is.null(filterGroup)){
    plotSignalMatrix = plotSignalMatrix
    plotBoxStats = plotBoxStats
  } else {
    if (length(filterGroup) == 1){
      if (filterGroup %in% levels(factor(metadata[,get(filterGroupColumn)]))){
        plotSignalMatrix = plotSignalMatrix[get(filterGroupColumn) == filterGroup]
        plotBoxStats = plotBoxStats[get(filterGroupColumn) == filterGroup]
      } else {
        stop("The input cell group is not in predefined list of options.")
      }
    } else if (all(filterGroup %in% levels(factor(metadata[,get(filterGroupColumn)])))) {
      plotSignalMatrix = plotSignalMatrix[get(filterGroupColumn) %in% filterGroup]
      plotBoxStats = plotBoxStats[get(filterGroupColumn) %in% filterGroup]
    } else {
      stop("At least one of the input cell groups is not in predefined list of 
           options.")
    }
  }
  
  plotTables = list(plotSignalMatrix = plotSignalMatrix, plotBoxStats = plotBoxStats)
  return(plotTables)
  }

# Internal helper function to plot jitter plot
# 
# @param plotSignalMatrix Result from \code{reshapeDataToPlot} function followed 
#        by \code{filterGroups} function
# @param myLabels Labels created to group corresponding tissue type together
# @param colorColumn Column in 'plotSignalMatrix' used for coloring groups
#
# @return ggplot object - jitter plot
SummarySignalJitterPlot = function(plotSignalMatrix, colorColumn){
  if("mixedVar" %in% names(plotSignalMatrix)){
    myLabels = unique(plotSignalMatrix[, .(mixedVar, colName)])
    myLabels[, spaceLabel := gsub("_", " ", myLabels[, colName])]
    
    if(!is.null(colorColumn)){
      jitterPlot = ggplot(plotSignalMatrix, aes(x=mixedVar, y=signal, 
                                                color=get(colorColumn)))+
        labs(color = colorColumn)
    } else{
      jitterPlot = ggplot(plotSignalMatrix, aes(x=mixedVar, y=signal))
    }
    jitterPlot = jitterPlot + scale_x_discrete(labels=myLabels$spaceLabel)
  } else {
    jitterPlot = ggplot(plotSignalMatrix, aes(x=colName, y=signal))
  }

  
  if ("name" %in% names(plotSignalMatrix)){
    jitterPlot = jitterPlot + facet_grid(name ~ .)
  }
  
  jitterPlot = jitterPlot + 
    geom_jitter(alpha=0.5, height=0, width=0.35) +
    geom_boxplot(outlier.colour=NA, fill=NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          text = element_text(size=10)) +
    xlab("") +
    ylab("signal") + 
    theme_blank_facet_label() +
    theme(strip.text.y.right = element_text(angle = 0))
  return(jitterPlot)
}

# Internal helper function to plot box plot
# 
# @param plotSignalMatrix Result from \code{reshapeDataToPlot} function followed 
#        by \code{filterGroups} function
# @param plotBoxStats Result from \code{reshapeDataToPlot} function followed 
#        by \code{filterGroups} function
# @param colorColumn Column in 'plotSignalMatrix' used for coloring groups
#
# @return ggplot object - box plot
SummarySignalBoxPlot = function(plotSignalMatrix, plotBoxStats, 
                                colorColumn){
  if("mixedVar" %in% names(plotSignalMatrix)){
    myLabels = unique(plotSignalMatrix[, .(mixedVar, colName)])
    myLabels[, spaceLabel := gsub("_", " ", myLabels[, colName])]
    
    if(!is.null(colorColumn)){
      boxPlot = ggplot(plotSignalMatrix, aes(x=mixedVar, y=signal, 
                                                fill=get(colorColumn)))+
        labs(fill = colorColumn)
    } else{
      boxPlot = ggplot(plotSignalMatrix, aes(x = mixedVar, y = signal))
    }
    boxPlot = boxPlot + scale_x_discrete(labels=myLabels$spaceLabel)
  } else {
    boxPlot = ggplot(plotSignalMatrix, aes(x=colName, y=signal))
  }
  
  if ("name" %in% names(plotSignalMatrix)){
    boxPlot = boxPlot + facet_grid(name ~ .)
  }
  boxPlot = boxPlot + geom_boxplot(outlier.colour=NA, alpha=0.9) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          text = element_text(size=10)) +
    xlab("") +
    ylab("normalized signal") + 
    ylim(min(plotBoxStats$value), max(plotBoxStats$value)) +
    theme_blank_facet_label() +
    theme(strip.text.y.right = element_text(angle = 0))
  return(boxPlot)
}

# Internal helper function to plot bar plot
# 
# @param plotBoxStats Result from \code{reshapeDataToPlot} function followed 
#        by \code{filterGroups} function
# @param myLabels Labels created to group corresponding tissue type together
# @param colorColumn Column in 'plotBoxStats' used for coloring groups
#
# @return ggplot object - bar plot

SummarySignalBarPlot = function(plotBoxStats,colorColumn){
  if("mixedVar" %in% names(plotBoxStats)){
    myLabels = unique(plotBoxStats[, .(mixedVar, colName)])
    myLabels[, spaceLabel := gsub("_", " ", myLabels[, colName])]
    
    if(!is.null(colorColumn)){
      barPlot = ggplot(plotBoxStats[plotBoxStats$boxStats == "median",], 
                       aes(x=mixedVar, 
                           y=value, 
                           fill=get(colorColumn))) +
        labs(fill=colorColumn)
    } else{
      barPlot = ggplot(plotBoxStats[plotBoxStats$boxStats == "median",], 
                       aes(x=mixedVar, 
                           y=value))
    }
    barPlot = barPlot + scale_x_discrete(labels=myLabels$spaceLabel)
  } else {
    barPlot = ggplot(plotBoxStats[plotBoxStats$boxStats == "median",], 
                     aes(x=colName, 
                         y=value))
  }
  
  if ("name" %in% names(plotBoxStats)){
    barPlot = barPlot + facet_grid(name ~ .)
  }
  barPlot = barPlot +
    geom_col(alpha=0.9)+
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          text = element_text(size=10)) +
    xlab("") +
    ylab("med(signal)") + 
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
# @param colorColumn Column in 'plotSignalMatrix' used for coloring groups
#
# @return ggplot object - violin plot

SummarySignalViolinPlot = function(plotSignalMatrix,plotBoxStats, 
                                   colorColumn){
  if("mixedVar" %in% names(plotSignalMatrix)){
    myLabels = unique(plotSignalMatrix[, .(mixedVar, colName)])
    myLabels[, spaceLabel := gsub("_", " ", myLabels[, colName])]
    
    if(!is.null(colorColumn)){
      violinPlot = ggplot(plotSignalMatrix, aes(x=mixedVar, y=signal, 
                                                fill=get(colorColumn))) +
        labs(fill = colorColumn)
    } else{
      violinPlot = ggplot(plotSignalMatrix, aes(x=mixedVar, y=signal))
    }
    violinPlot = violinPlot + scale_x_discrete(labels=myLabels$spaceLabel)
  } else {
    violinPlot = ggplot(plotSignalMatrix, aes(x=colName, y=signal))
  }
  
  
  if ("name" %in% names(plotSignalMatrix)){
    violinPlot = violinPlot + facet_grid(name ~ .)
  }
  
  violinPlot = violinPlot + 
    geom_violin(alpha=0.8, scale = "width", trim = TRUE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          text = element_text(size=10)) +
    xlab("") +
    ylab("signal") + 
    theme_blank_facet_label() +
    theme(strip.text.y.right = element_text(angle = 0))
  if("mixedVar" %in% names(plotSignalMatrix)){
    violinPlot = violinPlot +
      geom_point(data = plotBoxStats[plotBoxStats$boxStats == "median",], 
                                        aes(x = mixedVar, y = value), 
                 color = "black", size = 2)
  } else {
    violinPlot = violinPlot +
      geom_point(data = plotBoxStats[plotBoxStats$boxStats == "median",], 
                 aes(x = colName, y = value), 
                 color = "black", size = 2)
  }
  return(violinPlot)
}
