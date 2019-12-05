## Width calculation function
calcWidth = function(query) { 
  width(query)
}

## plotting function
plotQthist = function(df, EndBarColor = "maroon", MiddleBarColor = "gray57", quantile=NULL, bins=NULL) {
  # if user gives no bins or quantiles, we calculate everything based on size of dataset
  if(is.null(quantile) & is.null(bins)){
    q = calcQuantiles(df) 
    div = calcDivisions(df, q)
  }
  # if user gives both bins and quantiles, use their values entirely
  if(!is.null(quantile) & !is.null(bins)){
    q = quantile
    div = calcDivisions(df, q, bins=bins)
  }
  # if user gives quantiles but no bins, we use their quantile to calculate the divisions
  # causes nonsensical results when q is greater than 5
  if(!is.null(quantile) & is.null(bins)){
    q = quantile
    div = calcDivisions(df, q)
  }
  # if user gives bins but no quantiles, we calculate quantiles then insert bin number into the calcDivisions function
  if(is.null(quantile) & !is.null(bins)){
    q = calcQuantiles(df) 
    div = calcDivisions(df, q, bins= bins)
  }
  
  
  colors_vect = c( EndBarColor , rep(MiddleBarColor, (length(div)-3)), EndBarColor) # creates a vector for the colors
  
  df = cutDists(df, divisions= div) # calculating a frequency table with the specified divisions
  if ("name" %in% names(df)){
    # It has multiple regions
    g = ggplot(df, aes(x=cuts, y=Freq, fill=name)) + 
      facet_grid(. ~name)
  } else {
    g = ggplot(df, aes(x=cuts, y=Freq))
  }
  
  g = g +
    geom_bar(stat="identity", fill = colors_vect) + 
    theme_classic() + 
    theme(aspect.ratio=1) + 
    theme_blank_facet_label() +
    ylab("Frequency") +
    xlab("") +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + # vlab()
    theme(plot.title = element_text(hjust = 0.5)) + # Center title
    ggtitle("Quantile Trimmed Histogram") +
    theme(legend.position="bottom") +
    geom_text(aes(label= paste(q,"%", sep='')), data=df[c(1,length(df$Freq)),], vjust=-1)
  return(g)
}


######################################## helper functions

############# calculating the quantiles
calcQuantiles = function(df){
  n = length(df) # finding number of observations
  if (n > 1000) {n = 1000}
  if (n < 100) {n = 100}
  q = (11 - round(n/100)) # finding quantiles that should be used for end boxes
  return(q)
}

############# calulating the divisions
calcDivisions = function(df, q, bins=NULL){
  q = as.numeric(q)
  if(q>=50){
    stop("Quantile can not be larger than 50. Optimal size is under 10.")
  }
  if(!is.null(bins)){
    b = bins +1
  } 
  else {
    b = ((20-(2*q))/q) # FIX THIS!!! finding the number of bins based on the quantiles
  }
  quant = unname(quantile(df, probs = c((q/100), (1-(q/100)))))
  seq_10 = seq(quant[1], quant[2], length = b)
  div = c(-Inf, round(seq_10), Inf)
  return(div)
}

############ LABEL CUTS
# need to put the comma function in labelCuts because it takes a numeric vector and produces character labels
labelCuts = function(breakPoints, digits=1, collapse="-", infBins=FALSE) {
  labels = 
    apply(round(cbind( breakPoints[-length(breakPoints)],	
                       breakPoints[-1]),digits), 1, paste0, collapse=collapse) 
  
  if (infBins) {
    labels[1] = paste0("<", breakPoints[2])
    labels[length(labels)] = paste0(">", breakPoints[length(breakPoints)-1])
  }
  return(labels)
}

############## CUT DISTANCES
### This function calls label cuts
cutDists = function(dists, divisions = c(-Inf, -1e6, -1e4, -1000, -100, 0, 100, 1000, 10000, 1e6, Inf)) {
  if (is.list(dists)) {
    x = lapply(dists, cutDists)
    
    # To accommodate multiple lists, we'll need to introduce a new 'name'
    # column to distinguish them.
    nameList = names(dists)
    if(is.null(nameList)) {
      nameList = 1:length(query) # Fallback to sequential numbers
    }
    
    # Append names
    xb = rbindlist(x)
    xb$name = rep(nameList, sapply(x, nrow))
    
    return(xb)
  }
  
  labels = labelCuts(signif(divisions,3), collapse=" to ", infBins=TRUE)
  cuts = cut(dists, divisions, labels)
  df = as.data.frame(table(cuts))
  return(df)
}


######### CHANGING THEME
theme_blank_facet_label = function() {
  return(theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )
  )
}




