## Width calculation function
calcWidth = function(query) { 
  width(query)
}

## plotting function
plotQthist = function(df, EndBarColor = "gray57", MiddleBarColor = "gray27", quantile=NULL, bins=NULL) {
  
      output = calcDivisions(df, quantile=quantile, bins=bins)
      bins = output[["bins"]]
      quantile = output[["quantile"]]
      div = output[["divisions"]]
      
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
        geom_text(aes(label= paste(quantile,"%", sep='')), data=df[c(1,length(df$Freq)),], vjust=-1)
      return(g)
}

######################################## helper functions 

############# calulating the divisions for the graph
calcDivisions = function(df, bins=NULL, quantile = NULL){
      # calculating bins
      if(!is.null(bins)){
        b = bins
      } 
      else {
        n = length(df)
        if (n > 10000) {n = 10000}
        if (n < 500) {n = 500}
        b = round(n^.15 + (n/200) ) # finding number of bins based on the size of dataset
      }
      # calculating quantiles
      if(!is.null(quantile)){
        if(quantile >= 20){
          stop("Quantile can not be larger than 20. Optimal size is under 10.")
        }
        q = quantile
      }
      else{
        q = round(25/(b)) # finding the quantile on each size based on number of bins
        q = max(1,q) # minimum percentage on each side is 1%
      }
      quant = unname(quantile(df, probs = c((q/100), (1-(q/100)))))
      seq_10 = seq(quant[1], quant[2], length = b+1)
      div = c(-Inf, round(seq_10), Inf)
      listOutput <- list("bins"= b,"quantile"= q, "divisions" = div)
      return(listOutput)
}

############ LABEL CUTS
# I changed this helper function so that the x-axis would have commas
labelCuts = function(breakPoints, digits=1, collapse="-", infBins=FALSE) {
      roundedLabels = signif(round(cbind( breakPoints[-length(breakPoints)],breakPoints[-1]),digits), 3)
      is.na(roundedLabels) <- sapply(roundedLabels, is.infinite) #set the Inf values to NA so that formatC could add commas
      labelsWithCommas = formatC(roundedLabels, format="d", big.mark=",")
      labels = apply(labelsWithCommas, 1, paste0, collapse=collapse) 
      if (infBins) {
        labels[1] = paste0("<", formatC(breakPoints[2], format="d", big.mark=","))
        labels[length(labels)] = paste0(">", formatC(breakPoints[length(breakPoints)-1], format="d", big.mark=","))
      }
      return(labels)
}
