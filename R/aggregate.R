
# Quick way to count overlaps between a region set and one or more other 
# region sets. 
#
# Count how many regions
# from the first region set (queryRegionDT) overlap with each of the regions
# from the other region set/s.
# Uses only the midpoint of the first region set when finding overlaps. 
#
# @param queryRegionDT data.frame/data.table. Must have "chr" and "start"
# columns.
# @param regionsGRL GRangesList or GRanges. E.g. Binned chromosomes, with
# each chromosome as a GRanges object in the GRangesList.
# @return a data.table with the following columns: 
# regionID,  chr, start,  end, withinGroupID, regionGroupID, N
# the coordinates refer to the regions from regionsGRL.
# "regionGroupID" refers to which GRanges from regionsGRL the given
# region was a member of. "withinGroupID" refers to the index of the given
# region within its GRanges object.
# "regionID" has the index for the given region as if all the GRanges from
# regionsGRL were combined into a single GRanges object
# The "N" column has the counts for number of query regions 
# overlapping with that regionsGRL region
calcOLCount = function(queryRegionDT, regionsGRL) {
    jExpr = ".N"
    queryRegionDT = queryRegionDT
    
    # Assert that regionsGRL is a GRL.
    # If regionsGRL is given as a GRanges, we convert to GRL
    if(methods::is(regionsGRL,"GRanges")) {
        regionsGRL = GRangesList(regionsGRL)
    } else if (! methods::is(regionsGRL, "GRangesList")) {
        stop("regionsGRL is not a GRanges or GRangesList object")
    }
    
    # convert query regions to just the midpoint
    if ("end" %in% colnames(queryRegionDT)) {
        # assign to "start" because BSdtToGRanges only keeps the start coord
        queryRegionDT$start = round((queryRegionDT$start + queryRegionDT$end)/2) 
    }
    
    # only keeps start column
    bsgr = BSdtToGRanges(list(queryRegionDT))

    # It's required to do a findoverlaps on each region individually,
    # Not on a GRL, because of the way overlaps with GRLs work. So,
    # we must convert the GRL to a GR, but we must keep track of which
    # regions came from which group.
    regionsGR = unlist(regionsGRL)
    
    regionsGRL.length = sapply(regionsGRL, length)
    
    # Build a table to keep track of which regions belong to which group
    region2group = data.table(
        regionID=1:length(regionsGR), 
        chr=as.vector(seqnames(regionsGR)), 
        start=as.vector(start(regionsGR)), 
        end=as.vector(end(regionsGR)),
        withinGroupID= as.vector(unlist(sapply(regionsGRL.length, seq))),
        regionGroupID=rep(1:length(regionsGRL), regionsGRL.length))
    setkey(region2group, regionID)
    
    
    message("Finding overlaps...")
    fo = findOverlaps(bsgr[[1]], regionsGR)
    
    setkey(queryRegionDT, chr, start)
    # Gut check:
    # stopifnot(all(elementMetadata(bsgr[[1]])$readCount == queryRegionDT$readCount))
    
    message("Setting regionIDs...")
    queryRegionDT = queryRegionDT[queryHits(fo),] #restrict the table to CpGs in any region.
    
    if (NROW(queryRegionDT) < 1) {
        warning("No overlapping regions in the given region list; please expand your regionsGRL")
        return(NULL)
    }
    
    queryRegionDT[,regionID:=subjectHits(fo)] #record which region they overlapped.
    #queryRegionDT[queryHits(fo),regionID:=subjectHits(fo)]
    #if (!keep.na) {
    #	queryRegionDT = queryRegionDT[queryHits(fo),]
    #}
    
    # Build the by string
    byString = paste0("list(regionID)")
    
    # Now actually do the aggregate:
    message("Combining...")
    bsCombined = queryRegionDT[,eval(parse(text=jExpr)), by=eval(parse(text=byString))]
    setkey(bsCombined, regionID)
    
    e = region2group[bsCombined,]
    setkey(e, regionID)
    return(e)
}
