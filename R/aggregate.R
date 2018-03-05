
# Aggregate a BSDT across regions or region groups
# 
# This function is as BScombineByRegion, but can handle not only multiple
# samples in BSDT, but also simultaneously multiple region sets by passing
# a regionsGRL (GRangesList object).
# you can use jExpr to do other functions.
# 
# Given a bisulfite data table as input, with an identifier column for
# different samples; plus a GRanges objects with regions to aggregate.
#
# @param BSDT The bisulfite data.table (output from one of the parsing
# functions for methylation calls) that you wish to aggregate. It can
# be a combined table, with individual samples identified by column passed
# to splitFactor.
# @param regionsGRL Regions across which you want to aggregate.
# @param excludeGR A GenomicRanges object with regions you want to 
# exclude from the aggregation function. These regions will be eliminated
# from the input table and not counted.
# @param jExpr You can pass a custom command in the j slot to data.table
# specifying which columns to aggregate, and which functions to use. You
# can use buildJ() to build a jExpr argument easily.
# @param byRegionGroup You can aggregate by regionID or by regionGroupID; 
# this reflects the regionsGRL that you pass; by default, BSAggregate will
# aggregate each region individually -- scores will then be contiguous, and
# the output is 1 row per region.
# Turn on this flag to aggregate across all region groups, making the result
# uncontiguous, and resulting in 1 row per *region group*.
#
# @export
BSAggregate = function(BSDT, regionsGRL, excludeGR=NULL, regionsGRL.length=NULL,
	splitFactor=NULL, keepCols=NULL, sumCols=NULL, jExpr=NULL,
	byRegionGroup=FALSE, keep.na=FALSE) {

	# Assert that regionsGRL is a GRL.
	# If regionsGRL is given as a GRanges, we convert to GRL
	if( methods::is(regionsGRL,"GRanges")) {
		regionsGRL = GRangesList(regionsGRL)
	} else if (! methods::is(regionsGRL, "GRangesList")) {
		stop("regionsGRL is not a GRanges or GRangesList object")
	}

	if(! is.null(excludeGR)) {
		BSDT = BSFilter(BSDT, minReads=0, excludeGR)
	}

	bsgr = BSdtToGRanges(list(BSDT))

	additionalColNames = setdiff(colnames(BSDT), c("chr","start", "end","hitCount","readCount", splitFactor))

	colModes = sapply(BSDT,mode)
	if (is.null(sumCols)) {
		sumCols = setdiff(colnames(BSDT),c("chr", "start", "end", "strand", splitFactor, keepCols))
		# Restrict to numeric columns.		
		sumCols = intersect(sumCols, names(colModes[which(colModes == "numeric")]))

	}
	# It's required to do a findoverlaps on each region individually,
	# Not on a GRL, because of the way overlaps with GRLs work. So,
	# we must convert the GRL to a GR, but we must keep track of which
	# regions came from which group.
	regionsGR = unlist(regionsGRL)
	
	if(is.null(regionsGRL.length)) {
		if (length(regionsGRL) > 100) {
		message("BSAggregate: Calculating sizes. You can speed this up by supplying a regionsGRL.length vector...", appendLF=FALSE)
		}
		regionsGRL.length = sapply(regionsGRL, length)
		message("Done counting regionsGRL lengths.")
	}

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

	setkey(BSDT, chr, start)
	# Gut check:
	# stopifnot(all(elementMetadata(bsgr[[1]])$readCount == BSDT$readCount))

	message("Setting regionIDs...")
	BSDT = BSDT[queryHits(fo),] #restrict the table to CpGs in any region.

	if (NROW(BSDT) < 1) {
		warning("No BSDT sites in the given region list; please expand your regionsGRL")
		return(NULL)
	}

	BSDT[,regionID:=subjectHits(fo)] #record which region they overlapped.
	#BSDT[queryHits(fo),regionID:=subjectHits(fo)]
	#if (!keep.na) {
	#	BSDT = BSDT[queryHits(fo),]
	#}

	if (is.null(jExpr)) {
		cols=c(sumCols, keepCols)
		funcs = c(rep("sum", length(sumCols)), rep("unique", length(keepCols)))
		jExpr = buildJ(cols, funcs)
	}
	message("jExpr: ", jExpr)
	
	# Define aggregation column. aggregate by region or by region group?
	if (byRegionGroup) {
		agCol = "regionGroupID"
	} else {
		agCol = "regionID" # Default
	}

	# Build the by string
	if (is.null(splitFactor)) {
		byString = paste0("list(regionID)")
	} else {
		byString = paste0("list(", paste("regionID", paste0(splitFactor, ""), collapse=", ", sep=", "), ")")
	}

	# Now actually do the aggregate:
	message("Combining...")
	bsCombined = BSDT[,eval(parse(text=jExpr)), by=eval(parse(text=byString))]
	setkey(bsCombined, regionID)
	# Now aggregate across groups.
	# I do this in 2 steps to avoid assigning regions to groups,
	# which takes awhile. I think this preserve memory and is faster.

	# Define aggregation column. aggregate by region or by region group?
	if (byRegionGroup) {
		# must set allow=TRUE here in case there are multiple IDs (splitCol)
		bsCombined[region2group, regionGroupID:=regionGroupID, allow=TRUE]
		if (! is.null(splitFactor) ) { 
			byStringGroup = paste0("list(", paste("regionGroupID", paste0(splitFactor, collapse=", "), sep=", "), ")")
		} else {
			byStringGroup = "list(regionGroupID)"
		}
		bsCombined=bsCombined[,eval(parse(text=jExpr)), by=eval(parse(text=byStringGroup))]
		return(bsCombined)
	} else {
		e = region2group[bsCombined,]
		setkey(e, regionID)
		return(e)
	}
	# WARNING: There are now 2^2 ways to aggregate, sum vs mean
	# at each level: across regions, then across region sets. THis
	# doesn't give you a choice at this point. 
}
