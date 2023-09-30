#
# Generate input files for Himmeli. The network visualization
# is a streamlined depiction of the module enrichment in hub
# neighborhoods. For instance, only links that connect a key
# driver to another node are depicted.
#
# Input:
#   job               KDA data list as returned by kda.finish()
#
# Optional input:
#   modules           array of module names to be visualized
#   ndrivers          maximum number of drivers per module
#
# Written by Ville-Petteri Makinen 2013
#
kda2himmeli <- function(job, modules=NULL, ndrivers=5) {
  
    # Import node values.
    cat("\nImporting node data...\n")
    valdata <- tool.read(job$nodfile)
    
    # Set node sizes.
    z <- as.double(valdata$VALUE)
    z <- (z/quantile(z, 0.95) + rank(z)/length(z))
    valdata$SIZE <- pmin(4.0, z)
    
    # Select subset of genes.
    valdata <- kda2himmeli.identify(valdata, "NODE", job$graph$nodes)
    print(summary(valdata))
    
    # Select top scoring modules.
    cat("\nForwarding KDA results to Himmeli...\n")
    if(is.null(modules) & (is.null(job$ssearesults) == FALSE)) {
        tmp <- job$ssearesults
        tmp <- tmp[order(tmp$P),]
        modules <- tmp$modules
        if(length(modules) > 8) modules <- modules[1:8]
    }
    
    # Convert module names to indices.
    if(is.null(modules) == FALSE) {
        modules <- match(modules, job$modules)
        modules <- modules[which(modules > 0)]
        if(length(modules) < 1) stop("Unknown module names.")
    }
    
    # Select top key drivers from each module.
    drivers <- kda2himmeli.drivers(job$results, modules, ndrivers)
    mods <- unique(drivers$MODULE)
    palette <- kda2himmeli.colormap(length(mods))
    
    # Process each module separately.
    edgdata <- data.frame()
    noddata <- data.frame()
    for(i in 1:length(mods)) {
        rows <- which(drivers$MODULE == mods[i])
        tmp <- kda2himmeli.exec(job, valdata, drivers[rows,], mods, palette)
        edgdata <- rbind(edgdata, tmp$edat)
        noddata <- rbind(noddata, tmp$vdat)
    }
    
    # Create work folder.
    dpath <- file.path(job$folder, "himmeli")
    if(file.exists(dpath) == FALSE) dir.create(path=dpath, recursive=TRUE)
    
    # Save data files.
    edgfile <- file.path(job$folder, "kda2himmeli.edges.txt")
    nodfile <- file.path(job$folder, "kda2himmeli.nodes.txt")
    tool.save(frame=edgdata, file=edgfile)
    tool.save(frame=noddata, file=nodfile)
    
    # Configuration data.
    gname <- file.path(dpath, job$label)
    instr <- paste("GraphName", gname, sep="\t")
    instr <- c(instr, paste("EdgeFile", edgfile, sep="\t"))
    instr <- c(instr, paste("EdgeTailVariable", "TAIL", sep="\t"))
    instr <- c(instr, paste("EdgeHeadVariable", "HEAD", sep="\t"))
    instr <- c(instr, paste("EdgeWeightVariable", "WEIGHT", sep="\t"))
    instr <- c(instr, paste("EdgeColorVariable", "COLOR", sep="\t"))
    instr <- c(instr, paste("VertexFile", nodfile, sep="\t"))
    instr <- c(instr, paste("VertexNameVariable", "NODE", sep="\t"))
    instr <- c(instr, paste("VertexColorVariable", "COLOR", sep="\t"))
    instr <- c(instr, paste("VertexLabelVariable", "LABEL", sep="\t"))
    instr <- c(instr, paste("VertexShapeVariable", "SHAPE", sep="\t"))
    instr <- c(instr, paste("VertexSectorVariable", "SECTOR", sep="\t"))
    instr <- c(instr, paste("VertexSizeVariable", "SIZE", sep="\t"))
    instr <- c(instr, paste("DistanceUnit", "1.1", sep="\t"))
    instr <- c(instr, paste("ChassisMode", "on", "2.0", sep="\t"))
    
    # Color info.
    modnames <- job$modules[mods]
    for(j in 1:ncol(palette)) {
        c <- palette[,j]
        value <- sprintf("%02d%02d%02d", c[1], c[2], c[3])
        value <- paste("VertexColorInfo", modnames[j], value, sep="\t")
        instr <- c(instr, value)
    }
    
    # Save configuration data.
    fname <- file.path(job$folder, "kda2himmeli.config.txt")
    write.table(x=instr, sep="\t", file=fname, na="",
    row.names=FALSE, quote=FALSE)
    cat("\rSaved ", length(instr), " rows in '", fname, "'.\n", sep="")
    return(job)
}

#----------------------------------------------------------------------------

kda2himmeli.exec <- function(job, valdata, drivers, modpool, palette) {
    
  # Create star topology.
  edgdata <- data.frame()
  for(i in unique(drivers$NODE)) {
      tmp <- kda2himmeli.edges(job$graph, i, job$depth, job$direction)
      edgdata <- rbind(edgdata, tmp)
  }
    
    # Select affected nodes.
    tmp <- c(edgdata$TAIL, edgdata$HEAD)
    pos <- match(valdata$NODE, tmp)
    valdata <- valdata[which(pos > 0),]
    
    # Assign node shapes.
    valdata$SHAPE <- "circle"
    pos <- match(valdata$NODE, drivers$NODE)
    valdata[which(pos > 0),"SHAPE"] <- "star"
    
    # Trace module memberships.
    noddata <- kda2himmeli.colorize(valdata, job$moddata, modpool, palette)
    
    # Trim edge dataset.
    edgdata <- edgdata[which(edgdata$TAIL != edgdata$HEAD),]
    edgdata <- unique(edgdata[,c("TAIL", "HEAD", "WEIGHT")])
    edgdata$COLOR <- "cfcfcf"
    
    # Restore original identities.
    edgdata$TAIL <- job$graph$nodes[edgdata$TAIL]
    edgdata$HEAD <- job$graph$nodes[edgdata$HEAD]
    noddata$NODE <- job$graph$nodes[noddata$NODE]
    noddata$LABEL <- noddata$NODE
    
    # Make identities unique for the current module.
    modtag <- job$modules[drivers[1,"MODULE"]]
    for(i in 1:nrow(noddata)){
      noddata[i,"NODE"] <- paste(noddata[i,"NODE"], modtag, sep="@")
    }
    for(i in 1:nrow(edgdata)) {
        edgdata[i,"TAIL"] <- paste(edgdata[i,"TAIL"], modtag, sep="@")
        edgdata[i,"HEAD"] <- paste(edgdata[i,"HEAD"], modtag, sep="@")
    }
    
    # Return results.
    res <- list(edat=edgdata, vdat=noddata)
    return(res)
}

#----------------------------------------------------------------------------

kda2himmeli.drivers <- function(data, modules, ndriv) {
    nmods <- 8 # optimal number of colors
    
    # Select modules.
    if(is.null(modules) == FALSE) {
        nmods <- length(unique(modules))
        pos <- match(data$MODULE, modules)
        data <- data[which(pos > 0),]
        if(nrow(data) < 1) stop("No data on target modules.")
    }
    
    # Include only significant key drivers.
    rows <- which(data$FDR < 0.05)
    data <- data[rows,]
    if(nrow(data) < 1) stop("No key drivers for target modules.")
    
    # Separate modules.
    st <- tool.aggregate(data$MODULE)
    blocks <- st$blocks
    
    # Collect top drivers.
    peaks <- double()
    scores <- data$P
    for(k in 1:length(blocks)) {
        rows <- blocks[[k]]
        peaks[k] <- min(scores[rows])
    }
    
    # Collect drivers from top modules.
    ind <- integer()
    mask <- order(peaks)
    mask <- mask[1:min(length(mask),nmods)]
    for(k in mask) {
        rows <- blocks[[k]]
        rows <- rows[order(scores[rows])]
        rows <- rows[1:min(length(rows),ndriv)]
        ind <- c(ind, rows)
    }
    return(data[ind,c("MODULE","NODE")])
}

#----------------------------------------------------------------------------

kda2himmeli.edges <- function(graph, center, depth, direction) {
    g <- tool.subgraph.search(graph, center, depth, direction)
    if(direction >= 0) {
        g$HEAD <- g$RANK
        g$TAIL <- center
    } else {
        g$TAIL <- g$RANK
        g$HEAD <- center
    }
    g$WEIGHT <- (g$STRENG)/(1.0 + g$LEVEL)
    return(g)
}

#---------------------------------------------------------------------------

kda2himmeli.colorize <- function(noddata, moddata, modpool, palette) {
    
    # Collect module memberships.
    pos <- match(moddata$NODE, noddata$NODE)
    moddata <- moddata[which(pos > 0),c("MODULE","NODE")]
    moddata <- unique(moddata)
    
    # Merge duplicate rows.
    st <- tool.aggregate(moddata$NODE)
    blocks <- st$blocks
    colors <- rep(NA, length(blocks))
    sectors <- rep("", length(blocks))
    for(k in 1:length(blocks)) {
        rows <- blocks[[k]]
        mods <- moddata[rows,"MODULE"]
        mods <- intersect(mods, modpool)
        mods <- sort(unique(mods))
        n <- length(mods)
        if(n < 1) next
        if(n < 2) {
            c <- palette[,which(modpool == mods[1])]
            colors[k] <- sprintf("%02d%02d%02d", c[1], c[2], c[3])
            next;
        }
        for(i in 1:n) {
            c <- palette[,which(modpool == mods[i])]
            c <- sprintf("%02d%02d%02d", c[1], c[2], c[3])
            if(i < 2) {
                sectors[k] <- paste(sectors[k], "1:", c, sep="")
            } else {
                sectors[k] <- paste(sectors[k], ",1:", c, sep="")
            }
        }
    }
    
    # Combine results.
    res <- data.frame(NODE=as.integer(st$labels), COLOR=colors,
    SECTOR=sectors, stringsAsFactors=FALSE)
    res <- merge(noddata, res, all.x=TRUE)
    
    # Fill in missing values.
    rows <- which(is.na(res$COLOR))
    res[rows,"COLOR"] <- "cfcfcf"
    return(res)
}

#----------------------------------------------------------------------------

kda2himmeli.colormap <- function(ncolors) {
    palette <- rainbow(n=(ncolors + 2))
    palette <- (col2rgb(palette)/255)
    
    # Dampen raw colors.
    kappa <- pmax(0.25*palette[1,], 0.25*palette[2,], 0.4*palette[3,])
    palette[1,] <- ((1.0 - kappa)*palette[1,] + kappa)
    palette[2,] <- ((1.0 - kappa)*palette[2,] + kappa)
    palette <- round(99*palette)
    
    # Remove strongest greens (easy to confuse on screen).
    while(ncol(palette) > ncolors) {
        rb <- (palette[1,] + palette[3,])
        palette <- palette[,order(rb)]
        ind <- which.max(palette[,2])
        mask <- setdiff(1:ncol(palette), ind)
        palette <- palette[,mask]
    }
    
    # Sort by components.
    x <- (10000*palette[1,] + 100*palette[2,] + palette[3,])
    return(palette[,order(x)])
}

#----------------------------------------------------------------------------

kda2himmeli.identify <- function(dat, varname, labels) {
    if(nrow(dat) < 1) return(dat)
    
    # Find matching identities.
    pos <- match(dat[,varname], labels)
    rows <- which(pos > 0)
    
    # Select subset.
    dat[,varname] <- pos
    res <- dat[rows,]
    return(res)
}


#
# Generate input files for cytoscape. The network visualization
# is a streamlined depiction of the module enrichment in hub
# neighborhoods. Hence, only links that connect a key
# driver to another node are depicted within a particular depth.
#
# Input:
#   job               KDA data list as returned by kda.finish()
#
# Optional input:
#   modules           array of module names to be visualized
#   ndrivers          maximum number of drivers per module
#   depth             search depth of subgraphs
#
# Written by Zeyneb Kurt 2015, Modified by Jessica Ding 2021
#
kda2cytoscape <- function(job, node.list=NULL, modules=NULL, ndrivers=5,
depth=1) {
    # Select top scoring modules.
    cat("\nForwarding KDA results to Cytoscape...\n")
    
    if(is.null(modules))  # if module names were not provided by user,
      # take module names from kda results
      modules <- job$modules[unique(job$results$MODULE)]
    
    # Convert module names to indices.
    if(!is.null(modules)) {
        modules <- match(modules, job$modules)
        modules <- modules[which(modules > 0)]
        if(length(modules) < 1) stop("Unknown module names.")
    }
    
    # Select top key drivers from each module.
    drivers <- kda2cytoscape.drivers(job$results, modules, ndrivers)
    mods <- unique(drivers$MODULE)
    allmods <- c(mods, setdiff(unique(job$moddata$MODULE), mods))
    modnames <- job$modules[allmods]
    modnames[which(mods == 0)] <- "NON.MODULE"
    
    if(!is.null(node.list)){
      if(all(is.na(match(node.list, job$graph$nodes)))){
        stop("These nodes are not in the provided graph")
      }
      else{
        drivers <- c()
        nds <- which(!is.na(match(job$graph$nodes, node.list)))
        for(ii in 1:length(nds)){
            if (length(which(job$results$NODE == nds[ii])) > 0 ){
                mdls <- job$results$MODULE[which(job$results$NODE == 
                nds[ii])]
                drivers <- rbind(drivers, cbind(mdls, nds[ii]))
            }
            else drivers <- rbind(drivers, cbind(0, nds[ii]))
        }
      }
    }
    drivers <- as.data.frame(drivers)
    colnames(drivers) <- c("MODULE" , "NODE","FDR")
    
    # Create work folder.
    dpath <- file.path(job$folder, "cytoscape")
    if(file.exists(dpath) == FALSE)
      dir.create(path=dpath, recursive=TRUE)
    # Save top KDAs into file
    # make cols for all modules
    palette <- kda2cytoscape.colormap(length(allmods))
    drivers$MODNAMES <- modnames[match(drivers$MODULE, mods)]
    drivers$NODNAMES <- job$graph$nodes[drivers$NODE]
    for(i in 1:nrow(drivers)){
        k <- which(mods == drivers$MODULE[i])
        drivers$COLOR[i] <- palette[k]
    }
    
    kdafile <- file.path(dpath, paste(job$label, ".kda2cytoscape.top.kds.txt", sep =""))
    tool.save(frame=drivers, file=kdafile)
    
    # Process each module separately.
    edgdata <- data.frame()
    noddata <- data.frame()
    for(i in 1:length(mods)) {
        rows <- which(drivers$MODULE == mods[i])
        if(length(rows) > 0){
            tmp <- kda2cytoscape.exec(job, drivers[rows,], allmods, palette,
                                      depth)
            if (!is.null(tmp)){
                edgdata <- rbind(edgdata, tmp$edat)
                noddata <- rbind(noddata, tmp$vdat)
            }
        }
    }

    keep.list <- c()
    unq.nodes <- unique(noddata$NODE)
    for ( i in 1:length(unq.nodes) ) {
        rws <- which(noddata[,1] == unq.nodes[i])
        if (any(is.na(noddata[rws, "SECTOR"])) ) {
            nn <- which(as.integer(noddata[rws, 3]) == 200)
            if (length(nn) >0) keep.list <- c(keep.list, rws[nn[1]])
            else keep.list <- c(keep.list, rws[1])
        }
        else {
            kd <- which(noddata[rws, "SHAPE"] == "Diamond")
            if (length(kd) >0){
                keep.list <- c(keep.list, rws[kd[1]])
                lenlen <- length(strsplit(noddata[rws[kd[1]], "SECTOR"], "1:")[[1]])
                noddata[rws, "COLOR"] <-
                    strsplit(noddata[rws[kd[1]], "SECTOR"], "1:")[[1]][lenlen]
                noddata[rws, "COLOR"] <- paste0("#", noddata[rws[kd[1]], "COLOR"])
            }
            else{
                keep.list <- c(keep.list, rws[1])
                lenlen <- length(strsplit(noddata[rws[1], "SECTOR"], "1:")[[1]])
                noddata[rws, "COLOR"] <-
                    strsplit(noddata[rws[1], "SECTOR"], "1:")[[1]][lenlen]
                noddata[rws, "COLOR"] <- paste0("#", noddata[rws[1], "COLOR"])
            }
        }
    }
    
    noddata <- noddata[keep.list, ]
    noddata$SECTOR[is.na(noddata$SECTOR)] <- ""
    noddata$URL[is.na(noddata$URL)] <- ""
    
    # Save data files.
    edgfile <- file.path(dpath, paste(job$label, ".kda2cytoscape.edges.txt", sep =""))
    nodfile <- file.path(dpath, paste(job$label, ".kda2cytoscape.nodes.txt", sep =""))
    
    whole.network=cbind(job$graph$tails, job$graph$heads, job$graph$weight)
    TAIL.ids=match(edgdata$TAIL, job$graph$nodes)
    HEAD.ids=match(edgdata$HEAD, job$graph$nodes)
    edge.weights.to.be.found=cbind(TAIL.ids, HEAD.ids)
    direction=job$direction
    if(direction > 0)
      edge.weights.to.be.found=merge(edge.weights.to.be.found, whole.network,
                                     by.x=c(1,2), by.y=c(1,2))
    if(direction < 0)
      edge.weights.to.be.found=merge(edge.weights.to.be.found, whole.network,
                                     by.x=c(1,2), by.y=c(2,1))
    if(direction == 0)
      edge.weights.to.be.found=unique(rbind(merge(edge.weights.to.be.found, 
                                                  whole.network, by.x=c(1,2), 
                                                  by.y=c(1,2)), 
                                            merge(edge.weights.to.be.found, 
                                                  whole.network, by.x=c(1,2), 
                                                  by.y=c(2,1))))
    edge.weights.to.be.found[,1]=job$graph$nodes[edge.weights.to.be.found[,1]]
    edge.weights.to.be.found[,2]=job$graph$nodes[edge.weights.to.be.found[,2]]
    
    edgdata.wout.wgts=edgdata[,-3]
    edge.weights.to.be.found=merge(edgdata.wout.wgts, edge.weights.to.be.found, 
                                   by.x=c(1,2), by.y=c(1,2))
    edge.weights.to.be.found= edge.weights.to.be.found[,c(1,2,5,3,4)]
    
    colnames(edge.weights.to.be.found)[3] <- "WEIGHT"
    if(sum(duplicated(edge.weights.to.be.found[,c("TAIL","HEAD","COLOR","MODULE")]))){
      cat("\nDuplicate edges found!\n Check network file.\n\n")
    }
    if (!is.null(edgdata)) tool.save(frame=edge.weights.to.be.found, file=edgfile)
    
    num <- NULL
    color <- NULL
    mod <- NULL
    
    # Color info.
    for(j in 1:length(modnames)) {
      num[j] <- j
      color[j] <- palette[j]
      mod[j] <- modnames[j]
    }
    
    instr <- data.frame(NUMBER=num, MODULES=mod, COLOR=color)
    
    noddata$MODULE <- vapply(noddata$URL, function(x){
      mods <- c()
      for(col in instr$COLOR){
        if(grepl(gsub("#","",col), x)){
          mods <- c(mods, as.character(instr$MODULES[instr$COLOR==col]))
        }
      }
      if(length(mods)>0){
        return(do.call("paste",c(as.list(mods), list("sep"=","))))
      } else{
        return("")
      }
    }, "character")
    
    kd_res <- job$results
    kd_res$MODULE <- job$modules[kd_res$MODULE]
    kd_res$NODE <- job$graph$nodes[kd_res$NODE]
    kd_res <- kd_res[kd_res$FDR<0.05,]
    
    kd_info = t(vapply(noddata$NODE, function(x){
      info <- c()
      if(x %in% kd_res$NODE){
        # KD_MODULE
        info <- c(info, 
                  do.call("paste",
                          c(as.list(kd_res$MODULE[kd_res$NODE==x]),
                            list("sep"=","))))
        # KD_TOP_MODULE
        kdres <- kd_res[kd_res$NODE==x,]
        top_mod <- kdres$MODULE[which(min(kdres$P)==kdres$P)]
        info <- c(info, 
                  do.call("paste",c(kdres$MODULE[which(min(kdres$P)==kdres$P)],
                                    list("sep"=","))))
        # KD_BORDER_COLOR
        if(length(top_mod)>1) top_mod = top_mod[1]
        info <- c(info, instr$COLOR[instr$MODULES==top_mod])
        
        # KD_subnetwork_shown
        if(x %in% drivers$NODNAMES){
          info <- c(info, "YES")
        } else {
          info <- c(info, "NO")
        }
        
        return(info)
        
      } else{
        return(rep("Not a KD", 4))
      }
    }, FUN.VALUE = character(4)))
    
    colnames(kd_info) <- c("KD_MODULE","KD_TOP_MODULE",
                           "KD_BORDER_COLOR","KD_subnetwork_shown")
    noddata <- cbind(noddata, kd_info)
    
    modmember_nonkd <- which(noddata$URL!="" & noddata$SHAPE!="Diamond")
    noddata$SIZE[modmember_nonkd] <- 100
    kd_sn_not_shown <- which(noddata$KD_subnetwork_shown=="NO")
    noddata$SHAPE[kd_sn_not_shown] <- "Diamond"
    noddata$SIZE[kd_sn_not_shown] <- 100
    noddata$LABELSIZE[kd_sn_not_shown] <- 20
    
    for(col in 1:ncol(noddata)){
      noddata[,col] <- as.character(noddata[,col])
    }
    
    if (!is.null(noddata)){
      tool.save(frame=noddata[, c("NODE", "LABEL", "COLOR", 
                                  "SIZE", "SHAPE", "SECTOR",
                                  "URL", "MODULE","KD_MODULE",
                                  "KD_TOP_MODULE","KD_BORDER_COLOR",
                                  "KD_subnetwork_shown","LABELSIZE")], 
                file=nodfile)
    }
    
    fname <- file.path(dpath, paste(job$label, ".module.color.mapping.txt", sep =""))
    write.table(instr, file=fname, sep="\t", row.names=FALSE, quote=FALSE)
    
    cat("\rInformation is stored into relevant files.\n", sep="")
    return(job)
}

#----------------------------------------------------------------------------

kda2cytoscape.exec <- function(job, drivers, modpool, palette, graph.depth=1){

  # Create star topology.
  edgdata <- data.frame()
  for(j in unique(drivers$NODE)) {
    tmp <- kda2cytoscape.edges(job$graph, j, graph.depth, job$direction)
    edgdata <- rbind(edgdata, tmp)
  }
  
  # Select affected nodes.
  tmp <- data.frame(unique(c(edgdata$TAIL, edgdata$HEAD)),
  stringsAsFactors=FALSE)
  names(tmp) <- "NODE"
  # Assign node shapes.
  shapes <- rep("Ellipse", length(tmp))
  sizes <- rep(50, length(tmp))
  
  tmp$SHAPE <- shapes
  tmp$SIZE <- sizes
  
  pos <- match(tmp$NODE, drivers$NODE)
  tmp[which(pos > 0),"SHAPE"] <- "Diamond"
  tmp[which(pos > 0),"SIZE"] <- 200
  
  # Trace module memberships.
  noddata <- kda2cytoscape.colorize(tmp, job$moddata, modpool, palette)
  
  if(!is.null(noddata)){
    # Trim edge dataset.
    edgdata <- edgdata[which(edgdata$TAIL != edgdata$HEAD),]
    edgdata <- unique(edgdata[,c("TAIL", "HEAD", "WEIGHT")])
    edgdata$COLOR <- "#cfcfcf"
    
    # Restore original identities.
    edgdata$TAIL <- job$graph$nodes[edgdata$TAIL]
    edgdata$HEAD <- job$graph$nodes[edgdata$HEAD]
    modtag <- job$modules[drivers[1,"MODULE"]]
    if(length(modtag) == 0) modtag <- "NON.MODULE"
    
    noddata$NODE <- job$graph$nodes[as.numeric(noddata$NODE)]
    noddata$LABEL <- noddata$NODE
    
    # Make identities unique for the current module.
    
    noddata[1:nrow(noddata), "MODULE"] <- modtag
    
    edgdata[1:nrow(edgdata), "MODULE"] <- modtag
    
    # Return results.
    res <- list(edat=edgdata, vdat=noddata)
    return(res)
  }
  else res=NULL
  return(res)
}

#----------------------------------------------------------------------------

kda2cytoscape.drivers <- function(data, modules, ndriv) {

    nmods <- length(unique(data$MODULE))
    # Select modules.
    if(!is.null(modules)) {
        nmods <- length(unique(modules))
        pos <- match(data$MODULE, modules)
        data <- data[which(pos > 0),]
        if(nrow(data) < 1) stop("No data on target modules.")
    }
    
    # Include only significant key drivers.
    rows <- which(data$FDR < 0.05)
    data <- data[rows,]
    if(nrow(data) < 1) stop("No key drivers for target modules.")
    
    # Separate modules.
    st <- tool.aggregate(data$MODULE)
    blocks <- st$blocks
    
    # Collect top drivers.
    peaks <- double()
    scores <- data$P
    for(k in 1:length(blocks)) {
        rows <- blocks[[k]]
        peaks[k] <- min(scores[rows])
    }
    
    # Collect drivers from modules.
    ind <- integer()
    mask <- order(peaks)
    mask <- mask[1:min(length(mask),nmods)]
    for(k in mask) {
        rows <- blocks[[k]]
        rows <- rows[order(scores[rows])]
        rows <- rows[1:min(length(rows),ndriv)]
        ind <- c(ind, rows)
    }
    
    kds_include <- unique(data[ind,c("MODULE","NODE")]$NODE)
    return(data[data$NODE %in% kds_include,c("MODULE","NODE","FDR")])
}

#----------------------------------------------------------------------------

kda2cytoscape.edges <- function(graph, center, depth, direction) {
    global.HEAD <- c()
    global.TAIL <- c()
    global.WEIGHT <- c()
    
    centeri <- center
    for(i in 1:depth){
        new.cent <- c()
        for(j in centeri){
            gi <- tool.subgraph.search(graph, j, depth=1, direction)
            if(direction >= 0) {
                global.HEAD <- c(global.HEAD, gi$RANK)
                global.TAIL <- c(global.TAIL, rep(j, length(gi$RANK)) )
            } else {
                global.TAIL <- c(global.TAIL, gi$RANK)
                global.HEAD <- c(global.HEAD, rep(j, length(gi$RANK)) )
            }
            new.cent <- unique(c(new.cent, gi$RANK))
            global.WEIGHT <- c(global.WEIGHT,gi$STRENG)
        }
        centeri <- new.cent[which(!(new.cent %in% centeri))]
    }
    g <- unique(data.frame(cbind(HEAD=global.HEAD, TAIL=global.TAIL, 
                                 WEIGHT=global.WEIGHT)))
    
    return(g)
}

#---------------------------------------------------------------------------

kda2cytoscape.colorize <- function(noddata, moddata, modpool, palette) {

    # Google chart service.
    urlbase <- "http://chart.apis.google.com/chart?cht=p&chs=200x200"
    urlbase <- paste(urlbase, "chf=bg,s,00000000", sep="&")
    
    # Collect module memberships.
    pos <- match(moddata$NODE, noddata$NODE)
    moddata <- moddata[which(pos > 0),c("MODULE","NODE")]
    moddata <- unique(moddata)
    
    # Merge duplicate rows.
    if(length(moddata$NODE) > 0){
        st <- tool.aggregate(moddata$NODE)
        blocks <- st$blocks
        colors <- rep(NA, length(blocks))
        sectors <- rep("", length(blocks))
        urls <- rep("", length(blocks))
        
        for(k in 1:length(blocks)) {
            chd <- ""
            chco <- ""
            rows <- blocks[[k]]
            mods <- moddata[rows,"MODULE"]
            mods <- intersect(mods, modpool)
            mods <- sort(unique(mods))
            n <- length(mods)
            if(n < 1) next
            
            for(i in 1:n) {
                c <- palette[which(modpool == mods[i])]
                c <- gsub("#","",c)
                if(i < 2) {
                    chd <- "chd=t:1"
                    chco <- paste("chco=", c, sep="")
                    cc <- paste("#", c, sep="")
                    colors[k] <- cc
                    sectors[k] <- paste(sectors[k], "1:", c, sep="")
                } else {
                    chd <- paste(chd, 1, sep=",")
                    chco <- paste(chco, c, sep="|")
                    sectors[k] <- paste(sectors[k], ",1:", c, sep="")
                }
            }
            urls[k] <- paste(urlbase, chd, chco, sep="&")
        }
        
        # keep 1st clr in colors, keep all sectors in sectors
        res <- data.frame(NODE=as.integer(st$labels), COLOR=colors,
                          SECTOR=sectors, URL=urls, stringsAsFactors=FALSE)
        res <- merge(noddata, res, all.x=TRUE)
        res$LABELSIZE <- 12
        
        # Fill in missing values.
        rows <- which(is.na(res$COLOR))
        res[rows,"COLOR"] <- "#cfcfcf"
        res[which(res$SHAPE=="Ellipse" & is.na(res$URL)), "SIZE"] <- 50
        res[which(as.integer(res$SIZE) == 100), "LABELSIZE"] <- 20
        res[which(res$SHAPE == "Diamond"), "LABELSIZE"] <- 30
        return(res)
    }
    else res <- NULL
    return (res)
}

#----------------------------------------------------------------------------

kda2cytoscape.colormap <- function(ncolors) {

  palette <- rainbow(n=(ncolors + 2))
  palette <- vapply(palette, function(x){
    col <- col2rgb(x)
    col <- col + (255 - col)*0.5
    col <- rgb(t(col), maxColorValue=255)
    return(col)
  }, "character")
  palette <- (col2rgb(palette)/255)
  
  # Dampen raw colors.
  kappa <- pmax(0.25*palette[1,], 0.25*palette[2,], 0.4*palette[3,])
  palette[1,] <- ((1.0 - kappa)*palette[1,] + kappa)
  palette[2,] <- ((1.0 - kappa)*palette[2,] + kappa)
  palette <- round(255*palette)
  
  x <- (10000*palette[1,] + 100*palette[2,] + palette[3,])
  
  cols <- apply(palette, 2, function(x){
    return(rgb(x[1], x[2], x[3], maxColorValue = 255))
  })
  
  return(cols)
}

#----------------------------------------------------------------------------

kda2cytoscape.identify <- function(dat, varname, labels) {
    if(nrow(dat) < 1) return(dat)
    
    # Find matching identities.
    pos <- match(dat[,varname], labels)
    rows <- which(pos > 0)
    
    # Select subset.
    dat[,varname] <- pos
    res <- dat[rows,]
    return(res)
}

# Determine statistical significance of key driver genes.
#
# Input:
#   job$graph            see tool.graph()
#   job$graph$hubs       nodes considered hubs (indexed)
#   job$graph$hubnets    lists of neighboring nodes (indexed)
#   job$graph$cohubsets  lists of overlapping hubs (indexed)
#   job$module2nodes     lists of node indices for each module
#
# Output:
#   job$results          data frame of results (indexed)
#
# Written by Ville-Petteri Makinen 2013, Modified by Jessica Ding 2021
#
kda.analyze <- function(job) {
    set.seed(job$seed)
    
    cat("\nAnalyzing network...\n")
    nmods <- length(job$modules)
    
    # Analyze modules.
    res <- data.frame()
    hubs <- job$graph$hubs
    for(i in 1:nmods) {
        members <- job$module2nodes[[i]]
        p <- kda.analyze.exec(members, job$graph, job$nperm)
        mask <- which(p >= 0.0)
        if(length(mask) < 1) next
        
        # Find top hit.
        tmp <- data.frame(MODULE=i, NODE=hubs[mask], P=p[mask])
        pmin <- min(tmp$P)
        hit <- which(tmp$P == pmin)
        hit <- tmp$NODE[hit[1]]
        
        # Update results.
        nmemb <- length(members)
        name <- job$graph$nodes[hit]
        kd <- sprintf("%s, n=%d, p=%.2e", name, nmemb, pmin)
        cat(job$modules[i], ": ", kd, "\n", sep="")
        res <- rbind(res, tmp)
    }
    
    # Estimate false discovery rates.
    res$FDR <- p.adjust(res$P, method="fdr")
    job$results <- res
    return(job)
}

#----------------------------------------------------------------------------

kda.analyze.exec <- function(memb, graph, nsim) {
    hubs <- graph$hubs
    hubnets <- graph$hubnets
    nhubs <- length(hubs)
    nnodes <- length(graph$nodes)
    nmemb <- length(memb)
    
    # Observed enrichment scores.
    obs <- rep(NA, nhubs)
    for(k in 1:nhubs) {
        g <- hubnets[[hubs[k]]]
        obs[k] <- kda.analyze.test(g$RANK, g$STRENG, memb, nnodes)
    }
    
    # Estimate P-values.
    pvals <- rep(NA, nhubs)
    for(k in which(obs > 0)) {
        g <- hubnets[[hubs[k]]]
        
        # First pass.
        x <- kda.analyze.simulate(obs[k], g, nmemb, nnodes, 200)
        if(sum(is.na(x))==length(x)) next
        
        # Estimate preliminary P-value.
        param <- tool.normalize(x)
        z <- tool.normalize(obs[k], param)
        p <- pnorm(z, lower.tail=FALSE)
        if(p*nhubs > 2.0) next
        
        # Estimate final P-value.
        n <- min(nsim, abs(nsim - length(x))+1) # it was: (nsim - length(x))
        y <- kda.analyze.simulate(obs[k], g, nmemb, nnodes, n)
        param <- tool.normalize(c(x, y))
        z <- tool.normalize(obs[k], param)
        p <- pnorm(z, lower.tail=FALSE)
        p <- max(p, .Machine$double.xmin)
        if(p*nhubs > 1.0) next
        
        # Apply Bonferroni adjustment.
        pvals[k] <- p*nhubs
    }
    return(pvals)
}

#----------------------------------------------------------------------------

kda.analyze.simulate <- function(o, g, nmemb, nnodes, nsim) {
    neigh <- as.integer(g$RANK)
    w <- as.double(g$STRENG)
    
    # Simulate null distribution.
    nfalse <- 0
    x <- rep(NA, nsim)
    
    deviation = FALSE
    tries <- 0
    while(!deviation){
        if(tries>10) return(x) # unlikely to occur
        tries <- tries + 1
        for(n in 1:nsim) {
            if(nfalse > 20) break # was 10
            memb <- sample.int(nnodes, nmemb) 
            x[n] <- kda.analyze.test(neigh, w, memb, nnodes)
            if(is.na(x[n])) x[n] <- rnorm(1)
            nfalse <- (nfalse + as.integer(x[n] >= o))
        }
        # Trim results.
        x <- x[which(0*x == 0)]
        
        # check if values are sufficiently deviated
        if(length(x[x>min(x)])==1){
            deviation = TRUE
        } else if(sd(x[x>min(x)])!=0 & length(unique(x))>4){
            deviation = TRUE
        } else {
            x <- rep(NA, nsim)
            nfalse <- 0
        }
    }
    
    return(x)
}

#----------------------------------------------------------------------------

#kda.analyze.test <- function(ind, w, members, nnodes) {
#  shared <- which(match(ind, members) > 0)
#  obsmass <- sum(w[shared])
#  return(obsmass)
#}

#----------------------------------------------------------------------------

kda.analyze.test <- function(neigh, w, members, nnodes) {
    
    # Check if enough neighbors.
    nneigh <- length(neigh)
    if(nneigh < 3) return(0.0)
    
    # Find member nodes.
    nmemb <- length(members)
    pos <- match(neigh, members)
    shared <- which(pos > 0)
    
    # Background edge mass.
    totmass <- sum(w)
    
    # Edge mass captured by members.
    obsmass <- sum(w[shared])
    
    # Effective number of shared nodes.
    rho <- (obsmass/totmass)
    nobserv <- rho*nneigh
    
    # Expected number of shared nodes.
    nexpect <- (nmemb/nnodes)*nneigh
    
    # Calculate enrichment score.
    if(nobserv < 1.0) return(NA)
    z <- (nobserv - nexpect)/(sqrt(nexpect) + 1.0)
    return(z)
}
#
# Set parameters for a weighted key driver analysis
#
# Input:
#   plan$label        unique identifier for the analysis
#   plan$folder       parent folder for results
#   plan$netfile      path to network file
#                     columns: TAIL HEAD WEIGHT
#   plan$modfile      path to module file
#                     columns: MODULE NODE
#
# Optional:
#   plan$inffile      path to module info file
#                     columns: MODULE DESCR
#   plan$nodfile      path to node selection file
#                     columns: NODE
#   plan$depth        search depth for subgraph search
#   plan$direction    use zero for undirected, negative for
#                     downstream and positive for upstream
#   plan$maxoverlap   maximum allowed overlap between two
#                     key driver neighborhoods
#   plan$minsize      minimum module size
#   plan$mindegree    minimum node degree to qualify as a hub
#   plan$maxdegree    maximum node degree to include
#   plan$edgefactor   influence of node strengths:
#                     0.0 no influence, 1.0 full influence
#   plan$seed         seed for random number generator
#
# Output:
#   job               data structure for KDA
#
# Written by Ville-Petteri Makinen 2013
#
kda.configure <- function(plan) {
  if(is.null(plan$folder)) stop("No parent folder.")
  if(is.null(plan$label)) stop("No job label.")
  jdir <- file.path(plan$folder, "kda")
  
  if(!dir.exists(jdir)) dir.create(path=jdir, recursive=TRUE)
  if(file.access(jdir, 2) != 0) stop("Cannot access '" + plan$folder + "'.")
  
  cat("\nKDA Version:12.01.2021\n")
  cat("\nParameters:\n")
  plan$stamp <- Sys.time()
  if(is.null(plan$netfile)) stop("No network file.")
  if(is.null(plan$modfile)) stop("No module file.")
  
  if(is.null(plan$depth)) plan$depth <- 1
  plan$depth <- round(plan$depth)
  if(plan$depth < 1) stop("Unusable search depth.")
  cat("  Search depth: ", plan$depth, "\n", sep="")
  
  if(is.null(plan$direction)) plan$direction <- 0
  plan$direction <- round(plan$direction)
  if(abs(plan$direction) > 1) stop("Unusable search direction.")
  cat("  Search direction: ", plan$direction, "\n", sep="")
  
  if(is.null(plan$maxoverlap)) plan$maxoverlap <- 0.33
  if(plan$maxoverlap > 1) stop("Unusable overlap limit.")
  if(plan$maxoverlap < 0) stop("Unusable overlap limit.")
  cat("  Maximum overlap: ", plan$maxoverlap, "\n", sep="")
  
  if(is.null(plan$minsize)) plan$minsize <- 20
  if(plan$minsize < 1) stop("Unusable size limit.")
  cat("  Minimum module size: ", plan$minsize, "\n", sep="")
  
  if(is.null(plan$mindegree)) plan$mindegree <- "automatic"
  if(is.null(plan$maxdegree)) plan$maxdegree <- "automatic"
  cat("  Minimum degree: ", plan$mindegree, "\n", sep="")
  cat("  Maximum degree: ", plan$maxdegree, "\n", sep="")
  
  if(is.null(plan$edgefactor)) plan$edgefactor <- 0.5
  if(plan$edgefactor > 1) stop("Unusable edge factor.")
  if(plan$edgefactor < 0) stop("Unusable edge factor.")
  cat("  Edge factor: ", plan$edgefactor, "\n", sep="")
  
  if(is.null(plan$nperm)) plan$nperm <- 2000
  
  if(is.null(plan$seed)) plan$seed <- 1
  cat("  Random seed: ", plan$seed, "\n", sep="")
  return(plan)
}
# Organize and save results.
#
# Input:
#   job$label      unique identifier for the analysis
#   job$folder     output folder for results
#
# Results are also saved in tab-delimited text files.
#
# Written by Ville-Petteri Makinen 2013
#
kda.finish <- function(job) {
  cat("\nFinishing results...\n")
  if (nrow(job$results)==0){
    cat("No Key Driver Found!!!!")
  } else {
    
    # Estimate additional measures.
    res <- kda.finish.estimate(job)
    
    # Save full results.
    res <- kda.finish.save(res, job)
    
    # Create a simpler file for viewing.
    res <- kda.finish.trim(res, job)
    
    # Create a summary file of top hits.
    res <- kda.finish.summarize(res, job)
    
    # Create a inputs and parameters file
    job <- kda.finish.param(job)
    
    return(job)
  }
}

#---------------------------------------------------------------------------

kda.finish.estimate <- function(job) {
    res <- job$results
    
    # Collect module sizes.
    sizes <- rep(0, nrow(res))
    modnames <- res$MODULE
    mod2nod <- job$module2nodes
    for(i in 1:nrow(res)) {
        key <- modnames[i]
        sizes[i] <- length(mod2nod[[key]])
    }
    res$N.mod <- sizes
    
    # Collect overlaps with hub neighborhoods.
    nnodes <- length(job$graph$nodes)
    hubnets <- job$graph$hubnets
    sizes <- matrix(nrow=nrow(res), ncol=4)
    nodenames <- res$NODE
    for(i in 1:nrow(res)) {
        key <- modnames[i]
        node <- nodenames[i]
        memb <- mod2nod[[key]]
        g <- hubnets[[node]]
        sizes[i,1] <- nrow(g)
        sizes[i,2] <- length(intersect(g$RANK, memb))
        sizes[i,3] <- nrow(g)*length(memb)/nnodes
        sizes[i,4] <- sum(node == memb)
    }
    
    # Update data frame.
    res$N.neigh <- sizes[,1]
    res$N.obsrv <- sizes[,2]
    res$N.expct <- sizes[,3]
    res$MEMBER <- (sizes[,4] > 0)
    res$FILL <- (res$N.obsrv)/(res$N.neigh + 1e-20)
    res$FOLD <- (res$N.obsrv)/(res$N.expct + 1e-20)
    return(res)
}

#---------------------------------------------------------------------------

kda.finish.param <- function(job) {
  inputs <- data.frame("Input or parameter" = c("Network file",
                                                "Marker set file",
                                                "Search depth",
                                                "Search direction",
                                                "Maximum overlap",
                                                "Minimum module size",
                                                "Minimum degree",
                                                "Maximum degree",
                                                "Edge factor",
                                                "Number of permutations",
                                                "Random seed"),
                       "Value" = c(basename(job$netfile),
                                   basename(job$modfile),
                                   job$depth,
                                   job$direction,
                                   job$maxoverlap,
                                   job$minsize,
                                   job$mindegree,
                                   job$maxdegree,
                                   job$edgefactor,
                                   job$nperm,
                                   job$seed), 
                       check.names = FALSE)
  
  jdir <- file.path(job$folder, "kda")
  fname <- paste(job$label, ".param.txt", sep="")
  tool.save(frame=inputs, file=fname, directory=jdir)
  
  return(job)
}

#---------------------------------------------------------------------------

kda.finish.save <- function(res, job) {
    
    # Collect co-hubs.
    mtx <- matrix(nrow=0, ncol=0)
    nodes <- job$graph$nodes
    masters <- unique(res$NODE)
    cohubsets <- job$graph$cohubsets
    for(key in masters) {
        cohubs <- cohubsets[[key]]
        cohubs <- unique(c(key, cohubs))
        tmp <- matrix(key, nrow=length(cohubs), ncol=2)
        tmp[,2] <- cohubs
        if(ncol(mtx) < 2) mtx <- tmp
        else mtx <- rbind(mtx, tmp)
    }
    
    # Convert indices to identities.
    dat <- data.frame(HUB=mtx[,1], NODE=mtx[,2])
    dat$NODE <- nodes[dat$NODE]
    dat$HUB <- nodes[dat$HUB]
    
    # Save co-hub information.
    jdir <- file.path(job$folder, "kda")
    fname <- paste(job$label, ".hubs.txt", sep="")
    tool.save(frame=dat, file=fname, directory=jdir)
    
    # Merge with module info.
    if(nrow(job$modinfo) > 0) res <- merge(res, job$modinfo, all.x=TRUE)
    
    # Convert indices to identities.
    res$MODULE <- job$modules[res$MODULE]
    res$NODE <- job$graph$nodes[res$NODE]
    
    # Add description back for merged modules 
    # (will be DESCR for reported name)
    if(!is.null(job$inffile)){
      if(!is.null(res$DESCR)){
        modinfo <- read.delim(job$inffile, stringsAsFactors = FALSE)
        modinfo <- modinfo[,c("MODULE","DESCR")]
        modinfo <- rbind(modinfo, 
                         data.frame("MODULE"=paste0(modinfo$MODULE,",.."),
                                    "DESCR"=modinfo$DESCR, 
                                    stringsAsFactors = FALSE))
        for(iter in 1:nrow(res)){
          if(is.na(res$DESCR[iter]) & grepl(",..",res$MODULE[iter])){
            if(gsub(",..","",res$MODULE[iter]) %in% modinfo$MODULE){
              res$DESCR[iter] <- 
                modinfo$DESCR[modinfo$MODULE==gsub(",..","",res$MODULE[iter])]
            }
          }
        }
        res$DESCR[is.na(res$DESCR)] <- ""
      }
    }
    
    # Sort and save results.
    res <- res[order(res$P),]
    fname <- paste(job$label, ".results.txt", sep="")
    tool.save(frame=res, file=fname, directory=jdir)
    return(res)
}

#---------------------------------------------------------------------------

kda.finish.trim <- function(res, job) {
    
    # Select columns.
    header <- c("MODULE", "NODE", "P", "FDR", "FOLD")
    if(is.null(res$DESCR) == FALSE) header <- c(header, "DESCR")
    res <- res[,header]
    
    # Make numbers nicer to look at.
    preals <- res$P; pvals <- rep("", nrow(res))
    fdreals <- res$FDR; fdrates <- rep("", nrow(res))
    ldreals <- res$FOLD; folds <- rep("", nrow(res))
    for(i in 1:nrow(res)) {
      pvals[i] <- sprintf("%.2e", preals[ i])
      if(is.na(fdreals[i])) 
        fdrates[i] <- ""
      else
        fdrates[i] <- sprintf("%.4f", fdreals[i])
      folds[i] <- sprintf("%.2f", ldreals[i])
    }
    
    # Update results.
    res$P <- pvals
    res$FDR <- fdrates
    res$FOLD <- folds
    
    # Rename columns for post-processing.
    trimres <- res
    header[[3]] <- paste("P.", job$label, sep="")
    header[[4]] <- paste("FDR.", job$label, sep="")
    names(trimres) <- header
    
    # Save P-values.
    jdir <- file.path(job$folder, "kda")
    fname <- paste(job$label, ".pvalues.txt", sep="")
    tool.save(frame=trimres, file=fname, directory=jdir)
    return(res)
}

#---------------------------------------------------------------------------

kda.finish.summarize <- function(res, job) {
    
    # Determine ranking scores.
    nres <- nrow(res)
    rA <- rank(as.double(res$P))
    rB <- (nres - rank(as.double(res$FOLD)))
    scores <- (rA*nrow(res) + rB)
    
    # Determine blocks of modules.
    struct <- tool.aggregate(res$MODULE)
    blocks <- struct$blocks
    
    # Find the top node for each block.
    tops <- rep(0, length(blocks))
    for(k in 1:length(blocks)) {
        rows <- blocks[[k]]
        tmp <- scores[rows]
        pos <- which(tmp == min(tmp))
        tops[k] <- rows[pos]
    }
    
    # Select top drivers.
    scores <- scores[tops]
    res <- res[tops,]
    
    # Save P-values.
    res <- res[order(scores),]
    jdir <- file.path(job$folder, "kda")
    fname <- paste(job$label, ".tophits.txt", sep="")
    tool.save(frame=res, file=fname, directory=jdir)
    return(res)
}
#
# Prepare graph topology for weighted key driver analysis.
#
# Input:
#   job$graph         see tool.graph()
#   job$depth         search depth for subgraph search
#   job$direction     use zero for undirected, negative for
#                     downstream and positive for upstream
#   job$maxoverlap    maximum allowed overlap between two
#                     key driver neighborhoods
#   job$mindegree     minimum hub degree to include
#   job$edgefactor    influence of node strengths:
#                     0.0 no influence, 1.0 full influence
#
# Output:
#   job$graph$hubs       nodes considered hubs (indexed)
#   job$graph$hubnets    lists of neighboring nodes (indexed)
#   job$graph$cohubsets  lists of overlapping hubs (indexed)
#
# Written by Ville-Petteri Makinen 2013, Modified by Le Shu 2015
#
kda.prepare <- function(job) {
    
    # Determine minimum hub degree.
    nnodes <- length(job$graph$nodes)
    #if(job$mindegree == "automatic") {
    #  dmin <- nnodes/median(job$modulesizes)
    #  job$mindegree <- round(0.1*dmin)
    #}
    if (job$mindegree == "automatic") {
        dmin <- as.numeric(quantile(job$graph$stats$DEGREE,0.75))
        job$mindegree <- dmin
        cat("\nMinimum degree set to", dmin,"\n")
    }
    if (job$maxdegree == "automatic") {
        dmax <- as.numeric(quantile(job$graph$stats$DEGREE,1))
        job$maxdegree <- dmax
        cat("\nMaximum degree set to", dmax,"\n")
    }
    # Collect neighbors.
    cat("\nCollecting hubs...\n")
    job$graph <- kda.prepare.screen(job$graph, job$depth, job$direction,
    job$edgefactor, job$mindegree, job$maxdegree)
    if(length(job$graph$hubs) < 1) stop("No usable hubs detected.")
    
    # Collect overlapping co-hubs.
    job$graph <- kda.prepare.overlap(job$graph, job$direction,
    job$maxoverlap)
    
    # Print report.
    nhubs <- length(job$graph$hubs)
    nmem <- (object.size(job$graph))*(0.5^20)
    cat(sprintf("%d hubs (%.2f%%)\n", nhubs, 100*nhubs/nnodes))
    cat("Graph: ", nmem, " Mb\n", sep="")
    return(job)
}

#---------------------------------------------------------------------------

kda.prepare.screen <- function(graph, depth, direction, efactor, dmin, dmax) {
    stamp <- Sys.time()
    hubnets <- list()
    accepted <- integer()
    nnodes <- length(graph$nodes)
    
    # Determine strength cutoff.
    stren <- rep(0, nnodes)
    if(direction <= 0) stren <- (stren + graph$outstats$STRENG)
    if(direction >= 0) stren <- (stren + graph$instats$STRENG)
    slimit <- quantile(stren, 0.75)
    
    degr <- graph$stats$DEGREE
    
    # Select hubs.
    accepted <- integer()
    for(i in which(degr >= dmin & degr <= dmax)) {
        g <- tool.subgraph.search(graph, i, depth, direction)
        
        # Apply edge factor.
        g$STRENG <- (g$STRENG)^efactor
        
        # Use average strength for hub itself (by definition, the hub
        # typically has huge strength within its neighborhood).
        mask <- which(g$LEVEL < 1)
        g[mask,"STRENG"] <- median(g$STRENG)
        
        # Progress report.
        delta <- as.double(Sys.time() - stamp)
        if((delta >= 10.0) & (i > 1)) {
            cat(sprintf("%d hubs (%d nodes)\n", length(accepted), i))
            stamp <- Sys.time()
        }
        
        # Exclude extreme neighborhoods.
        #if(nrow(g) > nnodes/3) next
        #if(nrow(g) < dmin) next
        
        # Store subnetwork.
        hubnets[[i]] <- g[,c("RANK", "STRENG")]
        accepted <- c(accepted, i)
    }
    
    # Return results.
    graph$hubs <- accepted
    graph$hubnets <- hubnets
    return(graph)
}

#---------------------------------------------------------------------------

kda.prepare.overlap <- function(graph, direction, rmax) {
    hubs <- graph$hubs
    nhubs <- length(hubs)
    hubnets <- graph$hubnets
    stamp <- Sys.time()
    
    # Determine node strengths.
    nnodes <- length(graph$nodes)
    stren <- rep(0.0, nnodes)
    if(direction <= 0) stren <- (stren + graph$outstats$STRENG)
    if(direction >= 0) stren <- (stren + graph$instats$STRENG)
    
    # Sort hubs according to strength.
    mask <- order(stren[hubs], decreasing=TRUE)
    hubs <- hubs[mask]
    
    # Collect overlapping co-hubs.
    cohubsets <- list()
    for(i in 1:nhubs) {
        key <- hubs[i]
        
        # Progress report.
        delta <- as.double(Sys.time() - stamp)
        if((delta >= 10.0) & (i > 1)) {
            cat(sprintf("%d/%d co-hub sets\n", i, nhubs))
            stamp <- Sys.time()
        }
        
        # Neighborhood topology.
        g <- hubnets[[key]]
        neighbors <- g$RANK
        locals <- intersect(neighbors, hubs)
        locals <- setdiff(locals, key)
        strenA <- g$STRENG
        
        # Calculate overlaps.
        cohubs <- integer()
        overlaps <- double()
        for(k in locals) {
            w <- hubnets[[k]]
            strenB <- w$STRENG
            
            # Find overlapping nodes.
            posA <- match(neighbors, w$RANK)
            posB <- match(w$RANK, neighbors)
            sharedA <- which(posA > 0)
            sharedB <- which(posB > 0)
            uniqA <- which(is.na(posA))
            uniqB <- which(is.na(posB))
            
            # Calculate strength sums.
            wsharedA <- sum(strenA[sharedA])
            wsharedB <- sum(strenB[sharedB])
            wuniqA <- sum(strenA[uniqA])
            wuniqB <- sum(strenB[uniqB])
            
            # Average symmetric sum.
            wshared <- 0.5*(wsharedA + wsharedB)
            
            # Overlap ratio.
            r <- wshared/(wuniqA + wuniqB + wshared)
            if(r < rmax) next
            
            # Update data structures.
            cohubs <- c(cohubs, k)
        }
        
        # Store results.
        mask <- order(stren[cohubs], decreasing=TRUE)
        cohubsets[[key]] <- cohubs[mask]
    }
    
    # Return results.
    graph$cohubsets <- cohubsets
    return(graph)
}
#
# Import data for weighted key driver analysis.
#
# Input:
#   job               data structure for wKDA
#
# Output:
#   job$graph         indexed topology, see tool.graph()
#   job$modules       module identities
#   job$modinfo       module descriptions (indexed)
#   job$moddata       module data (indexed)
#                     columns: MODULE NODE
#   job$module2nodes  lists of node indices for each module
#   job$modulesizes   module sizes
#
# Written by Ville-Petteri Makinen 2013, Modified by Le Shu 2015
#
kda.start <- function(job) {
    
    # Import topology.
    edgdata <- kda.start.edges(job)
    moddata <- kda.start.modules(job, edgdata)
    
    # Import module descriptions.
    modinfo <- tool.read(job$inffile, c("MODULE", "DESCR"))
    modinfo <- unique(modinfo)
    if(nrow(modinfo) > 0) print(summary(modinfo))
    
    # Create an indexed graph structure.
    job$graph <- tool.graph(edgdata)
    nmem <- (object.size(job$graph))*(0.5^20)
    cat("Graph: ", nmem, " Mb\n", sep="")
    remove(edgdata)
    gc(FALSE)
    
    # Convert identities to indices.
    modules <- unique(moddata$MODULE)
    modinfo <- kda.start.identify(modinfo, "MODULE", modules)
    moddata <- kda.start.identify(moddata, "MODULE", modules)
    moddata <- kda.start.identify(moddata, "NODE", job$graph$nodes)
    
    # Collect module members.
    st <- tool.aggregate(moddata$MODULE)
    blocks <- st$blocks
    members <- as.integer(moddata$NODE)
    for(k in 1:length(blocks)) {
        mask <- blocks[[k]]
        blocks[[k]] <- members[mask]
    }
    
    # Finish results.
    job$modules <- modules
    job$modinfo <- modinfo
    job$moddata <- moddata
    job$module2nodes <- blocks
    job$modulesizes <- st$lengths
    return(job)
}

#----------------------------------------------------------------------------

kda.start.edges <- function(job) {
    
    # Import edge data.
    cat("\nImporting edges...\n")
    varnames <- c("TAIL", "HEAD", "WEIGHT")
    edgdata <- tool.read(file=job$netfile, vars=varnames)
    edgdata$WEIGHT <- suppressWarnings(as.double(edgdata$WEIGHT))
    edgdata <- edgdata[which(edgdata$WEIGHT > 0),]
    
    # Collect node names.
    nodes <- character()
    if(job$direction >= 0) nodes <- c(nodes, edgdata$TAIL)
    if(job$direction <= 0) nodes <- c(nodes, edgdata$HEAD)
    
    # Select nodes.
    if(is.null(job$nodfile) == FALSE) {
        cat("Selecting nodes...\n")
        symdata <- tool.read(file=job$nodfile, vars=c("NODE"))
        nodes <- intersect(as.character(symdata), nodes)
    }
    
    # Calculate node degrees.
    struct <- tool.aggregate(nodes)
    degrees <- struct$lengths
    nodes <- struct$labels
    
    # Automatic maximum degree.
    #if(job$maxdegree == "automatic") {
    #  dmax <- floor(0.05*length(nodes))
    #  mask <- which(degrees <= dmax)
    #  degrees <- degrees[mask]
    #  nodes <- nodes[mask]
    #  job$maxdegree <- dmax
    #}
    
    # Filter edges.
    posT <- match(edgdata$TAIL, nodes)
    posH <- match(edgdata$HEAD, nodes)
    rows <- which(posT*posH > 0)
    edgdata <- edgdata[rows,]
    
    # Print report.
    print(summary(edgdata))
    return(edgdata)
}

#----------------------------------------------------------------------------

kda.start.modules <- function(job, edgdata) {
    
    # Import module data.
    cat("\nImporting modules...\n")
    moddata <- tool.read(file=job$modfile, vars=c("MODULE", "NODE"))
    moddata <- unique(moddata)
    
    # Collect node names.
    nodes <- character()
    if(job$direction >= 0) nodes <- c(nodes, edgdata$TAIL)
    if(job$direction <= 0) nodes <- c(nodes, edgdata$HEAD)
    
    # Filter module members.
    pos <- match(moddata$NODE, nodes)
    moddata <- moddata[which(pos > 0),]
    
    # Remove small modules.
    st <- tool.aggregate(moddata$MODULE)
    mask <- which(st$lengths >= job$minsize)
    mods <- as.character(st$labels[mask])
    pos <- match(moddata$MODULE, mods)
    moddata <- moddata[which(pos > 0),]
    
    # Print report.
    print(summary(moddata))
    return(moddata)
}

#----------------------------------------------------------------------------

kda.start.identify <- function(dat, varname, labels) {
    if(nrow(dat) < 1) return(dat)
    
    # Find matching identities.
    pos <- match(dat[,varname], labels)
    rows <- which(pos > 0)
    
    # Select subset.
    dat[,varname] <- pos
    res <- dat[rows,]
    return(res)
}
#
# Generate inputs for weighted key driver analysis.
#
# Input:
#   job               MSEA data list as returned by ssea.finish().
#
# Optional input:
#   symbols           data.frame for translating gene symbols
#                     columns: FROM TO
#   rmax              maximum allowed overlap between gene sets
#
# Output:
#   plan$label        unique identifier for the analysis
#   plan$parent       parent folder for results
#   plan$modfile      path to module file
#                     columns: MODULE NODE
#   plan$inffile      path to module info file
#                     columns: MODULE DESCR
#   plan$nodfile      path to node selection file
#                     columns: NODE
#
# Written by Ville-Petteri Makinen 2013, Modified by Jessica Ding 2021
#
ssea2kda <- function(job, symbols=NULL, rmax=NULL, min.module.count=NULL, filter=NULL) {
    cat("\nForwarding MSEA results to KDA...\n")
    if(is.null(rmax)) rmax <- 0.33
    if(is.null(filter)) filter <- 0.25
    
    # Collect genes and top markers from original files.
    noddata <- ssea2kda.import(job$genfile, job$marfile)
    
    # Select candidate modules.
    res <- job$results
    res <- res[order(res$P),]
    rows <- which(res$FDR < filter)
    if(!is.null(min.module.count))
      if(length(rows) < min.module.count) rows <- (1:min.module.count)
    if(length(rows) < nrow(res)) res <- res[rows,]
    if(length(rows) < 3) stop("No modules to forward to MSEA with given filter!\n")
    
    # Collect member genes.
    moddata <- job$moddata
    pos <- match(moddata$MODULE, res$MODULE)
    moddata <- moddata[which(pos > 0),]
    
    # Restore original identities.
    modinfo <- job$modinfo
    modinfo$MODULE <- job$modules[modinfo$MODULE]
    moddata$MODULE <- job$modules[moddata$MODULE]
    moddata$GENE <- job$genes[moddata$GENE]
    
    moddata <- moddata[!grepl("_ctrl", moddata$MODULE),]
    moddata_orig <- moddata
    
    # Merge and trim overlapping modules.
    moddata$OVERLAP <- moddata$MODULE
    if(length(unique(moddata$MODULE))>1){
      moddata <- tool.coalesce(items=moddata$GENE, groups=moddata$MODULE,
      rcutoff=rmax)
      moddata$MODULE <- moddata$CLUSTER
      moddata$GENE <- moddata$ITEM
      moddata$OVERLAP <- moddata$GROUPS
      moddata <- moddata[,c("MODULE", "GENE", "OVERLAP")]
      moddata <- unique(moddata)
      
      # add back 'lost' genes
      if(length(setdiff(moddata_orig$GENE, moddata$GENE))>0){
        addBack <- data.frame("GENE"=unique(setdiff(moddata_orig$GENE, moddata$GENE)))
        modules <- c()
        for(gene in addBack$GENE){
          modules <- c(modules, 
                       do.call("paste",
                               c(as.list(moddata_orig$MODULE[moddata_orig$GENE==gene]), 
                                 list("sep"=", "))))
        }
        addBack$Module <- modules
        merged_modules <- unique(moddata$OVERLAP)
        names(merged_modules) <- unique(moddata$MODULE)
        for(iter in 1:nrow(addBack)){
          # get module name from merged
          first_match = unlist(strsplit(addBack$Module[iter], split = ", "))[1]
          mod <- names(merged_modules)[grepl(first_match, merged_modules)]
          if(length(mod)>1){
            for(set in mod){
              submods = unlist(strsplit(merged_modules[[set]], split = ","))
              if(first_match %in% submods){
                mod <- set
                break
              }
            }
          }
          temp <- data.frame("MODULE"=mod, 
                             "GENE"=addBack$GENE[iter], 
                             "OVERLAP"=merged_modules[names(merged_modules)==mod],
                             stringsAsFactors = FALSE)
          moddata <- rbind(moddata, temp)
        }
        moddata <- moddata[order(moddata$MODULE),]
      }
    }
    
    # Calculate enrichment scores for merged modules.
    tmp <- unique(moddata[,c("MODULE","OVERLAP")])
    if(length(unique(moddata$MODULE))>2){
        res <- ssea2kda.analyze(job, moddata)
    } else{
        cat("\nOnly one or two modules passing filter.\n")
        cat("Skipping enrichment score calculation for merged modules.\n\n")
        res <- job$results
        res$MODULE <- job$modules[res$MODULE]
        res <- res[res$MODULE %in% unique(moddata$MODULE),1:6]
    }
    res <- merge(res, tmp)
    res <- merge(res, modinfo, all.x=TRUE)
    
    # Mark modules with overlaps.
    for(i in which(moddata$MODULE != moddata$OVERLAP))
      moddata[i,"MODULE"] <- paste(moddata[i,"MODULE"], "..", sep=",")
    for(i in which(res$MODULE != res$OVERLAP))
      res[i,"MODULE"] <- paste(res[i,"MODULE"], "..", sep=",")
    
    # Separate merged genes.
    nodes <- character()
    genenames <- job$genes
    for(i in 1:length(genenames)) {
        segm <- strsplit(genenames[i], ",", fixed=TRUE)
        nodes <- c(nodes, segm[[1]])
    }
    
    # Expand rows with merged genes.
    moddata_a <- moddata
    tmp <- data.frame(stringsAsFactors=FALSE)
    for(i in 1:nrow(moddata)) {
        segm <- strsplit(moddata[i,"GENE"], ",", fixed=TRUE)
        segm <- segm[[1]]
        if(length(segm) < 2) next
        batch <- data.frame(MODULE=moddata[i,"MODULE"], GENE=segm,
        stringsAsFactors=FALSE)
        tmp <- rbind(tmp, batch)
        moddata[i,] <- NA
    }
    tmp$OVERLAP <- moddata_a$OVERLAP[match(tmp$MODULE, moddata_a$MODULE)]
    moddata <- na.omit(moddata[,c("MODULE","GENE","OVERLAP")])
    moddata <- unique(rbind(moddata, tmp))
    
    # Translate gene symbols.
    moddata$NODE <- moddata$GENE
    noddata$NODE <- noddata$GENE
    if(is.null(symbols) == FALSE) {
        moddata$NODE <- tool.translate(words=moddata$NODE, from=symbols$FROM,
        to=symbols$TO)
        noddata$NODE <- tool.translate(words=noddata$NODE, from=symbols$FROM,
        to=symbols$TO)
        moddata <- na.omit(moddata)
        noddata <- na.omit(noddata)
    }
    
    # Save module info for KDA.
    res <- res[order(res$P),]
    names(res)[5] = "NMARKER"
    inffile <- "msea2kda.info.txt"
    tool.save(frame=res, file=inffile, directory=job$folder)  
    
    # Save modules for KDA.
    modfile <- "msea2kda.modules.txt"
    tool.save(frame=unique(moddata[,c("MODULE", "OVERLAP","NODE", "GENE")]),
    file=modfile, directory=job$folder)  
    
    # Save nodes for KDA.
    nodfile <- "msea2kda.nodes.txt"
    out_node = unique(noddata[,c("NODE", "GENE", "MARKER", "VALUE")])
    names(out_node)[3] = "MARKER"
    tool.save(frame=out_node, file=nodfile, directory=job$folder)  
    
    # Return KDA plan template.
    plan <- list()
    plan$label <- job$label
    plan$folder <- job$folder
    plan$inffile <- file.path(job$folder, inffile)
    plan$modfile <- file.path(job$folder, modfile)
    plan$nodfile <- file.path(job$folder, nodfile)
    plan$ssearesults <- res
    return(plan)
}

#---------------------------------------------------------------------------

ssea2kda.import <- function(genfile, marfile) {
    
    # Import marker values.
    cat("\nImporting marker values...\n")
    locdata <- tool.read(marfile, c("MARKER", "VALUE"))
    locdata$VALUE <- as.double(locdata$VALUE)
    rows <- which(0*(locdata$VALUE) == 0)
    locdata <- unique(na.omit(locdata[rows,]))
    print(summary(locdata))
    
    # Import mapping data.
    cat("\nImporting mapping data...\n")
    gendata <- tool.read(genfile, c("GENE", "MARKER"))
    gendata <- unique(na.omit(gendata))
    print(summary(gendata))
    
    # Merge datasets.
    data <- merge(gendata, locdata)
    
    # Find top marker.
    mask <- integer()
    st <- tool.aggregate(data$GENE)
    blocks <- st$blocks
    for(k in 1:length(blocks)) {
        rows <- blocks[[k]]
        ind <- which.max(data[rows,"VALUE"])
        mask <- c(mask, rows[ind])
    }
    return(data[mask,])
}

#---------------------------------------------------------------------------

ssea2kda.analyze <- function(job, moddata) {
    
    # Convert identities to indices.
    moddata <- ssea.start.identify(moddata, "MODULE", job$modules)
    moddata <- ssea.start.identify(moddata, "GENE", job$genes)
    
    # Collect row indices for each module.
    st <- tool.aggregate(moddata$MODULE)
    keys <- as.integer(st$labels)
    blocks <- st$blocks
    nmods <- length(blocks)
    
    # Prepare gene lists.
    genlists <- list()
    for(k in 1:nmods) genlists[[k]] <- integer()
    
    # Collect gene sets.
    modsizes <- rep(0, nmods)
    genes <- as.integer(moddata$GENE)
    for(k in 1:length(blocks)) {
      key <- keys[k]
      rows <- blocks[[k]]
      members <- unique(genes[rows])
      genlists[[key]] <- as.integer(members)
      modsizes[[key]] <- length(members)
    }
    
    # Determine marker set sizes.
    modlengths <- rep(0, nmods)
    moddensities <- rep(0.0, nmods)
    loclists <- job$database$gene2loci
    for(k in 1:length(genlists)) {
      locset <- integer()
      for(i in genlists[[k]])
        locset <- c(locset, loclists[[i]])
      modlengths[[k]] <- length(unique(locset))
      moddensities[[k]] <- modlengths[[k]]/length(locset)
    }
    
    # Update database.
    job$database$modulesizes <- modsizes
    job$database$modulelengths <- modlengths
    job$database$moduledensities <- moddensities
    job$database$module2genes <- genlists
    
    # Run enrichment analysis.
    job <- ssea.analyze(job)
    res <- job$results
    
    # Restore module identities.
    res$NGENES <- modsizes[res$MODULE]
    res$NMARKER <- modlengths[res$MODULE]
    res$DENSITY <- moddensities[res$MODULE]
    res$MODULE <- job$modules[res$MODULE]
    return(res)
}
#
# Marker set enrichment analysis with gene-level permutations.
#
# Input:
#   job$seed       seed for random number generator
#   job$permtype   random permutation algorithm
#   job$nperm      maximum nubmer of permutations
#   job$database   see ssea.prepare()
#
# Output:
#   job$results    data frame:
#                  MODULE    module identity (indexed)
#                  P         enrichment P-values
#                  FREQ      enrichment P-values (raw frequencies)
#
# Written by Ville-Petteri Makinen 2013, Modified by Le Shu 2016, Modified by Jessica Ding 2021
#
ssea.analyze <- function(job) {
    cat("\nEstimating enrichment...\n")  
    set.seed(job$seed)
    
    # Observed enrichment scores.
    db <- job$database
    scores_int <- ssea.analyze.observe(db)
    scores <- scores_int[["scores"]]
    chis_sd <- scores_int[["chis_sd"]]
    chis_se <- scores_int[["chis_se"]]
    chis_variance <- scores_int[["chis_var"]]
    nmods <- length(scores)
    
    # Simulated scores.
    nperm <- job$nperm
    trim_start = job$trim
    trim_end = 1 - job$trim
    nullsets <- ssea.analyze.simulate(db, scores, nperm, job$permtype, trim_start, trim_end)
    
    # Estimate scores based on Gaussian distribution.
    cat("\nNormalizing scores...\n")
    znull <- double()
    zscores <- NA*scores
    for(i in which(0*scores == 0)) {
        x <- nullsets[[i]]    
        x <- x[which(0*x == 0)]
        param <- tool.normalize(x)
        z <- tool.normalize(x, param)
        zscores[i] <- tool.normalize(scores[i], param)
        znull <- c(znull, z)
    }
    
    # Estimate hit frequencies.
    freq <- NA*scores
    nnull <- length(znull)
    for(i in which(0*scores == 0))
        freq[i] <- sum(zscores[i] <= znull)/nnull
    
    # Estimate scores based on frequencies.
    hscores <- NA*freq
    rows <- which(freq > 5.0/nnull)
    hscores[rows] <- qnorm(freq[rows], lower.tail=FALSE)
    
    # Fill in scores for low frequencies.
    if(length(rows) < length(freq)) {
        omega <- which.max(zscores)
        hscores[omega] <- zscores[omega]
        rows <- c(rows, omega)
        pt <- approx(x=zscores[rows], y=hscores[rows], xout=zscores)
        hscores <- pt$y
    }
    
    # Estimate statistical significance.
    z <- 0.5*(zscores + hscores)
    pvalues <- pnorm(z, lower.tail=FALSE)
    
    # Collect results.
    res <- data.frame(MODULE=(1:nmods), stringsAsFactors=FALSE)
    res$P <- pvalues
    res$FREQ <- freq
    res$ZSCORE <- z
    res$CHI <- scores
    res$CHI_SD <- chis_sd
    res$CHI_SE <- chis_se
    
    job$zscores <- z
    job$var <- chis_variance
    
    # Remove missing scores.
    targets <- which(0*scores == 0)
    job$results <- res[targets,]
    return(job)
}

#----------------------------------------------------------------------------

ssea.analyze.simulate <- function(db, observ, nperm, permtype, trim_start, trim_end) {
  
  #############################################################################
  #####This is an additional process to trim genes with exceptionally high value####
  ###################################################################################
  gene2loci <- db$gene2loci
  locus2row <- db$locus2row
  observed <- db$observed
  #Calcuate individual gene enrichment score
  trim_scores <- rep(NA, length(gene2loci))
  
  for(k in 1:length(trim_scores)) {
    genes <- k
    
    # Collect markers.
    loci <- integer()
    for(i in genes) 
      loci <- c(loci, gene2loci[[i]])
    
    # Determine data rows.
    loci <- unique(loci)
    rows <- locus2row[loci]
    NMARKER <- length(rows)
    
    # Calculate total counts.
    e <- (NMARKER/length(locus2row))*colSums(observed)
    o <- observed[rows,]
    if(NMARKER > 1) o <- colSums(o)
    
    # Estimate enrichment.
    int <- ssea.analyze.statistic(o, e)
    trim_scores[k] <- int[["mean_z"]]
  }
  cutoff=as.numeric(quantile(trim_scores,probs=c(trim_start,trim_end)))
  gene_sel=which(trim_scores>cutoff[1]&trim_scores<cutoff[2])
  
  # Include only non-empty modules for simulation.
  nmods <- length(db$modulesizes)
  targets <- which(db$modulesizes > 0)
  hits <- rep(NA, nmods)
  hits[targets] <- 0
  
  # Prepare data structures to hold null samples.
  keys <- rep(0, nperm)
  scores <- rep(NA, nperm)
  scoresets <- list()
  for(i in 1:nmods)
    scoresets[[i]] <- double()
  
  # Simulate random scores.
  nelem <- 0
  snull <- double()
  stamp <- Sys.time()
  for(k in 1:nperm) {
    if(permtype == "gene") snull <- ssea.analyze.randgenes(db, targets, gene_sel)
    if(permtype == "marker") snull <- ssea.analyze.randloci(db, targets)
    if(length(snull) < 1) stop("Unknown permutation type.")
    
    # Check data capacity.
    ntarg <- length(targets)
    if((nelem + ntarg) >= length(keys)) {
      keys <- c(keys, rep(0, nelem))
      scores <- c(scores, rep(NA, nelem))
    }
    
    # Collect scores.
    for(i in 1:ntarg) {
      nelem <- (nelem + 1)
      t <- targets[i]
      keys[nelem] <- t
      scores[nelem] <- snull[i]
      hits[t] <- (hits[t] + (snull[i] > observ[t]))
    }
    
    # Drop less significant modules.
    hmax <- min(sqrt(nperm/k), 10)
    targets <- which(hits < hmax)
    if(length(targets) < 1) break
    
    # Progress report.
    delta <- as.double(Sys.time() - stamp)
    if((delta < 10.0) & (k < nperm)) next
    cat(sprintf("%d/%d cycles\n", k, nperm))
    stamp <- Sys.time()
  }
  
  # Remove missing entries.
  scores <- scores[1:nelem]
  keys <- keys[1:nelem]
  
  # Organize null scores into lists.
  st <- tool.aggregate(keys)
  labels <- as.integer(st$labels)
  blocks <- st$blocks
  for(i in 1:length(blocks)) {
    key <- labels[i]
    rows <- blocks[[i]]
    scoresets[[key]] <- scores[rows]
  }
  return(scoresets)
}

#----------------------------------------------------------------------------

ssea.analyze.observe <- function(db) {
  mod2gen <- db$module2genes
  gene2loci <- db$gene2loci
  locus2row <- db$locus2row
  observed <- db$observed
  expected <- db$expected
  nmods <- length(mod2gen)  
  
  # Test every module.
  scores <- rep(NA, nmods)
  chis_sd <- rep(NA, nmods)
  chis_se <- rep(NA, nmods)
  chis_variance <- rep(NA, nmods)
  for(k in 1:nmods) {
    genes <- mod2gen[[k]]
    
    # Collect markers.
    loci <- integer()
    for(i in genes) 
      loci <- c(loci, gene2loci[[i]])
    
    # Determine data rows.
    loci <- unique(loci)
    rows <- locus2row[loci]
    NMARKER <- length(rows)    
    
    # Calculate total counts.
    #e <- NMARKER*expected
    e <- (NMARKER/length(locus2row))*colSums(observed)
    o <- observed[rows,]
    if(NMARKER > 1) o <- colSums(o)
    
    # Estimate enrichment.
    # scores[k] <- ssea.analyze.statistic(o, e)
    int <- ssea.analyze.statistic(o, e)
    scores[k] <- int[["mean_z"]]
    chis_sd[k] <- int[["chis_sd"]]
    chis_se[k] <- int[["chis_se"]]
    if("matrix" %in% class(int[["chis_var"]])){
      chis_variance[k] <- NA
    } else{
      chis_variance[k] <- int[["chis_var"]]
    }
  }
  return(list("scores"=scores,
              "chis_sd"=chis_sd,
              "chis_se"=chis_se,
              "chis_var"=chis_variance))
}

#----------------------------------------------------------------------------

ssea.analyze.randgenes <- function(db, targets, gene_sel) {
  mod2gen <- db$module2genes
  modsizes <- db$modulesizes
  modlengths <- db$modulelengths
  gene2loci <- db$gene2loci
  locus2row <- db$locus2row
  observed <- db$observed
  expected <- db$expected
  
  # Test target modules.
  scores <- double()
  nrows <- length(locus2row)
  #npool <- length(gene2loci)
  for(k in targets) {
    msize <- modsizes[[k]]
    NMARKER <- modlengths[[k]]
    
    # Collect pre-defined number of markers from random genes.
    loci <- integer()
    #genes <- sample.int(npool, (msize + 10))
    genes <- sample(gene_sel, (msize + 10))
    while(length(loci) < NMARKER) {
      for(i in genes) {
        tmp <- gene2loci[[i]]
        loci <- c(loci, tmp)
      }
      loci <- unique(loci)
      genes <- sample(gene_sel, (msize + 10))
    }
    
    # Determine data rows.
    loci <- loci[1:NMARKER]
    rows <- locus2row[loci]
    
    # Calculate total counts.
    #e <- NMARKER*expected
    e <- (NMARKER/length(locus2row))*colSums(observed)
    o <- observed[rows,]
    if(NMARKER > 1) o <- colSums(o)
    
    # Estimate enrichment.
    # z <- ssea.analyze.statistic(o, e)
    int <- ssea.analyze.statistic(o, e)
    z <- int[["mean_z"]]
    scores <- c(scores, z)
  }
  return(scores)
}

#----------------------------------------------------------------------------

ssea.analyze.randloci <- function(db, targets) {
  modlengths <- db$modulelengths
  locus2row <- db$locus2row
  observed <- db$observed
  expected <- db$expected
  
  # Test target modules.
  scores <- double()
  nrows <- length(locus2row)
  for(k in targets) {
    NMARKER <- modlengths[[k]]
    
    # Determine data rows.
    loci <- sample.int(nrows, NMARKER)
    rows <- locus2row[loci]
    
    # Calculate total counts.
    #e <- NMARKER*expected
    e <- (NMARKER/length(locus2row))*colSums(observed)
    o <- observed[rows,]
    if(NMARKER > 1) o <- colSums(o)
    
    # Estimate enrichment.
    z <- ssea.analyze.statistic(o, e)
    scores <- c(scores, z[["mean_z"]])
  }
  
  # Return results.
  return(scores)
}

#----------------------------------------------------------------------------

ssea.analyze.statistic <- function(o, e) {
    z <- (o - e)/(sqrt(e) + 1.0)
    
    chis <- z
    chis_sd <- sd(chis)
    chis_se <- sd(chis)/sqrt(length(chis))
    chis_variance <- var(chis) # sd squared
    
    return(list("mean_z"=mean(z),
                "chis_sd"=chis_sd,
                "chis_se"=chis_se,
                "chis_var"=chis_variance))
}
#
# Add internal positive control modules.
#
# Input:
#   job$modules      module identities as characters
#   job$genes        gene identities as characters
#   job$moddata      preprocessed module data (indexed identities)
#   job$database     see ssea.prepare()
#
# Output:
#   job$modules      augmented module names
#   job$moddata      augmented module data
#   job$database     augmented database
#
# Written by Ville-Petteri Makinen 2013
#
ssea.control <- function(job) {
  cat("\nAdding positive controls...\n")
  db <- job$database
  gene2loci <- db$gene2loci
  locus2row <- db$locus2row
  observed <- db$observed
  expected <- db$expected
  
  # Find slots for control module.
  modules <- job$modules
  slotA <- which(modules == "_ctrlA")
  slotB <- which(modules == "_ctrlB")
  if(length(slotA) != 1) stop("No control slot A.")
  if(length(slotB) != 1) stop("No control slot B.")
  
  # Calculate gene scores.
  ngens <- length(gene2loci)
  scores <- rep(NA, ngens)
  for(k in 1:ngens) {
    loci <- gene2loci[[k]]
    NMARKER <- length(loci)
    if(NMARKER < 1) next
    
    # Calculate total counts.
    rows <- locus2row[loci]
    e <- NMARKER*expected
    o <- observed[rows,]
    if(NMARKER > 1) o <- colSums(o)
    
    # Estimate enrichment score.
    #scores[k] <- ssea.analyze.statistic(o, e)
    int <- ssea.analyze.statistic(o, e)
    scores[k] <- int[["mean_z"]]
  }
  
  # Select top genes.
  sizes <- db$modulesizes
  sizes <- sizes[which(sizes > 0)]
  ntop <- floor(median(sizes))
  genesA <- order(scores, decreasing=TRUE)
  genesA <- genesA[1:ntop]
  
  # Collect genes within modules.
  members <- integer()
  mod2gen <- db$module2genes
  for(k in 1:length(mod2gen))
    members <- c(members, mod2gen[[k]])
  members <- unique(members)
  
  # Select top genes among module members.
  genesB <- order(scores[members], decreasing=TRUE)
  genesB <- members[genesB[1:ntop]]
  
  # Collect markers.
  locsetA <- integer()
  locsetB <- integer()
  for(k in genesA)
    locsetA <- c(locsetA, gene2loci[[k]])
  for(k in genesB)
    locsetB <- c(locsetB, gene2loci[[k]])
  locsetA <- unique(locsetA)
  locsetB <- unique(locsetB)
  
  # Force matching number of markers.
  modlen <- min(length(locsetA), length(locsetB))
  
  # Create new modules.
  db$genescores <- scores
  db$modulesizes[[slotA]] <- ntop
  db$modulesizes[[slotB]] <- ntop
  db$modulelengths[[slotA]] <- modlen
  db$modulelengths[[slotB]] <- modlen
  db$moduledensities[[slotA]] <- length(locsetA)/sum(db$genesizes[genesA])
  db$moduledensities[[slotB]] <- length(locsetB)/sum(db$genesizes[genesB])
  db$module2genes[[slotA]] <- genesA
  db$module2genes[[slotB]] <- genesB
  
  # Update module data.
  tmpA <- data.frame(MODULE=slotA, GENE=genesA)
  tmpB <- data.frame(MODULE=slotB, GENE=genesB)
  job$moddata <- rbind(job$moddata, tmpA, tmpB)
  
  # Return results.
  job$database <- db
  remove(db)
  gc(FALSE)
  nmem <- (object.size(job))*(0.5^20)
  cat("Job: ", nmem, " Mb\n", sep="")
  return(job)
}
# Organize and save results.
#
# Input:
#   job$label      unique identifier for the analysis
#   job$folder     output folder for results
#   job$results    data frame:
#                  MODULE    module identity (indexed)
#                  P         enrichment P-values
#   job$database   see ssea.prepare()
#
# Output:
#   job$results       updated columns in data frame:
#                     NGENES   number of distinct member genes
#                     NMARKER    number of distinct member markers
#                     DENSITY  ratio of distinct to non-distinct markers
#                     FDR      false discovery rates
#   job$generesults   data frame:
#                     GENE     gene identity (indexed)
#                     NMARKER    gene size
#                     SCORE    unadjusted enrichment score
#                     MARKER    marker with maximum value
#                     VALUE    marker value
#
# Results are also saved in tab-delimited text files.
#
# Written by Ville-Petteri Makinen 2013
#
ssea.finish <- function(job) {
    cat("\nPostprocessing results...\n")
    job <- ssea.finish.fdr(job)
    job <- ssea.finish.genes(job)
    job <- ssea.finish.details(job)
    job <- ssea.finish.param(job)
    return(job)
}

#----------------------------------------------------------------------------

ssea.finish.genes <- function(job) {
    db <- job$database
    gen2loci <- db$gene2loci
    loc2row <- db$locus2row
    values <- job$locdata$VALUE
    
    # Collect top loci within genes.
    toploci <- integer()
    topvals <- double()
    ngenes <- length(gen2loci)
    for(k in 1:ngenes) {
        locset <- gen2loci[[k]]
        locvals <- values[loc2row[locset]]
        ind <- which.max(locvals)
        toploci[k] <- locset[ind]
        topvals[k] <- locvals[ind]
    }
    
    # Create data frame.
    res <- data.frame(GENE=(1:ngenes))
    res$SCORE <- db$genescores
    res$NMARKER <- db$genesizes
    res$MARKER <- toploci
    res$VALUE <- topvals
    job$generesults <- res
    
    # Restore original identities.
    res$GENE <- job$genes[res$GENE]
    res$MARKER <- job$loci[res$MARKER]
    
    # Save results.
    jdir <- file.path(job$folder, "msea")
    fname <- paste(job$label, ".genes.txt", sep="")
    names(res)[3:4] = c("NMARKER","MARKER")
    tool.save(frame=res, file=fname, directory=jdir)
    return(job)
}

#----------------------------------------------------------------------------

ssea.finish.fdr <- function(job, jobs=NULL) {
  res <- job$results
  
  # Add module statistics.
  db <- job$database
  if(is.null(jobs)){
    res$NGENES <- db$modulesizes[res$MODULE]
    res$NMARKER <- db$modulelengths[res$MODULE]
    res$DENSITY <- db$moduledensities[res$MODULE]
  } else{
    res$NGENES <- job$db_df$modulesizes[match(res$MODULE, job$db_df$MODULE)]
    res$NMARKER <- job$db_df$modulelengths[match(res$MODULE, job$db_df$MODULE)]
    res$DENSITY <- job$db_df$moduledensities[match(res$MODULE, job$db_df$MODULE)]
  }
  
  # Estimate false discovery rates.
  # Find number of modules that made it into the analysis.
  nmods <- length(unique(job$moddata$MODULE)) - 2 # minus 2 for ctrls
  res$FDR <- tool.fdr(res$P, nmods=nmods)
  # Merge with module info.
  res <- merge(res, job$modinfo, all.x=TRUE)
  # Sort according to significance.
  res <- res[order(res$P),]
  job$results <- res
  
  # Restore module names.
  if(is.null(jobs)){
    res$MODULE <- job$modules[res$MODULE]
  } else{
    if(nrow(job$modinfo)>2){
      job$modinfo$MODULE <- job$modules[job$modinfo$MODULE]
      res$DESCR <- job$modinfo$DESCR[match(res$MODULE, job$modinfo$MODULE)]
    } 
  }
  
  # Prepare results for post-processing.
  header <- rep("MODULE", 3)
  header[[2]] <- paste("P.", job$label, sep="")
  header[[3]] <- paste("FDR.", job$label, sep="")
  if(!is.null(jobs)){
    res <- res[,c("MODULE", "P", "FDR", "Cochran.Q","Cochran.P","I2","DESCR")]
    names(res) <- c(header,"Cochran.Q","Cochran.P","I2","DESCR")
  } else{
    res <- res[,c("MODULE", "P", "FDR","DESCR")]
    names(res) <- c(header,"DESCR")
  }
  
  # Make numbers nicer to look at.
  pvals <- character()
  fdrates <- character()
  qstat <- character()
  cpvals <- character()
  for(i in 1:nrow(res)) {
    pvals[i] <- sprintf("%.2e", res[i,2])
    fdrates[i] <- sprintf("%.4f", res[i,3])
    if(!is.null(jobs)){
      qstat[i] <- sprintf("%.2e", res[i,4])
      cpvals[i] <- sprintf("%.2e", res[i,5])
    }
  }
  res[,2] <- pvals
  res[,3] <- fdrates
  if(!is.null(jobs)){
    res[,4] <- qstat
    res[,5] <- cpvals
  }
  
  # Save P-values.
  if(!is.null(jobs)){
    jdir <- file.path(job$folder, "meta")
  } else {
    jdir <- file.path(job$folder, "msea")
  }
  fname <- paste(job$label, ".pvalues.txt", sep="")
  tool.save(frame=res, file=fname, directory=jdir)
  
  # Print individual significance results for 
  # different studies in meta
  if(!is.null(jobs)){
    combine_res <- data.frame("MODULE"=res$MODULE)
    for(i in 1:length(jobs)){
      p <- NULL
      fdr <- NULL
      ires <- jobs[[i]]$results
      # bring back module identities
      ires$MODULE <- jobs[[i]]$modules[ires$MODULE]
      for(iter in 1:length(res$MODULE)){
        if(res$MODULE[iter] %in% ires$MODULE){
          p[iter] <- ires$P[ires$MODULE==res$MODULE[iter]]
          fdr[iter] <- ires$FDR[ires$MODULE==res$MODULE[iter]]
        }
        else{
          p[iter] <- 1
          fdr[iter] <- 1
        }
      }
      if(paste0(jobs[[i]]$label,".P") %in% colnames(combine_res)){
        cat("Warning! Duplicate label names in individual studies.\n")
        cat("Please choose unique labels for each study.\n")
      }
      combine_res[,paste0(jobs[[i]]$label,".P")] = p
      combine_res[,paste0(jobs[[i]]$label,".FDR")] = fdr
      ires_stat <- data.frame("MODULE"=jobs[[i]]$modules,
                              "Zscore"=jobs[[i]]$zscores,
                              "Chi_Var"=jobs[[i]]$var)
      combine_res[,paste0(jobs[[i]]$label,".Z")] = 
        ires_stat$Zscore[match(combine_res$MODULE,
                               ires_stat$MODULE)]
      combine_res[,paste0(jobs[[i]]$label,".Chi_Var")] = 
        ires_stat$Chi_Var[match(combine_res$MODULE,
                                ires_stat$MODULE)]
    }
    
    combine_res$META.P = res[,paste("P.", job$label, sep="")]
    combine_res$META.FDR = res[,paste("FDR.", job$label, sep="")]
    combine_res$Cochran.Q = res[,"Cochran.Q"]
    combine_res$Cochran.P = res[,"Cochran.P"]
    combine_res$I2 = res[,"I2"]
    if("DESCR" %in% colnames(res)) combine_res$DESCR = res$DESCR
    if(sum(combine_res$DESCR!="")==0 | sum(!is.na(combine_res$DESCR))==0){
      combine_res$DESCR = NULL
    }
    for(i in 1:length(jobs)){
      for(col in c("TOPGENES","TOPMARKERS","TOPVALUES")){
        combine_res[,paste0(jobs[[i]]$label,".",col)] = 
          jobs[[i]]$topgenes[,col][match(combine_res$MODULE,
                                         jobs[[i]]$topgenes$MODULE)]
        na_val <- is.na(combine_res[,paste0(jobs[[i]]$label,".",col)])
        combine_res[,paste0(jobs[[i]]$label,".",col)][na_val] <- ""
      }
    }
    fname = paste(job$label, ".combined.results.txt", sep="")
    job$combined_results <- combine_res
    tool.save(frame=combine_res, file=fname, directory=jdir)
  }
  return(job)
}

#----------------------------------------------------------------------------

ssea.finish.details <- function(job, jobs=NULL) {
  
  # Find signficant modules.
  res <- job$results
  mask <- which(res$FDR < 1)
  if(length(mask) < 5) {
    mask <- order(res$P)
    mask <- mask[1:min(5,length(mask))]
  }
  
  # Collect gene members of top modules.
  dtl <- data.frame()
  mod2genes <- job$database$module2genes
  if(!is.null(jobs)){
    names(mod2genes) <- job$modules
  }
  for(k in res[mask,"MODULE"]) {
    genset <- mod2genes[[k]]
    tmp <- data.frame(MODULE=k, GENE=genset)
    dtl <- rbind(dtl, tmp)
  }
  
  # Merge with gene results.
  dtl <- merge(dtl, job$generesults, all.x=TRUE)
  
  # Merge with module info.
  if(nrow(job$modinfo) > 0)
    dtl <- merge(dtl, job$modinfo, all.x=TRUE)
  else
    dtl$DESCR = ""
  
  # Merge with module statistics.
  dtl <- merge(res[,c("MODULE", "P", "FDR")], dtl, all.y=TRUE)
  
  # Sort according to enrichment and marker value.
  scores <- 1000*rank(-(dtl$P))
  gscores <- tool.unify(dtl$VALUE)  
  rows <- order((scores + gscores), decreasing=TRUE)
  dtl <- dtl[rows,]
  
  # Restore names and sort columns.
  if(is.null(jobs)){
    dtl$MODULE <- job$modules[dtl$MODULE]
  }
  dtl$GENE <- job$genes[dtl$GENE]
  dtl$MARKER <- job$loci[dtl$MARKER]
  dtl <- dtl[,c("MODULE", "FDR", "GENE", "NMARKER", "MARKER", "VALUE", "DESCR")]
  
  # subset dtl to only module members
  moddata <- read.delim(job$modfile)
  mod_order <- unique(dtl$MODULE)
  comma <- dtl[grepl(",",dtl$GENE) & !grepl("_ctrlA|_ctrlB",dtl$MODULE),]
  dtl <- dtl[!(grepl(",",dtl$GENE) & !grepl("_ctrlA|_ctrlB",dtl$MODULE)),]
  
  addBack <- do.call("rbind", c(apply(comma, 1, function(x){
    genes <- unlist(strsplit(x["GENE"], split = ","))
    genes <- genes[genes %in% moddata$GENE[moddata$MODULE==x["MODULE"]]]
    return(data.frame("MODULE"=x["MODULE"],
                      "FDR"=as.numeric(x["FDR"]),
                      "GENE"=genes,
                      "NMARKER"=as.numeric(x["NMARKER"]),
                      "MARKER"=x["MARKER"],
                      "VALUE"=as.numeric(x["VALUE"]),
                      "DESCR"=x["DESCR"]))
  })))
  
  dtl <- rbind(dtl, addBack)
  dtl <- dtl[order(dtl$VALUE, decreasing = TRUE),]
  dtl <- dtl[order(match(dtl$MODULE, mod_order)),]
  
  # Make numbers look nicer.
  values <- character()
  fdrates <- character()
  for(i in 1:nrow(dtl)) {
    values[i] <- sprintf("%.2f", dtl[i,"VALUE"])
    fdrates[i] <- sprintf("%.2f%%", 100*dtl[i,"FDR"])
  }
  dtl$FDR <- fdrates
  dtl$VALUE <- values
  
  # Save contents.
  jdir <- file.path(job$folder, "msea")
  fname <- paste(job$label, ".details.txt", sep="")
  names(dtl) = c("MODULE","FDR","GENE","NMARKER","MARKER","VALUE","DESCR")
  tool.save(frame=dtl, file=fname, directory=jdir)
  
  # Add top five genes, markers, values details to result file
  # Restore module names.
  res$MODULE <- job$modules[res$MODULE]
  res$TOPGENES <- ""
  res$TOPMARKERS <- ""
  res$TOPVALUES <- ""
  for(i in 1:length(res$MODULE)){
    if(res$MODULE[i] %in% dtl$MODULE){
      if(length(dtl$GENE[dtl$MODULE==res$MODULE[i]])<10){
        inc <- 1:length(dtl$GENE[dtl$MODULE==res$MODULE[i]])
      } else{
        inc <- 1:10
      }
      res$TOPGENES[i] <- 
        do.call("paste",
                c(dtl$GENE[dtl$MODULE==res$MODULE[i]][inc],
                  list("sep"=";")))
      res$TOPMARKERS[i] <- 
        do.call("paste",
                c(dtl$MARKER[dtl$MODULE==res$MODULE[i]][inc],
                  list("sep"=";")))
      res$TOPVALUES[i] <- 
        do.call("paste",
                c(dtl$VALUE[dtl$MODULE==res$MODULE[i]][inc],
                  list("sep"=";")))
    } else{
      res$TOPGENES[i] <- ""
      res$TOPMARKERS[i] <- ""
      res$TOPVALUES[i] <- ""
    }
  }
  
  # Save full results.
  jdir <- file.path(job$folder, "msea")
  fname <- paste(job$label, ".results.txt", sep="")
  res <- res[,c("MODULE","P","FREQ","NGENES","NMARKER",
                "DENSITY","FDR","DESCR","ZSCORE","CHI",
                "CHI_SD","CHI_SE","TOPGENES","TOPMARKERS",
                "TOPVALUES")]
  job$topgenes <- res[,c("MODULE","TOPGENES","TOPMARKERS","TOPVALUES")]
  job$resultfile <- file.path(jdir, fname)
  job$msea_results <- res
  tool.save(frame=res, file=fname, directory=jdir)
  
  return(job)
}

ssea.finish.param <- function(job) {
  inputs <- data.frame(
    "Input or parameter" = c(
      "Marker association file",
      "Marker mapping file",
      "Marker set file",
      "Permutation type",
      "Number of permutations",
      "Random seed",
      "Trim",
      "Minimum gene count",
      "Maximum gene count",
      "Maximum overlap between genes"
    ),
    "Value" = c(
      basename(job$marfile),
      ifelse(is.null(job$genfile), NA, basename(job$genfile)),
      basename(job$modfile),
      job$permtype,
      job$nperm,
      job$seed,
      job$trim,
      job$mingenes,
      job$maxgenes,
      job$maxoverlap
    ), 
   check.names = FALSE
  )
  
  # Save contents.
  jdir <- file.path(job$folder, "msea")
  fname <- paste(job$label, ".param.txt", sep="")
  tool.save(frame=inputs, file=fname, directory=jdir)
  
  return(job)
}


#
# Merge multiple MSEA results into meta MSEA.
#
# Input:
#   jobs    MSEA list objects
#   label   label for meta job
#   folder  parent folder for meta job
#
# Output:
#   meta    MSEA list object
#
# Written by Ville-Petteri Makinen 2013, Modified by Le Shu 2015
#
ssea.meta <- function(jobs, label, folder) {
  
  # Create meta job.
  #cat("\nMerging jobs...\n")
  meta <- list()
  meta$label <- label
  meta$folder <- folder 
  meta$modfile <- "undefined"
  meta$genfile <- "undefined"
  meta$marfile <- "undefined"
  
  jdir <- file.path(folder, "meta")
  if(!dir.exists(jdir)) dir.create(path=jdir, recursive=TRUE)
  if(file.access(jdir, 2) != 0)
    stop("Cannot access '" + plan$folder + "'.")
  
  meta <- ssea.start.configure(meta)
  
  # Collect data.
  meta$results <- data.frame()
  meta$modinfo <- data.frame()
  meta$moddata <- data.frame()
  meta$gendata <- data.frame()
  meta$locdata <- data.frame()
  forhet <- list()
  for(k in 1:length(jobs)) {
    job <- jobs[[k]]
    results <- job$results
    modinfo <- job$modinfo
    moddata <- job$moddata
    gendata <- job$gendata
    locdata <- job$locdata
    
    # Restore original identities.
    results$MODULE <- job$modules[results$MODULE]
    modinfo$MODULE <- job$modules[modinfo$MODULE]
    moddata$MODULE <- job$modules[moddata$MODULE]
    moddata$GENE <- job$genes[moddata$GENE]
    gendata$GENE <- job$genes[gendata$GENE]
    gendata$MARKER <- job$loci[gendata$MARKER]
    locdata$MARKER <- job$loci[locdata$MARKER]
    
    # Update meta sets.
    meta$results <- rbind(meta$results, results)
    meta$modinfo <- rbind(meta$modinfo, modinfo)
    meta$moddata <- rbind(meta$moddata, moddata)
    meta$gendata <- rbind(meta$gendata, gendata)
    meta$locdata <- rbind(meta$locdata, locdata)
    
    # get data for het stats
    job$zscores[is.nan(job$zscores)] <- 0
    job$var[is.na(job$var)] <- 1
    forhet[[paste0("Job_",k)]] <- data.frame("MODULE"=job$modules,
                                             "Zscore"=job$zscore,
                                             "Var_Inv"=1/(job$var))
  }
  
  # Remove duplicate rows (non-numeric values only).
  meta$modinfo <- unique(meta$modinfo)
  meta$moddata <- unique(meta$moddata)
  meta$gendata <- unique(meta$gendata)
  
  # Determine identities.
  meta$modules <- unique(meta$moddata$MODULE)
  meta$genes <- unique(meta$gendata$GENE)
  meta$loci <- unique(meta$gendata$MARKER)
  
  # Convert identities to indices.
  meta$results <- ssea.start.identify(meta$results, "MODULE", meta$modules)
  meta$modinfo <- ssea.start.identify(meta$modinfo, "MODULE", meta$modules)
  meta$moddata <- ssea.start.identify(meta$moddata, "MODULE", meta$modules)
  meta$moddata <- ssea.start.identify(meta$moddata, "GENE", meta$genes)
  meta$gendata <- ssea.start.identify(meta$gendata, "GENE", meta$genes)
  meta$gendata <- ssea.start.identify(meta$gendata, "MARKER", meta$loci)
  meta$locdata <- ssea.start.identify(meta$locdata, "MARKER", meta$loci)
  
  # Convert marker values to z-scores.
  values <- meta$locdata$VALUE
  qvals <- tool.unify(values)
  qvals <- pmax(qvals, .Machine$double.xmin)
  qvals <- pmin(qvals, (1.0 - .Machine$double.eps))
  zvals <- qnorm(qvals)
  
  # Merge matching markers.
  st <- tool.aggregate(meta$locdata$MARKER)
  blocks <- st$blocks
  for(k in 1:length(blocks)) {
    rows <- blocks[[k]]
    z <- sum(zvals[rows])
    zvals[rows] <- z/sqrt(length(rows))
  }
  
  # Convert back to original data space.
  meta$locdata$VALUE <- quantile(values, pnorm(zvals))
  meta$locdata <- unique(meta$locdata)
  
  # Construct hierarchical representation.
  cat("\nPreparing data structures...\n")
  ngens <- length(meta$genes) 
  nmods <- length(meta$modules)
  meta$database <- ssea.prepare.structure(meta$moddata, meta$gendata,
                                          nmods, ngens)
  
  # Determine test cutoffs.
  lengths <- meta$database$modulelengths
  mu <- median(lengths[which(lengths > 0)])
  meta$quantiles <- seq(0.5, (1.0 - 1.0/mu), length.out=10)
  
  # Calculate hit counts.
  NMARKER <- length(meta$loci)
  hits <- ssea.prepare.counts(meta$locdata, NMARKER, meta$quantiles)
  meta$database <- c(meta$database, hits)
  
  # Check result values.
  meta$results <- meta$results[,c("MODULE", "P")]
  pvalues <- pmax(meta$results$P, .Machine$double.xmin)
  pvalues <- pmin(pvalues, (1.0 - .Machine$double.eps))
  meta$results$P <- pvalues
  
  # Calculate meta P-values.
  cat("\nPostprocessing meta results...\n")
  st <- tool.aggregate(meta$results$MODULE)
  blocks <- st$blocks
  for(k in 1:length(blocks)) {
    rows <- blocks[[k]]
    tmp <- meta$results[rows,]
    z <- qnorm(tmp$P)
    z <- sum(z)/sqrt(length(z))
    meta$results[rows,"P"] <- NA
    meta$results[rows[1],"P"] <- pnorm(z)
  }
  
  # restore orig identities
  meta$results$MODULE <- meta$modules[meta$results$MODULE]
  meta$db_df <- data.frame("MODULE"=meta$modules,
                           "modulesizes"=meta$database$modulesizes,
                           "modulelengths"=meta$database$modulelengths,
                           "moduledensities"=meta$database$moduledensities)
  
  # for each module get Z scores and inverse of variance of chis
  meta$results$Cochran.Q <- vapply(meta$results$MODULE, function(x){
    xval <- c()
    weightval <- c()
    for(j in names(forhet)){
      xval <- c(xval, forhet[[j]]$Zscore[forhet[[j]]$MODULE==x])
      weightval <- c(weightval, forhet[[j]]$Var_Inv[forhet[[j]]$MODULE==x])
    }
    weightval[weightval==1] <- 0
    cochran_res <- tool.cochranQ(x = xval, weights = weightval)
    if(is.nan(cochran_res["Q"])){
      return(NA)
    } else {
      return(cochran_res["Q"])
    }
  }, FUN.VALUE = numeric(1))
  meta$results$Cochran.P <- vapply(meta$results$MODULE, function(x){
    xval <- c()
    weightval <- c()
    for(j in names(forhet)){
      xval <- c(xval, forhet[[j]]$Zscore[forhet[[j]]$MODULE==x])
      weightval <- c(weightval, forhet[[j]]$Var_Inv[forhet[[j]]$MODULE==x])
    }
    weightval[weightval==1] <- 0
    cochran_res <- tool.cochranQ(x = xval, weights = weightval)
    return(cochran_res["p-value"])
  }, FUN.VALUE = numeric(1))
  degfrd <- length(jobs) - 1
  meta$results$Cochran.DF <- degfrd
  meta$results$I2 <- vapply(meta$results$Cochran.Q, function(x){
    if(is.na(x)){
      return(NA)
    } else{
      i2 = 100*((x-degfrd)/x)
      if(i2<0){
        return(0)
      } else{
        return(i2)
      }
    }
  }, FUN.VALUE=numeric(1))
  
  meta$results$Cochran.Q[is.nan(meta$results$Cochran.Q)] <- NA
  meta$results$Cochran.P[is.nan(meta$results$Cochran.P)] <- NA
  
  # Finish and save statistics.
  meta$results <- na.omit(meta$results)
  meta <- ssea.finish.fdr(meta, jobs = jobs)
  #meta <- ssea.finish.details(meta)
  return(meta)
}
#
# Prepare an indexed database for Marker set enrichment analysis.
#
# Input:
#   job$modules      module identities as characters
#   job$genes        gene identities as characters
#   job$loci         marker identities as characters
#   job$moddata      preprocessed module data (indexed identities)
#   job$gendata      preprocessed mapping data (indexed identities)
#   job$locdata      preprocessed marker data (indexed identities)
#   job$mingenes     minimum module size allowed
#   job$maxgenes     maximum module size allowed
#   job$maxoverlap   maximum module overlap allowed (use 1.0 to skip)
#
# Optional input:
#   job$quantiles    quantile points for test statistic
#
# Output:
#   job$modules      finalized module names
#   job$moddata      finalized module data
#   job$gendata      finalized mapping data
#   job$locdata      finalized marker data
#   job$quantiles    verified quantile points
#   job$database$modulesizes
#                    gene counts for modules
#   job$database$modulelengths
#                    distinct marker counts for modules
#   job$database$moduledensities
#                    ratio between distinct and non-distinct markers
#   job$database$genesizes
#                    marker count for each gene
#   job$database$module2genes
#                    gene lists for each module
#   job$database$gene2loci
#                    marker lists for each gene
#   job$database$locus2row
#                    row indices in the marker data frame for each marker
#   job$database$observed
#                    matrix of observed counts of values that exceed each
#                    quantile point for each marker
#   job$database$expected
#                    1.0 - quantile points
#
# The database uses indexed identities for modules, genes and marker.
# Output also includes all the other items from input list.
#
# Written by Ville-Petteri Makinen 2013
#
ssea.prepare <- function(job) {
  cat("\nPreparing data structures...\n")
  
  # Remove extreme modules.
  st <- tool.aggregate(job$moddata$MODULE)
  mask <- which((st$lengths >= job$mingenes) &
                  (st$lengths <= job$maxgenes))
  pos <- match(job$moddata$MODULE, st$labels[mask])
  job$moddata <- job$moddata[which(pos > 0),]
  removed <- setdiff(st$labels, unique(job$moddata$MODULE))
  cat("The following modules were removed because they did not\npass the mingenes and maxgenes parameters:\n", 
      do.call("paste", c(job$modules[as.numeric(removed)],list("sep"=", "))), "\n")
  if(nrow(job$moddata)==0) stop("No modules passing mingenes and maxgenes.")
  
  # Construct hierarchical representation.
  ngens <- length(job$genes)
  nmods <- length(job$modules)
  db <- ssea.prepare.structure(job$moddata, job$gendata, nmods, ngens)
  
  # Determine test cutoffs.
  if(is.null(job$quantiles)) {
    lengths <- db$modulelengths
    mu <- median(lengths[which(lengths > 0)])
    job$quantiles <- seq(0.5, (1.0 - 1.0/mu), length.out=10)
  }
  
  # Calculate hit counts.
  NMARKER <- length(job$loci)
  hits <- ssea.prepare.counts(job$locdata, NMARKER, job$quantiles)
  db <- c(db, hits)
  
  # Return results.
  job$database <- db
  remove(db)
  gc(FALSE)
  nmem <- (object.size(job))*(0.5^20)
  cat("Job: ", nmem, " Mb\n", sep="")
  return(job)
}

#----------------------------------------------------------------------------

ssea.prepare.structure <- function(moddata, gendata, nmods, ngens) {
  
  # Prepare list structures.
  genlists <- list()
  loclists <- list()
  for(k in 1:nmods) genlists[[k]] <- integer()
  for(k in 1:ngens) loclists[[k]] <- integer()
  modsizes <- rep(0, nmods)
  modlengths <- rep(0, nmods)
  moddensities <- rep(0.0, nmods)
  gensizes <- rep(0, ngens)
  
  # Collect row indices for each module.
  st <- tool.aggregate(moddata$MODULE)
  keys <- as.integer(st$labels)
  blocks <- st$blocks
  
  # Collect gene lists.
  genes <- as.integer(moddata$GENE)
  for(k in 1:length(blocks)) {
    key <- keys[k]
    rows <- blocks[[k]]
    members <- unique(genes[rows])
    genlists[[key]] <- as.integer(members)
    modsizes[[key]] <- length(members)
  }
  
  # Collect row indices for each gene.
  st <- tool.aggregate(gendata$GENE)
  keys <- as.integer(st$labels)
  blocks <- st$blocks
  
  # Collect marker lists.
  loci <- as.integer(gendata$MARKER)
  for(k in 1:length(blocks)) {
    key <- keys[k]
    rows <- blocks[[k]]
    members <- unique(loci[rows])
    loclists[[key]] <- as.integer(members)
    gensizes[[key]] <- length(members)
  }
  
  # Count distinct markers in each module.
  for(k in which(modsizes > 0)) {
    locset <- integer()
    genset <- genlists[[k]]
    for(i in genset)
      locset <- c(locset, loclists[[i]])
    modlengths[[k]] <- length(unique(locset))
    moddensities[[k]] <- modlengths[[k]]/length(locset)
  }
  
  # Check data integrity.
  if(sum(gensizes == 0) > 0) stop("Incomplete locus data.")
  if(length(gensizes) != ngens) stop("Inconsistent gene data.")
  if(length(modsizes) != nmods) stop("Inconsistent module data.")
  
  # Return results.
  res <- list()
  res$modulesizes <- modsizes
  res$modulelengths <- modlengths
  res$moduledensities <- moddensities
  res$genesizes <- gensizes
  res$module2genes <- genlists
  res$gene2loci <- loclists
  return(res)
}

#----------------------------------------------------------------------------

ssea.prepare.counts <- function(locdata, NMARKER, quantiles) {
    
    # Make sure there are at least two points to prevent
    # R automagic on matrices to mess things up.
    if(length(quantiles) < 2) quantiles <- rep(quantiles, 2)
    
    # Create mapping table.
    nrows <- nrow(locdata)
    locmap <- rep(0, NMARKER)
    locmap[locdata$MARKER] <- (1:nrows)
    
    # Convert values to standardized range.
    values <- tool.unify(locdata$VALUE)
    
    # Create bit matrix of values above quantiles.
    nquant <- length(quantiles)
    bits <- matrix(data=FALSE, nrow=nrows, ncol=nquant)
    for(i in 1:nrows)
      bits[i,] <- (values[i] > quantiles)
    
    # Return results.
    res <- list()
    res$locus2row <- locmap
    res$observed <- bits
    res$expected <- (1.0 - quantiles)
    return(res)
}
#
# Create a job for Marker set enrichment analysis
#
# Input:
#   plan$label       unique identifier for the analysis
#   plan$folder      output folder for results
#   plan$modfile     path to module file
#                    columns: MODULE GENE
#   plan$marfile     path to marker file
#                    columns: MARKER VALUE
#   plan$genfile     path to gene file
#                    columns: GENE MARKER
#
# Optional input:
#   plan$inffile     path to module info file
#                    columns: MODULE DESCR
#   plan$seed        seed for random number generator
#   plan$permtype    'gene' for gene-level, 'marker' for marker-level
#   plan$nperm       maximum number of random permutations
#   plan$mingenes    minimum number of genes per module (after merging)
#   plan$maxgenes    maximum number of genes per module
#   plan$quantiles   cutoffs for test statistic
#   plan$maxoverlap  maximum overlap allowed between genes
#   plan$trim        percentile of genes taken from beginning and end of 
#                    trait associations to avoid signal inflation of null
#                    background in gene permutation
#
# Output:
#   job$modules       module identities as characters
#   job$genes         gene identities as characters
#   job$loci          marker identities as characters
#   job$moddata       preprocessed module data (indexed identities)
#   job$modinfo       preprocessed module info (indexed identities)
#   job$gendata       preprocessed mapping data (indexed identities)
#   job$locdata       preprocessed marker data (indexed identities)
#   job$geneclusters  genes with shared markers
#
#  Output also includes the items from input list.
#
# Written by Ville-Petteri Makinen 2013, Modified by Le Shu 2015, Modified by Jessica Ding 2021
#
ssea.start <- function(plan) {
    
    jdir <- file.path(plan$folder, "msea")
    if(!dir.exists(jdir)) dir.create(path=jdir, recursive=TRUE)
    if(file.access(jdir, 2) != 0)
        stop("Cannot access '" + plan$folder + "'.")
    
    cat("\nMSEA Version:Oct2023\n")
    # Check parameters.
    job <- ssea.start.configure(plan)
    
    # Import gene sets.
    cat("\nImporting modules...\n")
    modinfo <- tool.read(plan$inffile, c("MODULE", "DESCR"))
    moddata <- tool.read(plan$modfile, c("MODULE", "GENE"))
    moddata <- unique(na.omit(moddata))
    modinfo <- unique(na.omit(modinfo))
    if(nrow(modinfo) > 0) print(summary(modinfo))
    print(summary(moddata))
    
    # Add slots for control modules.
    modules <- unique(moddata$MODULE)
    if((sum(modules == "_ctrlA") > 0) | (sum(modules == "_ctrlA") > 0))
      stop("Module names '_ctrlA' and '_ctrlB' are reserved.")
    modules <- c(modules, "_ctrlA", "_ctrlB")
    tmp <- data.frame(MODULE=c("_ctrlA", "_ctrlB"),
    DESCR=c("Top genes", "Top genes (module members)"))
    modinfo <- rbind(modinfo, tmp)
    
    # Import marker values.
    cat("\nImporting marker values...\n")
    if(!is.null(job$locfile)) job$marfile <- job$locfile
    locdata <- tool.read(job$marfile, c("MARKER", "VALUE"))
    locdata$VALUE <- as.double(locdata$VALUE)
    rows <- which(0*(locdata$VALUE) == 0)
    locdata <- unique(na.omit(locdata[rows,]))
    print(summary(locdata))
    
    # Import mapping data.
    cat("\nImporting mapping data...\n")
    gendata <- tool.read(job$genfile, c("GENE", "MARKER"))
    gendata <- unique(na.omit(gendata))
    print(summary(gendata))
    
    # Remove genes with no marker values.
    pos <- match(gendata$MARKER, locdata$MARKER)
    gendata <- gendata[which(pos > 0),]
    
    # Merge overlapping genes.
    cat("\nMerging genes containing shared markers...\n")
    gendata <- tool.coalesce(items=gendata$MARKER, groups=gendata$GENE,
                             rcutoff=job$maxoverlap)
    job$geneclusters <- gendata[,c("CLUSTER","GROUPS")]
    job$geneclusters <- unique(job$geneclusters)
    
    # Update gene symbols.
    moddata <- ssea.start.relabel(moddata, gendata)
    gendata <- unique(gendata[,c("GROUPS", "ITEM")])
    names(gendata) <- c("GENE", "MARKER")
    
    # Collect identities.
    job$modules <- modules
    job$loci <- intersect(gendata$MARKER, locdata$MARKER)
    pos <- match(gendata$MARKER, job$loci)
    job$genes <- gendata[which(pos > 0), "GENE"]
    job$genes <- unique(job$genes)
    
    # Exclude missing data and factorize identities.
    job$modinfo <- ssea.start.identify(modinfo, "MODULE", job$modules)
    job$moddata <- ssea.start.identify(moddata, "MODULE", job$modules)
    job$moddata <- ssea.start.identify(job$moddata, "GENE", job$genes)
    removed <- setdiff(1:(length(job$modules)-2), unique(job$moddata$MODULE))
    cat("The following modules were removed because they did not\nhave any genes from the association data:\n", 
        do.call("paste", c(job$modules[as.numeric(removed)],list("sep"=", "))), "\n")
    job$gendata <- ssea.start.identify(gendata, "GENE", job$genes)
    job$gendata <- ssea.start.identify(job$gendata, "MARKER", job$loci)
    job$locdata <- ssea.start.identify(locdata, "MARKER", job$loci)
    
    # Show job size.
    nmem <- (object.size(job))*(0.5^20)
    cat("Job: ", nmem, " Mb\n", sep="")
    
    # Clean-up.
    remove(modinfo)
    remove(moddata)
    remove(gendata)
    remove(locdata)
    gc(FALSE)
    return(job)
}

#----------------------------------------------------------------------------

ssea.start.configure <- function(plan) {
  #bypass if running Meta-MSEA
  if (!is.null(plan$folder) & 
      !is.null(plan$label) & 
      plan$marfile == "undefined" &
      plan$modfile == "undefined"){
    plan$permtype <- "gene"
    plan$nperm <- 20000
    plan$seed <- 1
    plan$mingenes <- 10
    plan$maxgenes <- 500
    plan$maxoverlap <- 0.33
    cat("\nRunning Meta-MSEA...\n")
  } else {
    plan$stamp <- Sys.time()
    
    if(is.null(plan$permtype)) plan$permtype <- "gene"
    if(is.null(plan$nperm)) plan$nperm <- 20000
    if(is.null(plan$seed)) plan$seed <- 1
    if(is.null(plan$mingenes)) plan$mingenes <- 10
    if(is.null(plan$maxgenes)) plan$maxgenes <- 500
    if(is.null(plan$maxoverlap)) plan$maxoverlap <- 0.33
    if(is.null(plan$trim)) plan$trim <- 0.002
    
    #Resolve inconsistencies
    if(plan$permtype == "locus") plan$permtype = "marker"
    if ("LOCUS" %in% names(read.delim(plan$marfile))){
      dir.create("tmp", showWarnings = FALSE)
      tmploci=tool.read(plan$marfile,c("LOCUS","VALUE"))
      tmpgene=tool.read(plan$genfile,c("GENE","LOCUS"))
      names(tmploci)=c("MARKER","VALUE")
      names(tmpgene)=c("GENE","MARKER")
      write.table(tmploci,"tmp/marker.txt",quote = FALSE,row.names = FALSE,
                  sep = "\t")
      write.table(tmpgene,"tmp/gene.txt",quote = FALSE,row.names = FALSE,
                  sep = "\t")
      plan$marfile = "tmp/marker.txt"
      plan$genfile = "tmp/gene.txt"
    }
    
    if(is.null(plan$folder)) stop("No output folder.")
    if(is.null(plan$label)) stop("No job label.")
    if(is.null(plan$modfile)) stop("No module file.")
    if(is.null(plan$genfile)){
      cat("\nNo gene file provided.\nAssuming non-GWAS enrichment.\n")
      # make 'fake' mappping file
      locdata <- tool.read(plan$marfile, c("MARKER", "VALUE"))
      locdata$VALUE <- as.double(locdata$VALUE)
      rows <- which(0*(locdata$VALUE) == 0)
      locdata <- unique(na.omit(locdata[rows,]))
      gendata <- data.frame("GENE"=locdata$MARKER,
                            "MARKER"=locdata$MARKER)
      dir.create("tmp", showWarnings = FALSE)
      write.table(gendata,"tmp/gene.txt",quote = FALSE,row.names = FALSE,
                  sep = "\t")
      plan$genfile = "tmp/gene.txt"
      
      # change to appropriate parameters
      plan$permtype <- "marker"
      plan$maxoverlap <- 1
    }
    if(is.null(plan$marfile)) stop("No marker file.")
    
    # accommodate utf converted files
    data <- read.delim(plan$marfile)
    if(grepl("X...", colnames(data)[1]) | grepl("X...", colnames(data)[2])){
      cat("Accommodating utf converted marker file...\n")
      colnames(data) <- gsub("X...","",colnames(data))
      system(paste0("chmod +x ",plan$marfile))
      write.table(data, plan$marfile, quote = F,row.names = F,sep = "\t")
    }
    data <- read.delim(plan$genfile)
    if(grepl("X...", colnames(data)[1]) | grepl("X...", colnames(data)[2])){
      cat("Accommodating utf converted gene file...\n")
      colnames(data) <- gsub("X...","",colnames(data))
      system(paste0("chmod +x ",plan$genfile))
      write.table(data, plan$genfile, quote = F,row.names = F,sep = "\t")
    }
    
    cat("\nParameters:\n")
    cat("  Permutation type: ", plan$permtype, "\n", sep="")  
    cat("  Permutations: ", plan$nperm, "\n", sep="")
    cat("  Random seed: ", plan$seed, "\n", sep="")
    cat("  Trim: ", plan$trim, "\n", sep="")
    cat("  Minimum gene count: ", plan$mingenes, "\n", sep="")
    cat("  Maximum gene count: ", plan$maxgenes, "\n", sep="")
    cat("  Maximum overlap between genes: ", plan$maxoverlap, "\n", 
        sep="")
    
    if(is.null(plan$quantiles) == FALSE) {
      cat("  Test quantiles:");
      for(q in plan$quantiles)
        cat(sprintf(" %.2f", 100*q), "%", sep="")
      cat("\n")
    }
  }
  
  return(plan)
}

#----------------------------------------------------------------------------

ssea.start.relabel <- function(dat, grp) {
    
    # New gene group symbols.
    oldgenes <- character()
    newgenes <- character()
    syms <- unique(grp[,c("CLUSTER","GROUPS")])
    rows <- which(syms$CLUSTER != syms$GROUPS)
    for(i in rows) {
        g <- syms[i,"GROUPS"]
        a <- strsplit(g, ",", fixed=TRUE)
        a <- a[[1]]
        b <- rep(g, length(a))
        oldgenes <- c(oldgenes, a)
        newgenes <- c(newgenes, b)
    }
    
    # Update dataset.
    if(length(newgenes) < 1) return(dat)
    pos <- match(dat$GENE, oldgenes)
    rows <- which(pos > 0)
    dat[rows,"GENE"] <- newgenes[pos[rows]]
    return(unique(dat))
}

#----------------------------------------------------------------------------

ssea.start.identify <- function(dat, varname, labels) {
    if(nrow(dat) < 1) return(dat)
    
    # Find matching identities.
    pos <- match(dat[,varname], labels)
    rows <- which(pos > 0)
    
    # Select subset.
    dat[,varname] <- pos
    res <- dat[rows,]
    return(res)
}

#
# Sort an array and find the indices of blocks with the same values.
# The second argument sets the minimum block size to be included.
#
# Output list:
#   res$labels     shared values within blocks
#   res$lengths    numbers of entries in blocks
#   res$blocks     integer arrays of entry positions within blocks
#   res$ranks      entry positions included in blocks
#
# Written by Ville-Petteri Makinen 2013
#
tool.aggregate <- function(entries, limit=1) {
    res <- list()
    res$blocks <- list()
    
    # Check input size.
    nelem <- length(entries)
    if(nelem < 1) stop("Unusable input.")
    
    # Factorize entries.
    entries <- as.factor(entries)
    
    # Convert factors to integers.
    elevels <- levels(entries)
    entries <- as.integer(entries)
    
    # Sort entries.
    mask <- order(entries)
    
    # Remove missing entries.
    rows <- which(entries > 0)
    mask <- intersect(mask, rows)
    if(length(mask) < 1) stop("Unusable input.")
    
    # Single entry.
    if(length(mask) < 2) {
        label <- na.omit(elevels)
        res$labels <- label
        res$lengths <- 1
        res$blocks[[1]] <- 1
        res$ranks <- 1
        return(res)
    }
    
    # Starting point.
    nstack <- 1
    stack <- integer(length=nelem)
    stack[nstack] <- mask[1]
    prev <- entries[mask[1]]
    mask <- mask[2:nelem]
    
    # Find segments of identical entries.
    buffer <- list()
    subsets <- list()
    for(i in mask) {
        if(entries[i] != prev) {
            
            # Clear stack.
            if(nstack >= limit) {
                ind <- (length(buffer) + 1)
                buffer[[ind]] <- stack[1:nstack]
            }
            nstack <- 0
            
            # Buffering for speed-up.
            if(length(buffer) > 120) {
                subsets <- c(subsets, buffer)
                buffer <- list()
            }
        }
        
        # Add item to stack.
        nstack <- (nstack + 1)
        stack[nstack] <- i
        prev <- entries[i]
    }
    
    # Clear last item(s).
    if(nstack >= limit) {
        ind <- (length(buffer) + 1)
        buffer[[ind]] <- stack[1:nstack]
    }
    
    # Clear buffer.
    if(length(buffer) > 0)
      subsets <- c(subsets, buffer)
    
    # Check if any subsets.
    nuniq <- length(subsets)
    if(nuniq < 1) return(NULL)
    
    # Determine additional attributes.
    loci <- rep(NA, nelem)
    sizes <- integer(nuniq)
    identities <- integer(nuniq)
    for(k in 1:nuniq) {
        mask <- as.integer(subsets[[k]])
        identities[k] <- entries[mask[1]]
        sizes[k] <- length(mask)
        loci[mask] <- mask
    }
    
    # Finish.
    res$labels <- elevels[identities]
    res$lengths <- sizes
    res$blocks <- subsets
    res$ranks <- na.omit(loci)
    return(res)
}
#
# Use hierarchical clustering to assign nodes into clusters.
#
# Input:
#   edges    data.frame:
#            A        item name
#            B        item name
#            POSa     item name rank
#            POSb     item name rank
#            R        overlap between A and B
#  cutoff    maximum overlap not considered clustered
#
# Output:
#   res      data frame
#            CLUSTER  cluster rank
#            NODE     item name
#
# Written by Ville-Petteri Makinen 2013
#
tool.cluster <- function(edges, cutoff=NULL) {
    
    # Default output.
    a <- edges$A
    b <- edges$B
    labels <- unique(c(a, b))
    res <- data.frame(CLUSTER=labels, stringsAsFactors=FALSE)
    res$NODE <- labels
    
    # Check if clustering is needed.
    r <- as.double(edges$R)
    posA <- as.integer(edges$POSa)
    posB <- as.integer(edges$POSb)
    ndim <- max(c(posA, posB))
    if(sum(posA != posB) < 1) return(res)
    if(max(r) <= 0.0) return(res)
    
    # Allocate distance matrix.
    mtx <- matrix(data=0.0, nrow=ndim, ncol=ndim)
    labels <- rep(NA, ndim)
    
    # Collect group labels.
    for(k in 1:nrow(edges)) {
        labels[posA[k]] <- a[k]
        labels[posB[k]] <- b[k]
    }
    
    # Recreate matrix form.
    for(k in 1:nrow(edges)) {
        i <- posA[k]
        j <- posB[k]
        mtx[i,j] <- r[k]
        mtx[j,i] <- r[k]
    }
    
    # Hierarchical clustering.
    d <- as.dist(1 - mtx)
    tree <- hclust(d)
    
    # Height cutoff.
    hlim <- max(tree$height)
    if(is.null(cutoff) == FALSE) hlim <- (1.0 - cutoff)
    
    # Find clusters.
    clusters <- tool.cluster.static(tree, hlim)
    
    # Enumerate clusters with singletons included.
    mask <- which(clusters == 0)
    clusters[mask] <- -mask
    clusters <- as.factor(clusters)
    clusters <- as.integer(clusters)
    
    # Create supergroups.
    res <- data.frame(CLUSTER=clusters, stringsAsFactors=FALSE)
    res$NODE <- labels[1:length(clusters)]
    return(res)
}

#---------------------------------------------------------------------------

tool.cluster.static <- function(dendro, hlim) {
    merged <- dendro$merge
    heights <- dendro$height
    ndim <- (length(heights) + 1)
    
    # Assign clusters by static cut.
    clusters <- rep(0, ndim)
    for(i in which(heights <= hlim)) {
        a <- as.integer(merged[i, 1])
        b <- as.integer(merged[i, 2])
        
        # Put previous cluster (if any) in A.
        if(a < 0) {
            tmp <- a
            a <- b
            b <- tmp
        }
        
        # De novo merging.
        if(a < 0) {
            clusters[-a] <- i
            clusters[-b] <- i
            next
        }
        
        # Merge with previous cluster.
        if(b < 0) {
            mask <- which(clusters == a)
            mask <- c(mask, -b)
            clusters[mask] <- i
            next;
        }
        
        # Merge two clusters.
        mask <- which((clusters == a) | (clusters == b))
        clusters[mask] <- i
    }
    return(clusters)
}
# Calculate overlaps between groups of items.
#
# Input:
#   items     array of item identities
#   groups    array of group identities for items
#
# Optional input:
#   rcutoff   maximum overlap not coalesced
#   ncore     minimum number of items required for trimming
#
# Output:
#   res       data frame:
#             CLUSTER  cluster identities
#             ITEM     item identities
#             GROUPS   comma separated group identities
#
# Due to trimming, output may contain fewer distinct items.
# Cluster identities are a subset of group identities.
#  
# Written by Ville-Petteri Makinen 2013, Modified by Jessica Ding 2021
#
tool.coalesce <- function(items, groups, rcutoff=0.0, ncore=NULL) {
    
    # Check arguments.
    if(length(items) != length(groups))
      stop("Incompatible inputs.")
    
    # Default output.
    res <- data.frame(CLUSTER=groups, GROUPS=groups,
                      ITEM=items, stringsAsFactors=FALSE)
    if(rcutoff >= 1.0) return(res)
    if(length(items)==length(unique(items))){
        cat("Genes are all unique - no need to merge modules\n")
        return(res)
    }
    
    # Check that group names are usable.
    grlabels <- unique(groups)
    if(is.character(grlabels)) {
        for(s in grlabels) {
            segm <- strsplit(s, ",", fixed=TRUE)
            if(length(segm[[1]]) < 2) next
            stop("Group labels must not contain ','.")
        }
    }
    
    # Determine core item set size.
    nitems <- length(items)
    if(is.null(ncore)) {
        ncore <- nitems/length(unique(groups))
        ncore <- round(ncore)
    }
    
    # Convert to integers.
    itemlev <- as.factor(items)
    grouplev <- as.factor(groups)
    members <- as.integer(itemlev)
    modules <- as.integer(grouplev)
    itemlev <- levels(itemlev)
    grouplev <- levels(grouplev)
    
    # Determine item freguencies.
    freq <- table(members)
    labels <- as.integer(names(freq))
    freq <- as.integer(freq)
    
    # Limit the number of comparisons.
    kappa <- 0
    while(TRUE) {
        kappa <- (kappa + 1)
        shared <- which(freq > kappa)
        shared <- labels[shared]
        
        # Determine groups with overlaps.
        pos <- match(members, shared)
        rows <- which(pos > 0)
        mods <- unique(modules[rows])
        if(length(mods) < 2000) break
    }
    
    # Show warning if overlaps too extensive.
    if(kappa > 1) {
        cat("WARNING! Limited overlap analysis due ")
        cat("to large number of groups.\n")
    }
    
    # Determine subset with shared items.
    pos <- match(modules, mods)
    incl <- which(pos > 0)
    excl <- setdiff((1:nitems), incl)
    
    # Find and trim clusters.
    if(length(incl)==0) return(res)
    res <- tool.coalesce.exec(members[incl], modules[incl], rcutoff, ncore)
    res <- rbind(res, tool.coalesce.exec(members[excl], modules[excl], 1.0))
    
    # Convert identities back to original.
    res$ITEM <- itemlev[res$ITEM]
    res$CLUSTER <- grouplev[res$CLUSTER]
    groupdat <- rep("", nrow(res))
    groupsets <- as.character(res$GROUPS)
    for(i in 1:nrow(res)) {
        gset <- strsplit(groupsets[i], ",", fixed=TRUE)
        gset <- as.integer(gset[[1]])
        groupdat[i] <- paste(grouplev[gset], collapse=",")
    }
    res$GROUPS <- groupdat
    return(res)
}

#---------------------------------------------------------------------------

tool.coalesce.exec <- function(items, groups, rcutoff, ncore) {
    if(is.numeric(items) == FALSE) stop("Unusable input.")
    if(is.numeric(groups) == FALSE) stop("Unusable input.")
    
    # Default output.
    res <- data.frame(CLUSTER=groups, GROUPS=groups, ITEM=items,
    stringsAsFactors=FALSE)
    if(rcutoff >= 1.0) return(res)
    
    # Iterative merging and trimming.
    res$COUNT <- 0.0
    while(TRUE) {
        clust <- tool.coalesce.find(res, rcutoff)
        if(is.null(clust)) break
        res <- tool.coalesce.merge(clust, ncore)
    }
    
    # Select columns.
    res$COUNT <- NULL
    res$CLUSTER <- 0
    itemdat <- res$ITEM
    clustdat <- res$CLUSTER
    groupdat <- res$GROUPS
    
    # Select representative label for clusters.
    st <- tool.aggregate(res$GROUP)
    blocks <- st$blocks
    labels <- st$labels
    for(k in 1:length(blocks)) {
        rows <- blocks[[k]]
        locals <- itemdat[rows]
        nodes <- strsplit(labels[[k]], ",", fixed=TRUE)
        nodes <- as.numeric(nodes[[1]])
        overlaps <- rep(0, length(nodes))
        for(i in 1:length(nodes)) {
            mask <- which(groups == nodes[i])
            shared <- intersect(items[mask], locals)
            overlaps[i] <- length(shared)
        }
        mask <- order(overlaps, decreasing=TRUE)
        clustdat[rows] <- nodes[mask[1]]
        groupdat[rows] <- paste(nodes[mask], collapse=",")
    }
    
    # Finish results.
    res$CLUSTER <- clustdat
    res$GROUPS <- groupdat
    return(res)
}

#---------------------------------------------------------------------------

tool.coalesce.find <- function(data, rmax) {
    
    # Harmonize column names.
    data <- data[,c("GROUPS", "ITEM", "COUNT")]
    names(data) <- c("NODE", "ITEM", "COUNT")
    
    # Find clusters.
    edges <- tool.overlap(items=data$ITEM, groups=data$NODE)
    clustdat <- tool.cluster(edges, cutoff=rmax)
    nclust <- length(unique(clustdat$CLUSTER))
    nnodes <- length(unique(clustdat$NODE))
    if(nclust >= nnodes) return(NULL)
    
    # Merge with original dataset.
    res <- merge(clustdat, data)
    return(res)
}

#---------------------------------------------------------------------------

tool.coalesce.merge <- function(data, ncore) {
    
    # Determine item clusters.
    st <- tool.aggregate(data$CLUSTER)
    blocks <- st$blocks
    
    # Trim clusters.
    res <- data.frame()
    for(k in 1:length(blocks)) {
        rows <- blocks[[k]]
        batch <- data[rows,]
        
        # Item hit counts.
        st <- tool.aggregate(batch$ITEM)
        counts <- as.integer(st$lengths)
        labels <- as.integer(st$labels)
        segm <- st$blocks
        
        # Add hit counts from previous round.
        for(j in 1:length(segm)) {
            mask <- segm[[j]]
            nj <- mean(batch[mask,"COUNT"])
            counts[j] <- sqrt(counts[j] + nj)
        }
        
        # Remove rarest items.
        levels <- sort(unique(counts))
        for(kappa in levels) {
            mask <- which(counts > kappa)
            nmask <- length(mask)
            if(nmask < ncore) break
            counts <- counts[mask]
            labels <- labels[mask]
        }
        
        # Collect groups.
        nodeset <- unique(batch$NODE)
        nodeset <- paste(nodeset, collapse=",")
        
        # Update results.
        tmp <- data.frame(GROUPS=nodeset, ITEM=labels, COUNT=counts, 
                          stringsAsFactors=FALSE)
        res <- rbind(res, tmp)
    }
    return(res)
}

#---------------------------------------------------------------------------

tool.cochranQ <-function(x, weights) {
  bar.x <- sum(weights*x)/sum(weights)
  Q <- sum(weights*((x-bar.x)^2))
  k <- length(x)
  p.value <- pchisq(Q, k-1, lower.tail = FALSE)
  output <- c(Q,p.value,as.integer(k-1))
  names(output) <- c("Q", "p-value", "df")
  return(output)
}

# 
# Written by Ville-Petteri Makinen 2013
# Modified by Jessica Ding 2023
#
tool.fdr <- function(p, nmods=NULL) {
  if(length(p) < 12){
    cat("WARNING! Number of pathways analyzed less than 10,\n")
    cat("running Bonferroni FDR correction.\n")
    fdrs <- p*nmods # bonferroni correction
    fdrs[fdrs>1] <- 1 # cap fdr at 1
    return(fdrs) 
  } else {
    return(tool.fdr.bh(p))
  }
  #return(tool.fdr.empirical(p, f))
}

#---------------------------------------------------------------------------

tool.fdr.bh <- function(p) {
    nelem <- length(p)
    
    # Revert to z-scores.
    pvals <- pmax(p, .Machine$double.xmin)
    pvals <- pmin(pvals, (1 - .Machine$double.eps))
    z <- qnorm(pvals)
    
    # Benjamini–Hochberg (1995) false discovery rate.
    fdr <- p.adjust(pvals, method="fdr")
    fdr <- pmax(fdr, .Machine$double.xmin)
    
    # Sort data points.
    mask <- order(z)
    xcoord <- z[mask]
    ycoord <- fdr[mask]
    
    # Remove steps in FDR values.
    prev <- 1
    for(i in 2:nelem) {
        if(ycoord[i] == ycoord[prev]) {
            ycoord[i] = NA
            next
        }
        window <- prev:(i-1)
        xcoord[prev] <- mean(xcoord[window])
        prev <- i
    }
    
    # Select distinct points.
    rows <- which(0*ycoord == 0)
    xcoord <- xcoord[rows]
    ycoord <- ycoord[rows]
    
    # Add sentinels.
    xcoord <- c((min(z) - 1.0), xcoord, (max(z) + 1.0))
    ycoord <- c(0.0, ycoord, 1.0)
    
    # Interpolate missing points.
    points <- approx(xcoord, ycoord, xout=z)
    return(points$y)
}

#---------------------------------------------------------------------------

tool.fdr.empirical <- function(p, f0) {
    
    # Estimate raw false discovery rates.
    f <- rank(p)/length(p)
    fdr <- pmin(pmax(f0/f, p), 1)
    
    # Pre-defined sections for sampling the FDR curve.
    p <- pmin(p, (1.0 - .Machine$double.eps))
    p <- pmax(p, .Machine$double.xmin)
    z <- qnorm(p, lower.tail=TRUE)
    alpha <- sort(unique(floor(z)))
    
    # Estimate pivot points.
    xcoord <- (min(z) - 1.0)
    ycoord <- 0.0
    for(a in alpha) {
        o <- (a + 1.0)
        elem <- which((z >= a) & (z < o))
        xcoord <- c(xcoord, mean(z[elem]))
        ycoord <- c(ycoord, mean(fdr[elem]))
    }
    
    # Interpolate back to the original resolution.
    points <- approx(xcoord, ycoord, xout=z)
    return(points$y)
}
#
# Convert an edge dataset into indexed graph representation.
# The input is a data frame with three columns TAIL, HEAD and WEIGHT
# for the edge end-points.
#
# Return value:
#   res$nodes      - N-element array of node names
#   res$tails      - K-element array of node indices
#   res$heads      - K-element array of node indices
#   res$weights    - K-element array of edge weights
#   res$tail2edge  - N-element list of adjacent edge indices
#   res$head2edge  - N-element list of adjacent edge indices
#   res$outstats   - N-row data frame of node statistics
#   res$instats    - N-row data frame of node statistics
#   res$stats      - N-row data frame of node statistics
#
# Written by Ville-Petteri Makinen 2013
#
tool.graph <- function(edges) {
    tails <- as.character(edges$TAIL)
    heads <- as.character(edges$HEAD)
    wdata <- as.double(edges$WEIGHT)
    
    # Remove empty end-points and non-positive weights.
    mask <- which((tails != "") & (heads != "") &
    (tails != heads) & (wdata > 0))
    tails <- tails[mask]
    heads <- heads[mask]
    wdata <- wdata[mask]
    nedges <- length(mask)
    
    # Create factorized representation.
    labels <- as.character(c(tails, heads))
    labels <- as.factor(labels)
    labelsT <- as.integer(labels[1:nedges])
    labelsH <- as.integer(labels[(nedges+1):(2*nedges)])
    
    # Create edge lists.
    nodnames <- levels(labels)
    nnodes <- length(nodnames)
    elistT <- tool.graph.list(labelsT, nnodes)
    elistH <- tool.graph.list(labelsH, nnodes)
    
    # Collect results.
    res <- list()
    res$nodes <- as.character(nodnames)
    res$outstats <- tool.graph.degree(elistT, wdata)
    res$instats <- tool.graph.degree(elistH, wdata)
    res$stats <- (res$outstats + res$instats)
    res$tail2edge <- elistT
    res$head2edge <- elistH
    res$tails <- as.integer(labelsT)
    res$heads <- as.integer(labelsH)
    res$weights <- wdata
    return(res)
}

#---------------------------------------------------------------------------

tool.graph.list <- function(entries, nnodes) {
    
    # Allocate list.
    groups <- list()
    for(i in 1:nnodes)
      groups[[i]] <- integer()
    
    # Find entry groups.
    st <- tool.aggregate(entries)
    labels <- as.integer(st$labels)
    blocks <- st$blocks
    
    # Reorganize according to node indices.
    nblocks <- length(blocks)
    for(k in 1:nblocks) {
        ind <- labels[[k]]
        groups[[ind]] <- blocks[[k]]
    }
    
    # Return edge lists.
    return(groups)
}

#---------------------------------------------------------------------------

tool.graph.degree <- function(node2edge, weights) {
    nnodes <- length(node2edge)
    stren <- rep(0.0, nnodes)
    degrees <- rep(0, nnodes)
    for(i in 1:nnodes) {
        rows <- node2edge[[i]]
        degrees[i] <- length(rows)
        if(degrees[i] < 1) next
        stren[i] <- sum(weights[rows])
    }
    res <- data.frame(DEGREE=degrees, STRENG=stren)
    return(res)
}
#
#
#
tool.metap <- function(datasets, idcolumn, pcolumn, weights=NULL) {
    
    # Collect identities.
    id <- character()
    nsets <- length(datasets)
    for(i in 1:nsets) {
        dat <- datasets[[i]]
        id <- as.character(dat[,idcolumn])
    }
    
    # Remove duplicates.
    id <- unique(id[which(id != "")])
    
    # Check weights.
    if(is.null(weights)) weights <- rep(1, nsets)
    if(length(weights) != nsets) stop("Incompatible weights.")
    if(sum(weights <= 0.0) > 0) stop("Non-positive weights.")
    weights <- nsets*(weights/sum(weights))
    
    # Calculate meta scores.
    h <- rep(0, length(id))
    z <- rep(0.0, length(id))
    for(i in 1:nsets) {
        dat <- datasets[[i]]
        pos <- match(dat[,idcolumn], id)
        rows <- which(pos > 0)
        pos <- pos[rows]
        p <- as.double(dat[rows,pcolumn])
        p <- pmax(p, .Machine$double.xmin)
        p <- pmin(p, (1 - .Machine$double.eps))
        z[pos] <- (z[pos] + weights[i]*qnorm(p))
        h[pos] <- (h[pos] + (weights[i])^2)
    }
    
    # Estimate meta P-values.
    rows <- which(h > 0.0)
    z[rows] <- z[rows]/sqrt(h[rows])
    pmeta <- pnorm(z)
    
    # Return results.
    res <- data.frame(ID=id, P=pmeta)
    names(res) <- c(idcolumn, pcolumn)
    res <- res[order(res$P),]
    return(res)
}
#
# One-sided.
#
tool.normalize <- function(x, prm=NULL, inverse=FALSE) {
    
    # Apply transform.
    if(is.null(prm) == FALSE) {
        if(inverse == TRUE) {
            if(is.na(prm$mu)) {
                y <- (x*(prm$scale) + prm$offset)
                return(y)
            }
            z <- (x*(prm$sigma) + prm$mu)
            z <- (exp(z) - 1.0)/(prm$scale)
            z <- (z + prm$offset)
            return(z)
        } else {
            z <- (x - prm$offset)*(prm$scale)
            if(is.na(prm$mu)) return(z)
            mask <- which(z < 0.0)
            z[mask] <- z[mask]/(1.0 - z[mask])
            z <- (log(z + 1.0) - prm$mu)/(prm$sigma)
        }
        return(z)
    }
    
    # Remove unusable values.
    x <- x[which(0*x == 0)]
    
    # Default parameters.
    prm$offset <- mean(x)
    prm$scale <- sd(x)
    prm$mu <- NA
    prm$sigma <- NA
    prm$quality <- 0.0
    if(length(x) < 10) return(prm)
    
    # Disable warnings.
    prev <- getOption("warn")
    options(warn=-1)
    
    # Ensure positivity and check size.
    xmin <- min(x)
    z <- (x[which(x > xmin)] - xmin)
    if(length(z) < 10) return(prm)
    
    # Scale by median.
    zmed <- median(z)
    z <- z/zmed
    
    # Find the best log transform.
    ctrl <- list(reltol=1e-3)
    gamma <- optim(par=1.0, fn=tool.normalize.quality, gr=NULL, z,
    lower=-9, upper=9, control=ctrl)
    
    # Apply transform.
    z <- log(exp(gamma$par)*z + 1.0)
    
    # Evaluate fit quality.
    mu <- mean(z)
    sigma <- sd(z)
    z <- (z - mu)/sigma
    kappa <- ks.test(x=z, y="pnorm", exact=FALSE)
    
    # Enable warnings.
    options(warn=prev)
    
    # Collect transformation parameters.
    prm$offset <- xmin
    prm$scale <- exp(gamma$par)/zmed
    prm$quality <- kappa$p.value
    prm$mu <- mu
    prm$sigma <- sigma
    return(prm)
}

#----------------------------------------------------------------------------

tool.normalize.quality <- function(g, z) {
    t <- log(exp(g)*z + 1.0)
    t <- (t - mean(t))/(sd(t) + 1e-20)
    if(length(t) < 1) return(0)
    suppressWarnings(res <- ks.test(x=t, y="pnorm", exact=FALSE))
    return(res$statistic)
}
#
# Calculate overlaps between groups of items.
#
# Input:
#   items     array of item identities
#   groups    array of group identities for items
#
#  Optional input:
#   nbackground   total number of items
#
#  Output:
#   res       data frame:
#             A      group name
#             B      group name
#             POSa   group name rank
#             POSb   group name rank
#             Na     group A size
#             Nb     group B size
#             Nab    shared items
#             R      overlap ratio
#             F      fold change to null expectation
#             P      overlap P-value (Fisher's test)
#
# Written by Ville-Petteri Makinen 2013
#
tool.overlap <- function(items, groups, nbackground=NULL) {
    
    # Check arguments.
    if(length(items) != length(groups))
      stop("Incompatible inputs.")
    
    # Remove duplicate entries.
    data <- data.frame(ITEM=items, GROUP=groups, stringAsFactors=FALSE)
    data <- unique(data)
    
    # Convert to integers.
    items <- as.integer(as.factor(data$ITEM))
    groups <- as.factor(data$GROUP)
    grouplevs <- levels(groups)
    groups <- as.integer(groups)
    
    # Determine block structure.
    modules <- tool.aggregate(groups)
    labels <- modules$labels
    blocks <- modules$blocks
    
    # Convert labels to original identities.
    labels <- as.integer(labels)
    labels <- grouplevs[labels]
    
    # Default background size.
    nitems <- length(unique(items))
    if(is.null(nbackground)) nbackground <- nitems
    if(nbackground < nitems)
      stop("tool.overlap: Invalid background size.")
    
    # Number of shared items for each pair of blocks.
    row <- 1
    stamp <- Sys.time()
    nblocks <- length(blocks)
    nrows <- (nblocks + 1)*nblocks/2
    mtx <- matrix(nrow=nrows, ncol=5)
    for(a in 1:nblocks) {
        maskA <- blocks[[a]]
        setA <- items[maskA]
        nA <- length(maskA)
        for(b in a:nblocks) {
            maskB <- blocks[[b]]
            setB <- items[maskB]
            nB <- length(setB)
            ind <- intersect(setA, setB)
            nAB <- length(ind)
            mtx[row,] <- c(a, b, nA, nB, nAB)
            row <- (row + 1)
            
            # Progress report.
            if(as.double(Sys.time() - stamp) < 10.0) next
            cat(sprintf("\r%d/%d ", row, nrows))
            stamp <- Sys.time()
        }
    }
    cat(sprintf("\r%d comparisons\n", nrows))
    
    # Calculate fold enrichment.
    numA <- mtx[,3]
    numB <- mtx[,4]
    numAB <- mtx[,5]
    fAB <- (sqrt(numAB)*sqrt(nbackground)/sqrt(numA)/sqrt(numB))^2
    
    # Estimate statistical significance.
    pAB <- phyper(numAB, numB, (nbackground - numB), numA, lower.tail=FALSE)
    pAB <- (pAB + .Machine$double.xmin)
    
    # Fix diagonal values.
    mask <- which(mtx[,1] == mtx[,2])
    fAB[mask] <- 1.0
    pAB[mask] <- 1.0
    
    # Fix P-values for small fold changes.
    mask <- which(fAB < 1.0)
    pAB[mask] <- 1.0
    
    # Collect results.
    res <- data.frame(A=labels[mtx[,1]], B=labels[mtx[,2]],
    POSa=as.integer(mtx[,1]), POSb=as.integer(mtx[,2]),
    Na=numA, Nb=numB, Nab=numAB,
    R=numAB/(numA + numB - numAB),
    F=fAB, P=pAB, stringsAsFactors=FALSE)
    return(res)
}
#
# Read a data frame from a file. All lines with NAs are excluded.
# The second argument selects a subset of columns.
#
# Written by Ville-Petteri Makinen 2013
#
tool.read <- function(file, vars=NULL) {
    if(is.null(file)) return(data.frame())
    if(file == "") return(data.frame())
    dat <- read.delim(file=file, header=TRUE,
    na.strings=c("NA", "NULL", "null", ""),
    colClasses="character", comment.char="",
    stringsAsFactors=FALSE)
    if(is.null(vars) == FALSE) dat <- dat[,vars]
    dat <- na.omit(dat)
    return(dat)
}
#
# Save a data frame in tab-delimited file. Compression
# only works on UNIX-family systems with gzip.
#
# Written by Ville-Petteri Makinen 2013
#
tool.save <- function(frame, file, directory=NULL, verbose=TRUE,
compression=FALSE) {
    if(verbose) cat("\rWriting to file... ")
    
    # Concatenate directory prefix.
    fname <- file
    if(is.null(directory) == FALSE) {
        if(file.exists(directory) == FALSE)
          dir.create(path=directory, recursive=TRUE)
        fname <- file.path(directory, file)
    }
    
    # Write data to file.
    write.table(x=frame, sep="\t", file=fname, na="",
    row.names=FALSE, quote=FALSE)
    
    # Compress file.
    if(compression) {
      if(verbose) cat("\rCompressing file... ")
      system(paste("gzip -f \"", fname, "\"", sep=""))
      fname <- paste(fname, ".gz", sep="")
    }
    
    # Print report.
    n <- nrow(frame)
    if(verbose) cat("\rSaved ", n, " rows in '", fname, "'.\n", sep="")
    return(fname)  
}
#
# Determine the network neighbors for a set of nodes. The first
# argument is the indexed graph structure from tool.graph().
# The second argument is the list of seed node names and the third
# indicates the maximum number of links to connect neighbors.
# The fourth input sets the directionality: use a negative value
# for downstream, positive for upstream or zero for undirected.
#
# Output:
#   res$RANK    - indices of neighboring nodes (including seeds)
#   res$LEVEL   - num of edges away from seed 
#   res$STRENG  - sum of adjacent edge weights within neighborhood
#   res$DEGREE  - number of adjacent edges within neighborhood
#
# Written by Ville-Petteri Makinen 2013
#
tool.subgraph <- function(graph, seeds, depth=1, direction=0) {
    depth <- as.integer(depth[[1]])
    direction <- as.integer(direction[[1]])
    
    # Convert seed names to indices.
    nodes <- graph$nodes
    ranks <- match(seeds, nodes)
    ranks <- ranks[which(ranks > 0)]
    ranks <- as.integer(ranks)
    
    # Find neighbors.
    res <- tool.subgraph.search(graph, ranks, depth[1], direction)
    return(res)
}

#---------------------------------------------------------------------------

tool.subgraph.search <- function(graph, seeds, depth, direction) {
    tails <- graph$tails
    heads <- graph$heads
    weights <- graph$weights
    tail2edge <- list()
    head2edge <- list()
    
    # Downstream and upstream searches.
    if(direction <= 0) tail2edge <- graph$tail2edge
    if(direction >= 0) head2edge <- graph$head2edge
    
    # Collect neighboring nodes.
    visited <- seeds
    levels <- 0*seeds
    for(i in 1:depth) {
        
        # Find edges to adjacent nodes.
        foundT <- tool.subgraph.find(seeds, tail2edge, heads, visited)
        foundH <- tool.subgraph.find(seeds, head2edge, tails, visited)
        
        # Expand neighborhood.
        seeds <- unique(c(foundT, foundH))
        visited <- c(visited, seeds)
        levels <- c(levels, (0*seeds + i))
        if(length(seeds) < 1) break
    }
    
    # Calculate node degrees and strengths.
    res <- data.frame(RANK=visited, LEVEL=levels, DEGREE=0,
    STRENG=0.0, stringsAsFactors=FALSE)
    res <- tool.subgraph.stats(res, tail2edge, heads, weights)
    res <- tool.subgraph.stats(res, head2edge, tails, weights)
    return(res)
}

#---------------------------------------------------------------------------

tool.subgraph.find <- function(seeds, edgemap, heads, visited) {
    if(length(edgemap) < 1) return(integer())
    
    # Collect all neighboring nodes.
    nneigh <- length(seeds)
    neighbors <- seeds
    for(i in seeds) {
        
        # Edges adjacent to the seed.
        mask <- edgemap[[i]]
        nmask <- length(mask)
        if(nmask < 1) next
        
        # Check capacity.
        if((length(neighbors) - nneigh) < nmask)
          neighbors <- c(neighbors, 0*neighbors, 0*mask)
        
        # Store node indices.
        neighbors[nneigh+(1:nmask)] <- heads[mask]
        nneigh <- (nneigh + nmask)
    }
    
    # Remove nodes already visited.
    neighbors <- unique(neighbors[1:nneigh])
    neighbors <- setdiff(neighbors, visited)
    return(neighbors)
}

#---------------------------------------------------------------------------

tool.subgraph.stats <- function(frame, edgemap, heads, weights) {
    if(length(edgemap) < 1) return(frame)
    nmemb <- nrow(frame)
    wcounts <- frame$DEGREE
    wsums <- frame$STRENG
    members <- frame$RANK
    for(i in members) {
        edges <- edgemap[[i]]
        pos <- match(heads[edges], members)
        edges <- edges[which(pos > 0)]
        if(length(edges) < 1) next
        rows <- match(heads[edges], members)
        wcounts[rows] <- (wcounts[rows] + 1)
        wsums[rows] <- (wsums[rows] + weights[edges])
    }
    frame$DEGREE <- wcounts
    # if weights are all the same, return flat weights
    # instead of weights determined by node degree
    if(sd(weights)==0) wsums <- rep(1, length(wcounts))
    frame$STRENG <- wsums
    return(frame)
}
#
# Translate words.
#
# Written by Ville-Petteri Makinen 2013
#
tool.translate <- function(words, from, to) {
    
    # Check translation table.
    if(length(from) != length(to)) stop("Incompatible inputs.")
    rows <- which((is.na(from) == FALSE) & (is.na(to) == FALSE))
    from <- as.character(from[rows])
    to <- as.character(to[rows])
    
    # Find words that can be translated.
    pos <- match(words, from)
    words[which(is.na(pos))] <- NA
    
    # Translate words.
    rows <- which(pos > 0)
    pos <- pos[rows]
    words[rows] <- to[pos]
    return(words)
}
#
# Convert a distribution to uniform ranks ]0,1[ with respect
# to a background distribution (or self if no
# background available).
#
# Written by Ville-Petteri Makinen 2013
#
tool.unify <- function(xtrait, xnull=NULL) {
    if(is.null(xnull)) xnull <- xtrait
    nt <- length(xtrait)
    n0 <- length(xnull)
    
    # Enforce distinct values.
    xmin <- min(xnull)
    sigma <- (max(xnull) - xmin)
    amp <- 0.001*sd(xnull)
    jitter <- seq(-amp/2, amp/2, length.out=n0)
    xnull <- (xnull + jitter)
    xnull <- (xnull - min(xnull))
    xnull <- xnull/max(xnull)
    xnull <- (xmin + xnull*sigma)
    
    # Cumulative reference distribution.
    f0 <- seq(0, 1, length.out=n0)
    q0 <- quantile(xnull, f0)
    q0 <- as.double(q0)
    
    # Make sure null range contains trait values.
    xmin <- min(xtrait)
    xmax <- max(xtrait)  
    if(xmin < q0[1]) q0[1] <- xmin
    if(xmax > q0[n0]) q0[n0] <- xmax
    
    # Map observations on null cumulative axis.
    out <- approx(x=q0, y=f0, xout=xtrait)
    return(out$y)
}
