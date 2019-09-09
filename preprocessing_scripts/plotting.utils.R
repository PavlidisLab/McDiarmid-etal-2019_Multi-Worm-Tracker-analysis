DEFAULT_FONT =   theme(text=element_text(family="Arial", face="bold")) #Times New Roman, 12pt, Bold
GWB_COLORS = config::get("colors")$`gwb-scheme`
HM_COLORS = colorRampPalette(colors = config::get("colors")$`heatmap-scheme`)(7)



allele.labels = read.csv(
  file = config::get("strain-genotype-labels"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)
feature.labels = read.csv(
  file = config::get('feature-subcategories-labels'),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)

feature.class.colors=NULL
feature.class = NULL
feature.class.groups = NULL
feature.class.colormap = NULL

loadHeatmapFeatureColors <-
  ## Load heatmap feature colors scheme into ***global*** environment.
  ## i.e, pay attention to the <<- operators.
  function(X) {
    feature.class <<-
      read.csv(
        config::get("feature-subcategories"),
        sep = "\t",
        stringsAsFactors = F
      )
    feature.class.groups <<-
      as.character(unlist(lapply(colnames(X), function(x)
        feature.class[feature.class[, 1] == x, 2])))
    feature.class.colormap <<- list(
      "Baseline Locomotion" = config::get("colors")$`baseline-locomotion`,
      "Morphology" = config::get("colors")$morphology,
      "Spontaneous Recovery/Learning Retention" = config::get("colors")$`spontaneous-recovery_learning-retention`,
      "Initial Sensitivity" = config::get("colors")$`initial-sensitivity`,
      "Habituation Learning" = config::get("colors")$`habituation-learning`
    )
    
    feature.class.colors <<- unique(as.character(unlist(
      lapply(feature.class.groups, function(x)
        as.character(feature.class.colormap[names(feature.class.colormap) == x]))
    )))
    
  }

iheatmapper <- function(X,
                        colorClip = NULL,
                        useGroups = NULL,
                        titleName = "Multiworm Tracker Features",
                        valuesNames = "Values",
                        ...) {
  # Generate a heatmap
  X <- as.matrix(X)
  loadHeatmapFeatureColors(X)
  
  colnames(X) <- translate(colnames(X), dictionary = "features")
  rownames(X) <- translate(rownames(X), dictionary = "alleles")
  
  if (is.null(colorClip)) {
    p <- main_heatmap(
      X,
      name = valuesNames,
      x_categorical = TRUE,
      y_categorical = TRUE,
      layout = list(margin = list(b = 120)),
      colors = HM_COLORS,
      
      ...
    )
  } else{
    p <- main_heatmap(
      X,
      name = valuesNames,
      x_categorical = TRUE,
      y_categorical = TRUE,
      colors = HM_COLORS,
      layout = list(margin = list(b = 120)),
      zmin = colorClip[1],
      zmid = colorClip[2],
      zmax = colorClip[3],
      ...
    )
  }
  
  
  p <- p %>%
    add_col_groups(
      groups = feature.class.groups,
      side = "bottom",
      name = "Features",
      title = "Type",
      colors = rev(feature.class.colors)[c(5, 2, 3, 4, 1)] # This is a kind of weird workaround, how do you reorder to color-factor pairings? Doesn't seem supported. #(as.vector(sapply(feature.class.colors, function(x) strtrim(x = adjustcolor(x), width = 7)))) # feature.class.colors
    ) %>%
    add_col_title(titleName,
                  side = "top",
                  font = list(size = 24, family = "Arial"))  %>%
    add_col_labels(
      ticktext = paste0("<b>", colnames(X), "</b>"),
      font = list(size = 10, family = "Arial")
    ) %>%
    add_row_labels(
      ticktext = paste0("<b>", row.names(X), "</b>"),
      size = 0.3,
      font = list(size = 10, family = "Arial")
    )
  
  if (is.null(useGroups)) {
    p <- p %>%
      add_row_clustering(method = "hclust",
                         side = "right",
                         name = "Variant clustering") %>%
      add_col_clustering(method = "hclust",
                         name = "Feature clustering",
                         side = "top")
  } else{
    useGroups <-
      useGroups[order(useGroups$group, as.character(useGroups$gene)), ]
    class <- setColorClass(X, clusters = useGroups)
    
    colourScheme <- getPreDefColoursForClass(class)
    
    names(colourScheme) <- class
    unique.colors.from.scheme <-
      unique(as.character(colourScheme[order(names(colourScheme))]))
    
    featureGrpVector <- c()
    for (col in colnames(X)) {
      if (col %in% colnames(reversals)) {
        featureGrpVector <- c(featureGrpVector, "reversal")
      } else {
        featureGrpVector <- c(featureGrpVector, "morphologic")
      }
    }
    
    p <- p %>%
      add_row_clustering(
        method = "groups",
        side = "right",
        groups = names(colourScheme),
        colors = unique.colors.from.scheme,
        name = "Variant clustering"
      ) %>%
      add_col_clustering(
        method = "hclust",
        name = "Feature type",
        groups = featureGrpVector,
        side = "top",
        show_colorbar = F
      )
    
  }
  
  return(p)
}

translate <-
  ## Translate a vector of values x using a specified dictionary (alleles, or features.)
  function(x, dictionary = "alleles", on_failure=c("fail", "usekey")) {
    
    on_failure = match.arg(on_failure)
    
    if (dictionary == "alleles") {
      dictionary = allele.labels
    } else if (dictionary == "features") {
      dictionary = feature.labels
    }
    
    ret =  as.character(unlist(sapply(
      X = c(x),
      FUN = function(key) {
        translation = dictionary[dictionary$Source == key, "Label"]
        if (on_failure == "fail") {
          translation = ifelse(is.na(translation), stop("Unknown translation"), translation)
        } else if (on_failure == "usekey") {
          translation = ifelse(is.na(translation), key, translation) 
        } else {
          stop("Unexpected behaviour in translate()" )  
        }
        
        translation
      }
    )))

    if (!assertthat::are_equal(length(c(x)),
                               length(ret))) {
      stop("Translation failed.")
    }

    return(ret)
  }


pvplot <-
  # Plot and save a pvclust object
  # X: a dataframe
  # breaks: Vector of alpha cutoff(s) for cluster boxes
  # filename: Path to save file
  # HIGHLIGHT : A vector of groups to highlight in colour.
  # HIGHLIGHT_v : A vector palette of colours
  function(X,
           breaks = c(0.95),
           path = NULL,
           filename = NULL,
           HIGHLIGHT = NULL,
           HIGHLIGHT_v = NULL) {
    
    if (is.null(path)){
      path=config::get("clustering")$plots  
    }
    
    if (is.null(HIGHLIGHT)) {
      # Don't colour highlight labels.
      HIGHLIGHT = c(0)
    } else if (is.null(HIGHLIGHT_v)) {
      # Highlight certain labels.
      HIGHLIGHT_v <-
        HIGHLIGHT$variable[match(
          gsub(
            pattern = " ",
            replacement = "",
            X$hclust$labels[X$hclust$order]
          ),
          gsub(
            pattern = " ",
            replacement = "",
            days.melted$day
          )
        )]
      for (color in unique(HIGHLIGHT_v)) {
        HIGHLIGHT_v[HIGHLIGHT_v == color]  <-
          brewer.pal(n = length(HIGHLIGHT_v), name = 'Paired')[color]
      }
      HIGHLIGHT_V <-
        brewer.pal(n = length(HIGHLIGHT_v), name = 'Paired')
    }
    
    
    # print("my breaks are")
    # print(breaks)
    colors = c('red',
               'blue',
               'green',
               'pink')
    while (!is.null(dev.list())) {
      dev.off()
    }
    
    if (!is.null(filename)) {
      png(
        filename = paste0(path, filename),
        width = 1920,
        height = 1080,
        units = "px"
      )
    }
    par(mar = c(10, 1, 1, 12))
    X %>% as.dendrogram %>%
      #set("branches_k_color", k = 2, value = c("purple", "orange")) %>%
      set("labels_col", HIGHLIGHT_v) %>%
      set("labels_cex", 0.75) %>%
      plot
    X %>% text
    
    
    for (i in 1:length(breaks)) {
      #p %>% pvrect(alpha=breaks[i], pv="au", type="geq", max.only=MAX.ONLY, border=colors[i])
      X %>% pvrect(
        alpha = breaks[i],
        pv = "au",
        type = "geq",
        max.only = MAX.ONLY,
        border = colors[i]
      )
    }
    
    if (!is.null(filename)) {
      # Print to device for fun.
      dev.off()
      pvplot(
        X = X,
        breaks = breaks,
        filename = NULL,
        HIGHLIGHT_v = HIGHLIGHT_v
      )
    }
    return()
  }


phenotypePlotByMean <- 
  
  # A point and error bar plot for the t.test results
  function(X=NULL,
           confidence.hi=NULL,
           confidence.lo=NULL,
           pvalues=NULL,
           metric=NULL,
           useClusters=NULL,
           showLabels=T,
           limitOutlier=NULL,
           seed = NA,
           repulsion = 30) {
    
    X<-data.frame(X)
    # x = cbind(confidence.hi, N2=numeric(0))
    isWT <- grep(x=rownames(X), pattern = "^N2") 
    
    X.WT.mean <- mean(X[isWT,c(metric)])
    X.WT.CI.hi.hi <- max(confidence.hi[grepl("^N2.", x = rownames(confidence.hi)),metric])
    X.WT.CI.hi.lo <- min(confidence.hi[grepl("^N2.", x = rownames(confidence.hi)),metric])
    
    X.WT.CI.lo.hi <- max(confidence.lo[grepl("^N2.", x = rownames(confidence.lo)),metric])
    X.WT.CI.lo.lo <- min(confidence.lo[grepl("^N2.", x = rownames(confidence.lo)),metric])
    # Remove wts and replace with mean
    X[isWT, metric] = X.WT.mean
    # confidence.hi[isWT, metric] = 0.0
    # confidence.lo[isWT, metric] = 0.0
    # pvalues[isWT, metric] = 0.0
    
    X <- X[!grepl("^N2.", x = rownames(X)),]
    confidence.hi <- confidence.hi[!grepl("^N2.", x = rownames(confidence.hi)),]
    confidence.lo <- confidence.lo[!grepl("^N2.", x = rownames(confidence.lo)),]
    
    confidence.hi[grepl("^N2$", x = rownames(confidence.hi)),metric] = 0
    confidence.lo[grepl("^N2$", x = rownames(confidence.lo)),metric] = 0
    
    isWT <- grep(x=rownames(X), pattern = "^N2") # Update 
    
    class <- setColorClass(X, clusters = useClusters)
    colourScheme <- getPreDefColoursForClass(class)
    names(colourScheme) <- class
    
    outlier <- translate(rownames(X)) # Convert labels first
    
    outlier[  ( X[,metric] >= X.WT.mean & confidence.lo[,metric] <= X.WT.mean ) |
                ( X[,metric] <= X.WT.mean & confidence.hi[,metric] >= X.WT.mean ) ] <- ""
    
    # class[  rownames(X) != "N2" & (( X[,metric] >= 0 & confidence.lo[,metric] <= 0 ) |
    #             ( X[,metric] <= 0 & confidence.hi[,metric] >= 0 )) ] <- "ASD (n.s.)"
    # colourScheme[  rownames(X) != "N2" & (( X[,metric] >= 0 & confidence.lo[,metric] <= 0 ) |
    #                                       ( X[,metric] <= 0 & confidence.hi[,metric] >= 0 )) ] <- "gray39"
    # names(colourScheme)[  rownames(X) != "N2" & (( X[,metric] >= 0 & confidence.lo[,metric] <= 0 ) |
    #                                                       ( X[,metric] <= 0 & confidence.hi[,metric] >= 0 )) ] <- "ASD (n.s.)"
    
    class[  rownames(X) != "N2" & (( X[,metric] >= X.WT.mean & confidence.lo[,metric] <= X.WT.mean ) |
                                     ( X[,metric] <= X.WT.mean & confidence.hi[,metric] >= X.WT.mean )) ] <- "ASD (n.s.)"
    colourScheme[  rownames(X) != "N2" & (( X[,metric] >= X.WT.mean & confidence.lo[,metric] <= X.WT.mean ) |
                                            ( X[,metric] <= X.WT.mean & confidence.hi[,metric] >= X.WT.mean )) ] <- "gray39"
    names(colourScheme)[  rownames(X) != "N2" & (( X[,metric] >= X.WT.mean & confidence.lo[,metric] <= X.WT.mean ) |
                                                   ( X[,metric] <= X.WT.mean & confidence.hi[,metric] >= X.WT.mean )) ] <- "ASD (n.s.)"
    
    ORDER = rank( X[,c(metric)], ties.method = "first")
    
    p <- ggplot(X, 
                aes(x=ORDER, 
                    y=X[,c(metric)]
                ))  
    
    if (!showLabels){
      ## WT hbar
      p = p +
        geom_hline(yintercept = X.WT.CI.hi.hi,
                   linetype = "solid") +
        geom_hline(yintercept = X.WT.CI.hi.lo,
                   linetype = "dotted") +
        geom_hline(yintercept = X.WT.CI.lo.hi,
                   linetype = "dotted") +
        geom_hline(yintercept = X.WT.CI.lo.lo,
                   linetype = "solid")
      
    }
    
    p = p + 
      geom_point(size = DEFAULT_PT_SIZE/1.8, 
                 mapping = aes(colour=class)
      ) + 
      theme_classic() +
      theme(text=element_text(family="Arial")) +
      # geom_hline(aes(yintercept=Z3up),  size=DEFAULT_LN_SIZE, linetype='longdash') +
      # geom_hline(aes(yintercept=Z3down),  size=DEFAULT_LN_SIZE, linetype='longdash') +
      TITLE_SIZES +
      labs(title=paste(paste(translate(metric, dictionary = "features"), "by Genotype")), 
           subtitle= "", 
           x="Genotypes", 
           y = paste("Sample mean distance from wild-type" ) ) + 
      scale_fill_manual(name="Group", values = colourScheme, guide = F)  +
      scale_color_manual(name="Group", values = colourScheme)
    
    #    if (!is.null(confidence)){
    if (T){
      p <- p +
        geom_errorbar(aes(ymin=confidence.lo[,c(metric)], ymax=confidence.hi[,c(metric)], colour=class), 
                      width=DEFAULT_LN_WIDTH,
                      size=DEFAULT_LN_SIZE,
                      position=position_dodge(.9)
        )
    }
    
    if (showLabels && length(outlier) != sum(as.numeric(outlier %in% ""))) { # Are there any outliers to show?
      if (!is.null(limitOutlier)){
        
        LABEL.CANDIDATES = outlier[order(X[,c(metric)])] # Labels in correct order
        WT.POS = grep(pattern = "^N2", row.names(X[order(X[,c(metric)]), ])) # Labels in correct order
        
        # ISLABEL = LABEL.CANDIDATES[grep(TRUE, (LABEL.CANDIDATES[1:WT.POS] != ""))]
        LABEL.LEFT = LABEL.CANDIDATES[grep(TRUE, (LABEL.CANDIDATES[1:WT.POS] != ""))][1:limitOutlier]
        LABEL.RIGHT = LABEL.CANDIDATES[WT.POS:length(LABEL.CANDIDATES)][grep(TRUE, (LABEL.CANDIDATES[WT.POS:length(LABEL.CANDIDATES)] != ""))]
        if (length(LABEL.RIGHT) > limitOutlier) {
          LABEL.RIGHT = LABEL.RIGHT[(length(LABEL.RIGHT) - limitOutlier + 1):(length(LABEL.RIGHT))]    
        }
        LABEL.SAVED = unique(c(LABEL.LEFT, LABEL.RIGHT))
        LABEL.PRUNE = setdiff(outlier, c(LABEL.SAVED, ""))
        outlier[outlier %in% LABEL.PRUNE] <- ""
      }
      textcolours <- class
      textcolours[class == "Wildtype"] <- "white"
      textcolours[class != "Wildtype"] <- "black"
      
      # NUDGER = max(c( confidence.hi[,c(metric)], confidence.lo[,c(metric)])) - quantile( c( confidence.hi[,c(metric)], confidence.lo[,c(metric)]), probs = 0.85)
      NUDGER = quantile( confidence.hi[,c(metric)] - X[,metric], probs = 0.80)
      p <- p +
        geom_label_repel(aes(ORDER,
                             X[,c(metric)],
                             label = outlier,
                             fill = class
                             
        ),
        color = 'white', seed = seed,
        size = 7,
        force = repulsion,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.5, "lines"),
        segment.color = 'black',
        colour = textcolours, 
        fontface = "bold",
        max.iter = 3e3,
        nudge_y = ifelse(ORDER %% 2 == 0, -1 * NUDGER, NUDGER)
        )  
      
    }
    
    p = p +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    
    return(p)
    
  }


setColorClass <- function(X, clusters = NULL){ # TODO: Default parameter here is useful for testing but dangerous for production.
  class <- rownames(X)
  WTs <- grepl(pattern="^N2$|N2-", x=class)
  
  class.orig <- class
  class[T] <- "Other"
  
  for (i in 1:length(class.orig)) {
    gene = clusters[i,2]
    group = clusters[i,1]
    class[class.orig %in% gene] <- group 
  }
  
  class[WTs] <- "Wildtype"
  class[!WTs] <- "ASD"
  
  return( class ) 
}



getPreDefColoursForClass <- 
  # Alternative method to get colours for larger numbers of colours provided in the palette.
  function(classOfX){
    classColours <- classOfX
    colourCount <- length(unique(classOfX))
    #getPalette <- viridis_pal(option = "C")(colourCount) # TODO: Remove me
    # getPalette <- select_col_vector[1:colourCount]
    # names(getPalette) <- levels(classOfX)
    for ( i in 1:colourCount ){
      k <- unique(classOfX)[i]
      if ( k == 'Wildtype' ){
        classColours[classOfX == k] <- 'black'
      }else{
        classColours[classOfX == k] <- GWB_COLORS[length(GWB_COLORS)]
      }
      
    }
    return(classColours)
  }


getColoursForClass <- function(X){
  colourCount <- length(unique(X))
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  names(getPalette) <- levels(X)
  for ( i in 1:length(unique(X)) ){
    k <- unique(X)[i]
    if (is.na(k) | k=='Wildtype'){
      X[X == k] <- 'black'
    }else{
      X[X == k] <- getPalette(colourCount)[i]  
    }
    
  }
  return(X)
}


correlationHeatmap <-
  # pheatmap based correlation heatmaps
  function(X, title, zmin=-1.0, zmax=1.0, border_color = NA, ...){
    
    loadHeatmapFeatureColors(X)
    feature.class.df = data.frame(Feature=  feature.class$Class) # Set globally.
    rownames(feature.class.df) <- feature.class$Feature
    p <- 
      pheatmap(X,
               labels_row = translate(rownames(X), dictionary = "features"),
               labels_col = translate(colnames(X), dictionary = "features"), 
               color = colorRampPalette(HM_COLORS)(20),
               breaks = seq(zmin, zmax, length = 21),
               na_col = 'black',
               annotation_col = feature.class.df,
               annotation_colors = list( Feature = unlist(feature.class.colormap) ),
               fontsize = 10,
               border_color = border_color,
               cellwidth = 9, 
               cellheight = 9,
               ...
      )
    return(p)
  }


tstatBarPerGene <- 
  # A bar plot for log p-values from the t.test results
  function(t.stat=NULL,
           strain=NULL,
           class=NULL,
           THRESHOLD=0.05) {
    
    
    reordering = (as.character(sapply(colnames(t.stat), function(x) feature.class$Class[feature.class$Feature == x][1] )))
    reordering[reordering == "Morphology"] <- "Baseline2"
    t.stat.tmp = t.stat[,order(reordering)]
    
    X<-melt(t(data.frame(t.stat.tmp)[strain,]))
    names(X) <- c("Metric", "Strain", "Value")
    
    grups = as.character(sapply(as.character(X$Metric), function(x) feature.class$Class[feature.class$Feature == x][1] ))
    kolors = as.character(sapply(as.character(sapply(as.character(X$Metric), function(x) feature.class$Class[feature.class$Feature == x][1] )), function(y) feature.class.colormap[y]))
    
    X$Groups = grups
    X$Colors = kolors
    
    class = feature.class.groups
    class.colors = class
    colourScheme = feature.class.colors
    colourScheme = names(class.colors)
    
    MINY <- quantile(t.stat, probs = 0.01) #min(t.stat)
    MAXY <- quantile(t.stat, probs = 0.99) #max(t.stat)
    X$Value <- fence(X$Value, UB=MAXY, LB=MINY)
    
    X = X[order(X$Groups),]
    X$Metric <- mapvalues(X$Metric, from=as.character(X$Metric), translate(as.character(X$Metric), dictionary = "features"))
    # levels(X$Metric) <- translate(as.character(X$Metric), dictionary = "features")
    
    p <- ggplot(X, 
                aes(
                  x = Metric,
                  y = Value ,
                  fill = X$Groups
                ) 
    ) + 
      theme_classic() + 
      geom_bar(stat="identity", position="dodge", width=0.75 ) +
      scale_y_continuous( limits = c(MINY,MAXY) ) +
      scale_fill_manual(name="Feature type", labels = unique(X$Groups), values = unique(X$Colors) ) +
      # scale_x_discrete(labels=translate(as.character(X$Metric), dictionary = "features"))+
      theme(axis.text.x = element_text(angle=60, hjust=1, size=10,face = "bold", colour = 'black')) +
      ylab(label = "T-Statistic") +
      ggtitle(paste0("Phenotypic profile for ", translate(strain, dictionary = "alleles") )) + 
      TITLE_SIZES + 
      DEFAULT_FONT
    
    return(p)
    
  }


multiplot <-
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  
function(...,
         plotlist = NULL,
         file,
         rows = 1,
         layout = NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, rows * ceiling(numPlots / rows)),
                     nrow = rows,
                     ncol = ceiling(numPlots / rows))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]],
            vp = viewport(
              layout.pos.row = matchidx$col, ## Swapped the rows and cols so they render left to right.
              layout.pos.col = matchidx$row
            ))
    }
  }
}


