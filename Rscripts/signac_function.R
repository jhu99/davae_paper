SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

ApplyMatrixByGroup <- function(
  mat,
  groups,
  fun,
  normalize = FALSE,
  group.scale.factors = NULL,
  scale.factor = NULL
) {
  if (normalize) {
    if (is.null(x = group.scale.factors) | is.null(x = scale.factor)) {
      stop("If normalizing counts, supply group scale factors")
    }
  }
  results <- list()
  all.groups <- unique(x = groups)
  for (i in seq_along(along.with = all.groups)) {
    pos.cells <- names(x = groups)[groups == all.groups[[i]]]
    if (length(x = pos.cells) > 1) {
      totals <- fun(x = mat[pos.cells, ])
    } else {
      totals <- mat[pos.cells, ]
    }
    results[[i]] <- data.frame(
      group = all.groups[[i]],
      count = totals,
      position = as.numeric(colnames(x = mat)),
      stringsAsFactors = FALSE
    )
  }
  coverages <- as.data.frame(
    x = do.call(what = rbind, args = results), stringsAsFactors = FALSE
  )
  if (normalize) {
    scale.factor <- SetIfNull(
      x = scale.factor, y = median(x = group.scale.factors)
    )
    coverages$norm.value <- coverages$count /
      group.scale.factors[coverages$group] * scale.factor
  } else {
    coverages$norm.value <- coverages$count
  }
  return(coverages)
}

GetGroups <- function(
  object,
  group.by,
  idents
) {
  if (is.null(x = group.by)) {
    obj.groups <- Idents(object = object)
  } else {
    obj.md <- object[[group.by]]
    obj.groups <- obj.md[, 1]
    names(obj.groups) <- rownames(x = obj.md)
  }
  if (!is.null(idents)) {
    obj.groups <- obj.groups[obj.groups %in% idents]
  }
  return(obj.groups)
}


SingleCoveragePlot <- function(
  object,
  region,
  annotation = NULL,
  peaks = NULL,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  window = 100,
  downsample = 0.1,
  height.tracks = 10,
  extend.upstream = 0,
  extend.downstream = 0,
  ymax = NULL,
  scale.factor = NULL,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-")
) {
  cells <- SetIfNull(x = cells, y = colnames(x = object))
  if (!is.null(x = idents)) {
    ident.cells <- WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
  }
  if (!is(object = region, class2 = 'GRanges')) {
    region <- StringToGRanges(regions = region, sep = sep)
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  reads.per.group <- AverageCounts(
    object = object,
    group.by = group.by,
    verbose = FALSE
  )
  cells.per.group <- CellsPerGroup(
    object = object,
    group.by = group.by
  )
  cutmat <- CutMatrix(
    object = object,
    region = region,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  group.scale.factors <- reads.per.group * cells.per.group
  scale.factor <- SetIfNull(
    x = scale.factor, y = median(x = group.scale.factors)
  )
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  coverages <- ApplyMatrixByGroup(
    mat = cutmat,
    fun = colSums,
    groups = obj.groups,
    group.scale.factors = group.scale.factors,
    scale.factor = scale.factor,
    normalize = FALSE
  )
  if (!is.na(x = window)) {
    coverages <- group_by(.data = coverages, group)
    coverages <- mutate(.data = coverages, coverage = rollapply(
      data = norm.value,
      width = window,
      FUN = mean,
      align = 'center',
      fill = NA
    ))
    coverages <- ungroup(x = coverages)
  } else {
    coverages$coverage <- coverages$norm.value
  }
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  stepsize <- 1 / downsample
  retain_positions <- seq(from = start.pos, to = end.pos, by = stepsize)
  downsampled_coverage <- coverages[coverages$position %in% retain_positions, ]
  ymax <- SetIfNull(x = ymax, y = signif(
    x = max(downsampled_coverage$coverage, na.rm = TRUE), digits = 2)
  )
  ymin <- 0
  downsampled_coverage <- downsampled_coverage[!is.na(
    x = downsampled_coverage$coverage
  ), ]
  
  gr <- GRanges(
    seqnames = chromosome,
    IRanges(start = start.pos, end = end.pos)
  )
  p <- ggplot(
    data = downsampled_coverage,
    mapping = aes(x = position, y = coverage, fill = group)
  ) +
    geom_area(stat = 'identity') +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
    xlab(label = paste0(chromosome, ' position (bp)')) +
    ylab(label = paste0(as.character(x = ymin), ' - ',
                        as.character(x = ymax), ')')) +
    ylim(c(ymin, ymax)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      legend.position = 'none',
     strip.text.y = element_blank(),
    )
  if (!is.null(x = peaks)) {
    # subset to covered range
    peak.intersect <- subsetByOverlaps(x = peaks, ranges = gr)
    peak.df <- as.data.frame(x = peak.intersect)
    peak.plot <- ggplot(data = peak.df, mapping = aes(color = 'darkgrey')) +
      geom_segment(aes(x = start, y = 0, xend = end, yend = 0, size = 2),
                   data = peak.df) +
      theme_classic() +
      ylab(label = "Peaks") +
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            legend.position = 'none') +
      xlab(label = paste0(chromosome, ' position (bp)')) +
      xlim(c(start.pos, end.pos)) +
      scale_color_identity()
    # remove axis from coverage plot
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank()
    )
  } else {
    peak.plot <- NULL
  }
  if (!is.null(x = annotation)) {
    if (!inherits(x = annotation, what = 'GRanges')) {
      stop("Annotation must be a GRanges object or EnsDb object.")
    }
    annotation.subset <- subsetByOverlaps(x = annotation, ranges = gr)
    annotation.df <- as.data.frame(x = annotation.subset)
    # adjust coordinates so within the plot
    annotation.df$start[annotation.df$start < start.pos] <- start.pos
    annotation.df$end[annotation.df$end > end.pos] <- end.pos
    annotation.df$direction <- ifelse(
      test = annotation.df$strand == "-", yes = -1, no = 1
    )
    if (nrow(x = annotation.df) > 0) {
      gene.plot <- ggplot(
        data = annotation.df,
        mapping = aes(
          xmin = start,
          xmax = end,
          y = strand,
          fill = strand,
          label = gene_name,
          forward = direction)
      ) +
        geom_gene_arrow(
          arrow_body_height = unit(x = 4, units = "mm"),
          arrowhead_height = unit(x = 4, units = "mm"),
          arrowhead_width = unit(x = 5, units = "mm"),
          forward = TRUE) +
        geom_gene_label(
          grow = TRUE,
          reflow = TRUE,
          height = unit(x = 4, units = "mm")
        ) +
        xlim(start.pos, end.pos) +
        xlab(label = paste0(chromosome, ' position (bp)')) +
        ylab("Genes") +
        theme_classic() +
        theme(legend.position = 'none',
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank())
      
      # remove axis from coverage plot
      p <- p + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank()
      )
      if (!is.null(x = peak.plot)) {
        peak.plot <- peak.plot + theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank()
        )
        p <- p +
          peak.plot +
          gene.plot +
          plot_layout(ncol = 1, heights = c(height.tracks, 1, 1))
      } else {
        p <- p +
          gene.plot +
          plot_layout(ncol = 1, heights = c(height.tracks, 1))
      }
    } else {
      if (!is.null(peak.plot)) {
        p <- p +
          peak.plot +
          plot_layout(ncol = 1, heights = c(height.tracks, 1))
      }
    }
  } else {
    if (!is.null(peak.plot)) {
      p <- p +
        peak.plot +
        plot_layout(ncol = 1, heights = c(height.tracks, 1))
    }
  }
  return(p)
}


CoveragePlot <- function(
  object,
  region,
  annotation = NULL,
  peaks = NULL,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  window = 100,
  downsample = 0.1,
  height.tracks = 10,
  extend.upstream = 0,
  extend.downstream = 0,
  scale.factor = NULL,
  ymax = NULL,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-"),
  ...
) {
  if (length(x = region) > 1) {
    plot.list <- lapply(
      X = seq_along(region),
      FUN = function(x) {
        SingleCoveragePlot(
          object = object,
          region = region[x],
          annotation = annotation,
          peaks = peaks,
          assay = assay,
          fragment.path = fragment.path,
          group.by = group.by,
          window = window,
          downsample = downsample,
          ymax = ymax,
          scale.factor = scale.factor,
          extend.upstream = extend.upstream,
          extend.downstream = extend.downstream,
          cells = cells,
          idents = idents,
          sep = sep
        )
      }
    )
    return(wrap_plots(plot.list, ...))
  } else {
    return(SingleCoveragePlot(
      object = object,
      region = region,
      annotation = annotation,
      peaks = peaks,
      assay = assay,
      fragment.path = fragment.path,
      group.by = group.by,
      window = window,
      downsample = downsample,
      height.tracks = height.tracks,
      extend.upstream = extend.upstream,
      extend.downstream = extend.downstream,
      ymax = ymax,
      scale.factor = scale.factor,
      cells = cells,
      idents = idents,
      sep = sep
    ))
  }
}