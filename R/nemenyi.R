nemenyi <- function (data, conf.level = 0.95, sort = c(TRUE, FALSE), 
          plottype = c("vline", "none", "mcb", "vmcb", "line", "matrix"), 
          select = NULL, labels = NULL, title = NULL, ...){
  sort <- sort[1]
  plottype <- match.arg(plottype, c("vline", "none", "mcb", 
                                    "vmcb", "line", "matrix"))
  if (length(dim(data)) != 2) {
    stop("Data must be organised as methods in columns and observations in rows.")
  }
  data <- as.matrix(data)
  data <- na.exclude(data)
  rows.number <- nrow(data)
  cols.number <- ncol(data)
  if (!is.null(select) && (select > cols.number)) {
    select <- NULL
  }
  if (plottype != "none") {
    sort <- TRUE
  }
  if (is.null(labels)) {
    labels <- colnames(data)
    if (is.null(labels)) {
      labels <- 1:cols.number
    }
  }
  else {
    labels <- labels[1:cols.number]
  }
  fried.pval <- stats::friedman.test(data)$p.value
  if (fried.pval <= 1 - conf.level) {
    fried.H <- "Ha: Different"
  }
  else {
    fried.H <- "H0: Identical"
  }
  r.stat <- stats::qtukey(conf.level, cols.number, Inf) * sqrt((cols.number * 
                                                                  (cols.number + 1))/(12 * rows.number))
  ranks.matrix <- t(apply(data, 1, function(x) {
    rank(x, na.last = "keep", ties.method = "average")
  }))
  ranks.means <- colMeans(ranks.matrix)
  ranks.intervals <- rbind(ranks.means - r.stat, ranks.means + 
                             r.stat)
  if (sort == TRUE) {
    order.idx <- order(ranks.means)
  }
  else {
    order.idx <- 1:cols.number
  }
  ranks.means <- ranks.means[order.idx]
  ranks.intervals <- ranks.intervals[, order.idx]
  labels <- labels[order.idx]
  if (!is.null(select)) {
    select <- which(order.idx == select)
  }
  if (plottype != "none") {
    args <- list(...)
    args.nms <- names(args)
    if (!is.null(title)){
      args$main <- title
    } else {
      if (!("main" %in% args.nms)) {
        args$main <- paste0("Friedman: ", format(round(fried.pval, 
                                                       3), nsmall = 3), " (", fried.H, ") \n Critical distance: ", 
                            format(round(r.stat, 3), nsmall = 3), sep = "")
      }
    }
    if (!("xaxs" %in% names(args))) {
      args$xaxs <- "i"
    }
    if (!("yaxs" %in% names(args))) {
      args$yaxs <- "i"
    }
    nc <- max(nchar(labels))
    nc <- nc/1.75 + 1
    nr <- nchar(sprintf("%1.2f", round(max(ranks.means), 
                                       2)))/1.75
    parmar.def <- parmar <- graphics::par()$mar
  }
  if ((plottype == "mcb") | (plottype == "vmcb")) {
    cmp <- RColorBrewer::brewer.pal(3, "Set1")[1:2]
    if (fried.pval > 1 - conf.level) {
      pcol <- "gray"
    }
    else {
      pcol <- cmp[2]
    }
    mnmx <- range(ranks.means) + c(-0.5, 0.5) * r.stat
    mnmx <- mnmx + diff(mnmx) * 0.04 * c(-1, 1)
    if (plottype == "mcb") {
      if (!("xlab" %in% names(args))) {
        args$xlab <- ""
      }
      if (!("ylab" %in% names(args))) {
        args$ylab <- "Mean ranks"
      }
      if (is.null(args$xlim)) {
        args$xlim <- c(0, cols.number + 1)
      }
      if (is.null(args$ylim)) {
        args$ylim <- mnmx
      }
    }
    else {
      if (!("ylab" %in% names(args))) {
        args$ylab <- ""
      }
      if (!("xlab" %in% names(args))) {
        args$xlab <- "Mean ranks"
      }
      if (is.null(args$ylim)) {
        args$ylim <- c(0, cols.number + 1)
      }
      if (is.null(args$xlim)) {
        args$xlim <- mnmx
      }
    }
    args$x <- args$y <- NA
    args$axes <- FALSE
    if ((plottype == "mcb") && (parmar[1] < (nc + nr))) {
      parmar[1] <- nc + nr
    }
    if ((plottype == "vmcb") && (parmar[2] < (nc + nr))) {
      parmar[2] <- nc + nr
      parmar[3] <- 1.5 # 
    }
    par(mar = parmar)
    if (is.null(select)) {
      select <- 1
    }
    do.call(plot, args)
    if (plottype == "mcb") {
      polygon(c(0, rep(cols.number + 1, 2), 0), rep(ranks.means[select], 
                                                    4) + r.stat/2 * c(1, 1, -1, -1), col = "gray90", 
              border = NA)
      points(1:cols.number, ranks.means, pch = 20, lwd = 3)
      axis(1, at = c(1:cols.number), labels = paste0(labels, 
                                                     " - ", sprintf("%1.2f", round(ranks.means, 2))), 
           las = 2)
      axis(2)
      for (i in 1:cols.number) {
        lines(rep(i, times = 2), ranks.means[i] + c(-1, 
                                                    1) * 0.5 * r.stat, type = "o", lwd = 1, col = pcol, 
              pch = 20)
      }
      idx <- abs(ranks.means[select] - ranks.means) < r.stat
      points((1:cols.number)[idx], ranks.means[idx], pch = 20, 
             lwd = 3, col = cmp[1])
    }
    else {
      polygon(rep(ranks.means[select], 4) + r.stat/2 * 
                c(1, 1, -1, -1), c(0, rep(cols.number + 1, 2), 
                                   0), col = "gray90", border = NA)
      points(ranks.means, 1:cols.number, pch = 20, lwd = 3)
      axis(2, at = c(1:cols.number), labels = paste0(labels, 
                                                     " - ", sprintf("%1.2f", round(ranks.means, 2))), 
           las = 2)
      axis(1)
      for (i in 1:cols.number) {
        lines(ranks.means[i] + c(-1, 1) * 0.5 * r.stat, 
              rep(i, times = 2), type = "o", lwd = 1, col = pcol, 
              pch = 20)
      }
      idx <- abs(ranks.means[select] - ranks.means) < r.stat
      points(ranks.means[idx], (1:cols.number)[idx], pch = 20, 
             lwd = 3, col = cmp[1])
    }
    box(which = "plot", col = "black")
  }
  if (plottype == "matrix") {
    rline <- array(NA, c(cols.number, 2))
    nem.mat <- array(0, c(cols.number, cols.number))
    for (i in 1:cols.number) {
      rline[i, ] <- c(which(ranks.means > ranks.intervals[1, 
                                                          i])[1], tail(which(ranks.means < ranks.intervals[2, 
                                                                                                           i]), 1))
      nem.mat[i, rline[i, 1]:rline[i, 2]] <- 1
    }
    diag(nem.mat) <- 2
    if (parmar[1] < nc) {
      parmar[1] <- nc
    }
    nr <- nchar(sprintf("%1.2f", round(max(ranks.means), 
                                       2)))/1.75
    if (parmar[2] < (nc + nr)) {
      parmar[2] <- nc + nr
    }
    par(mar = parmar)
    if (!("xlab" %in% names(args))) {
      args$xlab <- ""
    }
    if (!("ylab" %in% names(args))) {
      args$ylab <- ""
    }
    args$x <- 1:cols.number
    args$y <- 1:cols.number
    args$z <- nem.mat[order(order.idx), cols.number:1]
    args$axes <- FALSE
    cmp <- c("white", RColorBrewer::brewer.pal(3, "Set1")[2], 
             "black")
    if (fried.pval > 1 - conf.level) {
      cmp[3] <- "gray60"
      cmp[2] <- "gray70"
    }
    args$col <- cmp
    do.call(image, args)
    for (i in 1:cols.number) {
      abline(v = i + 0.5)
      abline(h = i + 0.5)
    }
    axis(1, at = 1:cols.number, labels = labels[order(order.idx)], 
         las = 2)
    axis(2, at = 1:cols.number, labels = paste0(rev(labels), 
                                                " - ", sprintf("%1.2f", round(rev(ranks.means), 2))), 
         las = 2)
    box()
    if (!is.null(select)) {
      polygon(order.idx[select] + c(-0.5, 0.5, 0.5, -0.5), 
              which((nem.mat[order(order.idx), cols.number:1])[order.idx[select], 
              ] == 2) + c(-0.5, -0.5, 0.5, 0.5), col = "red")
    }
  }
  if ((plottype == "line") | (plottype == "vline")) {
    rline <- matrix(NA, nrow = cols.number, ncol = 2)
    for (i in 1:cols.number) {
      tloc <- which((abs(ranks.means - ranks.means[i]) < 
                       r.stat) == TRUE)
      rline[i, ] <- c(min(tloc), max(tloc))
    }
    rline <- unique(rline)
    rline <- rline[apply(rline, 1, min) != apply(rline, 1, 
                                                 max), ]
    if (length(rline) == 2) {
      rline <- as.matrix(rline)
      rline <- t(rline)
    }
    k <- nrow(rline)
    cmp <- colorRampPalette(RColorBrewer::brewer.pal(12, 
                                                     "Paired"))(k)
    if (fried.pval > 1 - conf.level) {
      cmp <- rep("gray", times = k)
    }
    lbl <- paste0(labels, " - ", sprintf("%1.2f", round(ranks.means, 
                                                        2)))
    if (!("ylab" %in% names(args))) {
      args$ylab <- ""
    }
    if (!("xlab" %in% names(args))) {
      args$xlab <- ""
    }
    args$x <- args$y <- NA
    args$axes <- FALSE
    if (plottype == "line") {
      if (is.null(args$xlim)) {
        args$xlim <- c(1, cols.number)
      }
      if (is.null(args$ylim)) {
        args$ylim <- c(0, k + 1)
      }
    }
    else {
      if (is.null(args$xlim)) {
        args$xlim <- c(0, k + 1)
      }
      if (is.null(args$ylim)) {
        args$ylim <- c(1, cols.number)
      }
    }
    if ((plottype == "line") && (parmar[1] < (nc + nr))) {
      parmar[1] <- nc + nr
    }
    if ((plottype == "vline") && (parmar[2] < (nc + nr))) {
      parmar[2] <- nc + nr
    }
    par(mar = parmar)
    do.call(plot, args)
    if (plottype == "line") {
      points(1:cols.number, rep(0, cols.number), pch = 20, 
             lwd = 4)
      if (k > 0) {
        for (i in 1:k) {
          lines(rline[i, ], c(i, i), col = cmp[i], lwd = 4)
          lines(rep(rline[i, 1], times = 2), c(0, i), 
                col = "gray", lty = 2)
          lines(rep(rline[i, 2], times = 2), c(0, i), 
                col = "gray", lty = 2)
        }
      }
      axis(1, at = c(1:cols.number), labels = lbl, las = 2)
      if (!is.null(select)) {
        points(select, 0, pch = 20, col = RColorBrewer::brewer.pal(3, 
                                                                   "Set1")[1], cex = 2)
      }
    }
    else {
      points(rep(0, cols.number), 1:cols.number, pch = 20, 
             lwd = 4)
      if (k > 0) {
        for (i in 1:k) {
          lines(c(i, i), rline[i, ], col = cmp[i], lwd = 4)
          lines(c(0, i), rep(rline[i, 1], times = 2), 
                col = "gray", lty = 2)
          lines(c(0, i), rep(rline[i, 2], times = 2), 
                col = "gray", lty = 2)
        }
      }
      axis(2, at = c(1:cols.number), labels = lbl, las = 2)
      if (!is.null(select)) {
        points(0, select, pch = 20, col = RColorBrewer::brewer.pal(3, 
                                                                   "Set1")[1], cex = 2)
      }
    }
  }
  if (plottype != "none") {
    par(mar = parmar.def)
  }
  return(structure(list(means = ranks.means, intervals = ranks.intervals, 
                        fpval = fried.pval, fH = fried.H, cd = r.stat, conf.level = conf.level, 
                        k = cols.number, n = rows.number), class = "nemenyi"))
}