deseq.fisher = function(sigtab, FDR = 0.05) {

  library("dplyr")

  sigPos = dim(filter(sigtab, padj <= FDR, log2FoldChange > 0))[1]
  notPos = dim(filter(sigtab, padj >  FDR, log2FoldChange > 0))[1]
  sigNeg = dim(filter(sigtab, padj <= FDR, log2FoldChange < 0))[1]
  notNeg = dim(filter(sigtab, padj >  FDR, log2FoldChange < 0))[1]

  mat = matrix(c(sigNeg, sigPos, notNeg, notPos), 2, 2)
  colnames(mat) = c("Sig", "Not")
  rownames(mat) = c("Neg", "Pos")
  print(mat)
  fisher.test(mat)$p.value
}

##

pqPlot = function(pvalues, Q = 0.1, m = "BH", b = 0.01) {

  # adjust p-value to q-values
  df = data.frame(pvalues)
  colnames(df) <- c("p")
  df$q = p.adjust(df$p, method = m)

  qTemp = filter(df,q <= Q)
  pCut = max(qTemp$p)
  df$y = "red"
  df$y[df$p < pCut] = "blue"
  df$y[df$p > Q] = "black"

  pTemp = filter(df, p <= Q)
  qCut = max(pTemp$q)
  df$x = "red"
  df$x[df$q < Q] = "blue"
  df$x[df$q > qCut] = "black"

  points = ggplot(df, aes(x = p, y = q)) +
    geom_point(color = ifelse(df$p < Q, ifelse(df$q < Q, "blue", "red"), "black")) +
    geom_hline(yintercept = Q) +
    geom_vline(xintercept = Q)

  # marginal density of x - plot on top
  plot_top <- ggplot(df, aes(p)) +
    geom_histogram(binwidth = b, aes(fill = y)) +
    theme(legend.position = "none") +
    #geom_vline(xintercept = c(Q,pCut)) +
    scale_fill_identity()

  # marginal density of y - plot on the right
  plot_right <- ggplot(df, aes(q)) +
    geom_histogram(binwidth = b, aes(fill = x)) +
    coord_flip() +
    theme(legend.position = "none") +
    #geom_vline(xintercept = c(Q,qCut)) +
    scale_fill_identity()

  empty <- ggplot() + geom_point(aes(1, 1), colour = "white") + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank())

  grid.arrange(plot_top, empty, points, plot_right, ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
}

permanova = function(phyloObject, group = "TYPE", metric = "jsd", palettePick = FALSE) {

  phyloObject = zeroOTU(phyloObject)
  sd = data.frame(sample_data(phyloObject))
  sd = sd[group]
  colnames(sd) = "GROUP"

  palette = ""
  if (palettePick[1] == FALSE) {
    colorCount = length(unique(sd$GROUP))
    palette = getPalette(colorCount)
  } else{
    palette = palettePick
  }

  distJSD = phyloseq::distance(phyloObject, metric)
  adon.GROUP = adonis(distJSD ~ GROUP, sd)
  GROUP.p = adon.GROUP$aov.tab$`Pr(>F)`[1]

  jsd.pco = dudi.pco(cailliez(distJSD), scannf=F, nf=2)
  df = cbind(jsd.pco$li, sd)
  df = select(df, A1, A2, GROUP)
  df.x = aggregate(A1 ~ GROUP, data = df, mean)
  colnames(df.x) = c("GROUP", "X0")
  df.y = aggregate(A2 ~ GROUP, data = df, mean)
  colnames(df.y) = c("GROUP", "Y0")
  df = left_join(df, df.x, by = "GROUP")
  df = left_join(df, df.y, by = "GROUP")

  p1 = ggplot(df, aes(x = A1, y = A2)) +
    geom_segment(aes(x = X0, y = Y0, xend = A1, yend = A2, color = GROUP)) +
    geom_point(pch = 21, color = "black", size = 5, aes(fill = GROUP)) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette)

  list("plot" = p1,
       "pval" = GROUP.p,
       "group" = group,
       "metric" = metric)
}

makeNetworkPhyloseq = function(p, FDR = 0.05) {

  library("igraph")
  library("Hmisc")

  matrix = data.frame(otu_table(p),check.names = F)
  matrix = t(matrix)

  sequenceCor = cor(matrix)
  sequencePval = rcorr(as.matrix(matrix))$P

  sequenceCor.df = melt(sequenceCor)
  sequencePval.df = melt(sequencePval)

  colnames(sequenceCor.df) <- c("X", "Y", "value")
  colnames(sequencePval.df) <- c("X", "Y", "p")
  sequenceDF = merge(sequenceCor.df, sequencePval.df)
  sequenceDF = filter(sequenceDF, X != Y) # remove self correlations
  sequenceDF$BH = p.adjust(sequenceDF$p, method = "BH")
  correlation = filter(sequenceDF, BH <= FDR)

  el = select(correlation, X, Y, value)
  # el = filter(el, value >= 0.9) # Pearson Correlation Cutoff (or `value <= -0.2`)
  el[,1]=as.character(el[,1])
  el[,2]=as.character(el[,2])
  el=as.matrix(el)

  graph = graph.edgelist(el[,1:2], directed = F)
  E(graph)$weight=as.numeric(el[,3])
  graph = simplify(graph, edge.attr.comb = list("mean"))

  # save layout
  l <- layout.fruchterman.reingold(graph)
  V(graph)$x <- l[,1]
  V(graph)$y <- l[,2]

  E(graph)$color <- ifelse(E(graph)$weight >= 0, "black", "red")

  returnObject = list("graph" = graph, "cor" = correlation, "triangles" = tri)
  returnObject
}
