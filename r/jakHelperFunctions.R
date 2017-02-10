deseq.fisher = function(sigtab,FDR = 0.05) {

  library("dplyr")

  sigPos = dim(filter(sigtab,padj<=FDR,log2FoldChange >0))[1]
  notPos = dim(filter(sigtab,padj>FDR,log2FoldChange >0))[1]
  sigNeg = dim(filter(sigtab,padj<=FDR,log2FoldChange <0))[1]
  notNeg = dim(filter(sigtab,padj>FDR,log2FoldChange <0))[1]

  mat = matrix(c(sigNeg,sigPos,notNeg,notPos), 2,2)
  colnames(mat) = c("Sig","Not")
  rownames(mat) = c("Neg","Pos")
  print(mat)
  fisher.test(mat)$p.value
}

##

pqPlot = function(pvalues,Q=0.1,m = "BH",b = 0.01) {

  # adjust p-value to q-values
  df = data.frame(pvalues)
  colnames(df) <- c("p")
  df$q =p.adjust(df$p,method = m)

  qTemp = filter(df,q<=Q)
  pCut = max(qTemp$p)
  df$y = "red"
  df$y[df$p < pCut] = "blue"
  df$y[df$p > Q] = "black"

  pTemp = filter(df,p<=Q)
  qCut = max(pTemp$q)
  df$x = "red"
  df$x[df$q < Q] = "blue"
  df$x[df$q > qCut] = "black"

  points = ggplot(df, aes(x=p,y=q)) +
    geom_point(color=ifelse(df$p<Q,ifelse(df$q<Q,"blue","red"),"black")) +
    geom_hline(yintercept = Q) +
    geom_vline(xintercept = Q)

  # marginal density of x - plot on top
  plot_top <- ggplot(df, aes(p)) +
    geom_histogram(binwidth = b,aes(fill=y)) +
    theme(legend.position = "none") +
    #geom_vline(xintercept = c(Q,pCut)) +
    scale_fill_identity()

  # marginal density of y - plot on the right
  plot_right <- ggplot(df, aes(q)) +
    geom_histogram(binwidth = b,aes(fill=x)) +
    coord_flip() +
    theme(legend.position = "none") +
    #geom_vline(xintercept = c(Q,qCut)) +
    scale_fill_identity()

  empty <- ggplot() + geom_point(aes(1, 1), colour = "white") + theme(plot.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank())

  grid.arrange(plot_top, empty, points, plot_right, ncol = 2, nrow = 2, widths = c(4,1), heights = c(1, 4))
}
