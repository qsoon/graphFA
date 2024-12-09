# localization operator
localization.op <- function(evalues, evectors, g, i=NULL){ # support of T_i g is centered at node i
  if(is.null(i)){
    res <- evectors %*% (t(Conj(evectors))*g(evalues)) # ith column: T_i g
  } else{
    res <- as.vector(evectors %*% (Conj(evectors)[i,]*g(evalues))) # T_i g vector when i is specified
  }
  return(res)
}

# graph windowed Fourier transform
graph.window.FT <- function(x, S, g, M){
  # eigenres <- eigen(S)
  # evalues <- eigenres$values
  # evectors <- eigenres$vectors
  # lmax <- max(evalues)
  val <- eigensort(S)
  evalues <- val$evalues
  evectors <- val$evectors
  lmax <- max(evalues)
  
  C <- NULL
  normconst <- c()
  for(m in 1:M){
    gm <- function(lambda){
      return(g(lambda, sigma.sq=lmax*(M+1)/M^2, m, tau=lmax*(M+1)/M^2))
    }
    Tgm <- localization.op(evalues, evectors, gm)
    C <- cbind(C, t(Conj(Tgm))%*%x) # where C[i,m] = <x, Ti gm>
    normconst <- c(normconst, norm(Tgm, type="F")^2)
  }
  return(list(C=C, normconst=normconst))
}

# cpsd estimate
cpsd.graph <- function(x, y, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  if(is.null(method)){
    # x, y should be R realizations (N x R matrix)
    # eigenres <- eigen(S)
    # x.tilde <- t(Conj(eigenres$vectors))%*%x
    # y.tilde <- t(Conj(eigenres$vectors))%*%y
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    x.tilde <- t(Conj(evectors))%*%x
    y.tilde <- t(Conj(evectors))%*%y
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  } else if(method=="window"){
    # x, y should be one realization (N x 1 vector)
    C1 <- graph.window.FT(x, S, g, M)
    C2 <- graph.window.FT(y, S, g, M)
    cpsd <- colSums(C1$C * Conj(C2$C)) / C1$normconst
  } else if(method=="random"){
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    WB <- windowbank.random(N=length(evalues), M=M, V=evectors, sigma=sigma, seed=seed)
    x.tilde <- t(Conj(evectors))%*%(t(WB)*as.vector(x))
    y.tilde <- t(Conj(evectors))%*%(t(WB)*as.vector(y))
    cpsd <- rowMeans(x.tilde * Conj(y.tilde))
  }
  return(cpsd)
}

# psd estimate
psd.graph <- function(x, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  return(cpsd.graph(x=x, y=x, S=S, g=g, M=M, sigma=sigma, method=method, seed=seed))
}

# coherence estimate
coherence.graph <- function(x, y, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  if(is.null(method)){
    cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M, seed=seed)
    psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M, seed=seed)
    psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M, seed=seed)
  } else if(method=="window"){
    cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M, method=method, seed=seed)
    psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M, method=method, seed=seed)
    psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M, method=method, seed=seed)
  } else if(method=="random"){
    cpsd <- cpsd.graph(x=x, y=y, S=S, M=M, sigma=sigma, method=method, seed=seed)
    psd.x <- cpsd.graph(x=x, y=x, S=S, M=M, sigma=sigma, method=method, seed=seed)
    psd.y <- cpsd.graph(x=y, y=y, S=S, M=M, sigma=sigma, method=method, seed=seed)
  }
  return(cpsd*Conj(cpsd) / psd.x / psd.y)
}

# cross spectrum analysis All in one
cross.spectrum.graph <- function(x, y, S, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  cpsd <- cpsd.graph(x=x, y=y, S=S, g=g, M=M, sigma=sigma, method, seed=seed)
  psd.x <- cpsd.graph(x=x, y=x, S=S, g=g, M=M, sigma=sigma, method, seed=seed)
  psd.y <- cpsd.graph(x=y, y=y, S=S, g=g, M=M, sigma=sigma, method, seed=seed)
  return(list(cpsd=cpsd, psd.x = psd.x, psd.y = psd.y, coherence=cpsd*Conj(cpsd) / psd.x / psd.y))
}

graph.stationary.level <- function(cov, S){
  val <- eigensort(S)
  evalues <- val$evalues
  evectors <- val$evectors
  gamma <- t(Conj(evectors))%*%cov%*%evectors
  return(norm(diag(gamma), type=2) / norm(gamma, type="F"))
}

# robust cpsd estimate 
huberloss <- function(t, c){
  res <- c()
  for(i in t){
    if(abs(i) <= c){
      res <- c(res, i^2)
    } else{
      res <- c(res, 2*c*(abs(i)-c/2))
    }
  }
  return(res)
}

PIRLS <- function(lambda, y, tau, alpha, B, P, gamma_init, max_iterations, tol, check){
  n <- length(y)
  gamma <- gamma_init
  for (iteration in 1:max_iterations){
    gamma.old <- gamma
    if(check=="Zheng"){
      W <- diag(grad_checkft_by_Zheng(as.numeric(y-B%*%gamma), tau, alpha) / as.numeric((2*(y-B%*%gamma))))
    } else if(check=="Oh"){
      W <- diag(grad_checkft_by_Oh(as.numeric(y-B%*%gamma), tau, alpha) / as.numeric((2*(y-B%*%gamma))))
    }
    
    tmp <- solve(t(B)%*%W%*%B + n*lambda*P)%*%t(B)%*%W
    gamma <- tmp%*%y
    gammahat <- gamma
    if (norm(gamma-gamma.old, type="2") < tol) {
      cat("Iteration:", iteration, "\n")
      cat("Converged!\n")
      break
    }
  }
  
  res <- gammahat
  
  return(res)
}

grad_checkft_by_Oh <- Vectorize(function(x, tau, c){
  if((-c<=x) & (x<c)){
    return(tau*x/c*as.numeric(x >= 0) + (1-tau)*x/c*as.numeric(x < 0))
  }
  else{
    return(tau - as.numeric(x < -c))
  }
}, vectorize.args = "x")

# random window generation
windowbank.random <- function(N, M, V, sigma, seed=1){
  res <- matrix(0, nrow=M, ncol=N)
  set.seed(seed)
  for(i in 1:M){
    W.tilde <- diag(N) + rnorm(N^2, 0, sigma)
    res[i,] <- diag(V %*% W.tilde %*% t(Conj(V)))
  }
  return(res)
}


# cpsd estimate
# robust.cpsd.graph.window <- function(x, y, S, g, M){
#   C1 <- graph.window.FT(x, S, g, M)
#   C2 <- graph.window.FT(y, S, g, M)
#   cpsd <- colSums(C1$C * Conj(C2$C)) / C1$normconst
#   return(cpsd)
# }



plot_graph_custom <- function (z, size = 0.75, edge_color, vertex_color=NULL) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$color <- factor(edge_color, levels = unique(edge_color))
  p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                  xend = y1, yend = y2, colour = color), 
                                              linewidth = 2, data = df2) +
    scale_color_manual(values=color.cand, labels = paste("Line", 1:8, sep=""),
                       name = "Line number") + 
    geom_point(aes(fill=vertex_color), size = size, shape=21) + 
    scale_fill_gradient(low="white", high="black", na.value = "yellow", name = "People") +
    theme_void() +
    theme(legend.margin = margin(10,10,10,10))
  print(p1)
}

plot_graph_custom3 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2),  linewidth=e.size, color = "gray", data = df2) + ggtitle(title) +
      geom_point(size = v.size) + theme_void()+theme(aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}

plot_graph_custom4 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="gray", high="gray", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}

plot_graph_custom5 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.3, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}


plot_graph_custom_directed <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30, directed=TRUE) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(directed==FALSE){
    if(signal==FALSE){
      # p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
      #                                                 xend = y1, yend = y2, color=w), linewidth=e.size, 
      #                                             data = df2) +
      #   scale_color_gradient(low="gray", high="gray", name="Edge Weight")+
      #   geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      #   scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      #   # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      #   theme_void() +
      #   guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      #   theme(legend.margin = margin(mg), plot.margin = margin(mg), 
      #         plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
      #         legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
      # print(p1)
      p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = y1, y = y2, 
                                                      xend = x, yend = y),  linewidth=e.size, color = "gray", data = df2) + ggtitle(title) +
        geom_point(size = v.size) + theme_void()+theme(aspect.ratio=ratio)
      print(p1)
    }
    
    else{
      p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                      xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                  data = df2) +
        scale_color_gradient(low="grey", high="black", name="Edge Weight")+
        geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
        scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
        # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
        theme_void() +
        guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
        theme(legend.margin = margin(mg), plot.margin = margin(mg), 
              plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
              legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
      print(p1)
    }
  } else{
    if(signal==FALSE){
      # p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = y1, y = y2, 
      #                                                 xend = x, yend = y, color=w), linewidth=e.size, 
      #                                             data = df2) +
      #   scale_color_gradient(low="gray", high="gray", name="Edge Weight")+
      #   geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      #   scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      #   # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      #   theme_void() +
      #   guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      #   theme(legend.margin = margin(mg), plot.margin = margin(mg), 
      #         plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
      #         legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
      # print(p1)
      p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = y1, y = y2, 
                                                      xend = x, yend = y), arrow=arrow(), linewidth=e.size, color = "gray", data = df2) + ggtitle(title) +
        geom_point(size = v.size) + theme_void()+theme(aspect.ratio=ratio)
      print(p1)
    }
    
    else{
      p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = y1, y = y2, 
                                                      xend = x, yend = y, color=w), arrow=arrow(), linewidth=e.size, 
                                                  data = df2) +
        scale_color_gradient(low="grey", high="black", name="Edge Weight")+
        geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
        scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
        # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
        theme_void() +
        guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
        theme(legend.margin = margin(mg), plot.margin = margin(mg), 
              plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
              legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
      print(p1)
    }
  }
  
}

# plot_graph_custom_directed <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30, directed=FALSE) 
# {
#   if (is(z$sA, "sparseMatrix")) {
#     z$sA <- summary(z$sA)
#   }
#   x <- z$xy[, 1]
#   y <- z$xy[, 2]
#   # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
#   w <- rownames(z$xy)
#   ind_i <- z$sA[, 1]
#   ind_j <- z$sA[, 2]
#   y1 <- x[ind_j]
#   y2 <- y[ind_j]
#   df1 <- data.frame(x = x, y = y, w=w)
#   df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
#   df2$w <- z$sA[,3]
#   if(is.null(min)){
#     min <- min(vertex_color)
#   }
#   if(is.null(max)){
#     max <- max(vertex_color)
#   }
#   
#   if(directed==FALSE){
#     if(signal==FALSE){
#       p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
#                                                       xend = y1, yend = y2),  linewidth=e.size, color = "gray", data = df2) + ggtitle(title) +
#         geom_point(size = v.size) + theme_void()+theme(aspect.ratio=ratio)
#       print(p1)
#     }
#     
#     else{
#       p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
#                                                       xend = y1, yend = y2, color=w), linewidth=e.size, 
#                                                   data = df2) +
#         scale_color_gradient(low="grey", high="black", name="Edge Weight")+
#         geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
#         scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
#         # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
#         theme_void() +
#         guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
#         theme(legend.margin = margin(mg), plot.margin = margin(mg), 
#               plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
#               legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
#       print(p1)
#     }
#   } else {
#     if(signal==FALSE){
#       p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = y1, y = y2, 
#                                                       xend = x, yend = y), arrow=arrow(),  linewidth=e.size, color = "gray", data = df2) + ggtitle(title) +
#         geom_point(size = v.size) + theme_void()+theme(aspect.ratio=ratio)
#       print(p1)
#     }
#     
#     else{
#       p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
#                                                       xend = y1, yend = y2, color=w), linewidth=e.size, 
#                                                   data = df2) +
#         scale_color_gradient(low="grey", high="black", name="Edge Weight")+
#         geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
#         scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
#         # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
#         theme_void() +
#         guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
#         theme(legend.margin = margin(mg), plot.margin = margin(mg), 
#               plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
#               legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
#       print(p1)
#     }
#   }
#   
# }


GFPCA <- function(X, S, q=NULL, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  if(is.null(method)){
    # X[[i]] should be R realizations (n x R matrix)
    # p-dimensional graph signal
    p <- length(X)
    if(nrow(X[[1]])!=nrow(S)){
      stop("The number of rows of list X's element matrices should be equal to the number of rows of GSO!")
    }
    n <- nrow(S)
    R <- ncol(X[[1]])
    
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    # evalues.rev <- rev(evalues)
    # evectors.rev <- val$evectors[,n:1]
    
    P_array <- array(NA, dim = c(n,p,p))
    for(i in 1:p){
      for(j in 1:p){
        P_array[,i,j] <- cpsd.graph(x=X[[i]], y=X[[j]], S=S) 
      }
    }
    
    H.hat <- array(NA, dim = c(n,p,p))
    G.hat <- array(NA, dim = c(n,p,p))
    tau.hat <- matrix(NA, nrow=n, ncol=p)
    
    for(l in 1:n){ # \lambda_1 <= ... <= \lambda_n
      tmp <- eigensort(P_array[l,,])
      G.hat[l,,] <- tmp$evectors[,p:1]
      H.hat[l,,] <- t(Conj(G.hat[l,,]))
      tau.hat[l,] <- rev(tmp$evalues)
    }
    
    if(is.null(q)){
      q <- p
    }
    
    Y <- list()
    # Y[[i]] : dimension-reduced graph signal for R realizations (n x R matrix) / length q
    for(i in 1:q){
      Y[[i]] <- matrix(0, nrow=n, ncol=R)
      for(j in 1:p){
        Y[[i]] <- Y[[i]] + evectors %*% diag(H.hat[,i,j]) %*% t(Conj(evectors)) %*% X[[j]]
      }
    }
    mu.x.hat <- do.call(cbind, lapply(X, FUN=rowMeans)) # n x p matrix
    
    mu.y.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
    
    for(j in 1:p){
      for(k in 1:p){
        mu.y.hat[,j] <- mu.y.hat[,j] + evectors %*% diag(H.hat[,j,k]) %*% t(Conj(evectors)) %*% mu.x.hat[,k]
      }
    }
    
    mu.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
    X.hat <- list()
    # X.hat[[i]] : reconstructed graph signal for R realizations (n x R matrix) / length p
    for(i in 1:p){
      mu.hat[,i] <- mu.x.hat[,i] 
      X.hat[[i]] <- matrix(0, nrow=n, ncol=R)
      for(j in 1:q){
        mu.hat[,i] <- mu.hat[,i] - evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% mu.y.hat[,j]
        X.hat[[i]] <- X.hat[[i]] + evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% Y[[j]]
      }
      X.hat[[i]] <- X.hat[[i]] + mu.hat[,i]
    }  
    
    res <- list()
    res$X <- X ; res$Y <- Y ; res$X.hat <- X.hat
    res$mu.hat <- mu.hat ; res$tau.hat <- tau.hat ; res$H.hat <- H.hat ; res$G.hat <- G.hat ; res$Px <- P_array
    return(res)
  } 
  # else if(method=="window"){
  #     # X[[i]] should be one realization (n x 1 vector)
  #     # p-dimensional graph signal
  #     if(is.vector(X[[1]])){
  #       X <- lapply(X, as.matrix)
  #     }
  #     p <- length(X)
  #     if(nrow(X[[1]])!=nrow(S)){
  #       stop("The length of list X's element vector OR the number of rows of list X's element matrices should be equal to the number of rows of GSO!")
  #     }
  #     n <- nrow(S)
  #     R <- ncol(X[[1]])
  #     
  #     val <- eigensort(S)
  #     evalues <- val$evalues
  #     evectors <- val$evectors
  #     # evalues.rev <- rev(evalues)
  #     # evectors.rev <- val$evectors[,n:1]
  #     
  #     P_array <- array(NA, dim = c(n,p,p))
  #     for(i in 1:p){
  #       for(j in 1:p){
  #         P_array[,i,j] <- cpsd.graph(x=X[[i]], y=X[[j]], S=S, g=g, M=M, method=method) 
  #       }
  #     }
  #     
  #     H.hat <- array(NA, dim = c(n,p,p))
  #     G.hat <- array(NA, dim = c(n,p,p))
  #     tau.hat <- matrix(NA, nrow=n, ncol=p)
  #     
  #     for(l in 1:n){ # \lambda_1 <= ... <= \lambda_n
  #       tmp <- eigensort(P_array[l,,])
  #       G.hat[l,,] <- tmp$evectors[,p:1]
  #       H.hat[l,,] <- t(Conj(G.hat[l,,]))
  #       tau.hat[l,] <- rev(tmp$evalues)
  #     }
  #     
  #     if(is.null(q)){
  #       q <- p
  #     }
  #     
  #     Y <- list()
  #     # Y[[i]] : dimension-reduced graph signal for R realizations (n x R matrix) / length q
  #     for(i in 1:q){
  #       Y[[i]] <- matrix(0, nrow=n, ncol=R)
  #       for(j in 1:p){
  #         Y[[i]] <- Y[[i]] + evectors %*% diag(H.hat[,i,j]) %*% t(Conj(evectors)) %*% X[[j]]
  #       }
  #     }
  #     mu.x.hat <- do.call(cbind, lapply(X, function(x) rep(mean(x), n))) # n x p matrix
  #     # mu.x.hat <- do.call(cbind, X) # n x p matrix
  #     # window로 하는방법 생각
  #     
  #     mu.y.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
  #     
  #     for(j in 1:p){
  #       for(k in 1:p){
  #         mu.y.hat[,j] <- mu.y.hat[,j] + evectors %*% diag(H.hat[,j,k]) %*% t(Conj(evectors)) %*% mu.x.hat[,k]
  #       }
  #     }
  #     
  #     mu.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
  #     X.hat <- list()
  #     # X.hat[[i]] : reconstructed graph signal for R realizations (n x R matrix) / length p
  #     for(i in 1:p){
  #       mu.hat[,i] <- mu.x.hat[,i] 
  #       X.hat[[i]] <- matrix(0, nrow=n, ncol=R)
  #       for(j in 1:q){
  #         mu.hat[,i] <- mu.hat[,i] - evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% mu.y.hat[,j]
  #         X.hat[[i]] <- X.hat[[i]] + evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% Y[[j]]
  #       }
  #       X.hat[[i]] <- X.hat[[i]] + mu.hat[,i]
  #     }  
  #     
  #     res <- list()
  #     res$X <- lapply(as.vector(X)) ; res$Y <- lapply(as.vector(Y)) ; res$X.hat <- lapply(as.vector(X.hat))
  #     res$mu.hat <- mu.hat ; res$tau.hat <- tau.hat ; res$H.hat <- H.hat ; res$G.hat <- G.hat
  #     return(res)
  #   } 
  else if(method=="random"){
    # X[[i]] should be one realization (n x 1 vector)
    # p-dimensional graph signal
    if(is.vector(X[[1]])){
      X <- lapply(X, as.matrix)
    }
    p <- length(X)
    if(nrow(X[[1]])!=nrow(S)){
      stop("The length of list X's element vector OR the number of rows of list X's element matrices should be equal to the number of rows of GSO!")
    }
    n <- nrow(S)
    R <- ncol(X[[1]])
    
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    # evalues.rev <- rev(evalues)
    # evectors.rev <- val$evectors[,n:1]
    
    P_array <- array(NA, dim = c(n,p,p))
    for(i in 1:p){
      for(j in 1:p){
        # WB <- windowbank.random(N=n, M=M, V=evectors, sigma=sigma) # M x n
        P_array[,i,j] <- cpsd.graph(x=X[[i]], y=X[[j]], S=S, g=g, M=M, sigma=sigma, method=method, seed=seed)
      }
    }
    
    H.hat <- array(NA, dim = c(n,p,p))
    G.hat <- array(NA, dim = c(n,p,p))
    tau.hat <- matrix(NA, nrow=n, ncol=p)
    
    for(l in 1:n){ # \lambda_1 <= ... <= \lambda_n
      tmp <- eigensort(P_array[l,,])
      G.hat[l,,] <- tmp$evectors[,p:1]
      H.hat[l,,] <- t(Conj(G.hat[l,,]))
      tau.hat[l,] <- rev(tmp$evalues)
    }
    
    if(is.null(q)){
      q <- p
    }
    
    Y <- list()
    # Y[[i]] : dimension-reduced graph signal for R realizations (n x R matrix) / length q
    for(i in 1:q){
      Y[[i]] <- matrix(0, nrow=n, ncol=R)
      for(j in 1:p){
        Y[[i]] <- Y[[i]] + evectors %*% diag(H.hat[,i,j]) %*% t(Conj(evectors)) %*% X[[j]]
      }
    }
    
    WB <- windowbank.random(N=n, M=M, V=evectors, sigma=sigma, seed=seed)
    mu.x.hat <- do.call(cbind, lapply(X, function(x) rowMeans(t(WB)*as.vector(x)))) # n x p matrix
    
    mu.y.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
    
    for(j in 1:p){
      for(k in 1:p){
        mu.y.hat[,j] <- mu.y.hat[,j] + evectors %*% diag(H.hat[,j,k]) %*% t(Conj(evectors)) %*% mu.x.hat[,k]
      }
    }
    
    mu.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
    X.hat <- list()
    # X.hat[[i]] : reconstructed graph signal for R realizations (n x R matrix) / length p
    for(i in 1:p){
      mu.hat[,i] <- mu.x.hat[,i] 
      X.hat[[i]] <- matrix(0, nrow=n, ncol=R)
      for(j in 1:q){
        mu.hat[,i] <- mu.hat[,i] - evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% mu.y.hat[,j]
        X.hat[[i]] <- X.hat[[i]] + evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% Y[[j]]
      }
      X.hat[[i]] <- X.hat[[i]] + mu.hat[,i]
    }  
    
    res <- list()
    res$X <- lapply(X, as.vector) ; res$Y <- lapply(Y, as.vector) ; res$X.hat <- lapply(X.hat, as.vector)
    res$mu.hat <- mu.hat ; res$tau.hat <- tau.hat ; res$H.hat <- H.hat ; res$G.hat <- G.hat ; res$Px <- P_array
    return(res)
  }
}

plot_graph_custom4 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="gray", high="gray", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}

plot_graph_custom5 <- function (z, e.size=1, v.size=3, vertex_color=NULL, min=NULL, max=NULL, value, label.title.size=15, label.text.size=10, ratio=1, signal=TRUE, mg=c(10,10,10,10), title="", main.title.size=30) 
{
  if (is(z$sA, "sparseMatrix")) {
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  # w <- paste(rownames(z$xy),"_",c(1:nrow(z$xy)), sep="")
  w <- rownames(z$xy)
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y, w=w)
  df2 <- data.frame(x = x[ind_i], y = y[ind_i], y1 = y1, y2 = y2)
  df2$w <- z$sA[,3]
  if(is.null(min)){
    min <- min(vertex_color)
  }
  if(is.null(max)){
    max <- max(vertex_color)
  }
  
  if(signal==FALSE){
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2) +
      scale_color_gradient(low="gray", high="gray", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = "black", limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(fill=guide_colourbar(order=1), color=guide_colourbar(order=2)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
  
  else{
    p1 <- ggplot(df1, aes(x, y)) + geom_segment(aes(x = x, y = y, 
                                                    xend = y1, yend = y2, color=w), linewidth=e.size, 
                                                data = df2, show.legend=FALSE) +
      scale_color_gradient(low="grey", high="black", name="Edge Weight")+
      geom_point(aes(fill=vertex_color), size = v.size, shape=21) + 
      scale_fill_gradientn(colors = colorRampPalette(c("blue", "skyblue", "green", "yellow", "orange", "red"))(500), limits=c(min, max), na.value = "gray", name = value) +
      # geom_text(aes(label=w), data=df1, hjust=0.6, vjust = -1.1, size = 4) +
      theme_void() +
      guides(color=guide_colourbar(order=1)) + ggtitle(title) +
      theme(legend.margin = margin(mg), plot.margin = margin(mg), 
            plot.title=element_text(size=main.title.size, face="bold", hjust = 0.5), legend.title=element_text(size=label.title.size),
            legend.text = element_text(size=label.text.size), aspect.ratio=ratio)
    print(p1)
  }
}

GFPCA <- function(X, S, q=NULL, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  if(is.null(method)){
    # X[[i]] should be R realizations (n x R matrix)
    # p-dimensional graph signal
    p <- length(X)
    if(nrow(X[[1]])!=nrow(S)){
      stop("The number of rows of list X's element matrices should be equal to the number of rows of GSO!")
    }
    n <- nrow(S)
    R <- ncol(X[[1]])
    
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    # evalues.rev <- rev(evalues)
    # evectors.rev <- val$evectors[,n:1]
    
    P_array <- array(NA, dim = c(n,p,p))
    for(i in 1:p){
      for(j in 1:p){
        P_array[,i,j] <- cpsd.graph(x=X[[i]], y=X[[j]], S=S) 
      }
    }
    
    H.hat <- array(NA, dim = c(n,p,p))
    G.hat <- array(NA, dim = c(n,p,p))
    tau.hat <- matrix(NA, nrow=n, ncol=p)
    
    for(l in 1:n){ # \lambda_1 <= ... <= \lambda_n
      tmp <- eigensort(P_array[l,,])
      G.hat[l,,] <- tmp$evectors[,p:1]
      H.hat[l,,] <- t(Conj(G.hat[l,,]))
      tau.hat[l,] <- rev(tmp$evalues)
    }
    
    if(is.null(q)){
      q <- p
    }
    
    Y <- list()
    # Y[[i]] : dimension-reduced graph signal for R realizations (n x R matrix) / length q
    for(i in 1:q){
      Y[[i]] <- matrix(0, nrow=n, ncol=R)
      for(j in 1:p){
        Y[[i]] <- Y[[i]] + evectors %*% diag(H.hat[,i,j]) %*% t(Conj(evectors)) %*% X[[j]]
      }
    }
    mu.x.hat <- do.call(cbind, lapply(X, FUN=rowMeans)) # n x p matrix
    
    mu.y.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
    
    for(j in 1:p){
      for(k in 1:p){
        mu.y.hat[,j] <- mu.y.hat[,j] + evectors %*% diag(H.hat[,j,k]) %*% t(Conj(evectors)) %*% mu.x.hat[,k]
      }
    }
    
    mu.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
    X.hat <- list()
    # X.hat[[i]] : reconstructed graph signal for R realizations (n x R matrix) / length p
    for(i in 1:p){
      mu.hat[,i] <- mu.x.hat[,i] 
      X.hat[[i]] <- matrix(0, nrow=n, ncol=R)
      for(j in 1:q){
        mu.hat[,i] <- mu.hat[,i] - evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% mu.y.hat[,j]
        X.hat[[i]] <- X.hat[[i]] + evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% Y[[j]]
      }
      X.hat[[i]] <- X.hat[[i]] + mu.hat[,i]
    }  
    
    res <- list()
    res$X <- X ; res$Y <- Y ; res$X.hat <- X.hat
    res$mu.hat <- mu.hat ; res$tau.hat <- tau.hat ; res$H.hat <- H.hat ; res$G.hat <- G.hat ; res$Px <- P_array
    return(res)
  } 
  # else if(method=="window"){
  #     # X[[i]] should be one realization (n x 1 vector)
  #     # p-dimensional graph signal
  #     if(is.vector(X[[1]])){
  #       X <- lapply(X, as.matrix)
  #     }
  #     p <- length(X)
  #     if(nrow(X[[1]])!=nrow(S)){
  #       stop("The length of list X's element vector OR the number of rows of list X's element matrices should be equal to the number of rows of GSO!")
  #     }
  #     n <- nrow(S)
  #     R <- ncol(X[[1]])
  #     
  #     val <- eigensort(S)
  #     evalues <- val$evalues
  #     evectors <- val$evectors
  #     # evalues.rev <- rev(evalues)
  #     # evectors.rev <- val$evectors[,n:1]
  #     
  #     P_array <- array(NA, dim = c(n,p,p))
  #     for(i in 1:p){
  #       for(j in 1:p){
  #         P_array[,i,j] <- cpsd.graph(x=X[[i]], y=X[[j]], S=S, g=g, M=M, method=method) 
  #       }
  #     }
  #     
  #     H.hat <- array(NA, dim = c(n,p,p))
  #     G.hat <- array(NA, dim = c(n,p,p))
  #     tau.hat <- matrix(NA, nrow=n, ncol=p)
  #     
  #     for(l in 1:n){ # \lambda_1 <= ... <= \lambda_n
  #       tmp <- eigensort(P_array[l,,])
  #       G.hat[l,,] <- tmp$evectors[,p:1]
  #       H.hat[l,,] <- t(Conj(G.hat[l,,]))
  #       tau.hat[l,] <- rev(tmp$evalues)
  #     }
  #     
  #     if(is.null(q)){
  #       q <- p
  #     }
  #     
  #     Y <- list()
  #     # Y[[i]] : dimension-reduced graph signal for R realizations (n x R matrix) / length q
  #     for(i in 1:q){
  #       Y[[i]] <- matrix(0, nrow=n, ncol=R)
  #       for(j in 1:p){
  #         Y[[i]] <- Y[[i]] + evectors %*% diag(H.hat[,i,j]) %*% t(Conj(evectors)) %*% X[[j]]
  #       }
  #     }
  #     mu.x.hat <- do.call(cbind, lapply(X, function(x) rep(mean(x), n))) # n x p matrix
  #     # mu.x.hat <- do.call(cbind, X) # n x p matrix
  #     # window로 하는방법 생각
  #     
  #     mu.y.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
  #     
  #     for(j in 1:p){
  #       for(k in 1:p){
  #         mu.y.hat[,j] <- mu.y.hat[,j] + evectors %*% diag(H.hat[,j,k]) %*% t(Conj(evectors)) %*% mu.x.hat[,k]
  #       }
  #     }
  #     
  #     mu.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
  #     X.hat <- list()
  #     # X.hat[[i]] : reconstructed graph signal for R realizations (n x R matrix) / length p
  #     for(i in 1:p){
  #       mu.hat[,i] <- mu.x.hat[,i] 
  #       X.hat[[i]] <- matrix(0, nrow=n, ncol=R)
  #       for(j in 1:q){
  #         mu.hat[,i] <- mu.hat[,i] - evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% mu.y.hat[,j]
  #         X.hat[[i]] <- X.hat[[i]] + evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% Y[[j]]
  #       }
  #       X.hat[[i]] <- X.hat[[i]] + mu.hat[,i]
  #     }  
  #     
  #     res <- list()
  #     res$X <- lapply(as.vector(X)) ; res$Y <- lapply(as.vector(Y)) ; res$X.hat <- lapply(as.vector(X.hat))
  #     res$mu.hat <- mu.hat ; res$tau.hat <- tau.hat ; res$H.hat <- H.hat ; res$G.hat <- G.hat
  #     return(res)
  #   } 
  else if(method=="random"){
    # X[[i]] should be one realization (n x 1 vector)
    # p-dimensional graph signal
    if(is.vector(X[[1]])){
      X <- lapply(X, as.matrix)
    }
    p <- length(X)
    if(nrow(X[[1]])!=nrow(S)){
      stop("The length of list X's element vector OR the number of rows of list X's element matrices should be equal to the number of rows of GSO!")
    }
    n <- nrow(S)
    R <- ncol(X[[1]])
    
    val <- eigensort(S)
    evalues <- val$evalues
    evectors <- val$evectors
    # evalues.rev <- rev(evalues)
    # evectors.rev <- val$evectors[,n:1]
    
    P_array <- array(NA, dim = c(n,p,p))
    for(i in 1:p){
      for(j in 1:p){
        # WB <- windowbank.random(N=n, M=M, V=evectors, sigma=sigma) # M x n
        P_array[,i,j] <- cpsd.graph(x=X[[i]], y=X[[j]], S=S, g=g, M=M, sigma=sigma, method=method, seed=seed)
      }
    }
    
    H.hat <- array(NA, dim = c(n,p,p))
    G.hat <- array(NA, dim = c(n,p,p))
    tau.hat <- matrix(NA, nrow=n, ncol=p)
    
    for(l in 1:n){ # \lambda_1 <= ... <= \lambda_n
      tmp <- eigensort(P_array[l,,])
      G.hat[l,,] <- tmp$evectors[,p:1]
      H.hat[l,,] <- t(Conj(G.hat[l,,]))
      tau.hat[l,] <- rev(tmp$evalues)
    }
    
    if(is.null(q)){
      q <- p
    }
    
    Y <- list()
    # Y[[i]] : dimension-reduced graph signal for R realizations (n x R matrix) / length q
    for(i in 1:q){
      Y[[i]] <- matrix(0, nrow=n, ncol=R)
      for(j in 1:p){
        Y[[i]] <- Y[[i]] + evectors %*% diag(H.hat[,i,j]) %*% t(Conj(evectors)) %*% X[[j]]
      }
    }
    
    WB <- windowbank.random(N=n, M=M, V=evectors, sigma=sigma, seed=seed)
    mu.x.hat <- do.call(cbind, lapply(X, function(x) rowMeans(t(WB)*as.vector(x)))) # n x p matrix
    
    mu.y.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
    
    for(j in 1:p){
      for(k in 1:p){
        mu.y.hat[,j] <- mu.y.hat[,j] + evectors %*% diag(H.hat[,j,k]) %*% t(Conj(evectors)) %*% mu.x.hat[,k]
      }
    }
    
    mu.hat <- matrix(0, nrow=n, ncol=p) # n x p matrix
    X.hat <- list()
    # X.hat[[i]] : reconstructed graph signal for R realizations (n x R matrix) / length p
    for(i in 1:p){
      mu.hat[,i] <- mu.x.hat[,i] 
      X.hat[[i]] <- matrix(0, nrow=n, ncol=R)
      for(j in 1:q){
        mu.hat[,i] <- mu.hat[,i] - evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% mu.y.hat[,j]
        X.hat[[i]] <- X.hat[[i]] + evectors %*% diag(G.hat[,i,j]) %*% t(Conj(evectors)) %*% Y[[j]]
      }
      X.hat[[i]] <- X.hat[[i]] + mu.hat[,i]
    }  
    
    res <- list()
    res$X <- lapply(X, as.vector) ; res$Y <- lapply(Y, as.vector) ; res$X.hat <- lapply(X.hat, as.vector)
    res$mu.hat <- mu.hat ; res$tau.hat <- tau.hat ; res$H.hat <- H.hat ; res$G.hat <- G.hat ; res$Px <- P_array
    return(res)
  }
}

# determine #(factors)
theta_func <- function(X, B, Ghat){
  # Ghat : R x k
  k <- ncol(Ghat)
  R <- nrow(Ghat)
  p <- nrow(B)
  return(norm(t(X) - B%*%t(Ghat), type="F")^2 / (p*R))
}

g_penality <- function(p,R,type=1){
  Cpr <- min(sqrt(p),sqrt(R))
  g1 <- (p+R)/(p*R) * log(p*R / (p+R))
  g2 <- (p+R)/(p*R) * log(Cpr^2)
  g3 <- log(Cpr^2) / (Cpr^2)
  if(type==2){
    return(g2)
  } else if(type==3){
    return(g3)
  }
  return(g1)
}

IC_criterion <- function(X, B, Ghat, g_penality, type=1){
  k <- ncol(Ghat)
  R <- nrow(Ghat) # Ghat: R x k
  p <- nrow(B) # B: p x k
  
  return(log(theta_func(X, B, Ghat)) + k*g_penality(p,R,type))
}


gfactor <- function(X, S, q=NULL, g=NULL, M=NULL, sigma=NULL, method=NULL, seed=1){
  # X: p x n x R
}