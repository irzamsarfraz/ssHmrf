#' @title createSpatialNetwork
#' @description Creates a KNN based spatial network from
#' input spatial coordinates.
#' @param coords Input a named list where X element contains
#' x coordinates and Y element contains y coordinates.
#' @param k Specify number of initial clusters to generate.
#' Default \code{4}.
#'
#' @return A data.table containg the spatial network.
#' @export
createSpatialNetwork <- function(coords, k = 4){
  spatial_locations <- data.table(sdimx = coords$X, sdimy = coords$Y, cell_ID = seq(1, nrow(coords)))
  temp_spatial_locations = spatial_locations[, grepl('sdim', colnames(spatial_locations)), with = FALSE]
  temp_spatial_locations = as.matrix(temp_spatial_locations)
  first_dimension = colnames(temp_spatial_locations)[[1]]
  second_dimension = colnames(temp_spatial_locations)[[2]]

  from = to = NULL
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))


  spatial_locations_matrix = as.matrix(spatial_locations[, c("sdimx", "sdimy"), with = F])
  knn_spatial <- dbscan::kNN(x = spatial_locations_matrix,
                             k = k)

  knn_sptial.norm = data.frame(from = rep(1:nrow(knn_spatial$id), k),
                               to = as.vector(knn_spatial$id),
                               weight = 1/(1 + as.vector(knn_spatial$dist)),
                               distance = as.vector(knn_spatial$dist))
  nw_sptial.norm = igraph::graph_from_data_frame(knn_sptial.norm, directed = FALSE)
  network_DT = data.table::as.data.table(knn_sptial.norm)

  xbegin_name = "sdimx_begin"
  ybegin_name = "sdimy_begin"
  xend_name = "sdimx_end"
  yend_name = "sdimy_end"

  spatial_network_DT = data.table::data.table(from = cell_ID_vec[network_DT$from],
                                              to = cell_ID_vec[network_DT$to],
                                              xbegin_name = spatial_locations[network_DT$from, "sdimx"],
                                              ybegin_name = spatial_locations[network_DT$from, "sdimy"],
                                              xend_name = spatial_locations[network_DT$to, "sdimx"],
                                              yend_name = spatial_locations[network_DT$to, "sdimy"],
                                              distance = network_DT$distance,
                                              weight = network_DT$weight)

  data.table::setnames(spatial_network_DT,
                       old = c('xbegin_name', 'ybegin_name', 'xend_name', 'yend_name'),
                       new = c(xbegin_name, ybegin_name, xend_name, yend_name))
  data.table::setorder(spatial_network_DT, from, to)

  print(head(spatial_network_DT))

  spatial_network <- spatial_network_DT

  return(spatial_network)
}

#' @title initializeHMRF
#' @description Initializes the parameters of the HMRF using
#' the expression matrix and the spatial network previously
#' computed.
#' @param y Input expression matrix where rows are cells
#' and columns are features.
#' @param spatial_network Input a spatial network previously
#' generated from \code{\link{createSpatialNetwork}} function.
#' @param spatial_genes A character vector of names of genes
#' to use with HMRF. It is recommended to use only the top
#' most spatially variable genes at this step.
#' @param k Specify a k value which defines the number of
#' spatial domains to identify.
#' @param hmrf_seed Set a seed value
#' @param nstart Set a nstart value
#' @param factor_step Set a factor_step value
#' @param tolerance Set a tolerance value
#'
#' @return A list of initialized parameters for HMRF.
#' @export
initializeHMRF <- function(y, spatial_network, spatial_genes, k = 3, hmrf_seed = 100, nstart = 1000, factor_step = 1.05, tolerance = 1e-5){
  numcell <- dim(y)[1]
  m<-dim(y)[2]
  ncol.nei = max(table(c(spatial_network$to,spatial_network$from)))
  nei = matrix(-1,ncol = ncol.nei,nrow = numcell)
  rownames(nei) = rownames(y)
  for(i in 1:numcell)
  {
    nei.i = c(spatial_network$from[spatial_network$to==rownames(nei)[i]],
              spatial_network$to[spatial_network$from==rownames(nei)[i]])
    if(length(nei.i)>0)nei[i,1:length(nei.i)] = sort(match(nei.i,rownames(y)))
  }
  numnei<-as.integer(rowSums(nei!=(-1)))
  nn<-nei
  numedge <- 0
  for(i in 1:numcell){
    numedge<-numedge + length(nn[i,nn[i,]!=-1])
  }
  edgelist <- matrix(0, nrow=numedge, ncol=2)
  edge_ind <- 1
  for(i in 1:numcell){
    neighbors <- nn[i, nn[i,]!=-1]
    for(j in 1:length(neighbors)){
      edgelist[edge_ind,] <- c(i, neighbors[j])
      edge_ind<-edge_ind+1
    }
  }

  pp<-tidygraph::tbl_graph(edges=as.data.frame(edgelist), directed=F)
  yy<-pp%>%mutate(color=as.factor(color_dsatur()))
  colors<-as.list(yy)$nodes$color
  cl_color <- sort(unique(colors))
  blocks<-lapply(cl_color, function(cl){which(colors==cl)})

  kk = smfishHmrf.generate.centroid(y=y,par_k = k,par_seed=hmrf_seed,nstart=nstart)
  mu<-t(kk$centers) #should be dimension (m,k)
  lclust<-lapply(1:k, function(x) which(kk$cluster == x))
  damp<-array(0, c(k));
  sigma<-array(0, c(m,m,k))
  for(i in 1:k){
    sigma[, , i] <- cov(as.matrix(y[lclust[[i]], ]))
    di<-findDampFactor(sigma[,,i], factor=factor_step, d_cutoff=tolerance, startValue=0.0001)
    damp[i]<-ifelse(is.null(di), 0, di)
  }

  spatial_genes_selected <- spatial_genes #genes # should run binspect here

  inputList <- list(y=y, nei=nei, numnei=numnei, blocks=blocks,
                    damp=damp, mu=mu, sigma=sigma, k=k, genes=spatial_genes_selected, edgelist=edgelist)

  return(inputList)
}

#' runHMRF
#'
#' @param y Input expression matrix where rows are cells
#' and columns are features.
#' @param nei Specify the neighbours against each cell.
#' @param numnei Specify the number of neighbours
#' against each cell.
#' @param blocks Specify the groups of clusters computed
#' by a graph coloring algorithm where each no two
#' adjcacent cells have same colors.
#' @param damp Set damp param
#' @param mu Set mu param
#' @param sigma Set sigma param
#' @param k Set k value that defines the number of
#' spatial domains to identify.
#' @param genes Specify list of genes
#' @param edgelist Specify all edges between cells.
#' @param beta_initial_value Specify initial beta value.
#' @param beta_increment_value Specify increment value.
#' @param beta_iterations Specify total iterations.
#'
#' @return A list containing results of HMRF.
#' @export
runHMRF <- function(y, nei, numnei, blocks,
                    damp, mu, sigma, k, genes, edgelist, beta_initial_value = 0, beta_increment_value = 10, beta_iterations = 5){
  y <- as.matrix(y)
  #betas <- c(0,10,5)
  #beta_init = betas[1]
  beta_init <- beta_initial_value
  #beta_increment = betas[2]
  beta_increment = beta_increment_value
  #beta_num_iter = betas[3]
  beta_num_iter = beta_iterations


  beta_seq = (1:beta_num_iter-1)*beta_increment+beta_init
  # beta_seq = sequence(beta_num_iter,beta_init,beta_increment)
  beta_seq = sort(unique(c(0,beta_seq)))
  res <- c()
  for(beta_current in beta_seq){
    print(sprintf("Doing beta=%.3f", beta_current))
    tc.hmrfem<-smfishHmrf.hmrfem.multi(y=y, neighbors=nei, beta=beta_current, numnei=numnei,
                                       blocks=blocks, mu=mu, sigma=sigma, verbose=T, err=1e-7, maxit=50, dampFactor=damp)
    #smfishHmrf.hmrfem.multi.save(name, outdir, beta_current, tc.hmrfem, k)
    #do_one(name, outdir, k, y, nei, beta_current, numnei, blocks, mu, sigma, damp)

    ### stop the loop if there is a samll maximum probablity (<1/k+0.05v) of any gene
    if(sum(apply(tc.hmrfem$prob,1,max)<(1/k+0.05))>0)
    {cat(paste0('\n HMRF is stopping at large beta >= ',beta_current,', numerical error occurs, results of smaller betas were stored\n'));
      break()}

    t_key <- sprintf("k=%d b=%.2f", k, beta_current)
    tc.hmrfem$sigma=NULL
    tc.hmrfem$mu=NULL
    rownames(tc.hmrfem$prob) = rownames(y)
    rownames(tc.hmrfem$unnormprob) = rownames(y)
    names(tc.hmrfem$class) = rownames(y)
    res[[t_key]] <- tc.hmrfem
    # beta_current <- beta_current + beta_increment
  }

  result.hmrf = res
  HMRFoutput <- NULL
  HMRFoutput$res <- result.hmrf
  HMRFoutput$beta_seq <- beta_seq
  HMRFoutput$k <- k
  return(HMRFoutput)
}

#' @title plotHMRF
#' @description Plots the spatial domains
#' identified previously on a 2D spatial plot.
#' @param spatial_network Input the spatial network.
#' @param HMRFoutput Input the HMRF results.
#'
#' @return A list of plots.
#' @export
plotHMRF <- function(spatial_network, HMRFoutput){
  plots <- list()
  t_key = sprintf("k=%d b=%.2f", HMRFoutput$k, HMRFoutput$beta_seq)
  for(kk in intersect(t_key,names(HMRFoutput$res)))
  {
    dt_kk = HMRFoutput$res[[kk]]$class
    spatial_network$color <- dt_kk[spatial_network$from]
    p <- ggplot(spatial_network, aes(sdimx_begin, sdimy_begin)) + geom_point(aes(colour = color))
    plots[[kk]] <- p
  }
  return(plots)
}
