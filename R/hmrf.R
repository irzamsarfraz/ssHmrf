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

  # Tabulate spatial locations with cell_ID
  spatial_locations <- data.table::data.table(sdimx = coords$sdimx,
                                              sdimy = coords$sdimy,
                                              cell_ID = seq(1, nrow(coords)))

  # Run kNN from dbscan (needs matrix)
  spatial_locations_matrix = as.matrix(spatial_locations[, c("sdimx", "sdimy"), with = F])
  knn_spatial <- dbscan::kNN(x = spatial_locations_matrix, k = k)

  # Output from dbscan is a list of values so converting into a data.frame
  network_DT <- data.frame(from = rep(1:nrow(knn_spatial$id), k),
                               to = as.vector(knn_spatial$id),
                               weight = 1/(1 + as.vector(knn_spatial$dist)),
                               distance = as.vector(knn_spatial$dist))

  # Adding x & y coordinates for cells (begin & end)
  network_DT = data.table::data.table("from" = spatial_locations$cell_ID[network_DT$from],
                                      "to" = spatial_locations$cell_ID[network_DT$to],
                                      "sdimx_begin" = spatial_locations[network_DT$from, "sdimx"],
                                      "sdimy_begin" = spatial_locations[network_DT$from, "sdimy"],
                                      "sdimx_end" = spatial_locations[network_DT$to, "sdimx"],
                                      "sdimy_end" = spatial_locations[network_DT$to, "sdimy"],
                                      "distance" = network_DT$distance,
                                      "weight" = network_DT$weight)

  return(network_DT)
}

#' @title initializeHMRF
#' @description Initializes the parameters of the HMRF using
#' the expression matrix and the spatial network previously
#' computed.
#' @param expr Input expression matrix where rows are cells
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
initializeHMRF <- function(expr,
                           spatial_network,
                           spatial_genes,
                           k = 3,
                           hmrf_seed = 100,
                           nstart = 1000,
                           factor_step = 1.05,
                           tolerance = 1e-5){

  # Get number of cells and features
  num_cells <- dim(expr)[1]
  num_features <-dim(expr)[2]

  # Find max number of neighbours for any cell
  max_neighbours = max(table(c(spatial_network$to, spatial_network$from)))

  # Create a matrix with columns equal to max number of neighbours with default values -1
  #   so neighbouring cell indices could be inserted later for connected neighbours of each cell
  neighbour_matrix = matrix(-1, ncol = max_neighbours, nrow = num_cells)
  rownames(neighbour_matrix) = rownames(expr)

  # Iterate over each cell and find which cells are neighbours of which cells
  for(i in 1:num_cells)
  {
    neighbour_matrix.i = c(spatial_network$from[spatial_network$to==rownames(neighbour_matrix)[i]],
              spatial_network$to[spatial_network$from==rownames(neighbour_matrix)[i]])
    if(length(neighbour_matrix.i)>0)neighbour_matrix[i,1:length(neighbour_matrix.i)] = sort(match(neighbour_matrix.i,rownames(expr)))
  }

  # For each cell count number of neighbours
  numnei<-as.integer(rowSums(neighbour_matrix!=(-1)))
  nn<-neighbour_matrix

  # Identify and store all edges between all cells
  numedge <- 0
  for(i in 1:num_cells){
    numedge<-numedge + length(nn[i,nn[i,]!=-1])
  }
  edgelist <- matrix(0, nrow=numedge, ncol=2)
  edge_ind <- 1
  for(i in 1:num_cells){
    neighbors <- nn[i, nn[i,]!=-1]
    for(j in 1:length(neighbors)){
      edgelist[edge_ind,] <- c(i, neighbors[j])
      edge_ind<-edge_ind+1
    }
  }

  # Convert edgeList into a graph structure for easier manipulation
  pp<-tidygraph::tbl_graph(edges=as.data.frame(edgelist), directed=F)

  # Use a graph coloring algorithm to assign colors to cells so no two adjacent cells have same color
  yy<-pp%>%mutate(color=as.factor(color_dsatur()))
  colors<-as.list(yy)$nodes$color
  cl_color <- sort(unique(colors))

  # Form a block for each color and includes all cells for that color
  blocks<-lapply(cl_color, function(cl){which(colors==cl)})

  # Generate initial centroids from expression data
  kk = smfishHmrf::smfishHmrf.generate.centroid(y=expr, par_k = k, par_seed=hmrf_seed, nstart=nstart)

  # Cluster centers for each feature in each cluster
  mu<-t(kk$centers) #should be dimension (num_features, k)
  lclust<-lapply(1:k, function(x) which(kk$cluster == x))


  # Damp values for each k cluster using findDampFactor() function
  #  which uses sigma where sigma is cov for each cluster using cells in that cluster
  damp<-array(0, c(k));
  sigma<-array(0, c(num_features,num_features,k))
  for(i in 1:k){
    sigma[, , i] <- cov(as.matrix(expr[lclust[[i]], ]))
    di<-findDampFactor(sigma[,,i], factor=factor_step, d_cutoff=tolerance, startValue=0.0001)
    damp[i]<-ifelse(is.null(di), 0, di)
  }

  # Spatially variable genes selected for downstream
  spatial_genes_selected <- spatial_genes #genes # should run binspect here .. TODO

  inputList <- list(expr=expr, neighbour_matrix=neighbour_matrix, numnei=numnei, blocks=blocks,
                    damp=damp, mu=mu, sigma=sigma, k=k, genes=spatial_genes_selected, edgelist=edgelist)

  return(inputList)
}

#' runHMRF
#'
#' @param expr Input expression matrix where rows are cells
#' and columns are features.
#' @param neighbour_matrix Specify the neighbours against each cell.
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
runHMRF <- function(expr,
                    neighbour_matrix,
                    numnei,
                    blocks,
                    damp,
                    mu,
                    sigma,
                    k,
                    genes,
                    edgelist,
                    beta_initial_value = 0,
                    beta_increment_value = 10,
                    beta_iterations = 5){

  expr <- as.matrix(expr)

  # Run HMRF for each beta value
  beta_seq = (1:beta_iterations-1)*beta_increment_value+beta_initial_value
  beta_seq = sort(unique(c(0,beta_seq)))
  res <- c()
  for(beta_current in beta_seq){
    print(sprintf("Doing beta=%.3f", beta_current))
    tc.hmrfem<- smfishHmrf::smfishHmrf.hmrfem.multi(y=expr, neighbors=neighbour_matrix, beta=beta_current, numnei=numnei,
                                       blocks=blocks, mu=mu, sigma=sigma, verbose=T, err=1e-7, maxit=50, dampFactor=damp)

    ### stop the loop if there is a small maximum probablity (<1/k+0.05v) of any gene
    if(sum(apply(tc.hmrfem$prob,1,max)<(1/k+0.05))>0)
    {cat(paste0('\n HMRF is stopping at large beta >= ',beta_current,', numerical error occurs, results of smaller betas were stored\n'));
      break()}

    t_key <- sprintf("k=%d b=%.2f", k, beta_current)
    tc.hmrfem$sigma=NULL
    tc.hmrfem$mu=NULL
    rownames(tc.hmrfem$prob) = rownames(expr)
    rownames(tc.hmrfem$unnormprob) = rownames(expr)
    names(tc.hmrfem$class) = rownames(expr)
    res[[t_key]] <- tc.hmrfem
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
