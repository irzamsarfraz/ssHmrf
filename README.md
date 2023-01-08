# ssHmrf

ssHMRF is an R package built on top of smFishHmrf to allow easy identification of spatial domains from spatial transcriptomics data using hidden markov random fields.

## Package Usage
Below we describe an example to use the package using a sample dataset provided with the package:

### Import Data
```
expression_filepath <- "fcortex.expression.txt"
coordinates_filepath <- "fcortex.coordinates.txt"
genes_filepath <- "genes"

expr <- readExpression(expression_filepath)
coords <- readCoordinates(coordinates_filepath)
fields <- coords$field
coords <- coords$coords
genes <- readGenes(genes_filepath)
```

### Create Spatial Network
```
spatial_network <- createSpatialNetwork(coords = coords, k = 4)
```

### Initialize HMRF Parameters
```
hmrfParams <- initializeHMRF(expr = expr, 
                             spatial_network = spatial_network, 
                             spatial_genes = genes, 
                             k = 3)
```                           

### Run HMRF
```
hmrfResults <- runHMRF(expr = expr, 
                       nei = hmrfParams$nei, 
                       numnei = hmrfParams$numnei, 
                       blocks = hmrfParams$blocks, 
                       damp = hmrfParams$damp, 
                       mu = hmrfParams$mu, 
                       sigma = hmrfParams$sigma, 
                       edgelist = hmrfParams$edgelist, 
                       k = 3, 
                       genes = genes, 
                       beta_initial_value = 0, 
                       beta_increment_value = 10, 
                       beta_iterations = 5)
```

### Visualize Results
```
hmrfPlots <- plotHMRF(spatial_network = spatial_network, HMRFoutput = hmrfResults)
hmrfPlots
```
