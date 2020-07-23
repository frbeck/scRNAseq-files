##### libraries #####

library(dplyr)
library(Seurat)
library(patchwork)
library(leidenbase)
library(monocle3)
library(scater)
library(loomR)
library(patchwork)
library(slingshot)
library(topGO)
library(PANTHER.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichR)
library(UniProt.ws)
library(DBI)
library(AnnotationDbi)
library(reshape2)
library(class)
library(caret)
library(rpart)
library(rpart.plot)
library(C50)
library(hdf5r)
library("Matrix")

##### PCA_run  #####

PCA_run = function(tissue){
  tissue = NormalizeData(tissue)
  tissue = FindVariableFeatures(tissue)
  tissue = ScaleData(tissue, features = all.genes)
  tissue = RunPCA(tissue, features = VariableFeatures(object = tissue))
  return(tissue)
}

##### umap_tsne #####

umap_tsne = function(tissue, dims, res = NULL){
  if(is.null(res)){
    res = 1
  }
  tissue <- FindNeighbors(tissue, dims = 1:dims)
  tissue <- FindClusters(tissue, resolution = res)
  tissue <- RunUMAP(tissue, dims = 1:dims)
  tissue = RunTSNE(tissue, dims = 1:dims)
  return(tissue)
}

##### process_data #####

process_data = function(tissue){
  tissue = umap_tsne(PCA_run(tissue), dims = determine_dims(tissue))
  return(tissue)
}

##### correlation_analysis_1_gene #####

cor_analysis_1_gene = function(tissue, feature, nhits) {
  tissue.mtx = tissue[["RNA"]]@data
  tissue.mtx.md = as.matrix(tissue.mtx)
  gene.names = rownames(tissue.mtx.md)
  goi.rownum = match(feature, gene.names, nomatch = NA_integer_, incomparables = NULL)
  goi.values = as.numeric(tissue.mtx.md[goi.rownum, colnames(tissue.mtx.md)])
  goi.cor = apply(tissue.mtx.md, 1, function(n){cor(goi.values, n)})
  goi.cor = as.matrix(goi.cor)
  goi.cor.sorted = goi.cor[order(-goi.cor), , drop = FALSE]
  tophits = rownames(head(goi.cor.sorted, nhits))
  return(tophits)
}

##### cor_analysis #####

cor_analysis = function(tissue, goi.names, nhits){
  tissue.mtx = tissue[["RNA"]]@data
  tissue.mtx.md = as.matrix(tissue.mtx)
  gene.names = rownames(tissue.mtx.md)
  cell.names = colnames(tissue.mtx.md)
  goi.rownums = na.omit(match(goi.names, gene.names, nomatch = NA_integer_, incomparables = NULL))
  goi.values.mtx = matrix(0, nrow = length(cell.names), ncol = length(goi.rownums))
  for(i in 1:length(goi.rownums)) {
    goi.values.mtx[,i] = as.numeric(tissue.mtx.md[goi.rownums[i], colnames(tissue.mtx.md)])
  }
  goi.cor = apply(tissue.mtx.md, 1, function(n){cor(goi.values.mtx, n)})
  goi.cor = na.omit(t(as.matrix(goi.cor)))
  colnames(goi.cor) = c(all.genes[goi.rownums])
  goi.cor.comb = cbind(goi.cor, sum = rowSums(goi.cor), product = apply(goi.cor, 1, prod), mean = rowMeans(goi.cor))
  goi.cor.sorted = goi.cor.comb[order(-goi.cor.comb[,"mean"]), , drop = FALSE]
  topnhits = rownames(head(goi.cor.sorted, nhits))
  features = c(topnhits)
  return(features)
}

##### find_zone_markers #####

find_zone_markers = function(tissue, ident.1 = NULL, ident.2 = NULL){
  if(is.null(ident.1))
    ident.1 = "Z3"
  if(is.null(ident.2))
    ident.2 = "Z1"
  Idents(tissue) = "orig.ident"
  zone3.markers <- FindMarkers(tissue, ident.1, ident.2, min.pct = 0.25)
  return(zone3.markers)
}

##### zone_avg #####

zone_avg = function(tissue, idents = NULL) {
  if(is.null(idents)) 
    Idents(tissue) = "orig.ident"
  orig.levels <- levels(tissue)
  Idents(tissue) <- gsub(pattern = " ", replacement = "_", x = Idents(tissue))
  orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
  levels(tissue) <- orig.levels
  zone.averages <- AverageExpression(tissue, return.seurat = TRUE)
  return(zone.averages)
}

##### cell_annotation #####

cell_annotation = function(tissue, cluster.mtx){
  cellype.metadata = sapply(as.matrix(tissue@meta.data$seurat_clusters), function(x){
    for(i in 1:ncol(cluster.mtx)) {
      if(x %in% cluster.mtx[,i]){
        return(colnames(cluster.mtx)[i])
      }
    }
  })
  tissue = AddMetaData(tissue, metadata = cellype.metadata[1:length(cellype.metadata)], col.name = "cell_type")
  return(tissue)
}

##### make_marker_list #####

make_marker_list = function(j, vsZone){
  for(m in 1:length(types.tissue)){
    maker.list = rownames(eval(as.name(paste(types.tissue[m], zones[1,j], "vs", zones[vsZone,j], sep = "."))))
    assign(paste("marker.list", m, sep = "."), maker.list)
  }
  return(list(marker.list.1, marker.list.2, marker.list.3, marker.list.4))
}

##### retrieve_gene_data #####

retrieve_gene_data = function(shared.markers, vsZone, mode = NULL){
  if(is.null(mode))
    mode = FALSE
  if(mode == FALSE){
    for(l in 1:length(shared.markers)){
      for(m in 1:length(types.tissue)){
        gene.tissue = eval(as.name(paste(types.tissue[m], zones[1,j], "vs", zones[vsZone,j], sep = ".")))
        gene.tissue = gene.tissue[shared.markers[l],]
        rownames(gene.tissue) = paste(types.tissue[m], shared.markers[l], sep = ".")
        assign(paste("gene", types.tissue[m], sep = "."), t(gene.tissue))
        rm(gene.tissue)
      }
      marker.mtx = cbind(gene.nephron, gene.ureteric, gene.immune, gene.vas)
      assign(paste("gene.data", l, sep = "."), marker.mtx)
      if(l==1){
        gene.data = gene.data.1
      }
      else{
        gene.data = cbind(gene.data, eval(as.name(paste("gene.data", l, sep = '.'))))
      }
    }
    return(gene.data)
  }
  if(mode == TRUE){
    data.tissue = eval(as.name(paste(types.tissue[m], zones[1,j], "vs", zones[vsZone,j], sep = ".")))
    gene.data = data.tissue[shared.markers,]
    return(gene.data)
  }
}

##### not_tissue ####

not.tissue = function(j,m,vsZone){
  not.tissue.names = types.tissue[-m]
  for(i in 1:length(not.tissue.names)){
    genes.not.tissue = rownames(eval(as.name(paste(not.tissue.names[i], zones[1,j], "vs", zones[vsZone,j], sep = "."))))
    assign(paste("genes.not.tissue", toString(i), sep = "."), genes.not.tissue)
  }
  genes.not.tissue = c(genes.not.tissue.1, genes.not.tissue.2, genes.not.tissue.3)
  return(genes.not.tissue)
}

##### determine_dims #####

determine_dims = function(tissue){
  c1 = log(2.5, base = exp(1))/7 
  c2 = (36*log(2, base =  exp(1))+7*log(3, base = exp(1))+13*log(5, base = exp(1)))/(log(2/5, base = exp(1)))
  x = nrow(tissue@meta.data)
  return(round((log(x, exp(1))/c1)+c2))
}

##### find_max_set #####

find_max_set = function(tissue, i){
  for(n in 1:ncol(marker.mtx)){
    marker.list = names(which(marker.mtx[,n] == 1))
    setAVG = mean(tissue.gene.mtx[c(marker.list),i+1])
    if(n == 1){
      setAVG.mtx = setAVG
    }
    if(n != 1){
      setAVG.mtx = cbind(setAVG.mtx, setAVG)
    }
  }
  colnames(setAVG.mtx) = c(1:n)
  return(max.col(setAVG.mtx))
}

##### slingshot_transform ####

slingshot_transform = function(tissue, 
                               largeData = NULL, 
                               label = NULL, 
                               dim.red = NULL, 
                               return.seurat = NULL) {
  
  if(is.null(largeData)){largeData = FALSE}
  if(is.null(label)){label = 'seurat_clusters'}
  if(is.null(dim.red)){dim.red = 'UMAP'}
  if(is.null(return.seurat)){return.seurat = FALSE}
  
  if(length(tissue@reductions) == 0){
    tissue = umap_tsne(PCA_run(tissue), dims = determine_dims(tissue))}
  
  if(!largeData){
    tissue.sce <- slingshot(as.SingleCellExperiment(tissue), 
                            clusterLabels = label, 
                            reducedDim = dim.red)}
  if(largeData){
    tissue.sce <- slingshot(as.SingleCellExperiment(tissue), 
                            clusterLabels = label, 
                            reducedDim = dim.red,
                            approx_points = 5)}
  
  if(!return.seurat){
    rd <- reducedDims(tissue.sce)[[dim.red]]
    attributes(tissue.sce)$rd = rd
    cl <- tissue.sce[[label]]
    attributes(tissue.sce)$cl = cl
    
    lin1 <- getLineages(rd, cl, start.clus = names(table(tissue[[label]])[1]))
    attributes(tissue.sce)$lin1 = lin1
    
    if(largeData)
    {crv1 <- getCurves(lin1)}#, approx_points = 5)}
    if(!largeData)
    {crv1 <- getCurves(lin1)}
    
    attributes(tissue.sce)$crv1 = crv1
    
    return(tissue.sce)
  }
  
  if(return.seurat){
    tissue$Pseudotime = tissue.sce$slingPseudotime_1
    return(tissue)
  }
}

##### retr_cd #####

retr_cd = function(tissue.sce, datatype){
  return(attributes(tissue.sce)[[dataype]])
}

##### find_tissue_markers #####

find_tissue_markers = function(tissue, ident.1 = NULL, ident.2 = NULL){
  if(is.null(ident.1))
    print("missing input: 'ident.1'")
  if(is.null(ident.2))
    print("missing input: 'ident.2'")
  Idents(tissue) = "Lineages"
  tissue.markers <- FindMarkers(tissue, ident.1, ident.2, min.pct = 0.25)
  return(tissue.markers)
}

##### find_cluster_markers #####

find_cluster_markers = function(tissue, ident.1 = NULL, ident.2 = NULL){
  if(is.null(ident.1))
    print("missing input: 'ident.1'")
  if(is.null(ident.2))
    print("missing input: 'ident.2'")
  Idents(tissue) = "seurat_clusters"
  tissue.markers <- FindMarkers(tissue, ident.1, ident.2, min.pct = 0.25)
  return(tissue.markers)
}

##### knn_pred_tissue #####

knn_tiss_pred = function(tissue, gene.mtx, markers = NULL, k = NULL){
  if(is.null(markers)){
    warning('no markers given - please use markers for more accurate processing')
  } else {
    gene.mtx = gene.mtx[,rownames(markers)]
  }
  if(is.null(k) || k == 1){
    warning('k = 1: no probability values will be assigned')
    prob = FALSE
    k = 1
  } else{
    prob = TRUE
  }
  
  tr <<- sample(nrow(gene.mtx), nrow(gene.mtx))
  nw <<- sample(nrow(gene.mtx), nrow(gene.mtx))
  n.pred = ncol(gene.mtx)
  knnres = knn(gene.mtx[tr, -n.pred], 
               gene.mtx[nw, -n.pred], 
               gene.mtx$Lineages[tr], 
               k = k,
               prob = prob)
  gene.mtx$pred_Lineages[nw] = knnres
  return(gene.mtx)
}

##### knn_zone_pred #####

knn_zone_pred = function(tissue, gene.mtx.input, markers = NULL, k = NULL){
  if(is.null(markers)){
    warning('no markers given - please use markers for more accurate processing')
  } else {
    gene.mtx.input = gene.mtx.input[, c(markers, "orig.ident")]
  }
  if(is.null(k) || k == 1){
    warning('k = 1: no probability values will be assigned')
    prob = FALSE
    k = 1
  } else {
    prob = TRUE
  }
  
  gene.mtx.input$Lineages = NULL
  message(paste('Running knn on', fullNames.tissue[i], sep = ' '))
  tr <<- sample(nrow(gene.mtx.input), round(nrow(gene.mtx.input)*0.6))
  nw <<- sample(nrow(gene.mtx.input), nrow(gene.mtx.input))
  n.pred <<- ncol(gene.mtx.input)
  message(paste('Using', (n.pred - 1), 'markers', sep = ' '))
  knnres = knn(gene.mtx.input[tr, -n.pred], 
               gene.mtx.input[nw, -n.pred], 
               gene.mtx.input$orig.ident[tr], 
               k = k,
               prob = prob)
  gene.mtx.input$pred_ident[nw] = knnres
  return(gene.mtx.input)
}

##### svm_pred #####

svm_pred = function(i, k = NULL, tuneLength = NULL){
  if(is.null(k)){
    k = 5
    warning(paste('was not given: autoset to', k, sep = ' '))
  }
  if(is.null(tuneLength)){
    tuneLength = 10
    warning(paste('was not given: autoset to', tuneLength, sep = ' '))
  }
  
  if(!exists(paste(types.tissue[i], 'ret.mtx', sep = '.'))){
    message(paste('Creating missing matrix', paste(types.tissue[i], 'ret.mtx', sep = '.', collapse = '\n'), sep = ' '))
    markers = eval(as.name(paste('markers', types.tissue[i], sep = '.')))
    topN.tiss.zone <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    topN = rbind(topN.zone, topN.tiss.zone)
    
    assign(paste(types.tissue[i], 'ret.mtx', sep = '.'),
           subset(gene.mtx, Lineages == fullNames.tissue[i])[,c(topN$gene, 'orig.ident')], 
           envir = .GlobalEnv)
  }
    
  message(paste('Running support vector machine model on', fullNames.tissue[i], '...', sep = ' ', collapse = '\n'))
  message(paste('Creating', k, 'folds', sep = ' '))
  myFolds <<- createFolds(eval(as.name(paste(types.tissue[i], 'ret.mtx', sep = '.')))$orig.ident, k = k)

  message(paste('Training control', sep = ' '))
  myControl <<- trainControl(
    classProb = TRUE,
    verboseIter = FALSE,
    savePredictions = TRUE,
    index = myFolds
  )
  
  message(paste('Training support vector machine model with tune length:', tuneLength, sep = ' '))
  svm_model <- train(orig.ident ~ .,
                     eval(as.name(paste(types.tissue[i], 'ret.mtx', sep = '.'))),
                     method = "svmRadial",
                     tuneLength = tuneLength,
                     trControl = myControl)
  
  rm(myFolds, myControl)
  return(svm_model)
}

##### draw_confusion_matrix #####

draw_confusion_matrix <- function(cm) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 435, 'Class1', cex=1.2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 435, 'Class2', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'Class1', cex=1.2, srt=90)
  text(140, 335, 'Class2', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  

##### BioMart #####

convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}