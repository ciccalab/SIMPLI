library(Seurat)
library(irlba)

RunMultiCCA_Changed <- function(object.list, genes.use, add.cell.ids = NULL, niter = 25, num.ccs = 1, standardize = TRUE) 
{
  set.seed(42)
  if (length(object.list) < 3){
    stop("Must give at least 3 objects/matrices for MultiCCA")
  }
  mat.list <- list()
  if (class(object.list[[1]]) == "seurat") {
    if (missing(x = genes.use)) {
      genes.use <- c()
      for (obj in object.list) {
        genes.use <- c(genes.use, obj@var.genes)
      }
      genes.use <- unique(genes.use)
      if (length(x = genes.use) == 0) {
        stop("No variable genes present. Run MeanVarPlot and retry")
      }
    }
    for (i in 1:length(object.list)) {
      mat.list[[i]] <- object.list[[i]]@scale.data[genes.use, ]
    }
  }
  else {
    stop("input data not Seurat objects")
  }
  if (!missing(add.cell.ids)) {
    if (length(add.cell.ids) != length(object.list)) {
      stop("add.cell.ids must have the same length as object.list")
    }
    object.list <- lapply(seq_along(object.list), function(i) {
      RenameCells(object = object.list[[i]], add.cell.id = add.cell.ids[i])
    })
  }
  names.list <- lapply(object.list, slot, name = "cell.names")
  names.intersect <- Reduce(intersect, names.list)
  if (length(names.intersect) > 0) {
    stop("duplicate cell names detected, please set 'add.cell.ids'")
  }
  num.sets <- length(mat.list)
  if (standardize) {
    for (i in 1:num.sets) {
      mat.list[[i]] <- Seurat:::Standardize(mat.list[[i]], display_progress = F)
    }
  }
  ws <- list()
  for (i in 1:num.sets) {
    ws[[i]] <- irlba(mat.list[[i]], nv = num.ccs)$v[, 1:num.ccs,  drop = F]
  }
  ws.init <- ws
  ws.final <- list()
  cors <- NULL
  for (i in 1:length(ws)) {
    ws.final[[i]] <- matrix(0, nrow = ncol(mat.list[[i]]),  ncol = num.ccs)
  }
  for (cc in 1:num.ccs) {
    print(paste0("Computing CC ", cc))
    ws <- list()
    for (i in 1:length(ws.init)) {
      ws[[i]] <- ws.init[[i]][, cc]
    }
    cur.iter <- 1
    crit.old <- -10
    crit <- -20
    storecrits <- NULL
    while (cur.iter <= niter && abs(crit.old - crit)/abs(crit.old) > 
           0.001 && crit.old != 0) {
      crit.old <- crit
      crit <- Seurat:::GetCrit(mat.list, ws, num.sets)
      storecrits <- c(storecrits, crit)
      cur.iter <- cur.iter + 1
      for (i in 1:num.sets) {
        ws[[i]] <- Seurat:::UpdateW(mat.list, i, num.sets, ws, ws.final)
      }
    }
    for (i in 1:length(ws)) {
      ws.final[[i]][, cc] <- ws[[i]]
    }
    cors <- c(cors, Seurat:::GetCors(mat.list, ws, num.sets))
  }
  results <- list(ws = ws.final, ws.init = ws.init, num.sets = num.sets, 
                  cors = cors)
  combined.object <- object.list[[1]]
  for (i in 2:length(object.list)) {
    combined.object <- MergeSeurat(object1 = combined.object, 
                                   object2 = object.list[[i]], do.scale = F, do.center = F, 
                                   do.normalize = F, min.cells = -1, min.genes = -1, is.expr = -1)
  }
  combined.object <- NormalizeData(combined.object)
  combined.object@meta.data$orig.ident <- combined.object@cell.names
  combined.object <- ScaleData(object = combined.object)
  combined.object@scale.data[is.na(x = combined.object@scale.data)] <- 0
  combined.object@var.genes <- genes.use
  cca.data <- results$ws[[1]]
  for (i in 2:length(object.list)) {
    cca.data <- rbind(cca.data, results$ws[[i]])
  }
  rownames(cca.data) <- colnames(combined.object@data)
  cca.data <- apply(cca.data, MARGIN = 2, function(x) {
    if (sign(x[1]) == -1) {
      x <- x * -1
    }
    return(x)
  })
  combined.object <- SetDimReduction(object = combined.object, reduction.type = "cca", slot = "cell.embeddings",
    new.data = cca.data)
  combined.object <- SetDimReduction(object = combined.object, reduction.type = "cca", slot = "key", new.data = "CC")
  combined.object <- ProjectDim(object = combined.object, reduction.type = "cca", do.print = FALSE)
  combined.object <- SetDimReduction(object = combined.object, reduction.type = "cca", slot = "gene.loadings",
    new.data = GetGeneLoadings(object = combined.object, reduction.type = "cca", use.full = TRUE, genes.use = genes.use))
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("RunMultiCCA"))]
  parameters.to.store$object.list <- NULL
  combined.object <- Seurat:::SetCalcParams(object = combined.object, calculation = "RunMultiCCA",
    ... = parameters.to.store)
  return(combined.object)
}