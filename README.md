# Single-cell-RNA-sequencing
# Step 1: Load required libraries
library(Seurat)

# Step 2: Load and preprocess the data
data <- Read10X("path/to/your/data")  # Replace with the path to your data
seurat_obj <- CreateSeuratObject(counts = data)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# Step 3: Perform dimensionality reduction
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)

# Step 4: Visualize the results
seurat_obj <- RunUMAP(seurat_obj)
DimPlot(seurat_obj)

# Step 5: Identify cluster-specific marker genes
seurat_obj <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Step 6: Explore cluster identities
Idents(seurat_obj) <- "cluster"
FeaturePlot(seurat_obj, features = c("gene1", "gene2", "gene3"), cols = c("blue", "red", "green"))

# Step 7: Perform differential gene expression analysis
seurat_obj <- FindMarkers(seurat_obj, ident.1 = 0, ident.2 = 1, min.pct = 0.1, logfc.threshold = 0.25)

# Step 8: Perform gene set enrichment analysis
enrichment_results <- RunGSVA(seurat_obj, genesets = your_genesets, verbose = FALSE)

# Step 9: Integrate multiple datasets (if applicable)
seurat_obj_integrated <- FindIntegrationAnchors(object.list = list(seurat_obj1, seurat_obj2))
seurat_obj_integrated <- IntegrateData(anchorset = seurat_obj_integrated)
seurat_obj_integrated <- ScaleData(seurat_obj_integrated)

# Step 10: Perform batch correction (if applicable)
seurat_obj_integrated <- FindVariableFeatures(seurat_obj_integrated)
seurat_obj_integrated <- ScaleData(seurat_obj_integrated)
seurat_obj_integrated <- RunPCA(seurat_obj_integrated)
seurat_obj_integrated <- FindNeighbors(seurat_obj_integrated)
seurat_obj_integrated <- FindClusters(seurat_obj_integrated)

# Step 11: Visualize integrated results
seurat_obj_integrated <- RunUMAP(seurat_obj_integrated)
DimPlot(seurat_obj_integrated, group.by = "dataset")

# Step 12: Identify cluster-specific marker genes in integrated data
seurat_obj_integrated <- FindAllMarkers(seurat_obj_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Step 13: Perform differential gene expression analysis in integrated data
seurat_obj_integrated <- FindMarkers(seurat_obj_integrated, ident.1 = 0, ident.2 = 1, min.pct = 0.1, logfc.threshold = 0.25)

# Step 14: Perform gene set enrichment analysis in integrated data
enrichment_results_integrated <- RunGSVA(seurat_obj_integrated, genesets = your_genesets, verbose = FALSE)
