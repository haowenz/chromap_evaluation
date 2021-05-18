# chromap_evaluation
Evaluation scripts for Chromap manuscript on different data are in their own folders.

### Simulated data
We used the [paftools](https://github.com/lh3/minimap2/tree/master/misc) to evaluate the alignment results from Chromap and other aligners.

### ChIP-seq 

### Hi-C 
All the Hi-C evaluation scripts are in its folder.

### scATAC-seq
We used the R command 

	rds<-readRDS("MAESTRO.rds")
	write.table(cbind(rds$ATAC@meta.data, Embeddings(rds$ATAC@reductions$umap)), "new.tsv", sep="\t", quote=FALSE) 

to extract clustering information from MAESTRO.
