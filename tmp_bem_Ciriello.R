library(pheatmap)

input_fold <- "/group/iorio/Datasets/tumour_binary_Ciriello/"
object_to_load <- list.files(input_fold)

obj <- list()
for(id_name in seq_len(length(object_to_load))) {
  
  obj[[id_name]] <- readRDS(paste0(input_fold, object_to_load[id_name]))
  obj[[id_name]]$M$M$total <- obj[[id_name]]$M$M$missense + obj[[id_name]]$M$M$truncating
  
}

id <- which(sapply(obj, function(x) "NSCLC" %in% x$sample.class))
samples_id <- lapply(id, function(x) names(which(obj[[x]]$sample.class == "NSCLC")))
bem <- lapply(id, function(x) obj[[x]]$M$M$total[,samples_id[[x]]])
bem <- do.call(cbind, bem)
bem <- bem[rowSums(bem) > 2, ]
pheatmap(t(bem)[, c("TP53", "KRAS")])
