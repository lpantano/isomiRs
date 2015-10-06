

.get_mirna_mrna_ma = function(.toy, min_total=0.1){
  .seen = vector()
  .sort_mirna = sort(colSums(.toy>0), decreasing = T)
  .coop = vector("list", length=length(.sort_mirna))
  # Join miRNAs that target common genes
  for (nc1 in 1:length(.sort_mirna)){
    for (nc2 in 1:length(.sort_mirna)){
      if (nc1 != nc2 & sum(.seen==nc2)==0){
        .name1 = names(.sort_mirna)[nc1]
        .name2 = names(.sort_mirna)[nc2]
        .these1 = rownames(.toy)[.toy[,.name1]>0]
        .these2 = rownames(.toy)[.toy[,.name2]>0]
        .min = min(length(.these1), length(.these2))
        .common = length(intersect(.these1, .these2))
        if ( (.common / .min) > 0.50 ){
          .seen = unique(c(.seen , nc2, nc1))
          .coop[[nc1]] = c(.coop[[nc1]], nc2)
        }
      }
    }
  }

  .others = c()
  # miRNA not paired with any other
  for (nc1 in 1:length(.sort_mirna)){
    if (sum(.seen==nc1)==0 ){
      .these1 = rownames(.toy)[.toy[,nc1]>0]
      if ( length(.these1) / nrow(.toy) >= min_total ){
        .coop[[nc1]] = c(.coop[[nc1]], nc1)
      }#else{
      #.others = c(.others, nc1)
      #}
      .seen = c(.seen , nc1)
    }
  }

  # Create binary matrix for clustering
  # Create reduce matrix with coop miRNAs in one column
  .toy_binary <- .toy
  .toy_binary_red = lapply(1:length(.sort_mirna), function(x){
    if ( ! is.null(.coop[[x]]) ){
      .join = unique(c(names(.sort_mirna)[x], names(.sort_mirna)[.coop[[x]]]))
      .name = gsub('hsa-', '', gsub("miR-", "", paste0(.join, collapse = ",")))
      .news = rowSums(.toy_binary[,.join, drop=F])
      if ( sum(.news) >= nrow(.toy_binary) * min_total )
        return(list(v=.news, n=.name, real=.join))
    }
  })
  .used = unique(unlist(lapply(.toy_binary_red, function(x){
    x$real
  })))
  .names = unlist(lapply(.toy_binary_red, function(x){
    x$n
  }))
  .toy_binary_red = do.call(cbind, lapply(.toy_binary_red, function(x){
    x$v
  }))
  colnames(.toy_binary_red) = .names
  .toy_binary_red = as.data.frame(.toy_binary_red)

  # Clustering columns
  .toy_hc = rownames(.toy_binary_red)
  if (ncol(.toy_binary_red) > 1){
    cat("columns", ncol(.toy_binary_red))
    .toy_hc = hclust(dist(t(.toy_binary_red), method='binary'), method="ward.D2")
  }

  # Create categorical matrix for viz porpuse
  .toy_red = .toy_binary_red
  .toy_red[.toy_binary_red==0] = "None"
  .toy_red[.toy_binary_red!=0] = "Target"

  .ann_colors = lapply(colnames(.toy_red), function(x) {c(None="azure", Target="grey")})
  names(.ann_colors) = colnames(.toy_red)
  return(list(binary=.toy_binary_red, category=.toy_red, cols=.ann_colors, clustering=.toy_hc, used =.used))
}



.create_interactome <- function(de_mirna, de_mrna, mapping){
  network_list = list()
  for (contrast in intersect(names(de_mirna, de_mrna)) ){
    genes_de = de_mrna[[contrast]]
    mirna_de = de_mirna[[contrast]]
    target_de = mapping %>% filter(gene %in% names(genes_de))  %>% filter(mirna %in% names(mirna_de))
    idx_mirna = match(target_de$mirna, names(genes_de))
    idx_mrna = match(target_de$gene, names(mirna_de))
    target_de[!is.na(idx_mirna),"mirna_fc"] = mirna_de[idx_mirna[!is.na(idx_mirna)]]
    target_de[!is.na(idx_mrna),"mrna_fc"] = genes_de[idx_mrna[!is.na(idx_mrna)]]
    network = target_de %>% filter(!is.na(mirna_fc) & !is.na(mrna_fc) & mirna_fc * mrna_fc < 0)
    if (dim(network)[1]>0){
      network_list[[contrast]] = network %>% mutate(contrast = contrast)
    }
  }
  do.call(rbind, network_list)
}