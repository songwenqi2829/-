library(dplyr)
library(limma)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# 读取并预处理临床信息数据
clinical <- read.csv(file = '~/pj/rawdata/clinical.tsv',sep = '\t')
rownames(clinical) <- clinical$attrib_name
clinical$attrib_name <- NULL
clinical <- t(clinical) %>% as.data.frame()
clinical$sample <- rownames(clinical)
clinical <- subset(clinical,Tumor.Status == 'With tumor')
group <- clinical$Integrated.Phenotype
names(group) <- clinical$sample
clinical$Mutation_Phenotype <- gsub('-','_',clinical$Mutation_Phenotype)

# 读取并预处理RNA测序数据
rna_exp <- read.csv(file = '~/pj/rawdata/RNAseq.tsv',sep = '\t')
rownames(rna_exp) <- rna_exp$attrib_name
rna_exp$attrib_name <- NULL

samples <- intersect(colnames(rna_exp),rownames(clinical))
rna_exp <- rna_exp[,samples]
clin_rna <- clinical[samples,c('Mutation_Phenotype','CIN','Integrated.Phenotype')]
colnames(clin_rna) <- c('MUT','CIN','EMT')

# 读取并预处理蛋白组数据
pro_exp <- read.csv(file = '~/pj/rawdata/proteomics.tsv',sep = '\t')
rownames(pro_exp) <- pro_exp$attrib_name
pro_exp$attrib_name <- NULL

samples <- intersect(colnames(pro_exp),rownames(clinical))
pro_exp <- pro_exp[,samples]
clin_pro <- clinical[samples,c('Mutation_Phenotype','CIN','Integrated.Phenotype')]
colnames(clin_pro) <- c('MUT','CIN','EMT')

# 定义差异分析函数
compute_diff <- function(exp, group, compare) {
  common_sample <- intersect(colnames(exp), names(group))
  exp   <- exp[, common_sample, drop = FALSE]
  group <- group[common_sample]
  
  comp <- unlist(strsplit(compare, "-"))
  if (length(comp) != 2) {stop("compare must be in the format 'GroupA-GroupB'")}
  group1 <- comp[1];group2 <- comp[2]
  
  keep <- group %in% c(group1, group2)
  exp   <- exp[, keep, drop = FALSE]
  group <- factor(group[keep], levels = c(group1, group2))
  
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(exp, design)
  contrast.matrix <- makeContrasts(contrasts = paste0(group1, "-", group2),levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2,coef = 1,number = Inf,adjust.method = "BH",sort.by = "P")
  
  res$gene <- rownames(res)
  res$comparison <- compare
  return(res)
}

# 差异分析
group <- clin_rna$MUT;names(group) <- rownames(clin_rna)
rna_Hyper_other <- compute_diff(rna_exp,group = group,compare = 'MSI_H-MSS')
group <- clin_rna$CIN;names(group) <- rownames(clin_rna)
rna_Epi_other <- compute_diff(rna_exp,group = group,compare = 'CIN.high-CIN.low')
group <- clin_rna$EMT;names(group) <- rownames(clin_rna)
rna_EMT_other <- compute_diff(rna_exp,group = group,compare = 'Epithelial-EMT')

group <- clin_pro$MUT;names(group) <- rownames(clin_pro)
pro_Hyper_other <- compute_diff(pro_exp,group = group,compare = 'MSI_H-MSS')
group <- clin_pro$CIN;names(group) <- rownames(clin_pro)
pro_Epi_other <- compute_diff(pro_exp,group = group,compare = 'CIN.high-CIN.low')
group <- clin_pro$EMT;names(group) <- rownames(clin_pro)
pro_EMT_other <- compute_diff(pro_exp,group = group,compare = 'Epithelial-EMT')

# 绘制图像
data <- get('rna_EMT_other')
name <- 'EMT'

data$diff_class <- "Not significant"
data$diff_class[data$P.Value < 0.05 & data$logFC > 1] <- paste0('Up in ',name)
data$diff_class[data$P.Value < 0.05 & data$logFC < -1] <- paste0('Down in ',name)

# 定义颜色
cols <- setNames(
  c("firebrick", "steelblue", "grey70"),
  c(paste0("Up in ", name),
    paste0("Down in ", name),
    "Not significant")
)

# 绘图
ggplot(data, aes(x = logFC, y = -log10(P.Value), color = diff_class)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_manual(values = cols) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  theme_classic() +
  labs(
    title = paste0("Differential expression: ", name, " vs Others"),
    x = paste0("log2 Fold Change (", name, "/Other)"),
    y = "-log10(P value)",
    color = ""
  ) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top")


# 富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

result_list <- c('rna_EMT_other','rna_Hyper_other','rna_Epi_other',
                 'pro_EMT_other','pro_Hyper_other','pro_Epi_other')

for (res_name in result_list) {
  
  # 获取数据
  res <- get(res_name)
  
  # 上调基因
  sig_up <- rownames(subset(res, P.Value < 0.05 & logFC > 1))
  if (length(sig_up) > 0) {
    gene_df_up <- bitr(sig_up, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    entrez_up <- gene_df_up$ENTREZID
    
    ego_up <- enrichGO(gene = entrez_up,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)
    
    # 赋值给变量，例如 rna_EMT_other_up_go
    assign(paste0(res_name, "_up_go"), ego_up)
    
    # 保存结果
    if (nrow(ego_up) > 0) {
      write.csv(as.data.frame(ego_up), file = paste0(res_name, "_up_BP_GO.csv"), row.names = FALSE)
    }
  }
  
  # 下调基因
  sig_down <- rownames(subset(res, P.Value < 0.05 & logFC < -1))
  if (length(sig_down) > 0) {
    gene_df_down <- bitr(sig_down, fromType = "SYMBOL",
                         toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    entrez_down <- gene_df_down$ENTREZID
    
    ego_down <- enrichGO(gene = entrez_down,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         readable = TRUE)
    
    # 赋值给变量，例如 rna_EMT_other_down_go
    assign(paste0(res_name, "_down_go"), ego_down)
    
    # 保存结果
    if (nrow(ego_down) > 0) {
      write.csv(as.data.frame(ego_down), file = paste0('~/pj/result/',res_name, "_down_BP_GO.csv"), row.names = FALSE)
    }
  }
  
  cat("Finished:", res_name, "\n")
}

# 绘制图像
barplot(rna_Hyper_other_up_go,showCategory = 7)



