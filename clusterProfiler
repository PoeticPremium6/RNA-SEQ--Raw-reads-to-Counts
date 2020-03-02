Args <- commandArgs(T)
de_file = Args[1]
library("clusterProfiler")
library(org.Mm.eg.db)

read.table(file=de_file,header=T,row.names=1)  ->de
up <- row.names(de)[de$logFC>= 1&de$FDR < 0.05]
#up <-up$Name
head(up)
down <- row.names(de)[de$logFC<= -1&de$FDR < 0.05]
#down <-down$Name
head(down)
ego_up_all <- enrichGO(gene          =up,
                universe      = row.names(de),
                     OrgDb="org.Mm.eg.db",
                ont           = "ALL",
				keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
ego_down_all <- enrichGO(gene          =down,
                universe      = row.names(de),
                     OrgDb="org.Mm.eg.db",
                ont           = "ALL",
				keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)




ego_up_cc <- enrichGO(gene          =up,
                universe      = row.names(de),
                     OrgDb="org.Mm.eg.db",
                ont           = "CC",
				keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
               
ego_up_bp <- enrichGO(gene          = up,
                universe      = row.names(de),
                OrgDb="org.Mm.eg.db",
                ont           = "BP",
				keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
               
        readable      = TRUE)
ego_up_mf <- enrichGO(gene          = up,
                universe      = row.names(de),
                    OrgDb="org.Mm.eg.db",
                ont           = "MF",
				keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)

ego_down_cc <- enrichGO(gene          =down,
                universe      = row.names(de),
                     OrgDb="org.Mm.eg.db",
                ont           = "CC",
				keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
ego_down_bp <- enrichGO(gene          = down,
                universe      = row.names(de),
                OrgDb="org.Mm.eg.db",
                ont           = "BP",
				keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
               
        readable      = TRUE)
ego_down_mf <- enrichGO(gene          = down,
                universe      = row.names(de),
                    OrgDb="org.Mm.eg.db",
                ont           = "MF",
				keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)

merge_result(list(Up_regulated=ego_up_bp, Down_regulated=ego_down_bp)) ->up_down_bp
dotplot(up_down_bp,showCategory=10,font.size=18, includeAll = TRUE) + scale_size(range=c(2, 15))
dotplot(up_down_bp,showCategory=10,font.size=18, includeAll = TRUE) + scale_size(range=c(2, 15))
dev.copy2pdf(device='x11',file="GO_bp_dot.pdf")
write.table(ego_up_all,file="GO_up_detail.xls",quote=F,sep="\t")
write.table(ego_down_all,file="GO_down_detail.xls",quote=F,sep="\t")
#kegg
up_entrez <- bitr(up, fromType = 'ENSEMBL', toType = c('SYMBOL', 'ENTREZID'), OrgDb = 'org.Mm.eg.db')
down_entrez <- bitr(down, fromType = 'ENSEMBL', toType = c('SYMBOL', 'ENTREZID'), OrgDb = 'org.Mm.eg.db')
#de_entrez <- bitr(de, fromType = 'ENSEMBL', toType = c('SYMBOL', 'ENTREZID'), OrgDb = 'org.Mm.eg.db')
de_entrez <- bitr(row.names(de), fromType = 'ENSEMBL', toType = c('SYMBOL', 'ENTREZID'), OrgDb = 'org.Mm.eg.db')
up_down_entrez <- rbind(up_entrez,down_entrez)
gene_matrix <- de$logFC
names(gene_matrix) <- de_entrez$ENTREZID 
kegg_enrich <- enrichKEGG(gene = up_down_entrez$ENTREZID,
                          organism = 'mouse',
 universe = de_entrez$ENTREZID,
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)

dotplot(kegg_enrich)
dev.copy2pdf(device='x11',file="KEGG_DE_dot.pdf")
dotplot(kegg_enrich,showCategory=37)
dev.copy2pdf(device='x11',file="KEGG_all_dot.pdf")
write.table(file="KEGG.summary.xls",kegg_enrich,quote=F,sep="\t")
kegg_enrich1 <-setReadable(kegg_enrich,OrgDb="org.Mm.eg.db",keytype="ENTREZID")
write.table(file="KEGG.summary1.xls",kegg_enrich1,quote=F,sep="\t")
###view
mmu04650 <- pathview(gene.data = gene_matrix, 
         pathway.id = "04650", 
         species = "mouse",
 limit      = list(gene=max(abs(gene_matrix)), cpd=1))
mmu05340 <- pathview(gene.data = gene_matrix, 
		 pathway.id = "05340", 
		 species = "mouse",
 limit      = list(gene=max(abs(gene_matrix)), cpd=1))

pathview(gene.data = gene_matrix,
                 pathway.id = "mmu05202",
                 species = "mouse",
 limit      = list(gene=max(abs(gene_matrix)), cpd=1))
 pathview(gene.data = gene_matrix,
                 pathway.id = "mmu05222",
                 species = "mouse",
 limit      = list(gene=max(abs(gene_matrix)), cpd=1))
 pathview(gene.data = gene_matrix,
                 pathway.id = "mmu05010",
                 species = "mouse",
 limit      = list(gene=max(abs(gene_matrix)), cpd=1))
 pathview(gene.data = gene_matrix,
                 pathway.id = "mmu04940",
                 species = "mouse",
 limit      = list(gene=max(abs(gene_matrix)), cpd=1))
 pathview(gene.data = gene_matrix,
                 pathway.id = "mmu04640",
                 species = "mouse",
 limit      = list(gene=max(abs(gene_matrix)), cpd=1))

######DOSE only for human
library(DOSE)
ego_up_DO <- enrichDO(gene          =up,
                universe      = row.names(de),
                     OrgDb="org.Mm.eg.db",
                ont           = "DO",
				keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
ego_down_DO <- enrichDO(gene          =down,
                universe      = row.names(de),
                     OrgDb="org.Mm.eg.db",
                ont           = "DO",
				keyType       = 'ENSEMBL',
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)


merge_result(list(Up_regulated=ego_up_DO, Down_regulated=ego_down_DO)) ->up_down_DO
dotplot(up_down_DO,showCategory=10,font.size=18, includeALL = TRUE) + scale_size(range=c(2, 15))
dotplot(up_down_DO,showCategory=10,font.size=18, includeALL = TRUE) + scale_size(range=c(2, 15))
dev.copy2pdf(device='x11',file="GO_DO_dot.pdf")
write.table(ego_up_DO,file="DO_up_detail.xls",quote=F,sep="\t")
write.table(ego_down_DO,file="DO_down_detail.xls",quote=F,sep="\t")

