if (!requireNamespace("enrichR", quietly = TRUE)) {
  install.packages("enrichR")
}
library(enrichR)

data <- read.csv("ORA_results.csv", stringsAsFactors = FALSE)
genes <- unique(data$level_1)
database <- c("KEGG_2021_Human", "GO_Biological_Process_2021")
enriched <- enrichr(genes, database)

write.csv(enriched[["KEGG_2021_Human"]], "KEGG_enrichment.csv", row.names = FALSE)
write.csv(enriched[["GO_Biological_Process_2021"]], "GO_BP_enrichment.csv", row.names = FALSE)
