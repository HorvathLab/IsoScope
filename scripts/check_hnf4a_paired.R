library(DBI)
library(RSQLite)

con <- dbConnect(SQLite(), "shiny_data/hcc/quants.db")
smp <- dbGetQuery(con, "SELECT * FROM samples")
q <- dbGetQuery(con, "SELECT * FROM quants WHERE gene = 'HNF4A'")
dbDisconnect(con)

p1 <- c("ENST00000316099.10", "ENST00000415691.2", "ENST00000443598.6")
p2 <- c("ENST00000316673.9", "ENST00000894460.1", "ENST00000457232.5", "ENST00000609795.5")
cols <- intersect(smp$sample_label, names(q))
tot <- colSums(q[, cols])
p1s <- 100 * colSums(q[q$transcript_id %in% p1, cols]) / tot
p2s <- 100 * colSums(q[q$transcript_id %in% p2, cols]) / tot

smp$grp <- ifelse(grepl("adj", smp$group_label, ignore.case = TRUE), "adj", "tum")
pairs <- intersect(smp$subject_id[smp$grp == "adj"], smp$subject_id[smp$grp == "tum"])
pm <- function(v, p, g) mean(v[intersect(smp$sample_label[smp$subject_id == p & smp$grp == g], cols)])

out <- data.frame(
  patient = pairs,
  gene_tpm_adj = sapply(pairs, function(p) round(pm(tot, p, "adj"), 1)),
  gene_tpm_tum = sapply(pairs, function(p) round(pm(tot, p, "tum"), 1)),
  P1_adj = sapply(pairs, function(p) round(pm(p1s, p, "adj"), 1)),
  P1_tum = sapply(pairs, function(p) round(pm(p1s, p, "tum"), 1)),
  P2_adj = sapply(pairs, function(p) round(pm(p2s, p, "adj"), 1)),
  P2_tum = sapply(pairs, function(p) round(pm(p2s, p, "tum"), 1))
)
print(out, row.names = FALSE)
cat("\nMedian P1 share adj -> tum:", median(out$P1_adj), "->", median(out$P1_tum),
    " down in", sum(out$P1_tum < out$P1_adj), "/", nrow(out),
    " p =", wilcox.test(out$P1_tum, out$P1_adj, paired = TRUE)$p.value, "\n")
cat("Median P2 share adj -> tum:", median(out$P2_adj), "->", median(out$P2_tum),
    " up in", sum(out$P2_tum > out$P2_adj), "/", nrow(out),
    " p =", wilcox.test(out$P2_tum, out$P2_adj, paired = TRUE)$p.value, "\n")
