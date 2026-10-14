suppressMessages({library(DBI); library(RSQLite); library(data.table)})
ga <- fread("shiny_data/ccle/group_avg.csv")
r <- ga[gene == "RBP1", .(transcript_id, PANCREAS, STOMACH)]
cat("RBP1 total TPM  pancreas", round(sum(r$PANCREAS), 1), " stomach", round(sum(r$STOMACH), 1), "\n")
r[, `:=`(pct_pancreas = round(100 * PANCREAS / sum(PANCREAS)), pct_stomach = round(100 * STOMACH / sum(STOMACH)))]
print(r[order(-PANCREAS)][1:4])
con <- dbConnect(SQLite(), "shiny_data/ccle/quants.db")
p <- dbGetQuery(con, 'SELECT transcript_id, "HCC1143_BREAST", "HCC1187_BREAST" FROM quants WHERE gene = "PARP1"')
dbDisconnect(con)
ent <- function(x) { x <- x[x > 0]; if (length(x) < 2) return(0); q <- x / sum(x); -sum(q * log2(q)) / log2(length(x)) }
for (cl in c("HCC1143_BREAST", "HCC1187_BREAST")) {
  v <- p[[cl]]; o <- order(-v)[1:3]
  cat(cl, " total", round(sum(v), 1), " entropy", round(ent(v), 2), "\n")
  print(data.frame(tx = p$transcript_id[o], pct = round(100 * v[o] / sum(v))))
}
