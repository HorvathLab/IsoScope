.libPaths(c("/home/azka.sajjad/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
library(DBI)
library(RSQLite)

db <- dbConnect(SQLite(), "out/ccle/quants.db")
smp <- dbReadTable(db, "samples")
hit <- apply(smp, 1, function(r) any(grepl("^(BT-?549|MDA-?MB-?231|MCF-?7)(_|$)", r, ignore.case = TRUE)))
smp <- smp[hit, ]
print(smp)

qcols <- dbListFields(db, "quants")
idcol <- if (all(smp$sample_id %in% qcols)) "sample_id" else "sample_label"
stopifnot(all(smp[[idcol]] %in% qcols))

tx <- c("ENST00000295956.9", "ENST00000358537.7", "ENST00000369329.8", "ENST00000369332.7")
sql <- sprintf("SELECT gene, transcript_id, %s FROM quants WHERE transcript_id IN (%s)",
               paste0('"', smp[[idcol]], '"', collapse = ", "),
               paste0("'", tx, "'", collapse = ", "))
q <- dbGetQuery(db, sql)
dbDisconnect(db)
stopifnot(nrow(q) == 4)

m <- t(as.matrix(q[, smp[[idcol]], drop = FALSE]))
colnames(m) <- q$transcript_id
out <- data.frame(sample = smp[[idcol]], round(m[, tx], 2), check.names = FALSE)
out$FLNB_skip_frac <- round(100 * out$ENST00000358537.7 / (out$ENST00000358537.7 + out$ENST00000295956.9), 1)
print(out, row.names = FALSE)
