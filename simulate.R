#! /usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

nfrags <- as.integer(args[1])
locs <- unlist(strsplit(args[2], ","))

library(Biostrings)
library(polyester)

stopifnot(nfrags %% length(locs) == 0)

params <- list(
    outdir="sims",
    fraglen=250, fragsd=25, readlen=75, paired=TRUE,
    error_model="illumina5", error_rate=0.005,
    seed=123456, total_reads=nfrags, num_reps=1
    )

countmat <- data.frame(
    sample_01=rep.int(nfrags/length(locs), length(locs)),
    row.names=locs
    )

write.table(countmat, file=file.path(params$outdir, 'countmat.tsv'), sep='\t', quote=F)

allseqs <- readDNAStringSet("refs/HML2_extracted.fna")
writeXStringSet(allseqs[rownames(countmat)], "tmp.fasta")

simulate_experiment_countmat(
                    fasta = "tmp.fasta",
                    readmat = as.matrix(countmat),
                    outdir = params$outdir,
                    fraglen = params$fraglen,
                    fragsd = params$fragsd,
                    readlen = params$readlen,
                    error_model = params$error_model,
                    error_rate = params$error_rate,
                    seed = params$seed,
                    )

unlink("tmp.fasta")
