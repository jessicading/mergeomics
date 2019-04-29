source("../Mergeomics.R")

job.ssea <- list()
job.ssea$label <- ""
job.ssea$folder <- ""
job.ssea$genfile <- ""
job.ssea$locfile <- ""		
job.ssea$modfile <- ""
job.ssea$inffile <- ""
job.ssea$permtype <- "gene"
job.ssea$nperm <- 10000
job.ssea <- ssea.start(job.ssea)
job.ssea <- ssea.prepare(job.ssea)
job.ssea <- ssea.control(job.ssea)
job.ssea <- ssea.analyze(job.ssea,trim_start=0.005,trim_end=0.995)
job.ssea <- ssea.finish(job.ssea)

