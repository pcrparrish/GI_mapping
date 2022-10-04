
#####################################################################
## Project overview
#####################################################################

## Goals:
## - Analyze James's and Qing's pgRNA-based poison exon KO screen.

## Notes:
## - User-specified data structures:
##     screener.comparisons: tibble specifying comparisons between samples
## - Automatically created data structures:
##     screener.models: tibble holding model parameters for data normalization
##     d.counts: tibble holding raw counts in wide format
##     d.pg: tibble holding per-sgRNA fold-changes in long format
##     d.targets: tibble holding per-target fold-changes in long format
##       This is the "master" data structure that should be used for downstream
##       analyses. In addition to the sample comparisons specified in
##       screener.comparisons, it contains data for the reserved replicate name "all,"
##       which is computed by pooling data across all replicate comparisons for a
##       given origin / time.reference / time.comparison grouping.
## - Pooling data across replicates as described above is appropriate if the replicates
##   are biologically similar and therefore have similar distributions of fold-changes.
##   However, if the replicates are different, such that one replicate has a much more
##   extreme fold-change distribution (e.g., with a much greater dispersion), then
##   pooling across replicate comparisons is a poor statistical choice, since the
##   resulting distribution will be bi- or multi-modal and therefore recalcitrant
##   to simple statistical analyses. In that scenario, it is best to either remove
##   the problematic replicate(s) or analyze the replicates separately rather than
##   pooling them.

#####################################################################
## Required directories
#####################################################################

dir.pigfarm = file.path (Sys.getenv ("HOME"), "projects", "pigfarm")
dir.gene_expression = file.path (Sys.getenv ("HOME"), "projects", "gene_expression")
dir.project = file.path (dir.gene_expression,
                         "R", "projects", "2018", "poison_exon_ko_screen", Sys.getenv ("USER"))
dir.figs = gsub (file.path ("R", "projects"), file.path ("figs"), dirname (dir.project))

#####################################################################
## Load image
#####################################################################

## load image if available
file.image = file.path (dir.project, ".saved.Rdata", "image.Rdata")
if (file.exists (file.image)) {

  ## don't load the image if it appears that we've already done so
  if (exists ("d.pg") && nrow (d.pg) > 0) {
    message ("Found master data structure in memory; not loading image.")
  } else {
    message (paste0 ("Loading image from file '", file.image, "'..."), appendLF = FALSE)
    load (file.image)
    message ("done.")
  }

  ## redefine required directories to use the current home directory
  ## this allows images to be used across multiple machines w/ different
  ## home directory locations
  dir.image = gsub ("/$", "", gsub (file.path ("projects", "pigfarm"), "", dir.pigfarm))
  if (Sys.getenv ("HOME") != dir.image) {
    dir.pigfarm = gsub (dir.image, Sys.getenv ("HOME"), dir.pigfarm)
    dir.gene_expression = gsub (dir.image, Sys.getenv ("HOME"), dir.gene_expression)
    dir.project = gsub (dir.image, Sys.getenv ("HOME"), dir.project)
    dir.figs = gsub (dir.image, Sys.getenv ("HOME"), dir.figs)
    file.image = gsub (dir.image, Sys.getenv ("HOME"), file.image)
  }

}

## keep track of whether I need to update the image
## (I set this variable to TRUE if a master data structure is modified)
file.image.update = FALSE

#####################################################################
## Analysis parameters
#####################################################################

## reduce the number of parallel jobs to help avoid memory problems
options (cores = 1)

## number of samples to draw when computing FDRs
screener.numsamples = 10^4

## suppress the dplyr progress bar, because it obscures my own logging
options (dplyr.show_progress = FALSE)

#####################################################################
## Datasets
#####################################################################

if (!exists ("datasets")) datasets = list()
datasets[["human"]] = c (
  "2017/muhlemann.smg6_smg7_upf1_kd",
  "2018/bradley.poison_exon_ko")
datasets = lapply (datasets, unique)

assembly = list()
assembly[["human"]] = "hg19"

#####################################################################
## Load required packages
#####################################################################

source (file.path (dir.gene_expression, "R/spliceosome/data_structures/load_data.R"))
source (file.path (dir.gene_expression, "R/spliceosome/data_structures/load_annotation.R"))

## packages that we want to install and load
packages = c ("stringr")
for (package in packages) {
  ## install package if necessary
  if (!package %in% installed.packages())
    biocLite (package, suppressUpdates = TRUE)
  suppressPackageStartupMessages (require (package, character.only = TRUE,
                                           warn.conflicts = FALSE, quietly = TRUE))
  ## detach and reattach plyr to ensure that its functions are in the current namespace
  if (package == "dplyr") {
    detach ("package:dplyr")
    suppressPackageStartupMessages (require (package, character.only = TRUE,
                                             warn.conflicts = FALSE, quietly = TRUE))
  }
}

#####################################################################
## Load data
#####################################################################

for (species in names (datasets)) {
  load.rnaSets (
    species = species, assembly = assembly[[species]],
    datasets = datasets[[species]])
}

#####################################################################
## Global variables
#####################################################################

species = "human"

#####################################################################
## Sample comparisons
#####################################################################

## Goals:
## - Specify samples to compare.

## Notes:
## - Replicate comparisons are specified individually for the reference and comparison
##   timepoints. For example, if replicates A and B are specified for the reference
##   timepoint and replicates B and C are specified for the comparison timepoint,
##   then the comparisons A vs. B and B vs. C are performed.
## - Replicates can have any names except for "all". "all" is a reserved name that is
##   used to indicate the concatenation of all replicate comparisons for a given
##   origin / time.reference / time.comparison grouping.
## - To accomodate comparisons to the plasmid pool, I am making a copy of the
##   plasmid.0d.2A.gRNA_1.bam and plasmid.0d.2A.gRNA_2.bam samples to effectively
##   rename them to PC9.00d.2A.gRNA_1.bam and PC9.00d.2A.gRNA_2.bam. See hack below.

screener.comparisons = tibble (
  "origin" = character(), "time.reference" = character(), "time.comparison" = character(),
  "replicates.reference" = list(), "replicates.comparison" = list())

screener.comparisons = rbind (
    screener.comparisons,
    ## HeLa cell viability screen
    tibble (
        "origin" = "HeLa", "time.reference" = "0d", "time.comparison" = c ("8d", "14d"),
        "replicates.reference" = list (c ("1A", "2A", "2B", "2C", "2D")),
        "replicates.comparison" = list (c ("1A", "2A", "2B", "2C", "2D"))),
    ## PC9 cell viability screen
    tibble (
        "origin" = "PC9", "time.reference" = "0d", "time.comparison" = c ("8d", "14d"),
        "replicates.reference" = list (c ("2A", "2B", "2C", "2D")),
        "replicates.comparison" = list (c ("2A", "2B", "2C", "2D"))),
    ## parental PC9 in vivo screen: early tumors
    tibble (
        "origin" = "PC9_parental", "time.reference" = "preinj", "time.comparison" = c ("xeno_early"),
        "replicates.reference" = list (c ("2A")), "replicates.comparison" = list (c ("2A"))),
    ## parental PC9 in vivo screen: late tumors
    tibble (
        "origin" = "PC9_parental", "time.reference" = "preinj", "time.comparison" = c ("xeno_late"),
        "replicates.reference" = list (c ("2A")), "replicates.comparison" = list (c ("2A"))),
    ## PC9-Cas9 in vivo screen: early tumors
    tibble (
        "origin" = "PC9", "time.reference" = "preinj", "time.comparison" = c ("xeno_early"),
        "replicates.reference" = list (c ("2A", "2A", "2A", "2A")),
        "replicates.comparison" = list (c ("2A", "2B", "2C", "2D"))),
    ## PC9-Cas9 in vivo screen: late tumors
    tibble (
        "origin" = "PC9", "time.reference" = "preinj", "time.comparison" = c ("xeno_late"),
        "replicates.reference" = list (c ("2A", "2A", "2A", "2A", "2A", "2A", "2A", "2A", "2A", "2A")),
        "replicates.comparison" = list (c ("2A", "2B", "2C", "2D", "2E", "2F", "2G", "2H", "2I", "2J"))),
    ## HACK: PC9 comparison to plasmid pool (named PC9.00d.2A.gRNA_#.bam)
    tibble (
        "origin" = "PC9", "time.reference" = "00d", "time.comparison" = c ("0d", "8d", "14d"),
        "replicates.reference" = list (c ("2A", "2A", "2A", "2A")),
        "replicates.comparison" = list (c ("2A", "2B", "2C", "2D"))))

## if FALSE, then don't re-create existing files
rebuild = TRUE

genes.coding = genes.noncoding = list()
for (species in names (rnaAnnotations)) {
  genes.coding[[species]] = genes.compute.coding (rnaAnnotations[[species]])
  genes.noncoding[[species]] = genes.compute.noncoding (rnaAnnotations[[species]])
}

############################################################
## Build annotations
############################################################

## to do: update the below annotation parsing code once the annotations are
## sanitized and the mapping is updated to use pgRNA IDs of the format seq_1|seq_2

## Goals:
## - Parse annotation file and initialize d.counts, which holds raw counts of reads
##   supporting each pgRNA.

## Notes:
## - Columns in d.counts have the following meaning:
##   target: unique ID for each target
##     This is hacked for the upstream exon.

if (!exists ("d.counts")) {

  basedir = file.path (Sys.getenv ("HOME"), "projects", "shortread_data", "screens",
                       "human", "2018/bradley.poison_exon_ko_screen")
  file = file.path (basedir, "annotations", "pgrna", "pgrna.txt")

  stopifnot (file.exists (file))
  d.counts = read_tsv (file)
  ## restrict to the minimal information needed and reorder columns
  d.counts = d.counts %>%
    select (id, seq_1, seq_2,
            gene, target_type)
  ## for non-targeting controls, set 'gene' to NA since it isn't defined
  d.counts$gene[d.counts$target_type == "NTC"] = NA

  ## create Boolean column 'control' specifying whether a pgRNA is non-targeting
  d.counts = d.counts %>% mutate ("control" = target_type == "NTC")

}

############################################################
## Parse screen data
############################################################

## Goals:
## - Compute read counts supporting each pgRNA and store in columns named
##   counts_sample in d.counts.

basedir = file.path (Sys.getenv ("HOME"), "projects", "shortread_data", "screens",
                     "human", "2018/bradley.poison_exon_ko_screen")
## filenames follow the convention plasmid.0d.1A.gRNA_1.bam
files.1 = Sys.glob (file.path (basedir, "bam", "*.gRNA_1.bam"))
files.2 = gsub ("gRNA_1", "gRNA_2", files.1)

if (!all (file.exists (files.2)))
  stop ("Failed to find properly matched BAM files for reads 1 and 2.")
for (i in 1:length (files.1)) {

  ## parse sample information from filename
  sample = gsub ("\\.gRNA_1.bam", "", basename (files.1[i]))

  ## don't re-compute counts if they're already present in d.counts
  if (paste ("counts", sample, sep = "_") %in% colnames (d.counts))
    next
  file.image.update = TRUE

  timing = Sys.time()
  messenger (paste0 ("Parsing reads for sample ", sample, " (", i, " / ", length (files.1), ")...\n"))

  ## read in BAM output files
  param = ScanBamParam (
    ## restrict to mapped reads
    flag = scanBamFlag (isUnmappedQuery = FALSE),
    ## only read in the necessary fields
    what = c ("qname", "rname"))
  bam.1 = scanBam (files.1[i], param = param)
  bam.2 = scanBam (files.2[i], param = param)

  ## convert to tibbles
  ## we could use an all-purpose conversion method as follows:
  ##   d.bam.1 = do.call ("data.frame", bam.1) %>% as_tibble()
  ##   d.bam.2 = do.call ("data.frame", bam.2) %>% as_tibble()
  ## however, it's much faster to manually construct the tibbles
  d.bam.1 = tibble ("qname" = bam.1[[1]]$qname,
                    "rname" = factor2character (bam.1[[1]]$rname))
  d.bam.2 = tibble ("qname" = bam.2[[1]]$qname,
                    "rname" = factor2character (bam.2[[1]]$rname))
  rm (bam.1, bam.2)

  ## perform inner join between alignments of read 1 and 2 in order to obtain
  ## all possible pairings implied by each read alignment
  d.bam = inner_join (d.bam.1 %>% select (qname, "rname.1" = rname),
                      d.bam.2 %>% select (qname, "rname.2" = rname),
                      by = "qname")
  rm (d.bam.1, d.bam.2)
  ## check whether each pairing is correct
  d.bam = d.bam %>% mutate ("paired" = rname.1 == rname.2)
  
  ## if a given set of reads have one or more correct pairings, then keep the
  ## correct pairings and discard all incorrect pairings for those reads.
  ## The simplest way to achieve this is by grouping over reads as follows:
  ##   d.bam = d.bam %>% group_by (qname) %>%
  ##     mutate ("any_paired" = any (paired)) %>% ungroup() %>%
  ##     filter (paired | !any_paired) %>%
  ##     select (-any_paired)
  ## However, grouping over millions of reads is very slow. Therefore, do it manually.
  qname2anypaired = sapply (split (d.bam$paired, f = d.bam$qname), any)
  d.bam = d.bam %>% left_join (tibble ("qname" = names (qname2anypaired),
                                       "any_paired" = qname2anypaired),
                               by = "qname") %>%
    filter (paired | !any_paired) %>%
    select (-any_paired)

  ## compute weights for each set of reads
  ## As above, the straightforward method is slow:
  ##   d.bam = d.bam %>% group_by (qname) %>%
  ##     mutate ("weight" = 1 / n()) %>% ungroup()
  ## Therefore, do it manually for speed.
  qname2n = table (d.bam$qname)
  d.bam = d.bam %>% left_join (tibble ("qname" = names (qname2n),
                                       "n" = as.integer (qname2n)),
                               by = "qname") %>%
    mutate ("weight" = 1 / n) %>%
    select (-n)

  ## display statistics on the numbers of correctly paired reads
  n.paired = d.bam %>% filter (paired) %>% summarize (sum (weight)) %>% collect() %>% .[[1]]
  n = d.bam %>% summarize (sum (weight)) %>% collect() %>% .[[1]]
  messenger (paste0 ("   ", prettynum (n.paired), " / ", prettynum (n), " (", signif (100 * n.paired / n, 3), "%) correctly paired reads\n"))

  ## compute counts for correctly paired reads
  d.bam = d.bam %>% filter (paired)
  stopifnot (all (d.bam$rname.1 == d.bam$rname.2))
  d.bam = d.bam %>%
    select ("id" = rname.1, weight) %>%
    group_by (id) %>% summarize ("counts" = sum (weight)) %>% ungroup()

  ## store counts
  d.counts = d.counts %>% left_join (d.bam, by = "id")
  rm (d.bam)

  ## some pgRNAs may have no counts; replace those NAs with 0
  d.counts$counts[is.na (d.counts$counts)] = 0

  ## rename counts as counts_sample
  ## the awful rename_ invocation is necesary because I can't, for the life of me,
  ## figure out how to use rlang's crazy syntax to rename by a variable
  d.counts = d.counts %>% rename_ (.dots = setNames ("counts", paste ("counts", sample, sep = "_")))

  messenger (paste0 ("done (", signif (as.numeric (difftime (Sys.time(), timing, units = "mins")), 2), "m elapsed).\n"))

}

############################################################
## HACK: update annotations
############################################################

## to do: update the below annotation parsing code once the annotations are
## sanitized and the mapping is updated to use pgRNA IDs of the format seq_1|seq_2

file = file.path (dir.project, "data", "pgRNA.pe_lib.ann.csv")
rob = read_csv (file) %>%
  select (-pgRNA)
rob = rob %>%
  mutate ("seq_1" = sapply (strsplit (rob$sgRNAs, ","), "[[", 1),
          "seq_2" = sapply (strsplit (rob$sgRNAs, ","), "[[", 2),
          "id" = paste (seq_1, seq_2, sep = "|"),
          "type" = gsub ("\\.", "-", pairInfo),
          "event" = gsub ("@inc1$", "", inc1_ID)) %>%
  select (id, seq_1, seq_2, id, type, event)
## a given pgRNA might target multiple distinct events; simply remove duplicate entries
## for simplicity
rob = rob[!duplicated (rob$id), ]

d.counts = d.counts %>%
  select (-id) %>%
  left_join (rob,
             by = c ("seq_1", "seq_2")) %>%
  mutate ("target" = event)
d.counts$target[d.counts$target_type %in% c ("conservedPE", "unconservedPE")] = paste (d.counts$target[d.counts$target_type %in% c ("conservedPE", "unconservedPE")], "inc", sep = "@")
d.counts$target[d.counts$target_type == "upstreamExon"] = paste (d.counts$target[d.counts$target_type == "upstreamExon"], "upstreamExon", sep = "@")

d.counts = d.counts %>% select (id, control, type,
                                event, target, target_type,
                                starts_with ("counts"))

#####################################################################
## Fold-changes: per-pgRNA
#####################################################################

## Goals:
## - Compute raw fold-changes for each pgRNA and store in d.pg.

## Notes:
## - The parameter pseudocountscale is used to compute a per-pgRNA pseudocount
##   for the reference and comparison samples, defined as:
##   reference: pseudocountscale * (counts for pgRNA in reference)
##   comparison: pseudocountscale * (counts for pgRNA in reference) /
##     (counts for all pgRNAs in comparison / counts for all pgRNAs in reference)
##   The second scaling factor for the comparison sample normalizes to total library size.
##   In principle, this pseudocount takes account of the fact that pgRNA representation
##   is not even in our library, so we should use a stronger pseudocount for pgRNAs
##   that are overrepresented in the library and a weaker pseudocount for pgRNAs
##   that are underrepresented in the library. A pseudocount >= minpseudocount
##   is always assigned, even for pgRNAs with 0 counts in the reference sample.

foo = function (d, d.pg,
                pseudocountscale=0.05, minpseudocount=5) {

  ## check that data has been properly grouped
  stopifnot (nrow (d) == 1)

  origin = d$origin
  time.reference = d$time.reference
  time.comparison = d$time.comparison
  replicates.reference = unlist (d$replicates.reference)
  replicates.comparison = unlist (d$replicates.comparison)
  stopifnot (length (replicates.reference) == length (replicates.comparison))

  d.pg.grouping = d.pg %>% filter (origin == !!origin & time.reference == !!time.reference & time.comparison == !!time.comparison)
  
  for (i in 1:length (replicates.reference)) {

    replicate.reference = replicates.reference[i]
    replicate.comparison = replicates.comparison[i]

    ## don't recompute
    if (nrow (d.pg.grouping %>% filter (replicate.reference == !!replicate.reference & replicate.comparison == !!replicate.comparison)))
      next

    ## create a temporary data frame with static variable names
    d.pg.i = d.counts %>%
      mutate_ ("counts.reference" = paste ("counts", collapse (origin, time.reference, replicate.reference, sep = "."), sep = "_"),
               "counts.comparison" = paste ("counts", collapse (origin, time.comparison, replicate.comparison, sep = "."), sep = "_")) %>%
      select (id, control, type,
              event, target, target_type,
              counts.reference, counts.comparison)

    ## compute per-pgRNA pseudocounts
    pseudocounts.reference = pseudocountscale * d.pg.i$counts.reference
    pseudocounts.reference[pseudocounts.reference < minpseudocount] = minpseudocount
    pseudocounts.comparison = pseudocounts.reference * (sum (d.pg.i$counts.comparison) / sum (d.pg.i$counts.reference))
    pseudocounts.comparison[pseudocounts.comparison < minpseudocount] = minpseudocount
    d.pg.i = d.pg.i %>% mutate ("pseudocounts.reference" = pseudocounts.reference,
                                "pseudocounts.comparison" = pseudocounts.comparison)

    ## compute fold-change, remembering to normalize to library size
    ## we compute the library size explicitly, rather than within dplyr, because
    ## we need to sum over rows rather than perform a vectorized computation
    n.reference = sum (d.pg.i$counts.reference) + sum (d.pg.i$pseudocounts.reference)
    n.comparison = sum (d.pg.i$counts.comparison) + sum (d.pg.i$pseudocounts.comparison)
    d.pg.i = d.pg.i %>%
      mutate ("fc" = ((counts.comparison + pseudocounts.comparison) / n.comparison) / ((counts.reference + pseudocounts.reference) / n.reference))
    
    ## add information specifying the current comparison
    d.pg.i = d.pg.i %>%
      mutate ("origin" = origin, "time.reference" = time.reference, "time.comparison" = time.comparison, "replicate.reference" = replicate.reference, "replicate.comparison" = replicate.comparison)

    ## plot counts before and after pseudocount addition

    ## create a matrix of panels
    p = list()

    ######### panel 1: scatter plot of raw counts for non-targeting pgRNAs #########
    panel = 1 ## panel index
    p[[panel]] <- ggplot (mapping = aes (x = counts.reference, y = counts.comparison),
                          data = d.pg.i %>% filter (control))
    p[[panel]] <- p[[panel]] + geom_abline (intercept = 0, slope = 1,
                                            color = "black", alpha = 0.25,
                                            linetype = "solid", size = 0.5)
    p[[panel]] <- p[[panel]] + geom_point (shape = 16, alpha = 0.25,
                                           size = 0.5, color = "black")
    xlim = ylim = range (pretty (c (0, max (75, c (quantile (d.pg.i$counts.reference, 0.95), quantile (d.pg.i$counts.comparison, 0.95))))))
    p[[panel]] <- p[[panel]] + coord_cartesian (xlim = xlim, ylim = ylim)
    p[[panel]] <- p[[panel]] + xlab (bquote ("counts /" ~ 10^6 ~ "(" * .(time.reference) * ")"))
    p[[panel]] <- p[[panel]] + ylab (bquote ("counts /" ~ 10^6 ~ "(" * .(time.comparison) * ")"))
    p[[panel]] <- p[[panel]] + ggtitle ("non-targeting pgRNAs")

    ######### panel 2: scatter plot of raw counts for targeting pgRNAs #########
    panel = 2 ## panel index
    p[[panel]] <- ggplot (mapping = aes (x = counts.reference, y = counts.comparison),
                          data = d.pg.i %>% filter (!control))
    p[[panel]] <- p[[panel]] + geom_abline (intercept = 0, slope = 1,
                                            color = "black", alpha = 0.25,
                                            linetype = "solid", size = 0.5)
    p[[panel]] <- p[[panel]] + geom_point (shape = 16, alpha = 0.25,
                                           size = 0.5, color = "black")
    xlim = ylim = range (pretty (c (0, max (75, c (quantile (d.pg.i$counts.reference, 0.95), quantile (d.pg.i$counts.comparison, 0.95))))))
    p[[panel]] <- p[[panel]] + coord_cartesian (xlim = xlim, ylim = ylim)
    p[[panel]] <- p[[panel]] + xlab (bquote ("counts /" ~ 10^6 ~ "(" * .(time.reference) * ")"))
    p[[panel]] <- p[[panel]] + ylab (bquote ("counts /" ~ 10^6 ~ "(" * .(time.comparison) * ")"))
    p[[panel]] <- p[[panel]] + ggtitle ("targeting pgRNAs")

    ######### panel 3: scatter plot of raw counts + pseudocounts for non-targeting pgRNAs #########
    panel = 3 ## panel index
    p[[panel]] <- ggplot (mapping = aes (x = counts_with_pseudocounts.reference, y = counts_with_pseudocounts.comparison),
                          data = d.pg.i %>% filter (control) %>%
                            mutate ("counts_with_pseudocounts.reference" = counts.reference + pseudocounts.reference,
                                    "counts_with_pseudocounts.comparison" = counts.comparison + pseudocounts.comparison))
    p[[panel]] <- p[[panel]] + geom_abline (intercept = 0, slope = 1,
                                            color = "black", alpha = 0.25,
                                            linetype = "solid", size = 0.5)
    p[[panel]] <- p[[panel]] + geom_point (shape = 16, alpha = 0.25,
                                           size = 0.5, color = "black")
    xlim = ylim = range (pretty (c (0, max (75, c (quantile (d.pg.i$counts.reference, 0.95), quantile (d.pg.i$counts.comparison, 0.95))))))
    p[[panel]] <- p[[panel]] + coord_cartesian (xlim = xlim, ylim = ylim)
    p[[panel]] <- p[[panel]] + xlab (bquote ("counts + pseudocounts /" ~ 10^6 ~ "(" * .(time.reference) * ")"))
    p[[panel]] <- p[[panel]] + ylab (bquote ("counts + pseudocounts /" ~ 10^6 ~ "(" * .(time.comparison) * ")"))
    p[[panel]] <- p[[panel]] + ggtitle ("non-targeting pgRNAs")

    ######### panel 4: scatter plot of raw counts + pseudocounts for targeting pgRNAs #########
    panel = 4 ## panel index
    p[[panel]] <- ggplot (mapping = aes (x = counts_with_pseudocounts.reference, y = counts_with_pseudocounts.comparison),
                          data = d.pg.i %>% filter (!control) %>%
                            mutate ("counts_with_pseudocounts.reference" = counts.reference + pseudocounts.reference,
                                    "counts_with_pseudocounts.comparison" = counts.comparison + pseudocounts.comparison))
    p[[panel]] <- p[[panel]] + geom_abline (intercept = 0, slope = 1,
                                            color = "black", alpha = 0.25,
                                            linetype = "solid", size = 0.5)
    p[[panel]] <- p[[panel]] + geom_point (shape = 16, alpha = 0.25,
                                           size = 0.5, color = "black")
    xlim = ylim = range (pretty (c (0, max (75, c (quantile (d.pg.i$counts.reference, 0.95), quantile (d.pg.i$counts.comparison, 0.95))))))
    p[[panel]] <- p[[panel]] + coord_cartesian (xlim = xlim, ylim = ylim)
    p[[panel]] <- p[[panel]] + xlab (bquote ("counts + pseudocounts /" ~ 10^6 ~ "(" * .(time.reference) * ")"))
    p[[panel]] <- p[[panel]] + ylab (bquote ("counts + pseudocounts /" ~ 10^6 ~ "(" * .(time.comparison) * ")"))
    p[[panel]] <- p[[panel]] + ggtitle ("targeting pgRNAs")

    ## create plot
    file = file.path (dir.figs, origin, replicate.reference, replicate.comparison,
                      collapse (time.reference, time.comparison, sep = "-vs-"),
                      paste ("scatter", "pseudocount_effect", sep = "."))
    if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
      file = device.open (plotdim = c (2 * 3, 2 * 3), format = "pdf", file)
      multiplot (plotlist = p, cols = 2, byrow = TRUE)
      device.close (file)
    }
    
    ## remove temporary variables
    d.pg.i = d.pg.i %>% select (-counts.reference, -counts.comparison, -pseudocounts.reference, -pseudocounts.comparison)

    ## store
    d.pg.grouping = rbind (d.pg.grouping,
                           d.pg.i)

  }

  ## reorder columns to fit master tibble
  if (nrow (d.pg)) {
    stopifnot (all (colnames (d.pg) %in% colnames (d.pg.grouping)))
    d.pg.grouping = d.pg.grouping %>% select (colnames (d.pg))
    stopifnot (identical (colnames (d.pg), colnames (d.pg.grouping)))
  }

  return (d.pg.grouping)

}

if (!exists ("d.pg")) {
  d.pg = tibble (
    "origin" = character(), "time.reference" = character(), "time.comparison" = character(), "replicate.reference" = character(), "replicate.comparison" = character(),
    "id" = character(), "control" = logical(), "type" = character(),
    "event" = character(), "target" = character(), "target_type" = character(),
    "fc" = numeric())
}
d.pg.new = screener.comparisons %>%
  group_by (origin, time.reference, time.comparison) %>%
  do (foo (., d.pg,
            pseudocountscale = 0.05, minpseudocount = 5)) %>%
  ungroup()
if (!identical (d.pg, d.pg.new)) {
  file.image.update = TRUE
  d.pg = d.pg.new
}
rm (d.pg.new)

#####################################################################
## Target annotation: per-pgRNA
#####################################################################

## Goals:
## - Compute useful information about targets (cassette exons), including parent gene
##   expression, exon inclusion, etc.

## Notes:
## - Although the target annotation is stored in d.pg, the annotation is defined
##   at a per-target, not per-pgRNA, level.
## - Entries are set to NA for non-targeting pgRNAs.

## get junctions upstream and downstream of the cassette exon, which we need in order
## to access the splicing event annotations
## we define these over all rows of d.pg, with NAs representing non-targeting pgRNAs,
## for easier row indexing below
junctions.inc1 = junctions.inc2 = rep (NA, nrow (d.pg))
junctions.inc1[!d.pg$control] = events.event2junction (events = d.pg$event[!d.pg$control], isoforms = 1,
                                                       identifyingJunctions = 1)
junctions.inc2[!d.pg$control] = events.event2junction (events = d.pg$event[!d.pg$control], isoforms = 1,
                                                       identifyingJunctions = 2)

## gene and gene name
if (!all (c ("gene", "geneName") %in% colnames (d.pg))) {
  genes = geneNames = rep (NA, nrow (d.pg))
  ## the as.vector calls strips away the names to prevent future warnings when
  ## joining with a tibble lacking names for these columns
  genes[!d.pg$control] = as.vector (get.splicing.annotation (rnaAnnotations[[species]], "gene")[junctions.inc1[!d.pg$control]])
  geneNames[!d.pg$control] = as.vector (genes.id2name (rnaAnnotations[[species]], genes[!d.pg$control]))
  d.pg = d.pg %>%
    mutate ("gene" = genes,
            "geneName" = geneNames)
  file.image.update = TRUE
}

## sequence conservation
## conservation is defined based on:
## - the upstream 3' and downstream 5' splice sites of the cassette exon for targets
##   of type unconservedPE and conservedPE
## - the downstream 5' splice site of the upstream exon for targets of type upstreamExon
if (! "conservation" %in% colnames (d.pg)) {
  conservation = matrix (
    NA,
    nrow = nrow (d.pg), ncol = 300)
  ## - compute for the cassette exon
  conservation[!d.pg$control, ] = cbind (
    ## upstream 3' splice site: [-100, +50]
    get.splicing.annotation (rnaAnnotations[[species]], "3ss_conservation")[junctions.inc1[!d.pg$control], 1:150],
    ## downstream 5' splice site: [-50, +100]
    get.splicing.annotation (rnaAnnotations[[species]], "5ss_conservation")[junctions.inc2[!d.pg$control], 51:200])
  conservation = rowMeans (conservation)
  ## - compute for the upstream exon
  ##   note that the below indexing implicitly requires !d.pg$control to be TRUE,
  ##   so we don't need to enforce that when accessing junctions.inc1
  indices = d.pg$target_type == "upstreamExon"
  ##   We really should compute this over the 5' end of the upstream exon as well.
  ##   However, it's a pain to perform the necessary coordinate matching, so we'll
  ##   just go with the 3' end.
  conservation[indices] = rowMeans (
    ## upstream 5' splice site: [-50, +100]
    get.splicing.annotation (rnaAnnotations[[species]], "5ss_conservation")[junctions.inc1[indices], 51:200])
  ## - store
  d.pg = d.pg %>%
    mutate ("conservation" = conservation)
  file.image.update = TRUE
}

## parent gene expression
if (! "expr" %in% colnames (d.pg)) {
  dataset = "2018/bradley.poison_exon_ko"
  ## compute normalized gene expression
  f = tmm.compute.normalization (
    rnaSets[[species]][[dataset]],
    genes = genes.compute.coding (rnaAnnotations[[species]]))
  expr = .genes.compute.expression (
    rnaSets[[species]][[dataset]], "geneExpr",
    f = f)
  ## use the first RNA-seq dataset to compute gene expression for HeLa and PC9 cells,
  ## since it was performed closer in time to the screen
  expr = cbind (expr,
                "HeLa" = expr[, "HeLa.0d.2A"],
                "PC9" = expr[, "PC9.0d.2A"])
  ## create a column corresponding to PC9_parental (fake the gene expression data for
  ## those cells by assuming that it's identical to that for PC9 cells)
  expr = cbind (expr, "PC9_parental" = expr[, "PC9", drop = TRUE])
  ## confirm that we have RNA-seq samples for each origin
  stopifnot (all (unique (screener.comparisons$origin) %in% colnames (expr)))
  ## convert to long format for subsequent join with d.pg
  d.expr = rownames_to_column (expr %>% as.data.frame, var = "gene") %>%
    gather (key = "origin", value = "expr", ... = -gene) %>% as_tibble()
  ## warn if genes targeted by pgRNAs aren't present
  indices = d.pg$control | d.pg$gene %in% d.expr$gene
  if (!all (indices))
    warning (paste0 ("Failed to find gene expression estimates for ", prettynum (sum (!indices)), " genes."))
  ## store
  d.pg = d.pg %>% left_join (d.expr, by = c ("origin", "gene"))
  file.image.update = TRUE
}

## cassette exon inclusion
## - compute psi values for both basal and NMD-inhibited states
if (! all (c ("psi.basal", "psi.nmd_inhibited",
              "psi_muhlemann.basal", "psi_muhlemann.nmd_inhibited") %in% colnames (d.pg))) {

  dataset = "2018/bradley.poison_exon_ko"

  ## compute psi for basal (unperturbed) state
  psi.basal = get.expr (get.spliceEvent (rnaSets[[species]][[dataset]]))
  ## average across replicates for HeLa and PC9 cells
  psi.basal = cbind (psi.basal,
                     "HeLa" = rowMeans (psi.basal[, paste ("HeLa.ctrl.0d", c ("A", "B", "C"), sep = ".")]),
                     "PC9" = rowMeans (psi.basal[, paste ("PC9.ctrl.0d", c ("A", "B", "C"), sep = ".")]))
  ## create a column corresponding to PC9_parental (fake the splicing data for
  ## those cells by assuming that it's identical to that for PC9 cells)
  psi.basal = cbind (psi.basal, "PC9_parental" = psi.basal[, "PC9", drop = TRUE])
  ## confirm that we have RNA-seq samples for each origin
  stopifnot (all (unique (screener.comparisons$origin) %in% colnames (psi.basal)))
  ## restrict to those samples
  psi.basal = psi.basal[, unique (screener.comparisons$origin)]
  ## convert to long format for subsequent join with d.pg
  d.psi.basal = rownames_to_column (psi.basal %>% as.data.frame, var = "event.isoform") %>%
    gather (key = "origin", value = "psi.basal", ... = -event.isoform) %>% as_tibble()
  ## warn if event.isoform IDs not present
  d.pg = d.pg %>%
    mutate ("event.isoform" = paste (d.pg$event, "1", sep = "."))
  indices = d.pg$control | d.pg$event.isoform %in% d.psi.basal$event.isoform
  if (!all (indices))
    warning (paste0 ("Failed to find isoform estimates for ", prettynum (sum (!indices)), " events."))
  ## store
  d.pg = d.pg %>%
    left_join (d.psi.basal, by = c ("origin", "event.isoform")) %>%
    select (-event.isoform)

  ## compute psi for NMD-inhibited state
  psi.nmd_inhibited = get.expr (get.spliceEvent (rnaSets[[species]][[dataset]]))
  ## average across replicates for HeLa and PC9 cells
  psi.nmd_inhibited = cbind (psi.nmd_inhibited,
                     "HeLa" = rowMeans (psi.nmd_inhibited[, paste ("HeLa.chx.0d", c ("A", "B", "C"), sep = ".")]),
                     "PC9" = rowMeans (psi.nmd_inhibited[, paste ("PC9.chx.0d", c ("A", "B", "C"), sep = ".")]))
  ## create a column corresponding to PC9_parental (fake the splicing data for
  ## those cells by assuming that it's identical to that for PC9 cells)
  psi.nmd_inhibited = cbind (psi.nmd_inhibited, "PC9_parental" = psi.nmd_inhibited[, "PC9", drop = TRUE])
  ## confirm that we have RNA-seq samples for each origin
  stopifnot (all (unique (screener.comparisons$origin) %in% colnames (psi.nmd_inhibited)))
  ## restrict to those samples
  psi.nmd_inhibited = psi.nmd_inhibited[, unique (screener.comparisons$origin)]
  ## convert to long format for subsequent join with d.pg
  d.psi.nmd_inhibited = rownames_to_column (psi.nmd_inhibited %>% as.data.frame, var = "event.isoform") %>%
    gather (key = "origin", value = "psi.nmd_inhibited", ... = -event.isoform) %>% as_tibble()
  ## warn if event.isoform IDs not present
  d.pg = d.pg %>%
    mutate ("event.isoform" = paste (d.pg$event, "1", sep = "."))
  indices = d.pg$control | d.pg$event.isoform %in% d.psi.nmd_inhibited$event.isoform
  if (!all (indices))
    warning (paste0 ("Failed to find isoform estimates for ", prettynum (sum (!indices)), " events."))
  ## store
  d.pg = d.pg %>%
    left_join (d.psi.nmd_inhibited, by = c ("origin", "event.isoform")) %>%
    select (-event.isoform)

  ## now compute psi values for the basal and NMD-inhibited states using the Muhlemann dataset
  dataset = "2017/muhlemann.smg6_smg7_upf1_kd"

  ## compute psi for basal (unperturbed) state
  psi.basal = get.expr (get.spliceEvent (rnaSets[[species]][[dataset]]))
  ## average across replicates for HeLa and PC9 cells
  psi.basal = cbind (psi.basal,
                     "HeLa" = psi.basal[, "HeLa.control_kd.B"],
                     "PC9" = NA)
  ## create a column corresponding to PC9_parental (fake the splicing data for
  ## those cells by assuming that it's identical to that for PC9 cells)
  psi.basal = cbind (psi.basal, "PC9_parental" = psi.basal[, "PC9", drop = TRUE])
  ## confirm that we have RNA-seq samples for each origin
  stopifnot (all (unique (screener.comparisons$origin) %in% colnames (psi.basal)))
  ## restrict to those samples
  psi.basal = psi.basal[, unique (screener.comparisons$origin)]
  ## convert to long format for subsequent join with d.pg
  d.psi.basal = rownames_to_column (psi.basal %>% as.data.frame, var = "event.isoform") %>%
    gather (key = "origin", value = "psi_muhlemann.basal", ... = -event.isoform) %>% as_tibble()
  ## warn if event.isoform IDs not present
  d.pg = d.pg %>%
    mutate ("event.isoform" = paste (d.pg$event, "1", sep = "."))
  indices = d.pg$control | d.pg$event.isoform %in% d.psi.basal$event.isoform
  if (!all (indices))
    warning (paste0 ("Failed to find isoform estimates for ", prettynum (sum (!indices)), " events."))
  ## store
  d.pg = d.pg %>%
    left_join (d.psi.basal, by = c ("origin", "event.isoform")) %>%
    select (-event.isoform)

  ## compute psi for NMD-inhibited state
  ## we only have data for HeLa cells treated with shRNAs against SMG6 and SMG7
  psi.nmd_inhibited = psi.basal
  psi.nmd_inhibited[, "HeLa"] = get.expr (get.spliceEvent (rnaSets[[species]][[dataset]]))[, "HeLa.smg6_kd-smg7_kd.B"]
  ## for other cell types, set all psi values to NA since we don't have NMD inhibition
  ## data
  psi.nmd_inhibited[, setdiff (colnames (psi.nmd_inhibited), "HeLa")] = NA
  ## convert to long format for subsequent join with d.pg
  d.psi.nmd_inhibited = rownames_to_column (psi.nmd_inhibited %>% as.data.frame, var = "event.isoform") %>%
    gather (key = "origin", value = "psi_muhlemann.nmd_inhibited", ... = -event.isoform) %>% as_tibble()
  ## warn if event.isoform IDs not present
  d.pg = d.pg %>%
    mutate ("event.isoform" = paste (d.pg$event, "1", sep = "."))
  indices = d.pg$control | d.pg$event.isoform %in% d.psi.nmd_inhibited$event.isoform
  if (!all (indices))
    warning (paste0 ("Failed to find isoform estimates for ", prettynum (sum (!indices)), " events."))
  ## store
  d.pg = d.pg %>%
    left_join (d.psi.nmd_inhibited, by = c ("origin", "event.isoform")) %>%
    select (-event.isoform)

  file.image.update = TRUE

}

## reorder columns
d.pg = d.pg %>%
  select (origin, time.reference, time.comparison, replicate.reference, replicate.comparison,
          id, control, type,
          event, target, target_type,
          gene, geneName, conservation, expr, starts_with ("psi"),
          fc)

#####################################################################
## Normalization to unexpressed genes: per-pgRNA
#####################################################################

## Goals:
## - Normalize fold-changes s.t. the median fold-change for pgRNAs targeting
##   unexpressed genes is 1.

## Notes:
## - The normalized fold-changes in fc_norm are computed by dividing unnormalized
##   fold-changes by the median for all pgRNAs targeting unexpressed genes.

## threshold for defining unexpressed genes
maxExpression = 1

foo = function (d, d.pg) {

  ## check that data has been properly grouped
  stopifnot (nrow (d) == nrow (d %>% select (id) %>% distinct()))

  ## compute and store normalization
  normalization = d %>% filter (!control & expr <= maxExpression) %>% select (fc) %>% collect() %>% .[[1]] %>% median()
  d = d %>% mutate ("normalization" = normalization)

  ## normalize fold-changes
  d = d %>% mutate ("fc_norm" = fc / normalization) %>% select (-normalization)

  return (d)

}

d.pg.new = d.pg %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison) %>%
  do (foo (., d.pg)) %>%
  ungroup()
if (!identical (d.pg, d.pg.new)) {
  file.image.update = TRUE
  d.pg = d.pg.new
}
rm (d.pg.new)

#####################################################################
## p-values and FDRs: per-target
#####################################################################

## Goals:
## - Compute a p-value for a difference in distribution of fold-changes for all
##   pgRNAs designed to hit a given target relative to all pgRNAs targeting
##   unexpressed genes.
## - Compute per-target FDRs by converting per-target p-values into empirical FDRs
##   by sampling from the fold-changes of pgRNAs targeting unexpressed genes
##   to generate fake targets and then computing associated p-values.

## Notes:
## - p-values and FDRs are computed using normalized fold-changes (fc_norm).
## - Use a two-sided Mann-Whitney test to compute p-values.
## - If the reserved replicate name "all" is specified for both replicate.reference
##   and replicate.comparison, then duplicate entries are allowed for a given pgRNA
##   ID. The code assumes that replicate "all" corresponds to pooling data across
##   replicate comparisons, and so it's appropriate to have multiple entries for a
##   given pgRNA ID.

## threshold for defining unexpressed genes
maxExpression = 1

## numsamples: number of samples to draw when computing FDRs
##   Intuitively, I think that this should be >= the number of targets.
foo = function (d, d.targets, numsamples) {

  origin = d$origin[1]
  time.reference = d$time.reference[1]
  time.comparison = d$time.comparison[1]
  replicate.reference = d$replicate.reference[1]
  replicate.comparison = d$replicate.comparison[1]

  ## handle the case of the reserved replicate name "all"
  if (xor (replicate.reference == "all", replicate.comparison == "all"))
    stop ("Either both or neither of 'replicate.reference' and 'replicate.comparison' must be the reserved replicate name 'all'.")

  ## check that data has been properly grouped
  if (replicate.reference != "all")
    stopifnot (nrow (d) == nrow (d %>% select (id) %>% distinct()))

  ## don't recompute
  d.targets.grouping = d.targets %>%
    filter (origin == !!origin & time.reference == !!time.reference & time.comparison == !!time.comparison &
            replicate.reference == !!replicate.reference & replicate.comparison == !!replicate.comparison)
  if (nrow (d.targets.grouping))
    return (d.targets.grouping)

  timing = Sys.time()
  messenger (paste0 ("Computing p-values for comparison ",
                     collapse (c (origin,
                                  collapse (time.reference, time.comparison, sep = " vs. "),
                                  collapse (replicate.reference, replicate.comparison, sep = " vs. ")),
                               sep = " / "), "..."))
  
  ######### compute p-values #########

  ## vector of normalized fold-changes for pgRNAs targeting unexpressed genes
  fc_norm.control = d %>% filter (!control & expr <= maxExpression) %>% select (fc_norm) %>% collect() %>% .[[1]]
  ## list mapping targets to vectors of unnormalized and normalized fold-changes for
  ## targeting pgRNAs
  target2fc = split (d %>% filter (!control) %>% select (fc) %>% collect() %>% .[[1]],
                     f = d %>% filter (!control) %>% select (target) %>% collect() %>% .[[1]])
  target2fc_norm = split (d %>% filter (!control) %>% select (fc_norm) %>% collect() %>% .[[1]],
                          f = d %>% filter (!control) %>% select (target) %>% collect() %>% .[[1]])

  ## compute p-values for a difference in distribution in normalized fold-changes
  target2pval = mapply (function (x, y) wilcox.test (x, y, alternative = "two.sided")$p.value,
                        listify (fc_norm.control, length (target2fc_norm)),
                        target2fc_norm)
  names (target2pval) = names (target2fc_norm)

  ## store mean fold-change, unnormalized and normalized fold-changes for all
  ## targeting pgRNAs, and p-value
  d.targets.grouping = tibble (
    "origin" = origin, "time.reference" = time.reference, "time.comparison" = time.comparison, "replicate.reference" = replicate.reference, "replicate.comparison" = replicate.comparison,
    "target" = names (target2fc_norm),
    "fc_mean" = sapply (target2fc, mean.geometric), "fc_norm_mean" = sapply (target2fc_norm, mean.geometric),
    "fc" = target2fc, "fc_norm" = target2fc_norm, 
    "pval" = target2pval) %>%
    left_join (d %>% select (event, target, target_type) %>% distinct(),
               by = "target")

  ######### compute FDRs #########

  ## compute the median number of pgRNAs associated with each target
  numguides = median (d %>% filter (!control) %>% group_by (target) %>% summarize ("n" = n()) %>% ungroup() %>% select (n) %>% collect() %>% .[[1]])
  messenger (paste0 ("estimating FDRs with ", prettynum (numsamples), " groups of ", numguides, " pgRNAs targeting unexpressed genes..."))
  ## list of vectors of fold-changes for groups of numguides randomly selected pgRNAs
  ## targeting unexpressed genes (to mimic the numguides targeting pgRNAs present for most targets)
  ## Iteratively re-sampling is straightfoward:
  ##   faketarget2fc_norm = lapply (1:numsamples, function (x) sample (fc_norm.control, numguides))
  ## but approximately 100X slower than sampling once to get all values and then
  ## splitting as necessary. The straightforward but slow way to get a vector for
  ## splitting is this:
  ##   sort (rep (1:(numguides * numsamples), numguides))
  ## The non-intuitive matrix cast below is much faster.
  faketarget2fc_norm = split (sample (fc_norm.control, numsamples * numguides, replace = TRUE),
                              f = as.vector (matrix (rep (1:numsamples, numguides), nrow = 18, byrow = T)))

  ## compute p-values for a difference in distribution
  faketarget2pval = mapply (function (y) wilcox.test (x = fc_norm.control, y = y,
                                                      alternative = "two.sided")$p.value,
                            faketarget2fc_norm)

  ## use the distribution of p-values to estimate FDRS for each gene
  fdr = rep (NA, nrow (d.targets.grouping))
  for (i in 1:length (fdr)) {
    fdr[i] = sum (faketarget2pval <= d.targets.grouping$pval[i]) / length (faketarget2pval)
  }

  ## store FDRs for all targets
  d.targets.grouping = d.targets.grouping %>% mutate ("fdr" = fdr)

  messenger (paste0 ("done (", signif (as.numeric (difftime (Sys.time(), timing, units = "mins")), 2), "m elapsed).\n"))

  return (d.targets.grouping)

}

if (!exists ("d.targets")) {
  d.targets = tibble (
    "origin" = character(), "time.reference" = character(), "time.comparison" = character(), "replicate.reference" = character(), "replicate.comparison" = character(),
    "event" = character(), "target" = character(), "target_type" = character(),
    "fc_mean" = numeric(), "fc_norm_mean" = numeric(),
    "fc" = list(), "fc_norm" = list(),
    "pval" = numeric(), "fdr" = numeric())
}
## compute for all replicate comparisons
d.targets.new = d.pg %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison) %>%
  do (foo (., d.targets, numsamples = screener.numsamples)) %>%
  ungroup()
## compute for the reserved replicate comparison "all" vs. "all" by pooling across available replicates
d.targets.new = rbind (
  d.targets.new,
  d.pg %>%
  group_by (origin, time.reference, time.comparison) %>%
  mutate ("replicate.reference" = "all", "replicate.comparison" = "all") %>%
  do (foo (., d.targets, numsamples = screener.numsamples)) %>%
  ungroup())
## store
if (!identical (d.targets, d.targets.new)) {
  file.image.update = TRUE
  d.targets = d.targets.new
}
rm (d.targets.new)

## store target annotation
d.targets = d.targets %>%
  left_join (d.pg %>% filter (!control) %>%
             select (origin, event, target_type,
                     gene, geneName, conservation, expr, starts_with ("psi")) %>% distinct(),
             by = c ("origin", "event", "target_type"))

## reorder columns
d.targets = d.targets %>%
  select (origin, time.reference, time.comparison, replicate.reference, replicate.comparison,
          event, target, target_type,
          gene, geneName, conservation, expr, starts_with ("psi"),
          fc_mean, fc_norm_mean, fc, fc_norm, pval, fdr)

#####################################################################
## Save image
#####################################################################

if (file.image.update) {
  messenger (paste0 ("Saving image to file '", file.image, "'..."))
  dir.create (dirname (file.image), showWarnings = FALSE, recursive = TRUE)
  ## save everything except for core data structures, which are bulky
  ## and frequently need to be reloaded anyways due to database peculiarities
  save (list = setdiff (ls (all.names = TRUE), c ("rnaSet", "rnaSets", "rnaAnnotations", "genomeAnnotations")),
        file = file.image,
        envir = .GlobalEnv)
  file.image.update = FALSE
  messenger ("done.\n")
}

#####################################################################
## Table of fold-changes for all pgRNAs and targets
#####################################################################

## Goals:
## - Create wide-format tables containing:
##   * fold-changes for all pgRNAs
##   * mean fold-changes, FDRs, and p-values for all targets

## Notes:
## - These tables are intended to be supplementary information for the paper.

## Rename column fc to fc|origin|time.reference|time.comparison|replicate.reference|replicate.comparison.
##
## d: tibble restricted to a single origin/time.reference/time.comparison/replicate.reference/replicate.comparison
foo = function (d) {

  ## if d is derived from d.pg
  if ("id" %in% colnames (d)) {
    stopifnot (nrow (d) == nrow (d %>% select (id) %>% distinct()))
  }
  ## if d is derived from d.targets
  else {
    stopifnot (nrow (d) == nrow (d %>% select (event, target_type) %>% distinct()))
  }
  
  ## rename columns
  for (colname in c ("fc", "fc_norm", "fc_mean", "fc_norm_mean", "pval", "fdr")) {
    if (! colname %in% colnames (d))
      next
    d = d %>%
      rename (!!paste (colname, d$origin[1], d$time.reference[1], d$time.comparison[1], d$replicate.reference[1], d$replicate.comparison[1], sep = "|") := !!colname)
  }

  ## strip away now-unnecessary columns
  d = d %>% select (-origin, -time.reference, -time.comparison, -replicate.reference, -replicate.comparison)

  return (d)

}

######### per-pgRNA #########

## convert from long to wide format:
## - split long tibble into a list of sub-tibbles and rename columns for each
##   sub-tibble based on the grouping
d = lapply (split (d.pg %>% select (-target, -gene, -geneName, -conservation, -expr, -starts_with ("psi")),
                   f = paste0 (d.pg$origin, d.pg$time.reference, d.pg$time.comparison, d.pg$replicate.reference, d.pg$replicate.comparison)),
            foo)
## - join the sub-tibbles to get a single wide tibble
d = d %>% Reduce (function (d1, d2) { left_join (d1, d2,
                                                 by = c ("id", "control", "type", "event", "target_type")) }, .)

## round numeric columns
d = d %>% mutate_if (is.numeric, signif, 5)

## write
file = file.path (dir.figs,
                  collapse ("fc", "per-pgRNA", "table", sep = "."))
dir.create (dirname (file), showWarnings = FALSE, recursive = TRUE)
suppressWarnings (write.table (d, file = file,
                               quote = FALSE, sep = "\t",
                               row.names = FALSE, col.names = TRUE))

######### per-target #########

## convert from long to wide format:
## - split long tibble into a list of sub-tibbles and rename columns for each
##   sub-tibble based on the grouping
d = lapply (split (d.targets %>% select (-target, -conservation, -expr, -starts_with ("psi"), -fc, -fc_norm),
                   f = paste0 (d.targets$origin, d.targets$time.reference, d.targets$time.comparison, d.targets$replicate.reference, d.targets$replicate.comparison)),
            foo)
## - join the sub-tibbles to get a single wide tibble
d = d %>% Reduce (function (d1, d2) { left_join (d1, d2,
                                                 by = c ("event", "target_type", "gene", "geneName")) }, .)

## round numeric columns
d = d %>% mutate_if (is.numeric, signif, 5)

## write
file = file.path (dir.figs,
                  collapse ("fc", "per-target", "table", sep = "."))
dir.create (dirname (file), showWarnings = FALSE, recursive = TRUE)
suppressWarnings (write.table (d, file = file,
                               quote = FALSE, sep = "\t",
                               row.names = FALSE, col.names = TRUE))

#####################################################################
## Basic global analyses: per-pgRNA
#####################################################################

## Goals:
## - Create histograms comparing fold-changes for different sets of pgRNAs or targets:
##     plot 1: control vs. targeting pgRNAs (plot of per-pgRNA fold-changes)
##     plot 2: unexpressed vs. expressed genes (plot of per-target fold-changes,
##       estimated by taking the mean over all pgRNAs targeting each target)
##     plots 3-4: excluded vs. included cassette exons (plots of per-target fold-changes,
##       estimated by taking the mean over all pgRNAs targeting each target),
##       with psi values computed for basal or NMD-inhibited states.
##       Plots restricted to cassette exons within expressed genes.
##   Each of the above plots is (manually) faceted by target_type.

## Conclusions:
## - 

## Notes:
## - The data typically have a long left tail and so do not follow a normal
##   distribution. Therefore, instead of fitting and plotting a normal distribution
##   as I normally do for such histograms, I simply illustrate the median. The mean
##   isn't a useful statistic to illustrate given the asymmetry of the distributions.

## thresholds for defining unexpressed and expressed genes
expressionLimits = c (1, 10)

## thresholds for defining excluded and included cassette exons
psiLimits = c (0.05, 0.25)

## for all histograms, I create a matrix of panels, where the columns are:
## - 1: all
## - 2: upstreamExon
## - 3: all poison exons
## - 4: conserved poison exon
## - 5: unconserved poison exon
column2target_types = list (
  c ("upstreamExon", "conservedPE", "unconservedPE"),
  "upstreamExon",
  c ("conservedPE", "unconservedPE"),
  "conservedPE",
  "unconservedPE")

## plot 1
foo = function (d) {
  
  ## check that data has been properly grouped
  stopifnot (nrow (d) == nrow (d %>% select (id) %>% distinct()))

  origin = d$origin[1]
  time.reference = d$time.reference[1]
  time.comparison = d$time.comparison[1]
  replicate.reference = d$replicate.reference[1]
  replicate.comparison = d$replicate.comparison[1]

  for (statistic in c ("fc", "fc_norm")) {

    ## plot 1: histogram of control vs. targeting pgRNAs
    rows = 1
    columns = 1:length (column2target_types)
    p = list()
    for (row in rows) {
      for (column in columns) {

        ## panel index
        panel = length (p) + 1

        ## restrict to control pgRNAs and those targeting the current target_type
        d.plot = d %>% filter (control | target_type %in% column2target_types[[column]])

        d.plot = d.plot %>%
          select (control, "statistic" = !!statistic) %>%
          mutate ("value" = log2 (statistic),
                  "category" = factor (control)) ## I need a character vector instead of Boolean in order to use categoryLabels below
        d.plot$category = factor (d.plot$category, levels = c ("TRUE", "FALSE"))
        
        xlim = c (-5, 2.5)
        ylim = c (0, 0.7)
        
        p[[panel]] <- .plot.histogram.helper (
          d = d.plot,
          binwidth = 0.1, alpha = 0.75,
          xlim = xlim, xlabel = bquote (log[2]~("fold-change")),
          categoryLabels = c ("TRUE" = "control", "FALSE" = "targeting"))

        ## figure out the y coordinate limits that ggplot2 chooses, so that we can place
        ## our text annotations below accordingly
        ymax = min (max (ylim), sapply (ggplot_build (p[[panel]])$layout$panel_params, function (x) max (x$y.range)))
        
        ## illustrate the medians of the distributions
        ## - control
        M = as.numeric (d.plot %>% filter (control) %>% summarize (median (value)))
        p[[panel]] <- p[[panel]] +
          geom_vline (xintercept = M,
                      linetype = "dashed", color = palette.cb[which (levels (d.plot$category) == TRUE)]) +
          annotate (geom = "text",
                    x = M + 0.01 * (xlim[2] - xlim[1]), y = 0.9 * ymax,
                    label = paste0 ("mu[1/2] == ", signif (M, 2)), parse = TRUE,
                    hjust = 0, vjust = 1,
                    color = palette.cb[which (levels (d.plot$category) == TRUE)], size = 3)
        ## - targeting
        M = as.numeric (d.plot %>% filter (!control) %>% summarize (median (value)))
        p[[panel]] <- p[[panel]] +
          geom_vline (xintercept = M,
                      linetype = "dashed", color = palette.cb[which (levels (d.plot$category) == FALSE)]) +
          annotate (geom = "text",
                    x = M - 0.01 * (xlim[2] - xlim[1]), y = 0.9 * ymax,
                    label = paste0 ("mu[1/2] == ", signif (M, 2)), parse = TRUE,
                    hjust = 1, vjust = 1,
                    color = palette.cb[which (levels (d.plot$category) == FALSE)], size = 3)

        ## remove legend title
        p[[panel]] <- p[[panel]] + theme (legend.title = element_blank())

        ## add title stating the category of analyzed events
        label = switch (as.character (column),
                        "1" = "all exons",
                        "2" = "upstream exons",
                        "3" = "poison exons",
                        "4" = "conserved poison exons",
                        "5" = "unconserved poison exons")
        p[[panel]] <- p[[panel]] + ggtitle (label)
        p[[panel]] <- p[[panel]] + theme (plot.title = element_text (size = rel (0.8)))

      }
    }

    file = file.path (dir.figs, origin, replicate.reference, replicate.comparison,
                      collapse (time.reference, time.comparison, sep = "-vs-"),
                      paste ("histogram", statistic,
                             collapse ("control", "targeting", sep = "-vs-"), sep = "."))
    if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
      dir.create (dirname (file), showWarnings = FALSE, recursive = TRUE)
      file = device.open (plotdim = c (length (columns) * 3, length (rows) * 3), format = "pdf", file)
      multiplot (plotlist = p, cols = length (columns), byrow = TRUE)
      device.close (file)
    }

  }

  return (tibble())

}

d.pg %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison) %>%
  do (foo (.)) %>% invisible()

## plots 2, 3, and 4
foo = function (d) {

  ## check that data has been properly grouped
  stopifnot (nrow (d) == nrow (d %>% select (event, target_type) %>% distinct()))

  origin = d$origin[1]
  time.reference = d$time.reference[1]
  time.comparison = d$time.comparison[1]
  replicate.reference = d$replicate.reference[1]
  replicate.comparison = d$replicate.comparison[1]

  for (statistic in c ("fc_mean", "fc_norm_mean")) {

    ## plot 2: density plots of unexpressed vs. expressed genes
    rows = 1
    columns = 1:length (column2target_types)
    p = list()
    for (row in rows) {
      for (column in columns) {

        ## panel index
        panel = length (p) + 1

        ## restrict to the current target_type
        d.plot = d %>% filter (target_type %in% column2target_types[[column]])

        ## categorize by gene expression
        d.plot = d.plot %>% filter (!is.na (gene))
        expressed = rep (NA, nrow (d.plot))
        expressed[d.plot$expr <= expressionLimits[1]] = FALSE
        expressed[d.plot$expr >= expressionLimits[2]] = TRUE
        d.plot = d.plot %>% mutate ("expressed" = expressed) %>%
          filter (!is.na (expressed))

        d.plot = d.plot %>%
          select (expressed, "statistic" = !!statistic) %>%
          mutate ("value" = log2 (statistic),
                  "category" = factor (expressed)) ## I need a character vector instead of Boolean in order to use categoryLabels below
        d.plot$category = factor (d.plot$category, levels = c ("FALSE", "TRUE"))

        xlim = quantile (d.plot$value, c (0.025, 0.975))
        
        p[[panel]] <- .plot.histogram.helper (
          d = d.plot,
          histogram = FALSE, density = TRUE,
          adjust = 0.75,
          alpha = 0.75,
          xlim = xlim, xlabel = bquote (log[2]~("fold-change")),
          categoryLabels = c ("FALSE" = bquote (expr. <= .(expressionLimits[1])~TPM),
                              "TRUE" = bquote (expr. >= .(expressionLimits[2])~TPM)),
          parse = TRUE)

        ## figure out the y coordinate limits that ggplot2 chooses, so that we can place
        ## our text annotations below accordingly
        ymax = sapply (ggplot_build (p[[panel]])$layout$panel_params, function (x) max (x$y.range))
        
        ## illustrate the medians of the distributions
        ## - expressed
        M = as.numeric (d.plot %>% filter (expressed) %>% summarize (median (value)))
        p[[panel]] <- p[[panel]] +
          geom_vline (xintercept = M,
                      linetype = "dashed", color = palette.cb[which (levels (d.plot$category) == TRUE)]) +
          annotate (geom = "text",
                    x = M - 0.01 * (xlim[2] - xlim[1]), y = 0.8 * ymax,
                    label = paste0 ("mu[1/2] == ", signif (M, 2)), parse = TRUE,
                    hjust = 1, vjust = 1,
                    color = palette.cb[which (levels (d.plot$category) == TRUE)], size = 3)
        ## - unexpressed
        M = as.numeric (d.plot %>% filter (!expressed) %>% summarize (median (value)))
        p[[panel]] <- p[[panel]] +
          geom_vline (xintercept = M,
                      linetype = "dashed", color = palette.cb[which (levels (d.plot$category) == FALSE)]) +
          annotate (geom = "text",
                    x = M + 0.01 * (xlim[2] - xlim[1]), y = 0.8 * ymax,
                    label = paste0 ("mu[1/2] == ", signif (M, 2)), parse = TRUE,
                    hjust = 0, vjust = 1,
                    color = palette.cb[which (levels (d.plot$category) == FALSE)], size = 3)

        ## illustrate a p-value for a difference
        if (length (unique (d.plot$expressed)) > 1) {
          pval = wilcox.test (x = d.plot %>% filter (expressed) %>% pull (value),
                              y = d.plot %>% filter (!expressed) %>% pull (value))$p.value
        } else {
          pval = NA
        }
        p[[panel]] <- p[[panel]] +
          annotate (geom = "text",
                    x = xlim[1], y = 0.95 * ymax,
                    label = paste0 ("italic (p) == ", signif (pval, 3)), parse = TRUE,
                    hjust = 0, vjust = 1,
                    color = palette.cb[which (levels (d.plot$category) == TRUE)], size = 3)

        ## remove legend title
        p[[panel]] <- p[[panel]] + theme (legend.title = element_blank())

        ## add title stating the category of analyzed events
        label = switch (as.character (column),
                        "1" = "all exons",
                        "2" = "upstream exons",
                        "3" = "poison exons",
                        "4" = "conserved poison exons",
                        "5" = "unconserved poison exons")
        p[[panel]] <- p[[panel]] + ggtitle (label)
        p[[panel]] <- p[[panel]] + theme (plot.title = element_text (size = rel (0.8)))

      }
    }

    file = file.path (dir.figs, origin, replicate.reference, replicate.comparison,
                      collapse (time.reference, time.comparison, sep = "-vs-"),
                      paste ("density", statistic,
                             collapse ("unexpressed", "expressed", sep = "-vs-"), sep = "."))
    if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
      dir.create (dirname (file), showWarnings = FALSE, recursive = TRUE)
      file = device.open (plotdim = c (length (columns) * 3, length (rows) * 3), format = "pdf", file)
      multiplot (plotlist = p, cols = length (columns), byrow = TRUE)
      device.close (file)
    }

    ## plots 3 and 4: density plots of excluded vs. included exons
    ## create a matrix of panels, where columns are target_type
    for (useMuhlemann in c (TRUE, FALSE)) {
      for (perturbation in c ("basal", "nmd_inhibited")) {
        rows = 1
        columns = 1:length (column2target_types)
        p = list()
        for (row in rows) {
          for (column in columns) {

            ## panel index
            panel = length (p) + 1

            ## restrict to the current target_type
            d.plot = d %>% filter (target_type %in% column2target_types[[column]])

            ## categorize by cassette exon inclusion
            d.plot = d.plot %>% filter (!is.na (gene))
            included = rep (NA, nrow (d.plot))
            included[(d.plot$expr >= expressionLimits[2]) & (d.plot %>% pull (paste (if (useMuhlemann) "psi_muhlemann" else "psi", perturbation, sep = ".")) <= psiLimits[1])] = FALSE
            included[(d.plot$expr >= expressionLimits[2]) & (d.plot %>% pull (paste (if (useMuhlemann) "psi_muhlemann" else "psi", perturbation, sep = ".")) >= psiLimits[2])] = TRUE
            d.plot = d.plot %>% mutate ("included" = included) %>%
              filter (!is.na (included))

            d.plot = d.plot %>%
              select (included, "statistic" = !!statistic) %>%
              mutate ("value" = log2 (statistic),
                      "category" = factor (included)) ## I need a character vector instead of Boolean in order to use categoryLabels below
            d.plot$category = factor (d.plot$category, levels = c ("FALSE", "TRUE"))

            ## catch the case of no data
            if (!nrow (d.plot)) {
              p[[panel]] <- ggplot()
              next
            }

            xlim = quantile (d.plot$value, c (0.025, 0.975))

            p[[panel]] <- .plot.histogram.helper (
              d = d.plot,
              histogram = FALSE, density = TRUE,
              adjust = 0.75,
              alpha = 0.75,
              xlim = xlim, xlabel = bquote (log[2]~("fold-change")),
              categoryLabels = c ("FALSE" = bquote (psi <= .(psiLimits[1])),
                                  "TRUE" = bquote (psi >= .(psiLimits[2]))),
              parse = TRUE)
            
            ## figure out the y coordinate limits that ggplot2 chooses, so that we can place
            ## our text annotations below accordingly
            ymax = sapply (ggplot_build (p[[panel]])$layout$panel_params, function (x) max (x$y.range))
            
            ## illustrate the medians of the distributions
            ## - included
            M = as.numeric (d.plot %>% filter (included) %>% summarize (median (value)))
            p[[panel]] <- p[[panel]] +
              geom_vline (xintercept = M,
                          linetype = "dashed", color = palette.cb[which (levels (d.plot$category) == TRUE)]) +
              annotate (geom = "text",
                        x = M - 0.01 * (xlim[2] - xlim[1]), y = 0.8 * ymax,
                        label = paste0 ("mu[1/2] == ", signif (M, 2)), parse = TRUE,
                        hjust = 1, vjust = 1,
                        color = palette.cb[which (levels (d.plot$category) == TRUE)], size = 3)
            ## - excluded
            M = as.numeric (d.plot %>% filter (!included) %>% summarize (median (value)))
            p[[panel]] <- p[[panel]] +
              geom_vline (xintercept = M,
                          linetype = "dashed", color = palette.cb[which (levels (d.plot$category) == FALSE)]) +
              annotate (geom = "text",
                        x = M + 0.01 * (xlim[2] - xlim[1]), y = 0.8 * ymax,
                        label = paste0 ("mu[1/2] == ", signif (M, 2)), parse = TRUE,
                        hjust = 0, vjust = 1,
                        color = palette.cb[which (levels (d.plot$category) == FALSE)], size = 3)

            ## illustrate a p-value for a difference
            if (length (unique (d.plot$included)) > 1) {
              pval = wilcox.test (x = d.plot %>% filter (included) %>% pull (value),
                                  y = d.plot %>% filter (!included) %>% pull (value))$p.value
            } else {
              pval = NA
            }
            p[[panel]] <- p[[panel]] +
              annotate (geom = "text",
                        x = xlim[1], y = 0.95 * ymax,
                        label = paste0 ("italic (p) == ", signif (pval, 3)), parse = TRUE,
                        hjust = 0, vjust = 1,
                        color = palette.cb[which (levels (d.plot$category) == TRUE)], size = 3)

            ## remove legend title
            p[[panel]] <- p[[panel]] + theme (legend.title = element_blank())

            ## add title stating the category of analyzed events
            label = switch (as.character (column),
                            "1" = "all exons",
                            "2" = "upstream exons",
                            "3" = "poison exons",
                            "4" = "conserved poison exons",
                            "5" = "unconserved poison exons")
            p[[panel]] <- p[[panel]] + ggtitle (label)
            p[[panel]] <- p[[panel]] + theme (plot.title = element_text (size = rel (0.8)))

          }
        }

        file = file.path (dir.figs, origin, replicate.reference, replicate.comparison,
                          collapse (time.reference, time.comparison, sep = "-vs-"),
                          paste ("density", statistic,
                                 collapse ("excluded", "included", sep = "-vs-"),
                                 collapse ("perturbation", perturbation, sep = "_"),
                                 collapse ("useMuhlemann", useMuhlemann, sep = "_"), sep = "."))
        if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
          dir.create (dirname (file), showWarnings = FALSE, recursive = TRUE)
          file = device.open (plotdim = c (length (columns) * 3, length (rows) * 3), format = "pdf", file)
          multiplot (plotlist = p, cols = length (columns), byrow = TRUE)
          device.close (file)
        }
      }
    }

  }

  return (tibble())

}

d.targets %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison) %>%
  do (foo (.)) %>% invisible()

#####################################################################
## Heat maps: per-target
#####################################################################

## Goals:
## - Compute heat map using the most variable targets.

## Conclusions:
## -

for (statistic in c ("pval", "fc_mean", "fc_norm_mean")) {
  for (topn in c (100, 250, 500)) {

    ## we need the data in wide rather than long format
    ## because the conversion is a pain, simply read in the table that we wrote
    ## to disk earlier
    file = file.path (dir.figs,
                      collapse ("fc", "per-target", "table", sep = "."))
    d.plot = read_tsv (file) %>%
      select (starts_with (statistic))
    colnames (d.plot) = gsub (paste0 (statistic, "\\|"), "", colnames (d.plot))
    ## remove columns for the reserved replicate "all"
    indices = grep ("all|all", colnames (d.plot))
    if (length (indices))
      d.plot = d.plot[, -indices]

    ## restrict to the n-most variable genes
    sds = apply (d.plot %>% as.matrix, 1, sd)
    names (sds) = 1:length (sds)
    sds = rev (sort (sds))
    indices = as.integer (names (sds)[1:topn])
    d.plot = d.plot[indices, ]
    
    file = file.path (dir.figs,
                      paste ("heatmap", "targets",
                             collapse ("statistic", statistic, sep = "_"),
                             collapse ("topn", topn, sep = "_"), sep = "."))
    if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
      .plot.heatmap (d = as.matrix (d.plot),
                     method = "correlation", hclustMethod = "ward.D2",
                     zscoreNormalization = TRUE,
                     plotdim = c (5, 5), format = "pdf", file = file)
    }

  }
}

#####################################################################
## Box plots of pgRNAs targeting events of interest
#####################################################################

## Goals:
## - Create box plots showing fold-changes for control pgRNAs vs.
##   pgRNAs targeting targets of interest.

# d.targets %>% filter (origin == "HeLa" & replicate.reference == "all") %>% filter (target_type == "conservedPE") %>% arrange (fdr) %>% print (n = 100)

events = c (
  ## SMG1
  "se@16:18911366:18937272:-|16:18908278:18937272:-",
  ## CPSF4
  "se@7:99050063:99051227:+|7:99050063:99051589:+",
  ## SRSF3
  "se@6:36566760:36567598:+|6:36566760:36568929:+",
  ## SNRNP70
  "se@19:49604728:49605371:+|19:49604728:49607891:+",
  ## SMNDC1
  "se@10:112060494:112063226:-|10:112058548:112063226:-",
  ## TRA2B
  "se@3:185649640:185655613:-|3:185644522:185655613:-",
  ## SRSF7
  "se@2:38976315:38976671:-|2:38975795:38976671:-",
  "se@2:38976488:38976671:-|2:38975795:38976671:-",
  ## MAX
  "se@14:65544146:65544631:-|14:65543381:65544631:-",
  ## SIRT1
  "se@10:69651312:69665920:+|10:69651312:69666547:+")

for (event in events) {

  for (statistic in c ("fc", "fc_norm")) {

    p = list()
    
    ## get available target_type values for this event
    target_types = d.pg %>% filter (event == !!event) %>% select (target_type) %>% distinct() %>% collect() %>% .[[1]]
    ## sort in defined order so that panels are consistent between figures
    target_types = target_types[na.omit (match (c ("upstreamExon", "conservedPE", "unconservedPE"), target_types))]
    for (target_type in target_types) {

      ## panel index
      panel = length (p) + 1

      ## get fold-changes for control pgRNAs and pgRNAs for the target of interest
      d.plot = d.pg %>% filter (control | (event == !!event & target_type == !!target_type))

      ## use geneName as category (get it from d.targets)
      geneName = d.targets %>%
        filter (event == !!event & target_type == !!target_type) %>%
        select (geneName) %>% distinct() %>% collect() %>% .[[1]]
      d.plot = d.plot %>% mutate ("geneName" = geneName)
      d.plot$geneName[d.plot$control == TRUE] = "control"

      ## create column with compact specification of the comparison
      d.plot = d.plot %>%
        mutate ("comparison" = paste (origin, time.reference, time.comparison, replicate.reference, replicate.comparison, sep = "|"))
      d.plot$geneName = factor (d.plot$geneName, levels = c ("control", geneName))

      d.plot = d.plot %>%
        select (control, target_type,
                comparison, "statistic" = !!statistic, geneName)

      ## box plot
      p[[panel]] <- ggplot (
        data = d.plot,
        mapping = aes (x = factor (comparison), y = log2 (statistic),
                       fill = geneName))
      p[[panel]] <- p[[panel]] + geom_boxplot (outlier.shape = NA)

      ylim = c (-3, 3)
      p[[panel]] <- p[[panel]] + coord_cartesian (ylim = ylim)
      
      ## use custom colors
      p[[panel]] <- p[[panel]] + scale_fill_manual (values = as.vector (palette.cb[1:2]))

      ## axis labels
      p[[panel]] <- p[[panel]] + theme (axis.text.x = element_text (angle = 45, hjust = 1, vjust = 1))
      p[[panel]] <- p[[panel]] + xlab (NULL) + ylab (bquote (log[2]~("fold-change")))
      ## plot title
      p[[panel]] <- p[[panel]] + ggtitle (bquote (italic (.(geneName)) ~
                                                    .(switch ((d.plot %>% filter (!control))$target_type[1],
                                                              "upstreamExon" = "(upstream exon)",
                                                              "conservedPE" = "(conserved poison exon)",
                                                              "unconservedPE" = "(unconserved poison exon)"))))

    }

    file = file.path (dir.figs,
                      paste ("box", statistic,
                             geneName, str_replace_all (event, "\\||\\:", "_"), sep = "."))
    if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
      file = device.open (plotdim = c (length (p) * (2 + (length (unique (d.plot$comparison)) / 4)), 3.5), format = "pdf", file)
      multiplot (plotlist = p, cols = 2)
      device.close (file)
    }

  }

}

#####################################################################
## Ranked plot of p-values, FDRs, and fold-change
#####################################################################

## Goals:
## - Create ranked plots of p-values, FDRs, or fold-change that are optionally
##   stratified by target_type.

## note samples to highlight
targets_to_highlight = c(
    ## U2AF1 poison exon
    "se@21:44521542:44524425:-|21:44520629:44524425:-@inc",
    ## U2AF1 upstream exon
    "se@21:44521542:44524425:-|21:44520629:44524425:-@upstreamExon",
    ## SRSF3 poison exon
    "se@6:36566760:36567598:+|6:36566760:36568929:+@inc",
    ## SRSF3 upstream exon
    "se@6:36566760:36567598:+|6:36566760:36568929:+@upstreamExon",
    ## SNRNP70 poison exon
    "se@19:49604728:49605371:+|19:49604728:49607891:+@inc",
    ## SNRNP70 upstream exon
    "se@19:49604728:49605371:+|19:49604728:49607891:+@upstreamExon",
    ## MBNL1 upstream exon 1
    "se@3:152163328:152173331:+|3:152163328:152174056:+@upstreamExon",
    ## MBNL1 upstream exon 2
    "se@3:152174150:152175798:+|3:152174150:152177060:+@upstreamExon",
    ## SMNDC1 upstream exon
    "se@10:112060494:112063226:-|10:112058548:112063226:-@inc",
    ## SMNDC1 poison exon
    "se@10:112060494:112063226:-|10:112058548:112063226:-@upstreamExon"
    ## SIRT1 upstream exon
    ##"se@10:69651312:69665920:+|10:69651312:69666547:+@upstreamExon",
    ## SIRT1 poison exon
    ##"se@10:69651312:69665920:+|10:69651312:69666547:+@inc",
)

## target_type.grouping: list of vectors of target_type to group together
##   A separate ranked plot is created for each entry in target_type.grouping.
## targets.highlight: targets to highlight in red
foo = function (d,
                target_type.grouping, targets.highlight=NULL,
                ## collapse all p-values/FDRs/fold-changes below minpval/minfdr/minfc
                minpval=1e-6, minfdr=1e-6, minfc=1e-1) {
  
  ## check that data has been properly grouped
  stopifnot (nrow (d) == nrow (d %>% select (event, target_type) %>% distinct()))

  origin = d$origin[1]
  time.reference = d$time.reference[1]
  time.comparison = d$time.comparison[1]
  replicate.reference = d$replicate.reference[1]
  replicate.comparison = d$replicate.comparison[1]

  for (statistic in c ("pval", "fdr", "fc_mean", "fc_norm_mean")) {

    d.plot = d %>%
      select (geneName, target, target_type, "statistic" = statistic)

    ## categorize data by target_type
    d.plot = d.plot %>% filter (target_type %in% unlist (target_type.grouping))
    category = rep (NA, nrow (d.plot))
    for (i in 1:length (target_type.grouping))
      category[d.plot$target_type %in% target_type.grouping[[i]]] = names (target_type.grouping)[i]
    stopifnot (!any.na (d.plot$target_type))
    d.plot = d.plot %>% mutate ("category" = category)
    d.plot$category = factor (d.plot$category, levels = names (target_type.grouping))

    ## compute relative rank for each grouping
    d.plot = d.plot %>%
      group_by (category) %>%
      arrange (statistic) %>% mutate ("x" = 100 * (1:n() / n())) %>%
      ungroup()

    ## keep the y axis bounded
    if (statistic == "pval") {
      d.plot$statistic[d.plot$statistic < minpval] = minpval
    } else if (statistic == "fdr") {
      d.plot$statistic[d.plot$statistic < minfdr] = minfdr
    } else if (statistic %in% c ("fc_mean", "fc_norm_mean")) {
      d.plot$statistic[d.plot$statistic < minfc] = minfc
    }

    if (statistic %in% c ("pval", "fdr")) {
      d.plot = d.plot %>% mutate ("y" = -log10 (statistic))
    } else if (statistic %in% c ("fc_mean", "fc_norm_mean")) {
      d.plot = d.plot %>% mutate ("y" = -log2 (statistic))
    }

    p <- ggplot (
      data = d.plot,
      mapping = aes (x = x, y = y,
                     color = category))

    ## plot all targets
    p <- p + geom_point (size = 0.25)

    ## highlight requested targets with parent gene name
    p <- p + geom_point (
               data = d.plot %>% filter (target %in% targets.highlight),
               size = 0.75, color = "red")
    p <- p + geom_text (
               data = d.plot %>% filter (target %in% targets.highlight),
               aes (label = geneName),
               hjust = 1, nudge_x = -0.01,
               size = 2, color = "red")

    ## use custom colors
    p <- p + scale_color_manual (values = as.vector (palette.cb))

    ## remove legend title
    p <- p + theme (legend.title = element_blank())

    ## axis labels
    p <- p + theme (axis.text.x = element_text (angle = 45, hjust = 1, vjust = 1))
    p <- p + xlab ("relative target rank (%)")
    p <- p + ylab (switch (statistic,
                           "pval" = bquote (-log[10]~(italic (p)*'-'*value)),
                           "fdr" = bquote (-log[10]~(FDR)),
                           "fc_mean" = bquote (-log[2]~("fold-change")),
                           "fc_norm_mean" = bquote (-log[2]~("normalized fold-change"))))

    file = file.path (dir.figs, origin, replicate.reference, replicate.comparison,
                      collapse (time.reference, time.comparison, sep = "-vs-"),
                      paste ("rank", statistic, sep = "."))
    if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
      file = device.open (plotdim = c (3, 3), format = "pdf", file)
      print (p)
      device.close (file)
    }

  }

  return (tibble())

}

d.targets %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison) %>%
  do (foo (., target_type.grouping = list (
                "upstream exon" = "upstreamExon",
                "conserved poison exon" = "conservedPE",
                "unconserved poison exon" = "unconservedPE"),
           targets.highlight = targets_to_highlight)) %>%
  invisible()

#####################################################################
## Volcano plots (fold-change vs. p-value or FDR)
#####################################################################

## Goals:
## - Create volcano plots of fold-change vs. p-value or FDR, stratified by target_type.

## Notes:
## - Volcano plots based on unnormalized fold-changes are not very meaningful,
##   as p-values and FDRs are computed based on normalized fold-changes.

## use the same panel layout as for the histograms
stopifnot (exists ("column2target_types"))

## targets.highlight: targets to highlight in red
foo = function (d,
                targets.highlight=NULL,
                ## collapse all p-values/FDRs below minpval/minfdr
                minpval=1e-4, minfdr=1e-4) {

  ## check that data has been properly grouped
  stopifnot (nrow (d) == nrow (d %>% select (event, target_type) %>% distinct()))

  origin = d$origin[1]
  time.reference = d$time.reference[1]
  time.comparison = d$time.comparison[1]
  replicate.reference = d$replicate.reference[1]
  replicate.comparison = d$replicate.comparison[1]

  for (statistic1 in c ("fc_mean", "fc_norm_mean")) {
    for (statistic2 in c ("pval", "fdr")) {

      rows = 1
      columns = 1:length (column2target_types)
      p = list()
      for (row in rows) {
        for (column in columns) {

          ## panel index
          panel = length (p) + 1

          d.plot = d %>%
            select (geneName, target, target_type,
                    "statistic1" = !!statistic1, "statistic2" = !!statistic2)

          ## categorize data by target_type
          d.plot = d.plot %>% filter (target_type %in% column2target_types[[column]])

          ## keep the y axis bounded
          if (statistic2 == "pval") {
            d.plot$statistic2[d.plot$statistic2 < minpval] = minpval
          } else if (statistic2 == "fdr") {
            d.plot$statistic2[d.plot$statistic2 < minfdr] = minfdr
          }
          
          d.plot = d.plot %>% mutate ("x" = statistic1,
                                      "y" = statistic2)

          p[[panel]] <- ggplot (
            data = d.plot,
            mapping = aes (x = log2 (x), y = -log10 (y)))

          ## plot all targets
          p[[panel]] <- p[[panel]] + geom_point (size = 0.1, color = "black")

          ## highlight requested targets with parent gene name
          p[[panel]] <- p[[panel]] + geom_point (
                                       data = d.plot %>% filter (target %in% targets.highlight),
                                       size = 0.75, color = "red")
          p[[panel]] <- p[[panel]] + geom_text (
                                       data = d.plot %>% filter (target %in% targets.highlight),
                                       aes (label = geneName),
                                       hjust = 1, nudge_x = -0.01,
                                       size = 2, color = "red")

          p[[panel]] <- p[[panel]] + coord_cartesian (xlim = c (-3, 1))

          ## axis labels
          p[[panel]] <- p[[panel]] + xlab (bquote (log[2]~("fold-change")))
          p[[panel]] <- p[[panel]] + ylab (switch (statistic2,
                                                   "pval" = bquote (-log[10]~(italic (p)*'-'*value)),
                                                   "fdr" = bquote (-log[10]~(FDR))))

          ## add title stating the category of analyzed events
          label = switch (as.character (column),
                          "1" = "all exons",
                          "2" = "upstream exons",
                          "3" = "poison exons",
                          "4" = "conserved poison exons",
                          "5" = "unconserved poison exons")
          p[[panel]] <- p[[panel]] + ggtitle (label)
          p[[panel]] <- p[[panel]] + theme (plot.title = element_text (size = rel (0.8)))

        }
      }
      
      file = file.path (dir.figs, origin, replicate.reference, replicate.comparison,
                        collapse (time.reference, time.comparison, sep = "-vs-"),
                        paste ("volcano", collapse (statistic1, statistic2, sep = "-vs-"), sep = "."))
      if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
        file = device.open (plotdim = c (length (columns) * 3, length (rows) * 3), format = "pdf", file)
        multiplot (plotlist = p, cols = length (columns), byrow = TRUE)
        device.close (file)
      }

    }
  }

  return (tibble())

}

d.targets %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison) %>%
  do (foo (.,
           targets.highlight = c ("se@10:69651312:69665920:+|10:69651312:69666547:+@inc",
                                  "se@10:112060494:112063226:-|10:112058548:112063226:-@inc",
                                  "se@10:112060494:112063226:-|10:112058548:112063226:-@upstreamExon"))) %>%
  invisible()

#####################################################################
## CDF plot of fractions of depleted targets
#####################################################################

## Goals:
## - Create CDF plots of fold-change for a given p-value threshold that are optionally
##   stratified by target_type.

## Notes:
## - Analysis is restricted to targets with pval <= maxpval.

## target_type.grouping: list of vectors of target_type to group together
##   A separate CDF plot is created for each entry in target_type.grouping.
foo = function (d,
                target_type.grouping) {

  ## check that data has been properly grouped
  stopifnot (nrow (d) == nrow (d %>% select (event, target_type) %>% distinct()))

  origin = d$origin[1]
  time.reference = d$time.reference[1]
  time.comparison = d$time.comparison[1]
  replicate.reference = d$replicate.reference[1]
  replicate.comparison = d$replicate.comparison[1]

  for (statistic in c ("fc_mean", "fc_norm_mean")) {
    for (maxpval in c (1e-2, 1e-3)) {

      d.plot = d %>%
        select (geneName, target, target_type,
                pval, "statistic" = !!statistic)

      ## categorize data by target_type
      d.plot = d.plot %>% filter (target_type %in% unlist (target_type.grouping))
      category = rep (NA, nrow (d.plot))
      for (i in 1:length (target_type.grouping))
        category[d.plot$target_type %in% target_type.grouping[[i]]] = names (target_type.grouping)[i]
      stopifnot (!any.na (d.plot$target_type))
      d.plot = d.plot %>% mutate ("category" = category)
      d.plot$category = factor (d.plot$category, levels = names (target_type.grouping))

      ## enforce maxpval
      d.plot = d.plot %>% filter (pval <= maxpval)

      p <- ggplot (
        data = d.plot,
        mapping = aes (x = -log2 (statistic),
                       color = category))

      p <- p + stat_ecdf (geom = "step")

      p <- p + coord_cartesian (xlim = c (0, 2))

      ## use custom colors
      p <- p + scale_color_manual (values = as.vector (palette.cb))

      ## remove legend title
      p <- p + theme (legend.title = element_blank())

      ## axis labels
      p <- p + theme (axis.text.x = element_text (angle = 45, hjust = 1, vjust = 1))
      p <- p + xlab (bquote (-log[2]~("fold-change")))
      p <- p + ylab ("cumulative density")

      file = file.path (dir.figs, origin, replicate.reference, replicate.comparison,
                        collapse (time.reference, time.comparison, sep = "-vs-"),
                        paste ("cdf", statistic,
                               collapse ("maxpval", maxpval, sep = "_"), sep = "."))
      if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
        file = device.open (plotdim = c (3, 3), format = "pdf", file)
        print (p)
        device.close (file)
      }

    }
  }

  return (tibble())

}

d.targets %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison) %>%
  do (foo (., target_type.grouping = list (
                "upstream exon" = "upstreamExon",
                "conserved poison exon" = "conservedPE",
                "unconserved poison exon" = "unconservedPE"))) %>%
  invisible()

#####################################################################
## Scatter plots comparing upstream and poison exon depletion
#####################################################################

## Goals:
## - Create scatter plots comparing fold-change for matched upstream and poison exons.

## Notes:
## - I think that matched upstreamExon targets were chosen for conservedPE targets,
##   but not for unconservedPE targets.

## fc.highlight: highlight poison exon targets with abs (log2 (fc)) >= abs (log2 (fc.highlight)) in blue
## events.highlight: highlight targets matched to these events in red
foo = function (d,
                fc.highlight=4, events.highlight=NULL,
                ## collapse all fold-changes below minfc
                minfc=1e-1) {

  ## check that data has been properly grouped
  stopifnot (nrow (d) == nrow (d %>% select (event, target_type) %>% distinct()))

  origin = d$origin[1]
  time.reference = d$time.reference[1]
  time.comparison = d$time.comparison[1]
  replicate.reference = d$replicate.reference[1]
  replicate.comparison = d$replicate.comparison[1]

  for (statistic in c ("fc_mean", "fc_norm_mean")) {

    d.plot = d %>%
      select (geneName, event, target_type, "statistic" = statistic)

    ## restrict to events with upstreamExon as well as either conservedPE and/or unconservedPE
    events = intersect (
      d.plot %>% filter (target_type == "upstreamExon") %>% select (event) %>% distinct() %>% collect() %>% .[[1]],
      d.plot %>% filter (target_type %in% c ("unconservedPE", "conservedPE")) %>% select (event) %>% distinct() %>% collect() %>% .[[1]])
    d.plot = d.plot %>% filter (event %in% events)

    d.plot = inner_join (
      d.plot %>% filter (target_type == "upstreamExon") %>% select (geneName, event, "x" = statistic),
      d.plot %>% filter (target_type %in% c ("unconservedPE", "conservedPE")) %>% select (geneName, event, "y" = statistic, target_type),
      by = c ("geneName", "event"))

    ## keep the axes bounded
    d.plot$x[d.plot$x < minfc] = minfc
    d.plot$y[d.plot$y < minfc] = minfc

    ## compute log fold-change
    d.plot = d.plot %>%
      mutate ("x" = log2 (x), "y" = log2 (y))

    p <- ggplot (
      data = d.plot,
      mapping = aes (x = x, y = y))

    ## plot all targets
    p <- p + geom_point (size = 0.25, color = "black")

    ## highlight poison exon targets with abs (fold-change) >= fc.highlight
    p <- p + geom_point (
               data = d.plot %>% filter (abs (y) >= abs (log2 (fc.highlight))),
               size = 0.75, color = "blue")
    p <- p + geom_text_repel (
               data = d.plot %>% filter (abs (y) >= abs (log2 (fc.highlight))),
               aes (label = geneName),
               hjust = 0,
               size = 2, color = "blue")

    ## highlight requested events with parent gene name
    p <- p + geom_point (
               data = d.plot %>% filter (event %in% events.highlight),
               size = 0.75, color = "red")
    p <- p + geom_text_repel (
               data = d.plot %>% filter (event %in% events.highlight),
               aes (label = geneName),
               hjust = 0,
               size = 2, color = "red")

    ## remove legend title
    p <- p + theme (legend.title = element_blank())

    ## axis labels
    p <- p + xlab (switch (statistic,
                           "fc_mean" = bquote (log[2]~("fc for upstream exon")),
                           "fc_norm_mean" = bquote (log[2]~("norm. fc for upstream exon"))))
    p <- p + ylab (switch (statistic,
                           "fc_mean" = bquote (log[2]~("fc for poison exon")),
                           "fc_norm_mean" = bquote (log[2]~("norm. fc for poison exon"))))

    file = file.path (dir.figs, origin, replicate.reference, replicate.comparison,
                      collapse (time.reference, time.comparison, sep = "-vs-"),
                      paste ("scatter", collapse (c ("upstream", "poison"), sep = "-vs-"),
                             statistic, sep = "."))
    if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
      file = device.open (plotdim = c (3, 3), format = "pdf", file)
      print (p)
      device.close (file)
    }

  }

  return (tibble())

}

d.targets %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison) %>%
  do (foo (.,
           events.highlight = c ("se@10:69651312:69665920:+|10:69651312:69666547:+"))) %>%
  invisible()

#####################################################################
## Library-specific plots: Scatter plots comparing fold-change between cell lines
#####################################################################

## Goals:
## - Create scatter plot illustrating fold-changes for targets called as hits in PC9
##   cells in PC9 and HeLa cells.

## Conclusions:
## - There is an excellent correlation.

## Notes:
## - Targets are restricted to those in genes that are expressed in both cell lines.
## - This plot is based solely on hits called in PC9 cells because of the smaller
##   dynamic range observed in that cell line.

## call targets as hits if they have abs (log2 (fc_norm_mean)) >= abs (log2 (minfc))
## and fdr <= maxfdr
minfc = 1.25
maxfdr = 1e-2

## restrict analysis to targets in genes w/ at least this expression level
minExpression = 10

time.reference = "0d"
time.comparison = "14d"
replicate.reference = replicate.comparison = "all"

for (statistic in c ("fc_mean", "fc_norm_mean")) {

  file = file.path (dir.figs, collapse ("HeLa", "PC9", sep = "-vs-"), replicate.reference, replicate.comparison,
                    collapse (time.reference, time.comparison, sep = "-vs-"),
                    paste ("scatter", statistic, sep = "."))
  if (file.exists (paste (file, "pdf", sep = ".")) && !rebuild) next

  ## build tibbles for samples that we want to compare
  d.hela = d.targets %>%
    filter (origin == "HeLa" & time.reference == !!time.reference & time.comparison == !!time.comparison & replicate.reference == !!replicate.reference & replicate.comparison == !!replicate.comparison) %>%
    select (target, target_type, expr, fdr, "statistic" = statistic)
  d.pc9 = d.targets %>%
    filter (origin == "PC9" & time.reference == !!time.reference & time.comparison == !!time.comparison & replicate.reference == !!replicate.reference & replicate.comparison == !!replicate.comparison) %>%
    select (target, target_type, expr, fdr, "statistic" = statistic)

  ## identify targets in genes that are expressed in both cell lines and called as
  ## hits in PC9 cells
  targets = unique (intersect (
    d.hela %>% filter (expr >= minExpression) %>% select (target) %>% collect() %>% .[[1]],
    d.pc9 %>% filter (expr >= minExpression & abs (log2 (statistic)) >= abs (log2 (minfc)) & fdr <= maxfdr) %>% select (target) %>% collect() %>% .[[1]]))

  ## compute log fold-change
  d.plot = left_join (
    d.hela %>% filter (target %in% targets) %>% select (target, target_type, "x" = statistic) %>% mutate ("x" = log2 (x)),
    d.pc9 %>% filter (target %in% targets) %>% select (target, target_type, "y" = statistic) %>% mutate ("y" = log2 (y)),
    by = c ("target", "target_type"))

  ## categorize data by target_type
  d.plot = d.plot %>% mutate ("category" = factor (d.plot$target_type, levels = c ("upstreamExon", "conservedPE", "unconservedPE")))
 
  p <- ggplot (
    data = d.plot,
    mapping = aes (x = x, y = y,
                   color = category))

  ## plot all targets
  p <- p + geom_point (size = 0.5)

  ## show linear fit
  p <- p + stat_smooth (method = "lm", formula = y ~ x,
                        color = "blue", fill = "gray", size = 0.5)

  xlim = range (pretty (d.plot$x))
  ylim = range (pretty (d.plot$y))
  p <- p + coord_cartesian (xlim = xlim, ylim = ylim)

  ## show Pearson correlation
  r = cor (d.plot$x, d.plot$y, method = "pearson")
  p <- p + annotate ("text",
                     x = xlim[1] + 0.05 * (xlim[2] - xlim[1]), y = ylim[2] - 0.05 * (ylim[2] - ylim[1]),
                     label = paste0 ("italic (r) == ", signif (r, 2)), parse = TRUE,
                     hjust = 0, size = 3, color = "black")

  ## use custom colors
  p <- p + scale_color_manual (values = as.vector (palette.cb))

  ## remove legend title and put legend in bottom-right
  p <- p + theme (legend.title = element_blank(),
                  legend.justification = c (1, 0), legend.position = c (1, 0),
                  legend.margin = margin (0.25, 0.25, 0.25, 0.25, unit = "lines"))
  
  ## axis labels
  p <- p + xlab (switch (statistic,
                         "fc_mean" = bquote (log[2]~("fold-change for HeLa")),
                         "fc_norm_mean" = bquote (log[2]~("norm. fold-change for HeLa"))))
  p <- p + ylab (switch (statistic,
                         "fc_mean" = bquote (log[2]~("fold-change for PC9")),
                         "fc_norm_mean" = bquote (log[2]~("norm. fold-change for PC9"))))
  
  file = device.open (plotdim = c (3, 3), format = "pdf", file)
  print (p)
  device.close (file)

}

#####################################################################
## Basic statistics of hits
#####################################################################

## Goals:
## - Report core statistics for depleted and enriched hits.

## Conclusions:
## - 

## Notes:
## - 

## call targets as hits if they have abs (log2 (fc_norm_mean)) >= abs (log2 (minfc))
## and fdr <= maxfdr
minfc = 1.25
maxfdr = 1e-2

minExpression = 10

d = d.targets

## group as poison exon or upstream constitutive
grouping = rep (NA, nrow (d))
grouping[d$target_type == "upstreamExon"] = "upstream exon"
grouping[d$target_type %in% c ("conservedPE", "unconservedPE")] = "poison exon"
d = d %>% mutate ("grouping" = grouping)

## compute statistics on the numbers of targeted exons in expressed genes
d.targeted = d %>%
  filter (expr >= minExpression) %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison,
            grouping) %>% mutate ("n.targeted" = n()) %>% ungroup() %>%
  select (origin, time.reference, time.comparison, replicate.reference, replicate.comparison,
          grouping, n.targeted) %>% distinct()

## identify hits
d = d %>% filter (abs (log2 (fc_norm_mean)) >= abs (log2 (minfc)) & fdr <= maxfdr)

## classify each hit as enriched or depleted
category = rep (NA, nrow (d))
category[log2 (d$fc_norm_mean) < 0] = "depleted"
category[log2 (d$fc_norm_mean) > 0] = "enriched"
stopifnot (!any.na (category))
d = d %>% mutate ("category" = category)

## store information on numbers of targeted exons
d = d %>% left_join (d.targeted,
                     by = setdiff (colnames (d.targeted), "n.targeted"))

## compute statistics for the numbers of hits and fractions of targeted exons that are called as hits
d = d %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison,
            grouping, category) %>%
  mutate ("n" = n(), "f.targeted" = n / n.targeted) %>% ungroup() %>%
  select (origin, time.reference, time.comparison, replicate.reference, replicate.comparison,
          grouping, category, n, n.targeted, f.targeted) %>% distinct() %>%
  arrange (origin, time.reference, time.comparison, replicate.reference, replicate.comparison,
           grouping, category)

d %>% filter (time.reference == "0d" & time.comparison == "14d" & replicate.comparison %in% c ("all"))

d %>% filter (time.reference == "preinj" & time.comparison == "xeno_early" & replicate.comparison %in% c ("all"))

d %>% filter (time.reference == "preinj" & time.comparison == "xeno_late" & replicate.comparison %in% c ("all"))

#####################################################################
## Density plots of psi values w/ or w/o NMD inhibition
## Association between psi values and fold changes w/ or w/o NMD inhibition
#####################################################################

## Goals:
## - Create plots exploring the association between exon inclusion and fold-change:
##     plot 1: density plot of psi values in basal and NMD-inhibited states
##     plot 2: box plots comparing fold changes for targets stratified by psi value

## Conclusions:
## - There is a large and significant increase in poison exon inclusion following
##   NMD inhibition.
## - Poison exons that are frequently included exhibit significantly more dramatic
##   depletion than do poison exons that are preferentially excluded. However, this
##   difference is /only/ readily visible and significant when the Muhlemann data,
##   (2017/muhlemann.smg6_smg7_upf1_kd), rather than our own data, is used to query
##   the NMD-inhibited state, likely because the Muhlemann data has much better NMD
##   inhibition.

## Notes:
## - The density plots (plot 1) are restricted to targets satisfying these criteria:
##   * in expressed genes, defined as expression >= minExpression
## - The box plots (plot 2) are restricted to targets satisfying these criteria:
##   * in expressed genes, defined as expression >= minExpression
##   * with psi values > 0 and delta psi > 0 following NMD inhibition
##     I've empirically observed that targets with 0 psi values exhibit unusual
##     behavior, perhaps because they aren't bona fide poison exons.
##   Targets are grouped by psi values in the NMD-inhibited state.
## - I experimented with comparing psi values for targets that were or weren't called
##   as hits. I didn't observe clear differences in the psi value distribuions unless
##   I required that hits exhibit large fold-changes. This is consistent with
##   the association between psi value and fold-change being quantitative,
##   rather than qualitative.

## threshold for defining expressed genes
minExpression = 10

foo = function (d) {

  ## check that data has been properly grouped
  stopifnot (nrow (d) == nrow (d %>% select (event, target_type) %>% distinct()))

  origin = d$origin[1]
  time.reference = d$time.reference[1]
  time.comparison = d$time.comparison[1]
  replicate.reference = d$replicate.reference[1]
  replicate.comparison = d$replicate.comparison[1]

  ## plot 1: density plots of psi (basal vs. NMD-inhibited cells)
  for (useMuhlemann in c (TRUE, FALSE)) {

    rows = 1
    columns = 1:length (column2target_types)
    p = list()
    for (row in rows) {
      for (column in columns) {

        ## panel index
        panel = length (p) + 1

        d.plot = d %>%
          ## restrict to the current target_type
          filter (target_type %in% column2target_types[[column]]) %>%
          ## restrict to targets in expressed genes
          filter (expr >= minExpression)

        d.plot = d.plot %>%
          mutate ("basal" = if (useMuhlemann) psi_muhlemann.basal else psi.basal,
                  "nmd_inhibited" = if (useMuhlemann) psi_muhlemann.nmd_inhibited else psi.nmd_inhibited) %>%
          select (basal, nmd_inhibited) %>%
          gather (key = "category", value = "value") %>%
          mutate ("value" = 100 * value)
        d.plot$category = factor (d.plot$category, levels = c ("basal", "nmd_inhibited"))

        ## restrict to targets with defined psi values
        d.plot = d.plot %>% filter (!is.na (value))
        if (!nrow (d.plot)) {
          p[[panel]] <- ggplot()
          next
        }

        xlim = c (0, 100)
        
        ## plot density plots
        p[[panel]] <- .plot.histogram.helper (
          d = d.plot,
          histogram = FALSE, density = TRUE,
          adjust = 0.5,
          alpha = 0.75,
          xlim = xlim, xlabel = bquote (psi['NMD+']~"(%)"),
          categoryLabels = c ("basal" = "basal",
                              "nmd_inhibited" = "NMD-inhibited"))

        ## figure out the y coordinate limits that ggplot2 chooses, so that we can place
        ## our text annotations below accordingly
        ymax = sapply (ggplot_build (p[[panel]])$layout$panel_params, function (x) max (x$y.range))

        ## illustrate the medians of the distributions
        ## - NMD-inhibited
        M = as.numeric (d.plot %>% filter (category == "nmd_inhibited") %>% summarize (median (value)))
        p[[panel]] <- p[[panel]] +
          geom_vline (xintercept = M,
                      linetype = "dashed", color = palette.cb[which (levels (d.plot$category) == "nmd_inhibited")]) +
          annotate (geom = "text",
                    x = M - 0.01 * (xlim[2] - xlim[1]), y = 0.85 * ymax,
                    label = paste0 ("mu[1/2] == ", signif (M, 2)), parse = TRUE,
                    hjust = 0, vjust = 1,
                    color = palette.cb[which (levels (d.plot$category) == "nmd_inhibited")], size = 3)
        ## - basal
        M = as.numeric (d.plot %>% filter (category == "basal") %>% summarize (median (value)))
        p[[panel]] <- p[[panel]] +
          geom_vline (xintercept = M,
                      linetype = "dashed", color = palette.cb[which (levels (d.plot$category) == "basal")]) +
          annotate (geom = "text",
                    x = M + 0.01 * (xlim[2] - xlim[1]), y = 0.75 * ymax,
                    label = paste0 ("mu[1/2] == ", signif (M, 2)), parse = TRUE,
                    hjust = 0, vjust = 1,
                    color = palette.cb[which (levels (d.plot$category) == "basal")], size = 3)

        ## illustrate a p-value for a difference
        if (length (unique (d.plot$category)) > 1) {
          pval = wilcox.test (x = d.plot %>% filter (category == "basal") %>% pull (value),
                              y = d.plot %>% filter (category == "nmd_inhibited") %>% pull (value))$p.value
        } else {
          pval = NA
        }
        p[[panel]] <- p[[panel]] +
          annotate (geom = "text",
                    x = xlim[2], y = 0.5 * ymax,
                    label = paste0 ("italic (p) == ", signif (pval, 3)), parse = TRUE,
                    hjust = 1, vjust = 0,
                    color = palette.cb[which (levels (d.plot$category) == "nmd_inhibited")], size = 3)

        ## remove legend title
        p[[panel]] <- p[[panel]] + theme (legend.title = element_blank())

        ## add title stating the category of analyzed events
        label = switch (as.character (column),
                        "1" = "all exons",
                        "2" = "upstream exons",
                        "3" = "poison exons",
                        "4" = "conserved poison exons",
                        "5" = "unconserved poison exons")
        p[[panel]] <- p[[panel]] + ggtitle (label)
        p[[panel]] <- p[[panel]] + theme (plot.title = element_text (size = rel (0.8)))

      }
    }

    file = file.path (dir.figs, origin, replicate.reference, replicate.comparison,
                      collapse (time.reference, time.comparison, sep = "-vs-"),
                      paste ("density", "psi",
                             collapse ("basal", "nmd_inhibited", sep = "-vs-"),
                             collapse ("useMuhlemann", useMuhlemann, sep = "_"), sep = "."))
    if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
      dir.create (dirname (file), showWarnings = FALSE, recursive = TRUE)
      file = device.open (plotdim = c (length (columns) * 3, length (rows) * 3), format = "pdf", file)
      multiplot (plotlist = p, cols = length (columns), byrow = TRUE)
      device.close (file)
    }

  }

  ## plot 2: box plots comparing fold-change for targets binned by psi value
  ## rows: psi computed from basal or NMD-inhibited cells
  perturbations = c ("basal", "nmd_inhibited")
  for (statistic in c ("fc_mean", "fc_norm_mean")) {
    for (useMuhlemann in c (TRUE, FALSE)) {

      rows = 1:length (perturbations)
      columns = 1:length (column2target_types)
      p = list()
      for (row in rows) {
        for (column in columns) {

          ## panel index
          panel = length (p) + 1

          d.plot = d %>%
            ## restrict to the current target_type
            filter (target_type %in% column2target_types[[column]]) %>%
            ## restrict to targets in expressed genes
            filter (expr >= minExpression) %>%
            mutate ("basal" = if (useMuhlemann) psi_muhlemann.basal else psi.basal,
                    "nmd_inhibited" = if (useMuhlemann) psi_muhlemann.nmd_inhibited else psi.nmd_inhibited) %>%
            ## restrict to targets with defined psi values
            filter (!is.na (basal) & !is.na (nmd_inhibited)) %>%
            ## restrict to targets that are detectably included following NMD inhibition
            ## and exhibit increased inclusion following NMD inhibition
            filter (nmd_inhibited > 0 & nmd_inhibited > basal)
          
          d.plot = d.plot %>%
            select (target_type,
                    "statistic" = !!statistic,
                    "psi" = !!perturbations[row])

          ## restrict to targets with defined psi values
          d.plot = d.plot %>% filter (!is.na (psi))
          ## categorize by psi
          breaks = c (0, 0.05, 0.25, 1)
          d.plot = d.plot %>%
            mutate ("category" = cut (psi, breaks = breaks,
                                      include.lowest = TRUE)) %>%
            filter (!is.na (category))
          stopifnot (!any.na (d.plot$category))

          if (!nrow (d.plot) || (length (unique (d.plot$category)) < length (levels (d.plot$category)))) {
            p[[panel]] <- ggplot()
            next
          }

          p[[panel]] <- ggplot (
            data = d.plot,
            mapping = aes (x = category, y = log2 (statistic)))
          p[[panel]] <- p[[panel]] + geom_boxplot (outlier.shape = NA, notch = TRUE)

          ylim = quantile (log2 (d.plot$statistic), c (0.05, 0.95))
          p[[panel]] <- p[[panel]] + coord_cartesian (ylim = ylim)

          ## axis labels
          p[[panel]] <- p[[panel]] + theme (axis.text.x = element_text (angle = 45, hjust = 1, vjust = 1))
          p[[panel]] <- p[[panel]] +
            xlab (switch (perturbations[row],
                          "basal" = bquote (psi["NMD+"]^"basal"),
                          "nmd_inhibited" = bquote (psi["NMD+"]^"NMD-inhibited"))) +
            ylab (switch (statistic,
                          "fc_mean" = bquote (log[2]~("fold-change")),
                          "fc_norm_mean" = bquote (log[2]~("norm. fold-change"))))

          ## figure out the x coordinate limits that ggplot2 chooses, so that we can place
          ## our text annotations below accordingly
          xmax = sapply (ggplot_build (p[[panel]])$layout$panel_params, function (x) max (x$x.range))

          ## illustrate a p-value for a difference
          if (length (unique (d.plot$category)) > 1) {
            pval = wilcox.test (x = d.plot %>% filter (category == levels (d.plot$category)[1]) %>% pull (statistic),
                                y = d.plot %>% filter (category == levels (d.plot$category)[length (levels (d.plot$category))]) %>% pull (statistic))$p.value
          } else {
            pval = NA
          }
          p[[panel]] <- p[[panel]] +
            annotate (geom = "text",
                      x = xmax, y = min (ylim),
                      label = paste0 ("italic (p) == ", signif (pval, 3)), parse = TRUE,
                      hjust = 1, vjust = 0,
                      size = 3)

          ## add title stating the category of analyzed events
          label = switch (as.character (column),
                          "1" = "all exons",
                          "2" = "upstream exons",
                          "3" = "poison exons",
                          "4" = "conserved poison exons",
                          "5" = "unconserved poison exons")
          p[[panel]] <- p[[panel]] + ggtitle (label)
          p[[panel]] <- p[[panel]] + theme (plot.title = element_text (size = rel (0.8)))

        }
      }

      file = file.path (dir.figs, origin, replicate.reference, replicate.comparison,
                        collapse (time.reference, time.comparison, sep = "-vs-"),
                        paste ("box", statistic,
                               "by_psi",
                               collapse ("useMuhlemann", useMuhlemann, sep = "_"), sep = "."))
      if (!file.exists (paste (file, "pdf", sep = ".")) || rebuild) {
        dir.create (dirname (file), showWarnings = FALSE, recursive = TRUE)
        file = device.open (plotdim = c (length (columns) * 1.75, length (rows) * 3), format = "pdf", file)
        multiplot (plotlist = p, cols = length (columns), byrow = TRUE)
        device.close (file)
      }

    }
  }

  return (tibble())

}

d.targets %>%
  group_by (origin, time.reference, time.comparison, replicate.reference, replicate.comparison) %>%
  do (foo (.)) %>% invisible()
