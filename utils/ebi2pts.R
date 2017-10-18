#!/usr/bin/env Rscript

AUTHOR = "Mark Pinese <m.pinese@garvan.org.au>"
VERSION = "0.0.2 (18 Nov 2017)"


suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(RMySQL))


read_population_associations = function(assoc_file, ancestry_file, pop_regex = "european")
{
    ancestry = read.table(ancestry_file, sep = "\t", header = TRUE, comment.char = "", quote = "", row.names = NULL, stringsAsFactors = FALSE)
    colnames(ancestry) = c(colnames(ancestry)[-1], "dummy")
    ancestry = ancestry[,-ncol(ancestry)]

    count.popn.initial = daply(ancestry, .(STUDY.ACCCESSION), function(a) sum(a$NUMBER.OF.INDIVDUALS[a$STAGE == "initial" & grepl(pop_regex, tolower(a$BROAD.ANCESTRAL.CATEGORY))]))
    count.popn.replication = daply(ancestry, .(STUDY.ACCCESSION), function(a) sum(a$NUMBER.OF.INDIVDUALS[a$STAGE == "replication" & grepl(pop_regex, tolower(a$BROAD.ANCESTRAL.CATEGORY))]))

    assoc = read.table(assoc_file, sep = "\t", header = TRUE, comment.char = "", quote = "", stringsAsFactors = FALSE)
    assoc = assoc[,c("PUBMEDID", "STRONGEST.SNP.RISK.ALLELE", "PVALUE_MLOG", "P.VALUE..TEXT.", "OR.or.BETA", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "STUDY.ACCESSION")]
    assoc$pop = pop_regex
    assoc$pop.n.initial = count.popn.initial[assoc$STUDY.ACCESSION]
    assoc$pop.n.replication = count.popn.replication[assoc$STUDY.ACCESSION]

    # Calculate ln(Beta), discarding invalid values
    assoc$log_beta = log(assoc$OR.or.BETA)
    assoc = assoc[!(is.na(assoc$log_beta)) & is.finite(assoc$log_beta),,drop=FALSE]

    # Discard non-SNP associations, or rows with incomplete data (missing at
    # least one of rsid, effect allele, or beta)
    temp = strsplit(assoc$STRONGEST.SNP.RISK.ALLELE, "-")
    assoc$rsid = sapply(temp, function(l) l[1])
    assoc$effect_allele = toupper(sapply(temp, function(l) l[2]))
    assoc = assoc[!is.na(assoc$effect_allele) & !is.na(assoc$OR.or.BETA) & nchar(assoc$rsid) > 2 & substr(assoc$rsid, 1, 2) == "rs" & grepl("^[ACGT]$", assoc$effect_allele),]
    assoc = assoc[,colnames(assoc) != "STRONGEST.SNP.RISK.ALLELE"]

    assoc$id = paste(assoc$STUDY.ACCESSION, gsub(":", ";", assoc$MAPPED_TRAIT), gsub("^ *\\(", "", gsub(" *\\)$", "", assoc$P.VALUE..TEXT)), sep = ":")

    assoc = assoc[,c("id", "STUDY.ACCESSION", "PUBMEDID", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "pop", "pop.n.initial", "pop.n.replication", "rsid", "effect_allele", "log_beta", "PVALUE_MLOG")]
    colnames(assoc) = c(    "id", "accession",       "pubmed",   "trait",        "trait_uri",        "pop", "n.initial",     "n.replication",     "rsid", "effect_allele", "log_beta", "log10_p")

    assoc
}


rebase_associations = function(assoc, dbname = "hg19", host = "genome-mysql.cse.ucsc.edu", user = "genome", snp_table = "snp147Common")
{
    # Use the rsid and alternate allele information to rebase the scores on the
    # reference sequence.  This may require inversion of the betas, and associated
    # updating of an offset term for each polygenic score.
    # 
    # Note that this code uses a connection to a UCSC-style MySQL DB.
    message(sprintf("%d rsIDs in input", length(unique(assoc$rsid))))

    # Fetch dbSNP data for the rsids in assoc from a remote (eg UCSC) database.
    ucsc = src_mysql(dbname = dbname, host = host, user = user)
    ucsc.snps = tbl(ucsc, snp_table)
    temp.rsids = unique(assoc$rsid)
    ucsc.snps.matching_rsids = as.data.frame(select(filter(ucsc.snps, name %in% temp.rsids & locType == "exact"), 
        chrom, chromEnd, name, refNCBI, strand, observed, alleles, alleleFreqs))
    colnames(ucsc.snps.matching_rsids) = c("chrom", "pos", "rsid", "ref", "strand", "observed", "alleles", "alleleFreqs")
    message(sprintf("%d rsIDs found in database", length(unique(ucsc.snps.matching_rsids$rsid))))

    # Filter down to biallelic unambiguous strand SNPs.
    ucsc.snps.matching_rsids = ucsc.snps.matching_rsids[ucsc.snps.matching_rsids$observed %in% c("A/C", "A/G", "C/T", "G/T"),,drop=FALSE]
    message(sprintf("%d strand-specific SNPs", length(unique(ucsc.snps.matching_rsids$rsid))))

    # Filter down to canonical chromosomes
    ucsc.snps.matching_rsids$chrom = sub("^chr", "", ucsc.snps.matching_rsids$chrom)
    ucsc.snps.matching_rsids = ucsc.snps.matching_rsids[grepl("^([0-9]+|X|Y)$", ucsc.snps.matching_rsids$chrom),,drop=FALSE]
    message(sprintf("%d rsIDs in canonical chromosomes", length(unique(ucsc.snps.matching_rsids$rsid))))

    # Extract alleles and allele frequencies
    temp = t(simplify2array(strsplit(ucsc.snps.matching_rsids$alleles, ",")))
    colnames(temp) = c("allele1", "allele2")
    ucsc.snps.matching_rsids = cbind(ucsc.snps.matching_rsids, temp)
    ucsc.snps.matching_rsids$allele1 = as.character(ucsc.snps.matching_rsids$allele1)
    ucsc.snps.matching_rsids$allele2 = as.character(ucsc.snps.matching_rsids$allele2)
    temp = t(simplify2array(strsplit(ucsc.snps.matching_rsids$alleleFreqs, ",")))
    colnames(temp) = c("af1", "af2")
    ucsc.snps.matching_rsids = cbind(ucsc.snps.matching_rsids, temp)
    ucsc.snps.matching_rsids$af1 = as.numeric(as.character(ucsc.snps.matching_rsids$af1))
    ucsc.snps.matching_rsids$af2 = as.numeric(as.character(ucsc.snps.matching_rsids$af2))
    ucsc.snps.matching_rsids = ucsc.snps.matching_rsids[,!(colnames(ucsc.snps.matching_rsids) %in% c("alleles", "alleleFreqs"))]
    stopifnot(max(abs(ucsc.snps.matching_rsids$af1 + ucsc.snps.matching_rsids$af2 - 1)) < 2e-6)

    # Convert all alleles to the positive strand
    COMPLEMENT = c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
    ucsc.snps.matching_rsids$allele1posstrand = ifelse(ucsc.snps.matching_rsids$strand == "+", ucsc.snps.matching_rsids$allele1, COMPLEMENT[ucsc.snps.matching_rsids$allele1])
    ucsc.snps.matching_rsids$allele2posstrand = ifelse(ucsc.snps.matching_rsids$strand == "+", ucsc.snps.matching_rsids$allele2, COMPLEMENT[ucsc.snps.matching_rsids$allele2])
    stopifnot(ucsc.snps.matching_rsids$allele1posstrand == ucsc.snps.matching_rsids$ref | ucsc.snps.matching_rsids$allele2posstrand == ucsc.snps.matching_rsids$ref)

    # Identify which allele is the "alternate" (ie not genomic reference).
    # Store this allele in the alt field, and its frequency in the aaf field.
    ucsc.snps.matching_rsids$alt = ifelse(ucsc.snps.matching_rsids$allele2posstrand == ucsc.snps.matching_rsids$ref, ucsc.snps.matching_rsids$allele1posstrand, ucsc.snps.matching_rsids$allele2posstrand)
    ucsc.snps.matching_rsids$aaf = ifelse(ucsc.snps.matching_rsids$allele2posstrand == ucsc.snps.matching_rsids$ref, ucsc.snps.matching_rsids$af1, ucsc.snps.matching_rsids$af2)
    ucsc.snps.matching_rsids = ucsc.snps.matching_rsids[,c("rsid", "strand", "chrom", "pos", "ref", "alt", "aaf")]

    # Augment the association data with the dbSNP fields.  Note that 
    # some rows may be lost (all = FALSE) if they did not have a matching
    # rsid field in the filtered ucsc.snps.matching_rsids.
    assoc = merge(assoc, ucsc.snps.matching_rsids, by = "rsid", all = FALSE)
    message(sprintf("%d rsIDs following merge", length(unique(assoc$rsid))))

    # In most cases the effect allele is reported according to the SNP strand.  However, 
    # in about 20% of studies it's reported according to the positive (ie reference) strand,
    # regardless of the SNP orientation in dbSNP.  Therefore we really can't take much for
    # granted when interpreting the effect_allele field, and need to robustly identify which 
    # of the ref or alt alleles equals effect_allele, regardless of strand convention.
    # Fortunately we subset to strand-unambiguous biallelic SNPs in the above, so this 
    # conversion is possible.
    assoc$effect_allele_pos = ifelse(assoc$effect_allele == assoc$ref | assoc$effect_allele == assoc$alt, assoc$effect_allele, COMPLEMENT[assoc$effect_allele])

    # (Partially) verify that conversion worked
    stopifnot(assoc$effect_allele_pos == assoc$ref | assoc$effect_allele_pos == assoc$alt)

    # Now we can finally rebase the associations to the reference.
    assoc$logBeta = ifelse(assoc$effect_allele_pos == assoc$alt, assoc$log_beta, -assoc$log_beta)
    offsets = daply(assoc, .(id), function(d) sum(d$log_beta[d$effect_allele_pos == d$ref]))

    list(rebased = assoc[,c("id", "rsid", "chrom", "pos", "ref", "alt", "aaf", "logBeta")], offsets = offsets)
}


drop_small_studies = function(assoc, min.n = 5000)
{
    selected_ids = unique(assoc$id[assoc$n.initial >= min.n & assoc$n.replication >= min.n])

    assoc[assoc$id %in% selected_ids,,drop=FALSE]
}


drop_small_models = function(assoc, min.loci = 5)
{
    loci_per_id = daply(assoc, .(id), nrow)
    selected_ids = names(loci_per_id)[loci_per_id >= min.loci]

    assoc[assoc$id %in% selected_ids,,drop=FALSE]
}


estimate_score_variance = function(assoc)
{
    daply(assoc, .(id), function(a) sum(a$logBeta^2 * a$aaf*(1-a$aaf) * 2))
}


export_rebased_scores = function(assoc.rebased, output_path)
{
    offsets = assoc.rebased$offsets
    coefficients = assoc.rebased$rebased
    coefficients$vid = paste(coefficients$chrom, coefficients$pos, coefficients$ref, coefficients$alt, sep = ":")

    offset_df = data.frame(id = names(offsets), rsid = NA, chrom = NA, pos = NA, ref = NA, alt = NA, aaf = NA, logBeta = offsets, vid = "OFFSET")
    offset_df = offset_df[offset_df$id %in% unique(coefficients$id),]

    augmented_coefficients = rbind(coefficients, offset_df)
    augmented_coefficients = augmented_coefficients[order(augmented_coefficients$id, augmented_coefficients$chrom, augmented_coefficients$pos),]

    write.table(augmented_coefficients[,c("id", "vid", "rsid", "aaf", "logBeta")], output_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}



library(docopt)

sprintf('Convert EBI GWAS summary tables to POPSTAR PTS format.  Requries a connection to a UCSC-style MySQL DB.

Usage:
  ebi2pts.R [options] <associations> <ancestry>
  ebi2pts.R -h

Options: 
  -n <MIN_N>      Minimum GWAS sample and replication cohort size [default: 5000]
  -m <MIN_M>      Minimum number of loci in each PTS [default: 10]
  -b <POPN>       Target population regular expression [default: european]
  -o <OUT>        Path to output PTS file [default: /dev/stdout]
  -host <DBHOST>  MySQL genome database host [default: genome-mysql.cse.ucsc.edu]
  -db <DBNAME>    MySQL genome database name [default: hg19]
  -user <DBUSER>  MySQL genome database user [default: genome]
  -tbl <DBTABLE>  MySQL genome database table [default: snp147Common]
  -h              Print this help

Arguments: 
  associations  Associations TSV from EBI GWAS Catalog (gwas_catalog_v1.0.1-associations_*)
  ancestry      Ancestry TSV from EBI GWAS Catalog (gwas_catalog-ancestry_*)

%s
%s
', AUTHOR, VERSION) -> doc

opts = docopt(doc)

if (opts$h)
{
    print(doc)
    stop(0)
}

associations.raw = read_population_associations(opts$associations, opts$ancestry, opts$b)
message(sprintf("Loaded %d polygenic scores, %d loci total", length(unique(associations.raw$id)), length(unique(associations.raw$rsid))))
if (as.integer(opts$n) > 0)
{
    message("Filtering by cohort size... ", appendLF = FALSE)
    associations.raw = drop_small_studies(associations.raw, min.n = as.integer(opts$n))
    message(sprintf("%d polygenic scores, %d loci total remain", length(unique(associations.raw$id)), length(unique(associations.raw$rsid))))
}
if (as.integer(opts$m) > 1)
{
    message("Filtering by score size... ", appendLF = FALSE)
    associations.raw = drop_small_models(associations.raw, min.loci = as.integer(opts$m))
    message(sprintf("%d polygenic scores, %d loci total remain", length(unique(associations.raw$id)), length(unique(associations.raw$rsid))))
}
message("Rebasing associations... ")
associations.rebased = rebase_associations(associations.raw, dbname = opts$db, host = opts$host, user = opts$user, snp_table = opts$tbl)
if (as.integer(opts$m) > 1)
    associations.rebased$rebased = drop_small_models(associations.rebased$rebased, min.loci = as.integer(opts$m))
message(sprintf("Done.  Final count: %d polygenic scores, %d loci total", length(unique(associations.raw$id)), length(unique(associations.raw$rsid))))
export_rebased_scores(associations.rebased, opts$o)
