library("pacman")
pacman::p_load(
  "shiny", "shinydashboard", "ggplot2",
  "MetaboAnalystR", "KEGGREST", "tidyverse", "ggfortify",
  "httr", "RColorBrewer", "XML",
  "methods", "rvest", "data.table", "DT"
)


lookup_chem_id <- function(cpd_names_vec) {
  "Look up various IDs on KEGG, through MetaboAnalystR.
    Keeps only KEGG and HMDB IDs.
    Requirements: MetaboAnalystR
    
    PARAMS
    ------
    cpd_names_vec: vector of characters; common chemical names of compounds
    
    RETURNS
    -------
    cpd_names_tbl: tibble of IDs.
    "
  # Call to Kegg
  mset <- InitDataObjects("NA", "utils", FALSE)
  mset <- Setup.MapData(mset, cpd_names_vec)
  mset <- CrossReferencing(mset, "name", T, T, T, T, T)
  mset <- CreateMappingResultTable(mset)

  # print warnings
  print(mset$msgset$nmcheck.msg[2])
  cpd_names_tbl <- as_tibble(mset$dataSet$map.table)
  cpd_names_tbl <- cpd_names_tbl %>% rename("Sample" = "Query")
  cpd_names_tbl <- cpd_names_tbl %>% select("Sample", "KEGG", "HMDB")

  # replace NA, string "NA", or empty cell with string "undef"
  cpd_names_tbl[cpd_names_tbl == ""] <- "undef"
  cpd_names_tbl[cpd_names_tbl == "NA"] <- "undef"
  cpd_names_tbl <- cpd_names_tbl %>% replace(., is.na(.), "undef")

  return(cpd_names_tbl)
}


clean_and_transform_tibble <- function(tbl, log_bool, z_transform_bool) {
  "Origrinal clean_and_transform function.  
	"
  # replace NA with zero
  tbl <- tbl %>% replace(., is.na(.), 0)

  # Replace 0 with half of smallest nonzero value
  metab_names <- colnames(tbl)[3:length(colnames(tbl))]
  min.nonzero.val <- sort(unique(as.vector(as.matrix(tbl %>% select(metab_names)))))[2] / 2
  tbl[tbl == 0] <- min.nonzero.val

  # Take log (ugh), transforming in-place
  if (log_bool) {
    for (nm in metab_names) {
      tbl[nm] <- log2(unlist(tbl[nm]))
    }
  }
  # Do z-transform
  if (z_transform_bool) {
    for (nm in metab_names) {
      tmp <- unlist(tbl[nm])
      tbl[nm] <- (tmp - mean(tmp)) / sd(tmp)
    }
  }
  return(tbl)
}


get_results_tibble <- function(tbl, input_alpha, grp_numerator, grp_denominator, input_fdr) {
  "Takes a tibble of cleaned, transformed (if necessary) values as input, and does:
	1. t-tests at input_alpha
	2. BH correction at input_fdr
	3. Computes FC of averages
	4. Gets FC colours

    PARAMS
    ------
    tbl: input tibble of metabolite abundances. 
    input_alpha: float; alpha at which t-tests are applied.
    grp_numerator: numerator of class/group for fold change calculations
    grp_denominator: denominator of class/group for fold change calculations
    input_fdr: float; fdr at which the benjamini-hochberg method is applied.

    OUTPUT
    ------
    tbl: tibble with the following column names:
        'Sample': str; compound name
        'fc': float; fold change of log2-abundances. 
        'KEGG': str; KEGG Id, retrieved through MetaboAnalystR. 
        'HMDB': str; HMDB Id, retrieved through MetaboAnalystR.
        'raw_p_val': float; raw p-value, output through the t-test.
        'adj_p_val': float; adjusted p-value, computed through the Benjamini-Hochberg procedure.
        'fc_colour': str; colour of fold changes, in hex colour code. 

    NOTES
    -----
     * Unknown KEGG and HMDB Ids get replaced with string `undef`. 
	"
  group.name.col <- names(tbl[, 2])
  metab_names <- colnames(tbl)[3:length(colnames(tbl))]
  # Init vec of t.stats p-vals, and FCs
  # Compute t-stats
  p.vals.ls <- vector(mode = "numeric", length = length(metab_names))
  fc.ls <- vector(mode = "numeric", length = length(metab_names))
  for (i in 1:length(metab_names)) {
    g1.vec <- as.vector(unlist(tbl %>% filter(!!sym(group.name.col) == grp_numerator) %>% select(metab_names[i])))
    g2.vec <- as.vector(unlist(tbl %>% filter(!!sym(group.name.col) == grp_denominator) %>% select(metab_names[i])))

    x <- t.test(g1.vec, g2.vec, conf.level = 1 - input_alpha)
    p.vals.ls[i] <- x$p.value

    mu1 <- mean(g1.vec)
    mu2 <- mean(g2.vec)
    fc.ls[i] <- mu1 / mu2
  }


  names(p.vals.ls) <- metab_names
  names(fc.ls) <- metab_names

  # adjust: BH correction
  adj.p.vals.ls <- p.adjust(p.vals.ls, method = "hochberg", n = length(p.vals.ls))

  # Compute ipath colours: non-significant, FC-up, FC-down, FC-neutral
  fc.colour.ls <- rep("#ACACAC", ncol(tbl)) # default gray (non-significant)
  names(fc.colour.ls) <- metab_names
  for (nm in metab_names) {
    if (adj.p.vals.ls[nm] < input_fdr) {
      if (fc.ls[nm] > 1.2) {
        fc.colour.ls[nm] <- "#0571b0" # blue
      } else if (fc.ls[nm] < 0.8) {
        fc.colour.ls[nm] <- "#ca0020" # red
      } else {
        fc.colour.ls[nm] <- "#FFFFFF" # black
      }
    }
  }

  # Get all KEGG IDs
  chem.id.tbl <- lookup_chem_id(metab_names)
  kegg.id.vec <- as.vector(unlist(chem.id.tbl %>% select("KEGG")))

  # enframe and merge all named lists
  fc.tbl <- tibble::enframe(fc.ls) %>% rename("fc" = value, "Sample" = name)
  p.val.tbl <- tibble::enframe(p.vals.ls) %>% rename("raw_p_val" = value, "Sample" = name)
  adj.p.val.tbl <- tibble::enframe(adj.p.vals.ls) %>% rename("adj_p_val" = value, "Sample" = name)
  fc.colour.tbl <- tibble::enframe(fc.colour.ls) %>% rename("fc_colour" = value, "Sample" = name)
  tbl <- list(fc.tbl, chem.id.tbl, p.val.tbl, adj.p.val.tbl, fc.colour.tbl) %>% reduce(inner_join, by = "Sample")

  return(tbl)
}


get_ipath_selection_str <- function(tbl, metab_name_colname, kegg_colname, fc_colour_colname, node_width_text) {
  "Get the ipath selection string from an input tibble. Usually goes after get_results_tibble() after the appropriate 
    filter()-ing or select()-ing. 

    PARAMS
    ------
    tbl: input tibble. 
    metab_name_colname: str; name of column containing metabolite names. 
    kegg_colname: str;name of column containing KEGG IDs
    fc_colour_colname: str; name of column containing colours in hex colour code format. 
    node_width_text: str; width of node, for ipath entry. Usually 'W20'

    RETURNS
    -------
    selection_str
    "

  # Extract columns as named lists, named by metab names
  kegg.vec <- as.vector(unlist(tbl[kegg_colname]))
  names(kegg.vec) <- as.vector(unlist(tbl[metab_name_colname]))
  fc.colour.ls <- as.vector(unlist(tbl[fc_colour_colname]))
  names(fc.colour.ls) <- as.vector(unlist(tbl[metab_name_colname]))

  selection_str_ls <- c()
  for (nm in names(kegg.vec)) {
    selection_str_ls <- c(selection_str_ls, paste(kegg.vec[nm], node_width_text, fc.colour.ls[nm], sep = " "))
  }

  selection_str <- paste0(selection_str_ls, collapse = "\n")
  return(selection_str)
}


call_maca_pw_analysis <- function(fn_auc_csv, kegg_species_id) {
  "Calls the pathway enrichment analysis module from MetaboAnalystR. 
    Does row-wise median-normalization and log-transforms the data.
    P-values of pathway enrichment are calculated using the `globaltest` algorithm, and pathway impact scores
    computed using the pathway centrality option. But impact should be disregarded as an overly-abstract
    graph theoretic notion that doesn't necessarily have any biological relevance. 

    PARAMS
    ------
    fn_auc_csv: str; filename of input run summary table as a csv file, with AUCs as values. 
    rownames are the sample names, column names are the metabolite names. Column 1 are the 
    experimental groupings. Because of the way this module works, only 2 groups are supported.

    RETURNS
    -------
    list of two outputs:
    tbl.out: output tibble of pathway analysis enrichment. columns:
        metabolite (compound common name), total cmpd, Hits, raw p (raw p value), -log p, 
        Holm adjust(ed p value), FDR, Impact. 
    pw.dict: named list of lists; each key is the pathway ID, and each value is a list of 
        compounds from the input data which appear in that particular pathway. 
    "

  mSet <- InitDataObjects("conc", "pathqea", FALSE)
  mSet <- Read.TextData(mSet, fn_auc_csv, "rowu", "disc")
  mSet <- CrossReferencing(mSet, "name")
  mSet <- CreateMappingResultTable(mSet)
  mSet <- SanityCheckData(mSet)
  mSet <- ReplaceMin(mSet)
  mSet <- PreparePrenormData(mSet)
  mSet <- Normalization(mSet, "MedianNorm", "LogNorm", "NULL", ratio = FALSE, ratioNum = 20)
  mSet <- SetKEGG.PathLib(mSet, kegg_species_id)
  mSet <- SetMetabolomeFilter(mSet, F)
  mSet <- CalculateQeaScore(mSet, "rbc", "gt")

  tbl.out <- as_tibble(mSet$analSet$qea.mat, rownames = "pw_name")
  pw.dict <- mSet$analSet$qea.hits

  return(list(tbl.out, pw.dict))
}

# This used to be in output$ipath
# Due for deletion soon

'
tbl0 <- read_csv(input$file1$datapath)
        
        group.name.col <- names(tbl0[,2])
        groupnames <- as.vector(unlist(unique(tbl0[[group.name.col]])))
        control.group.name <- "Plus" # name of denominator group. 
        alpha <- input$fdr_num
        kegg_species_id <- input$species_code
        
        #  ==================== START ==================== 
        cnames <- colnames(tbl0)[3:length(colnames(tbl0))]
        # Init vec of t.stats p-vals, and fold change of averages
        # Compute t-stats
        p.vals.ls <- vector(mode="numeric", length = length(cnames))
        fc.ls <- vector(mode="numeric", length = length(cnames))
        for (i in 1:length(cnames)) {
            g1.vec <- as.vector(unlist(tbl0 %>% filter(!!sym(group.name.col)==groupnames[1]) %>% select(cnames[i])))
            g2.vec <- as.vector(unlist(tbl0 %>% filter(!!sym(group.name.col)==groupnames[2]) %>% select(cnames[i])))
            x <- t.test(g1.vec, g2.vec)
            p.vals.ls[i] <- x$p.value
            
            mu1 <- mean(g1.vec)
            mu2 <- mean(g2.vec)
            fc.ls[i] <- mu1/mu2
        }
        names(p.vals.ls) <- colnames(tbl0)[3:length(colnames(tbl0))]
        names(fc.ls) <- colnames(tbl0)[3:length(colnames(tbl0))]
        
        # adjust: BH correction
        # pick out names where p < alpha
        p.vals.ls <- p.adjust(p.vals.ls, method = "hochberg", n = length(p.vals.ls))
        significant.metabs <- c()
        for (nm in names(p.vals.ls)) {
            if (p.vals.ls[[nm]] < alpha) {
                significant.metabs <- c(significant.metabs, nm)
            }
        }
        print(paste0(length(significant.metabs), " significant metabolites found at alpha = ", alpha))
        
        # Filter for significant metabs
        p.vals.ls2 <- p.vals.ls[significant.metabs]
        fc.ls2 <- fc.ls[significant.metabs]
        
        # Colour by FC
        th.lower <- 0.8
        th.upper <- 1.2
        fc.colour.ls <- rep("#000000", length(significant.metabs))
        names(fc.colour.ls) <- significant.metabs
        for (i in 1:length(significant.metabs)) {
            nm <- significant.metabs[i]
            if (fc.ls[nm] <= th.lower) {
                fc.colour.ls[i] <- "#d95f02" #orange
            }
            if (fc.ls[nm] > th.upper) {
                fc.colour.ls[i] <- "#1b9e77" #green
            }
        }
        print("check")
        fc.tbl <- tibble::enframe(fc.ls2) %>% rename("fc"=value, "Sample"=name)
        fc.colour.tbl <- tibble::enframe(fc.colour.ls) %>% rename("fc_colour"=value, "Sample"=name)
        fc.tbl <- inner_join(fc.tbl, fc.colour.tbl, by="Sample")
        p.val.tbl <- tibble::enframe(p.vals.ls2) %>% rename("FDR"=value, "Sample"=name)
        
        #  ==================== GET CPD IDS ==================== 
        kegg.id.tbl <- lookup_chem_id(names(fc.ls2))
        kegg.id.vec <- as.vector(unlist(kegg.id.tbl %>% select("KEGG")))
        
        tbl1 <- inner_join(fc.tbl, kegg.id.tbl, by="Sample")
        tbl1 <- inner_join(tbl1, p.val.tbl, by="Sample") %>% select(Sample, fc, fc_colour, KEGG, FDR)
    '
