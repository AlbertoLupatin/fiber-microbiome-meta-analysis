#' Merge a list of metadata into a new dataframe. It also mutate Study and Invididual columns as factor
#'
#' @param Metadata_list a list of metadata with unique individuals
#'
#' @return A dataframe with all the metadata in the list

merge_metadata <- function(Metadata_list) {

  metadata_merged <- purrr::reduce(.x = Metadata_list,
                                   merge, all = T) %>%
    column_to_rownames("Run") %>%
    mutate(Study = as.factor(Study),
           Individual = as.factor(Individual))

  return(metadata_merged)
}


#' Retrieve Specie Name from Full Taxonomy
#'
#' @description
#' The function retrieves the specie name and eventually its identificative number.
#' It substitutes `Specie incertae sedis` with `sp.`
#'
#' @param otu_table an OTU table having the full taxonomy string in the columns
#'
#' @return an OTU table having the specie name and its identificative number as column names

clearNames <- function(specie_names) {
  specie_names <- gsub(".*s__", "s__", specie_names)
  specie_names <- gsub("(?<!s__)\\[.*?\\]", "", specie_names, perl = TRUE)
  specie_names <- gsub("\\/.*\\|", "|", specie_names, perl = TRUE)
  specie_names <- gsub("s__", "", specie_names)
  specie_names <- gsub("\\|ext_mOTU_v31_", "_e", specie_names)
  specie_names <- gsub("\\|ref_mOTU_v31_", "_r", specie_names)
  specie_names <- gsub("\\|meta_mOTU_v31_", "_r", specie_names)
  specie_names <- gsub("Specie incertae sedis", "sp.", specie_names)
  specie_names <- gsub("species incertae sedis", "sp.", specie_names)
  specie_names <- gsub("\\b(\\w+)\\s+\\1\\b", "\\1", specie_names, perl = TRUE)
  specie_names <- gsub(" CAG:\\d+", "", specie_names, perl = TRUE)
  specie_names <- gsub("_", " ", specie_names, perl = TRUE)
  specie_names <- gsub("\\[|\\]", "", specie_names)
  specie_names <- str_replace(specie_names, "\\s{2,}", " ")
  specie_names <- str_replace(specie_names, "^\\s+", "")
  specie_names <- str_replace(specie_names, "\\b(\\w+)\\s+\\1\\b", "\\1")

  return(specie_names)
}


#' Save a list of plots to a pdf
#'
#' @description
#' Save a list of plots to a pdf
#'
#' @param plots A list of plots
#' @param filename String with the complete path for where do you want the file to be saved
#' @param width plot width
#' @param height plot height
#'
#' @return A pdf file with one plot for each page

savePlots <- function(plots, filename, width, height) {

  if (missing(width) & missing(height)) {pdf(filename)}
  else if (missing(width) & !(missing(height))) {pdf(filename, height = height)}
  else if (!(missing(width)) & missing(height)) {pdf(filename, width = width)}

  for (i in 1:length(plots)){
    print(plots[[i]])
  }
  dev.off()

}

#' Run SIAMcat Linar Mixed Model
#'
#' @description
#' From the phyloseq and metadata inputs, runs SIAMcat's lmm through the before and after timepoints in the metadata
#'
#' @param phyloseq Phyloseq object
#' @param metadata Dataframe with the timepoints before and after and, eventually, study variable
#' @param study_covariate include the study variable as a random effect in the linear model
#'
#' @return A siamcat object with the associations stored in it

runSIAMCAT <- function(phyloseq, metadata, study_coveriate = TRUE) {

  label <- create.label(meta = metadata, label = "Timepoint", case = "After", verbose = 0)

  sc <- siamcat(feat = phyloseq@otu_table, meta = metadata, label = label, verbose = 0)

  sc <- filter.features(sc, filter.method = "prevalence", cutoff = 0.1, verbose = 0)
  sc <- normalize.features(sc, norm.method = "log.std")

  if (study_coveriate & length(unique(metadata$Study)) > 1) {
    sc <- check.associations(sc, formula = "feat~label+(1|Study)+(1|Individual)", verbose = 0)
  } else {
    sc <- check.associations(sc, formula = "feat~label+(1|Individual)", verbose = 0)
  }

  return(sc)
}

