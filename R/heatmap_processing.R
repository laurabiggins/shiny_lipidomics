library(pheatmap)
library(magrittr)

#==============
# functions

# this assumes that the first column contains names
remove_0_variance <- function(dataset, name_column = 1){
  variances <- apply(dataset[, -1], 1, var)
  dataset[variances > 0, ]
}

# determine which group the sample name is
get_group <- function(sample_name, groups) {
  matches <- sapply(groups, grepl, sample_name, ignore.case = TRUE)
  names(matches)[matches]
}

#' create_paired_seq
#'
#' output can be used for creating a constrasting colour palette
#' create_paired_seq() would return a vector containing
#' 1  6  2  7  3  8  4  9  5 10
#'
#' @param seq_length sequence length, integer value
#'
#' @return integer vector
#' @export
#'
#' @examples
#' create_paired_seq()
create_paired_seq <- function(seq_length = 10){
  
  if(!seq_length %% 2 == 0) {
    warning("seq_length in paired seq function should be an even number, rounding up.")
  }  
  half_seq <- ceiling(seq_length/2)
  two_no <- function(x, half_seq) c(x, x + half_seq)
  as.vector(sapply(1:half_seq, two_no, half_seq))
} 

#===========================================
# sorting out lipid classes

parse_lipid_classes <- function(raw_dataset, lipid_column){
  
  lipid_classes <- raw_dataset %>%
    dplyr::select(tidyselect::all_of(lipid_column)) %>%
    tidyr::separate(
      col = tidyselect::all_of(lipid_column), 
      into = c("class", "info"), 
      sep = " "
    ) %>%
    dplyr::pull(class)
  
  tibble::tibble(
    class = forcats::fct_inorder(lipid_classes),
    names = dplyr::pull(raw_dataset, lipid_column)
  )
}

#' get_lipid_summary
#'
#'summarise lipid information
#'
#' @param all_lipids 
#' @param lipid_class 
#'
#' @return 2 column tibble, first containing the lipid classes, the 2nd containing 
#' the counts. Column names "class" and "count"
#' @export
#'
#' @examples
get_lipid_summary <- function(all_lipids, lipid_class = "class"){
  
  lipid_tibble <- all_lipids %>%
    dplyr::group_by(.data[[lipid_class]]) %>%
    dplyr::count() 
  
  colnames(lipid_tibble) <- c("class", "count")
  
  lipid_tibble
}  

# We only want row annotations once per lipid group, in the middle
# of the row, so we need to pass a vector with a load of blanks among the class names. 
#' Title
#'
#' @param lipid_summary input should be a lipid summary object that comes from 
#' the function get_lipid_summary
#'
#' @return
#' @export
#'
#' @examples
create_lipid_class_labels <- function(lipid_summary, total_no_lipids){
  
  ordered_counts <- dplyr::pull(lipid_summary, count)
  ordered_classes <- as.vector(dplyr::pull(lipid_summary, class))
  
  cumulative_counts <- head(cumsum(ordered_counts), -1) # remove final one
  pos <- c(0, cumulative_counts) + ceiling(ordered_counts/2)
  
  lipid_labels <- rep("", total_no_lipids)
  lipid_labels[pos] <- ordered_classes
  lipid_labels
}

#' create_lipid_colours
#' 
#' for creating a colour for each lipid class
#'
#' @param class_names 
#' @param no_colours number of different colours/shades
#' @param colours  2 colours to create a palette from
#'
#' @return named vector of length class names. The names are the class names. The 
#' colours will repeat 
#' @export
#'
#' @examples
#' create_lipid_colours(letters[1:20])
create_lipid_colours <- function(class_names, no_colours = 10, colours = c("lightgrey", "navy")){
  
  my_seq <- rep_len(create_paired_seq(no_colours), length.out = length(class_names))
  stats::setNames(
    object = colorRampPalette(colours)(no_colours)[my_seq], 
    class_names
  )
}

