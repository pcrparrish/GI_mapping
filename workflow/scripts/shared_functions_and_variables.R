## FOR USE WITH GI MAPPING PIPELINE


## ggplot themes
## see: https://www.rdocumentation.org/packages/ggplot2/versions/2.1.0/topics/theme_update
## and https://stackoverflow.com/questions/23173915/can-ggplot-theme-formatting-be-saved-as-an-object
plot_theme <- theme(axis.text = element_text(colour="black"),
                    axis.ticks = element_line(color="black"))

plot_options <- list(theme_bw(base_size = 12))


## set default aspect ratios
wide_ar <- 0.75
square_ar <- 1

print_kbl <- function(tbl) {
  kbl(tbl) %>%
    kable_styling(full_width = FALSE, 
                  position = "left",
                  bootstrap_options = c("striped", "hover", "responsive"))
}

save_tbl <- function(tbl){
  tbl_str <- deparse(substitute(tbl))
  tbl_name <- str_split(tbl_str, pattern = "\\.")[[1]][2]
  write_tsv(tbl, file.path(out_dir, "tables", "tsv", paste0(params$cell_line, "_", tbl_name, ".txt")))
  write_rds(tbl, file.path(out_dir, "tables", "rds", paste0("d.", params$cell_line, "_", tbl_name)))
}


save_plot <- function(plt){
  plt_str <- deparse(substitute(plt))
  ggsave(plot = plt, 
         filename = file.path(out_dir, "plots", "pdf", paste0(params$cell_line, "_", plt_str, ".pdf")))
  ggsave(plot = plt, 
         filename = file.path(out_dir, "plots", "png", paste0(params$cell_line, "_", plt_str, ".png")))
}

make_out_dir <- function(out_dir){
  if(dir.exists(out_dir)){
    ## print a message with the output directory location
    print(paste("Output directory already exists at:", out_dir, sep = " "))
  } else{
    ## make output dirs
    dir.create(file.path(out_dir, "tables", "rds"), recursive = TRUE)
    dir.create(file.path(out_dir, "tables", "tsv"), recursive = TRUE)
    dir.create(file.path(out_dir, "plots", "png"), recursive = TRUE)
    dir.create(file.path(out_dir, "plots", "pdf"), recursive = TRUE)
    print(paste("Output directory created at:", out_dir, sep = " "))
  }
}