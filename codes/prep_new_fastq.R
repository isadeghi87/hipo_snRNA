setwd(dir = '/omics/odcf/analysis/hipo2/hipo_K35R')
setwd('~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/hipo_k35/')

library(stringr)

## read sample sheet 

dat = readxl::read_excel(path = "./datasets/KATI-HIPO-K35 .xls",skip = 5,sheet = 2)
dat = subset(dat,Rerun == 'Yes')
unique(dat$Species)

# List of IDs to search for
# human <- c("K35R-7P6VFX", "K35R-4DZ4CT", "K35R-VFHB3C", "K35R-BNTVHX", "K35R-BZ8JPB")
# patterns <- "^([^-]+-[^-]+)"
# dat$id =str_extract(dat$`HIPO ID`, patterns)

# Base directory containing the original files
base_dir <- '/omics/odcf/project/hipo2/hipo_K35R/sequencing/10x_scRNA_sequencing/view-by-pid/'

# Directory to store the symbolic links

# Find directories containing the sequences
dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
sequence_dirs <- dirs[grep('sequence', dirs)]

## iterate over species


# Function to create symbolic links with the new naming convention
create_symlink_with_rename <- function(original_dir, new_dir_base) {
  files <- list.files(path = original_dir, pattern = "fastq.gz$", full.names = TRUE)
  
  # Extract sample ID from the directory path
  sample_id <- str_extract(original_dir, paste(ids, collapse = "|"))
  
  # Extract different parts of the filename (e.g., sample, lane, read type)
  for (i in seq_along(files)) {
    original_path <- files[i]
    
    # Extract read type, lane, and sample number
    file_name <- basename(original_path)
    read_type <- ifelse(grepl("_R1", file_name), "R1", "R2")
    
    # Modify lane and sample number based on the filename structure
    lane_number <- str_extract(file_name, "L00[1-9]")
    sample_number <- str_extract(file_name, "S[0-9]+")
    
    # Extract unique file prefix (for cases like AS-723277 etc.)
    file_prefix <- str_extract(file_name, "^[^-]+-[^-]+")
    
    # Construct the new filename
    new_file_name <- paste0(file_prefix, "_", sample_number, "_", lane_number, "_", read_type, "_001.fastq.gz")
    
    # Full path for the new symbolic link
    full_new_path <- file.path(new_dir_base, file_prefix, new_file_name)
    
    # Create the directory for the sample if it doesn't exist
    dir.create(dirname(full_new_path), recursive = TRUE, showWarnings = FALSE)
    
    # Create the symbolic link
    command <- paste("ln -s", shQuote(original_path), shQuote(full_new_path))
    system(command)
  }
}

species <- na.omit(unique(dat$Species))
for (sp in species) {
  data <- dat[dat$Species == sp,]
  ids <- unique(na.omit(data$`HIPO ID`))
  print(sp)
  print(ids)
  
  # Link base directory for species
  link_base_dir <- paste0("/omics/odcf/analysis/hipo2/hipo_K35R/data/linked_fastqs/", sp)
  
  # Filter directories by matching ids
  filtered_dirs <- sequence_dirs[grep(paste(ids, collapse = "|"), sequence_dirs)]
  
  # Process each filtered directory
  for (dir_path in filtered_dirs) {
    create_symlink_with_rename(dir_path, link_base_dir)
  }
}