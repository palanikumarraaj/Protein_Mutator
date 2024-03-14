library(Biostrings)
args <- commandArgs(trailingOnly = TRUE)
user_ID <- args[1]
file_path <- paste0("Protein_Normal/", user_ID, ".fasta")

# Step 2: Load mutations from user-provided ID.tsv file
#mutation_file <- user_ID
mutation_file <- paste0(user_ID, ".csv")
mutations <- read.csv(mutation_file, header = TRUE, stringsAsFactors = FALSE)
#mutations$replacement <- as.character(mutations$replacement)
#mutations$V3[mutations$V3 == "TRUE"] <- "T"
#head(mutations)

data <- mutations
# Function to convert TRUE/FALSE to T/F

convert_logical <- function(x) {
  ifelse(x == "TRUE", "T", ifelse(x == "FALSE", "F", x))
}

for (col in names(data)) {
  data[, col] <- sapply(data[, col], convert_logical)
}
#head(data)
mutations <- data
# Step 3: Convert loaded data to list format
mutations <- lapply(1:nrow(mutations), function(i) {
  list(position = mutations[i, 1], original = mutations[i, 2], replacement = mutations[i, 3])
})

# Load the FASTA file
fasta_seq <- readAAStringSet(file_path)
##############
# Function to create a modified FASTA file
create_modified_fasta <- function(fasta_seq, mutation, output_file_path) {
  modified_seq <- replaceAt(fasta_seq, mutation$position, mutation$replacement)
  header <- paste0(names(fasta_seq), "_", mutation$original, mutation$position, mutation$replacement)
  writeXStringSet(modified_seq, output_file_path, format = "fasta")
}
# Create a file with all 3 mutations
all_mutations_fasta <- fasta_seq
for (mutation in mutations) {
  all_mutations_fasta <- replaceAt(all_mutations_fasta, mutation$position, mutation$replacement)
}
# Extract the FASTA prefix from the file path
fasta_prefix <- tools::file_path_sans_ext(basename(file_path))

# Function to create a modified FASTA filename
create_modified_filename <- function(mutation) {
  return(paste0(fasta_prefix, "_", mutation$original, mutation$position, mutation$replacement, ".fasta"))
}

###
# Create the Protein_Mutant folder if it doesn't exist
if (!dir.exists("Protein_Mutant")) {
  dir.create("Protein_Mutant")
}

for (mutation in mutations) {
  # Create modified filename with path
  output_filename <- file.path("Protein_Mutant", create_modified_filename(mutation))
  
  # Create the modified FASTA file
  create_modified_fasta(fasta_seq, mutation, output_filename)
}

# Combine all mutations into a single string
all_mutations_str <- paste(sapply(mutations, function(m) paste0(m$original, m$position, m$replacement)), collapse = "_")

# Create filename for combined mutations with path
combined_filename <- file.path("Protein_Mutant", paste0(fasta_prefix, "_", all_mutations_str, ".fasta"))

# Create the file with all mutationsstr
create_modified_fasta(all_mutations_fasta, list(), combined_filename)
