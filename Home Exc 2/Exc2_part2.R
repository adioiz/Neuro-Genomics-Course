#################
#### PART 2 #####
#################

# Load the Circadian RNA-seq data into a matrix -> data:
Matrix_of_part_5 <- read.csv("/Users/adioiz/Documents/Learnings/Neuro-Genomics Course/Home Exc 1/CircadianRNAseq.csv")
data <- as.matrix(Matrix_of_part_5)

## Set constants
circadian_frequency <- 1 / 24
sampling_interval <- 4  # 4-hour sampling interval
N <- length(data[1, 2:13])

# Function to calculate G-factors
calculate_G_factors <- function(data) {
  circadian_powers <- numeric(nrow(data))
  
  for (i in 1:nrow(data)) {
    expression_levels <- as.numeric(data[i, 2:13])
    
    # Skip genes with all zero expression levels
    if (all(expression_levels == 0, na.rm = TRUE)) {
      circadian_powers[i] <- 0
      next
    }
    
    # Perform FFT and compute the power spectrum
    circadian_fft <- fft(expression_levels)
    circadian_power <- Mod(circadian_fft)^2
    unique_circadian_powers <- circadian_power[2:7]
    
    # Safeguard against NAs and zero sums
    sum_circadian_powers <- sum(unique_circadian_powers, na.rm = TRUE)
    if (is.na(sum_circadian_powers) || sum_circadian_powers == 0) {
      circadian_powers[i] <- 0
      next
    }
    
    # Normalize and find the circadian frequency index
    circadian_power_normalized <- unique_circadian_powers / sum_circadian_powers
    frequencies <- seq(1 / (N * sampling_interval), 1 / (2 * sampling_interval), 
                       length.out = length(circadian_power_normalized))
    circadian_index <- which.min(abs(frequencies - circadian_frequency))
    circadian_powers[i] <- circadian_power_normalized[circadian_index]
  }
  
  return(sort(circadian_powers, decreasing = TRUE, na.last = TRUE))
}

## 1. Calculate G-factors for the original data
circadian_powers <- calculate_G_factors(data)

# Count G factors for specific cutoffs
vector_for_G_factor_cutoffs_amount <- sapply(0:100, function(i) {
  cutoff <- i / 100
  sum(circadian_powers >= cutoff, na.rm = TRUE)
})
temp_for_backup <- vector_for_G_factor_cutoffs_amount

## 2. Shuffle the columns to remove biological temporal order
averaged_vector_result <- numeric(101)

set.seed(42)  # For reproducibility
for (averaging_index in 1:100) {
  # Shuffle columns (samples)
  shuffled_data <- data
  shuffled_data[, 2:13] <- shuffled_data[, sample(2:13)]
  
  # Calculate G-factors for shuffled data
  shuffled_powers <- calculate_G_factors(shuffled_data)
  
  # Count G factors for specific cutoffs
  shuffled_cutoffs <- sapply(0:100, function(i) {
    cutoff <- i / 100
    sum(shuffled_powers >= cutoff, na.rm = TRUE)
  })
  
  averaged_vector_result <- averaged_vector_result + shuffled_cutoffs
}

# Average the shuffled results over 100 iterations
averaged_vector_result <- averaged_vector_result / 100

## 3. Plot the results
relevant_x_values <- seq(0, 1, by = 0.01)

# Plot 1: Original vs Shuffled G factor counts
plot(relevant_x_values, temp_for_backup, type = "o", col = "blue", pch = 16, lwd = 2,
     xlab = "G factor used as a cutoff", ylab = "The vector values", 
     main = "The result for cutoff of shuffled data")
lines(relevant_x_values, averaged_vector_result, col = "green", pch = 16, lwd = 2)
legend("topright", legend = c("Original", "Shuffled"), col = c("blue", "green"), pch = 16, lwd = 2)

# Plot 2: True Positive Fractions
True_positive_fractions <- (temp_for_backup - averaged_vector_result) / temp_for_backup



plot(relevant_x_values, True_positive_fractions, type = "o", col = "blue", pch = 16, lwd = 2,
     xlab = "G factor used as a cutoff", ylab = "True positive fraction",
     main = "True positive fractions by the G factor cutoffs")

## 4. Find G-factor threshold for True Positive Rate of 0.8
threshold_index <- which.min(abs(True_positive_fractions - 0.8))  # Closest to 0.8
G_factor_threshold <- relevant_x_values[threshold_index]
cat("G-factor cutoff corresponding to 0.8 True Positive Rate:", G_factor_threshold, "\n")

# Extract the genes detected as circadian above this threshold
genes_names_above_threshold <- data[which(circadian_powers >= G_factor_threshold), ncol(data)]
genes_above_threshold <- data[which(circadian_powers >= G_factor_threshold), 2:ncol(data)]
cat("Number of circadian genes detected:", length(genes_above_threshold), "\n")

# Print the gene names
cat("List of detected circadian genes:\n")
print(genes_names_above_threshold)


write.csv(genes_above_threshold, "/Users/adioiz/Documents/Learnings/Neuro-Genomics Course/Home Exc 2/genes_above_threshold.csv", row.names=FALSE)




