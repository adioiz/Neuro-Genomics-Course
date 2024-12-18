##################
##### Part 2 #####
##################

if (!require("BiocManager", quietly = TRUE))
  {install.packages("BiocManager")
    BiocManager::install("pasilla")
    }
library("pasilla")

## Load the actual data from the pasilla package
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
rownames(coldata) <- sub("fb", "", rownames(coldata))
cts <- cts[, rownames(coldata)]

## Examine the first 10 lines of this matrix - what kind of information is in the matrix?

# Display the first 10 rows of the 'cts' matrix
first10lines <- head(cts, 10)
print(first10lines)

## Define a variable to hold the dimensions of the cts matrix
cts_dimensions <- dim(cts)
print(cts_dimensions)

## Calculate the sum of reads for each sample
sample_sums <- colSums(cts)

# Check if all values in the vector are the same
all_equal <- all(sample_sums == sample_sums[1])

# Print the result
print(all_equal)

## Create a normalized version of the cts matrix
# Set the target read count to the total of the first column
target_sum <- sample_sums[1]

# Calculate scaling factors for each column
scaling_factors <- target_sum / sample_sums

# Create a normalized version of the cts matrix
cts_normalized <- cts

# Apply scaling factors to each column except the first one
cts_normalized[, -1] <- sweep(cts[, -1], 2, scaling_factors[-1], "*")

# Print the normalized matrix
print(cts_normalized)

## Make sure that in the normalized matrix the sum of reads is the same in all samples
normalized_sums <- colSums(cts_normalized)
# Round the normalized sums to the nearest integer
normalized_sums <- round(normalized_sums)
# Verify if all column sums are equal to the target sum

if (all(normalized_sums == normalized_sums[1])) {
  print("All columns have the same total read count.")
} else {
  print("The columns do not have the same total read count.")
}
##################
##### Part 3 #####
##################

## Does the data fit a Poisson distribution?
# Select columns for untreated & treated conditions
cts_untreated <- cts_normalized[, grep("untreated", colnames(cts_normalized))]
cts_treated <- cts_normalized[, grep("treated", colnames(cts_normalized))]

# Calculate the mean expression for each gene across samples
mean_expression <- rowMeans(cts_untreated)
mean_expression_for_treated <- rowMeans(cts_treated)
# Calculate the variance of expression for each gene across samples
variance_expression <- apply(cts_untreated, 1, var)
variance_expression_for_treated <- apply(cts_treated, 1, var)

# Add 1 to avoid log of 0, then log-transform the mean and variance
log_mean_expression <- log(mean_expression + 1)
log_variance_expression <- log(variance_expression + 1)
log_mean_expression_for_treated <- log(mean_expression_for_treated + 1)
log_variance_expression_for_treated <- log(variance_expression_for_treated + 1)
# Plot the log-transformed variance vs. log-transformed mean
plot(log_mean_expression, log_variance_expression, 
     xlab = "Log Mean Expression", 
     ylab = "Log Variance Expression", 
     main = "Variance vs Mean for Untreated Conditions",
     pch = 16, col = "blue")

# Add the line y = x to the plot for reference
abline(0, 1, col = "red", lty = 2)  # y = x line (slope = 1, intercept = 0)

## Does the data fit a Negative Binomial distribution?

# Fit the nonlinear model variance = mean + a * mean^2:
curved_fit <- nls(log_variance_expression ~ log_mean_expression + a * (log_mean_expression^2), 
           start = list(a = 0))
curved_fit_for_treated <- nls(log_variance_expression_for_treated ~ log_mean_expression_for_treated + a * (log_mean_expression_for_treated^2), 
                  start = list(a = 0))
# Extract the fitted coefficient 'a'
a <- coef(curved_fit)
a_treated <- coef(curved_fit_for_treated)
# Generate points for the fitted curve
curve_points <- log_mean_expression + a * (log_mean_expression^2)
curve__points_treated <- log_mean_expression_for_treated + a_treated * (log_mean_expression_for_treated^2)
# Plot the log-transformed mean vs log-transformed variance untreated
plot(log_mean_expression, log_variance_expression, 
     xlab = "Log Mean Expression", 
     ylab = "Log Variance Expression", 
     main = "Variance vs Mean with Fitted Negative Binomial Curve - untreated",
     pch = 16, col = "blue")

# Add the fitted curve to the plot
points(log_mean_expression, curve_points, col = "red", lwd = 2, pch=16)

# Plot the log-transformed mean vs log-transformed variance treated
plot(log_mean_expression_for_treated, log_variance_expression_for_treated, 
     xlab = "Log Mean Expression", 
     ylab = "Log Variance Expression", 
     main = "Variance vs Mean with Fitted Negative Binomial Curve - treated",
     pch = 16, col = "blue")

# Add the fitted curve to the plot
points(log_mean_expression_for_treated, curve__points_treated, col = "red", lwd = 2, pch=16)

##################
##### Part 4 #####
##################

## Detect at least one gene that has different expression levels in the first sample

log_first_treated_sample <- log(cts_normalized[, "treated1"] + 1)
log_first_untreated_sample <- log(cts_normalized[, "untreated1"] + 1)
# Plot to compare the expression levels
plot(log_first_untreated_sample,
     xlab = "Index", 
     ylab = "Log Expression", 
     main="Genes expression levels - treated v.s. untreated, in log scale", 
     col = 'red')
points(log_first_treated_sample, col= 'blue')
# a great way to separate some of the genes here is by using the following line that is set at e^10, for visualization:
visualization_assisting_line <- vector("numeric", length = 15000)
visualization_assisting_line <- visualization_assisting_line+11 # the height can change, if it is needed after normalizing cts.
lines(visualization_assisting_line, col ='green', lwd=2)
legend("topright", legend = c("untreated1", "treated1"), fill = c("red", "blue"), border = c("red", "blue"), bty = "n")
# finding the name and the index of the rows in cts where the expression levels of the genes are on different sides of that line:
threshold_untreated = exp(10.9)
threshold_treated = exp(11.1)
selected_indices <- which(cts_normalized[, "treated1"] > threshold_treated & 
                            cts_normalized[, "untreated1"] < threshold_untreated)
selected_indices_reverse_way <- which(cts_normalized[, "treated1"] < threshold_treated & 
                            cts_normalized[, "untreated1"] > threshold_untreated)
selected_indices <- c(selected_indices, selected_indices_reverse_way)
selected_genes_names <- rownames(cts)[selected_indices]
selected_genes_indices <- unname(selected_indices)
print(paste("Indices of the selected genes:", paste(selected_indices, collapse = ", ")))
print(paste("The selected genes:", paste(selected_genes_names, collapse = ", ")))

selected_gene <- selected_genes_indices[1]
expression_levels_of_the_detected_gene <- cts_normalized[selected_gene, ]

plot(expression_levels_of_the_detected_gene,
     main = "Expression Levels of the detected Gene", 
     xlab = "Samples", ylab = "Expression Level")

## install the package 'DESeq2' 
# copy and paste the lines below into R 
if (!requireNamespace("BiocManager", quietly = TRUE))
{install.packages("BiocManager") 
  BiocManager::install("DESeq2")} 
## use the installed library 
# copy and paste the line below into R
library("DESeq2") 

## Detect the probability that each one of the genes is differentially expressed between the two conditions 
# copy and paste the following lines into R
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition) 
dds <- DESeq(dds) 
res <- results(dds) 
res
## Sort the genes according to the p-value, such that lower p-values will appear first

res_without_na <- res[!is.na(res$pvalue), ] # remove NA values
ordered_res <- res_without_na[order(res_without_na$pvalue, decreasing = FALSE), ]
p_values <- ordered_res$pvalue
# Top 10 differ genes
ten_most_divide_genes <- head(ordered_res, 10)

log2FoldChange_vector <- ten_most_divide_genes$log2FoldChange
padj_vector <- ten_most_divide_genes$padj

# Check if the visually inspected gene is among the top 10 most differentially expressed genes
is_in_top_10 <- selected_genes_names %in% rownames(ten_most_divide_genes)

# Display the result
#if (is_in_top_10) {
  #print(paste("The visually inspected gene (", selected_genes_names, ") is among the 10 genes."))
#} else {
  #print(paste("The visually inspected gene (", selected_genes_names, ") is NOT among the 10 genes."))
#}


##################
##### Part 5 #####
##################
# Load the Circadian RNA-seq data into a matrix -> data :
file = "/Users/adioiz/Documents/Learnings/Neuro-Genomics Course/Home Exc 1/CircadianRNAseq.csv"
data <- as.matrix(read.csv(file, header = TRUE))

per1a_data <- data[data[, "GeneSymbol"] == "per1a", ]
# Extract the time-course data (columns 2 to 13)
expression_levels_of_per1a_gene <- as.numeric(per1a_data[2:13])

# Plot the expression levels of gene per1a
plot(expression_levels_of_per1a_gene, type = "o", 
     xlab = "Time Points", 
     ylab = "Expression Level", 
     main = "Expression of 'per1a' in time points", 
     xaxt = "n", 
     col = "blue", pch = 16, lwd = 2)

# Customize x-axis with time point:
axis(1, at = 1:length(colnames(data)[2:13]), labels = colnames(data)[2:13], las = 2)

# Compute the FFT:
per1a_fft<- fft(expression_levels_of_per1a_gene)
N = length(expression_levels_of_per1a_gene) # number of samples
# Compute the power spectrum (magnitude squared of FFT)
per1a_fft_power <- Mod(per1a_fft)^2
unique_per1a_powers <- per1a_fft_power[2:7] # Remove the zero-frequency and take only the first six
# Normalize the power spectrum:
per1a_power_normalized <- unique_per1a_powers / sum(unique_per1a_powers)

# Generate frequency labels (1/hour for a 24-hour cycle)
sampling_interval <- 4  # 4-hour sampling interval
frequencies <- seq(1 / (N * sampling_interval), 
                   1 / (2 * sampling_interval), 
                   length.out = length(per1a_power_normalized))

# Plot the normalized FFT power
plot(frequencies, per1a_power_normalized, 
     xlab = "Frequency [1/hour]", 
     ylab = "Normalized FFT Power", 
     main = "Normalized FFT Power Spectrum")
abline(v = 1/24, col = "green", lwd = 2, lty = 2)
# Add text below the x-line
text(x = 1/24 + 0.0008, y = 0.05, labels = "Circadian Frequency", pos = 1, col = "green")

## Process all the genes in the dataset and sort them according to the normalized FFT power in the circadian frequency of 1/24

circadian_index <- 2 # fs/N = 1/48 so the circadian frequency is in the second sample (1/24)
circadian_powers <- numeric(nrow(data))

for (i in 1:nrow(data)) {
  expression_levels <- as.numeric(data[i, 2:13])
  
  # Skip genes with all zero expression levels
  if (all(expression_levels == 0)) {
    circadian_powers[i] <- NA
    next
  }
  
  # Perform FFT and Compute the power spectrum:
  circadian_fft <- fft(expression_levels)
  circadian_power <- Mod(circadian_fft)^2
  unique_circadian_powers <- circadian_power[2:7]
  circadian_power_normalized <- unique_circadian_powers / sum(unique_circadian_powers)
  
  frequencies <- seq(1 / (N * sampling_interval), 
                     1 / (2 * sampling_interval), 
                     length.out = length(circadian_power_normalized))
  
  circadian_powers[i] <- circadian_power_normalized[circadian_index]
}

# Sort the circadian powers and get the indices of the top 10
sorted_indices <- order(circadian_powers, decreasing = TRUE, na.last = TRUE)

top_10_circadian_genes <- data[sorted_indices[1:10], "GeneSymbol"]
top_40_circadian_genes <- data[sorted_indices[1:40], "GeneSymbol"]
# Print the top 10 genes
print("Top 10 genes with the highest normalized FFT power at 1/24:")
print(top_10_circadian_genes)

write.csv(circadian_powers, "/Users/adioiz/Documents/Learnings/Neuro-Genomics Course/Home Exc 2/g_factors.csv", row.names=FALSE)



## sanity check
for (i in sorted_indices[1:10]){
  expression_levels <- as.numeric(data[i, 2:13])
  # Perform FFT and Compute the power spectrum:
  fft <- fft(expression_levels)
  circadian_power <- Mod(fft)^2
  unique_powers <- circadian_power[2:7]
  power_normalized <- unique_powers / sum(unique_powers)
  
  frequencies <- seq(1 / (N * sampling_interval), 
                     1 / (2 * sampling_interval), 
                     length.out = length(power_normalized))
  plot(frequencies, power_normalized, 
       xlab = "Frequency [1/hour]", 
       ylab = "Normalized FFT Power", 
       main = "Normalized FFT Power Spectrum")
  abline(v = 1/24, col = "green", lwd = 2, lty = 2)
}


##################
##### Part 6 #####
##################

# Create a numerical matrix from all the count data
data <- read.csv(file, header = TRUE)
data_numeric_vector <- as.numeric(as.matrix(data[, 2:13]))
data_numeric <- matrix(data_numeric_vector, nrow = nrow(data), ncol = 12)
colnames(data_numeric) <- colnames(data[2:13])

# Calculate the variance and mean for each gene in all the time points
mean_data <- rowMeans(data_numeric, na.rm=TRUE)
variance_data <- apply(data_numeric, 1, var, na.rm=TRUE) # Calculate the variance

# Log-transform the mean and variance
log_mean_data <- log(mean_data + 1)
log_variance_data <- log(variance_data + 1)

# Sort all the vectors expressions in ascending order:
sorted_indices <- order(mean_data, decreasing = FALSE)

sorted_mean_data <- mean_data[sorted_indices]
sorted_variance_data <- variance_data[sorted_indices]

## Bin the data
# Define the bins and range
num_bins <- 20
min_value <- 3
max_value <- max(mean_data)
bins <- seq(min_value, max_value, length.out = num_bins + 1)
bin_indices <- cut(mean_data, breaks = bins, include.lowest = TRUE, labels = FALSE)

# Calculate z-scores for each bin:
z_scores <- numeric(length(mean_data)) # vector to store z-scores
for (bin in 1:num_bins) {
  # Get the indices of genes in the current bin
  genes_in_bin <- which(bin_indices == bin)
  
  # Skip bins with no genes
  if (length(genes_in_bin) == 0) next
  
  variances_in_bin <- sorted_variance_data[genes_in_bin]
  mean_variance_bin <- mean(variances_in_bin)
  sd_variance_bin <- sd(variances_in_bin) # standard deviation of the variances
  
  z_scores[genes_in_bin] <- (variances_in_bin - mean_variance_bin) / sd_variance_bin
}

# Sort z-scores in descending order
sorted_indices <- order(z_scores, decreasing = TRUE, na.last = TRUE)
z_scores <- z_scores[sorted_indices]

top_40_z_scores <-z_scores[1:40]
data_sorted_by_z_scores_top_40 <- data[sorted_indices[1:40], "GeneSymbol"]

common_genes <- intersect(top_40_circadian_genes, data_sorted_by_z_scores_top_40)






