##################
##### Part 2 #####
##################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pasilla")
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
cts_untreated <- cts[, grep("untreated", colnames(cts))]
cts_treated <- cts[, grep("treated", colnames(cts))]

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

