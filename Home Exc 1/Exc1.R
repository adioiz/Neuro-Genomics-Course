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





## Part 4:

plot.new()
log_of_genes_expression_levels_in_the_first_treated_sample <- log(Normalized_cts[,1]+1)
log_of_genes_expression_levels_in_the_first_untreated_sample <- log(Normalized_cts[,4]+1)
plot(log_of_genes_expression_levels_in_the_first_untreated_sample,main="Genes expression levels - treated v.s. untreated, in log scale", col = 'red')
points(log_of_genes_expression_levels_in_the_first_treated_sample, col= 'blue')
# a great way to separate some of the genes here is by using the following line that is set at e^10, for visualization:
visualization_assisting_line <- vector("numeric", length = 15000)
visualization_assisting_line <- visualization_assisting_line+11 # the height can change, if it is needed after normalizing cts.
points(visualization_assisting_line, col ='green')
# finding the name and the index of the rows in cts where the expression levels of the genes are on different sides of that line:
for(index in 1:the_dimensions_of_cts_matrix[1]) #accessing the amount of rows in cts.
{if (Normalized_cts[index,1]>exp(11.1) & Normalized_cts[index,4] < exp(10.9)) # checking by values that are close to 10.
{
  print(index)
  print(rownames(cts)[index])
}
}

# plotting the expression levels of the last gene that was detected in the loop.

plot.new()
expression_levels_of_the_detected_gene <- vector("numeric", length = the_dimensions_of_cts_matrix[2])

for(j in 1:the_dimensions_of_cts_matrix[2]) #accessing the amount of columns in cts.
{expression_levels_of_the_detected_gene[j]<-Normalized_cts[index,j]}
plot(expression_levels_of_the_detected_gene)




## install the package 'DESeq2' 
# copy and paste the lines below into R 
if (!requireNamespace("BiocManager", quietly = TRUE))
{install.packages("BiocManager") 
  BiocManager::install("DESeq2")} 
## use the installed library 
# copy and paste the line below into R
library("DESeq2") 
## use DESeq to detect the probability that each one of the genes is differentially expressed between the two conditions 
# copy and paste the following lines into R
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition) 
dds <- DESeq(dds) 
res <- results(dds) 
res


## not sure about the following lines yet !

#verifying that pvalues are on the 5th column: 
colnames(res)[5]
#creating_a_temp_version:
ordered_res = res
#ordering by the pvalues in increasing order
indices_for_ordering_res<-order(res[,5], decreasing = FALSE)
ordered_res<-res[indices_for_ordering_res,]
#Examining the ordered_res matrix helps understanding that the last relevant row is 12358.
rownames(ordered_res)[12349:12358]
#presenting ‘log2FoldChange’ values:
ordered_res[12349:12358,2]
#presenting ‘padj’ values:
ordered_res[12349:12358,6]

# Part 5:

# Importing the matrix:

Matrix_of_part_5 <- read.csv("/Users/Yuval/Downloads/CircadianRNAseq.csv")
Matrix_of_part_5 <- as.matrix(Matrix_of_part_5)

x_axis_values <-c(2:13)
# After examining the matrix, per1a is discovered at row 1248.

# Next is plotting the data for that gene over the hours:
plot(Matrix_of_part_5[1248,], xaxt='n')
axis(1,at=x_axis_values,labels=c('11 [pM]','3 [AM]','7 [AM]','11 [AM]','3 [PM]','7 [PM]','11 [PM]','3 [AM]','7 [AM]','11 [AM]','3 [PM]','7 [PM]'))

# Converting to the frequency domain:

numeric_version_for_the_matrix <- as.numeric(Matrix_of_part_5)

fft(numeric_version_for_the_matrix)
