# BF591_Final_Project
R Shiny application that features multiple bioinformatics processes implemented in R.


#### Data Description:
  - Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls

#### 12.5.1 Sample Information Exploration
The distinct values and distributions of sample information are important to understand before conducting analysis of corresponding sample data. This component allows the user to load and examine a sample information matrix.

Inputs:
  - Sample information matrix in CSV format
 
Shiny Functionalities:

  - [ ] Tab with a summary of the table that includes a summary of the type and values in each column, e.g.: Number of rows: X Number of columns: Y
	
  - [x] Tab with a data table displaying the sample information, with sortable columns

  - [ ] Tab with histograms, density plots, or violin plots of continuous variables.

If you want to make it fancy, allow the user to choose which column to plot and another column to group by!

#### 12.5.2 Counts Matrix Exploration
Exploring and visualizing counts matrices can aid in selecting gene count filtering strategies and understanding counts data structure. This component allows the user to choose different gene filtering thresholds and assess their effects using diagnostic plots of the counts matrix.

Inputs:
  - Normalized counts matrix, by some method or other, in CSV format
  - Input controls that filter out genes based on their statistical properties:
  - Slider to include genes with at least X percentile of variance
  - Slider to include genes with at least X samples that are non-zero
 
Shiny Functionalities:
  - [ ] Tab with text or a table summarizing the effect of the filtering, including: number of samples, total number of genes, number and % of genes passing current filter, number and % of genes not passing current filter
  - [ ] Tab with diagnostic scatter plots, where genes passing filters are marked in a darker color, and genes filtered out are lighter: median count vs variance (consider log scale for plot), median count vs number of zeros
  - [ ] Tab with a clustered heatmap of counts remaining after filtering - consider enabling log-transforming counts for visualization be sure to include a color bar in the legend
  - [ ] Tab with a scatter plot of principal component analysis projections. You may either: allow the user to select which principal components to plot in a scatter plot (e.g. PC1 vs PC2), allow the user to plot the top N principal components as a beeswarm plot, be sure to include the % variance explained by each component in the plot labels

#### 12.5.3 Differential Expression
Differential expression identifies which genes, if any, are implicated in a specific biological comparison. This component allows the user to load and explore a differential expression dataset.

Inputs:
  - Results of a differential expression analysis in CSV format.
  - If results are already made available, you may use those
  - Otherwise perform a differential expression analysis using DESeq2, limma, or edgeR from the provided counts file
 
Shiny Functionalities:
  - [ ] Tab with sortable table displaying differential expression results
Optional: enable gene name search functionality to filter rows of table
  - [ ] Tab with content similar to that described in Assignment 7

#### 12.6 Choose-your-own Adventure
Implement one or more of the following components as part of your app.

