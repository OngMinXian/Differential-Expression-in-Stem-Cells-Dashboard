# Differential Expression in Stem Cells Dashboard



https://user-images.githubusercontent.com/95210624/230481990-ab814381-e5a7-4448-8cd0-4fe44c3cb3f9.mp4


 
### About the project 
--- 
Stem cell differentiation is the process whereby a stem cell develops into a more specialized cell. In this project, human embryonic stem cells and their differentiated stem cells, namely Mesodermal, Ectodermal and Endodermal, will be compared using differential expression analysis in order to quantify the changes in expression levels between them and find genes of interest. By implementing a dashboard using Shiny (R package), I aim to provide users with multiple ways to visualize the results of the differential expression analysis.

###  Features
--- 
**Sidebar.** Users can select a differentiated stem cell type to compare with undifferentiated stem cells. Users can select the significance value for the false discovery rate (FDR) and the minimum log2(fold change) value which in turn will result in varying genes to be identified as significant. Significant genes can be considered up-regulated or down-regulated and users can choose from either option or both. Genomic track options allow different plot types to be plotted to allow for different visualization methods.

**MA plot.** log2(fold change) is plotted against mean expression (average of the normalized count values) of every gene. Points with adjusted p-values lesser than the significance level and absolute log2(fold change) greater than the stated minimum were considered significant and would be colored (Blue for up-regulated genes, green for down-regulated genes and red for both was selected). Running DESeq2 repeatedly upon each user input was observed to result in an inefficient application. Hence, the results of the DESeq2 analysis were saved into a csv file for quicker plotting.

**Volcano plot.** -log10(adjusted p-value) was plotted against log2(fold change). Likewise, significant genes were colored. Using R’s default plotting functions, a scatter plot was used.

**Significant gene data table.** Every significant gene can be found in the paginated data table with a default sorting of the most significant (lowest adjusted p-value) gene. Other information such as gene id, HGNC symbol, base mean, log2(fold change) and adjusted p-value. There is a search functionality as well as changing the number of entries per page. The HGNC symbol was queried using biomaRt. In order to prevent multiple queries which can be quite slow on the application, every gene was queried as a pre-processing step and saved as a csv file for quick access.

**Heatmap.** The regularized logarithmic transformed counts of the top 20 genes up-regulated and top 20 down-regulated genes were plotted in a heatmap (plot shows only up-regulated or down-regulated genes or both depending on selection). Each row represents a gene which is labeled with their HGNC symbol (genes with no HGNC symbol found were labeled with their gene id instead). Like the results of the DESeq2 analysis, the counts were saved as a CSV. Luckily, DESeq2’s rlog and pheatmap worked with matrices and not just DESeq2 result objects. Gene information. Information of the selected gene is displayed, namely Ensembl gene id, HGNC symbol, chromosome name, start and end position and the strand (1 means plus and -1 means minus). As mentioned above about the HGNC symbol, all this information was similarly stored in the same CSV so as to allow for a more efficient application.

**Genomic track.** Using the bigwig files found on the Encode project as well as the gene information retrieved using biomaRt, the transcription expression level of the selected genes can be plotted using GViz. Note that the plot type can be changed in the sidebar to allow for different visualization methods (histogram, step, line and heatmap).

###  Future work
--- 
While there was an attempt to optimize the application, load times are still significantly slow and loading animations were added in to prevent users from feeling that the application broke while it was loading. A suggestion would be to convert the bigwig files into smaller file formats as only the areas of the significant genes are required to plot the genomic tracks or use cloud hosting to host the large files. Different users might encounter different formatting issues due to the application’s UI/UX not customized to the size of each user’s screen dimension but to the plot/table generated. This can be fixed by detecting the user’s screen dimension and calculating an optimal length and width of each element.

###  Installation
--- 
Unforunately due to the large sizes of the bigwig files, downloading the repository would not allow a user to open the application locally. Found in samples.csv contains the experiments in Ensembl where the files are located. If you have any troubles setting it up, feel free to create an issue.
