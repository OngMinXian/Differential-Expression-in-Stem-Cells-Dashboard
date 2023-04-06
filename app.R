###################################################
# Import libraries
###################################################
library(shiny)
library(Gviz)
library(rtracklayer)
library(DT)
library(shinyjs)
library(shinydashboard)
library(DESeq2)
library(pheatmap)
library(shinycssloaders)

############################################
# Reads stored results and define functions
############################################
metadata <- read.csv("samples.csv")

counts <- read.csv("counts.csv")
rownames(counts) <- counts$X

res.ecto <- read.csv("Results/ecto.csv")
res.endo <- read.csv("Results/endo.csv")
res.meso <- read.csv("Results/meso.csv")

res.sig.global <- 'None'

plot_MA <- function(res, a, choice, l2fc) {
  a <- as.numeric(a)
  if (choice == 'Both') {
    return(geneplotter::plotMA(data.frame(res$baseMean, res$log2FoldChange , (res$padj < a) & (abs(res$log2FoldChange) > log2(l2fc))), colSig = 'red'))
  } else if (choice == 'Up') {
    return(geneplotter::plotMA(data.frame(res$baseMean, res$log2FoldChange , (res$padj < a) & (res$log2FoldChange > log2(l2fc))), colSig = 'blue'))
  } else {
    return(geneplotter::plotMA(data.frame(res$baseMean, res$log2FoldChange , (res$padj < a) & (res$log2FoldChange < (-log2(l2fc)))), colSig = 'green'))
  }
  
}

plot_volcano <- function(res, a, choice, l2fc) {
  a <- as.numeric(a)
  c <- rep('black', nrow(res))
  if (choice == 'Both') {
    c[(res$padj < a) & (abs(res$log2FoldChange) > log2(l2fc))] <- 'red'
  } else if (choice == 'Up') {
    c[(res$padj < a) & (res$log2FoldChange > log2(l2fc))] <- 'blue'
  } else {
    c[(res$padj < a) & (res$log2FoldChange < (-log2(l2fc)))] <- 'green'
  }
  
  plot(res$log2FoldChange, -log10(res$padj), col = c, pch = 16, cex = 0.3, xlim = c(-4, 4), ylim = c(0, 50))
  
}

ensembl_gene_id <- 'info$ensembl_gene_id'
hgnc_symbol <- 'info$hgnc_symbol'
chromosome_name <- 'info$chromosome_name'
start_position <- 'info$start_position'
end_position <- 'info$end_position'
strand <- 'info$strand'

info_ecto <- read.csv('Gene_info/info_ecto.csv')
rownames(info_ecto) <- info_ecto$ensembl_gene_id 

info_endo <- read.csv('Gene_info/info_endo.csv')
rownames(info_endo) <- info_endo$ensembl_gene_id 

info_meso <- read.csv('Gene_info/info_meso.csv')
rownames(info_meso) <- info_meso$ensembl_gene_id 

############################################
# User Interface
############################################

ui <- dashboardPage(
  
  dashboardHeader(title = "Differentiated HSC expression"),

  dashboardSidebar(
    
    selectInput('cell_type',
                label = "Choose a differentiated stem cell type",
                choices = list('None',
                               'Mesodermal',
                               'Ectodermal',
                               'Endodermal'),
                selected = 'none'),
    
    conditionalPanel(
      condition = "input.cell_type != 'None'",
      
      radioButtons('up_down_choice',
                   label = "Choose up or down regulated genes",
                   choices = c('Both', 'Up', 'Down'),
                   selected = 'Both'),
      
      selectInput('p_value',
                  label = "Choose a significance level for the false discover rate",
                  choices = list(0.1,
                                 0.05,
                                 0.01),
                  selected = 0.1),
      
      numericInput('l2fc',
                   label = "Choose a minimal log2 fold change (abs)",
                   min = 0,
                   max = 10,
                   step = 0.2,
                   value = 1),
      
      selectInput('Gene_id',
                  label = "Select or search for a gene",
                  choices = list('None'),
                  selected = 'None'),
      
      conditionalPanel(
        condition = "input.Gene_id != 'None'",
      
        selectInput('plot_choice',
                    label = 'Choice a plot type for genomic tracks',
                    choices = list('histogram', 'step', 'line', 'heatmap'),
                    selected = 'histogram'))
      
      )
    
    ),
  
  dashboardBody(
    
    conditionalPanel(
      condition = "input.cell_type == 'None'",
      
      fluidRow(
        box(title = 'Dashboard Description', width = 12, solidHeader = TRUE, status = "primary", collapsible = TRUE,
            uiOutput('Description'))
      )
      
    ),
    
    conditionalPanel(
      condition = "input.cell_type != 'None'",
      
      fluidRow(
        column(width = 6,
          box(title = "MA plot", width = NULL, solidHeader = TRUE, status = "primary", collapsible = TRUE,
              withSpinner(plotOutput('MA_plot')))),
        
        column(width = 6,
          box(title = "Volcano plot", width = NULL, solidHeader = TRUE, status = "primary", collapsible = TRUE,
              withSpinner(plotOutput('Volcano_plot'))))
          
      ),
      
      fluidRow(
        column(width = 6,
               box(title = "Significant genes", width = NULL, solidHeader = TRUE, status = "primary", collapsible = TRUE,
                   withSpinner(dataTableOutput('Sig_genes'))),
               
               conditionalPanel(
                 condition = "input.Gene_id != 'None'",
                 box(title = "Selected gene information", width = NULL, solidHeader = TRUE, status = "primary", collapsible = TRUE,
                     withSpinner(tableOutput('info')))
                )
               ),
        
        
        column(width = 6,
               box(title = "Heatmap", width = NULL, solidHeader = TRUE, status = "primary", collapsible = TRUE,
                   withSpinner(plotOutput('heat_map'))))
        
      ),
      
      conditionalPanel(
        condition = "input.Gene_id != 'None'",
        
        fluidRow(
          box(title = "Genomic tracks", width = 12, solidHeader = TRUE, status = "primary", collapsible = TRUE,
              withSpinner(plotOutput('stem_1')), withSpinner(plotOutput('diff_1')))
        )
        
      )
  
    )
    
  )
    
)
  

############################################
# Server Logic
############################################
server <- function(input, output, session) {
  
  output$Description <- renderUI({
    
    updateSelectInput(session, 'Gene_id', select = 'None')
    
    description_text <- '<p><span>For this assignment, our group opted to create this dashboard such that users can better visualise the effects of various parameters on discovery of significantly different genes. The analysis is performed on stem cells as the reference cell line, along with 2 additional differentiated cell lines, mesodermal and ectodermal cells, on top of the endodermal cells explored in the practical.</span></p>
<p><span>This dashboard has multiple components, namely:</span></p>
<ul>
    <li><span>MA Plot</span></li>
    <li>Volcano Plot</li>
    <li>Table containing data of significant genes</li>
    <li>Heat Map plot</li>
    <li>Comparison between BigWig data between stem cell and selected cell type based on selected gene</li>
</ul>
<p><strong><span>Guide on usage of the dashboard</span></strong><strong><span>&nbsp;</span></strong></p>
<ul>
    <li><span>In the sidebar of this dashboard, select the desired cell type to compare against stem cells.</span></li>
    <li><span>After selecting cell type, additional options include:</span>
        <ul>
            <li><span>Choice for up-regulated or down-regulated genes</span></li>
            <li><span>Significance value for false discovery rate</span></li>
            <li><span>Minimum log2 fold change</span></li>
            <li><span>Selection to search for specific gene</span></li>
        </ul>
    </li>
    <li><span>Changes to the additional options will be reflected in the various plots</span></li>
    <li><span>To observe the difference in BigWig data between the various cell lines, do make a specific gene selection in the sidebar under the additional options. Note that the plot will require some time to load.</span></li>
</ul>'
    
    if (input$cell_type == 'None') {
      shinyjs::hide(id = 'Gene_id')
      HTML(description_text)
    }
  })
  
  output$MA_plot <- renderPlot({
    if (input$cell_type == 'Mesodermal') {
      res <- res.meso
    } else if (input$cell_type == 'Ectodermal') {
      res <- res.ecto
    } else if (input$cell_type == 'Endodermal') {
      res <- res.endo
    } else {
      return(NULL)
    }
    return(plot_MA(res, input$p_value, input$up_down_choice, input$l2fc))
  })
  
  output$Volcano_plot <- renderPlot({
    if (input$cell_type == 'Mesodermal') {
      res <- res.meso
    } else if (input$cell_type == 'Ectodermal') {
      res <- res.ecto
    } else if (input$cell_type == 'Endodermal') {
      res <- res.endo
    } else {
      return(NULL)
    }
    return(plot_volcano(res, input$p_value, input$up_down_choice, input$l2fc))
  })
  
  output$Sig_genes <- renderDataTable({
    if (input$cell_type == 'Mesodermal') {
      res <- res.meso
    } else if (input$cell_type == 'Ectodermal') {
      res <- res.ecto
    } else if (input$cell_type == 'Endodermal') {
      res <- res.endo
    } else {
      return(NULL)
    }
    
    res.sig <- res[(res$padj < as.numeric(input$p_value)) & (!is.na(res$padj)),]
    
    if (input$up_down_choice == 'Up') {
      res.sig <- res.sig[res.sig$log2FoldChange > input$l2fc,]
    } else if (input$up_down_choice == 'Down') {
      res.sig <- res.sig[res.sig$log2FoldChange < (-input$l2fc),]
    } else {
      res.sig <- res.sig[abs(res.sig$log2FoldChange) > input$l2fc,]
    }
    
    res.sig.global <<- res.sig[order(res.sig$padj, decreasing = FALSE),]
    res.sig <- res.sig.global
    
    choice <- as.list(c('None', res.sig$gene_id))
    choice <- sapply(choice, function(x) strsplit(x, "[.]")[[1]][1])
    updateSelectInput(session, 'Gene_id', label = "Select or search for a gene", choices = choice)
    shinyjs::show(id = 'Gene_id')
    
    if (input$cell_type == 'Mesodermal') {
      res.sig$hgnc_symbol <- info_meso[choice[-1], 'hgnc_symbol']
    } else if (input$cell_type == 'Ectodermal') {
      res.sig$hgnc_symbol <- info_ecto[choice[-1], 'hgnc_symbol']
    } else if (input$cell_type == 'Endodermal') {
      res.sig$hgnc_symbol <- info_endo[choice[-1], 'hgnc_symbol']
    }
    
    res.sig.ordered <- res.sig[, c('gene_id', 'hgnc_symbol', 'baseMean', 'log2FoldChange', 'padj')]
    res.sig.ordered <- res.sig.ordered[order(res.sig.ordered$padj, decreasing = FALSE),]
    res.sig.ordered$padj <- as.character(signif(res.sig.ordered$padj, 3))
    return(DT::datatable(res.sig.ordered, options = list(lengthMenu = c(5, 10, 20), pageLength = 5)))
  })
  
  output$heat_map <- renderPlot({
    
    if (input$cell_type == 'None') { return() }
    
    c <- input$cell_type
    res.sig <- res.sig.global
    ids <- sapply(res.sig$gene_id, function(x) strsplit(x, "[.]")[[1]][1])

    res.sig.up <- res.sig[res.sig$log2FoldChange > input$l2fc,]
    res.sig.up <- res.sig.up[order(res.sig.up$padj, decreasing = FALSE),][c(1:20),]
    res.sig.down <- res.sig[res.sig$log2FoldChange < (-input$l2fc),]
    res.sig.down <- res.sig.down[order(res.sig.down$padj, decreasing = FALSE),][c(1:20),]
    
    if (input$up_down_choice == 'Up') {
      mat <- counts[res.sig.up$gene_id,]
    } else if (input$up_down_choice == 'Down') {
      mat <- counts[res.sig.down$gene_id,]
    } else {
      mat <- rbind(counts[res.sig.up$gene_id,], counts[res.sig.down$gene_id,])
    }
    
    mat <- mat[,-1]
    
    if (input$cell_type == 'Mesodermal') {
      mat <- mat[,c('stemcell1', 'stemcell2', 'mesodermal1', 'mesodermal2')]
    } else if (input$cell_type == 'Ectodermal') {
      mat <- mat[,c('stemcell1', 'stemcell2', 'ectodermal1', 'ectodermal2')]
    } else if (input$cell_type == 'Endodermal') {
      mat <- mat[,c('stemcell1', 'stemcell2', 'endodermal1', 'endodermal2')]
    }
    
    ids <- sapply(rownames(mat), function(x) strsplit(x, "[.]")[[1]][1])
    if (input$cell_type == 'Mesodermal') {
      symbols <- info_meso[ids, 'hgnc_symbol']
      symbols[is.na(symbols)] <- rownames(mat)[is.na(symbols)]
      rownames(mat) <- symbols
    } else if (input$cell_type == 'Ectodermal') {
      symbols <- info_ecto[ids, 'hgnc_symbol']
      symbols[is.na(symbols)] <- rownames(mat)[is.na(symbols)]
      rownames(mat) <- symbols
    } else if (input$cell_type == 'Endodermal') {
      symbols <- info_endo[ids, 'hgnc_symbol']
      symbols[is.na(symbols)] <- rownames(mat)[is.na(symbols)]
      rownames(mat) <- symbols
    }
    
    mat <- rlog(as.matrix(mat))
    return(pheatmap(mat, scale = 'row'))
    
  })
  
  output$info <- renderTable({
    
    if (input$cell_type == 'Mesodermal') {
      info <- info_meso[info_meso$ensembl_gene_id == input$Gene_id, ]
    } else if (input$cell_type == 'Ectodermal') {
      info <- info_ecto[info_ecto$ensembl_gene_id == input$Gene_id, ]
    } else if (input$cell_type == 'Endodermal') {
      info <- info_endo[info_endo$ensembl_gene_id == input$Gene_id, ]
    } else {
      return()
    }
   
    ensembl_gene_id <<- info$ensembl_gene_id
    hgnc_symbol <<- info$hgnc_symbol
    chromosome_name <<- info$chromosome_name
    start_position <<- info$start_position
    end_position <<- info$end_position
    strand <<- info$strand

    info
  })
  
  output$stem_1 <- renderPlot({
    
    if (input$cell_type == 'None') { return() }
    c <- input$Gene_id
    
    if (strand == '1') {
      stem_1 <- import.bw("ENCSR348EFG/ENCFF022QBF_plus.bigWig", as = "GRanges")
    } else {
      stem_1 <- import.bw("ENCSR348EFG/ENCFF381GKG_minus.bigWig", as = "GRanges")
    }
    
    if (input$plot_choice == 'line') {
      c <- 'l'
    } else if (input$plot_choice == 'step') {
      c <- 's'
    } else if (input$plot_choice == 'histogram') {
      c <- 'histogram'
    } else if (input$plot_choice == 'heatmap') {
      c <- 'heatmap'
    }
    
    track <- DataTrack(stem_1, chromosome = chromosome_name, name = 'Stem cell 1')
    return(plotTracks(track, from = start_position, to = end_position, type = c))
    
  })
  
  output$diff_1 <- renderPlot({
    
    if (input$cell_type == 'None') { return() }
    c <- input$Gene_id
    
    if (strand == '1') {
      if (input$cell_type == 'Mesodermal') {
        diff_1 <- import.bw("ENCSR500UOD/ENCFF388QFN_plus.bigWig", as = "GRanges")
      } else if (input$cell_type == 'Ectodermal') {
        diff_1 <- import.bw("ENCSR851BRK/ENCFF255CRM_plus.bigWig", as = "GRanges")
      } else if (input$cell_type == 'Endodermal') {
        diff_1 <- import.bw("ENCSR002CTR/ENCFF775NPV_plus.bigWig", as = "GRanges")
      }
    } else {
      if (input$cell_type == 'Mesodermal') {
        diff_1 <- import.bw("ENCSR500UOD/ENCFF553BUT_minus.bigWig", as = "GRanges")
      } else if (input$cell_type == 'Ectodermal') {
        diff_1 <- import.bw("ENCSR851BRK/ENCFF247PDX_minus.bigWig", as = "GRanges")
      } else if (input$cell_type == 'Endodermal') {
        diff_1 <- import.bw("ENCSR002CTR/ENCFF462HCZ_minus.bigWig", as = "GRanges")
      }
    }
    
    if (input$plot_choice == 'line') {
      c <- 'l'
    } else if (input$plot_choice == 'step') {
      c <- 's'
    } else if (input$plot_choice == 'histogram') {
      c <- 'histogram'
    } else if (input$plot_choice == 'heatmap') {
      c <- 'heatmap'
    }
    
    track <- DataTrack(diff_1, chromosome = chromosome_name, name = paste(input$cell_type, 'cell', '1'))
    return(plotTracks(track, from = start_position, to = end_position, type = c))
    
  })
  
}

############################################
# Run the app 
############################################
shinyApp(ui = ui, server = server)

