#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)

# Define UI for application that draws a histogram
ui <- fluidPage(
    navbarPage("GOVisual", tabPanel("Kegg",
                                           titlePanel("Upload Entrezid File"),
                                           sidebarLayout(
                                               sidebarPanel(
                                                   fileInput("file1", "Choose txt file"),
                                                   checkboxInput("header", "Header", TRUE),
                                                   downloadButton('downloadPlot', 'Download Plot')
                                                   
                                               ),
                                               mainPanel(
                                                   plotOutput("dotplot1")
                                               )
                                           )
    ),
    tabPanel("Biological Process",
             titlePanel("Upload Entrezid File"),
             sidebarLayout(
                 sidebarPanel(
                     fileInput("file2", "Choose txt file"),
                     checkboxInput("header", "Header", TRUE),
                     downloadButton('downloadPlot2', 'Download Plot')
                     
                 ),
                 mainPanel(
                     plotOutput("dotplot2")
                 )
             )
    ),
    tabPanel("Cellular Components",
             titlePanel("Upload Entrezid File"),
             sidebarLayout(
                 sidebarPanel(
                     fileInput("file3", "Choose txt file"),
                     checkboxInput("header", "Header", TRUE),
                     downloadButton('downloadPlot3', 'Download Plot')
                     
                 ),
                 mainPanel(
                     plotOutput("dotplot3")
                 )
             )
    ),
    tabPanel("Molecular Function",
             titlePanel("Upload Entrezid File"),
             sidebarLayout(
                 sidebarPanel(
                     fileInput("file4", "Choose txt file"),
                     checkboxInput("header", "Header", TRUE),
                     downloadButton('downloadPlot4', 'Download Plot')
                     
                 ),
                 mainPanel(
                     plotOutput("dotplot4")
                 )
             )
    )
    )
)
# Define server logic required to draw a histogram
server <- function(input, output){
    ###########Kegg plot input setting for download#################
    plotInput1 <- reactive({
        file_to_read = input$file1
        if(is.null(file_to_read)){
            return()
        }
        data1 = read.table(file_to_read$datapath, header = input$header)
        enri = enrichKEGG(data1$Entrezid)
        barpl <- barplot(enri, width = 1)
    })
    ###########BP plot input setting for download#################
    plotInput2 <- reactive({
        file_to_read = input$file2
        if(is.null(file_to_read)){
            return()
        }
        data1 = read.table(file_to_read$datapath, header = input$header)
        enri = enrichGO(
            data1$Entrezid,
            org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.2,
            minGSSize = 10,
            maxGSSize = 500,
            readable = FALSE,
            pool = FALSE
        )
        dotpl <- dotplot(enri)
    })
    ###########CC plot input setting for download#################
    plotInput3 <- reactive({
        file_to_read = input$file3
        if(is.null(file_to_read)){
            return()
        }
        data1 = read.table(file_to_read$datapath, header = input$header)
        enri = enrichGO(
            data1$Entrezid,
            org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "CC",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.2,
            minGSSize = 10,
            maxGSSize = 500,
            readable = FALSE,
            pool = FALSE
        )
        dotpl <- dotplot(enri)
    })
    ###########MF plot input setting for download#################
    plotInput4 <- reactive({
        file_to_read = input$file4
        if(is.null(file_to_read)){
            return()
        }
        data1 = read.table(file_to_read$datapath, header = input$header)
        enri = enrichGO(
            data1$Entrezid,
            org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "MF",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.2,
            minGSSize = 10,
            maxGSSize = 500,
            readable = FALSE,
            pool = FALSE
        )
        dotpl <- dotplot(enri)
    })
    output$dotplot1 <- renderPlot({
        file_to_read = input$file1
        if(is.null(file_to_read)){
            return()
        }
        data1 = read.table(file_to_read$datapath, header = input$header)
        enri = enrichKEGG(data1$Entrezid)
        barplot(enri)
        
    })
    ########### Download Kegg pathway plot#####################
    output$downloadPlot <- downloadHandler(
        filename = function() { paste(input$file1, '.pdf', sep='') },
        content = function(file) {
            ggsave(file,plotInput(), width = 8)
        }
    )
    output$dotplot2 <- renderPlot({
        file_to_read = input$file2
        if(is.null(file_to_read)){
            return()
        }
        data1 = read.table(file_to_read$datapath, header = input$header)
        enri = enrichGO(
            data1$Entrezid,
            org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.2,
            minGSSize = 10,
            maxGSSize = 500,
            readable = FALSE,
            pool = FALSE
        )
        dotplot(enri)
        
    })
    ########### Download BP pathway plot#####################
    output$downloadPlot2 <- downloadHandler(
        filename = function() { paste(input$file2, '.pdf', sep='') },
        content = function(file) {
            ggsave(file,plotInput2(), width = 8)
        }
    )
    output$dotplot3 <- renderPlot({
        file_to_read = input$file3
        if(is.null(file_to_read)){
            return()
        }
        data1 = read.table(file_to_read$datapath, header = input$header)
        enri = enrichGO(
            data1$Entrezid,
            org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "CC",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.2,
            minGSSize = 10,
            maxGSSize = 500,
            readable = FALSE,
            pool = FALSE
        )
        dotplot(enri)
        
    })
    ########### Download CC pathway plot#####################
    output$downloadPlot3 <- downloadHandler(
        filename = function() { paste(input$file3, '.pdf', sep='') },
        content = function(file) {
            ggsave(file,plotInput3(), width = 8)
        }
    )
    output$dotplot4 <- renderPlot({
        file_to_read = input$file4
        if(is.null(file_to_read)){
            return()
        }
        data1 = read.table(file_to_read$datapath, header = input$header)
        enri = enrichGO(
            data1$Entrezid,
            org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "MF",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.2,
            minGSSize = 10,
            maxGSSize = 500,
            readable = FALSE,
            pool = FALSE
        )
        dotplot(enri)
        
    })
    ########### Download MF pathway plot#####################
    output$downloadPlot4 <- downloadHandler(
        filename = function() { paste(input$file4, '.pdf', sep='') },
        content = function(file) {
            ggsave(file,plotInput4(), width = 10)
        }
    )
}
# Run the application 
shinyApp(ui = ui, server = server)
