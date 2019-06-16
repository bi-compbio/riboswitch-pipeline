library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(formattable)

#' UI of the Riboswitches explorer
#'
#' Provides the controls needed to input a gene id, mechanism, some network
#' distribution algorithms, some community detection algorithms, and a network3D
#' pane to visualize the network. It can also be exported to Cytoscape or PNG.
ui <- dashboardPage(
  dashboardHeader(title="Riboswitches"),
  dashboardSidebar(disable=TRUE),
  dashboardBody(
    fluidRow(
      column(width=9,
        tabBox(width=NULL,
          tabPanel("Summary stats",
            div(style='overflow-x: scroll',
              tableOutput("summaryStatsFDR"),
              tableOutput("summaryStatsFC"),
              tableOutput("summaryStatsFDRFC")
            )
          ),
          tabPanel("Summary differential abundance",
            div(style='overflow-x: scroll',
              h3(textOutput("summaryDETitle")),
              formattableOutput("summaryDE")
            )
          ),
          tabPanel("Differential abundance",
            fluidRow(
              plotlyOutput("MA.fdr"),
              plotlyOutput("MA.fc")
            ),
            DT::dataTableOutput("MA")
          ),
          tabPanel("Motifs", plotOutput("motifs")),
          tabPanel("Networks", plotOutput("networks"))
        )
      ),
      column(width=3,
        box(width=NULL, status="warning", collapsible=TRUE,
          selectInput("project", "Project:", choices=list()),
          selectInput("switch", "Switch:", choices=list()),
          selectInput("contrast", "Contrast:", choices=list()),
          numericInput("FDR", "FDR cutoff:", 0.01),
          numericInput("FC", "log2 FC cutoff (abs):", 0.5),
          hr(),
          fluidRow(column(12, align="center", actionButton("go", "Go!")))
        )
      )
    )
  )
)

#' Server part of the project selection interface
#'
#' This module creates all the dynamic controls of the UI, fills out the UI and
#' provides the reactive functions that give access to the data.
#'
#' @param input Input object from the shiny app
#' @param output Output object from the shiny app
#' @param session Session object from the shiny app
#' @return nothing
server <- function(input, output, session) {

  source("report.helpers.R")
  session$onSessionEnded(stopApp)

  ##
  ## fill the dinamic content of the UI
  ##
  ## read the project locations
  projects <- reactive({
    validate(need(file.exists("app.projects.csv"), "Projects database not found"))
    read.csv("app.projects.csv", comment.char="#")
  })

  ## read project file
  res <- reactive({
    withProgress(message="Loading edgeR results", value=0, {
      p <- projects()$File[projects()$Name == input$project]
      validate(need(file.exists(p), "Project data not found"))
      local({
        load(p)
        res
      })
    })
  })

  # feed the project names, switches and contrasts
  observe({
    updateSelectInput(session, "project", choices=sort(projects()$Name))
  })

  observe ({
    updateSelectInput(session, "switch", choices=sort(names(res())))
    updateSelectInput(session, "contrast", choices=sort(names(res()[[input$switch]]$LRT)))
  })

  ##
  ## Summary stats
  ##
  output$summaryStatsFDR <- renderTable({
    input$go
    FDR <- isolate(input$FDR)
    FC  <- isolate(input$FC)
    switch <- isolate(input$switch)

    if(!is.null(switch)) {
      x <- res()[[switch]]

      sapply(x$LRT, function(x) {
        c(`Survivors`      =sum(x$table$FDR < FDR),
          `(logFC >0) =ON` =sum(x$table$FDR < FDR & x$table$logFC > 0),
          `(logFC <0) =OFF`=sum(x$table$FDR < FDR & x$table$logFC < 0))
      })
    }
  }, rownames=TRUE)

  output$summaryStatsFC <- renderTable({
    input$go
    FDR <- isolate(input$FDR)
    FC  <- isolate(input$FC)
    switch <- isolate(input$switch)

    if(!is.null(switch)) {
      x <- res()[[switch]]

      sapply(x$LRT, function(x) {
        c(`Survivors`      =sum(abs(x$table$logFC) > FC),
          `(logFC >0) =ON` =sum(x$table$logFC >  FC),
          `(logFC <0) =OFF`=sum(x$table$logFC < -1 * FC))
      })
    }
  }, rownames=TRUE)

  output$summaryStatsFDRFC <- renderTable({
    input$go
    FDR <- isolate(input$FDR)
    FC  <- isolate(input$FC)
    switch <- isolate(input$switch)

    if(!is.null(switch)) {
      x <- res()[[switch]]

      sapply(x$LRT, function(x) {
        c(`Survivors`      =sum(x$table$FDR < FDR & abs(x$table$logFC) > FC),
          `(logFC >0) =ON` =sum(x$table$FDR < FDR & x$table$logFC >  FC),
          `(logFC <0) =OFF`=sum(x$table$FDR < FDR & x$table$logFC < -1 * FC))
      })
    }
  }, rownames=TRUE)

  ##
  ## Summary DE
  ##
  output$summaryDETitle <- renderText({
    input$go
    switch   <- isolate(input$switch)
    contrast <- isolate(input$contrast)

    paste(switch, contrast, sep="-")
  })

  output$summaryDE <- renderFormattable({
    input$go
    FDR <- isolate(input$FDR)
    FC  <- isolate(input$FC)
    switch   <- isolate(input$switch)
    contrast <- isolate(input$contrast)

    if(!is.null(switch) && !is.null(contrast)) {
      x <- res()[[switch]]
      y <- res()[[switch]]$LRT[[contrast]]

      # get the top 10 patterns
      i <- order(y$table$FDR)
      i <- i[y$table$FDR[i] < FDR & abs(y$table$logFC) > FC]       # from the ordered table, discard any not DE

      # make a table with information about FC on other contrasts
      top.patterns <- rownames(y$table)[i]
      df <- as.data.frame(do.call(cbind, lapply(x$LRT, function(x) x$table[top.patterns, "logFC"])))
      df$pattern_ <- top.patterns
      df <- df[, c(grep("^pattern_", colnames(df)), grep("^pattern_", colnames(df), invert=TRUE))]
      colnames(df) <- paste(sub("^pattern_", switch, colnames(df)), "switch")

      # format the table
      my_color_bar <- function() {  # function based on formattable::color_bar
        formatter("span", style=function(x) style(display="inline-block",
                                                  direction="rtl",
                                                  `background-color`=ifelse(x > 0, csscolor("lightgreen"), csscolor("lightpink")),
                                                  `border-radius`="4px",
                                                  `padding-right`="2px",
                                                  width=percent(abs(x) / max(abs(x)))))
      }

      if(nrow(df) > 0)
        formattable(df, list(formattable::area(col=2:ncol(df)) ~ my_color_bar()))
      else
        formattable(df)
    } else {
      NULL
    }
  })

  ##
  ## MA plots
  ##
  # FDR
  output$MA.fdr <- renderPlotly({
    input$go
    FDR <- isolate(input$FDR)
    FC  <- isolate(input$FC)
    switch   <- isolate(input$switch)
    contrast <- isolate(input$contrast)

    withProgress(message="Plotting...", value=0, {
      if(!is.null(switch) && !is.null(contrast)) {
        x <- res()[[switch]]$LRT[[contrast]]

        # MA-plots data structure
        df <- {
          out <- x$table[, c("logFC", "FDR", "logCPM")]
          out$logFDR  <- -log10(out$FDR)
          out$label   <- rownames(out)
          out$density <- get.density(out$logCPM, out$logFDR)
          out
        }
        rownames(df) <- NULL

        # do the FC plot
        p <- ggplot(df, aes(x=logCPM, y=logFDR, label=label, color=density)) +
          geom_point(alpha=.5) +
          scale_color_viridis() +
          ggtitle(paste(switch, contrast, sep=" - ")) +
          xlab("mean of normalized counts") +
          ylab("log2 FC") +
          theme_bw()
      } else {
        p <- ggplot() + theme_bw()
      }

      ggplotly(p)
    })
  })

  # FC
  output$MA.fc <- renderPlotly({
    input$go
    FDR <- isolate(input$FDR)
    FC  <- isolate(input$FC)
    switch   <- isolate(input$switch)
    contrast <- isolate(input$contrast)

    withProgress(message="Plotting...", value=0, {
      if(!is.null(switch) && !is.null(contrast)) {
        x <- res()[[switch]]$LRT[[contrast]]

        # MA-plots data structure
        df <- {
          out <- x$table[, c("logFC", "FDR", "logCPM")]
          out$logFDR  <- -log10(out$FDR)
          out$label   <- rownames(out)
          out$density <- get.density(out$logCPM, out$logFC)
          out
        }
        rownames(df) <- NULL

        # do the FC plot
        p <- ggplot(df, aes(x=logCPM, y=logFC, label=label, color=density)) +
          geom_point(alpha=.5) +
          scale_color_viridis() +
          ggtitle(paste(switch, contrast, sep=" - ")) +
          xlab("mean of normalized counts") +
          ylab("log2 FC") +
          theme_bw()
      } else {
        p <- ggplot() + theme_bw()
      }

      ggplotly(p)
    })
  })

  ##
  ## Motifs
  ##
  output$motifs <- renderPlot({
    input$go
    FDR <- isolate(input$FDR)
    FC  <- isolate(input$FC)
    switch   <- isolate(input$switch)
    contrast <- isolate(input$contrast)

    withProgress(message="Plotting...", value=0, {
      if(!is.null(switch) && !is.null(contrast)) {

        x <- res()[[switch]]$LRT[[contrast]]
        x.sig <- rownames(x$table)[x$table$FDR < FDR & abs(x$table$logFC) > FC]

        # calculate the pwm and plot the motif
        if(length(x.sig) < 1) return(NULL) # nothing to cluster if 1 motifs or less

        # calculate hamming distances between significant sequences and cluster
        d <- hamming(x.sig)
        hc <- hclust(as.dist(d))
        height <- round(nchar(gsub("_", "", x.sig[1])) * .75)  # expected similarity between sequences of a cluster
        motif.clusters <- cutree(hc, h=height)

        # new fancy seqlogo plot and cluster
        pfms <- sapply(1:max(motif.clusters), function(i) {
          my.pcm    <- countsMat(names(motif.clusters)[motif.clusters == i])
          my.pcm.o  <- new("pcm", mat=my.pcm, name=paste("CL", i))
          pcm2pfm(my.pcm.o)
        })

        # calculate the max FC per cluster
        fc <- sapply(1:max(motif.clusters), function(i) {
          j <- names(motif.clusters)[motif.clusters == i]
          k <- which.max(abs(x$table[j, ]$logFC))
          x$table[j, ]$logFC[k]
        })

        # align and plot motifs
        plot.piled.motifs(pfms, fc, switch, contrast)
      }
    })
  })

  ##
  ## Networks
  output$networks <- renderPlot({
    input$go
    FDR <- isolate(input$FDR)
    FC  <- isolate(input$FC)
    switch   <- isolate(input$switch)
    contrast <- isolate(input$contrast)

    # calculate the network based on sequence distances and plot
    withProgress(message="Plotting...", value=0, {
      if(!is.null(switch) && !is.null(contrast)) {
        x <- res()[[switch]]$LRT[[contrast]]
        i <- x$table$FDR < FDR & abs(x$table$logFC) > FC
        plot.net.edger.results(rownames(x$table)[i], x$table$logFC[i], switch, contrast)
      }
    })
  })
}

# Run the application
shinyApp(ui=ui, server=server)
