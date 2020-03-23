# INFEST app: insect feeding statistics

# --------------------------------------------------
# checar pacote
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}
packages <- c("shiny")
check.packages(packages)

# --------------------------------------------------
# A function to process the insect feeding events from the EPG system

epg <- function(filesIn) 
{
   files <- filesIn$name
   nfiles <- length(files)

   ins <- list()
   for(i in 1:nfiles) {
      ins[[i]] <- scan(filesIn$datapath[i], 
        what = "double", sep = ",", flush = TRUE, quiet = TRUE)
      ins[[i]] <- gsub(" ", "", ins[[i]], fixed = TRUE)
   }
   names(ins) <- files
   events <- sapply(ins, length)/2

   # aux function to get odd indexes from elements of ins list
   fodd <- function(x) {
      n <- length(x)
      aux <- seq.int(1, n, by = 1)
      x[aux %% 2 != 0]
   }
   waves <- lapply(ins, fodd)
   levels <- lapply(waves, unique)
   omax <- which.max(sapply(levels, length))
   lev <- levels[[omax]]
   waves <- lapply(waves, factor, levels = lev)
   waveevents <- t(sapply(waves, table))   # output and export csv

   # aux funct using even and odd indexes at each element of the list
   evenodd <- function(x) {
      n <- length(x)
      aux <- seq.int(1, n, by = 1)
      cumtimes <- as.numeric(x[aux %% 2 == 0]) / 60   # in min
      waves <- factor(x[aux %% 2 != 0], levels = lev)
      if ( any(c("Z", "z", "NP", "np", "Np", "nP") == waves[1]) ) {
         times <- c(0, diff(cumtimes))
      } else {
         times <- c(cumtimes[1], diff(cumtimes))
      }
      tapply(times, waves, FUN = sum)
   }
   times <- t(sapply(ins, evenodd))
   times[is.na(times)] <- 0   # output and export csv

   # output
   out <- list(events = waveevents, duration = times)
   return(out)
}


# server ---------------------------------------------------
server <- shinyServer(function(input, output){
   filesIn <- reactive({
      validate(
         need(input$files != "", "...waiting for the input files")
      )
      input$files
   })
   tab <- reactive({
      epg(filesIn())
   })
   output$duration <- renderTable({ 
      tab()$duration
   }, digits = 2, spacing = "xs", rownames = TRUE)
   output$events <- renderTable({ 
      tab()$events
   }, digits = 0, spacing = "xs", rownames = TRUE)
   output$downloadData1 <- downloadHandler(
         filename = function() {
            paste("duration", ".csv", sep = "")
         },
         content = function(file) {
            write.csv(tab()$duration, file)
   })
   output$downloadData2 <- downloadHandler(
         filename = function() {
            paste("events", ".csv", sep = "")
         },
         content = function(file) {
            write.csv(tab()$events, file)
   })
   output$graph1 <- renderPlot({
      data1 <- tab()$duration
      if (nrow(data1) < ncol(data1)) {
         showNotification("n is too small for biplot!", type = "error")
      } else {
         acp1 <- prcomp(scale(data1))
         imp <- paste0(100 * round(summary(acp1)$importance[2, 1:2], 4), "%")
         par(mfrow = c(1, 1), las = 1, cex = 0.9, mar = c(4.5, 4.5, 3, 3))
         biplot(acp1, scale = 0, cex = 0.8, 
            xlab = paste0("PC1 (", imp[1], ")"),
            ylab = paste0("PC2 (", imp[2], ")"))
         abline(h = 0, v = 0, lty = 2, col = "gray")
      }
   })
   output$graph2 <- renderPlot({
      data2 <- tab()$events
      if (nrow(data2) < ncol(data2)) {
         showNotification("n is too small for biplot!", type = "error")
      } else {
         acp2 <- prcomp(scale(data2))
         imp2 <- paste0(100 * round(summary(acp2)$importance[2, 1:2], 4), "%")
         par(mfrow = c(1, 1), las = 1, cex = 0.9, mar = c(4.5, 4.5, 3, 3))
         biplot(acp2, scale = 0, cex = 0.8, 
            xlab = paste0("PC1 (", imp2[1], ")"),
            ylab = paste0("PC2 (", imp2[2], ")"))
         abline(h = 0, v = 0, lty = 2, col = "gray")
      }
   })
   output$graph3 <- renderPlot({
      w <- tab()$events
      par(mfrow = c(ceiling(nrow(w)/4), 4), cex = 0.8, mar = c(1, 1, 2, 1))
      for(i in 1:nrow(w)) {
         pie(w[i, ], radius = 0.9, 
            col = topo.colors(ncol(w)), main = rownames(w)[i], cex.main = 0.8)
      }
   })
})


# ----------------------------------------------------------------
ui <- fluidPage(
   # App title
   titlePanel(title=h3(img(src="infest_logo.png", height = 90), 
      "Insect feeding statistics"), 
       windowTitle="INFEST version beta"), 
   helpText("Processing of insect feeding data from the EPG system"),
   # Menu here
   sidebarLayout(
      sidebarPanel(width = 3,
         fileInput("files", "Select one or more EPG files", multiple = TRUE),
         submitButton("Run"),
         tags$hr(),
         h6("Powered by:",
         helpText("Integrated Pest Management Lab."),
         tags$img(src = "agrometrics.jpg", heigth = 110, width = 110),
         tags$img(src = "logoIF.png", heigth = 90, width = 90),
         helpText("Developed by: Da Silva, A. R. (C) Copyright 2020\n
                  <anderson.silva[at]ifgoiano.edu.br>"),
         helpText("Da Silva, A.R. et al. (2020). INFEST: an R package for processing insect feeding data from the EPG system."),
	 helpText("GNU General Public Licence, version 3.0")
         )
      ),
   # Ouput here
   mainPanel(
      tabsetPanel(
         tabPanel("Duration of events",
            column(width = 6,
               downloadButton("downloadData1", "Export table (.csv)"),
               h5("Duration (in minutes) of events per waveform"),
               tableOutput("duration")
            ),
            column(width = 5, br(), br(),
               h5("Principal component biplot"),
               plotOutput("graph1")
            )
         ),
         tabPanel("Events count",
            column(width = 6,
               downloadButton("downloadData2", "Export table (.csv)"),
               h5("Number of events per waveform"),
               tableOutput("events")
            ),
            column(width = 5, br(), br(),
               h5("Principal component biplot"),
               plotOutput("graph2")
            ),
            column(width = 10, br(),
               h5("Pie chart: proportion of events per waveform"),
               plotOutput("graph3")
            )
         ))
   ))
)

# Run the app
# runApp()   # to run locally
shinyApp(ui = ui, server = server)  # to run from a server

