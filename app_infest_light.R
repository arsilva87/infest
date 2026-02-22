
# The INFEST app: insect feeding behavior statistics, the light version
# shinylive::export(appdir = getwd(), destdir = "docs")
# httpuv::runStaticServer("docs")

library(shiny)
library(htmltools)
library(gamlss)

# --------------------------------------------------
# A function to process the insect feeding events from the EPG system

epg <- function(filesIn, time_trim = NULL)
{
  files <- filesIn$name
  nfiles <- length(files)

  # aux function for reading .ANA files
  stylet_codes <- c("1" = "np", "2" = "C", "3" = "E1e",
                    "4" = "E1", "5" = "E2", "6" = "F",
                    "7" = "G", "8" = "pd", "9" = "II-2",
                    "10" = "II-3", "11" = "w11", "12" = "w12",
                    "99" = "T")
  read_ana <- function(f) {
    d <- try(read.table(f, sep = "\t", fileEncoding = "UTF-16"),
             silent = TRUE)
    if(inherits(d, "try-error")) d <- read.table(f, sep = "\t",
                                                 fileEncoding = "UTF-8")
    d$V1 <- stylet_codes[as.character(d$V1)] # waveform
    d$V2 <- as.double(gsub(",", ".", d$V2))  # time
    nr <- nrow(d)
    dl <- lapply(1:nr, function(i) {
      rbind(d$V1[i], d$V2[i])
    })
    unlist(dl)
  }

  f_ext <- tools::file_ext( tolower(filesIn$datapath) )
  ins  <- list()
  for(i in 1:nfiles) {
    ins[[i]] <- try(read_ana(filesIn$datapath[i]), silent = TRUE)
    if(inherits(ins[[i]], "try-error")) {
      ins[[i]] <- scan(filesIn$datapath[i],
                       what = "double", sep = ",", flush = TRUE, quiet = TRUE)
      ins[[i]] <- gsub(" ", "", ins[[i]], fixed = TRUE)
    }

  }
  names(ins) <- files
  events <- sapply(ins, length)/2

  check <- (2*events) %% 2 != 0  # is odd = maybe missing lines
  o <- which(check == TRUE)
  if (any(check == TRUE)) {
    ins[o] <- NULL
    showNotification(paste(" Error (please check line labels) in file(s):",
                           paste(files[o], collapse = ", ")),
                     type = "error", duration = 10)
  }

  fcheck2 <- function(x) {
    aux <- seq.int(1, length(x), by = 1)
    w <- x[aux %% 2 != 0]
    n <- length(w)
    any(w[-n] == w[-1])
  }
  check2 <- sapply(ins, fcheck2)
  o2 <- which(check2 == TRUE)
  if (any(check2 == TRUE)) {
    showNotification(paste("The following file(s) seem(s) to have labels consecutively repeated :",
                           paste(files[o2], collapse = ", ")),
                     type = "warning", duration = 10)
  }

  # time trim (optional)
  if(!is.null(time_trim)) {
    ins <- lapply(ins, function(x) {
      time_rec <- as.double(x)
      o <- which(time_rec >= time_trim[1] & time_rec <= time_trim[2])
      if(length(o) < 2) {
        showNotification(" Time interval is too short! ",
                         type = "error", duration = 5)
      }
      w <- seq.int(min(o)-1, max(o), by = 1)
      x[w]
    })
  }

  # aux function to get odd indexes from elements of ins list
  fodd <- function(x) {
    n <- length(x)
    aux <- seq.int(1, n, by = 1)
    y <- x[aux %% 2 != 0]
    y
  }
  waves <- lapply(ins, fodd)
  levels <- lapply(waves, unique)
  lev <- unique(unlist(levels))
  waves <- lapply(waves, factor, levels = lev)
  table_func <- function(x) {
    out <- table(x)
    names(out) <- paste0("NWEi_", names(out))
    out
  }
  waveevents <- t(sapply(waves, table_func))   # output and export csv

  # aux function to compute the seq variables
  countseq <- function(x) {
    n <- length(x)
    tabpair <- table(x[-n], x[-1], dnn = list("from", "to"))
    diag(tabpair) <- NA
    dfpair <- data.frame(tabpair)
    dfpair <- subset(dfpair, from != to)
    pares <- paste0(dfpair$from, "-", dfpair$to)
    count <- dfpair$Freq
    names(count) <- pares
    count
  }
  seqevents <- t(sapply(waves, countseq))   # output and export csv

  # get the sequential events as a 2-way table count for network graph
  tableseq <- function(x) {
    n <- length(x)
    tabseq <- table(x[-n], x[-1], dnn = c("from", "to"))
    diag(tabseq) <- table(x)
    tabseq
  }
  tabseq <- lapply(waves, tableseq)
  sumtabseq <- Reduce("+", tabseq)

  # transition matrix
  countseq_m <- function(x) {
    n <- length(x)
    tabpair <- table(x[-n], x[-1], dnn = list("from", "to"))
    tabpair
  }
  tran_mats <- lapply(waves, countseq_m)
  tran_mat <- Reduce("+", tran_mats)

  # DURATION - aux funct using even and odd indexes at each element of the list
  #nonprobe <- c("Z", "z", "np", "NP", "nP", "Np")
  w_nonprobe <- regexpr(c("np", "z"), lev, ignore.case = TRUE) > 0
  nonprobe <- lev[w_nonprobe]
  probe_waves <- lev[!(lev %in% nonprobe)]

  evenodd <- function(x, arq) {
    n <- length(x)
    aux <- seq.int(1, n, by = 1)
    waves <- factor(x[aux %% 2 != 0], levels = lev)
    cumtimes <- as.numeric(x[aux %% 2 == 0])
    if(arq == "ana") {
      times <- diff(cumtimes)
      n2 <- n/2
      out <- tapply(times, waves[-n2], FUN = sum)
    } else {
      times <- c(cumtimes[1], diff(cumtimes))
      out <- tapply(times, waves, FUN = sum)
    }
    names(out) <- paste0("WDi_", names(out))
    out
  }
  times <- t(sapply(ins, evenodd, arq = f_ext[1]))
  times[is.na(times)] <- 0

  # duration per event
  de <- times/waveevents
  de[is.na(de)] <- 0
  if(nrow(de) > 1 | !is.null(dim(de))) {
    colnames(de) <- gsub("WDi", "WDEi", colnames(de))
  } else {
    names(de) <- gsub("WDi", "WDEi", names(de))
  }

  # duration of first or second event 'W'
  DurFS_func <- function(x, pos, arq) {
    n <- length(x)
    aux <- seq.int(1, n, by = 1)
    waves <- factor(x[aux %% 2 != 0], levels = lev)
    cumtimes <- as.numeric(x[aux %% 2 == 0])
    if(arq == "ana") {
      times <- diff(cumtimes)
    } else {
      times <- c(cumtimes[1], diff(cumtimes))
    }
    DurFSW <- sapply(lev, function(w) {
      o <- which(waves %in% w)[pos]
      times[o]
    })
    names(DurFSW) <- paste0("DurFrs_", names(DurFSW))
    DurFSW
  }
  DurFrsW <- t(sapply(ins, DurFS_func, pos = 1, arq = f_ext[1]))
  DurFrsW[is.na(DurFrsW)] <- 0

  # number of events to first event W
  NumEvFrsW_fun <- function(x) {
    out <- sapply(probe_waves, function(w) {
      r <- which(w == x)[1] - 1
      names(r) <- paste0("NumEvFrs_", w)
      r
    }, USE.NAMES = FALSE)
    out
  }
  NumEvFrsW <- t(sapply(waves, NumEvFrsW_fun))

  final_data <- data.frame(times,
                           waveevents,
                           de,
                           DurFrsW,
                           NumEvFrsW)
  # output
  out <- list(duration = final_data,
              seqevents = seqevents,
              tabseqevents = sumtabseq,
              transmatrix = tran_mat)
  return(out)
}

# ----------------------------------------------------------
# aux function

f_rmse <<- function(model) {
  res = model[["y"]] - fitted(model)
  sqrt(mean(res^2))
}

zerofams <<- c("NO", "PARETO2", "ZAGA", "ZAIG",
               "PO", "NBII", "ZIP", "ZINBI", "PIG")

# ----------------------------------------------------------
# aux function to automatically fit all models and gen. report

lrt_func <- function(mod) {
  l <- logLik(mod)
  mod0 <- update(mod, formula = ~1, sigma.formula = ~1,
                 family = mod$family[1])
  l0 <- logLik(mod0)
  X2 <- 2*c(l - l0)
  grau <- attr(l, "df") - attr(l0, "df")
  pval <- pchisq(X2, grau, lower = F)
  if(inherits(mod, "try-error")) {
    out <- list(X2 = NA, dfr = NA, pvalue = NA)
  } else {
    out <- list(X2 = X2, dfr = grau, pvalue = pval)
  }
  class(out) <- "lrt_func"
  out
}

best_family_full <- function(y, Group, conf = 0.95, adj = "none") {
  fams <- c("NO", "EXP", "GA", "IGAMMA", "IG",
            "WEI", "PARETO2", "ZAGA", "ZAIG",
            "PO", "NBII", "ZIP", "ZINBI", "PIG")
  m1 <- m2 <- list()
  nf <- length(fams)
  d <- data.frame(y = y, Group = Group)
  d <- d[complete.cases(d), ]
  assign("ddd", d, envir = .GlobalEnv)
  for(i in 1:nf) {
    dat <- get("ddd", envir = .GlobalEnv)
    m1[[i]] <- try(gamlss(dat$y ~ dat$Group, sigma.formula = ~1,
                          family = fams[i],
                          #data = get("ddd", envir = .GlobalEnv),
                          trace=FALSE),
                   silent = TRUE)
    m2[[i]] <- try(gamlss(dat$y ~ dat$Group, sigma.formula = ~Group,
                          family = fams[i],
                          #data = get("ddd", envir = .GlobalEnv),
                          trace=FALSE),
                   silent = TRUE)
  }
  gams <- c(m1, m2)
  aics <- sapply(gams, function(w) {
    ifelse(inherits(w, "try-error"), NA, AIC(w))
  })
  o_aic <- which.min(aics)
  if(length(o_aic) == 0) o_aic = 1
  mod <- gams[[o_aic]]
  if(!inherits(mod, "try-error")) {
    mod$family <- gams[[o_aic]]$family
    mod$call[4] <- mod$family[1]
    lrt <- lrt_func(mod)
    dfr <- mod$df.residual
    aic <- mod$aic
    rmse <- f_rmse(mod)
  } else {
    mod$family <- c(NA, NA)
    mod$call[4] <- NA
    lrt <- list(X2 = NA, dfr = NA, pvalue = NA)
    dfr <- NA
    aic <- NA
    rmse <- NA
  }
  out <- list(dtf = d,
              fit = mod,
              lrt = lrt,
              conf = conf,
              adj = adj,
              dfr = dfr,
              family = mod$family,
              aic = aic,
              rmse = rmse)
  class(out) <- "infest"
  return(out)
}


# server ---------------------------------------------------
server <- function(input, output, session){
  filesIn <- reactive({
    validate(
      need(input$files != "", "...waiting for the input files")
    )
    input$files
  })

  observeEvent(input$showSidebar, {
    shinyjs::show(id = "Sidebar")
  })
  observeEvent(input$hideSidebar, {
    shinyjs::hide(id = "Sidebar")
  })

  tab <- eventReactive(input$run, {
    epg(filesIn(), time_trim = input$time_range * 60)
  })

  # ------------------------------------- duration
  dur <- reactive({
    df <- tab()$duration
    n <- nrow(df)
    df <- data.frame(Group = as.integer(rep(1, n)), round(df, 2))
    df
  })

  output$duration <- rhandsontable::renderRHandsontable({
    df_den <- data.frame(Group = dur()[, 1], dur()[, -1])
    updateSelectInput(session, inputId = 'y', choices = names(df_den)[-1])
    rhandsontable::rhandsontable(df_den, rowHeaderWidth=180)
  })
  mt <- reactive({
    DF = rhandsontable::hot_to_r(input$duration)
    if(!is.null(input$tab$changes$changes)) {
      row.no <- unlist(input$tab$changes$changes)[1]
      col.no <- unlist(input$tab$changes$changes)[2]
      new.val <- unlist(input$tab$changes$changes)[4]
    }
    DF
  })

  # --------------------------------------------------------------
  # statistical report
  observeEvent(input$report, {
    grupos <- as.factor(mt()[, 1])
    req(nlevels(grupos) > 1)
    variables <- tab()$duration
    dtf <- data.frame(Group = grupos, variables)
    pval_adj <- ifelse(input$tukey == TRUE, "tukey", "none")
    conf_lev <- input$conf1
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    nvar <- ncol(dtf)-1
    gams_l <- lapply(2:ncol(dtf), function(i) {
      out<- try(best_family_full(y = dtf[, i], Group = dtf$Group,
                                 conf = conf_lev, adj = pval_adj))
      progress$set(message = "Fitting models...", value = i/nvar,
                   detail = paste("variable", colnames(dtf)[i]))
      out
    })
    names(gams_l) <- colnames(dtf[, -1])
    fit_l <- lapply(gams_l, "[[", "fit")
    means_l <- lapply(fit_l, function(x) {
      try(emmeans::emmeans(x, specs = "Group", type = "response",
                           level = gams_l[[1]]$conf),
          silent = TRUE)
    })
    comp._l <- lapply(means_l, function(x) {
      try(multcomp::cld(x, Letters = letters,
                        adjust = pval_adj,
                        level = conf_lev,
                        alpha = 1-conf_lev),
          silent=TRUE)
    })
    ff <- function(x) {
      if(!inherits(x, "try-error")) {
        x <- x[order(x$Group), ]
        colnames(x)[2] <- "emmean"
        m <- ifelse(x$emmean > 100, round(x$emmean), round(x$emmean, 2))
        s <- ifelse(x$SE > 100, round(x$SE), round(x$SE, 2))
        g <- gsub(" ", "", x$".group")
        out <- paste(m, "±", s, g)
        names(out) <- paste0("Group_", x$Group)
      } else {
        out <- rep(NA, nlevels(dtf$Group))
      }
      out
    }
    all_mc. <- as.data.frame(t(sapply(comp._l, ff)))
    all_mc <- data.frame(variable = rownames(all_mc.), all_mc.)

    comp_l <- lapply(means_l, function(x) {
      try(pairs(x, adjust = pval_adj),
          silent=TRUE)
    })
    ff2 <- function(x) {
      if(!inherits(x, "try-error")) {
        d <- data.frame(x)
        w_test <- gsub(".ratio", "", colnames(d)[ncol(d)-1])
        out <- d[, c(1, 2, 4, 4, ncol(d))]
        out[, 4] <- w_test
      } else {
        out <- data.frame(NA, NA, NA, NA, NA)
        w_test <- NA
      }
      rownames(out) <- NULL
      colnames(out) <- c("contrast", "estimate/ratio", "df", "test", "p_value")
      out
    }
    comp_l2 <- lapply(comp_l, ff2)
    all_mc2. <- do.call(rbind, comp_l2)
    all_mc2 <- data.frame(variable = rownames(all_mc2.), all_mc2.)
    showNotification("DONE!", duration = 2, type = "message")

    # output report
    lrt_l <- lapply(gams_l, "[[", "lrt")
    out_d <- data.frame(variable = names(gams_l),
                        LRT = sapply(lrt_l, "[[", "X2"),
                        Df = sapply(lrt_l, "[[", "dfr"),
                        p_value = sapply(lrt_l, "[[", "pvalue"),
                        distribution = sapply(lapply(gams_l, "[[", "family"), "[", 2),
                        AIC = sapply(gams_l, "[[", "aic"),
                        RMSE = sapply(gams_l, "[[", "rmse")
    )
    rownames(out_d) <- NULL

    # A helper function to turn a data frame into a pure HTML table string
    df_to_html <- function(df) {
      thead <- paste0("<th>", names(df), "</th>", collapse = "")
      tbody <- apply(df, 1, function(x) {
        paste0("<td>", x, "</td>", collapse = "")
      })
      tbody <- paste0("<tr>", tbody, "</tr>", collapse = "")
      return(paste0("<table border='1'><thead>", thead, "</thead><tbody>", tbody, "</tbody></table>"))
    }

    output$report_display <- renderUI({
      tagList(
        div(style = "padding: 20px; border: 1px solid #ddd; background-color: #f9f9f9; border-radius: 8px;",
            h2("Statistical Analysis Report", style = "color: #2c3e50;"),
            hr(),
            h3("Likelihood Ratio Tests for the 'Group' factor"),
            HTML(df_to_html(out_d)),
            h3("Summary of multiple comparisons of means (± SE) of 'Group' levels"),
            HTML(df_to_html(all_mc)),
            h3("Pairwise comparisons of means"),
            HTML(df_to_html(all_mc2)),
            br(),
            tags$footer(
              p(em("This document was automatically generated with INFEST Light")),
              p(strong("Date and time:"), Sys.time())
            )
        )
      )
    })
  })
}


# ----------------------------------------------------------------
ui = navbarPage(title = tags$head(img(src = "https://github.com/arsilva87/infest/blob/5b234a69a2a96aa83323d17643c2d8fff98b128c/inst/infest/www/infest_2_0.png?raw=true", height = 65),
                                  "Insect Feeding Behavior Statistics - light version"),
                windowTitle = "INFEST Light",
                shinyjs::useShinyjs(),
                div(style = "margin-top:-20px"),
                # response variables  ------------------------------------------------------
                tabPanel("Variables/Home", icon = icon("bug"),
                         actionButton("hideSidebar", "Hide sidebar",
                                      style='padding:4px; font-size:80%;  background-color: lightyellow'),
                         actionButton("showSidebar", "Show sidebar",
                                      style='padding:4px; font-size:80%;  background-color: lightgreen'),
                         sidebarLayout(
                           div( id ="Sidebar",
                                sidebarPanel(width = 3,
                                             fileInput("files", "Select one or more EPG files", multiple = TRUE),
                                             sliderInput("time_range", "Timeline trim (min.)",
                                                         min = 0, max = 7200, value = c(0, 7200)),
                                             fluidRow(
                                               column(3,
                                                      tags$br(),
                                                      actionButton("run", "Run", icon = icon("r-project"),
                                                                   style="color: white; background-color: #2e6da4; border-color: white")
                                               )),
                                             tags$small("Press F5/Refresh page to restart the app")
                                )
                           ),
                           mainPanel(width = 9,
                                     h5("Response variables by insect"),
                                     tags$code("WARNING: only the column 'Group' should be edited!"),
                                     rhandsontable::rHandsontableOutput("duration")
                           )
                         )
                ),
                # Report --------------------------------------------
                tabPanel("Stats report", icon = icon("person-running"),
                         h5("NO SWEAT!"),
                         h5("Find the best-fitting model for each response variable
                            as a function of the factor 'Group', perform multiple comparisons of means
                            and create a statistical analysis report."),
                         tags$br(),
                         numericInput("conf1", "Confidence level",
                                      value = 0.95, min = 0, max = 1, step = 0.01),
                         checkboxInput("tukey", "Tukey's p-value adjustment", value = TRUE),
                         tags$code("WARNING: this is a time-consuming task"),
                         tags$br(),
                         actionButton("report", "Fit models", icon = icon("r-project"),
                                      style="color: white; background-color: #2e6da4; border-color: white"),
                         tags$br(),
                         uiOutput("report_display")
                ),
                tabPanel("About", icon = icon("info-circle"),
                         tags$p("Discover the full version of INFEST and know more from: ",
                                tags$a(
                                  "arsilva87.github.io/infest/",
                                  href = "https://arsilva87.github.io/infest/",
                                  target = "_blank")
                         )
                )
)

# Run the app
# runApp()   # to run locally
shinyApp(ui = ui, server = server)


