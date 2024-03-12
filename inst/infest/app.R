
# The INFEST app: insect feeding behavior statistics

#options(repos = c("http://cran.rstudio.com/",
#                  "http://github.com/kassambara/",
#                  "https://cran.r-project.org/"))

library(htmltools)
library(shinythemes)
library(igraph)
library(factoextra)
library(DT)
library(rhandsontable)
library(emmeans)
library(shinyjs)
library(gamlss)
library(ggeffects)
library(ggplot2)
library(reshape2)
library(plotly)
library(rmarkdown)
library(multcompView)

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
    if(inherits(d, "try-error")) d <- read.table(f, sep = "\t")
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
  nonprobe <- c("Z", "z", "np", "NP", "nP", "Np")
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

  # number of probes to first event W
  NumPrFrsW_func <- function(x) {
    out <- sapply(lev, function(w) {
      which(x %in% w)[1] - 1
    })
    onp <- which(names(out) %in% nonprobe)[1]
    out <- out[-onp]
    names(out) <- paste0("NumPrFrs_", names(out))
    out
  }
  NumPrFrsW <- t(sapply(waves, NumPrFrsW_func))

  final_data <- data.frame(times,
                           waveevents,
                           de,
                           DurFrsW,
                           NumPrFrsW)
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
# axu function for multiple binomial tests for transitions
transmatrix_test <- function(x) {
  x <- as.matrix(x)
  nwaves <- ncol(x)
  ntrans <- rowSums(x)
  invnt <- diag(ntrans^-1)
  prop_trans <- invnt %*% x
  pvals <- sapply(1:nwaves, function(i) {
    sapply(x[i,], function(y) {
      nn <- ifelse(ntrans[i] == 0, 1, ntrans[i])
      binom.test(y, nn, p = 1/(nwaves-1),
                 alternative = "greater")$p.value
    })
  })
  pvals_m <- round(t(pvals), 4)
  dimnames(pvals_m) <- dimnames(prop_trans) <- dimnames(x)
  out <- list(pvals = pvals_m, prop = prop_trans)
  class(out) <- "transmatrix_test"
  out
}

print.transmatrix_test <- function(x) {
  nwaves <- ncol(x$pvals)
  cat("\n  Multiple exact binomial tests for event transitions\n",
      "\n  Matrix of p-values:\n")
  print(x$pvals)
  cat("\nAlternative hypothesis: true probability of transition is greater than ",
      1/(nwaves-1), ", meaning: transition is not random\n", sep = "")
}

# ----------------------------------------------------------
# aux function to automatically fit all models and gen. report

boxplot_func <- function(df_den) {
  ggplot(df_den, aes(x = Group, y = y)) +
    geom_boxplot() +
    geom_jitter(color = "black", size = 1, shape = 1,
                alpha = 0.7, width = 0.2) +
    stat_summary(fun=mean, geom="point", shape=5, size = 3,
                 show.legend = FALSE, color='black') +
    coord_flip() + theme_bw()
}

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
  gams <- list()
  nm <- 2*length(fams)
  d <<- data.frame(y = y, Group = Group)
  for(i in 1:nm) {
    gams[[i]] <- try(gamlss(y ~ Group, sigma.formula = ~1,
                            family = fams[i], data = na.omit(d),
                            trace=FALSE),
                     silent = TRUE)
    gams[[i+1]] <- try(gamlss(y ~ Group, sigma.formula = ~Group,
                              family = fams[i], data = na.omit(d),
                              trace=FALSE),
                       silent = TRUE)
  }
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
  out <- list(dtf = data.frame(y, Group),
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

  output$duration <- renderRHandsontable({
    df_den <- data.frame(Group = dur()[, 1], dur()[, -1])
    updateSelectInput(session, inputId = 'y', choices = names(df_den)[-1])
    rhandsontable(df_den, rowHeaderWidth=180)
  })
  mt <- reactive({
    DF = hot_to_r(input$duration)
    if(!is.null(input$tab$changes$changes)) {
      row.no <- unlist(input$tab$changes$changes)[1]
      col.no <- unlist(input$tab$changes$changes)[2]
      new.val <- unlist(input$tab$changes$changes)[4]
    }
    DF
  })

  output$graph1 <- renderPlotly({
    dtf_vars <- tab()$duration
    cols <- grep(input$biplot_data, colnames(dtf_vars))
    dtf_vars <- dtf_vars[, cols]
    dtf. <- data.frame(lev = mt()[, 1], dtf_vars)
    dtf <- na.omit(dtf.)
    if( length(unique(dtf$lev)) < 2 ) lev <- NA else lev <- as.factor(dtf$lev)
    if (nrow(dtf) < ncol(dtf)) {
      showNotification("n is too small for a biplot!", type = "error")
    } else {
      acp <- prcomp(scale(dtf[, -1]))
      p <- fviz_pca_biplot(acp,
                           geom.var = c("arrow", "text"),
                           geom.ind = "text",
                           col.ind = "black",
                           addEllipses = TRUE, ellipse.level = input$conf1,
                           fill.ind = lev,
                           palette = topo.colors(nlevels(lev)),
                           col.var = "contrib", fill.var = "contrib",
                           legend.title = list(fill = "Group", color = "Contrib (%)"),
                           title = "")
      ggplotly(p)
    }
  })

  output$graph12 <- renderPlotly({
    dtf <- tab()$duration
    df_den <- data.frame(X = rownames(dtf), Group = as.factor(mt()[, 1]), dtf)
    updateSelectInput(session, inputId = 'barplot_file',
                      choices = c("All", rownames(dtf)),
                      selected = input$barplot_file)
    if(input$barplot_file != "All") {
      df_den <- subset(df_den, X == input$barplot_file)
    }
    dados_m <- reshape2::melt(df_den, id.vars = 1:2)
    p <- ggplot(dados_m, aes(y = variable, x = value, fill = Group)) +
      geom_col() + xlab("") +
      facet_wrap(~X) +
      theme_bw()
    plotly::ggplotly(p)
  })

  observeEvent(input$bxp1, {
    df_den <<- data.frame(Group = as.factor(mt()[, 1]),
                          tab()$duration)
    y <- df_den$y <- df_den[, input$y]
    yname <- input$y
    output$boxplot1 <- renderPlotly({
      p <- ggplot(df_den, aes(x = Group, y = y)) +
        geom_boxplot() + ylab(yname) +
        geom_jitter(color = "black", size = 0.9, shape = 1,
                    alpha = 0.7, width = 0.2) +
        stat_summary(fun=mean, geom="point", shape=5, size = 3,
                     show.legend = FALSE, color='black') +
        coord_flip() + theme_bw()
      ggplotly(p)
    })
  })

  observeEvent(input$fit_glm1, {
    pval_adj <- ifelse(input$tukey == TRUE, "tukey", "none")
    df_den <<- data.frame(Group = as.factor(mt()[, 1]), tab()$duration)
    y <- df_den$y <- df_den[, input$y]
    Group <- df_den$Group
    if(nlevels(Group) < 2) {
      showNotification("'Group' has only one level.
                      Go back to the data table and specify more levels.",
                       type = "error", duration = 5)
    } else {
      if( any(y == 0) & !(input$family1 %in% zerofams) ) {
        showNotification("Only positive values are accepted by this
                        distribution model. See the 'About' menu.",
                         type = "warning", duration = 5)
      }
      assign("familia1", input$family1, envir = .GlobalEnv)
      sigma1 <- ifelse(input$sigma1 == TRUE, expression(~Group), expression(~1))
      mod1 <- try(gamlss(y ~ Group, sigma.formula = eval(sigma1),
                         data = na.omit(df_den),
                         family = get("familia1", envir = .GlobalEnv),
                         control = gamlss.control(trace = F)))
      mod0 <- try(gamlss(y ~ 1, sigma.formula = eval(sigma1),
                         data = na.omit(df_den),
                         family = get("familia1", envir = .GlobalEnv),
                         control = gamlss.control(trace = F)))
      med1 <- try(emmeans::emmeans(mod1, "Group", type = "response",
                                   level = input$conf1))
      comp_med1 <- try(pairs(med1, adjust = pval_adj))
      comp_med1. <- try(multcomp::cld(med1, Letters = letters,
                        level = input$conf1, alpha = 1-input$conf1,
                        adjust = pval_adj)[])
      if(inherits(mod1, "try-error") |
         inherits(mod0, "try-error") |
         inherits(med1, "try-error")) {
        showNotification("Algorithm has not converged. Try again.",
                         type = "error", duration = 5)
      } else {
        res1 <- mod1[["y"]] - fitted(mod1)
        output$aic1 <- renderPrint({ cat("AIC =", AIC(mod1),
                                         "  RMSE =", sqrt(mean(res1^2)), "\n") })
        output$model1 <- renderPrint({
          l <- logLik(mod1)
          l0 <- logLik(mod0)
          X2 <- 2*c(l - l0)
          grau <- attr(l, "df") - attr(l0, "df")
          pval <- pchisq(X2, grau, lower = F)
          cat("Likelihood Ratio Test for the 'Group' factor\n",
              "LRT =", X2, "   df =", grau, "   Pr(Chisq) =", pval)
        })
        output$med1 <- renderPrint({
          if(input$letters1) {
            cat("Pairwise comparisons of means\n")
            comp_med1.
          } else {
            cat("Estimated marginal means\n")
            med1
          }
        })
        output$contrasts1 <- renderPrint({
          if(input$letters1) {
            NULL
          } else {
            cat("Pairwise comparisons of means\n")
            comp_med1
          }
        })
        output$graph13 <- renderPlotly({
          ggplotly(plot(med1, xlab = input$y) +
                     theme_bw(base_size = 12) +
                     ggtitle(paste0(100*input$conf1,
                                    "% confidence intervals for means\n")))
        })
        observeEvent(input$wp_glm1, {
          output$graph_wp1 <- renderPlot( {
            gamlss::wp(mod1)
            title("Worm-plot")
            } )
        })
      }
    }
  })

  # ---------------------------------------- event transition
  output$seqevents <- DT::renderDT({
    df3 <- data.frame(Group = as.factor(mt()[, 1]), tab()$seqevents)
    DT::datatable(df3, editable = F,
                  extensions = 'Buttons',
                  options = list(dom = 'Bfrtip',
                                 digits = 0,
                                 buttons = c('copy', 'csv', 'print'),
                                 text = "Download",
                                 paging = FALSE, searching = TRUE))
  })

  output$transmatrix <- DT::renderDT({
    df_tm <- tab()$transmatrix
    lev <- rownames(df_tm)
    tm <- matrix(df_tm, nrow = nrow(df_tm),
           dimnames = list(paste0("from_", lev),
                           paste0("to_", lev)))
    output$trans_tests <- renderPrint(width = 400, {
      pv <- transmatrix_test(tm)
      output$prop_trans <- DT::renderDT({
        DT::datatable(pv$prop, editable = F,
                      extensions = 'Buttons',
                      options = list(dom = 'Bfrtip',
                                     digits = 4,
                                     buttons = c('copy', 'csv', 'print'),
                                     text = "Download",
                                     paging = FALSE, searching = TRUE))
      })
      print(pv)
    })
    DT::datatable(tm, editable = F,
                  extensions = 'Buttons',
                  options = list(dom = 'Bfrtip',
                                 digits = 0,
                                 buttons = c('copy', 'csv', 'print'),
                                 text = "Download",
                                 paging = FALSE, searching = TRUE))
  })

  output$graph4 <- renderPlot({
    tt <- tab()$tabseqevents
    par(mar = c(1, 1, 1, 1))
    net <- graph_from_adjacency_matrix(tt, weighted = TRUE)
    net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE)
    w <- E(net)$weight/max(E(net)$weight)
    E(net)$width <- 3 * w
    plot(net,
         vertex.color = "yellow",
         vertex.shape = "circle",
         vertex.size = 20 * diag(tt)/max(diag(tt)) + 13,
         edge.color = gray(1-w), edge.arrow.size = 1)
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
      out<- best_family_full(y = dtf[, i], Group = dtf$Group,
                             conf = conf_lev, adj = pval_adj)
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
    comp_l <- lapply(means_l, function(x) {
      try(pairs(x, adjust = pval_adj), silent=TRUE)
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
    all_mc <- as.data.frame(t(sapply(comp._l, ff)))
    showNotification(paste("DONE!"), duration = 3, type = "message")

    output$downloadData <- downloadHandler(
      filename = function() {
        "infest_data.csv"
      },
      content = function(file) {
        write.csv(dtf, file, row.names = TRUE)
      }
    )

    df_tm <- tab()$transmatrix
    lev <- rownames(df_tm)
    tm <- matrix(df_tm, nrow = nrow(df_tm),
                 dimnames = list(paste0("from_", lev),
                                 paste0("to_", lev)))
    trans <- transmatrix_test(tm)

    output$downloadReport <- downloadHandler(
      filename = "infest_report.html",
      content = function(file) {
        showNotification("Please wait", duration = 2, type = "warning")
        rmarkdown::render("./www/report.Rmd",
                          params = list(mods = gams_l,
                                        meds = means_l,
                                        comp1 = comp._l,
                                        all_mc = all_mc,
                                        comp2 = comp_l,
                                        trans = trans),
                          output_file = file,
                          encoding = "UTF-8")
      }
    )
  })
}


# ----------------------------------------------------------------
ui = navbarPage(title = tags$head(img(src="infest_2_0.png", height = 65),
                                  "Insect Feeding Behavior Statistics"),
                theme = shinytheme("yeti"),
                windowTitle = "INFEST 2.0",
                useShinyjs(),
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
                                             tags$hr(),
                                             fluidRow(
                                               column(6,
                                               numericInput("conf1", "Confidence level",
                                                            value = 0.95, min = 0, max = 1, step = 0.01,
                                                            width = 95)),
                                               column(3,
                                               tags$br(),
                                               actionButton("run", "Run", icon = icon("r-project"),
                                                            style="color: white; background-color: #2e6da4; border-color: white")
                                             )),
                                             tags$small("Press F5/Refresh page to restart the app")
                                )
                           ),
                           mainPanel(width = 9,
                                     tabsetPanel(
                                       tabPanel(icon("table"),
                                                h5("Response variables by insect"),
                                                tags$code("WARNING: only the column 'Group' should be edited!"),
                                                rHandsontableOutput("duration"),
                                                tags$h1()),
                                       tabPanel(icon("chart-line"),
                                                h5("Principal component biplot"),
                                                radioButtons("biplot_data", "Variable type",
                                                             choices = c("WDi", "NWEi", "WDEi"), inline = TRUE),
                                                plotlyOutput("graph1"),
                                                tags$h1()),
                                       tabPanel(icon("chart-bar"),
                                                h5("Bar-plots per insect"),
                                                selectInput('barplot_file', 'Files/Insects', choices = "All"),
                                                plotlyOutput("graph12"),
                                                tags$h1()),
                                       tabPanel(icon("calculator"),
                                                h5("Fit 1-way GAMLSS model"),
                                                column(width = 3,
                                                       tags$code("Group factor: 'Group'"),
                                                       selectInput('y', 'Response variable', 'WDi_NP'),
                                                       actionButton("bxp1", "Box-plot", icon = icon("r-project"),
                                                                    style="color: white; background-color: #2e6da4; border-color: white"),
                                                       tags$h2(),
                                                       selectInput('family1', 'Distribution model',
                                                                   choices = list("Gaussian" = "NO",
                                                                                  "Exponential" = "EXP",
                                                                                  "Gamma" = "GA",
                                                                                  "Inverse Gamma" = "IGAMMA",
                                                                                  "Inverse Gaussian" = "IG",
                                                                                  "Weibull" = "WEI",
                                                                                  "Pareto type 2" = "PARETO2",
                                                                                  "Zero adjusted Gamma" = "ZAGA",
                                                                                  "Zero adjusted Inverse Gaussian" = "ZAIG",
                                                                                  "Poisson" = "PO",
                                                                                  "Negative Binomial" = "NBII",
                                                                                  "Zero-Inflated Poisson" = "ZIP",
                                                                                  "Zero-Inflated Negative Binomial" = "ZINBI",
                                                                                  "Poisson-Inverse Gaussian" = "PIG")),
                                                       checkboxInput("sigma1", "Model variances of 'Group'", value = FALSE),
                                                       checkboxInput("tukey", "Tukey's p-value adjustment", value = TRUE),
                                                       checkboxInput("letters1", "Display compact letters", value = FALSE),
                                                       actionButton("fit_glm1", "Fit", icon = icon("r-project"),
                                                                    style="color: white; background-color: #2e6da4; border-color: white"),
                                                       actionButton("wp_glm1", "Worm-plot", icon = icon("r-project"))
                                                ),
                                                column(width = 8,
                                                       conditionalPanel(condition = "input.bxp1",
                                                                        plotlyOutput("boxplot1"),
                                                                        tags$h1()
                                                       ),
                                                       conditionalPanel(condition = "input.bxp1 == FALSE",
                                                                        verbatimTextOutput("aic1"),
                                                                        verbatimTextOutput("model1"),
                                                                        verbatimTextOutput("med1"),
                                                                        verbatimTextOutput("contrasts1"),
                                                                        plotlyOutput("graph13")
                                                       ),
                                                       plotOutput("graph_wp1"),
                                                       tags$h1()
                                                )
                                       )
                                     )
                           )
                         )
                ),
                # Report --------------------------------------------
                tabPanel("Stats report", icon = icon("person-running"),
                         h5("NO SWEAT!"),
                         h5("Find the best-fitting model for each response variable
                            as a function of the factor 'Group', then generate and
                            download a full statistical analysis report."),
                         tags$code("WARNING: this is a time-consuming task"),
                         tags$br(),
                         actionButton("report", "Fit models", icon = icon("r-project"),
                                      style="color: white; background-color: #2e6da4; border-color: white"),
                         downloadButton("downloadData", "Data"),
                         downloadButton("downloadReport", "Report")
                ),
                # Event transitions --------------------------------------------
                tabPanel("Event transitions", icon = icon("sort-alpha-up"),
                         tabsetPanel(tabPanel(icon("table"),
                                              h5("Event transition (from-to) count by insect"),
                                              DTOutput("seqevents"),
                                              tags$h1()),
                                     tabPanel(icon("table"),
                                              h5("Transition matrix"),
                                              DTOutput("transmatrix"),
                                              tags$br(),
                                              h5("Transition probabilities"),
                                              DTOutput("prop_trans"),
                                              tags$br(),
                                              verbatimTextOutput("trans_tests"),
                                              tags$h1()),
                                     tabPanel(icon("project-diagram"),
                                              h5("Network plot"),
                                              plotOutput("graph4"),
                                              tags$h1())
                         )
                ),
                # --------------------------------------------------------------
                tabPanel("Legend", icon = icon("book"),
                         includeHTML('www/legend.html')
                ),
                tabPanel("About", icon = icon("info-circle"),
                         includeHTML('www/about.html')
                ),
                tabPanel("News", icon = icon("flag"),
                         includeHTML('www/news.html')
                ),
                footer = includeHTML('www/footer.html')
)

# Run the app
# runApp()   # to run locally
shinyApp(ui = ui, server = server,
         options = list(launch.browser = TRUE))  # to run from a server
