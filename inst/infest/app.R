
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


# --------------------------------------------------
# A function to process the insect feeding events from the EPG system

epg <- function(filesIn, time_trim = NULL)
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

   check <- (2*events) %% 2 != 0  # is odd = maybe missing lines
   o <- which(check == TRUE)
   if (any(check == TRUE)) {
      ins[o] <- NULL
      #warning("Error in file(s) ", files[o], ". Please check the line labels.")
      showNotification(paste(" Error (please check line labels) in file(s):",
                          paste(files[o], collapse = ", ")),
                       type = "error", duration = 15)
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
      #warning("File(s) ", files[o2], " seem(s) to have labels sequentialy repeated.")
      showNotification(paste("The following file(s) seem(s) to have labels sequentially repeated :",
                          paste(files[o2], collapse = ", ")),
                       type = "warning", duration = 15)
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
      x[aux %% 2 != 0]
   }
   waves <- lapply(ins, fodd)
   levels <- lapply(waves, unique)
   lev <- unique(unlist(levels))
   waves <- lapply(waves, factor, levels = lev)
   waveevents <- t(sapply(waves, table))   # output and export csv

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

   # aux funct using even and odd indexes at each element of the list
   nonprobe <- c("Z", "z", "np", "NP", "nP", "Np")
   evenodd <- function(x) {
      n <- length(x)
      aux <- seq.int(1, n, by = 1)
      waves <- factor(x[aux %% 2 != 0], levels = lev)
      cumtimes <- as.numeric(x[aux %% 2 == 0])
      ost <- which(waves %in% nonprobe == TRUE)[1]
      FP <- cumtimes[ost]
      times <- c(cumtimes[1], diff(cumtimes))
      c(t_1Pr = FP, tapply(times, waves, FUN = sum))
   }
   times <- t(sapply(ins, evenodd))
   times[is.na(times)] <- 0   # output and export csv

   # Z plus Np
   z_np1 <- colnames(times) %in% nonprobe
   if(sum(z_np1) > 1) {
     times <- as.data.frame(times)
     times$Z_NP <- rowSums(times[, z_np1])
   }
   z_np2 <- colnames(waveevents) %in% nonprobe
   if(sum(z_np2) > 1) {
     waveevents <- as.data.frame(waveevents)
     waveevents$Z_NP <- rowSums(waveevents[, z_np2])
   }

   de <- times[, -1]/waveevents
   de[is.na(de)] <- 0

   # output
   out <- list(duration = times,
               events = waveevents,
               durevents = de,
               seqevents = seqevents, tabseqevents = sumtabseq)
   return(out)
}

# ----------------------------------------------------------
# aux function for finding the best distribution

f_rmse <<- function(model) {
  res = model[["y"]] - fitted(model)
  sqrt(mean(res^2))
}

best_family <<- function(y) {
  fams <- c("NO", "EXP", "GA", "IGAMMA", "IG", "LOGNO",
             "WEI", "PARETO2", "ZAGA", "ZAIG")
  gams <- list()
  for(i in 1:length(fams)) {
    gams[[i]] <- try(gamlss(y ~ 1, sigma.formula = ~1,
                            family = fams[i]), silent = TRUE)
  }
  aics <- sapply(gams, function(w) {
    ifelse(inherits(w, "try-error"), NA, AIC(w))
  })
  o_aic <- which.min(aics)
  out <- list(fit = gams[[o_aic]],
              dfr = gams[[o_aic]]$df.residual,
              family = gams[[o_aic]]$family[2],
              rmse = f_rmse(gams[[o_aic]]))
  return(out)
}

best_family2 <<- function(y) {
  fams <- c("NO", "EXP", "GA", "IGAMMA", "IG", "LOGNO",
             "WEI", "PARETO2", "ZAGA", "ZAIG",
             "PO", "NBII", "ZIP", "ZINBI", "PIG")
  gams <- list()
  for(i in 1:length(fams)) {
    gams[[i]] <- try(gamlss(y ~ 1, sigma.formula = ~1,
                            family = fams[i]), silent = TRUE)
  }
  aics <- sapply(gams, function(w) {
    ifelse(inherits(w, "try-error"), NA, AIC(w))
  })
  o_aic <- which.min(aics)
  out <- list(fit = gams[[o_aic]],
              dfr = gams[[o_aic]]$df.residual,
              family = gams[[o_aic]]$family[2],
              rmse = f_rmse(gams[[o_aic]]))
  return(out)
}

zerofams <<- c("NO", "PARETO2", "ZAGA", "ZAIG",
               "PO", "NBII", "ZIP", "ZINBI", "PIG")


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
     den <- ifelse(input$unit == TRUE, 60, 1)
     df_den <- data.frame(Group = dur()[, 1], round(dur()[, -1]/den, 2))
     updateSelectInput(session, inputId = 'y', choices = names(df_den),
                       selected = "t_1Pr")
     updateSelectInput(session, inputId = 'g', choices = names(df_den))
     rhandsontable(df_den)
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

   output$graph1 <- renderPlot({
     df <- tab()$duration
     if( length(unique(mt()[, 1])) == 1 ) lev <- NA else lev <- factor(mt()[, 1])
     if (nrow(df) < ncol(df)) {
       showNotification("n is too small for a biplot!", type = "error")
     } else {
       acp1 <- prcomp(scale(df))
       fviz_pca_biplot(acp1,
                       geom.var = c("arrow", "text"),
                       geom.ind = "text",
                       col.ind = "black",
                       addEllipses = TRUE, ellipse.level = input$conf1,
                       fill.ind = lev,
                       palette = topo.colors(nlevels(lev)),
                       col.var = "contrib", fill.var = "contrib",
                       legend.title = list(fill = "Group", color = "Contrib (%)"),
                       title = "")
     }
   }, width = 600, height = "auto", res = 80)

   output$graph12 <- renderPlot({
     w <- as.matrix(tab()$duration)
     par(mfrow = c(ceiling(nrow(w)/4), 4), cex = 0.9, mar = c(2, 3, 2, 3))
     for(i in 1:nrow(w)) {
       perc = round(100*w[i, ]/sum(w[i, ]), 1)
       lab = paste0(names(w[i, ]), " (", perc, "%)")
       pie(w[i, ], radius = 0.9,
           col = topo.colors(ncol(w)), main = rownames(w)[i],
           labels = lab, cex.main = 0.9)
     }
   })

   observeEvent(input$best_gamlss1, {
     options(width=100)
     output$best_model1 <- renderTable({
       den <- ifelse(input$unit == TRUE, 60, 1)
       df_den <<- data.frame(Group = mt()[, 1], round(tab()$duration/den, 2))
       Group <- as.factor(df_den[, input$g])
       fits <<- apply(df_den[, -1], 2, best_family)
       mods_gam <- lapply(fits, "[[", "fit")
       #LRTs <- lapply(mods_gam, drop1, data = df_den)
       best_fits <- data.frame(
         Variable = names(mods_gam),
         Distribution = sapply(fits, "[[", "family"),
         DFr = round(sapply(fits, "[[", "dfr"), 1),
         AIC = sapply(mods_gam, "[[", "aic"),
         RMSE = sapply(fits, "[[", "rmse")
         #'LRT(Group)' = sapply(LRTs, "[[", "LRT")[2, ],
         #'p-value(Group)' = sapply(LRTs, "[[", "Pr(Chi)")[2, ]
         )
       best_fits
     }, digits = 4)
   })

   observeEvent(input$fit_glm1, {
     pval_adj <- ifelse(input$tukey == TRUE, "tukey", "none")
     den <- ifelse(input$unit == TRUE, 60, 1)
     df_den <<- data.frame(Group = mt()[, 1], round(tab()$duration/den, 2))
     y <- df_den$y <- df_den[, input$y]
     Group <- df_den$Group <- as.factor(df_den[, input$g])
     if(nlevels(Group) < 2) {
       showNotification("'Group' has only one level.
                      Go back to the WDi data table and specify more levels.",
                                             type = "error", duration = 5)
     } else {
       if( any(y == 0) & !(input$family1 %in% zerofams) ) {
         y <- df_den$y <- y + 0.1
         showNotification("Positive values only are allowed for this
                        distribution model,
                        thus a constant 0.1 was added to all values.",
                          type = "warning", duration = 5)
       }
       assign("familia1", input$family1, envir = .GlobalEnv)
       sigma1 <- ifelse(input$sigma1 == TRUE, expression(~Group), expression(~1))
       mod1 <- try(gamlss(y ~ Group, sigma.formula = eval(sigma1), data = df_den,
                      family = get("familia1", envir = .GlobalEnv),
                      control = gamlss.control(trace = F)))
       mod0 <- try(gamlss(y ~ 1, sigma.formula = eval(sigma1), data = df_den,
                      family = get("familia1", envir = .GlobalEnv),
                      control = gamlss.control(trace = F)))
       med1 <- try(emmeans::emmeans(mod1, "Group", type = "response"))
       comp_med1 <- try(pairs(med1, adjust = pval_adj))
       if(inherits(mod1, "try-error") |
          inherits(mod0, "try-error") |
          inherits(med1, "try-error") |
          inherits(comp_med1, "try-error")) {
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
         output$med1 <- renderPrint({ med1 })
         output$contrasts1 <- renderPrint({ comp_med1 })
         output$graph13 <- renderPlot({ plot(med1, xlab = input$y) +
             theme_bw(base_size = 12) })
       }
     }
   })

   # ----------------------------------------- events
   output$events <-  DT::renderDT({
     df2 <- data.frame(Group = mt()[, 1], tab()$events)
     updateSelectInput(session, inputId = 'y2', choices = names(df2),
                       selected = "t_1Pr")
     updateSelectInput(session, inputId = 'g2', choices = names(df2))
     DT::datatable(df2, editable = F,
                   extensions = 'Buttons',
                   options = list(dom = 'Bfrtip',
                                  digits = 0,
                                  buttons = c('copy', 'csv', 'print'),
                                  text = "Download",
                                  paging = FALSE, searching = TRUE))
   })

   output$graph2 <- renderPlot({
     df2 <- data.frame(Group = mt()[, 1], tab()$events)
     if( length(unique(df2[, 1])) == 1 ) lev <- NA else lev <- factor(df2[, 1])
     if (nrow(df2) < ncol(df2[, -1])) {
       showNotification("n is too small for a biplot!", type = "error")
     } else {
       acp2 <- prcomp(scale(df2[, -1]))
       fviz_pca_biplot(acp2,
                       geom.var = c("arrow", "text"),
                       geom.ind = "text",
                       col.ind = "black",
                       addEllipses = TRUE, ellipse.level = input$conf1,
                       fill.ind = lev,
                       palette = topo.colors(nlevels(lev)),
                       col.var = "contrib", fill.var = "contrib",
                       legend.title = list(fill = "Group", color = "Contrib (%)"),
                       title = "")
     }
   }, width = 600, height = "auto", res = 80)

   output$graph3 <- renderPlot({
      w <- tab()$events
      par(mfrow = c(ceiling(nrow(w)/4), 4), cex = 0.9, mar = c(2, 3, 2, 3))
      for(i in 1:nrow(w)) {
         perc = round(100*w[i, ]/sum(w[i, ]), 1)
         lab = paste0(names(w[i, ]), " (", perc, "%)")
         pie(w[i, ], radius = 0.9,
            col = topo.colors(ncol(w)), main = rownames(w)[i],
            labels = lab, cex.main = 0.9)
      }
   })

   observeEvent(input$best_gamlss2, {
     options(width=100)
     output$best_model2 <- renderTable({
       df_den <<- data.frame(Group = mt()[, 1], tab()$events)
       Group <- as.factor(df_den[, input$g])
       fits <<- apply(df_den[, -1], 2, best_family2)
       mods_gam <- lapply(fits, "[[", "fit")
       best_fits <- data.frame(
         Variable = names(mods_gam),
         Distribution = sapply(fits, "[[", "family"),
         DFr = round(sapply(fits, "[[", "dfr"), 1),
         AIC = sapply(mods_gam, "[[", "aic"),
         RMSE = sapply(fits, "[[", "rmse")
       )
       best_fits
     }, digits = 4)
   })

   observeEvent(input$fit_glm2, {
     pval_adj2 <- ifelse(input$tukey2 == TRUE, "tukey", "none")
     df2 <<- data.frame(Group = mt()[, 1], tab()$events)
     y <- df2$y <- df2[, input$y2]
     Group <- df2$Group <- as.factor(df2[, input$g2])
     if(nlevels(Group) < 2) {
       showNotification("'Group' has only one level.
                      Go back to the WDi data table and specify more levels.",
                        type = "error", duration = 5)
     } else {
       if( any(y == 0) & !(input$family2 %in% zerofams)) {
         y <- df2$y <- y + 1
         showNotification("Positive values only are allowed for this
                        distribution model,
                        thus a constant 1 was added to all values.",
                          type = "warning", duration = 5)
       }
       assign("familia2", input$family2, envir = .GlobalEnv)
       sigma2 <- ifelse(input$sigma2 == TRUE, expression(~Group), expression(~1))
       mod2 <- try(gamlss(y ~ Group, sigma.formula = eval(sigma2), data = df2,
                      family = get("familia2", envir = .GlobalEnv),
                      control = gamlss.control(trace = F)))
       mod0 <- try(gamlss(y ~ 1, sigma.formula = eval(sigma2), data = df2,
                      family = get("familia2", envir = .GlobalEnv),
                      control = gamlss.control(trace = F)))
       med2 <- try(emmeans::emmeans(mod2, "Group", type = "response"))
       comp_med2 <- try(pairs(med2, adjust = pval_adj2))
       if(inherits(mod2, "try-error") |
          inherits(mod0, "try-error") |
          inherits(med2, "try-error") |
          inherits(comp_med2, "try-error")) {
          showNotification("Algorithm has not converged. Try again.",
                          type = "error", duration = 5)
       } else {
         res2 <- mod2[["y"]] - fitted(mod2)
         output$aic2 <- renderPrint({ cat("AIC =", AIC(mod2),
                                          "  RMSE =", sqrt(mean(res2^2)), "\n") })
         output$model2 <- renderPrint({
           l <- logLik(mod2)
           l0 <- logLik(mod0)
           X2 <- 2*c(l - l0)
           grau <- attr(l, "df") - attr(l0, "df")
           pval <- pchisq(X2, grau, lower = F)
           cat("Likelihood Ratio Test for the 'Group' factor\n",
               "LRT =", X2, "   df =", grau, "   Pr(Chisq) =", pval)
         })
         output$med2 <- renderPrint({ med2 })
         output$contrasts2 <- renderPrint({ comp_med2 })
         output$graph23 <- renderPlot({ plot(med2, xlab = input$y2) +
             theme_bw(base_size = 12) })
       }
     }
   })

   # ----------------------------------------- duration per events
   output$durevents <-  DT::renderDT({
     den <- ifelse(input$unit == TRUE, 60, 1)
     df_de <<- data.frame(Group = mt()[, 1], round(tab()$durevents/den, 2))
     updateSelectInput(session, inputId = 'y_de', choices = names(df_de))
     updateSelectInput(session, inputId = 'g_de', choices = names(df_de))
     DT::datatable(df_de, editable = F,
                   extensions = 'Buttons',
                   options = list(dom = 'Bfrtip',
                                  digits = 0,
                                  buttons = c('copy', 'csv', 'print'),
                                  text = "Download",
                                  paging = FALSE, searching = TRUE))
   })

   output$graph_de <- renderPlot({
     df_de <- data.frame(Group = mt()[, 1], tab()$durevents)
     if( length(unique(df_de[, 1])) == 1 ) lev <- NA else lev <- factor(df_de[, 1])
     if (nrow(df_de) < ncol(df_de[, -1])) {
       showNotification("n is too small for a biplot!", type = "error")
     } else {
       acp_de <- prcomp(scale(df_de[, -1]))
       fviz_pca_biplot(acp_de,
                       geom.var = c("arrow", "text"),
                       geom.ind = "text",
                       col.ind = "black",
                       addEllipses = TRUE, ellipse.level = input$conf1,
                       fill.ind = lev,
                       palette = topo.colors(nlevels(lev)),
                       col.var = "contrib", fill.var = "contrib",
                       legend.title = list(fill = "Group", color = "Contrib (%)"),
                       title = "")
     }
   }, width = 600, height = "auto", res = 80)

   output$graph_pie_de <- renderPlot({
     w <- tab()$durevents
     par(mfrow = c(ceiling(nrow(w)/4), 4), cex = 0.9, mar = c(2, 3, 2, 3))
     for(i in 1:nrow(w)) {
       perc = round(100*w[i, ]/sum(w[i, ]), 1)
       lab = paste0(names(w[i, ]), " (", perc, "%)")
       pie(w[i, ], radius = 0.9,
           col = topo.colors(ncol(w)), main = rownames(w)[i],
           labels = lab, cex.main = 0.9)
     }
   })

   observeEvent(input$best_gamlss_de, {
     options(width=100)
     output$best_model_de <- renderTable({
       den <- ifelse(input$unit == TRUE, 60, 1)
       df_de <<- data.frame(Group = mt()[, 1], round(tab()$durevents/den, 2))
       Group <- as.factor(df_de[, input$g])
       fits <<- apply(df_de[, -1], 2, best_family)
       mods_gam <- lapply(fits, "[[", "fit")
       best_fits <- data.frame(
         Variable = names(mods_gam),
         Distribution = sapply(fits, "[[", "family"),
         DFr = round(sapply(fits, "[[", "dfr"), 1),
         AIC = sapply(mods_gam, "[[", "aic"),
         RMSE = sapply(fits, "[[", "rmse")
       )
       best_fits
     }, digits = 4)
   })

   observeEvent(input$fit_glm_de, {
     den <- ifelse(input$unit == TRUE, 60, 1)
     pval_adj_de <- ifelse(input$tukey_de == TRUE, "tukey", "none")
     df_de <<- data.frame(Group = mt()[, 1], round(tab()$durevents/den, 2))
     y <- df_de$y <- df_de[, input$y_de]
     Group <- df_de$Group <- as.factor(df_de[, input$g_de])
     if(nlevels(Group) < 2) {
       showNotification("'Group' has only one level.
                      Go back to the WDi data table and specify more levels.",
                        type = "error", duration = 5)
     } else {
       if( any(y == 0) & !(input$family_de %in% zerofams)) {
         y <- df_de$y <- y + 0.1
         showNotification("Positive values only are allowed for this
                        distribution model,
                        thus a constant 0.1 was added to all values.",
                          type = "warning", duration = 5)
       }
       assign("familia_de", input$family_de, envir = .GlobalEnv)
       sigma_de <- ifelse(input$sigma_de == TRUE, expression(~Group), expression(~1))
       mod_de <- try(gamlss(y ~ Group, sigma.formula = eval(sigma_de), data = df_de,
                      family = get("familia_de", envir = .GlobalEnv),
                      control = gamlss.control(trace = F)))
       mod0 <- try(gamlss(y ~ 1, sigma.formula = eval(sigma_de), data = df_de,
                          family = get("familia_de", envir = .GlobalEnv),
                          control = gamlss.control(trace = F)))
       med_de <- try(emmeans::emmeans(mod_de, "Group", type = "response"))
       comp_med_de <- try(pairs(med_de, adjust = pval_adj_de))
       if(inherits(mod_de, "try-error") |
          inherits(mod0, "try-error") |
          inherits(med_de, "try-error") |
          inherits(comp_med_de, "try-error")) {
         showNotification("Algorithm has not converged. Try again.",
                          type = "error", duration = 5)
       } else {
         res_de <- mod_de[["y"]] - fitted(mod_de)
         output$aic_de <- renderPrint({ cat("AIC =", AIC(mod_de),
                                            "  RMSE =", sqrt(mean(res_de^2)), "\n") })
         output$model_de <- renderPrint({
           l <- logLik(mod_de)
           l0 <- logLik(mod0)
           X2 <- 2*c(l - l0)
           grau <- attr(l, "df") - attr(l0, "df")
           pval <- pchisq(X2, grau, lower = F)
           cat("Likelihood Ratio Test for the 'Group' factor\n",
               "LRT =", X2, "   df =", grau, "   Pr(Chisq) =", pval)
         })
         output$med_de <- renderPrint({ med_de })
         output$contrasts_de <- renderPrint({ comp_med_de })
         output$graph2_de <- renderPlot({ plot(med_de, xlab = input$y_de) +
             theme_bw(base_size = 12) })
       }
     }
   })

   # ---------------------------------------- seq events
   output$seqevents <- DT::renderDT({
     df3 <- data.frame(Group = mt()[, 1], tab()$seqevents)
     updateSelectInput(session, inputId = 'y3', choices = names(df3),
                       selected = "t_1Pr")
     updateSelectInput(session, inputId = 'g3', choices = names(df3))
     DT::datatable(df3, editable = F,
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
}


# ----------------------------------------------------------------
ui = navbarPage(title = tags$head(img(src="infest_1_2.png", height = 65),
                                  "Insect Feeding Behavior Statistics"),
                theme = shinytheme("yeti"),
                windowTitle = "INFEST 1.2",
                useShinyjs(),
                tabPanel("WDi", icon = icon("clock"),
                    actionButton("hideSidebar", "Hide sidebar",
                                 style='padding:4px; font-size:80%;  background-color: lightyellow'),
                    actionButton("showSidebar", "Show sidebar",
                                 style='padding:4px; font-size:80%;  background-color: lightgreen'),
                    sidebarLayout(
                      div( id ="Sidebar",
                      sidebarPanel(width = 3,
                        fileInput("files", "Select one or more EPG files", multiple = TRUE),
                        sliderInput("time_range", "Time trim (min.)",
                                    min = 0, max = 7200, value = c(0, 7200)),
                        checkboxInput("unit", "Convert sec. to min.", value = TRUE),
                        actionButton("run", "Run", icon = icon("r-project"),
                                     style="color: white; background-color: #2e6da4; border-color: white"),
                        tags$h1(),
                        numericInput("conf1", "Confidence level for Group ellipses",
                                     value = 0.95, min = 0, max = 1, step = 0.05,
                                     width = 100)
                      )
                      ),
                      mainPanel(width = 9,
                         tabsetPanel(
                         tabPanel(icon("table"),
                            h5("Waveform Duration by insect"),
                            h6("Warning: only column 'Group' should be edited!"),
                            rHandsontableOutput("duration")),
                         tabPanel(icon("line-chart"),
                             h5("Principal component biplot"),
                             plotOutput("graph1")),
                         tabPanel(icon("chart-pie"),
                                  h5("Proportion of Waveform Duration by insect"),
                                  plotOutput("graph12")),
                         tabPanel(icon("magnifying-glass"),
                                  h5("Find the best fitting distriution for each variable"),
                                  h6("WARNING: this is done using only the intercept as predictor"),
                                  actionButton("best_gamlss1", "Run models", icon = icon("r-project"),
                                               style="color: white; background-color: #2e6da4; border-color: white"),
                                  tableOutput("best_model1")),
                         tabPanel(icon("calculator"),
                                  h5("Fit 1-way GAMLSS model"),
                                  column(width = 3,
                                    selectInput('g', 'Group factor', 'Group'),
                                    selectInput('y', 'Response variable', 't_1Pr'),
                                    selectInput('family1', 'Distribution model',
                                                choices = list("Gaussian" = "NO",
                                                               "Exponential" = "EXP",
                                                               "Gamma" = "GA",
                                                               "Inverse Gamma" = "IGAMMA",
                                                              "Inverse Gaussian" = "IG",
                                                              "Lognormal" = "LOGNO",
                                                              "Weibull" = "WEI",
                                                              "Pareto type 2" = "PARETO2",
                                                              "Zero adjusted Gamma" = "ZAGA",
                                                              "Zero adjusted Inverse Gaussian" = "ZAIG")),
                                    checkboxInput("sigma1", "Model variances of 'Group'", value = FALSE),
                                    checkboxInput("tukey", "Tukey's p-value adjustment", value = TRUE),
                                    actionButton("fit_glm1", "Fit", icon = icon("r-project"),
                                                 style="color: white; background-color: #2e6da4; border-color: white")
                                    ),
                                  column(width = 8,
                                         verbatimTextOutput("aic1"),
                                         verbatimTextOutput("model1"),
                                         verbatimTextOutput("med1"),
                                         verbatimTextOutput("contrasts1"),
                                         plotOutput("graph13"),
                                         tags$h1(" "), tags$h1(" "))
                                  )
                         )
                      )
                  )
                ),
                tabPanel("NWEi", icon = icon("bug"),
                    tabsetPanel(
                      tabPanel(icon("table"),
                          h5("Number of Waveform Events by insect"),
                          DTOutput("events")),
                       tabPanel(icon("line-chart"),
                          h5("Principal component biplot"),
                          plotOutput("graph2")),
                       tabPanel(icon("chart-pie"),
                          h5("Proportion of Waveform Events by insect"),
                          plotOutput("graph3")),
                      tabPanel(icon("magnifying-glass"),
                               h5("Find the best fitting distriution for each variable"),
                               h6("WARNING: this is done using only the intercept as predictor"),
                               actionButton("best_gamlss2", "Run models", icon = icon("r-project"),
                                            style="color: white; background-color: #2e6da4; border-color: white"),
                               tableOutput("best_model2")),
                       tabPanel(icon("calculator"),
                                h5("Fit 1-way GAMLSS model"),
                                column(width = 3,
                                       selectInput('g2', 'Group factor', 'Group'),
                                       selectInput('y2', 'Response variable', 'Z'),
                                       selectInput('family2', 'Distribution model',
                                                   choices = list("Gaussian" = "NO",
                                                                  "Exponential" = "EXP",
                                                                  "Gamma" = "GA",
                                                                  "Inverse Gamma" = "IGAMMA",
                                                                  "Inverse Gaussian" = "IG",
                                                                  "Lognormal" = "LOGNO",
                                                                  "Weibull" = "WEI",
                                                                  "Pareto type 2" = "PARETO2",
                                                                  "Zero adjusted Gamma" = "ZAGA",
                                                                  "Zero adjusted Inverse Gaussian" = "ZAIG",
                                                                  "Poisson" = "PO",
                                                                  "Negative Binomial" = "NBII",
                                                                  "Zero-Inflated Poisson" = "ZIP",
                                                                  "Zero-Inflated Negative Binomial" = "ZINBI",
                                                                  "Poisson-Inverse Gaussian" = "PIG")),
                                       checkboxInput("sigma2", "Model variances of 'Group'", value = FALSE),
                                       checkboxInput("tukey2", "Tukey's p-value adjustment", value = TRUE),
                                       actionButton("fit_glm2", "Fit", icon = icon("r-project"),
                                                    style="color: white; background-color: #2e6da4; border-color: white")
                                       ),
                                column(width = 8,
                                       verbatimTextOutput("aic2"),
                                       verbatimTextOutput("model2"),
                                       verbatimTextOutput("med2"),
                                       verbatimTextOutput("contrasts2"),
                                       plotOutput("graph23"),
                                       tags$h1(" "), tags$h1(" "))
                       )
                    )
                    ),
                tabPanel("WDEi", icon = icon("clock-rotate-left"),
                             tabsetPanel(
                               tabPanel(icon("table"),
                                        h5("Waveform Duration per Event by insect"),
                                        DTOutput("durevents")),
                               tabPanel(icon("line-chart"),
                                        h5("Principal component biplot"),
                                        plotOutput("graph_de")),
                               tabPanel(icon("chart-pie"),
                                        h5("Proportion of Waveform Duration per Event by insect"),
                                        plotOutput("graph_pie_de")),
                               tabPanel(icon("magnifying-glass"),
                                        h5("Find the best fitting distriution for each variable"),
                                        h6("WARNING: this is done using only the intercept as predictor"),
                                        actionButton("best_gamlss_de", "Run models", icon = icon("r-project"),
                                                     style="color: white; background-color: #2e6da4; border-color: white"),
                                        tableOutput("best_model_de")),
                               tabPanel(icon("calculator"),
                                        h5("Fit 1-way GAMLSS model"),
                                        column(width = 3,
                                               selectInput('g_de', 'Group factor', 'Group'),
                                               selectInput('y_de', 'Response variable', 'Z'),
                                               selectInput('family_de', 'Distribution model',
                                                           choices = list("Gaussian" = "NO",
                                                                          "Exponential" = "EXP",
                                                                          "Gamma" = "GA",
                                                                          "Inverse Gamma" = "IGAMMA",
                                                                          "Inverse Gaussian" = "IG",
                                                                          "Lognormal" = "LOGNO",
                                                                          "Weibull" = "WEI",
                                                                          "Pareto type 2" = "PARETO2",
                                                                          "Zero adjusted Gamma" = "ZAGA",
                                                                          "Zero adjusted Inverse Gaussian" = "ZAIG")),
                                               checkboxInput("sigma_de", "Model variances of 'Group'", value = FALSE),
                                               checkboxInput("tukey_de", "Tukey's p-value adjustment", value = TRUE),
                                               actionButton("fit_glm_de", "Fit", icon = icon("r-project"),
                                                            style="color: white; background-color: #2e6da4; border-color: white")
                                               ),
                                        column(width = 8,
                                               verbatimTextOutput("aic_de"),
                                               verbatimTextOutput("model_de"),
                                               verbatimTextOutput("med_de"),
                                               verbatimTextOutput("contrasts_de"),
                                               plotOutput("graph2_de"),
                                               tags$h1(" "), tags$h1(" "))
                               )
                             )
                ),
                tabPanel("Sequential events", icon = icon("sort-alpha-up"),
                         tabsetPanel(tabPanel(icon("table"),
                                  h5("Sequential (from-to) event count"),
                                  DTOutput("seqevents")),
                         tabPanel(icon("project-diagram"),
                                  h5("Network plot"),
                                  plotOutput("graph4"))
                         )
                ),
                tabPanel("Legend", icon = icon("book"),
                         includeHTML('www/legend.html')
                ),
                tabPanel("About", icon = icon("info-circle"),
                         includeHTML('www/about.html')
                ),
                footer = includeHTML('www/footer.html')
)

# Run the app
# runApp()   # to run locally
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))  # to run from a server
