library(shiny)
library(rstan)
library(rgenoud)
library(mclm)
library(tidyverse)
library(survminer)
library(survival)
library(shinyvalidate)
library(DT)
library(shinydashboard)
library(rjson)
library(readxl)

knitr::opts_chunk$set(echo = FALSE)
options(mc.cores = parallel::detectCores())


#TODO:
#celltypes by number instead of name?

#clean up UI (by putting everything in boxes?)
#implement rough AIC/BIC before sampling for each model?
  #use optim instead of genoud

#survival plots
#accept seed/random for cross-validation partitioning

#display crossvalidation status message

ui <- #fluidPage(style = "flex-direction: column",
  dashboardPage(
    dashboardHeader(disable = TRUE),
    dashboardSidebar(disable = TRUE),
    dashboardBody(
                #titlePanel("Bayesian Weibull Regression Survival Analysis"),
                
                fluidPage(
                  fluidRow(style="display:flex;",
                           
                           column(3, wellPanel(
                             radioButtons("fileType", "Select the type of input file",
                                          list(".csv", ".json", ".rds", ".xlsx", ".xls"))
                           )),
                           
                           column(3, wellPanel(
                             uiOutput("fileSelection")
                           ))
                  ),
                  
                  actionButton("fileSelected", "CONTINUE")
                  
                ),
                
                uiOutput("describeData"),
                
                uiOutput("giveCols"),
                
                uiOutput("modelOptions"),
                
                uiOutput("stan"),
                
                verbatimTextOutput("diag"),
                
                uiOutput("crossValidate"),
                
                verbatimTextOutput("crossValidationStatus"),
                
                verbatimTextOutput("CVResults"),
                
                plotOutput("predictionPlot"),
                
                uiOutput("predictions"),
                
                uiOutput("comparisonPlots")
    )
)

server <- function(input, output) {
  
  # vars <- c(1, 2)
  # varsNames <- c("one", "two")
  # names(vars) <- varsNames
  # 
  # output$predictionComparisonPlot <- renderUI(
  #   box(collapsible = TRUE, width = 12,
  #     dashboardPage(
  #       dashboardHeader(title = "Comparison Plots"),
  #       dashboardSidebar(menuItem(radioButtons(inputId = "test", label = "test", choices = vars, selected = vars[2]))),
  #       dashboardBody(
  #         plotOutput("tempPlot")
  #       )
  #       
  #     )
  #   )
  # )
  
  output$tempPlot <- renderPlot(hist(mtcars[, as.numeric(input$test)]))
  # {print(input$test)
  #   print(typeof(input$test))
  # if (input$test == 1) hist(mtcars$mpg) else if (input$test == 2) hist(mtcars$cyl)
  # })
  
  iv <- InputValidator$new()
  
  iv$add_rule("file", sv_required())
  
  
  iv$enable()
  
  
  
  values <- reactiveValues(continuousX = c(), discreteX = c(), T = c(), censoring = c(), 
                           P = 0, N = 0, X = NULL, QR = FALSE, norm = FALSE, fit = NULL,
                           X_norm = NULL, model = NULL, fold = 0)
  
  output$fileSelection <- renderUI(
    fileInput("file", "Input the file", multiple = FALSE, accept = input$fileType)
  )
  
  observeEvent(input$fileSelected, {
    
    req(iv$is_valid())
    
    file <- switch(input$fileType,
                   ".rds" = readRDS(input$file$datapath),
                   ".csv" = read.csv(input$file$datapath),
                   ".json" = fromJSON(file = input$file$datapath),
                   ".xlsx" = read_excel(input$file$datapath),
                   ".xls" = read_excel(input$file$datapath)
                   )
    
    values$file <- na.omit(file)
    
    file <- values$file
    
    names(file) <- sprintf("%s (%d)", names(file), 1:(length(names(file))))
    
    output$dataTable <- renderDataTable({
      datatable(file, 
                options = list(scrollX = TRUE, paging = TRUE, pageLength = 5),
                rownames = FALSE)
    })
    
    output$describeData <- renderUI({
    
      #shiny::validate(need(input$file, "No file given"))
      
      fluidPage(
        
        #fluidRow(
          dashboardPage(
            dashboardHeader(disable = TRUE),
            dashboardSidebar(disable = TRUE),
            dashboardBody(
              #fluidRow(
                box(title = "Data File", width = 20, collapsible = TRUE,
                  dataTableOutput("dataTable")
                )
              #)
            )
          #)
        ),
        
        fluidRow(
          
          column(3, wellPanel(
            sliderInput("numObs", "Select the number of observations:", 
                        min = 1, max = 100, value = 1)
          )),
          
          column(3, wellPanel(
            sliderInput("numContinuousX", "Select the number of continuous explanatory variables:",
                        min = 0, max = 10, value = 0)
          )),
          
          column(3, wellPanel(
            sliderInput("numDiscreteX", "Select the number of discrete explanatory variables:",
                        min = 0, max = 10, value = 0)
          )),
        ),
        
        actionButton("dataDescribed", "CONTINUE")
      )
    
    })
  })
  
  observeEvent(input$dataDescribed, {output$giveCols <- renderUI({
    
    fluidPage(
      fluidRow(
        column(3, wellPanel(
          numericInput("timeCol", "Please assign a column to failure time", value = 0, min = 1, step = 1)
        )),
        
        column(3, wellPanel(
          numericInput("censorCol", "Please assign a column to censoring status", value = 0, min = 1, step = 1)
        )),
        
        column(3, wellPanel(
          selectizeInput("continuousX", "Enter the columns with continuous explanatory variables",
                         choices = NULL, multiple = TRUE, options = list(create = TRUE))
        )),
        
        column(3, wellPanel(
          selectizeInput("discreteX", "Enter the columns with discrete explanatory variables",
                         choices = NULL, multiple = TRUE, options = list(create = TRUE))
        )),
      ),
      
      actionButton("doneWithData", "CONTINUE")
    )
  })
  
  iv$add_rule("timeCol", sv_required())
  iv$add_rule("timeCol", function(value) if (!is.integer(value) || value <= 0) "Must be a positive integer")
  iv$add_rule("censorCol", sv_required())
  iv$add_rule("censorCol", function(value) if (!is.integer(value) || value <= 0) "Must be a positive integer") 
  
  iv$add_rule("continuousX", function(value) if (input$numContinuousX != 0 && is.null(value)) "Required"
              else if (input$numContinuousX != length(value)) "Incorrect number of continuous covariates provided")
  iv$add_rule("discreteX", function(value) if (input$numDiscreteX != 0 && is.null(value)) "Required"
              else if (input$numDiscreteX != length(value)) "Incorrect number of discrete covariates provided")  
  })
  
  
  
  observeEvent(input$doneWithData, output$modelOptions <- renderUI({
    
    req(iv$is_valid())
    
    # #either no continuous variables or columns are given
    # shiny::validate(need(input$numContinuousX == 0 || !is.null(input$continuousX),
    #                      sprintf("No column numbers for continuous variables were provided")))
    # 
    # #correct number of continuous variables given
    # shiny::validate(need(input$numContinuousX == length(input$continuousX),
    #                      sprintf("The number of columns provided for continuous variables (%d) did not match the number of continuous variables expected (%d)",
    #                              length(input$continuousX), input$numContinuousX)))
    
    values$continuousX <- as.numeric(input$continuousX)
    
    
    #either no discrete variables or columns are given
    shiny::validate(need(input$numContinuousX == 0 || !is.null(input$continuousX),
                         sprintf("No column numbers for continuous variables were provided")))
    
    #correct number of discrete variables given
    shiny::validate(need(input$numDiscreteX == length(input$discreteX),
                         sprintf("The number of columns provided for discrete variables (%d) did not match the number of discrete variables expected (%d)",
                                 length(input$discreteX), input$numDiscreteX)))
    
    values$discreteX <- as.numeric(input$discreteX)
    
    #require input for time
    shiny::validate(need(input$timeCol > 0, sprintf("Invalid column number provided for survival time: %d", input$timeCol)))
    
    #require input for censoring
    shiny::validate(need(input$censorCol > 0, sprintf("Invalid column number provided for censoring status: %d", input$censorCol)))
    
    #No duplicate columns used !!!!!!!!!!!!NOT WORKING
    shiny::validate(need(!any(duplicated(c(input$timeCol, input$censorCol, values$discreteX, values$continuousX))),
                         "No duplicate columns allowed"))
    
    req(iv$is_valid())
    
    output$modelOptions <- renderUI(
      fluidPage(
        fluidRow(
          column(3, wellPanel(
            radioButtons("model", "Choose the model",
                         list("Exponentiated Weibull" = "e", "Generalized Weibull"= "g", "Weibull" = "w"))
          )),
          
          column(3, wellPanel(
            checkboxGroupInput("modelOpts", "Select the model options",
                               list("QR Decomp" = "q", "Normalized Data" = "n"))
          )),
          
          column(3, wellPanel(
            numericInput("chains", "Enter the number of MCMC chains", value = 4),
            
            numericInput("iter", "Enter the number of iterations to run each chain for", value = 1000)
          )),
        ),
        
        actionButton("sample", "SAMPLE")
      ))
    
  }))
  
  observeEvent(input$sample, {
    
    data <- values$file #data.matrix(values$file)
    
    # view(data)
    
    values$T <- as.numeric(data[, input$timeCol])
    values$N <- length(values$T)
    
    values$censoring <- as.numeric(data[, input$censorCol]) == 0
    
    #X <- matrix(0, nrow = values$N, ncol = 0)
    X <- data.matrix(data[, values$continuousX])
    #X <- as.matrix(sapply(data[, values$continuousX], as.numeric))
    #X <- data.matrix(data[, values$continuousX])
    
    # view(X)
    
    names <- list()
    
    if (input$numDiscreteX > 0)
      for (i in 1:input$numDiscreteX) {
        index <- values$discreteX[i]
        print(i)
        col <- data[, index]
        print(index)
        for (factor in unique(col)[-1]) {
          print(factor)
          #print(col == factor)
          X <- cbind(X, col == factor)
          names <- append(names, c(names(values$file)[index], factor))
        }
      }
    
    values$names <- append(names(values$file)[values$continuousX], names)
    
    values$X <- X
    values$P <- ncol(X)
    
    if (values$P == 0)
      values$X <- matrix(0, nrow = values$N, ncol = 1)
    
    if (length(input$modelOpts) == 0) {
      values$QR <- FALSE
      values$norm <- FALSE
    }
    else if (length(input$modelOpts) == 1) {
      values$QR <- input$modelOpts == "q"
      values$norm <- !values$QR
    }
    else {
      values$QR <- TRUE
      values$norm <- TRUE
    }
    
    X_norm <- matrix(0, nrow = nrow(X), ncol = ncol(X))
    #Apply scaling recommended by Gelman
    if (values$P != 0) {
      X_norm[, myseq(1, input$numContinuousX)] <- (values$X[, myseq(1, input$numContinuousX)] - mean(values$X[, myseq(1, input$numContinuousX)])) * .5 / sd(values$X[, 1:input$numContinuousX])
      X_norm[, myseq(input$numContinuousX + 1, values$P)] <- values$X[, myseq(input$numContinuousX + 1, values$P)] - mean(values$X[, myseq(input$numContinuousX + 1, values$P)])
    }
    
    values$X_norm <- X_norm
    
    if (values$norm)
      values$X <- values$X_norm
    
    # print(sprintf("values$QR is %s", values$QR))
    # print(sprintf("values$norm is %s", values$norm))
    
    values$model <- stan_model(switch(input$model,
                                      "e" = "exponentiated_regression.stan",
                                      "g" = "general_regression.stan",
                                      "w" = "weibull_regression.stan"))
    
    # print(list(generate_samples = FALSE, bayesian = TRUE,
    #            QR = values$QR, norm = values$norm,
    #            N = values$N, P = values$P, X = values$X,
    #            T = values$T, trun = values$censoring,
    #            N_new = values$N, X_new = values$X))
    
    # for (elt in list(generate_samples = FALSE, bayesian = TRUE,
    #                  QR = values$QR, norm = values$norm,
    #                  N = values$N, P = values$P, X = values$X,
    #                  T = values$T, trun = values$censoring,
    #                  N_new = values$N, X_new = values$X))
    #   print(is.numeric(elt))
    
    
    values$fit <- sampling(values$model, list(generate_samples = FALSE, bayesian = TRUE,
                                              QR = values$QR, norm = values$norm,
                                              N = values$N, P = values$P, 
                                              X = values$X,
                                              T = values$T, trun = values$censoring,
                                              N_new = values$N, X_new = values$X), #Placeholder values for now until generating predictions
                           iter = input$iter, chains = input$chains)
    
    values$X <- X
    
    print(values$fit)
    
    output$stan <- renderUI(
      actionButton("finishedSampling", "Sampling Complete. Diagnose")
    )
    
  })
  
  observeEvent(input$finishedSampling, {
    
    ll <- switch(input$model,
                 "e" = exp_weib_ll,
                 "g" = gen_weib_ll,
                 "w" = weib_ll)
    
    domains <- switch(input$model,
                      "e" = cbind(c(0, 0), c(Inf, Inf)),
                      "g" = cbind(c(0, -Inf), c(Inf, Inf)),
                      "w" = cbind(0, Inf))
    
    if (values$P != 0) {
      start <- get_posterior_mean(values$fit)[c(1, 
                                                (2 + values$P + nrow(domains)):(1 + 2 * values$P + nrow(domains)),
                                                (2 + values$P):(1 + values$P + nrow(domains))), input$chains + 1]
    } else
      start <- get_posterior_mean(values$fit)[c(1, 3:(2 + nrow(domains)))]
    
    # view(start) #for debugging
    # 
    # view(rbind(
    #   cbind(c(0, rep(-Inf, values$P)), rep(Inf, 1 + values$P)),
    #   domains))
    
    
    MLE <- genoud(ll, nvars = length(start), starting.values = start, Domains = rbind(
      cbind(c(0, rep(-Inf, values$P)), rep(Inf, 1 + values$P)),
      domains), optim.method = "Nelder-Mead", max = TRUE
    )
    
    # MLE <- list()
    # MLE$value <- 0
    
    aic <- -2 * MLE$value + 2 * length(start)
    bic <- -2 * MLE$value + length(start) * log(length(values$T))
    
    output$diag <- renderText(sprintf("AIC is %f\nBIC is %f\n\n%s", aic, bic, print_fit()))
    
    if (values$P != 0) {
      
      output$crossValidate <- renderUI(
        fluidPage(
          fluidRow(
            column(3, wellPanel(
              sliderInput("fold", "Select the number of folds for cross-validation", min = 2,
                          max = values$N, value = 5)
            ))
            
          ),
          
          actionButton("crossValidate", label = 
                         sprintf("%s", if (values$fold == 0) "Cross Validate" else if (values$fold > input$fold) "Cross-validation complete" else sprintf("Cross-validation in progress, computing fold %d", values$fold)),
                       )
        )
      )
    }
    
  })
  
  # observeEvent(input$crossValidate, 
  #              output$crossValidationStatus <- renderText(sprintf("Cross-validation in progress, computing fold %d", values$fold))
  #              )
  
  observeEvent(input$crossValidate, {
    
    output$crossValidationStatus <- renderText(sprintf("Cross-validation in progress, computing fold %d", values$fold))
    
    values$CVResults <- cross_validate(input$fold, values$QR, values$norm, values$model, values$T, 
                                       values$X, values$N, values$P, values$censoring, switch(input$model, "e" = 2, "g" = 2, "w" = 1),
                                       num_chains = input$chains, iterations = input$iter)
    
    output$crossValidationStatus <- NULL
    
    output$predictions <- renderUI({
      dashboardPage(
        dashboardHeader(disable = TRUE),
        dashboardSidebar(disable = TRUE),
        dashboardBody(
          box(title = "Predictions", width = 20, collapsible = TRUE,
              dataTableOutput("predictionsTable")
          )
        )
      )
    })
    
    predictionTable <- values$CVResults$samples_T
    colnames(predictionTable) <- sprintf("T_new[%d]", 1:values$N)
    rownames(predictionTable) <- c("Mean", sprintf("Sample %d", 1:(input$iter * input$chains / 2)))
    
    output$predictionsTable <- renderDataTable({
      datatable(as.data.frame(predictionTable), 
                extensions = "Buttons", 
                options = list(scrollX = TRUE, paging = TRUE, buttons = c("copy", "csv", "pdf", "excel"), dom = 'tB', server = FALSE), 
                rownames = TRUE)
    })
    
    output$CVResults <- renderText(sprintf("%s\n\n", analyze_results(values$CVResults)))
    
    survplot_comp <- data.frame(time = c(values$T, values$CVResults$generated_T), status = !c(values$censoring, numeric(values$N)), comp = c(rep("data", values$N), rep("prediction", values$N)))
    
    fit <- survfit(Surv(time, status) ~ comp, data = survplot_comp)
    
    #plot(ggsurvplot(fit, data = survplot_comp))
    
    output$predictionPlot <- renderPlot(ggsurvplot(fit, data = survplot_comp, title = "Comparison of Predicted and Actual Failure Times"))
    
    # output$predictionsBox <- renderUI(
    #   fluidRow(
    #     box(title = "Prediction Values", collapsible = TRUE, width = 12,
    #         renderDataTable(mtcars, options = list(scrollX = TRUE), rownames = FALSE))
    #   ))
    
    if (input$numContinuousX > 0) {
      print("test1")
      varies <- list()
      for (i in myseq(1, values$P - input$numContinuousX))
        varies <- append(varies, sprintf("%s %s", values$names[[2 * i + input$numContinuousX - 1]], values$names[[2 * i + input$numContinuousX]]))
      
      plots1 <- list()
      plots2 <- list()
      
      for (i in myseq(1, values$P - input$numContinuousX)) {
        survplot_cur <- data.frame(time = values$T, status = !values$censoring, category = values$X[, i + input$numContinuousX])
        fit_cur <- survfit(Surv(time, status) ~ category, data = survplot_cur)
        plots1 <- append(plots1, ggsurvplot(fit_cur, data = survplot_cur, title = sprintf("Comparison of Actual Failure Times Based on Whether the Observation was %s", varies[i])))
        survplot_cur <- data.frame(time = values$CVResults$generated_T, status = !numeric(values$N), category = values$X[, i + input$numContinuousX])
        fit_cur <- survfit(Surv(time, status) ~ category, data = survplot_cur)
        plots2 <- append(plots2, ggsurvplot(fit_cur, data = survplot_cur, title = sprintf("Comparison of Predicted Failure Times Based on Whether the Observation was %s", varies[i])))
      }
      
      print("test2")
      
      names(plots1) <- varies
      names(plots2) <- varies
      names(varies) <- varies
      
      print("test3")
      
      output$comparisonPlot <- renderPlot(
        plots1[[input$var]])
      
      print("test4")
      
      output$predictionComparisonPlot <- renderPlot(
        plots2[[input$var]])
      
      print("test5")
      
      #output$test <- renderText(sprintf("%s", input$var))
    
      output$comparisonPlots <- renderUI(

        box(collapsible = TRUE, width = 12, title = "Survival Plots", 
            fluidPage(
              fluidRow(
                radioButtons(inputId = "var", label = "Choose which categorical variable to compare:", choices = varies, selected = varies[1])
              ),
              
              fluidRow(
                plotOutput("comparisonPlot"), plotOutput("predictionComparisonPlot")
              )
            )
            
            # dashboardPage(
            #   dashboardHeader(title = "Survival Plots"),
            #   dashboardSidebar(menuItem(radioButtons(inputId = "var", label = "", choices = varies, selected = varies[1]))),
            #   dashboardBody(
            #     fluidRow(
            #       plotOutput("comparisonPlot"),
            #       plotOutput("predictionComparisonPlot")
            #     )
            #   )
            # )
        )
      )
    }
    
    

  })
  
  exp_weib_ll <- function(params) {
    params <- unlist(params)
    #ll <- numeric(N);
    P <- values$P
    Y <- log(values$T)
    beta0 <- params[1]; betas <- params[myseq(2, P + 1)]; 
    kappa <- params[P + 2]; theta <- params[P + 3];
    
    w <- (Y - beta0 - values$X_norm %*% betas) / kappa
    a <- 1 - exp(-exp(w))
    ll <- ifelse(values$censoring, log(1 - a ^ theta), (theta - 1) * log(a) + w - exp(w))
    sum(ll) + sum(!values$censoring) * log(theta / kappa)
  }
  
  gen_weib_ll <- function(params) {
    params <- unlist(params)
    #ll <- numeric(N);
    P <- values$P
    Y <- log(values$T)
    beta0 <- params[1]; betas <- params[myseq(2, P + 1)]; 
    alpha <- params[values$P + 2]; lambda <- params[values$P + 3];
    
    w <- (Y - beta0 - values$X_norm %*% betas) / alpha 
    ll <- ifelse(values$censoring, 1 / lambda * log(1 - lambda * exp(w)), (1 / lambda - 1) * log(1 - lambda * exp(w)) + w)
    sum(ll) - sum(!values$censoring) * log(alpha)
  }
  
  weib_ll <- function(params) {
    params <- unlist(params)
    #ll <- numeric(N);
    P <- values$P
    Y <- log(values$T)
    beta0 <- params[1]; betas <- params[myseq(2, P + 1)]; 
    b <- params[values$P + 2] 
    
    w <- (Y - beta0 - values$X_norm %*% betas) / b
    ll <- ifelse(values$censoring, -exp(w), w - exp(w))
    sum(ll) - sum(!values$censoring) * log(b)
  }
  
  print_fit <- function() {
    num_model_params <- switch(input$model,
                               "e" = 2,
                               "g" = 2,
                               "w" = 1)
    
    sink(file = "file_for_text_output.txt")
    print(values$fit)
    sink(file = NULL)
    full <- read_txt(file = "file_for_text_output.txt")
    
    P <- values$P
    if (P != 0) {
      keep <- full[c(5:6, 
                     (7 + P + num_model_params):(6 + 2 * P + num_model_params),
                     (6 + P + 1):(6 + P + num_model_params))]
    } else 
      keep <- full[c(5:6, 8:(7 + num_model_params))]
    
    
    write_txt(keep, file = "file_for_text_output.txt")
    
    if (P != 0) {
      sink(file = "file_for_text_output.txt", append = TRUE)
      
      for (i in myseq(1, input$numContinuousX)) 
        cat(sprintf("betas[%d] models the effect of %s\n", i, values$names[[i]]))
      
      for (i in myseq(1, (P - input$numContinuousX)))
        cat(sprintf("betas[%d] models the effect of %s %s\n", i + input$numContinuousX, values$names[[2 * i + input$numContinuousX - 1]], values$names[[2 * i + input$numContinuousX]]))
      
      sink(file = NULL)
    }
    
    read_file("file_for_text_output.txt")
  }
  
  analyze_results <- function(CVResults) {
    generated_T <- CVResults$generated_T[!values$censoring]
    T <- values$T[!values$censoring]
    sink("file_for_text_output.txt")
    SSE <- sum((T - generated_T) ^ 2)
    cat("Note: these estimates only compare predictions to non-censored data\n")
    cat(sprintf("The sum of squares due to error was %f\n\n", signif(SSE, 4)))
    cat(sprintf("The average squared error was %f\n\n", signif(SSE / (values$N - values$P), 4)))
    ASE <- sum(abs(T - generated_T)) / values$N
    cat(sprintf("The average absolute error was %f\n\n", signif(ASE, 4)))
    sink(file = NULL)
    
    read_file("file_for_text_output.txt")
  }
  
  cross_validate <- function(K, QR, norm, model, T, X, N, P, censoring, num_model_parameters
                             = 2, num_chains = 8, iterations = 5000, init_r = 2) {
    
    generated_Y <- numeric(N);
    generated_T <- numeric(N);
    samples_T <- matrix(0, ncol = N, nrow = 1 + num_chains * iterations / 2)
    fits <- list()
    
    #############################################################################
    # Partition the Data
    #############################################################################
    set.seed(123)
    num_parts <- K
    size_part <- trunc(1 / num_parts * N)
    indices <- sample(1:N, N, replace = FALSE)
    parts <- list()
    for (i in 1:(num_parts - 1)) {
      parts[[i]] <- indices[seq(size_part * (i - 1), size_part * i - 1)]
    }
    parts[[num_parts]] <- indices[seq(size_part * (num_parts - 1), N)]
    
    set.seed(Sys.time())
    
    #############################################################################
    # Run the Model With Each Partition as a Testing Set
    #############################################################################
    for (i in 1:num_parts) {
      values$fold <- i
      print(sprintf("################################ Part %d ################################", i))
      training_parts <- (1:num_parts)[-c(i)]
      testing_indices <- parts[[i]]
      training_indices <- unlist(parts[training_parts])
      X_train <- as.matrix(X[training_indices, ])
      X_test <- as.matrix(X[testing_indices, ])
      Ti <- as.numeric(T[training_indices])
      N_train <- nrow(X_train)
      print(sprintf("N_train is %d", N_train))
      N_test <- nrow(X_test)
      censoring_train <- censoring[training_indices]
      
      fit <- sampling(model, list(generate_samples = TRUE, bayesian = TRUE, QR = QR, norm = norm, N = N_train, X = X_train, T = Ti, P = P, N_new = N_test, X_new = X_test, trun = censoring_train), iter = iterations, chains = num_chains, init_r = init_r)
      
      fits[[i]] <- fit
      print(fit)
      means <- get_posterior_mean(fit, pars = "T_new")
      samples <- rstan::extract(fit, pars = "T_new")
      # view(means)
      # view(samples)
      
      generated_T[testing_indices] <- means[, num_chains + 1]
      # view(generated_T[testing_indices])
      # view(unlist(samples))
      # view(rbind(generated_T[testing_indices], matrix(unlist(samples), byrow = FALSE, ncol = N_test)))
      samples_T[, testing_indices] <- rbind(generated_T[testing_indices], matrix(unlist(samples), byrow = FALSE, ncol = N_test))
      
      # for (j in 1:N_test) {
      #   #generated_Y[testing_indices[j]] <- means[P + num_model_parameters + j, num_chains + 1]
      #   generated_T[testing_indices[j]] <- means[1 + 2 * P + num_model_parameters + N_test + j, num_chains + 1]
      #   samples_T[testing_indices[j], ] <- c(generated_T[testing_indices[j]],
      #                                        samples[[1 + 2 * P + num_model_parameters + N_test + j]])
      # }
    }
    values$fold <- values$fold + 1
    return(list(fits = fits, generated_T = generated_T, samples_T = samples_T))
  }
  
  myseq <- function(m, n) {
    if (m > n)
      return(integer())
    return(seq(m, n))
  }
}

shinyApp(ui = ui, server = server)