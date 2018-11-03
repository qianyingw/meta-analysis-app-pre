# This is a Shiny web application for meta-analysis
# by Qianying Wang @CAMARADES
# server.R
# https://camarades.shinyapps.io/meta-analysis/


library(shiny)
library(metafor)
library(meta)
library(shinythemes)
library(plotly)
library(ggplot2)
library(colourpicker)
library(dplyr)
library(RCurl)



shinyServer(function(input, output, session) {
  
  # ---------- Read csv file ------------------
  dataset <- reactive({
    inFile = input$file
    if (is.null(inFile)) { 
      stop("No file uploaded.")
    }
    dataset <- read.csv(inFile$datapath, header = T, row.names = NULL, stringsAsFactors = F)
    # dataset <- na.omit(dataset)
  })
  
  
  # ------------------ Calculate yi, vi --------------------
  # yi: effect size
  # vi: variance of effect size = SE^2
  yivi <- reactive({
    
    if (is.null(dataset()))
      stop("No file uploaded.")
    else
      dat <- dataset()
    # dat<- read.csv("C:/Users/wqy/Desktop/DevelopmentalAllSeparate_MODEL.csv",stringsAsFactors=F,header=T)
    
    
    
    # --------------------------- GEN ------------------------------------- 
    # If the csv file has included effect size and SE 
    # (does not need app to calculate effect size and SE)
    if (input$EffectMeasure == "GEN") {
      yi = dat[,"Effect.size"]
      vi = dat[,"SE"]^2
      ni = dat[,"N"]  # No.True.C + No.T
      yivi = cbind(dat, yi, vi, ni)
    } 
    
    # --------------------------- NMD -------------------------------------     
    if (input$EffectMeasure == "NMD") {
      # Sham 
      m.s = dat[,"Mean.in.Sham"]
      n.s = dat[,"Number.in.Sham"]
      A = dat[,"Standard.Deviation.in.Sham"]
      B = dat[,"Standard.Error.in.Sham"] * sqrt(n.s)
      X = A; X[is.na(A)] = B[is.na(A)]; X[is.na(B)] = A[is.na(B)]
      sd.s = X
      # Control
      m.c = dat[,"Mean.in.Control.Group"]
      n.c = dat[,"Number.in.Control.Group"]   
      A = dat[,"Standard.Deviation..C."]
      B = dat[,"Standard.Error..C."] * sqrt(n.c)
      X = A; X[is.na(A)] = B[is.na(A)]; X[is.na(B)] = A[is.na(B)]
      sd.c = X
      # Treatment
      m.t = dat[,"Mean.in.Treatment.Group"]
      n.t = dat[,"Number.in.Treatment.Group"]
      A = dat[,"Standard.Deviation..Rx."]
      B = dat[,"Standard.Error..Rx."] * sqrt(n.t)
      X = A; X[is.na(A)] = B[is.na(A)]; X[is.na(B)] = A[is.na(B)]
      sd.t = X
      
      if (input$CompareType == "SC") { # sham + control
        if (sum(names(dat) %in% "Control.Groups.per.Sham") == 0){
          n.s.true = rep(NA, nrow(dat))
        } else {
          n.s.true = n.s/dat[,"Control.Groups.per.Sham"]
        }
        yi = 100*(m.s-m.c)/m.s
        vi = (100*(sd.s/m.s))^2/n.s.true + (100*(sd.c/m.s))^2/n.c
        ni = n.s.true + n.c
      } 
      
      if (input$CompareType == "SCT"){ # sham + control + treatment
        n.c.true = n.c/dat[,"Treatment.Groups.per.Control"]
        # n.c.true = dat[,"True.No.of.C"]
        m.s[is.na(m.s)] = 0
        yi = 100*(m.c-m.t)/(m.c-m.s)
        vi = (100*(sd.c/(m.c-m.s)))^2/n.c.true + (100*(sd.t/(m.c-m.s)))^2/n.t
        ni = n.c.true + n.t
        
      }
      
      yivi = cbind(dat, yi, vi, ni)
    } # NMD
    
    # --------------------------- SMD -------------------------------------   
    if (input$EffectMeasure == "SMD") {
      
      dat[dat[,"Higher.score.is"]=="better", "Higher.score.is"]=-1
      dat[dat[,"Higher.score.is"]=="worse", "Higher.score.is"]=1
      direction=as.numeric(dat[,"Higher.score.is"])
      
      
      # Sham 
      m.s = dat[,"Mean.in.Sham"]
      n.s = dat[,"Number.in.Sham"]
      A = dat[,"Standard.Deviation.in.Sham"]
      B = dat[,"Standard.Error.in.Sham"] * sqrt(n.s)
      X = A; X[is.na(A)] = B[is.na(A)]; X[is.na(B)] = A[is.na(B)]
      sd.s = X
      # Control
      m.c = dat[,"Mean.in.Control.Group"]
      n.c = dat[,"Number.in.Control.Group"]   
      A = dat[,"Standard.Deviation..C."]
      B = dat[,"Standard.Error..C."] * sqrt(n.c)
      X = A; X[is.na(A)] = B[is.na(A)]; X[is.na(B)] = A[is.na(B)]
      sd.c = X
      # Treatment
      m.t = dat[,"Mean.in.Treatment.Group"]
      n.t = dat[,"Number.in.Treatment.Group"]
      A = dat[,"Standard.Deviation..Rx."]
      B = dat[,"Standard.Error..Rx."] * sqrt(n.t)
      X = A; X[is.na(A)] = B[is.na(A)]; X[is.na(B)] = A[is.na(B)]
      sd.t = X
      
      if (input$CompareType == "SC") { # sham + control
        if (sum(names(dat) %in% "Control.Groups.per.Sham") == 0){
          n.s.true = rep(NA, nrow(dat))
        } else {
          n.s.true = n.s/dat[,"Control.Groups.per.Sham"]
        }
        Sp = sqrt((n.s.true-1)*sd.s^2/(n.s.true+n.c-2) + (n.c-1)*sd.c^2/(n.s.true+n.c-2))
        yi = (m.s-m.c)/Sp*(1-3/(4*(n.s.true+n.c)-9))*direction
        vi = (n.s.true+n.c)/(n.s.true*n.c) + yi^2/2/(n.s.true+n.c-3.94)
        ni = n.s.true + n.c
      } 
      
      if (input$CompareType == "SCT"){ # sham + control + treatment
        n.c.true = n.c/dat[,"Treatment.Groups.per.Control"] 
        Sp = sqrt((n.c.true-1)*sd.c^2/(n.c.true+n.t-2) + (n.t-1)*sd.t^2/(n.c.true+n.t-2))
        yi = (m.c-m.t)/Sp*(1-3/(4*(n.c.true+n.t)-9))*direction
        vi = (n.c.true+n.t)/(n.c.true*n.t) + yi^2/2/(n.c.true+n.t-3.94)
        ni = n.c.true + n.t
      }
      
      yivi = cbind(dat, yi, vi, ni)
    } # SMD
    
    
    
    # --------------------------- Odds Ratio -------------------------------------         
    if (input$EffectMeasure == "OR") {
      
      n.s = dat[,"Number.in.Sham"]
      n.c = dat[,"Number.in.Control.Group"]   
      n.t = dat[,"Number.in.Treatment.Group"]
      
      if (input$CompareType == "SC") { # sham + control
        
        num.c = dat[,"Number.Affectd.by.Outcome.Measure.in.Sham.Group"]
        n.c = dat[,"Number.in.Sham"]
        
        num.t = dat[,"Number.Affectd.by.Outcome.Measure.in.Control.Group"]
        n.t = dat[,"Number.in.Control.Group"]
        
        ai = num.t
        bi = n.t - num.t
        ci = rep(NA, nrow(dat))
        di = rep(NA, nrow(dat))
        
        if (sum(names(dat) %in% "Number.Affectd.by.Outcome.Measure.in.Sham.Group") != 0){
          ci = num.c
          di = n.c - num.c
        }  
        if (sum(names(dat) %in% "Control.Groups.per.Sham") == 0){
          n.s.true = rep(NA, nrow(dat))
        } else {
          n.s.true = n.s/dat[,"Control.Groups.per.Sham"]
        }
        ni = n.s.true + n.c
      }
      
      if (input$CompareType == "SCT") { # sham + control + treatment
        
        num.c = dat[,"Number.Affectd.by.Outcome.Measure.in.Control.Group"]
        n.c = dat[,"Number.in.Control.Group"]
        
        num.t = dat[,"Number.Affectd.by.Outcome.Measure.in.Treatment.Group"]
        n.t = dat[,"Number.in.Treatment.Group"]
        
        ai = num.t
        bi = n.t - num.t
        ci = num.c
        di = n.c - num.c
        
        ni = n.c.true + n.t
      }
      
      
      # When any of ai, bi, ci and di is 0,
      # 0.5 is added to avoid computation problems
      for (j in 1:length(ai)) {
        if (ai[j] == 0 | bi[j] == 0 | ci[j] == 0 | di[j] == 0) {
          ai[j] = ai[j] + 0.5
          bi[j] = bi[j] + 0.5
          ci[j] = ci[j] + 0.5
          di[j] = di[j] + 0.5
        }
      }
      
      res <-  metabin(event.e = ai, n.e = ai + bi, 
                      event.c = ci, n.c = ci + di,
                      sm = "OR", method = "MH", method.tau = "DL")
      
      lnORi = res$TE     # log(ai*di/(bi*ci))
      SElnORi = res$seTE # sqrt(1/ai + 1/bi + 1/ci + 1/di)
      
      yivi = cbind(dat, ni, ai, bi, ci, di, lnORi, SElnORi)
    } # OR
    
    return(yivi)
  })
  
  
  # ----------------- Data Table ----------------------
  # Display data table
  output$DT <- renderDataTable({
    if (input$DataType == "n0"){ DT <- yivi() }
    if (input$DataType == "n1"){ DT <- Nest() }
    return (DT)
  }, options = list(pageLength = 10))
  
  # Download data table
  output$DownTable <- downloadHandler(
    filename = function() {
      paste("Table", input$FileType, sep = ".")
    },
    content = function(file) {
      if (input$DataType == "n0"){ DT <- yivi() }
      if (input$DataType == "n1"){ DT <- Nest() }
      sep <- switch(input$FileType, "csv" = ",", "tsv" = "\t")
      write.table(DT, file, sep = sep, row.names = F)
    }
  )
  
  
  
  
  # ----------------- Meta-analysis ----------------------
  output$GlobalOutput <- renderPrint({
    
    yv <- yivi() 
    if (is.null(yv)) { stop("No file uploaded.") }
    
    if (input$EffectMeasure == "OR") {
      res = summary(metabin(event.e = ai, n.e = ai + bi,
                            event.c = ci, n.c = ci + di,
                            data = yv,
                            sm = "OR", method = "MH",
                            method.tau = input$HetEstimator))
      
    } else {
      res = summary(metagen(TE = yi, 
                            seTE = sqrt(vi), 
                            data = yv, 
                            comb.fixed = F,
                            method.tau = input$HetEstimator, 
                            hakn = T))
    }
    
    res
    
  })
  
  # ----------------- Forest plot ----------------------
  Fplot <- function(){
    
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop("No file uploaded.") }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop("Cannot calculate effect size. Please change model or method.") } 
    }
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop("Your csv file doesn't contain Pub ID or Group letter. Studies can't be nested.")
      }
    }
    
    # meta-analysis
    if (input$EffectMeasure == "OR") {
      res = metabin(event.e = ai, n.e = ai + bi, 
                    event.c = ci, n.c = ci + di,
                    data = yv,
                    sm = "OR", method = "MH", 
                    comb.fixed = F,
                    method.tau = input$HetEstimator)
      
    } else {
      res = metagen(TE = yi, 
                    seTE = sqrt(vi), 
                    data = yv, 
                    comb.fixed = F,
                    method.tau = input$HetEstimator, 
                    hakn = T)
    }
    
    # Order of individual studies displayed on forest plot 
    if (input$ForOrder == "") { sortlab = NULL }             # original order
    if (input$ForOrder == "ine") { sortlab = res$TE }        # incresing effect size
    if (input$ForOrder == "inw") { sortlab = res$w.random }  # incresting weight
    if (input$ForOrder == "iny") { sortlab = yv$Year }       # increasing year
    if (input$ForOrder == "dee") { sortlab = -res$TE }       # decreasing effect size
    if (input$ForOrder == "dew") { sortlab = -res$w.random } # decreasing weight
    if (input$ForOrder == "dey") { sortlab = -yv$Year }      # decreasing year
    
    
    # whether plot weight on the right side of the forest plot
    if (input$ShowWeight == T) {
      right_cols = c("effect", "ci", "w.random")
      right_labs = c("Effect", "95% CI", "Weight")
    } else {
      right_cols = c("effect", "ci")
      right_labs = c("Effect", "95% CI")
    }
    
    
    # x-axis limit
    if (input$ForXmin == "") {
      x_min = min(res$lower, na.rm = T)
      if (input$EffectMeasure == "OR") {
        x_min = min(exp(res$lower), na.rm = T)
      }
    } else {
      x_min = as.numeric(input$ForXmin)
    }
    if (input$ForXmax == "") {
      x_max = max(res$upper, na.rm = T)
      if (input$EffectMeasure == "OR") {
        x_max = max(exp(res$upper), na.rm = T)
      }
    } else {
      x_max = as.numeric(input$ForXmax)
    }
    
    
    forest(res, # metagen object
           hetstat = F, # whether print results for heterogeneity measures
           
           studlab = paste0(yv$Surname, ", ", yv$Year),
           leftcols = c("studlab"), 
           rightcols = right_cols, 
           rightlabs = right_labs,
           
           sortvar = sortlab, # Year, w.random, TE
           
           xlim = c(x_min, x_max),
           xlab = input$ForXlab,
           
           col.square = input$ForSqCol, 
           col.diamond = input$ForDiaCol, 
           
           plotwidth = paste0(input$ForWidth, "cm"), # x-axis width,
           colgap.forest.left = paste0(input$GapLeft, "cm"),
           colgap.forest.right = paste0(input$GapRight, "cm"),
           
           col.square.lines = "black", 
           col.diamond.lines = "black", 
           col.i = "black",
           
           digits = 3,
           digits.se = 3,
           digits.weight = 3
    )
  }
  
  # Width and height of plot window
  FWinW <- reactive({ input$ForWinWidth })
  FWinH <- reactive({ input$ForWinHeight })
  
  # Show forest plot
  output$ForestPlot <- renderPlot({
    Fplot()
  }, width = FWinW, height = FWinH)
  
  # Download forest plot
  output$DownForest <- downloadHandler(
    filename = function() {
      paste("forest-", Sys.time(), ".png", sep = "")
    }, 
    content = function(file) {
      png(file, width = FWinW(), height = FWinW()) 
      Fplot()
      dev.off() 
    }
  )
  
  # ----------------- Select variables for Het analysis ----------------------
  # for stratified (one discrete variable)
  output$SubVar <- renderUI({
    
    yv <- dataset()
    if (is.null(yv)) { stop() }
    if (sum(is.na(yv$yi)) == nrow(yv)) { stop() } 
    
    use <- yv[, -which(names(yv) %in% c("Pub.ID","Year"))]
    
    selectInput(inputId = "subvar", 
                label = "Select discrete variable for stratified meta-analysis",
                choices = as.list(c(Choose='', names(use))), 
                multiple=F, selectize=F)
  })
  
  
  # for meta-regression (multi continuous variables)
  output$RegContVar <- renderUI({
    
    yv <- dataset()
    if (is.null(yv)) { stop() }
    if (sum(is.na(yv$yi)) == nrow(yv)) { stop() } 
    
    use <- yv[, -which(names(yv) %in% c("Pub.ID","Year"))]
    selectInput(inputId = "cbox", 
                label = "Select continuous variables for meta-regression",
                choices = as.list(names(use)), 
                multiple=T, selectize=F) 
  })
  
  # for meta-regression (multi discrete variables)
  output$RegDiscVar <- renderUI({
    
    yv <- dataset()
    if (is.null(yv)) { stop() }
    if (sum(is.na(yv$yi)) == nrow(yv)) { stop() } 
    
    use <- yv[, -which(names(yv) %in% c("Pub.ID","Year"))]
    selectInput(inputId = "dbox", 
                label = "Select discrete variables for meta-regression",
                choices = as.list(names(use)), 
                multiple=T, selectize=F)
  })
  
  
  
  # ---------------------- Het (stratified meta-analysis) ------------------------
  # Summary table of global meta-analysis shown on "Het" tab
  output$AllOutput <- renderPrint({
    
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop() }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop() } 
    }
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop()
      }
    }
    
    if(length(input$subvar) == 0) { stop() }
    
    # meta-analysis summary information
    res = summary(metagen(TE = yi, 
                          seTE = sqrt(vi), 
                          data = yv, 
                          comb.fixed = F,
                          method.tau = input$HetEstimator,
                          hakn = T))
    
    
    
    
    
    # create summary data frame of global meta-analysis 
    df = data.frame(effect = round(res$random$TE, 3),
                    se = round(res$random$seTE, 3),
                    
                    # confidence interval
                    lower = round(res$random$lower, 3),
                    upper = round(res$random$upper, 3),
                    
                    # heterogeneity statistic
                    Q = round(res$Q, 3),
                    tau2 = round((res$tau)^2, 3),
                    I2 = paste0(round(res$I2$TE*100,3),"%"),
                    
                    pvalue = res$random$p,
                    No.comparisons = round(res$k, 3),
                    No.animals = round(sum(yv$ni), 3),
                    No.pub = length(unique(yv$Pub.ID)), # number of publications
                    row.names = "study"
    )
    
    print(df, print.gap=3)
    
  }) 
  
  
  # Summary table of heterogeneity test shown on "Het" tab
  output$SubTest <- renderPrint({
    
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop("No file uploaded.") }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop("Cannot calculate effect size. Please change model or method.") } 
    }
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop("Your csv file doesn't contain Pub ID or Group letter. Studies can't be nested.")
      }
    }
    
    if(length(input$subvar)==0) { stop("No discrete variable selected") }
    
    # stratified meta-analysis
    
    
    # res = metabin(event.e = ai, n.e = ai + bi, 
    #               event.c = ci, n.c = ci + di,
    #               data = yv,
    #               sm = "OR", method = "MH", 
    #               comb.fixed = F,
    #               method.tau = input$HetEstimator)
    
    # res = metabin(event.e = ai, n.e = ai + bi,
    #               event.c = ci, n.c = ci + di,
    #               data = yv,
    #               byvar = yv[,"Random.Allocation.to.Group"],
    #               sm = "OR", method = "MH",
    #               comb.fixed = F,
    #               method.tau = "DL")
    
    res=summary(metagen(TE = yi, 
                        seTE = sqrt(vi), 
                        byvar = yv[,input$subvar], # grouping variable for stratified meta-analysis
                        data = yv, 
                        comb.fixed = F,
                        method.tau = input$HetEstimator,
                        hakn = T)) 
    
    df2 = data.frame(
      # residual heterogeneity Q
      Q = round(res$Q - sum(res$Q.w),3), 
      # degree of freedom
      df = res$df.Q.b,
      # chi-square statistic with significance level 0.05
      chi2 = round(qchisq(0.05, df=res$df.Q-res$df.Q.w, ncp=0, lower.tail=F, log.p=F), 3),
      pvalue = round(pchisq(res$Q - sum(res$Q.w), df=res$df.Q-res$df.Q.w, ncp = 0, lower.tail = F, log.p = F), 8),
      row.names = "Between groups")
    print(df2, print.gap=3)
  })
  
  
  
  
  # Summary table of stratified meta-analysis shown on "Het" tab
  output$SubOutput <- renderPrint({
    
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop() }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop() } 
    }
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop()
      }
    }
    
    if(length(input$subvar) == 0) { stop() }
    
    # stratified meta-analysis
    res = summary(metagen(TE = yi, 
                          seTE = sqrt(vi), 
                          byvar = yv[,input$subvar],
                          data = yv, 
                          comb.fixed = F,
                          method.tau = input$HetEstimator,
                          hakn = T))
    
    
    mod <- metagen(TE = yi, 
                   seTE = sqrt(vi), 
                   byvar = yv[,input$subvar],
                   data = yv, 
                   comb.fixed = F,
                   method.tau = input$HetEstimator,
                   hakn = T)
    
    # calculate number of animals in each group
    tem = aggregate(yv$ni, by = list(gro=yv[,input$subvar]), FUN = sum)
    temp = tem
    tlab = as.matrix(mod$bylevs)
    
    
    
    # to gurantee that variables in stratified meta-analysis output have same order
    # to variables in temp
    
    for (t in 1:nrow(temp)) {
      temp[t,] = tem[tem[,1] == tlab[t],]
    }
    
    
    
    df = data.frame(effect = round(res$within.random$TE, 3),
                    se = round(res$within.random$seTE, 3),
                    
                    # confidence interval
                    lower = round(res$within.random$lower, 3),
                    upper = round(res$within.random$upper, 3),
                    
                    # heterogeneity statistic
                    Q = round(res$Q.w, 3),
                    tau2 = round((res$tau.w)^2, 3),
                    I2 = paste0(round(res$I2.w$TE*100,3),"%"),
                    
                    pvalue = res$within.random$p,
                    No.comparisons = round(res$k.w, 3),
                    No.animals = round(temp[,2], 3),
                    row.names = res$bylevs)
    
    print(df, print.gap=3)
    
  }) 
  
  # ---------------------- Het (meta-regression) ------------------------
  # summary info of meta-regression
  output$RegOutput <- renderPrint({
    
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop("No file uploaded.") }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop("Cannot calculate effect size. Please change model or method.") } 
    }
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop("Your csv file doesn't contain Pub ID or Group letter. Studies can't be nested.")
      }
    }
    
    dlab = input$dbox   # discrete variables selected
    clab = input$cbox   # continuous variables selected
    
    
    # regression formula for different type of variables
    if(length(clab)==0 & length(dlab)==0) { stop("No variable selected") }
    if(length(clab)==0 & length(dlab)!=0) { 
      # only discrete variables in the regression model
      # mo = ~ factor(d1) + factor(d2) + ...
      mo = paste("~",paste(paste("factor(", paste(dlab, collapse=")+factor(")),")"),collapse="")
    }
    if(length(clab)!=0 & length(dlab)==0){ 
      # only continuous variables in the regression model
      # mo = ~ c1 + c2 + ...
      mo = paste("~",paste(clab, collapse="+"), collapse="")
    }
    if(length(clab)!=0 & length(dlab)!=0){ 
      # both discrete and continuous variables exist in the regression model
      # mo = ~ factor(d1) + factor(d2) + ... + c1 + c2 + ...
      form.d = paste(paste("+factor(", paste(dlab, collapse=")+factor(")),")")
      form.c = paste("~",paste(clab, collapse="+"), collapse="")
      mo = paste(form.c, form.d, collapse="+")
    }
    
    # meta-regression
    rma(yi = yi, 
        vi = vi, 
        mods = as.formula(mo), 
        data = yv, 
        method = input$HetEstimator, test="knha")
  })
  
  
  # ------------------- Select variables for Het plot --------------------
  # variables for sub forest plot (one discrete variable)
  output$SubPlotVar <- renderUI({
    
    yv <- dataset()
    if (is.null(yv)) { stop() }
    if (sum(is.na(yv$yi)) == nrow(yv)) { stop() } 
    
    use <- yv[, -which(names(yv) %in% c("Pub.ID","Year"))]
    
    selectInput(inputId = "dis.forest", 
                label = "Select variables for forest plot",
                choices = as.list(c(Choose='', names(use))), 
                multiple=F, selectize=F)
  })
  
  
  # variables for meta-regression plot
  output$RegPlotVar <- renderUI({
    
    yv <- dataset()
    if (is.null(yv)) { stop() }
    if (sum(is.na(yv$yi)) == nrow(yv)) { stop() } 
    
    use <- yv[, -which(names(yv) %in% c("Pub.ID","Year"))]
    
    selectInput(inputId = "con.regplot", 
                label = "Select one variable for meta-regression",
                choices = as.list(c(Choose='', names(use))), 
                selectize=F)
  })
  
  
  # ---------------------- Het plot (subforest plot) -----------------------
  Fplot2 <- function(){
    
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop("No file uploaded.") }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop("Cannot calculate effect size. Please change model or method.") } 
    }
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop("Your csv file doesn't contain Pub ID or Group letter. Studies can't be nested.")
      }
    }
    
    if (length(input$dis.forest)==0) { stop("No discrete variable selected.") }
    
    # stratified meta-analysis
    res <- metagen(TE = yi, 
                   seTE = sqrt(vi), 
                   byvar = yv[,input$dis.forest],
                   data = yv, 
                   method.tau = input$HetEstimator,
                   hakn = T)
    
    # columns to be plotted on the right side of the forest plot
    # show weight of fixed-effect model or not
    # show weight of random-effects model or not
    if (input$ShowFixWeight == T & input$ShowRandWeight == T) {
      right_cols = c("effect", "ci", "w.fixed", "w.random")
      right_labs = c("Effect", "95% CI", "W-fixed", "W-random")
    }
    if (input$ShowFixWeight == T & input$ShowRandWeight == F) {
      right_cols = c("effect", "ci", "w.fixed")
      right_labs = c("Effect", "95% CI", "W-fixed")
    }
    if (input$ShowFixWeight == F & input$ShowRandWeight == T) {
      right_cols = c("effect", "ci", "w.random")
      right_labs = c("Effect", "95% CI", "W-random")
    }
    if (input$ShowFixWeight == F & input$ShowRandWeight == F) {
      right_cols = c("effect", "ci")
      right_labs = c("Effect", "95% CI")
    }
    
    
    forest(res, 
           hetstat = F, 
           bylab = input$dis.forest, # variable selected for sub forest
           print.subgroup.labels = T,
           
           studlab = paste0(yv$Surname, ", ", yv$Year),
           leftcols = c("studlab"), 
           rightcols = right_cols,
           rightlabs = right_labs,
           
           xlab = input$SubXlab,
           
           col.square = input$SubSqCol, 
           col.diamond = input$SubDiaCol,
           
           plotwidth = paste0(input$SubWidth, "cm"),
           colgap.forest.left = paste0(input$SubGapLeft, "cm"),
           colgap.forest.right = paste0(input$SubGapRight, "cm"),
           
           col.square.lines = "black", 
           col.diamond.lines = "black", 
           col.i = "black"
    )
  }
  
  # width and height of plot window
  SWinW<- reactive({ input$SubWinWidth })
  SWinH <- reactive({ input$SubWinHeight })
  
  # display sub forest plot
  output$SubForest <- renderPlot({
    Fplot2()
  }, width = SWinW, height = SWinH)
  
  # download forest plot
  output$DownSubForest <- downloadHandler(
    filename = function() {
      paste("forest-", Sys.time(), ".png", sep="")
    }, 
    content = function(file) {
      png(file, width = SWinW(), height = SWinH()) 
      Fplot2()
      dev.off()  
    }
  )
  
  # ----------------- Het plot (regression plot) ----------------------
  
  output$RegPlot <- renderPlotly({
    
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop("No file uploaded.") }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop("Cannot calculate effect size. Please change model or method.") } 
    }
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop("Your csv file doesn't contain Pub ID or Group letter. Studies can't be nested.")
      }
    }
    
    if (length(input$con.regplot)==0) { stop("No continuous variable selected.") }
    
    
    lab = input$con.regplot          # Name of the continuous variable
    Cont = as.numeric(yv[,lab])
    
    # Univariate meta-regression
    res <- rma(yi = yi, 
               vi = vi, 
               mods = ~Cont, 
               data = yv, 
               method = input$HetEstimator,
               test="knha")
    
    Effect.Size = yv[,"yi"]
    Size = 1/yv[,"vi"]
    xnew = seq(from = min(Cont,na.rm = T) - (max(Cont,na.rm = T) - min(Cont, na.rm = T))*0.1, 
               to = max(Cont,na.rm = T) + (max(Cont,na.rm = T) - min(Cont, na.rm = T))*0.1, 
               len = length(Cont))
    yfit = xnew*res$b[2]+res$b[1]
    
    if (input$RegXlab == "type the x-axis label") {
      reg_xlab = lab
    } else { reg_xlab = input$RegXlab }
    reg_ylab = input$RegYlab
    
    
    plot_ly() %>%
      add_trace(data = yv, 
                x = ~Cont, 
                y = ~Effect.Size, 
                size = ~Size, 
                mode = "markers", 
                showlegend = F, 
                name = "Observed",
                
                hoverinfo = 'text',
                text = ~paste('Effect: ', round(Effect.Size,2), 
                              '</br>SE: ', round(sqrt(yv$vi),2),
                              '</br>Author: ', Surname,
                              '</br>Year: ', Year
                ),
                marker = list(color = input$RegPtCol)
      ) %>% 
      
      add_trace(data = yv, 
                x = ~xnew, 
                y = yfit, 
                mode = "lines", 
                showlegend = F, 
                name = "Meta-regression",
                line = list(color = input$RegLCol)
      ) %>%
      layout(                      
        xaxis = list(title = reg_xlab,
                     showgrid = F, 
                     showline = T, 
                     zeroline = T, 
                     zerolinecolor = toRGB("gray90"), 
                     zerolinewidth = 1,
                     linecolor = toRGB("black"),
                     linewidth = 1.2),       
        
        yaxis = list(title = reg_ylab, 
                     showgrid = F,
                     showline = T, 
                     zeroline = T, 
                     zerolinecolor = toRGB("gray90"), 
                     zerolinewidth = 1,
                     linecolor = toRGB("black"),
                     linewidth = 1.2),
        autosize = F, 
        width = input$RegWidth, 
        height = input$RegHeight
      ) # layout
  })
  
  
  
  # ------------------- Select variables for bar plot --------------------
  # variables for stratified bar plot
  output$BarVar <- renderUI({
    
    yv <- dataset()
    if (is.null(yv)) { stop() }
    if (sum(is.na(yv$yi)) == nrow(yv)) { stop() } 
    
    use <- yv[, -which(names(yv) %in% c("Pub.ID","Year"))]
    
    selectInput(inputId = "dis.bar", 
                label = "Select a variable for bar plot",
                choices = as.list(names(use)), 
                multiple = F, selectize = F)
  })
  
  
  # ----------------- Bar plot ----------------------
  
  # function for composing multiple plots 
  # (copy from website...) no need to change 
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    if (is.null(layout)) {
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      for (i in 1:numPlots) {
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  # bar plot
  Bplot <- function(){
    
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop("No file uploaded.") }
      if (sum(is.na(yv$yi)+is.infinite(yv$yi)) == nrow(yv)) { stop("Cannot calculate effect size. Please change model or method.") } 
    }
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop("Your csv file doesn't contain Pub ID or Group letter. Studies can't be nested.")
      }
    }
    
    lab = input$dis.bar 
    if (is.null(lab)) { stop("No discrete variable selected") }
    nbar = length(lab)
    
    # stratified bar plot
    if (input$HetMethod == "sub"){
      plots <- list()  # new plots list
      for (j in 1:nbar){
        
        # stratified meta-analysis
        mod <- metagen(TE = yi, 
                       seTE = sqrt(vi), 
                       byvar = yv[,lab[[j]]],
                       data = yv, 
                       comb.fixed = F, 
                       method.tau = input$HetEstimator,
                       hakn = T
        )
        
        
        tem = aggregate(yv$ni, by = list(gro=yv[,lab[[j]]]), FUN = sum)
        temp = tem
        tlab = as.matrix(mod$bylevs)
        
        # to gurantee that variables in stratified meta-analysis output have same order
        # to variables in temp
        for (t in 1:nrow(temp)) {
          temp[t, ] = tem[tem[,1]==tlab[t],]
        }
        
        
        df = data.frame(xlabel = temp[,1],
                        effect = mod$TE.random.w,
                        # width of each bar 
                        w = sqrt(temp[,2])/sum(sqrt(temp[,2])),
                        # confidence interval
                        c_low = mod$lower.random.w, 
                        c_up = mod$upper.random.w)
        
        
        
        # y-axis label
        bar_ylab = input$BarYlab 
        # title
        if (input$BarTitle == "type the title") {
          bar_title = lab[[j]]
        } else { 
          bar_title = input$BarTitle
        }
        
        # y-axis limit
        if (input$BarYmin == "") {
          y_min = min(0,df$c_low, mod$lower.random)
        } else {
          y_min = as.numeric(input$BarYmin)
        }
        if (input$BarYmax == "") {
          y_max = max(0, df$c_up, mod$upper.random)
        } else {
          y_max = as.numeric(input$BarYmax)
        }
        
        
        
        
        p = ggplot(df, aes(x = xlabel, 
                           y = effect, 
                           width = w)) +
          
          geom_col(fill = "white", colour = "black") +
          # rectangles area: 95% CI of global estimate of efficacy
          geom_rect(aes(xmin = -Inf, 
                        xmax = Inf, 
                        ymin = mod$lower.random, 
                        ymax = mod$upper.random), 
                    fill = "grey90", alpha = 0.2) +
          # vertical error bars: 95% CI for individual estimates
          geom_linerange(aes(ymin = c_low, 
                             ymax = c_up), data=df) +
          
          geom_hline(yintercept = 0,linetype = 2) +
          
          scale_y_continuous(expand = c(0, 0), limits = c(y_min, y_max)) + 
          
          labs(x = "",y = bar_ylab) +
          ggtitle(bar_title) +
          theme(axis.line = element_line(colour = "black"),
                
                axis.text = element_text(size = 12),
                
                
                plot.title = element_text(size = input$BarTitleSize,
                                          hjust = 0.5,
                                          vjust = 0.5),
                
                axis.title = element_text(size = input$BarYlabSize,
                                          hjust = 0.5,
                                          vjust = 0.5), # y
                
                axis.text.x = element_text(size = input$BarLabSize,
                                           hjust = input$BarLabPos,
                                           angle = input$BarLabAngle),  # bar label
                
                
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank() 
          ) 
        
        
        plots[[j]] = p
      }  # for (j in 1:nbar)
      
      multiplot(plotlist = plots, cols = 2)
    } # if (input$het == "sub")
    
    
    # meta-regression bar plot
    if (input$HetMethod == "reg"){
      plots <- list()  # new plots list
      for (j in 1:nbar){
        
        mylist <- split(yv, yv[,lab[[j]]], drop = T)
        
        # Do meta-regression
        res <- rma(yi = yi, 
                   vi = vi, 
                   data = yv, 
                   mods = ~factor(yv[,lab[[j]]])-1, 
                   method = input$HetEstimator,
                   test="knha")
        
        # global meta-analysis
        glob <- rma(yi = yi, 
                    vi = vi, 
                    data = yv, 
                    method = input$HetEstimator,
                    test="knha")
        
        lev = rownames(res$b) # names for all the levels of one variable
        tem = aggregate(yv$ni, by = list(gro=yv[,lab[[j]]]), FUN = sum)
        
        
        
        df = data.frame(xlabel = substr(lev, 23, nchar(lev)),
                        y = res$b,
                        low = res$ci.lb,
                        up = res$ci.ub,
                        n = sqrt(tem[,2]))
        
        # y-axis label
        bar_ylab = input$BarYlab
        # title
        if (input$BarTitle == "type the title") {
          bar_title = lab[[j]]
        } else { 
          bar_title = input$BarTitle 
        }
        
        # y-axis limit
        if (input$BarYmin == "") {
          y_min = min(0, df$low, glob$ci.lb)
        } else {
          y_min = as.numeric(input$BarYmin)
        }
        if (input$BarYmax == "") {
          y_max = max(0, df$up, glob$ci.ub)
        } else {
          y_max = as.numeric(input$BarYmax)
        }
        
        
        p = ggplot(df, aes(x = xlabel, 
                           y = y, 
                           width = n/sum(n))) +
          
          geom_col(fill = "white", colour = "black") +
          
          geom_rect(aes(xmin = -Inf, 
                        xmax = Inf, 
                        ymin = glob$ci.lb, 
                        ymax = glob$ci.ub), 
                    fill = "grey90", alpha = 0.2) +
          
          geom_linerange(aes(ymin = low, ymax = up), data = df) +
          
          geom_hline(yintercept = 0, linetype = 2) +
          
          scale_y_continuous(expand = c(0, 0), limits = c(y_min, y_max)) + 
          
          labs(x = "",y = bar_ylab) +
          ggtitle(bar_title) +
          theme(axis.line = element_line(colour = "black"),
                
                axis.text = element_text(size = 12),
                
                
                plot.title = element_text(size = input$BarTitleSize,
                                          hjust = 0.5,
                                          vjust = 0.5),
                
                axis.title = element_text(size = input$BarYlabSize,
                                          hjust = 0.5,
                                          vjust = 0.5), # y
                
                axis.text.x = element_text(size = input$BarLabSize,
                                           hjust = input$BarLabPos,
                                           angle = input$BarLabAngle),  # bar label
                
                
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank() 
          )
        plots[[j]] = p
      }  # for (j in 1:nbar)
      multiplot(plotlist = plots, cols = 2)
    } # if (input$het == "reg")
  }
  
  
  # width and height of plot window
  BarW <- reactive({ input$BarWidth })
  BarH <- reactive({ input$BarHeight })
  
  # display bar plot
  output$BarPlot <- renderPlot({
    Bplot()
  }, width = BarW, height = BarH)
  
  # download bar plot
  output$DownBar <- downloadHandler(
    filename = function() {
      paste("Bar-", Sys.time(), ".png", sep = "")
    },
    content = function(file) {
      png(file, width = BarW(), height = BarH()) # open the png device
      Bplot()
      dev.off()  # turn the device off
    }
  )
  
  
  
  Binfo <- function(){
    
    lab = input$dis.bar
    if (is.null(lab)) { stop("No variable selected") }
    nbar = length(lab)
    
    if (input$DataType == "n0") { yv <- yivi() }
    if (input$DataType == "n1") { yv <- Nest() }
    
    # stratified bar info
    if (input$HetMethod == "sub"){
      
      # stratified meta-analysis
      mod <- metagen(TE = yi, 
                     seTE = sqrt(vi), 
                     byvar = yv[,lab],
                     data = yv, 
                     comb.fixed = F,
                     method.tau = input$HetEstimator,
                     hakn = T)
      
      
      tem = aggregate(yv$ni, by = list(gro=yv[,lab]), FUN = sum)
      temp = tem
      tlab = as.matrix(mod$bylevs)
      
      # to gurantee that variables in stratified meta-analysis output have same order
      # to variables in temp
      for (t in 1:nrow(temp)) {
        temp[t,] = tem[tem[,1] == tlab[t],]
      }
      
      ndig = 2
      
      df = data.frame(Group = temp[,1],
                      ES = round(mod$TE.random.w, ndig),
                      
                      # confidence interval
                      ci.low = round(mod$lower.random.w, ndig),
                      ci.up = round(mod$upper.random.w, ndig),
                      
                      No.Animals = temp[,2],
                      No.Studies = mod$k.w,
                      stringsAsFactors = F)
      df = rbind(df,
                 c("Global", 
                   round(mod$TE.random, ndig),
                   round(mod$lower.random, ndig),
                   round(mod$upper.random, ndig), 
                   sum(yv[,"ni"]), mod$k)
      )
    } # if (input$HetMethod == "sub")
    
    # meta-regression bar info
    if (input$HetMethod == "reg"){
      
      # Split into subgroups
      mylist <- split(yv, yv[,lab], drop = T)
      # Do meta-regression
      res <- rma(yi = yi, 
                 vi = vi, 
                 mods = ~factor(yv[,lab])-1,
                 data = yv, 
                 method = input$HetEstimator,
                 test="knha")
      # global meta-analysis
      glob <- rma(yi = yi, 
                  vi = vi, 
                  data = yv,
                  method = input$HetEstimator,
                  test="knha")
      
      lev = rownames(res$b) # names for all the levels of one variable
      
      ndig = 2
      df = data.frame(Group = substr(lev, 18, nchar(lev)),
                      ES = round(res$b, ndig),
                      ci.low = round(res$ci.lb ,ndig),
                      ci.up = round(res$ci.ub, ndig),
                      stringsAsFactors = F)
      df = rbind(df,
                 c("Global", 
                   round(glob$b, ndig), 
                   round(glob$ci.lb, ndig), 
                   round(glob$ci.ub, ndig))
      )
    } # if (input$het == "reg")
    return(df)
  } 
  
  output$BarInfoTable <- renderTable({
    Binfo()
  })
  
  
  # ---------------- Trim-and-Fill ------------------------- 
  output$TafOutput <- renderPrint({
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop("No file uploaded.") }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop("Cannot calculate effect size. Please change model or method.") } 
    }
    
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop("Your csv file doesn't contain Pub ID or Group letter. Studies can't be nested.")
      }
    }
    
    
    res <- rma(yi = yi, vi = vi, ni = ni, 
               data = yv, 
               method = input$HetEstimator,
               test="knha")
    trimfill(res, estimator = "L0", side = input$TafSide)
    
    # res <- metagen(TE=yi, seTE=sqrt(vi), data=yv, comb.fixed=F, method.tau=input$meth)
    # summary(trimfill(res))
  })
  
  # ----------------- Funnel plot ----------------------
  Funplot <- function(){
    
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop("No file uploaded.") }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop("Cannot calculate effect size. Please change model or method.") } 
    }
    
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop("Your csv file doesn't contain Pub ID or Group letter. Studies can't be nested.")
      }
    }
    
    
    res <- rma(yi = yi, vi = vi, ni = ni, 
               data = yv, 
               method = input$HetEstimator,
               test="knha")
    taf <- trimfill(res, estimator = "L0", side = input$TafSide)
    
    xmin = min(min(taf$yi)-0.1*(max(taf$yi)-min(taf$yi)), 0)
    xmax = max(max(taf$yi)+0.1*(max(taf$yi)-min(taf$yi)), 0)
    
    
    # y-axis is inverse of standard error    
    if (input$TafYaxis == "seinv"){
      ymax = max(max(1/sqrt(taf$vi))+0.1*(max(1/sqrt(taf$vi))-min(1/sqrt(taf$vi))), 0)
      if (input$TafFill == "Yes"){
        # show imputed studies on funnel plot
        funnel(taf, xlim = c(xmin, xmax),  ylim = c(0.000000001,ymax),
               back = "white", 
               refline = res$b, 
               yaxis = input$TafYaxis,
               xlab = input$TafXlab, 
               ylab = input$TafYlab)
      } else {
        # not show imputed studies on funnel plot
        funnel(res, xlim = c(xmin, xmax), ylim = c(0.000000001,ymax),
               back = "white",
               yaxis = input$TafYaxis,  
               xlab = input$TafXlab, 
               ylab = input$TafYlab)
      }
    }  
    
    # y-axis is inverse of the square-root sample size
    if (input$TafYaxis == "sqrtninv"){
      ymax = max(max(1/sqrt(taf$ni))+0.1*(max(1/sqrt(taf$ni))-min(1/sqrt(taf$ni))), 0)
      # show imputed studies on funnel plot
      if (input$TafFill == "Yes"){
        funnel(taf, xlim = c(xmin, xmax), ylim = c(0.000000001, ymax), 
               back = "white", 
               refline = res$b, 
               yaxis = input$TafYaxis,
               xlab = input$TafXlab, 
               ylab = input$TafYlab)
      } else {
        # not show imputed studies on funnel plot
        funnel(res, xlim = c(xmin, xmax), ylim = c(0.000000001, ymax), 
               back = "white", 
               yaxis = input$TafYaxis,
               xlab = input$TafXlab, 
               ylab = input$TafYlab
        )
      }
    }
    abline(v = taf$b, lty = 2)
  }
  
  
  # width and height of plot window
  FnlW <- reactive({ input$FunnelWidth })
  FnlH <- reactive({ input$FunnelHeight })
  
  # displat funnel plot
  output$FunnelPlot <- renderPlot({
    Funplot()
  },  width = FnlW, height = FnlH)
  
  # download funnel plot
  output$DownFunnel <- downloadHandler(
    filename = function() {
      paste("funnel-", Sys.time(), ".png", sep = "")
    },
    content = function(file) {
      png(file, width = FnlW(), height = FnlH()) 
      Funplot()
      dev.off()  
    }
  )
  
  # ---------------- Egger's regression ------------------------- 
  output$EggerOutput <- renderPrint({
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop("No file uploaded.") }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop("Cannot calculate effect size. Please change model or method.") } 
    }
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop("Your csv file doesn't contain Pub ID or Group letter. Studies can't be nested.")
      }
    }
    
    Standardised.effect <- (yv$yi)/sqrt(yv$vi)
    Precision <- 1/sqrt(yv$vi)
    # egger regression
    egger <- lm(Standardised.effect ~ Precision)
    summary(egger)
    
  })
  
  
  
  # ----------------- Egger's regression plot ----------------------
  output$EggerPlot <- renderPlotly({
    if (input$DataType == "n0"){ 
      yv <- yivi() 
      if (is.null(yv)) { stop("No file uploaded.") }
      if (sum(is.na(yv$yi)) == nrow(yv)) { stop("Cannot calculate effect size. Please change model or method.") } 
    }
    if (input$DataType == "n1"){ 
      yv <- Nest() 
      if (is.null(yv)) {
        stop("Your csv file doesn't contain Pub ID or Group letter. Studies can't be nested.")
      }
    }
    
    Standardised.effect <- (yv$yi)/sqrt(yv$vi)
    Precision <- 1/sqrt(yv$vi)
    
    egger <- lm(Standardised.effect ~ Precision)
    
    # calculate 95% confidence band
    xnew = seq(0, to = max(Precision)+(max(Precision)-min(Precision))*0.1, 
               len = length(Precision))
    ynew = egger$coefficients[1] + xnew*egger$coefficients[2]
    pnew <- predict(egger, newdata = data.frame(Precision=xnew),
                    interval = "confidence", level = 0.95)
    
    
    
    
    plot_ly() %>%
      
      # 95% CI ribbon
      add_ribbons(x = ~xnew, 
                  ymin = pnew[,"lwr"], 
                  ymax = pnew[,"upr"], 
                  line = list(color = input$EggerRibCol),
                  fillcolor = input$EggerRibCol, 
                  showlegend = F, name = "95% confidence") %>%
      # regression line
      add_trace(x = ~xnew, 
                y = ynew, 
                mode = "lines", 
                line = list(color = input$EggerLCol),
                showlegend = F, name = "Egger's fit") %>%
      
      # scatter plot
      add_trace(x = ~Precision, 
                y = ~Standardised.effect, 
                mode = "markers", 
                marker = list(color = input$EggerPtCol),
                showlegend = F, name = "Observed") %>%
      
      layout(                        
        xaxis = list(title = input$EggerXlab, showgrid = F, 
                     showline = T, linewidth = 1,
                     rangemode="tozero",
                     zeroline = T, zerolinecolor = toRGB("gray90"), zerolinewidth = 1),       
        yaxis = list(title = input$EggerYlab, showgrid = F, 
                     showline = T, linewidth = 1, 
                     rangemode="tozero",
                     zeroline = T, zerolinecolor = toRGB("gray90"), zerolinewidth = 1),
        autosize = F, 
        width = input$EggerWidth, 
        height = input$EggerHeight
      )
  })
  
  
})
