# This is a Shiny web application for meta-analysis
# by Qianying Wang @CAMARADES
# ui.R
# https://camarades.shinyapps.io/meta-analysis-app/


library(shiny)
library(metafor)
library(meta)
library(shinythemes)
library(plotly)
library(ggplot2)
library(colourpicker)
library(dplyr)
library(RCurl)

shinyUI(fluidPage(
  
  theme = shinytheme("spacelab"),
  
  headerPanel(
    list(HTML('<img src="images.jpg" height=90 width=60 /> '), 
         HTML('<font size="5"> Camarades meta-analysis tool </font>')),
           #"Camarades meta-analysis tool"),
    windowTitle="Camarades meta-analysis tool"
  ),
  
  sidebarLayout(
    sidebarPanel(

      fileInput(inputId = "file", label = "Upload CSV File",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
      ),
      
      # helpText(a("User guide", href="http://www.dcn.ed.ac.uk/camarades/metaAnalysisUserGuide.html", 
      #            target="_blank")
      # ),
      
      radioButtons(inputId = "CompareType", label = "Select comparison",
                   choices = list("Analysis of model" = "SC", 
                                  "Analysis of intervention" = "SCT"), 
                   selected = "SCT"),
    
      hr(),
      
      radioButtons(inputId = "EffectMeasure", label = "Select effect size measure",
                   choices = list("Have calculated effect size and SE" = "GEN",
                                  "Normalised mean difference" = "NMD", 
                                  "Standardised mean difference" = "SMD",
                                  "Odds ratio" = "OR"), 
                   selected = "NMD"),
      
      
      hr(),
      
      
      radioButtons(inputId = "HetMethod", label = "Select heterogeneity analysis method",
                   choices = list("Stratified meta-analysis" = "sub", 
                                  "Meta-regression" = "reg"),
                   selected = "sub"),
      
      hr(),
      selectInput(inputId = "HetEstimator", label = "Select heterogeneity estimator",
                  choices = list("DerSimonian-Laird" = "DL", 
                                 "Hedges" = "HE", 
                                 "Hunter-Schmidt" = "HS", 
                                 "Sidik-Jonkman" = "SJ", 
                                 "Maximum-likelihood" = "ML", 
                                 "Restricted maximum-likelihood" = "REML", 
                                 "Empirical Bayes" = "EB"),
                  selected = "REML"),
      
      hr(),
      checkboxInput(inputId = "KnHaTest", label = "Fit with Knapp and Hartung method", value = F)
    ),
    
    mainPanel(
      tabsetPanel(
        
        tabPanel("Select data", 
                 verticalLayout(
                   wellPanel(
                     fluidRow(
                       column(6,
                              radioButtons(inputId = "DataType", label = "Select data you want to use for analysis",
                                           choices = list("Raw data" = "n0"), 
                                                          # "Data nested by group letter" = "n1"),
                                                          # "Data nested by group letter and another variable" = "n2"),
                                           selected = "n0")
                       ),
                       column(6,
                              radioButtons("FileType", "File type",
                                           choices = c("csv", "tsv")),
                              downloadButton('DownTable', 'Download table')
                       )
                     )
                     
                     
                     # radioButtons(inputId = "DataType", label = "Select data you want to use for analysis",
                     #              choices = list("Raw data" = "n0", 
                     #                             "Data nested by group letter" = "n1",
                     #                             "Data nested by group letter and another variable" = "n2"),
                     #              selected = "n0")
                   ),
                   
                   # conditionalPanel(condition = "input.DataType=='n2'",
                   #                  fluidRow(
                   #                    column(6,
                   #                           uiOutput("NestVar")
                   #                    ),
                   #                    column(6,
                   #                           uiOutput("NumLevel")
                   #                    )
                   #                  ),
                   #                  uiOutput("grbox")
                   #  ),

                   dataTableOutput("DT")
                   
                 )
        ),
        
        
        tabPanel("Meta-analysis",
                 verticalLayout(
                   verbatimTextOutput("GlobalOutput"),
                   wellPanel(
                     fluidRow(
                       column(6,
                              textInput(inputId = "ForXlab", label = "Label for x-axis", value = "Effect & REML"),
                              sliderInput(inputId = "GapLeft", label = "Left gap", min = 0, max = 10, value = 1),
                              sliderInput(inputId = "ForWinWidth", label = "Figure width (px)", min = 600, max = 2000, value = 600),
                              sliderInput(inputId = "ForWidth", label = "Plot width", min = 4, max = 50, value = 8)
                       ),
                       column(6,
                              fluidRow(
                                column(4, 
                                       selectInput(inputId = "ForOrder", label = "Order of studies",
                                                   choices = list("original order" = "",
                                                                  "increasing effect size" = "ine", 
                                                                  "increasing weight" = "inw", 
                                                                  "increasing year" = "iny", 
                                                                  "decreasing effect size" = "dee", 
                                                                  "decreasing weight" = "dew", 
                                                                  "decreasing year" = "dey"),
                                                   selected = "")                                       
                                ),
                                column(4, 
                                       colourInput(inputId = "ForSqCol", label = "Square color", value = "#E69F00")
                                ),
                                column(4,
                                       colourInput(inputId = "ForDiaCol", label = "Diamond color", value = "#56B4E9")
                                )
                              ),

                              sliderInput(inputId = "GapRight", label = "Right gap", min = 0, max = 10, value = 1),
                              numericInput(inputId = "ForWinHeight", label = "Figure height (px)", value = 400, min = NA, max = NA, step = NA),
                              # sliderInput(inputId = "ForWinHeight", label = "Figure height (px)", min=400, max=2000, value = 400),
                              
                              # HTML("<br><br>"),
                              
                              fluidRow(
                                column(6, 
                                       textInput(inputId = "ForXmin", label = "x-min", value = "")
                                ),
                                column(6,
                                       textInput(inputId = "ForXmax", label = "x-max", value = "")
                                )
                              ),
                              
                              fluidRow(
                                column(6, 
                                       checkboxInput("ShowWeight", label = "Show weight", value = T)
                                ),
                                column(6,
                                       downloadButton("DownForest", "Download forest plot")
                                )
                              )
                       )
                     )
                   ),
                   plotOutput(outputId = "ForestPlot")
                 )
        ),
        
        tabPanel("Heterogeneity",
                 verticalLayout(
                   wellPanel(
                     conditionalPanel(condition = "input.HetMethod == 'sub'",
                                      uiOutput("SubVar")
                     ),
                     conditionalPanel(condition = "input.HetMethod == 'reg'",
                                      fluidRow(
                                        column(6, uiOutput("RegContVar")),
                                        column(6, uiOutput("RegDiscVar"))
                                      )
                     ) 
                   ),
                   # stratified meta-analysis
                   conditionalPanel(condition = "input.HetMethod == 'sub'",
                                    helpText("Overall meta-analysis:"),
                                    verbatimTextOutput("AllOutput", placeholder = T),
                                    
                                    helpText("Results for subgroups (random effect model):"),
                                    verbatimTextOutput("SubOutput", placeholder = T),
                                    
                                    helpText("Test for subgroup difference (chi-square test):"),
                                    verbatimTextOutput("SubTest", placeholder = T)
                                    ),
                   # meta-regression
                   conditionalPanel(condition = "input.HetMethod == 'reg'",
                                    verbatimTextOutput("RegOutput", placeholder = T) ) 
                 )
        ),
        
        tabPanel("Heterogeneity plot",
                 verticalLayout(
                   
                     # subgroup forest plot
                     conditionalPanel(condition = "input.HetMethod == 'sub'",
                                      wellPanel(
                                        
                                        fluidRow(
                                          column(6,
                                                 uiOutput("SubPlotVar"),
                                                 sliderInput(inputId = "SubGapLeft", label = "Left gap", min = 0, max = 10, value = 1),
                                                 sliderInput(inputId = "SubWinWidth", label = "Width (px)", min = 600, max = 2000, value = 600),
                                                 sliderInput(inputId = "SubWidth", label = "Plot width (cm)", min = 4, max = 50, value = 8)
                                        
                                          ),
                                          column(6,
                                                 
                                                 fluidRow(
                                                   column(4, 
                                                          textInput(inputId = "SubXlab", label = "Label for x-axis", value = "Effect & REML")
                                                   ),
                                                   column(4, 
                                                          colourInput(inputId = "SubSqCol", label = "Square color", value = "#E69F00")
                                                   ),
                                                   column(4,
                                                          colourInput(inputId = "SubDiaCol", label = "Diamond color", value = "#56B4E9")
                                                   )
                                                 ),
                                                 
                                                 sliderInput(inputId = "SubGapRight", label = "Right gap", min = 0, max = 10, value = 1),
                                                 numericInput(inputId = "SubWinHeight", label = "Height (px)", value = 400, min = NA, max = NA, step = NA),
                                                 # sliderInput(inputId = "SubWinHeight", label = "Height (px)", min = 400, max = 2000, value = 400),
                                                
                                                 HTML("<br><br>"),
                                                 
                                                 fluidRow(
                                                   column(3, 
                                                          checkboxInput("ShowFixWeight", label = "Show w-fixed", value = T)
                                                   ),
                                                   column(3, 
                                                          checkboxInput("ShowRandWeight", label = "Show w-random", value = T)
                                                   ),
                                                   column(6, 
                                                          downloadButton("DownSubForest", "Download forest plot")
                                                   )
                                                 )
                                          )
                                        ) # fluidRow
                                      ) # wellPanel
                                      
                                      
                     ), # conditionPanel
                     
                     HTML("<br>"),
                     # meta-regression plot
                     conditionalPanel(condition = "input.HetMethod == 'reg'",
                                      wellPanel(
                                        fluidRow(
                                          column(6,
                                                 uiOutput("RegPlotVar"),
                                                 fluidRow(
                                                   column(6, 
                                                          textInput(inputId = "RegXlab", label = "Label for x-axis", value = "type the x-axis label")
                                                   ),
                                                   column(6,
                                                          textInput(inputId = "RegYlab", label = "Label for y-axis", value = "type the y-axis label")
                                                   )
                                                 ),
                                                 fluidRow(
                                                   column(6, 
                                                          colourInput(inputId = "RegPtCol", label = "Point color", value = "#131414")
                                                   ),
                                                   column(6,
                                                          colourInput(inputId = "RegLCol", label = "Line color", value = "#FC130F")
                                                   )
                                                 )
                                          ),
                                          column(6,
                                                 sliderInput(inputId = "RegWidth", label = "Width (px)", min = 300, max = 2000, value = 700),
                                                 HTML("<br><br>"),
                                                 sliderInput(inputId = "RegHeight", label = "Height (px)", min = 300, max = 2000, value = 400)
                                          )
                                        ) # fluidRow
                                      ) # wellPanel
                     ), # conditionPanel
                     HTML("<br>"),
                   
                     conditionalPanel(condition = "input.HetMethod == 'sub'",
                                      plotOutput("SubForest") ),
                     conditionalPanel(condition = "input.HetMethod == 'reg'",
                                      plotlyOutput("RegPlot") ) 
                 )
        ),
        
        
        
        tabPanel("Bar plot",
                 verticalLayout(
                   wellPanel(
                     fluidRow(
                       column(6, 
                              uiOutput("BarVar"),
                              sliderInput(inputId = "BarHeight", label = "Height (px)", min = 400, max = 2000, value = 400),
                              sliderInput(inputId = "BarTitleSize", label = "Title size", min = 8, max = 30, value = 12),
                              sliderInput(inputId = "BarYlabSize", label = "Y-axis label size", min = 8, max = 30, value = 12),
                              sliderInput(inputId = "BarLabAngle", label = "Bar label angle", min = 0, max = 90, step = 1, value = 0)
                              ),
                       column(6,
                              fluidRow(
                                column(6, 
                                       textInput(inputId = "BarYlab", label = "Label for y-axis", value = "Effect size")
                                ),       
                                column(6,
                                       textInput(inputId = "BarTitle", label = "Title", value = "type the title")
                                )
                              ),
                              
                              fluidRow(
                                column(6, 
                                       textInput(inputId = "BarYmin", label = "y-min", value = "")
                                ),
                                column(6,
                                       textInput(inputId = "BarYmax", label = "y-max", value = "")
                                )
                              ),
                              
                              sliderInput(inputId = "BarWidth", label = "Width (px)", min = 100, max = 2000, value = 400),
                              sliderInput(inputId = "BarLabSize", label = "Bar label size", min = 8, max = 30, value = 12),
                              sliderInput(inputId = "BarLabPos", label = "Bar label position", min = 0.1, max = 3, value = 0.5),
                              HTML("<br>"),
                              downloadButton("DownBar", "Download bar plot")
                       )
                     )
                      )
                 ),
                 tableOutput("BarInfoTable"),
                 plotOutput("BarPlot")
                 
        ),
        
        
        
        tabPanel("Trim and Fill",
                 verticalLayout(
                   verbatimTextOutput("TafOutput"),
                   wellPanel(
                     fluidRow(
                       column(6, 
                              radioButtons(inputId = "TafYaxis", label = "Select y-axis",
                                           choices = list("Inverse of the standard error" = "seinv",
                                                          "Inverse of the square-root sample size" = "sqrtninv"),
                                           selected = "seinv"),
                              
                              radioButtons(inputId = "TafFill", label = "Select funnel plot type",
                                           choices = list("Show imputed studies" = "Yes",
                                                          "Show only published studies" = "No"),
                                           selected = "Yes"),
                              
                              radioButtons(inputId = "TafSide", label = "Side of funnel plot the missing studies imputed",
                                           choices = list("Left" = "left",
                                                          "Right" = "right"),
                                           selected = "left"),
                              downloadButton("DownFunnel", "Download funnel plot")
                       ),
                       column(6,
                              fluidRow(
                                column(6, 
                                       textInput(inputId = "TafXlab", label = "Label for x-axis", value = "Effect size")
                                ),
                                column(6,
                                       textInput(inputId = "TafYlab", label = "Label for y-axis", value = "type the y-axis label")
                                )
                              ),
                              sliderInput(inputId = "FunnelWidth", label = "Width (px)", min = 300, max = 1200, value = 600),
                              sliderInput(inputId = "FunnelHeight", label = "Height (px)", min = 300, max = 1200, value = 400)
                              # HTML("<br><br>"),
                              
                       )
                     )
                   ),
                   plotOutput("FunnelPlot")
                 )
        ),
        
        tabPanel("Egger's Regression", 
               
                 wellPanel(
                   fluidRow(
                     column(6,
                            fluidRow(
                              column(6, 
                                     textInput(inputId = "EggerXlab", label = "Label for x-axis", value = "Precision (1/SE)")
                              ),
                              column(6,
                                     textInput(inputId = "EggerYlab", label = "Label for y-axis", value = "Effect size/SE")
                              )
                            ),
                            HTML("<br>"),
                            fluidRow(
                              column(4, 
                                     colourInput(inputId = "EggerPtCol", label = "Point color", value = "#131414")
                              ),
                              column(4,
                                     colourInput(inputId = "EggerLCol", label = "Line color", value = "#3F51B5")
                              ),
                              column(4,
                                     colourInput(inputId = "EggerRibCol", label = "Ribbon color", value = "#C8F0E7")
                              )
                            )
                     ),
                     column(6,
                            sliderInput(inputId = "EggerWidth", label = "Width (px)", min = 300, max = 2000, value = 700),
                            sliderInput(inputId = "EggerHeight", label = "Height (px)", min = 300, max = 2000, value = 400)
                     )
                   ) # fluidRow
                 ), # wellPanel
                 
                 plotlyOutput("EggerPlot"),
                 HTML("<br>"),
                 verbatimTextOutput("EggerOutput")
        )
        
        
      )
    )
  )
))




