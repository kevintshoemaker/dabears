###
### Shiny App
###

# Load necessary packages
library(shiny)
library(MuMIn)
library(openCR)
library(readxl)
library(splines)
library(xtable)
library(httr)

# Define UI
ui <- fluidPage(
  titlePanel("Nevada Black Bear Abundance Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("datafile", "Choose Excel File", 
                accept = c(".xlsx")),
      numericInput("basis_functions", "Number of Basis Functions", 5, min = 4, max = 6),
      selectInput("data_set", "Data Set to Use", choices = c("marked bears", "all")),
      checkboxInput("remove_cubs", "Remove Cubs That Were Never Seen Again", TRUE),
      
      selectInput("phi_vars", "Variables for phi", multiple = TRUE, choices = c("1","sex", "year","basis functions")),
      selectInput("p_vars", "Variables for p", multiple = TRUE, choices = c("1","sex", "year","basis functions")),
      selectInput("N_vars", "Variables for N", multiple = TRUE, choices = c("1","sex", "year","basis functions")),
      numericInput("start_year", "Start Year", 1997, min = 1997, max = 2022),
      numericInput("end_year", "End Year", 2022, min = 1997, max = 2022),
      actionButton("analyze", "Analyze Data")
    ),
    
    mainPanel(
      tableOutput("aic_table"),
      verbatimTextOutput("model_summary"),
      plotOutput("population_plot")
    )
  )
)


setwd("~/Dropbox/GitHub/Nevada-Black-Bear-Abundance")

#####################
input$datafile$datapath = "~/Dropbox/GitHub/Nevada-Black-Bear-Abundance/Data/CaptureHistories.xlsx" 
input$start_year=1993
input$end_eqr=2022
input$remove_cubs=TRUE
input$basis_functions=5
input$phi_vars=c("sex", "basis functions")
input$p_vars=c("sex")
input$N_vars=c("sex")
#####################

# Define server logic
server <- function(input, output, session) {
  
  observeEvent(input$analyze, {
    req(input$datafile)
    
    # Load the data
    file_path <- input$datafile$datapath
    
    data1 <- read_xlsx(file_path, sheet = "marked bears")
    
    start.row <- 3  
    
    # Update years based on user input
    years <- input$start_year:input$end_year  
    T <- length(years)
    seasons <- 1:4 
    J <- length(seasons)
    
    remove.cubs <- input$remove_cubs
    sex <- "Both"
    primary.sessions <- "annual"
    #? data.set <- input$data_set
    basis.functions <- input$basis_functions
    
    basis_vector <- character(basis.functions)
    
    # Loop to generate the basis functions
    for (i in 1:basis.functions) {
      basis_vector[i] <- paste0("bs", i)
    }
    
    # Combine the basis functions into a single character string
    basis_string <- paste(basis_vector, collapse = " + ")
    
    all.years.of.data <- 1997:((ncol(data1) - 16) / 4 + 1996)
    all.years.ind <- rep(all.years.of.data, each = length(seasons))
    all.years.ind <- all.years.ind[-length(all.years.of.data) * J]
    
    n1.total <- length(start.row:nrow(data1)) 
    ch1.all.tmp <- data1[start.row:n1.total, 17:(ncol(data1) - 1)]
    ch1.all.m.tmp <- as.matrix(ch1.all.tmp, nrow(ch1.all.tmp), ncol(ch1.all.tmp))
    ch1.all.m <- cbind(ch1.all.m.tmp)
    
    sex1.all <- data1[start.row:n1.total, 8]
    sex1.all.f <- factor(ifelse(sex1.all == "M", "Male", "Female"))
    
    age1.all.tmp <- data.frame(data1[start.row:n1.total, 7])
    age1.all.tmp[age1.all.tmp == "7 WK"] <- 7 / 52
    age1.all.tmp[age1.all.tmp == "8 WEEKS"] <- 8 / 52
    age1.all.tmp[age1.all.tmp == "4 MO"] <- 4 / 12
    age1.all.tmp[age1.all.tmp == "4 MOS"] <- 4 / 12
    age1.all.tmp[age1.all.tmp == "5 MO"] <- 5 / 12
    age1.all.tmp[age1.all.tmp == "6 MO"] <- 6 / 12
    age1.all.tmp[age1.all.tmp == "7 MO"] <- 7 / 12
    age1.all.tmp[age1.all.tmp == "7 MOS"] <- 7 / 12
    age1.all.tmp[age1.all.tmp == "8 MO"] <- 8 / 12
    age1.all.tmp[age1.all.tmp == "8 MOS"] <- 8 / 12
    age1.all.tmp[age1.all.tmp == "9 MO"] <- 9 / 12
    age1.all.tmp[age1.all.tmp == "9 MOS"] <- 9 / 12
    age1.all.tmp[age1.all.tmp == "10 MO"] <- 10 / 12
    age1.all.tmp[age1.all.tmp == "10 MOS"] <- 10 / 12
    age1.all.tmp[age1.all.tmp == "11 MO"] <- 11 / 12
    age1.all.tmp[age1.all.tmp == "11 MOS"] <- 11 / 12
    age1.all.tmp[age1.all.tmp == "15 MOS"] <- 15 / 12
    names(age1.all.tmp) <- "age"
    age1.all <- as.numeric(age1.all.tmp$age)
    
    col.ind <- all.years.ind %in% years
    row.ind1 <- apply(ch1.all.m[, col.ind], 1, sum) > 0
    
    if (remove.cubs) {
      keep <- (apply(ch1.all.m, 1, sum) > 1 | age1.all > (16 / 12))
    }
    if (!remove.cubs) keep <- rep(TRUE, nrow(ch1.all.m))
    
    ch1.m <- ch1.all.m[row.ind1 & keep, col.ind]
    n <- nrow(ch1.m)
    
    sex1 <- sex1.all.f[row.ind1 & keep]
    age1 <- age1.all[row.ind1 & keep]
    
    season <- factor(rep(1:4, T))
    season <- season[-(T * J)]
    year <- rep(1:T, each = J)
    year <- year[-(T * J)]
    year.f <- factor(year)
    year2 <- year^2
    bs <- bs(year, basis.functions, intercept = TRUE)
    
    sescov <- data.frame(season = season,
                         year = year,
                         year.f = year.f,
                         year2 = year2)
    
    for (i in 1:basis.functions) {
      new <- bs[, i]
      sescov[, ncol(sescov) + 1] <- new
      colnames(sescov)[ncol(sescov)] <- paste0("bs", i)
    }
    
    
    sex=sex1
    age=age1
    
    bear.df=data.frame(ch=pasty(ch1.m),sex=sex1)
    bearCH=suppressWarnings(unRMarkInput(bear.df))
    
    
    
    
    
    
    
    
    # Select variables for phi, p, and N
    phi_vars <- input$phi_vars
    if("basis functions" %in% phi_vars){
      phi_vars=phi_vars[phi_vars!="basis functions"]
      phi_formula.tmp <- paste(c(basis_string, phi_vars),collapse=" + ")
      phi_forumla=as.formula(paste("phi~", phi_formula.tmp))
    }else{
      phi_formula.tmp <- paste(c(basis_string, phi_vars),collapse=" + ")
      phi_forumla=as.formula(paste("phi~", phi_formula.tmp))
    }
    
    p_vars <- input$p_vars
    if("basis functions" %in% p_vars){
      p_vars=p_vars[p_vars!="basis functions"]
      p_formula.tmp <- paste(c(basis_string, p_vars),collapse=" + ")
      p_forumla=as.formula(paste("p~", p_formula.tmp))
    }else{
      p_formula.tmp <- paste(c(basis_string, p_vars),collapse=" + ")
      p_forumla=as.formula(paste("p~", p_formula.tmp))
    }
    N_vars <- input$N_vars
    if("basis functions" %in% N_vars){
      N_vars=N_vars[N_vars!="basis functions"]
      N_formula.tmp <- paste(c(basis_string, N_vars),collapse=" + ")
      N_forumla=as.formula(paste("N~", N_formula.tmp))
    }else{
      N_formula.tmp <- paste(c(basis_string, N_vars),collapse=" + ")
      N_forumla=as.formula(paste("N~", N_formula.tmp))
    }
    
    
    
    # Generate model based on the number of basis functions selected
      mod <- openCR.fit(bearCH, type = "JSSAb", 
                      model =  list(p_formula,
                                    phi_formula,
                                    N_formula),
                      distributions = c('binomial'),
                      sessioncov=sescov)
    
  })
}


# Run the application 
shinyApp(ui = ui, server = server)


