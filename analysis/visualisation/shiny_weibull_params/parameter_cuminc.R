##******************************************##
## Visualise chosen Weibull parametrisation ##
##******************************************##


library(shiny)
library(tidyverse)



# UI ----------------------------------------------------------------------


ui <- fluidPage(
  
  # App title ----
  titlePanel("Parameters"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar to demonstrate various slider options ----
    sidebarPanel(
      
      h3("Weibull params REL (event 1)", align = "center"),
      sliderInput("shape_1", "Shape Cause 1:",
                  min = 0, max = 6,
                  value = 1, step = 0.01),
      
      # Input:  ----
      sliderInput("base_1", "BaseHaz Cause 1:",
                  min = 0, max = 3,
                  value = 1, step = 0.01),
      
      # Input:  ----
      sliderInput("beta_1", "Beta Cause 1:",
                  min = -2, max = 2,
                  value = 0, step = 0.01),
      
      # Input:  ----
      sliderInput("gamma_1", "Gamma Cause 1:",
                  min = -2, max = 2,
                  value = 0, step = 0.01),
      
      h3("Weibull params NRM (event 2)", align = "center"),
      sliderInput("shape_2", "Shape Cause 2:",
                  min = 0, max = 6,
                  value = 1, step = 0.01),
      
      # Input:  ----
      sliderInput("base_2", "BaseHaz Cause 2:",
                  min = 0, max = 3,
                  value = 1, step = 0.01),
      
      # Input:  ----
      sliderInput("beta_2", "Beta Cause 2:",
                  min = -2, max = 2,
                  value = 0, step = 0.01),
      
      # Input:  ----
      sliderInput("gamma_2", "Gamma Cause 2:",
                  min = -2, max = 2,
                  value = 0, step = 0.01)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(width = 6,
              fluidRow(
                tabsetPanel(type = "tabs",
                            tabPanel(title = "Cumulative incidence",
                                     plotOutput("plot_CI")),
                            tabPanel(title = "Condit. cumul. Haz",
                                     plotOutput("plot_haz")),
                            tabPanel(title = "Condit. Haz",
                                     plotOutput("plot_condhaz")),
                            tabPanel(title = "Reverse KM.",
                                     plotOutput("rev_KM")),
                            tabPanel(title = "Prop Events",
                                     verbatimTextOutput("descrip"))
              ),
              fluidRow(align = "center",
                       
                hr(),
                       
                h3("Censoring and covariates", align = "center"),
                
                sliderInput("rate_cens", "Exp censoring rate:",
                            min = 0, max = 5,
                            value = 0.5, step = 0.01),
                
                # Input:  ----
                sliderInput("X", "Value of covariate X",
                            min = -3, max = 3,
                            value = 0, step = 0.1),
                
                # Input:  ----
                sliderInput("Z", "Value of covariate Z",
                            min = -3, max = 3,
                            value = 0, step = 0.1)
              )
        
      )
    )
  )
)


# Server  -----------------------------------------------------------------


server <- function(input, output) {
  
  theme_set(theme_bw(base_size = 16) +
              theme(plot.title = element_text(hjust = 0.5)))
  
  # Reactive expression to create data frame of all input values ----
  dat_cumin <- reactive({
    
    generate_dat(n = 500,
                 X_type = "continous", 
                 r = .4, 
                 ev1_pars = list("a1" = input$shape_1, 
                                 "h1_0" = input$base_1,
                                 "b1" = input$beta_1, 
                                 "gamm1" = input$gamma_1),
                 ev2_pars = list("a2" = input$shape_2, 
                                 "h2_0" = input$base_2, 
                                 "b2" = input$beta_2, 
                                 "gamm2" = input$gamma_2), 
                 rate_cens = input$rate_cens)
    
  })
  
  
  # Show the values in an HTML table ----
  output$plot_CI <- renderPlot({
    
    cumincs_plot_truepred(cox_long = setup_mstate(dat_cumin()), 
                          combo = data.frame(
                            "val_X" = input$X, 
                            "val_Z" = input$Z
                          ),
                          ev1_pars = list(
                            "a1" = input$shape_1, 
                            "h1_0" = input$base_1,
                            "b1" = input$beta_1, 
                            "gamm1" = input$gamma_1
                          ),
                          ev2_pars = list(
                            "a2" = input$shape_2, 
                            "h2_0" = input$base_2, 
                            "b2" = input$beta_2, 
                            "gamm2" = input$gamma_2
                          ), 
                          dat = dat_cumin())
  })
  
  
  
  output$descrip <- renderPrint({
    table(dat_cumin()$eps) / nrow(dat_cumin())
  })
  
  output$plot_haz <- renderPlot({
    
    t <- seq(0.01, 5, by = 0.01)
    
    lam1 <- input$shape_1 * exp((input$beta_1 * input$X + 
                                   input$gamma_1 * input$Z))
    lam2 <- input$shape_2 * exp((input$beta_2 * input$X + 
                                   input$gamma_2 * input$Z))
    
    cumhaz1 <- cumhaz_weib(alph = input$shape_1, lam = lam1, t = t) 
    cumhaz2 <- cumhaz_weib(alph = input$shape_2, lam = lam2, t = t) 
    
    cbind.data.frame(
      "t" = t,
      "haz1" = cumhaz1,
      "haz2" = cumhaz2
    ) %>% 
      gather("event", "cumhaz", .data$haz1, .data$haz2) %>% 
      ggplot(aes(.data$t, .data$cumhaz, col = .data$event)) +
      geom_line(size = 1.5) +
      ggtitle(paste0("Cumulative hazard for X = ", input$X,
                     " and Z = ", input$Z)) +
      xlab("Time") +
      ylab("Cumulative hazard") +
      scale_color_manual(
        "State", 
        values = c(2, 3),
        labels = c("REL", "NRM")
      ) +
      theme(legend.position = "bottom")
    
  })
  
  output$plot_condhaz <- renderPlot({

    t <- seq(0.01, 5, by = 0.01)
    
    lam1 <- input$shape_1 * exp((input$beta_1 * input$X + 
                                   input$gamma_1 * input$Z))
    lam2 <- input$shape_2 * exp((input$beta_2 * input$X + 
                                   input$gamma_2 * input$Z))
    
    haz1 <- haz_weib(alph = input$shape_1, lam = lam1, t = t) 
    haz2 <- haz_weib(alph = input$shape_2, lam = lam2, t = t) 
    
    cbind.data.frame(
      "t" = t,
      "haz1" = haz1,
      "haz2" = haz2
    ) %>% 
      gather("event", "hazard", .data$haz1, .data$haz2) %>% 
      ggplot(aes(.data$t, .data$hazard, col = .data$event)) +
      geom_line(size = 1.5) +
      ggtitle(paste0("Hazard for X = ", input$X,
                     " and Z = ", input$Z)) +
      xlab("Time") +
      ylab("Conditional hazard") +
      scale_color_manual(
        "State", 
        values = c(2, 3),
        labels = c("REL", "NRM")
      ) +
      theme(legend.position = "bottom")
  })
  
  output$rev_KM <- renderPlot({
    
    # Run reverse KM
    fit <- survfit(Surv(t, eps == 0) ~ 1, data = dat_cumin())
    
    cbind.data.frame("t" = fit$time,
                     "surv" = fit$surv) %>% 
      ggplot(aes(.data$t, .data$surv)) +
      geom_line(size = 1.5, linetype = "solid") +
      ggtitle("Reverse Kaplan-Meier") +
      xlab("Time") + ylab("Probability")
    
  })
}

shinyApp(ui = ui, server = server)
