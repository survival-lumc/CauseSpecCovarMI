# Load packages - MASS before tidyverse because of select
pacman::p_load(MASS, cmprsk, survival, 
               tidyverse, mice, shiny)

source("data_generation/dat_generation_MRR_KMweib.R")





# Make this into app
ui <- fluidPage(
  
  # App title ----
  titlePanel("Parameters"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar to demonstrate various slider options ----
    sidebarPanel(
      
      # Input:  ----
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
      
      # Input:  ----
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
                  value = 0, step = 0.01),
      
      sliderInput("rate_cens", "Exp censoring rate:",
                  min = 0, max = 5,
                  value = 0.5, step = 0.01)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(width = 6,
      tabsetPanel(type = "tabs",
        tabPanel(title = "Cumulative incidence",
                 plotOutput("plot_CI")),
        tabPanel(title = "Theoretical CI",
                 plotOutput("theor_CI")),
        tabPanel(title = "Hazards",
                 plotOutput("plot_haz")),
        tabPanel(title = "Density time",
                 plotOutput("plot_denstime")),
        tabPanel(title = "Condit. Haz X = Z = 0",
                 plotOutput("plot_condhaz")),
        tabPanel(title = "Prop Events",
                 verbatimTextOutput("descrip"))
      )
    )
  )
)

server <- function(input, output) {
  
  # Reactive expression to create data frame of all input values ----
  dat_cumin <- reactive({
    
    dat <- dat_gener_MRR_KM(N = 50000,
                            X_type = "contin",
                            mus = c(0, 0), 
                            covmat = matrix(c(1, 0.25, 
                                              0.25, 1), nrow = 2), 
                            mech = "MCAR", 
                            p = 0, 
                            cause2 = "weib", 
                            vals_t1 = c(input$shape_1, input$base_1, 
                                        input$beta_1, input$gamma_1), 
                            vals_t2 = c(input$shape_2, input$base_2, 
                                        input$beta_2, input$gamma_2),
                            rate_cens = input$rate_cens)
    
  })
  
  
  # Show the values in an HTML table ----
  output$plot_CI <- renderPlot({
    CI <- cuminc(ftime = dat_cumin()$t, fstatus = dat_cumin()$eps)
    plot(CI, curvlab = c("Cause 1", "Cause 2"), 
         color = c("blue", "black"), lwd = c(3, 3), xlim = c(0, 5))
  })
  
  output$theor_CI <- renderPlot({
    
    # Pick 500 times points so integration is shorted
    t <- dat_cumin()$t
    t <- sort(c(min(t), sample(t, size = 500, replace = F), max(t)))
    
    ci1 <- cuminc_weib(alph_ev = input$shape_1, lam_ev = input$base_1, 
                       alph_comp = input$shape_2, lam_comp = input$base_2,
                       t = t)
    
    ci2 <- cuminc_weib(alph_ev = input$shape_2, lam_ev = input$base_2, 
                       alph_comp = input$shape_1, lam_comp = input$base_1,
                       t = t)
    
    plot(t, ci1, col = "blue", type = "l", ylim = c(0, 1), lwd = 3,
         xlim = c(0, 5))
    lines(t, ci2, col = "black", lwd = 3, lty = 2)
    legend("topleft", c("Cause 1", "Cause 2"), col = c("blue", "black"),
           lty = c(1, 2), lwd = c(3, 3))
  })
  
  output$plot_denstime <- renderPlot({
    dens1 <- dat_cumin() %>% 
      filter(eps == 1) %>% 
      select(t)
    
    dens2 <- dat_cumin() %>% 
      filter(eps == 2) %>% 
      select(t)
    
    plot(density(dens1$t), col = "blue", xlim = c(0, 5), lwd = 3) 
    lines(density(dens2$t), col = "black", lwd = 3)
    legend("topright", c("Cause 1", "Cause 2"), col = c("blue", "black"),
           lty = c(1, 1), lwd = c(3, 3))
  })
  
  output$descrip <- renderPrint({
    table(dat_cumin()$eps) / nrow(dat_cumin())
  })
  
  output$plot_haz <- renderPlot({
    dat <- dat_cumin()
    
    plot(dat$t, dat$H1, type = "l", 
         ylim = c(0, max(dat$H1, dat$H2)), col = "blue",
         xlim = c(0, 5), lwd = 3)
    lines(dat$t, dat$H2, lwd = 3)
    legend("topleft", c("Cause 1", "Cause 2"), col = c("blue", "black"),
           lty = c(1, 1), lwd = c(3, 3))
  })
  
  output$plot_condhaz <- renderPlot({
    dat <- dat_cumin()
    
    haz1 <- haz_weib(alph = input$shape_1, lam = input$base_1, t = dat$t) 
    haz2 <- haz_weib(alph = input$shape_2, lam = input$base_2, t = dat$t) 
    
    plot(dat$t, haz1, type = "l", 
         ylim = c(0, max(haz1, haz2)), col = "blue",
         xlim = c(0, 5), lwd = 3)
    lines(dat$t, haz2, lwd = 3)
    legend("topleft", c("Cause 1", "Cause 2"), col = c("blue", "black"),
           lty = c(1, 1), lwd = c(3, 3))
  })
}

shinyApp(ui = ui, server = server)
