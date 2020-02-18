    #######################
## ui file table app ##
#######################


# Source help files
source("load_packages.R")

# Load/install packages
load_packages(c("shiny", "DT", "shinyWidgets", "shinydashboard", "tidyverse"))

# Load other useful objects
source("support_vars.R")


# Define UI for application that draws a histogram
shinyUI(
    fluidPage(
        # Application title
        titlePanel("Simulation results"),
        
        # Sidebar with a slider input for number of bins
        sidebarLayout(
            sidebarPanel(width = 3,
                         
                         # Button to run based on inputs
                         selectInput("select_scen", 
                                     h4("Select scenario"), 
                                     choices = scen_names,
                                     selectize = T),
                         
                         # Choice of m - uncheck all to have everything
                         checkboxGroupInput("m_imps", 
                                            h4("Choice of m"), 
                                            choices = list("m = 1" = 1, 
                                                           "m = 10" = 10, 
                                                           "m = 20" = 20),
                                            selected = 20),
                         
                         # Columns to display
                         checkboxGroupInput("show_vars", 
                                            h4("Variables to view:"),
                                            names(final_dat)[-1], 
                                            selected = c("var", "analy", "coef",
                                                         "se", "true")),
                         
                         selectInput("select_scen_comp", 
                                     h4("Select scenario to compare"), 
                                     choices = scen_names,
                                     selectize = T),
                         
                         # Run it
                         actionButton("apply", "Apply")
                         
            ),
            
            # Show a plot of the generated distribution
            mainPanel(
                fluidRow(
                    tabBox(
                        width = 13,
                        tabPanel(title = "Table - single",
                                 dataTableOutput("simulations_table")
                                 ),
                        tabPanel(title = "Table - compare",
                                 fluidRow(
                                     box(title = h4(uiOutput("title_scen")),
                                         width = 6,
                                         dataTableOutput("sim1")
                                     ),
                                     box(title = h4(uiOutput("title_comp")),
                                         width = 6,
                                         dataTableOutput("simulations_table_comp")
                                     )
                                 )
                        ),
                        tabPanel(title = "Plot - single",
                                 plotOutput("plot")
                                 ),
                        tabPanel(title = "Plot - compare",
                                 plotOutput("plot_compare1"),
                                 plotOutput("plot_compare2")
                                 )
                    )
                )
            )
        )
    )
)
