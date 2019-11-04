###########################
## Server file table app ##
###########################


# Load
source("support_vars.R")

# Set ggplot theme
theme_set(theme_MIsim_app())


# Define server function
shinyServer(function(input, output) {
    
    # Read in results
    load("final_dat.RData") 
    
    # Throw away seed
    final_dat <- final_dat %>% 
        select(-seed) %>%
        mutate_if(is.numeric, ~ round(., 3))
    
    dat <- eventReactive(input$apply, {
        
        if (is.null(input$m_imps)) {
            tabs <- final_dat %>% 
                filter(scen_name == input$select_scen)
        } else {
            tabs <- final_dat %>% 
                filter(m %in% c(input$m_imps, 0),
                       scen_name == input$select_scen)
        }
        
        tabs <- tabs %>% 
            select(-scen_name) %>% 
            
            # Order columns
            select(input$show_vars)
        
    })
    
    
    dat_comp <- eventReactive(input$apply, {
        
        if (is.null(input$m_imps)) {
            tabs <- final_dat %>% 
                filter(scen_name == input$select_scen_comp)
        } else {
            tabs <- final_dat %>% 
                filter(m %in% c(input$m_imps, 0),
                       scen_name == input$select_scen_comp)
        }
        
        tabs <- tabs %>% 
            select(-scen_name) %>% 
            
            # Order columns
            select(input$show_vars)
            
            # Remove var name and analysis in second table
    })
    
    

    output$simulations_table <- renderDataTable({ 
        datatable(dat(), options = list(
            pageLength = 28
            
        ),
        rownames = FALSE) %>% 
            formatStyle(
                'var',
                target = "row",
                fontWeight = styleEqual(c("X1", "X2"), c('', 'bold'))
            )
    }) 
    
    output$sim1 <- renderDataTable({ 
        datatable(dat(), options = list(
            pageLength = 28
            
        ),
        rownames = FALSE) %>% 
            formatStyle(
                'var',
                target = "row",
                fontWeight = styleEqual(c("X1", "X2"), c('', 'bold'))
            )
    }) 
    
    # For comparison
    output$simulations_table_comp <- renderDataTable({ 
        datatable(dat_comp(), options = list(
            pageLength = 28
            
        ),
        rownames = FALSE) %>% 
            formatStyle(
                'var',
                target = "row",
                fontWeight = styleEqual(c("X1", "X2"), c('', 'bold'))
            )
    })
    
    # Make titles for comparative plots
    output$title_scen <- renderText({
        HTML(paste0("<b>", scen_labels[which(scen_names == input$select_scen)]))
    })
    
    output$title_comp <- renderText({
        HTML(paste0("<b>", scen_labels[which(scen_names == input$select_scen_comp)]))
    })
    
    output$plot <- renderPlot({
        final_dat %>%
            filter(scen_name == input$select_scen,
                   m %in% c(input$m_imps, 0)) %>% 
            ggplot(aes(analy, coef)) +
            geom_hline(aes(yintercept = true), linetype = "dashed") +
            geom_errorbar(aes(ymin = coef - 1.96 * se_emp,
                              ymax = coef + 1.96 * se_emp,
                              colour = var), 
                          width = 0.2,
                          size = 1) + 
            geom_point(size = 2) +
            coord_flip() +
            ggtitle(scen_labels[which(scen_names == input$select_scen)]) +
            facet_grid(. ~ var) +
            scale_y_continuous(limits = c(-.75, 1.25), breaks = c(-.5, 0, .5, 1)) +
            xlab("Analysis") + ylab("Coefficient value") +
            theme(legend.position = "none")
            
    })
    
    output$plot_compare1 <- renderPlot({
        final_dat %>%
            filter(scen_name == input$select_scen,
                   m %in% c(input$m_imps, 0)) %>% 
            ggplot(aes(analy, coef)) +
            geom_hline(aes(yintercept = true), linetype = "dashed") +
            geom_errorbar(aes(ymin = coef - 1.96 * se_emp,
                              ymax = coef + 1.96 * se_emp,
                              colour = var), width = 0.2,
                          size = 1) + 
            geom_point(size = 2) +
            coord_flip() +
            ggtitle(scen_labels[which(scen_names == input$select_scen)]) +
            facet_grid(. ~ var, scales = "fixed") +
            scale_y_continuous(limits = c(-.75, 1.25), breaks = c(-.5, 0, .5, 1)) +
            xlab("Analysis") + ylab("Coefficient value") 
    })
    
    output$plot_compare2 <- renderPlot({
        final_dat %>%
            filter(scen_name == input$select_scen_comp,
                   m %in% c(input$m_imps, 0)) %>%  
            ggplot(aes(analy, coef)) +
            geom_hline(aes(yintercept = true), linetype = "dashed") +
            geom_errorbar(aes(ymin = coef - 1.96 * se_emp,
                              ymax = coef + 1.96 * se_emp,
                              colour = var), width = 0.2,
                          size = 1) + 
            geom_point(size = 2) +
            coord_flip() +
            ggtitle(scen_labels[which(scen_names == input$select_scen_comp)]) +
            facet_grid(. ~ var, scales = "fixed") +
            scale_y_continuous(limits = c(-.75, 1.25), breaks = c(-.5, 0, .5, 1)) +
            xlab("Analysis") + ylab("Coefficient value") 
    })
})


