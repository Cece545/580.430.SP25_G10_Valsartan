library(shiny)
library(plotly)
library(dplyr)
library(readr)
library(tidyr)
library(R.matlab)

# Define helper functions to load simulation data based on user input
timecourse_path <- function(a, k) {
  if (k == 0) {
    paste0("missed_dose_", a, "_recovery_0_timeseries.csv")
  } else {
    paste0("missed_dose_", a, "_recovery_", k, "m_5_timeseries.csv")
  }
}

metrics_path <- function(a) {
  paste0("missed_dose_", a, "_metrics.csv")
}

# Load AUC Data for 3D Plot
load_auc_data <- function() {
  case1_data <- read_csv("valsartan_case1_healthy.csv")
  case2_data <- read_csv("valsartan_case2_diseased.csv")
  case1_data$Population <- "Healthy"
  case2_data$Population <- "Diseased"
  list(case1 = case1_data, case2 = case2_data)
}

# Define UI
ui <- fluidPage(
  titlePanel("Valsartan Missed Dose Pharmacokinetics Viewer"),
  sidebarLayout(
    sidebarPanel(
      numericInput("dose_number", "Which dose was missed? (2 to 6):", value = 2, min = 2, max = 6, step = 1),
      numericInput("delay_factor", "How much delay? (k = 0 to 5):", value = 0, min = 0, max = 5, step = 1),
      sliderInput("time_window", "Select Time Range (hrs):", min = 0, max = 168, value = c(0, 168)),
      tags$hr(),
      tags$h4("Figure Caption:"),
      helpText(
        "This app visualizes PK/PD effects of missed or delayed valsartan doses.",
        "You can inspect free valsartan, AngII-Receptor complex concentration over time,",
        "and summary metrics (Cmax, Cmin, AUC).",
        "This app also includes a 3D plot for", 
        "relationship between [AngII-AT1R complex] and weight and baseline [receptor]."
      )
    ),
    mainPanel(
      plotlyOutput("valsartan_plot"),
      plotlyOutput("complex_plot"),
      plotlyOutput("metrics_plot"),
      tags$h4("3D plots: [AngII-AT1R receptor complex] vs. weight and baseline [receptor]",  style = "text-align: center; font-weight: normal;"),
      fluidRow(
        column(6, plotlyOutput("auc_3d_plot_healthy", height = "280px")),
        column(6, plotlyOutput("auc_3d_plot_diseased", height = "280px"))
      ),
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  timeseries_data <- reactive({
    dfs <- list()
    for (k in 0:5) {
      path <- file.path(timecourse_path(input$dose_number, k))
      if (file.exists(path)) {
        df <- read_csv(path, show_col_types = FALSE) %>%
          rename(Time = Time_hr, Valsartan = Free_Valsartan_nM, Complex = AngII_ATR1_Complex_nM) %>%
          filter(Time >= input$time_window[1], Time <= input$time_window[2]) %>%
          mutate(ScenarioK = k)
        dfs[[length(dfs) + 1]] <- df
      }
    }
    all_df <- bind_rows(dfs)
    return(all_df)
  })

  metrics_data <- reactive({
    path <- file.path(metrics_path(input$dose_number))
    req(file.exists(path))
    df <- read_csv(path)

    if (input$delay_factor == 0) {
      scenario_name <- "0"
    } else if (input$delay_factor == 5) {
      scenario_name <- "m"
    } else {
      scenario_name <- c("m/5", "2m/5", "3m/5", "4m/5")[input$delay_factor]
    }

    df_filtered <- df %>% filter(Scenario == scenario_name)
    return(df_filtered)
  })

  output$valsartan_plot <- renderPlotly({
    colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
    df <- timeseries_data()
    if (nrow(df) == 0) {
      return(plot_ly() %>% layout(title = "No data available for selected input"))
    }

    p <- plot_ly()
    scenarios <- sort(unique(df$ScenarioK))

    for (k in scenarios) {
      sub_df <- df %>% filter(ScenarioK == k)
      opacity_val <- ifelse(k == input$delay_factor, 1, 0.3)
      width_val <- ifelse(k == input$delay_factor, 3, 1)

      p <- p %>% add_trace(
        data = sub_df,
        x = ~Time, y = ~Valsartan,
        type = "scatter", mode = "lines",
        name = paste("Delay =", k),
        opacity = opacity_val,
        line = list(width = width_val, color = colors[k + 1]),
        showlegend = TRUE
      )
    }

    p %>% layout(
      title = "Free Valsartan Concentration Over Time",
      xaxis = list(title = "Time (hr)"),
      yaxis = list(title = "[Valsartan] (nM)")
    )
  })

  output$complex_plot <- renderPlotly({
    colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
    df <- timeseries_data()
    if (nrow(df) == 0) {
      return(plot_ly() %>% layout(title = "No data available for selected input"))
    }

    p <- plot_ly()
    scenarios <- sort(unique(df$ScenarioK))

    for (k in scenarios) {
      sub_df <- df %>% filter(ScenarioK == k)
      opacity_val <- ifelse(k == input$delay_factor, 1, 0.3) # Use full opacity for all
      width_val <- ifelse(k == input$delay_factor, 3, 1)

      p <- p %>% add_trace(
        data = sub_df,
        x = ~Time, y = ~Complex,
        type = "scatter", mode = "lines",
        name = paste("Delay =", k),
        opacity = opacity_val,
        line = list(width = width_val, color = colors[k + 1]),
        showlegend = TRUE
      )
    }

    p %>% layout(
      title = "AngII-Receptor Complex Concentration Over Time",
      xaxis = list(title = "Time (hr)"),
      yaxis = list(title = "[AngII-Receptor Complex] (nM)")
    )
  })

  output$metrics_plot <- renderPlotly({
    colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
    path <- file.path(metrics_path(input$dose_number))
    req(file.exists(path))
    df_all <- read_csv(path)

    # Define mapping from Scenario to numeric delay_factor
    scenario_map <- c("0" = 0, "m/5" = 1, "2m/5" = 2, "3m/5" = 3, "4m/5" = 4, "m" = 5)
    df_all <- df_all %>%
      mutate(delay_factor = scenario_map[Scenario]) %>%
      filter(!is.na(delay_factor))

      # Separate plots for AUC, Cmax, Cmin
    plots <- list()
    for (i in seq_along(c("AUC", "Cmax", "Cmin"))) {
      metric <- c("AUC", "Cmax", "Cmin")[i]
    
      df_metric <- df_all %>% 
        select(delay_factor, !!sym(metric)) %>%
        rename(Value = !!sym(metric)) %>%
        mutate(opacity = ifelse(delay_factor == input$delay_factor, 1, 0.3))
    
      # Show legend only for the first plot (AUC), suppress for Cmax and Cmin
      show_legend <- ifelse(i == 1, TRUE, FALSE)
    
      p <- plot_ly(
        data = df_metric,
        x = ~delay_factor,
        y = ~Value,
        type = "bar",
        color = ~factor(delay_factor, labels = c("Delay = 0", "Delay = 1", "Delay = 2", "Delay = 3", "Delay = 4", "Delay = 5")),
        colors = colors,
        opacity = ~opacity,
        text = ~round(Value, 2),
        textposition = "auto",
        showlegend = show_legend
      ) %>% layout(
        title = NULL,
        yaxis = list(title = ifelse(metric == "AUC", "AUC (nM*hr)", paste(metric, "(nM)"))),
        xaxis = list(
          title = "Delayed Time (hr)",
          tickvals = c(0, 1, 2, 3, 4, 5),
          ticktext = c("0", "4.8", "9.6", "14.4", "19.2", "24")
        )
      )

      plots[[length(plots) + 1]] <- p
    }

    # Combine subplots and adjust layout
    subplot(plots, nrows = 1, shareX = FALSE, titleX = TRUE, widths = c(0.25, 0.3, 0.25), margin = 0.05) %>% layout(
      yaxis = list(title = paste("AUC of [AngII-AT1R complex] (nM*hr) for Dose", input$dose_number), titlefont = list(size = 12)),
      yaxis2 = list(title = paste("Cmax of [AngII-AT1R complex] (nM) for Dose", input$dose_number), titlefont = list(size = 12)),
      yaxis3 = list(title = paste("Cmin of [AngII-AT1R complex] (nM) for Dose", input$dose_number), titlefont = list(size = 12)),
      title = "Metrics plot",
      legend = list(
        x = 1, y = 1,
        orientation = 'v',
        traceorder = 'normal',
        itemsizing = 'constant',
        font = list(size = 12)
      )
    )

  })
  # Load AUC Data
  auc_data <- reactive({
    load_auc_data()
  })

# Determine Axis Ranges
  y_range <- reactive({
    data <- auc_data()
    range(c(data$case1$Receptor_Healthy_uM, data$case2$Receptor_Diseased_uM))
  })

  z_range <- reactive({
    data <- auc_data()
    range(c(data$case1$AUC_AngII_AT1R_Healthy, data$case2$AUC_AngII_AT1R_Diseased))
  })

  output$auc_3d_plot_healthy <- renderPlotly({
    data <- auc_data()$case1
    plot_ly(
      data = data,
      x = ~Weight_kg, y = ~Receptor_Healthy_uM, z = ~AUC_AngII_AT1R_Healthy,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 4, color = "#377EB8")
    ) %>% layout(
      scene = list(
        xaxis = list(title = "Weight (kg)"),
        yaxis = list(title = "Receptor (uM)", range = y_range()),
        zaxis = list(title = "AUC (nM*hr)", range = z_range()),
        bgcolor = 'rgba(0, 0, 0, 0)'
      ),
      title = "Healthy Population"
    )
  })

  output$auc_3d_plot_diseased <- renderPlotly({
    data <- auc_data()$case2
    plot_ly(
      data = data,
      x = ~Weight_kg, y = ~Receptor_Diseased_uM, z = ~AUC_AngII_AT1R_Diseased,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 4, color = "#E41A1C")
    ) %>% layout(
      scene = list(
        xaxis = list(title = "Weight (kg)"),
        yaxis = list(title = "Receptor (uM)", range = y_range()),
        zaxis = list(title = "AUC (nM*hr)", range = z_range()),
        bgcolor = 'rgba(0, 0, 0, 0)'
      ),
      title = "Diseased Population"
    )
  })
}

# Run the app
shinyApp(ui = ui, server = server)
