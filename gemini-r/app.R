# Load necessary library
# Ensure you have 'shiny' installed: install.packages("shiny")
library(shiny)
# library(scales) # Used for alpha transparency - install if needed: install.packages("scales")

# --- Discrete Event Simulation Model Code ---
# Includes parameters, helper functions, and simulation logic.

# --- Simulation Parameters (Inspired by Python Example) ---
# Define the parameters governing the simulation model.
# These include time horizon, discount rate, patient numbers,
# event rates (stroke, mortality), quality of life utilities, and costs.
# Default parameters are set here.
params <- list(
  horizon = 30,                    # Simulation horizon in years
  discount_rate = 0.035,           # Annual discount rate for costs and QALYs
  n_patients = 1000,               # Number of simulated patients per arm (can be adjusted in UI)
  start_age = 50,                  # Starting age of patients

  # Event Rates (lambda for exponential distribution - events per year)
  baseline_stroke_rate = 0.05,     # Annual stroke risk for dupimab (in Healthy state)
  hr_stroke_supimab = 0.7,         # Hazard Ratio for stroke on Supimab vs Dupimab
  mortality_rate_healthy = 0.02,   # Background annual mortality rate (example, could be age-dependent)
  mortality_rate_post_stroke = 0.1,# Increased annual mortality rate after stroke

  # Utilities (Quality of Life weights)
  utility_healthy = 0.85,          # Utility QoL weight for Healthy state
  utility_post_stroke = 0.70,      # Utility QoL weight for Post-Stroke state

  # Costs (£)
  drug_cost_dupimab = 10000,       # Annual cost for dupimab
  drug_cost_supimab = 15000,       # Annual cost for supimab
  stroke_cost = 20000,             # One-off cost at stroke event
  annual_poststroke_cost = 5000    # Additional annual cost after stroke
)


# --- Function to Simulate One Patient ---
# This function simulates the life history of a single patient through the model.
# It tracks health state transitions, event occurrences, and accumulates costs and QALYs.
simulate_patient_r <- function(patient_id, treatment, current_params) {
  # Use current_params passed to the function, allowing UI adjustments
  # Initialize patient state variables
  current_time <- 0
  current_age <- current_params$start_age
  health_state <- "Healthy" # Possible states: Healthy, Post-Stroke, Dead
  total_discounted_cost <- 0
  total_discounted_qaly <- 0
  time_of_stroke <- Inf # Use Inf to indicate stroke hasn't happened yet

  # Determine treatment-specific parameters based on the assigned treatment arm
  if (treatment == "Supimab") {
    stroke_rate <- current_params$baseline_stroke_rate * current_params$hr_stroke_supimab
    drug_cost <- current_params$drug_cost_supimab
  } else { # Assuming Dupimab is the alternative
    stroke_rate <- current_params$baseline_stroke_rate
    drug_cost <- current_params$drug_cost_dupimab
  }

  # --- Simulation Loop (Event-driven) ---
  # The loop continues as long as the patient is not dead and the simulation time is within the horizon.
  while (current_time < current_params$horizon && health_state != "Dead") {

    # Determine current event rates based on the patient's current health state
    # Stroke can only occur if the patient is currently in the 'Healthy' state.
    current_stroke_rate <- ifelse(health_state == "Healthy", stroke_rate, 0)
    # Mortality rate depends on whether the patient has had a stroke.
    current_mortality_rate <- ifelse(health_state == "Healthy",
                                     current_params$mortality_rate_healthy,
                                     current_params$mortality_rate_post_stroke)

    # Sample time to next potential events from the *current* time using exponential distribution.
    # If a rate is 0, the time to that event is effectively infinity.
    time_to_next_stroke <- ifelse(current_stroke_rate > 0, rexp(1, rate = current_stroke_rate), Inf)
    time_to_next_death <- rexp(1, rate = current_mortality_rate)

    # Determine the time until the *next* event occurs (the minimum of the sampled times)
    time_to_next_event <- min(time_to_next_stroke, time_to_next_death)

    # Determine the time step for the current interval.
    # It's the time to the next event, unless that pushes the simulation beyond the horizon.
    time_step <- min(time_to_next_event, current_params$horizon - current_time)

    # Safety check: If time_step is effectively zero or negative, break the loop.
    if (time_step <= 1e-9) {
      # print(paste("Warning: Zero or negative time step for patient", patient_id, "at time", current_time))
      break
    }

    # --- Accrue Costs and QALYs for the Interval ---
    # Determine the utility and ongoing costs associated with the patient's state *during* this interval.
    current_utility <- ifelse(health_state == "Healthy", current_params$utility_healthy, current_params$utility_post_stroke)
    current_annual_ongoing_cost <- drug_cost + ifelse(health_state == "Post-Stroke", current_params$annual_poststroke_cost, 0)

    # Calculate discounted QALYs and Costs accrued during this time_step using continuous discounting integral.
    if (current_params$discount_rate > 0) {
        discount_integral_factor <- (1 / current_params$discount_rate) *
                                    (exp(-current_params$discount_rate * current_time) - exp(-current_params$discount_rate * (current_time + time_step)))
    } else { # Handle the case of zero discount rate (no discounting)
        discount_integral_factor <- time_step
    }

    interval_qaly <- current_utility * discount_integral_factor
    interval_ongoing_cost <- current_annual_ongoing_cost * discount_integral_factor

    # Add interval results to patient totals
    total_discounted_qaly <- total_discounted_qaly + interval_qaly
    total_discounted_cost <- total_discounted_cost + interval_ongoing_cost

    # --- Advance Time and Process Event ---
    # Calculate the absolute time at which the next event *would* occur if horizon allowed
    next_event_absolute_time <- current_time + time_to_next_event

    # Advance the simulation clock to the end of the current interval
    current_time <- current_time + time_step
    current_age <- current_age + time_step # Update age as well

    # --- Process the Event that Occurred (if it happened before the horizon) ---
    # Check if an event is scheduled to happen exactly at the new current_time
    # (meaning the time_step was determined by time_to_next_event, not the horizon)
    if (abs(current_time - next_event_absolute_time) < 1e-9 && current_time < current_params$horizon) {
        # Use a small tolerance (epsilon) for comparing floating point numbers
        epsilon <- 1e-9
        if (abs(time_to_next_stroke - time_to_next_event) < epsilon && time_to_next_stroke != Inf) {
            # Stroke event occurred
            health_state <- "Post-Stroke"
            time_of_stroke <- current_time
            # Add the discounted one-off cost of the stroke event (discounted to the time of event)
            total_discounted_cost <- total_discounted_cost + current_params$stroke_cost * exp(-current_params$discount_rate * current_time)

        } else if (abs(time_to_next_death - time_to_next_event) < epsilon) {
            # Death event occurred
            health_state <- "Dead"
            # The while loop condition (health_state != "Dead") will now cause termination
        }
        # If neither stroke nor death was the minimum time, it implies an issue (e.g., rates were zero).
        # The loop continues if state is not Dead and time < horizon.
    }
    # If current_time >= params$horizon, the loop terminates based on the while condition.

  } # End while loop for the patient

  # Return a list containing the final results for this patient
  return(list(
    patient_id = patient_id,
    treatment = treatment,
    total_cost = total_discounted_cost,
    total_qaly = total_discounted_qaly,
    time_of_stroke = ifelse(is.infinite(time_of_stroke), NA, time_of_stroke), # Store NA if no stroke
    time_of_death = ifelse(health_state == "Dead", current_time, NA)          # Store time of death or NA
  ))
}

# --- Main Simulation Runner Function (Modified for Shiny) ---
# This function orchestrates the simulation for both treatment arms and processes results.
run_des_simulation_shiny <- function(current_params) {
  # Initialize an empty list to store results from individual patient simulations
  results_list <- list()
  patient_counter <- 1 # Counter to index the results_list

  # Use shiny::withProgress to show simulation progress in the UI
  withProgress(message = 'Running Simulation', value = 0, {
      total_sims = 2 * current_params$n_patients

      # --- Simulate Dupimab Arm ---
      for (i in 1:current_params$n_patients) {
        results_list[[patient_counter]] <- simulate_patient_r(i, "Dupimab", current_params)
        patient_counter <- patient_counter + 1
        # Update progress bar
        incProgress(1/total_sims, detail = paste("Dupimab Patient", i))
      }

      # --- Simulate Supimab Arm ---
      for (i in 1:current_params$n_patients) {
        results_list[[patient_counter]] <- simulate_patient_r(i, "Supimab", current_params)
        patient_counter <- patient_counter + 1
        # Update progress bar
        incProgress(1/total_sims, detail = paste("Supimab Patient", i))
      }
  }) # End withProgress

  # --- Process Results ---
  # Convert the list of lists into a data frame.
  col_names <- names(results_list[[1]])
  results_matrix <- matrix(unlist(results_list), ncol = length(col_names), byrow = TRUE)
  colnames(results_matrix) <- col_names
  results_df <- data.frame(results_matrix, stringsAsFactors = FALSE)

  # Convert relevant columns back to numeric type.
  numeric_cols <- c("total_cost", "total_qaly", "time_of_stroke", "time_of_death")
  for (col in numeric_cols) {
     results_df[[col]] <- suppressWarnings(as.numeric(as.character(results_df[[col]])))
  }
  results_df$patient_id <- as.integer(as.character(results_df$patient_id))

  # --- Aggregate Results ---
  # Calculate mean cost and QALY per treatment arm using tapply.
  mean_costs <- tapply(results_df$total_cost, results_df$treatment, mean, na.rm = TRUE)
  mean_qalys <- tapply(results_df$total_qaly, results_df$treatment, mean, na.rm = TRUE)

  # Extract the calculated means, handling potential NAs if a treatment arm had no valid results
  mean_cost_dupimab <- ifelse("Dupimab" %in% names(mean_costs), mean_costs["Dupimab"], NA)
  mean_qaly_dupimab <- ifelse("Dupimab" %in% names(mean_qalys), mean_qalys["Dupimab"], NA)
  mean_cost_supimab <- ifelse("Supimab" %in% names(mean_costs), mean_costs["Supimab"], NA)
  mean_qaly_supimab <- ifelse("Supimab" %in% names(mean_qalys), mean_qalys["Supimab"], NA)

  # --- Calculate Incremental Results & ICER ---
  # Ensure means are valid numbers before calculating differences
  if (is.na(mean_cost_supimab) || is.na(mean_cost_dupimab) || is.na(mean_qaly_supimab) || is.na(mean_qaly_dupimab)) {
      inc_cost <- NA
      inc_qaly <- NA
      icer <- NA
      icer_text <- "N/A (Missing Data)"
  } else {
      inc_cost <- mean_cost_supimab - mean_cost_dupimab
      inc_qaly <- mean_qaly_supimab - mean_qaly_dupimab

      # Calculate the ICER text based on dominance or calculated value
      if (inc_qaly <= 1e-9) { # Use tolerance for floating point comparison
          icer <- NA
          icer_text <- ifelse(inc_cost < 0, "Supimab Dominant", "Supimab Dominated or Equal")
      } else {
          icer <- inc_cost / inc_qaly
          icer_text <- sprintf("£%.0f per QALY", icer)
      }
  }

  # Return a list containing the summary statistics and the full patient-level data frame
  return(list(
    summary = list(
      Dupimab = list(mean_cost = mean_cost_dupimab, mean_qaly = mean_qaly_dupimab),
      Supimab = list(mean_cost = mean_cost_supimab, mean_qaly = mean_qaly_supimab),
      Incremental = list(cost = inc_cost, qaly = inc_qaly, icer = icer, icer_text = icer_text)
    ),
    patient_data = results_df
  ))
}


# --- Shiny UI Definition ---
# Defines the user interface layout and elements.
ui <- fluidPage(
  theme = shinythemes::shinytheme("lumen"), # Apply a theme for better appearance
  titlePanel("DES Model: Supimab vs Dupimab Cost-Effectiveness"),
  sidebarLayout(
    sidebarPanel(
      width = 3, # Adjust sidebar width
      h4("Simulation Controls"),
      # Input for number of patients
      numericInput("n_patients_input", "Number of Patients per Arm:",
                   value = params$n_patients, min = 100, max = 10000, step = 100),
       # Input for discount rate
      sliderInput("discount_rate_input", "Annual Discount Rate (%):",
                  min = 0, max = 10, value = params$discount_rate * 100, step = 0.1),
      # Add more inputs here later if needed (e.g., costs, HR)
      # --- Action Button to trigger simulation ---
      actionButton("run_sim", "Run Simulation", icon = icon("play"), class = "btn-primary btn-lg"), # Style button
      hr(),
      # --- Display Summary Results ---
      h4("Summary Results"),
      tags$div(style="font-size: 0.9em; line-height: 1.6;", # Style summary text
          tags$strong("Dupimab:"),
          textOutput("summary_dup_cost"),
          textOutput("summary_dup_qaly"),
          tags$br(),
          tags$strong("Supimab:"),
          textOutput("summary_sup_cost"),
          textOutput("summary_sup_qaly"),
          tags$br(),
          tags$strong("Incremental (Supimab vs Dupimab):"),
          textOutput("summary_inc_cost"),
          textOutput("summary_inc_qaly"),
          tags$strong(textOutput("summary_icer")) # Make ICER bold
      )
    ),
    # --- Main Panel for Plots ---
    mainPanel(
      width = 9, # Adjust main panel width
      tabsetPanel( # Use tabs for better organization
          tabPanel("Summary Plots",
                   fluidRow(
                       column(6, h4("Mean Cost Comparison"), plotOutput("plot_mean_cost", height = "350px")),
                       column(6, h4("Mean QALY Comparison"), plotOutput("plot_mean_qaly", height = "350px"))
                   ),
                   fluidRow(
                       column(12, h4("Cost-Effectiveness Plane"), plotOutput("plot_ce_plane", height = "500px"))
                   )
          ),
          tabPanel("Distributions",
                   fluidRow(
                       column(6, h4("Cost Distribution (Dupimab)"), plotOutput("plot_hist_cost_dup", height = "300px")),
                       column(6, h4("Cost Distribution (Supimab)"), plotOutput("plot_hist_cost_sup", height = "300px"))
                   ),
                   fluidRow(
                       column(6, h4("QALY Distribution (Dupimab)"), plotOutput("plot_hist_qaly_dup", height = "300px")),
                       column(6, h4("QALY Distribution (Supimab)"), plotOutput("plot_hist_qaly_sup", height = "300px"))
                   )
          ),
          tabPanel("Patient Data",
                   h4("Raw Patient-Level Results"),
                   DT::dataTableOutput("patient_table") # Use DT for interactive table
          )
      )
    )
  )
)

# --- Shiny Server Logic ---
# Defines the reactive logic that runs the simulation and generates outputs.
server <- function(input, output, session) {

  # Reactive value to store simulation results (initially NULL)
  sim_results <- reactiveVal(NULL)

  # --- Run Simulation on Button Click ---
  # observeEvent listens for the button click.
  observeEvent(input$run_sim, {
    # Update parameters based on current UI inputs
    current_sim_params <- params # Start with default params
    current_sim_params$n_patients <- as.integer(input$n_patients_input) # Ensure integer
    current_sim_params$discount_rate <- input$discount_rate_input / 100 # Convert percentage

    # Input validation (optional but good practice)
    if(is.na(current_sim_params$n_patients) || current_sim_params$n_patients < 10) {
        showNotification("Please enter a valid number of patients (>= 10).", type = "error")
        return() # Stop execution
    }
    if(is.na(current_sim_params$discount_rate) || current_sim_params$discount_rate < 0) {
        showNotification("Please enter a valid discount rate (>= 0).", type = "error")
        return() # Stop execution
    }

    # Run the simulation function (defined above)
    # This function now includes the withProgress indicator
    results <- run_des_simulation_shiny(current_sim_params)

    # Store results in the reactive value, triggering updates in outputs
    sim_results(results)

    # Confirmation notification
    showNotification("Simulation complete!", type = "message", duration = 5)
  })

  # --- Render Summary Text Outputs ---
  # These renderText functions update the text in the sidebar when sim_results() changes.
  # req(sim_results()) ensures that the code only runs after results are available.
  output$summary_dup_cost <- renderText({
    req(sim_results())
    res_sum <- sim_results()$summary$Dupimab
    sprintf("Cost = £%.0f", ifelse(is.na(res_sum$mean_cost), 0, res_sum$mean_cost))
  })
  output$summary_dup_qaly <- renderText({
    req(sim_results())
    res_sum <- sim_results()$summary$Dupimab
    sprintf("QALYs = %.4f", ifelse(is.na(res_sum$mean_qaly), 0, res_sum$mean_qaly))
  })
   output$summary_sup_cost <- renderText({
    req(sim_results())
    res_sum <- sim_results()$summary$Supimab
    sprintf("Cost = £%.0f", ifelse(is.na(res_sum$mean_cost), 0, res_sum$mean_cost))
  })
  output$summary_sup_qaly <- renderText({
    req(sim_results())
     res_sum <- sim_results()$summary$Supimab
    sprintf("QALYs = %.4f", ifelse(is.na(res_sum$mean_qaly), 0, res_sum$mean_qaly))
  })
   output$summary_inc_cost <- renderText({
    req(sim_results())
    res_inc <- sim_results()$summary$Incremental
    sprintf("Incr. Cost = £%.0f", ifelse(is.na(res_inc$cost), 0, res_inc$cost))
  })
  output$summary_inc_qaly <- renderText({
    req(sim_results())
    res_inc <- sim_results()$summary$Incremental
    sprintf("Incr. QALYs = %.4f", ifelse(is.na(res_inc$qaly), 0, res_inc$qaly))
  })
   output$summary_icer <- renderText({
    req(sim_results())
    sprintf("ICER = %s", sim_results()$summary$Incremental$icer_text)
  })

  # --- Render Plots (using Base R) ---
  # These renderPlot functions generate the plots when sim_results() changes.

  # Bar Plot for Mean Costs
  output$plot_mean_cost <- renderPlot({
    req(sim_results())
    res_summary <- sim_results()$summary
    costs <- c(Dupimab = res_summary$Dupimab$mean_cost, Supimab = res_summary$Supimab$mean_cost)
    costs[is.na(costs) | !is.finite(costs)] <- 0 # Handle NA/Inf
    max_y <- max(costs, 0, na.rm = TRUE) * 1.15 # Ensure max_y is at least 0, add buffer

    # Set plot margins
    par(mar = c(4, 5, 2, 1)) # bottom, left, top, right

    bp <- barplot(costs,
            names.arg = c("Dupimab", "Supimab"),
            col = c("#a6cee3", "#fdbf6f"), # ColorBrewer Paired adjusted
            main = NULL, # Remove main title from plot itself
            ylab = "Cost (£)",
            ylim = c(0, max_y),
            las = 1, # Rotate y-axis labels
            cex.names = 1.1,
            cex.axis = 0.9,
            cex.lab = 1.1,
            border = "grey20" # Darker border
           )
     # Add text labels on bars
     text(x = bp, y = costs, label = sprintf("£%.0f", costs), pos = 3, cex = 0.9, col = "black", font=2)
     grid(nx=NA, ny=NULL, col="lightgray", lty="dotted") # Add horizontal grid lines
     box(col = "grey60") # Add box around plot
  })

  # Bar Plot for Mean QALYs
  output$plot_mean_qaly <- renderPlot({
    req(sim_results())
    res_summary <- sim_results()$summary
    qalys <- c(Dupimab = res_summary$Dupimab$mean_qaly, Supimab = res_summary$Supimab$mean_qaly)
    qalys[is.na(qalys) | !is.finite(qalys)] <- 0 # Handle NA/Inf
    max_y <- max(qalys, 0, na.rm = TRUE) * 1.15 # Ensure max_y is at least 0, add buffer

    par(mar = c(4, 5, 2, 1)) # bottom, left, top, right

    bp <- barplot(qalys,
            names.arg = c("Dupimab", "Supimab"),
            col = c("#a6cee3", "#fdbf6f"),
            main = NULL,
            ylab = "QALYs",
            ylim = c(0, max_y),
            las = 1,
            cex.names = 1.1,
            cex.axis = 0.9,
            cex.lab = 1.1,
            border = "grey20"
           )
     # Add text labels on bars
     text(x = bp, y = qalys, label = sprintf("%.3f", qalys), pos = 3, cex = 0.9, col = "black", font=2)
     grid(nx=NA, ny=NULL, col="lightgray", lty="dotted")
     box(col = "grey60")
  })

  # Cost-Effectiveness Plane
  output$plot_ce_plane <- renderPlot({
    req(sim_results())
    df <- sim_results()$patient_data
    # Filter out potential NA/Inf values that would cause plotting errors
    df_plot <- df[is.finite(df$total_cost) & is.finite(df$total_qaly), ]
    if (nrow(df_plot) == 0) return(NULL) # Don't plot if no valid data

    # Define colors and symbols
    plot_cols <- ifelse(df_plot$treatment == "Dupimab", "#a6cee3", "#fdbf6f")
    plot_pch <- ifelse(df_plot$treatment == "Dupimab", 16, 17) # Solid circle, solid triangle

    # Determine plot limits dynamically
    xlim_range <- range(df_plot$total_qaly, na.rm = TRUE, finite = TRUE)
    ylim_range <- range(df_plot$total_cost, na.rm = TRUE, finite = TRUE)
    xlim_buffer <- (xlim_range[2] - xlim_range[1]) * 0.05
    ylim_buffer <- (ylim_range[2] - ylim_range[1]) * 0.05
    
    # Ensure limits are sensible if range is zero
    if(xlim_buffer == 0) xlim_buffer = 1 
    if(ylim_buffer == 0) ylim_buffer = 1000

    par(mar = c(5, 5, 2, 2)) # Adjust margins

    plot(df_plot$total_qaly, df_plot$total_cost,
         col = scales::alpha(plot_cols, 0.6), # Use alpha transparency
         pch = plot_pch,
         cex = 1.0, # Adjust point size
         xlab = "Total Discounted QALYs",
         ylab = "Total Discounted Cost (£)",
         main = NULL,
         xlim = c(xlim_range[1] - xlim_buffer, xlim_range[2] + xlim_buffer),
         ylim = c(ylim_range[1] - ylim_buffer, ylim_range[2] + ylim_buffer),
         las = 1,
         cex.lab = 1.1,
         cex.axis = 0.9,
         bty = "l" # Box type 'L'
         )
    grid(col="lightgray", lty="dotted")
    legend("topleft",
           legend = c("Dupimab", "Supimab"),
           col = c("#a6cee3", "#fdbf6f"),
           pch = c(16, 17),
           pt.cex = 1.2,
           bty = "n", # No box around legend
           cex = 1.0)
     box(col = "grey60")
  })

  # Histogram for Costs (Dupimab)
  output$plot_hist_cost_dup <- renderPlot({
      req(sim_results())
      costs <- sim_results()$patient_data$total_cost[sim_results()$patient_data$treatment == "Dupimab"]
      costs_finite <- costs[is.finite(costs)] # Remove NA/Inf
      if (length(costs_finite) == 0) return(NULL)

      par(mar = c(4, 5, 2, 1))
      hist(costs_finite,
           main = NULL,
           xlab = "Total Discounted Cost (£)",
           ylab = "Frequency",
           col = "#a6cee3",
           border = "grey20",
           las = 1,
           cex.lab = 1.0,
           cex.axis = 0.9,
           breaks = 30)
      rug(costs_finite, col=scales::alpha("black", 0.3)) # Add rug plot
      box(col = "grey60")
  })

  # Histogram for Costs (Supimab)
  output$plot_hist_cost_sup <- renderPlot({
      req(sim_results())
      costs <- sim_results()$patient_data$total_cost[sim_results()$patient_data$treatment == "Supimab"]
      costs_finite <- costs[is.finite(costs)]
      if (length(costs_finite) == 0) return(NULL)

      par(mar = c(4, 5, 2, 1))
      hist(costs_finite,
           main = NULL,
           xlab = "Total Discounted Cost (£)",
           ylab = "Frequency",
           col = "#fdbf6f",
           border = "grey20",
           las = 1,
           cex.lab = 1.0,
           cex.axis = 0.9,
           breaks = 30)
      rug(costs_finite, col=scales::alpha("black", 0.3))
      box(col = "grey60")
  })

   # Histogram for QALYs (Dupimab)
  output$plot_hist_qaly_dup <- renderPlot({
      req(sim_results())
      qalys <- sim_results()$patient_data$total_qaly[sim_results()$patient_data$treatment == "Dupimab"]
      qalys_finite <- qalys[is.finite(qalys)]
      if (length(qalys_finite) == 0) return(NULL)

      par(mar = c(4, 5, 2, 1))
      hist(qalys_finite,
           main = NULL,
           xlab = "Total Discounted QALYs",
           ylab = "Frequency",
           col = "#a6cee3",
           border = "grey20",
           las = 1,
           cex.lab = 1.0,
           cex.axis = 0.9,
           breaks = 30)
      rug(qalys_finite, col=scales::alpha("black", 0.3))
      box(col = "grey60")
  })

  # Histogram for QALYs (Supimab)
  output$plot_hist_qaly_sup <- renderPlot({
      req(sim_results())
      qalys <- sim_results()$patient_data$total_qaly[sim_results()$patient_data$treatment == "Supimab"]
      qalys_finite <- qalys[is.finite(qalys)]
      if (length(qalys_finite) == 0) return(NULL)

      par(mar = c(4, 5, 2, 1))
      hist(qalys_finite,
           main = NULL,
           xlab = "Total Discounted QALYs",
           ylab = "Frequency",
           col = "#fdbf6f",
           border = "grey20",
           las = 1,
           cex.lab = 1.0,
           cex.axis = 0.9,
           breaks = 30)
      rug(qalys_finite, col=scales::alpha("black", 0.3))
      box(col = "grey60")
  })

  # Render DataTable
   output$patient_table <- DT::renderDataTable({
       req(sim_results())
       df <- sim_results()$patient_data
       # Round numeric columns for display
       numeric_cols_dt <- c("total_cost", "total_qaly", "time_of_stroke", "time_of_death")
       df[numeric_cols_dt] <- lapply(df[numeric_cols_dt], function(x) if(is.numeric(x)) round(x, 2) else x)

       DT::datatable(df,
                     options = list(pageLength = 10, scrollX = TRUE), # Sensible defaults
                     rownames = FALSE,
                     filter = 'top', # Add column filters
                     class = 'cell-border stripe hover') # Apply some styling
   })

}

# --- Run the Shiny App ---
# This command starts the Shiny application.
shinyApp(ui = ui, server = server)

