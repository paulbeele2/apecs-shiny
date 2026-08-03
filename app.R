library(shiny)
library(shinyWidgets)
library(readr)
library(dplyr)
library(bslib)
library(plotly)

required_files <- c(
  "output_ppv_als_grid.csv",
  "output_ppv_alsftd_grid.csv",
  "varying_param_als_grid.csv",
  "varying_param_alsftd_grid.csv",
  "www/APECS_logo.png",
  "www/APECS_relative_count.svg"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required files: ",
    paste(missing_files, collapse = ", ")
  )
}

als_grid <- read_csv("output_ppv_als_grid.csv", show_col_types = FALSE)
alsftd_grid <- read_csv("output_ppv_alsftd_grid.csv", show_col_types = FALSE)
varying_als_grid <- read_csv("varying_param_als_grid.csv", show_col_types = FALSE)
varying_alsftd_grid <- read_csv("varying_param_alsftd_grid.csv", show_col_types = FALSE)

common_vals <- sort(unique(varying_als_grid$ALSFTD_common))
rare_vals   <- sort(unique(varying_als_grid$ALSFTD_rare))
h2_vals     <- sort(unique(varying_als_grid$h2als))

format_count <- function(x) {
  format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE, trim = TRUE)
}

format_prob <- function(x) {
  paste0(
    format(
      round(x * 100, 1),
      nsmall = 1,
      big.mark = ",",
      decimal.mark = ".",
      scientific = FALSE,
      trim = TRUE
    ),
    "%"
  )
}

format_param <- function(x) {
  format(
    x,
    nsmall = 2,
    big.mark = ",",
    decimal.mark = ".",
    scientific = FALSE,
    trim = TRUE
  )
}

format_prevalence <- function(x) {
  matches <- gregexpr("\\d+", x)
  extracted <- regmatches(x, matches)

  formatted <- lapply(extracted, function(vals) {
    sapply(vals, function(v) {
      format(
        as.numeric(v),
        big.mark = ",",
        decimal.mark = ".",
        scientific = FALSE,
        trim = TRUE
      )
    })
  })

  regmatches(x, matches) <- formatted
  x
}

format_als_history <- function(row) {
  paste0(
    "1st degree ALS = ", ifelse(is.na(row$relatives_1st_als), "Unknown", row$relatives_1st_als),
    "<br>2nd degree ALS = ", ifelse(is.na(row$relatives_2nd_als), "Unknown", row$relatives_2nd_als),
    "<br>3rd degree ALS = ", ifelse(is.na(row$relatives_3rd_als), "Unknown", row$relatives_3rd_als)
  )
}

format_ftd_history <- function(row) {
  paste0(
    "1st degree FTD = ", ifelse(is.na(row$relatives_1st_ftd_unique), "Unknown", row$relatives_1st_ftd_unique),
    "<br>2nd degree FTD = ", ifelse(is.na(row$relatives_2nd_ftd_unique), "Unknown", row$relatives_2nd_ftd_unique),
    "<br>3rd degree FTD = ", ifelse(is.na(row$relatives_3rd_ftd_unique), "Unknown", row$relatives_3rd_ftd_unique)
  )
}

wilson_ci <- function(x, n, conf_level = 0.95) {
  if (n == 0 || is.na(n) || is.na(x) || x < 0 || x > n) {
    return(c(estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }

  z <- qnorm(1 - (1 - conf_level) / 2)
  p_hat <- x / n
  z2 <- z^2

  denom <- 1 + z2 / n
  center <- (p_hat + z2 / (2 * n)) / denom
  margin <- (z * sqrt((p_hat * (1 - p_hat) / n) + (z2 / (4 * n^2)))) / denom

  c(
    estimate = p_hat,
    ci_lower = max(0, center - margin),
    ci_upper = min(1, center + margin)
  )
}

build_plot <- function(row, prior_mendelian_n, prior_total_n) {
  prior_monogenic <- prior_mendelian_n / prior_total_n
  prior_polygenic <- 1 - prior_monogenic
  post_monogenic <- row$n_mendelian / row$n
  post_polygenic <- row$n_non_mendelian / row$n

  plot_df <- data.frame(
    scenario = factor(c("Prior", "Family history"), levels = c("Family history", "Prior")),
    monogenic = c(prior_monogenic, post_monogenic),
    polygenic = c(prior_polygenic, post_polygenic),
    monogenic_n = c(prior_mendelian_n, row$n_mendelian),
    polygenic_n = c(prior_total_n - prior_mendelian_n, row$n_non_mendelian),
    total_n = c(prior_total_n, row$n),
    monogenic_label = c(
      format_prob(prior_monogenic),
      format_prob(post_monogenic)
    ),
    polygenic_label = c(
      format_prob(prior_polygenic),
      format_prob(post_polygenic)
    ),
    ci_low_label = c("—", format_prob(row$PPV_CI_low)),
    ci_high_label = c("—", format_prob(row$PPV_CI_high))
  )

  plot_df$monogenic_hover <- ifelse(
    plot_df$scenario == "Family history",
    paste0(
      "<b>", plot_df$scenario, "</b><br>",
      "Monogenic: ", plot_df$monogenic_label, "<br>",
      "95% CI: ", plot_df$ci_low_label, " – ", plot_df$ci_high_label, "<br>",
      "n = ", format_count(plot_df$monogenic_n), " of ", format_count(plot_df$total_n),
      "<extra></extra>"
    ),
    paste0(
      "<b>", plot_df$scenario, "</b><br>",
      "Monogenic: ", plot_df$monogenic_label, "<br>",
      "n = ", format_count(plot_df$monogenic_n), " of ", format_count(plot_df$total_n),
      "<extra></extra>"
    )
  )

  plot_df$polygenic_hover <- ifelse(
    plot_df$scenario == "Family history",
    paste0(
      "<b>", plot_df$scenario, "</b><br>",
      "Polygenic: ", plot_df$polygenic_label, "<br>",
      "95% CI: ", plot_df$ci_low_label, " – ", plot_df$ci_high_label, "<br>",
      "n = ", format_count(plot_df$polygenic_n), " of ", format_count(plot_df$total_n),
      "<extra></extra>"
    ),
    paste0(
      "<b>", plot_df$scenario, "</b><br>",
      "Polygenic: ", plot_df$polygenic_label, "<br>",
      "n = ", format_count(plot_df$polygenic_n), " of ", format_count(plot_df$total_n),
      "<extra></extra>"
    )
  )

  plot_ly(plot_df) %>%
    add_trace(
      x = ~monogenic,
      y = ~scenario,
      name = "Monogenic",
      type = "bar",
      orientation = "h",
      marker = list(
        color = "#ed2024",
        line = list(color = "white", width = 1)
      ),
      text = ~monogenic_label,
      textposition = "inside",
      insidetextanchor = "middle",
      textfont = list(size = 11, color = "white"),
      hovertemplate = ~monogenic_hover
    ) %>%
    add_trace(
      x = ~polygenic,
      y = ~scenario,
      name = "Polygenic",
      type = "bar",
      orientation = "h",
      marker = list(
        color = "#5cbcd6",
        line = list(color = "white", width = 1)
      ),
      text = ~polygenic_label,
      textposition = "inside",
      insidetextanchor = "middle",
      textfont = list(size = 11, color = "white"),
      hovertemplate = ~polygenic_hover
    ) %>%
    layout(
      barmode = "stack",
      bargap = 0.25,
      font = list(size = 13),
      uniformtext = list(minsize = 9, mode = "show"),
      hoverlabel = list(
        font = list(
          size = 12,
          color = "white"
        )
      ),
      xaxis = list(
        title = "",
        tickformat = ".0%",
        range = c(0, 1),
        showgrid = TRUE,
        gridcolor = "#e9ecef",
        zeroline = FALSE,
        fixedrange = TRUE
      ),
      yaxis = list(
        title = "",
        automargin = TRUE,
        fixedrange = TRUE
      ),
      legend = list(
        orientation = "h",
        xanchor = "center",
        x = 0.5,
        yanchor = "bottom",
        y = 1.01,
        font = list(size = 12),
        traceorder = "normal"
      ),
      margin = list(l = 8, r = 8, t = 8, b = 8),
      dragmode = FALSE
    ) %>%
    config(
      displayModeBar = FALSE,
      staticPlot = FALSE,
      scrollZoom = FALSE,
      doubleClick = FALSE,
      showTips = FALSE
    )
}

make_match_table <- function(row, mode) {
  if (mode == "ALS only") {
    tags$table(
      class = "table table-striped table-bordered table-sm",
      tags$thead(
        tags$tr(
          tags$th("ALS family history"),
          tags$th("Number of monogenic index patients"),
          tags$th("Number of polygenic index patients"),
          tags$th("Prevalence of this specific family history")
        )
      ),
      tags$tbody(
        tags$tr(
          tags$td(HTML(format_als_history(row))),
          tags$td(format_count(row$n_mendelian)),
          tags$td(format_count(row$n_non_mendelian)),
          tags$td(format_prevalence(row$prevalence))
        )
      )
    )
  } else {
    tags$table(
      class = "table table-striped table-bordered table-sm",
      tags$thead(
        tags$tr(
          tags$th("ALS family history"),
          tags$th("FTD family history"),
          tags$th("Number of monogenic index patients"),
          tags$th("Number of polygenic index patients"),
          tags$th("Family history prevalence")
        )
      ),
      tags$tbody(
        tags$tr(
          tags$td(HTML(format_als_history(row))),
          tags$td(HTML(format_ftd_history(row))),
          tags$td(format_count(row$n_mendelian)),
          tags$td(format_count(row$n_non_mendelian)),
          tags$td(format_prevalence(row$prevalence))
        )
      )
    )
  }
}

ui <- page_fluid(
  theme = bs_theme(version = 5),

  tags$head(
    tags$style(HTML("
      .app-shell {
        max-width: 1600px;
        margin: 0 auto;
        padding-top: 10px;
        padding-bottom: 20px;
      }

      .app-header {
        display: flex;
        align-items: center;
        justify-content: flex-start;
        gap: 16px;
        flex-wrap: wrap;
        margin-bottom: 18px;
        padding: 8px 4px 2px 4px;
      }

      .app-header-logo {
        height: 82px;
        width: auto;
        object-fit: contain;
        flex-shrink: 0;
      }

      .app-header-text {
        min-width: 280px;
      }

      .app-title {
        font-size: 24px;
        font-weight: 600;
        line-height: 1.2;
      }

      .app-subtitle {
        font-size: 18px;
        font-weight: 400;
        line-height: 1.2;
      }

      .app-fullname {
        font-size: 15px;
        color: #555;
        margin-top: 4px;
      }

      .card-header {
        font-size: 15px !important;
        font-weight: 600;
      }

      .card-body,
      .card-body p,
      .card-body div,
      .card-body span,
      .sidebar,
      .sidebar .form-label,
      .sidebar .selectize-input,
      .sidebar .selectize-dropdown,
      .sidebar .form-select,
      .sidebar .control-label,
      .sidebar label,
      .table,
      .table th,
      .table td,
      .shiny-input-container,
      .shiny-input-container label,
      .form-control,
      .form-select {
        font-size: 13px !important;
        line-height: 1.5;
      }

      .sidebar .form-label,
      .sidebar label,
      .sidebar .control-label,
      .table th {
        font-weight: 600;
      }

      .shiny-html-output p:last-child,
      .card-body p:last-child {
        margin-bottom: 0;
      }

      .nav-tabs .nav-link {
        font-size: 14px;
        font-weight: 600;
      }

      .irs-grid-text {
        font-size: 11px !important;
      }

      @media (max-width: 767px) {
        .app-header {
          align-items: flex-start;
        }

        .app-header-logo {
          height: 64px;
        }

        .app-title {
          font-size: 20px;
        }

        .app-subtitle {
          font-size: 15px;
        }

        .app-fullname {
          font-size: 13px;
        }
      }
    "))
  ),

  div(
    class = "app-shell",

    div(
      class = "app-header",
      tags$img(
        src = "APECS_logo.png",
        class = "app-header-logo"
      ),
      div(
        class = "app-header-text",
        div(
          class = "app-title",
          HTML(
            "ALS Family History - Monogenic Probability Calculator<br>
            <span class='app-subtitle'>Based on Mendelian and complex inheritance theory</span>"
          )
        ),
        div(
          class = "app-fullname",
          tags$em("APECS - ALS PEdigree simulations under a Complex and Simple disease model")
        )
      )
    ),

    navset_card_tab(
      id = "analysis_tab",

      nav_panel(
        "Main parameters simulation",
        layout_sidebar(
          sidebar = sidebar(
            width = 260,
            position = "left",
            open = list(
              desktop = "open",
              mobile = "closed"
            ),

            selectInput("mode_main", "Model", choices = c("ALS only", "ALS + FTD")),

            selectInput(
              "als1_main",
              "1st-degree relatives with ALS",
              choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
              selected = "0"
            ),
            selectInput(
              "als2_main",
              "2nd-degree relatives with ALS",
              choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
              selected = "0"
            ),
            selectInput(
              "als3_main",
              "3rd-degree relatives with ALS",
              choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
              selected = "0"
            ),

            conditionalPanel(
              condition = "input.mode_main == 'ALS + FTD'",
              selectInput(
                "ftd1_main",
                "1st-degree relatives with FTD",
                choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
                selected = "0"
              ),
              selectInput(
                "ftd2_main",
                "2nd-degree relatives with FTD",
                choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
                selected = "0"
              ),
              selectInput(
                "ftd3_main",
                "3rd-degree relatives with FTD",
                choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
                selected = "0"
              )
            )
          ),

          layout_columns(
            col_widths = c(8, 4),

            card(
              full_screen = FALSE,
              card_header("How to count relatives"),
              div(
                style = "padding: 4px; height: 100%; overflow: hidden;",
                tags$img(
                  src = "APECS_relative_count.svg",
                  style = "width: 100%; height: 100%; object-fit: contain; display: block;"
                )
              )
            ),

            card(
              full_screen = FALSE,
              card_header("Counting affected relatives"),
              div(
                style = "padding: 5px;",
                p(
                  "Index patient (individual A) is marked by the black arrow. ",
                  "The degree of relatives to individual A is illustrated by the number in each individual."
                ),
                p(
                  "Comorbid ALS-FTD (individual B) is only counted as ALS once."
                ),
                p(
                  "Note that for ‘any dementia’-affected relative, both FTD- (individual D) ",
                  "and other dementia-affected (individual C) relatives are considered."
                ),
                p(
                  "APECS was benchmarked in 3-generation pedigrees, limiting ",
                  "predictions involving more distant affected relatives."
                )
              )
            )
          ),

          layout_columns(
            col_widths = c(8, 4),

            card(
              full_screen = FALSE,
              card_header("Prior Probability vs. Family History Probability for Monogenic Disease"),
              plotlyOutput("prob_bar_main", height = "400px")
            ),

            card(
              full_screen = FALSE,
              card_header("Estimated Probability of Monogenic Disease"),
              uiOutput("ppv_box_main")
            )
          ),

          card(
            card_header("Simulated Pedigrees matching Family History"),
            uiOutput("match_tbl_main")
          )
        )
      ),

      nav_panel(
        "Variable parameters",
        layout_sidebar(
          sidebar = sidebar(
            width = 260,
            position = "left",
            open = list(
              desktop = "open",
              mobile = "closed"
            ),

            selectInput("mode_var", "Model", choices = c("ALS only", "ALS + FTD")),

            selectInput(
              "als1_var",
              "1st-degree relatives with ALS",
              choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
              selected = "0"
            ),
            selectInput(
              "als2_var",
              "2nd-degree relatives with ALS",
              choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
              selected = "0"
            ),
            selectInput(
              "als3_var",
              "3rd-degree relatives with ALS",
              choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
              selected = "0"
            ),

            conditionalPanel(
              condition = "input.mode_var == 'ALS + FTD'",
              selectInput(
                "ftd1_var",
                "1st-degree relatives with FTD",
                choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
                selected = "0"
              ),
              selectInput(
                "ftd2_var",
                "2nd-degree relatives with FTD",
                choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
                selected = "0"
              ),
              selectInput(
                "ftd3_var",
                "3rd-degree relatives with FTD",
                choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
                selected = "0"
              )
            ),

            tags$hr(),

            sliderTextInput(
              inputId = "common_var",
              label = "ALS/FTD common variant penetrance",
              choices = format_param(common_vals),
              selected = format_param(0.20),
              grid = TRUE,
              force_edges = TRUE
            ),
            sliderTextInput(
              inputId = "rare_var",
              label = "ALS/FTD rare variant penetrance",
              choices = format_param(rare_vals),
              selected = format_param(0.50),
              grid = TRUE,
              force_edges = TRUE
            ),
            sliderTextInput(
              inputId = "h2_var",
              label = "ALS heritability",
              choices = format_param(h2_vals),
              selected = format_param(0.40),
              grid = TRUE,
              force_edges = TRUE
            )
          ),

          layout_columns(
            col_widths = c(8, 4),

            card(
              full_screen = FALSE,
              card_header("How to count relatives"),
              div(
                style = "padding: 4px; height: 100%; overflow: hidden;",
                tags$img(
                  src = "APECS_relative_count.svg",
                  style = "width: 100%; height: 100%; object-fit: contain; display: block;"
                )
              )
            ),

            card(
              full_screen = FALSE,
              card_header("Selected simulation parameters"),
              uiOutput("param_summary_var")
            )
          ),

          layout_columns(
            col_widths = c(8, 4),

            card(
              full_screen = FALSE,
              card_header("Prior Probability vs. Family History Probability for Monogenic Disease"),
              plotlyOutput("prob_bar_var", height = "400px")
            ),

            card(
              full_screen = FALSE,
              card_header("Estimated Probability of Monogenic Disease"),
              uiOutput("ppv_box_var")
            )
          ),

          card(
            card_header("Simulated Pedigrees matching Family History"),
            uiOutput("match_tbl_var")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {

  prior_main <- reactive({
    df <- if (input$mode_main == "ALS only") als_grid else alsftd_grid

    df %>%
      summarise(
        prior_total_n = sum(n),
        prior_mendelian_n = sum(n_mendelian)
      )
  })

  prior_var <- reactive({
    df <- if (input$mode_var == "ALS only") varying_als_grid else varying_alsftd_grid

    df %>%
      filter(
        ALSFTD_common == as.numeric(input$common_var),
        ALSFTD_rare == as.numeric(input$rare_var),
        h2als == as.numeric(input$h2_var)
      ) %>%
      summarise(
        prior_total_n = sum(n),
        prior_mendelian_n = sum(n_mendelian)
      )
  })

  selected_row_main <- reactive({
    df <- if (input$mode_main == "ALS only") als_grid else alsftd_grid

    if (input$als1_main != "Unknown") {
      df <- df %>% filter(relatives_1st_als == as.numeric(input$als1_main))
    }
    if (input$als2_main != "Unknown") {
      df <- df %>% filter(relatives_2nd_als == as.numeric(input$als2_main))
    }
    if (input$als3_main != "Unknown") {
      df <- df %>% filter(relatives_3rd_als == as.numeric(input$als3_main))
    }

    if (input$mode_main == "ALS + FTD") {
      if (input$ftd1_main != "Unknown") {
        df <- df %>% filter(relatives_1st_ftd_unique == as.numeric(input$ftd1_main))
      }
      if (input$ftd2_main != "Unknown") {
        df <- df %>% filter(relatives_2nd_ftd_unique == as.numeric(input$ftd2_main))
      }
      if (input$ftd3_main != "Unknown") {
        df <- df %>% filter(relatives_3rd_ftd_unique == as.numeric(input$ftd3_main))
      }
    }

    validate(
      need(nrow(df) > 0, "No matching pedigree pattern found in the lookup table.")
    )

    prior_total_n <- prior_main()$prior_total_n[[1]]

    agg <- df %>%
      summarise(
        relatives_1st_als = if (input$als1_main == "Unknown") NA_real_ else first(relatives_1st_als),
        relatives_2nd_als = if (input$als2_main == "Unknown") NA_real_ else first(relatives_2nd_als),
        relatives_3rd_als = if (input$als3_main == "Unknown") NA_real_ else first(relatives_3rd_als),
        relatives_1st_ftd_unique = if ("relatives_1st_ftd_unique" %in% names(df)) {
          if (input$ftd1_main == "Unknown") NA_real_ else first(relatives_1st_ftd_unique)
        } else {
          NA_real_
        },
        relatives_2nd_ftd_unique = if ("relatives_2nd_ftd_unique" %in% names(df)) {
          if (input$ftd2_main == "Unknown") NA_real_ else first(relatives_2nd_ftd_unique)
        } else {
          NA_real_
        },
        relatives_3rd_ftd_unique = if ("relatives_3rd_ftd_unique" %in% names(df)) {
          if (input$ftd3_main == "Unknown") NA_real_ else first(relatives_3rd_ftd_unique)
        } else {
          NA_real_
        },
        n = sum(n),
        n_mendelian = sum(n_mendelian),
        n_non_mendelian = sum(n_non_mendelian)
      )

    ci <- wilson_ci(agg$n_mendelian, agg$n)

    agg %>%
      mutate(
        PPV = unname(ci["estimate"]),
        PPV_CI_low = unname(ci["ci_lower"]),
        PPV_CI_high = unname(ci["ci_upper"]),
        prevalence = paste0(n, " out of ", prior_total_n, " patients")
      )
  })

  selected_row_var <- reactive({
    df_all <- if (input$mode_var == "ALS only") varying_als_grid else varying_alsftd_grid

    df_param <- df_all %>%
      filter(
        ALSFTD_common == as.numeric(input$common_var),
        ALSFTD_rare == as.numeric(input$rare_var),
        h2als == as.numeric(input$h2_var)
      )

    validate(
      need(nrow(df_param) > 0, "No pedigrees found for this parameter combination.")
    )

    total_selected_param_n <- sum(df_param$n)

    df <- df_param

    if (input$als1_var != "Unknown") {
      df <- df %>% filter(relatives_1st_als == as.numeric(input$als1_var))
    }
    if (input$als2_var != "Unknown") {
      df <- df %>% filter(relatives_2nd_als == as.numeric(input$als2_var))
    }
    if (input$als3_var != "Unknown") {
      df <- df %>% filter(relatives_3rd_als == as.numeric(input$als3_var))
    }

    if (input$mode_var == "ALS + FTD") {
      if (input$ftd1_var != "Unknown") {
        df <- df %>% filter(relatives_1st_ftd_unique == as.numeric(input$ftd1_var))
      }
      if (input$ftd2_var != "Unknown") {
        df <- df %>% filter(relatives_2nd_ftd_unique == as.numeric(input$ftd2_var))
      }
      if (input$ftd3_var != "Unknown") {
        df <- df %>% filter(relatives_3rd_ftd_unique == as.numeric(input$ftd3_var))
      }
    }

    validate(
      need(nrow(df) > 0, "No matching pedigree pattern found for this family history and parameter combination.")
    )

    agg <- df %>%
      summarise(
        relatives_1st_als = if (input$als1_var == "Unknown") NA_real_ else first(relatives_1st_als),
        relatives_2nd_als = if (input$als2_var == "Unknown") NA_real_ else first(relatives_2nd_als),
        relatives_3rd_als = if (input$als3_var == "Unknown") NA_real_ else first(relatives_3rd_als),
        relatives_1st_ftd_unique = if ("relatives_1st_ftd_unique" %in% names(df)) {
          if (input$ftd1_var == "Unknown") NA_real_ else first(relatives_1st_ftd_unique)
        } else {
          NA_real_
        },
        relatives_2nd_ftd_unique = if ("relatives_2nd_ftd_unique" %in% names(df)) {
          if (input$ftd2_var == "Unknown") NA_real_ else first(relatives_2nd_ftd_unique)
        } else {
          NA_real_
        },
        relatives_3rd_ftd_unique = if ("relatives_3rd_ftd_unique" %in% names(df)) {
          if (input$ftd3_var == "Unknown") NA_real_ else first(relatives_3rd_ftd_unique)
        } else {
          NA_real_
        },
        n = sum(n),
        n_mendelian = sum(n_mendelian),
        n_non_mendelian = sum(n_non_mendelian)
      )

    ci <- wilson_ci(agg$n_mendelian, agg$n)

    agg %>%
      mutate(
        ALSFTD_common = as.numeric(input$common_var),
        ALSFTD_rare = as.numeric(input$rare_var),
        h2als = as.numeric(input$h2_var),
        PPV = unname(ci["estimate"]),
        PPV_CI_low = unname(ci["ci_lower"]),
        PPV_CI_high = unname(ci["ci_upper"]),
        prevalence = paste0(n, " out of ", total_selected_param_n, " patients")
      )
  })

  output$ppv_box_main <- renderUI({
    row <- selected_row_main()

    div(
      style = "padding: 5px;",
      p(
        HTML(paste0(
          "Probability of monogenic disease: <strong>",
          format_prob(row$PPV),
          "</strong> (95% CI <strong>",
          format_prob(row$PPV_CI_low),
          "</strong> – <strong>",
          format_prob(row$PPV_CI_high),
          "</strong>);<br><br>",
          "Based on ",
          format_count(row$n),
          " matching simulated pedigrees."
        ))
      )
    )
  })

  output$prob_bar_main <- renderPlotly({
    row <- selected_row_main()
    prior <- prior_main()
    build_plot(row, prior$prior_mendelian_n[[1]], prior$prior_total_n[[1]])
  })

  output$match_tbl_main <- renderUI({
    row <- selected_row_main()
    make_match_table(row, input$mode_main)
  })

  output$param_summary_var <- renderUI({
    div(
      style = "padding: 5px;",
      p(
        HTML(paste0(
          "<strong>ALS/FTD moderate penetrance allele:</strong> ",
          format_param(as.numeric(input$common_var)), "<br>",
          "<strong>ALS/FTD high penetrance allele:</strong> ",
          format_param(as.numeric(input$rare_var)), "<br>",
          "<strong>ALS heritability:</strong> ",
          format_param(as.numeric(input$h2_var))
        ))
      )
    )
  })

  output$ppv_box_var <- renderUI({
    row <- selected_row_var()

    div(
      style = "padding: 5px;",
      p(
        HTML(paste0(
          "Probability of monogenic disease: <strong>",
          format_prob(row$PPV),
          "</strong> (95% CI <strong>",
          format_prob(row$PPV_CI_low),
          "</strong> – <strong>",
          format_prob(row$PPV_CI_high),
          "</strong>);<br><br>",
          "Based on ",
          format_count(row$n),
          " matching simulated pedigrees."
        ))
      )
    )
  })

  output$prob_bar_var <- renderPlotly({
    row <- selected_row_var()
    prior <- prior_var()
    build_plot(row, prior$prior_mendelian_n[[1]], prior$prior_total_n[[1]])
  })

  output$match_tbl_var <- renderUI({
    row <- selected_row_var()
    make_match_table(row, input$mode_var)
  })
}

shinyApp(ui = ui, server = server)