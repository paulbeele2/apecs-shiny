library(shiny)
library(readr)
library(dplyr)
library(bslib)
library(plotly)
library(wesanderson)

required_files <- c(
  "output_ppv_als_grid.csv",
  "output_ppv_alsftd_grid.csv",
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

prior_mendelian_n <- 44321
prior_total_n <- 272651

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

ui <- page_sidebar(
  title = div(
    class = "title-wrap",
    tags$img(
      src = "APECS_logo.png",
      class = "apecs-logo"
    ),
    div(
      class = "title-text",
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
  theme = bs_theme(version = 5),

  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      .title-wrap {
        display: flex;
        align-items: center;
        justify-content: flex-start;
        width: 100%;
        gap: 20px;
      }

      .title-text {
        min-width: 0;
      }

      .apecs-logo {
        height: 100px;
        width: auto;
        object-fit: contain;
        flex-shrink: 0;
      }

      .app-title {
        font-size: 22px;
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

      .svg-wrap {
        padding: 10px;
        min-height: 260px;
        max-height: 360px;
        overflow: hidden;
        display: flex;
        align-items: center;
        justify-content: center;
      }

      .svg-card-img {
        width: 100%;
        height: 100%;
        max-height: 340px;
        object-fit: contain;
        display: block;
      }

      .ppv-box {
        padding: 15px;
      }

      .table-wrap {
        overflow-x: auto;
        -webkit-overflow-scrolling: touch;
      }

      .table {
        margin-bottom: 0;
      }

      .table th,
      .table td {
        vertical-align: middle;
      }

      .table th:nth-child(n+2),
      .table td:nth-child(n+2) {
        text-align: center;
        white-space: nowrap;
      }

      @media (max-width: 768px) {
        .title-wrap {
          flex-direction: column;
          align-items: flex-start;
          gap: 12px;
        }

        .apecs-logo {
          height: 64px;
        }

        .app-title {
          font-size: 18px;
        }

        .app-subtitle {
          font-size: 14px;
        }

        .app-fullname {
          font-size: 12px;
        }

        .card-header {
          font-size: 14px !important;
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
          font-size: 12px !important;
          line-height: 1.4;
        }

        .svg-wrap {
          min-height: 180px;
          max-height: 240px;
          padding: 8px;
        }

        .svg-card-img {
          max-height: 220px;
        }

        .ppv-box {
          padding: 12px;
        }
      }
    "))
  ),

  sidebar = sidebar(
    width = 260,
    open = "desktop",
    selectInput("mode", "Model", choices = c("ALS only", "ALS + FTD")),

    selectInput(
      "als1",
      "1st-degree relatives with ALS",
      choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
      selected = "0"
    ),
    selectInput(
      "als2",
      "2nd-degree relatives with ALS",
      choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
      selected = "0"
    ),
    selectInput(
      "als3",
      "3rd-degree relatives with ALS",
      choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
      selected = "0"
    ),

    conditionalPanel(
      condition = "input.mode == 'ALS + FTD'",
      selectInput(
        "ftd1",
        "1st-degree relatives with FTD",
        choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
        selected = "0"
      ),
      selectInput(
        "ftd2",
        "2nd-degree relatives with FTD",
        choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
        selected = "0"
      ),
      selectInput(
        "ftd3",
        "3rd-degree relatives with FTD",
        choices = c("0", "1", "2", "3", "4", "5", "Unknown"),
        selected = "0"
      )
    )
  ),

  layout_columns(
    col_widths = breakpoints(
      sm = c(12, 12),
      lg = c(8, 4)
    ),

    card(
      full_screen = FALSE,
      card_header("How to count relatives"),
      div(
        style = "padding: 10px; height: 280px; overflow: hidden;",
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
        style = "padding: 15px;",
        p(
          "Index patient (individual A) is marked by the black arrow. ",
          "The degree of relatives to individual A is illustrated by the number in each individual."
        ),
        p(
          "Note that for ‘any dementia’-affected relative, both FTD- (individual D) ",
          "and other dementia-affected (individual C) relatives are considered."
        ),
        p(
          "Comorbid ALS-FTD (individual B) is only counted as ALS once."
        )
      )
    )
  ),

  layout_columns(
    col_widths = breakpoints(
      sm = c(12, 12),
      lg = c(8, 4)
    ),

    card(
      full_screen = FALSE,
      card_header("Prior Probability vs. Family History Probability for Monogenic Disease"),
      plotlyOutput("prob_bar", width = "100%", height = "260px")
    ),

    card(
      full_screen = FALSE,
      card_header("Estimated Probability of Monogenic Disease"),
      uiOutput("ppv_box")
    )
  ),

  card(
    card_header("Simulated Pedigrees matching Family History"),
    uiOutput("match_tbl")
  )
)

server <- function(input, output, session) {

  selected_row <- reactive({
    df <- if (input$mode == "ALS only") als_grid else alsftd_grid

    if (input$als1 != "Unknown") {
      df <- df %>% filter(relatives_1st_als == as.numeric(input$als1))
    }
    if (input$als2 != "Unknown") {
      df <- df %>% filter(relatives_2nd_als == as.numeric(input$als2))
    }
    if (input$als3 != "Unknown") {
      df <- df %>% filter(relatives_3rd_als == as.numeric(input$als3))
    }

    if (input$mode == "ALS + FTD") {
      if (input$ftd1 != "Unknown") {
        df <- df %>% filter(relatives_1st_ftd_unique == as.numeric(input$ftd1))
      }
      if (input$ftd2 != "Unknown") {
        df <- df %>% filter(relatives_2nd_ftd_unique == as.numeric(input$ftd2))
      }
      if (input$ftd3 != "Unknown") {
        df <- df %>% filter(relatives_3rd_ftd_unique == as.numeric(input$ftd3))
      }
    }

    validate(
      need(nrow(df) > 0, "No matching pedigree pattern found in the lookup table.")
    )

    agg <- df %>%
      summarise(
        relatives_1st_als = if (input$als1 == "Unknown") NA_real_ else first(relatives_1st_als),
        relatives_2nd_als = if (input$als2 == "Unknown") NA_real_ else first(relatives_2nd_als),
        relatives_3rd_als = if (input$als3 == "Unknown") NA_real_ else first(relatives_3rd_als),
        relatives_1st_ftd_unique = if ("relatives_1st_ftd_unique" %in% names(df)) {
          if (input$ftd1 == "Unknown") NA_real_ else first(relatives_1st_ftd_unique)
        } else {
          NA_real_
        },
        relatives_2nd_ftd_unique = if ("relatives_2nd_ftd_unique" %in% names(df)) {
          if (input$ftd2 == "Unknown") NA_real_ else first(relatives_2nd_ftd_unique)
        } else {
          NA_real_
        },
        relatives_3rd_ftd_unique = if ("relatives_3rd_ftd_unique" %in% names(df)) {
          if (input$ftd3 == "Unknown") NA_real_ else first(relatives_3rd_ftd_unique)
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

  output$ppv_box <- renderUI({
    row <- selected_row()

    div(
      class = "ppv-box",
      p(
        HTML(paste0(
          "Probability of monogenic disease: <strong>",
          format_prob(row$PPV),
          "</strong> (95% CI <strong>",
          format_prob(row$PPV_CI_low),
          "</strong> – <strong>",
          format_prob(row$PPV_CI_high),
          "</strong>)<br><br>",
          "Based on ",
          format_count(row$n),
          " matching simulated pedigrees."
        ))
      )
    )
  })

  output$prob_bar <- renderPlotly({
    row <- selected_row()

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

    darjeeling_cols <- grDevices::adjustcolor(
      wes_palette("Darjeeling1", 5, type = "discrete"),
      alpha.f = 0.9
    )

    plot_ly(plot_df) %>%
      add_trace(
        x = ~monogenic,
        y = ~scenario,
        name = "Monogenic",
        type = "bar",
        orientation = "h",
        marker = list(
          color = darjeeling_cols[1],
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
          color = darjeeling_cols[2],
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
        hoverlabel = list(font = list(size = 12)),
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
          font = list(size = 12)
        ),
        margin = list(l = 85, r = 8, t = 6, b = 24),
        dragmode = FALSE
      ) %>%
      config(
        displayModeBar = FALSE,
        staticPlot = FALSE,
        scrollZoom = FALSE,
        doubleClick = FALSE,
        showTips = FALSE
      )
  })

  output$match_tbl <- renderUI({
    row <- selected_row()

    if (input$mode == "ALS only") {
      div(
        class = "table-wrap",
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
      )
    } else {
      div(
        class = "table-wrap",
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
      )
    }
  })
}

shinyApp(ui = ui, server = server)