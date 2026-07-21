library(shiny)
library(readr)
library(dplyr)
library(bslib)
library(plotly)
library(wesanderson)

required_files <- c(
  "output_ppv_als_grid.csv",
  "output_ppv_alsftd_grid.csv",
  "www/APECS.png"
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
  format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
}

format_prob <- function(x) {
  paste0(
    format(
      round(x * 100, 1),
      nsmall = 1,
      big.mark = ".",
      decimal.mark = ",",
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
        big.mark = ".",
        decimal.mark = ",",
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
    style = "display: flex; align-items: center; justify-content: flex-start; width: 100%; gap: 20px;",
    tags$img(
      src = "APECS.png",
      height = "100px",
      style = "object-fit: contain;"
    ),
    div(
      div(
        style = "font-weight: 600; font-size: 22px; line-height: 1.2;",
        HTML(
          "ALS Family History - Monogenic Probability Calculator<br>
          <span style='font-weight: 400; font-size: 18px;'>Based on Mendelian and complex inheritance theory</span>"
        )
      ),
      div(
        style = "font-size: 16px; color: #555; margin-top: 4px;",
        tags$em("APECS - ALS PEdigree simulations under a Complex and Simple disease model")
      )
    )
  ),
  theme = bs_theme(version = 5),
  sidebar = sidebar(
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
  card(
    card_header("Estimated Probability of Monogenic Disease (PPV)"),
    verbatimTextOutput("ppv_text")
  ),
  card(
    card_header("Prior Probability vs Family History Probability"),
    plotlyOutput("prob_bar", height = "260px")
  ),
  card(
    card_header("Simulated Pedigrees Matching Input Family History"),
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

  output$ppv_text <- renderText({
    row <- selected_row()

    paste0(
      "Probability of monogenic disease: ", format_prob(row$PPV),
      " (95% CI ", format_prob(row$PPV_CI_low), " – ", format_prob(row$PPV_CI_high), ")\n",
      "Based on ", format_count(row$n), " matching simulated pedigrees."
    )
  })

  output$prob_bar <- renderPlotly({
    row <- selected_row()

    prior_monogenic <- prior_mendelian_n / prior_total_n
    prior_polygenic <- 1 - prior_monogenic
    post_monogenic <- row$n_mendelian / row$n
    post_polygenic <- row$n_non_mendelian / row$n

    plot_df <- data.frame(
      scenario = c("Prior probability", "Posterior probability"),
      monogenic = c(prior_monogenic, post_monogenic),
      polygenic = c(prior_polygenic, post_polygenic),
      monogenic_n = c(prior_mendelian_n, row$n_mendelian),
      polygenic_n = c(prior_total_n - prior_mendelian_n, row$n_non_mendelian),
      total_n = c(prior_total_n, row$n)
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
        marker = list(color = darjeeling_cols[1]),
        text = ~paste0(format_prob(monogenic), "<br>(n = ", format_count(monogenic_n), ")"),
        textposition = "inside",
        insidetextanchor = "middle",
        hovertemplate = paste(
          "<b>%{y}</b><br>",
          "Monogenic: %{text}<br>",
          "Total: %{customdata}<extra></extra>"
        ),
        customdata = ~format_count(total_n)
      ) %>%
      add_trace(
        x = ~polygenic,
        y = ~scenario,
        name = "Polygenic",
        type = "bar",
        orientation = "h",
        marker = list(color = darjeeling_cols[2]),
        text = ~paste0(format_prob(polygenic), "<br>(n = ", format_count(polygenic_n), ")"),
        textposition = "inside",
        insidetextanchor = "middle",
        hovertemplate = paste(
          "<b>%{y}</b><br>",
          "Polygenic: %{text}<br>",
          "Total: %{customdata}<extra></extra>"
        ),
        customdata = ~format_count(total_n)
      ) %>%
      layout(
        barmode = "stack",
        xaxis = list(
          title = "Proportion",
          tickformat = ".0%",
          range = c(0, 1)
        ),
        yaxis = list(
          title = "",
          automargin = TRUE
        ),
        legend = list(
          orientation = "h",
          xanchor = "center",
          x = 0.5,
          yanchor = "bottom",
          y = 1.10
        ),
        margin = list(l = 180, r = 30, t = 20, b = 50)
      )
  })

  output$match_tbl <- renderUI({
    row <- selected_row()

    if (input$mode == "ALS only") {
      tags$table(
        class = "table table-striped table-bordered table-sm",
        tags$thead(
          tags$tr(
            tags$th("ALS family history"),
            tags$th("Number of affected index patients"),
            tags$th("Number of monogenic index patients"),
            tags$th("Number of polygenic index patients"),
            tags$th("Prevalence of this specific family history")
          )
        ),
        tags$tbody(
          tags$tr(
            tags$td(HTML(format_als_history(row))),
            tags$td(format_count(row$n)),
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
            tags$th("Number of affected index patients"),
            tags$th("Number of monogenic index patients"),
            tags$th("Number of polygenic index patients"),
            tags$th("Prevalence of this specific family history")
          )
        ),
        tags$tbody(
          tags$tr(
            tags$td(HTML(format_als_history(row))),
            tags$td(HTML(format_ftd_history(row))),
            tags$td(format_count(row$n)),
            tags$td(format_count(row$n_mendelian)),
            tags$td(format_count(row$n_non_mendelian)),
            tags$td(format_prevalence(row$prevalence))
          )
        )
      )
    }
  })
}

shinyApp(ui = ui, server = server)