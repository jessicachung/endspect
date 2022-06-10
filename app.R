# TODO: need to batch norm GEO samples

library(shiny)
# library(shinyjs)
library(dplyr)
library(DT)
library(ggplot2)

debug <- TRUE
load("data/data.RData")


# Currently only use RWH data
sample_df <- sample_df %>% filter(study == "rwh")

cycle_plot <- function(dat, title="", subtitle="", plot_type="ggplot") {
  if (plot_type == "ggplot") {
    size <- 1
  } else {
    size <- 0.5
  }
  if (is.null(dat)) return(NULL)
  ggplot(dat,
         aes(x=time, y=exprs)) +
    geom_point(size=size) +
    scale_x_continuous(breaks=seq(0, 100, by=10)) +
    labs(title=title, subtitle=subtitle,
         x="Cycle Time", y="Normalised Expression") +
    theme_bw()
}



ui <- fluidPage(
  title="Endspect",
  # includeCSS("style.css"),

  br(), br(),

  fluidRow(
    column(4, wellPanel(
      h4("Settings:"),
      # https://stackoverflow.com/questions/47305703/how-to-color-the-slider-in-shiny-to-the-right-of-the-value-instead-of-to-the-lef

      # Input method:
      selectInput("select_input_method",
                  label=h5("Input Method:"),
                  choices=list(
                    "Slider" = "input_slider",
                    # "Slider Range" = "input_slider_range",
                    "Manual input" = "input_manual"
                  ),
                  selected="input_slider"
      ), # End selectInput

      hr(),

      h4("Cycle Region:"),

      conditionalPanel("input.select_input_method == 'input_slider'",

                       sliderInput("slider_cycle_region_midpoint",
                                   label=h5("Midpoint:"),
                                   min=0,
                                   max=100,
                                   value=60,
                                   post=" ± w",
                                   width="100%"),

                       numericInput("numeric_window_size",
                                    label=h5("Half window size (w):"),
                                    value=8, min=1, max=30)

      ), # End conditionalPanel (input_slider)

      conditionalPanel("input.select_input_method == 'input_manual'",
                       br(),
                       p(tags$b("Range 1")),
                       column(6,
                              numericInput("input_a1",
                                           label=h5("Start:"),
                                           value=60,
                                           min=0, max=100, step=0.5)
                       ),
                       column(6,
                              numericInput("input_a2",
                                           label=h5("End:"),
                                           value=62.5,
                                           min=0, max=100, step=0.5)
                       ),
                       p(tags$b("Range 2")),
                       column(6,
                              numericInput("input_b1",
                                           label=h5("Start:"),
                                           value=65,
                                           min=0, max=100, step=0.5)
                       ),
                       column(6,
                              numericInput("input_b2",
                                           label=h5("End:"),
                                           value=67.5,
                                           min=0, max=100, step=0.5)
                       ),
      ), # End conditionalPanel (input_manual)

      actionButton("submit_button", "Set region"),

      br(), hr(), br(),

      selectInput("advanced_options",
                  label="Advanced Options",
                  choices=list("Hide advanced options" = "hide_options",
                               "Show advanced options" = "show_options"),
                  selected="hide_options"),

      conditionalPanel("input.advanced_options == 'show_options'",
                       checkboxInput("show_p_values",
                                     "Show p-values for rapidly changing genes t-tests",
                                      value=FALSE)
      ) # End conditionalPanel (show_options)



    ) # End wellPanel
    ), # End column

    column(width=7, offset=0,
           tabsetPanel(type = "tabs",
                       tabPanel("README",
                                br(),
                                wellPanel("#TODO: Write documentation")
                       ), # End tabPanel README
                       tabPanel("Rapid Change",
                                br(),
                                wellPanel(textOutput("rapid_text")),
                                plotOutput("gg_rapid_change"),
                                br(),
                                dataTableOutput("dt_rapid_change")
                       ), # End tabPanel Rapid Change
                       tabPanel("Highest Expression",
                                br(),
                                wellPanel(textOutput("high_text")),
                                plotOutput("gg_high"),
                                br(),
                                dataTableOutput("dt_high")
                       ), # End tabPanel Highest Expression
                       tabPanel("Inspect Genes",
                                br(),
                                wellPanel(textOutput("inspect")),
                                plotOutput("gg_inspect"),
                                br(),
                                dataTableOutput("dt_inspect")
                       ), # End tabPanel Inspect Genes
                       selected="Rapid Change"
           ) # End tabsetPanel
    ) # End column
  ) # End fluidRow
) # End fluidPage




server <- function(input, output, session) {
  rv <- reactiveValues(
    a = c(51, 59),
    b = c(61, 69),
    a_samples = c(),
    b_samples = c(),
    ab_samples = c(),
    not_ab_samples = c(),
    rapid_df = NULL,
    high_df = NULL,
    inspect_df = gene_info %>% arrange(desc(R2)),
    rapid_selected_gene = "ENSG00000184502",
    high_selected_gene = "ENSG00000145832",
    inspect_selected_gene = "ENSG00000197442"
  )

  # observeEvent({input$select_input_method;
  #   input$slider_cycle_region_midpoint_1; input$slider_cycle_region_midpoint_2; input$slider_cycle_region_midpoint_w;
  #   input$numeric_window_size;

  observeEvent({input$submit_button
    }, {

    if (input$select_input_method == "input_slider") {
      midpoint <- input$slider_cycle_region_midpoint

      rv$a[1] <- (midpoint - input$numeric_window_size) %% 100
      rv$a[2] <- (midpoint - 1) %% 100
      rv$b[1] <- (midpoint + 1) %% 100
      rv$b[2] <- (midpoint + input$numeric_window_size) %% 100
      print(rv$a)
      print(rv$b)
    }

    if (input$select_input_method == "input_manual") {
      print(input$numeric_window_size)
      rv$a[1] <- input$input_a1 %% 100
      rv$a[2] <- input$input_a2 %% 100
      rv$b[1] <- input$input_b1 %% 100
      rv$b[2] <- input$input_b2 %% 100

      # If window a end and window b start are the same, remove the middle value for sample comparison
      if ((rv$a[2] == rv$b[1]) & (rv$a[1] == round(rv$a[1]))) {
        rv$a[2] <- (rv$a - 1) %% 100
        rv$b[1] <- (rv$b + 1) %% 100
      }

      if (debug) print(rv$a)
      if (debug) print(rv$b)
    }

  })

  # Get region samples
  observeEvent({rv$a; rv$b}, {
    if (rv$a[2] >= rv$a[1]) {
      rv$a_samples <- sample_df %>% filter(time %>% between(rv$a[1], rv$a[2])) %>% pull(sample_id)
    } else {
      rv$a_samples <- sample_df %>% filter(time > rv$a[1] | time < rv$a[2]) %>% pull(sample_id)
    }

    if (rv$b[2] >= rv$b[1]) {
      rv$b_samples <- sample_df %>% filter(time %>% between(rv$b[1], rv$b[2])) %>% pull(sample_id)
    } else {
      rv$b_samples <- sample_df %>% filter(time > rv$b[1] | time < rv$b[2]) %>% pull(sample_id)
    }

    if (rv$b[2] > rv$a[1]) {
      rv$ab_samples <- sample_df %>% filter(time %>% between(rv$a[1], rv$b[2])) %>% pull(sample_id)
    } else {
      rv$ab_samples <- sample_df %>% filter(time > rv$a[1] | time < rv$b[2]) %>% pull(sample_id)
    }
    rv$not_ab_samples <- sample_df %>% filter(! sample_id %in% rv$ab_samples) %>% pull(sample_id)

    if (debug) print(rv$a_samples)
    if (debug) print(rv$b_samples)

    if (length(rv$a_samples) == 0 | length(rv$b_samples) == 0) {
      #TODO: error
    }

    # if (debug) print(rv$ab_samples)
    # Make sure there's no overlap between a and b samples
    # stopifnot(length(intersect(rv$a_samples, rv$b_samples)) == 0)
  })


  observeEvent({rv$a_samples; rv$b_samples; rv$ab_samples}, {

    # Simple difference of means
    # rv$gene_df <- data.frame(ensembl_id=rownames(exprs),
    #                       mean_a=apply(exprs[,rv$a_samples], 1, mean),
    #                       mean_b=apply(exprs[,rv$b_samples], 1, mean),
    #                       mean_ab=apply(exprs[,rv$ab_samples], 1, mean),
    #                       mean_not_ab=apply(exprs[,rv$not_ab_samples], 1, mean)
    #                       ) %>%
    #   mutate(sub_diff=round(mean_b-mean_a, 4),
    #          abs_diff=abs(sub_diff),
    #          full_diff=round(mean_ab-mean_not_ab, 4))


    # TODO: Check samples are non-zero

    t1 <- lapply(1:nrow(exprs), function(i) {
      t <- t.test(exprs[i,rv$a_samples], exprs[i,rv$b_samples])
      return(c(t$estimate, t$p.value))
    })
    # t2 <- lapply(1:nrow(exprs), function(i) {
    #   t <- t.test(exprs[i,rv$ab_samples], exprs[i,rv$not_ab_samples])
    #   return(c(t$estimate, t$p.value))
    # })

    rv$rapid_df <- data.frame(ensembl_id=rownames(exprs),
                              mean_a=sapply(t1, function(x) x[1]),
                              mean_b=sapply(t1, function(x) x[2])) %>%
      mutate(diff=round(mean_b-mean_a, 4),
             p.value=sapply(t1, function(x) x[3])) %>%
      merge(gene_info, by="ensembl_id") %>%
      arrange(p.value)

    rv$high_df <- data.frame(ensembl_id=rownames(exprs),
                             mean_ab=apply(exprs[,rv$ab_samples], 1, mean),
                             mean_not_ab=apply(exprs[,rv$not_ab_samples], 1, mean)) %>%
      mutate(diff=round(mean_ab-mean_not_ab, 4)) %>%
      merge(gene_info, by="ensembl_id") %>%
      arrange(desc(diff))

    rownames(rv$rapid_df) <- NULL
    rownames(rv$high_df) <- NULL

    rv$rapid_selected_gene <- rv$rapid_df %>% pull(ensembl_id) %>% head(1) %>% as.character()
    rv$high_selected_gene <- rv$high_df %>% pull(ensembl_id) %>% head(1) %>% as.character()

  })

  # Update selected gene when row is clicked
  observeEvent({input$dt_rapid_change_cell_clicked}, {
    if (!is.null(input$dt_rapid_change_cell_clicked$row)) {
      rv$rapid_selected_gene <- as.character(rv$rapid_df[input$dt_rapid_change_cell_clicked$row, "ensembl_id"])
      if (debug) message(sprintf("Rapid df row clicked: %d (%s)",
                                 input$dt_rapid_change_cell_clicked$row,
                                 rv$rapid_selected_gene))
    }
  })

  # Update selected gene when row is clicked
  observeEvent({input$dt_high_cell_clicked}, {
    if (!is.null(input$dt_high_cell_clicked$row)) {
      rv$high_selected_gene <- as.character(rv$high_df[input$dt_high_cell_clicked$row, "ensembl_id"])
      if (debug) message(sprintf("High df row clicked: %d (%s)",
                                 input$dt_high_cell_clicked$row,
                                 rv$high_selected_gene))
    }
  })

  # Update selected gene when row is clicked
  observeEvent({input$dt_inspect_cell_clicked}, {
    if (!is.null(input$dt_inspect_cell_clicked$row)) {
      rv$inspect_selected_gene <- as.character(rv$inspect_df[input$dt_inspect_cell_clicked$row, "ensembl_id"])
      if (debug) message(sprintf("Inspect df row clicked: %d (%s)",
                                 input$dt_inspect_cell_clicked$row,
                                 rv$inspect_selected_gene))
    }
  })

  output$gg_rapid_change <- renderPlot({
    if (rv$b[2] > rv$a[1]) {
      midpoint <- ((rv$a[2] + rv$b[1])/2) %% 100
    } else {
      midpoint <- ((rv$a[2] + rv$b[1])) %% 100
    }
    dat <- sample_df %>%
      mutate(exprs=exprs[rv$rapid_selected_gene, sample_id])
    x <- gene_info %>% filter(ensembl_id == rv$rapid_selected_gene)
    cycle_plot(dat=dat, title=sprintf("%s (%s)", x$symbol, x$ensembl_id),
               subtitle=ifelse(is.na(x$gene_name), "", x$gene_name)) +
      geom_vline(xintercept=midpoint,
                 linetype="dashed",
                 color="red")
  })
  output$gg_high <- renderPlot({
    dat <- sample_df %>%
      mutate(exprs=exprs[rv$high_selected_gene, sample_id])
    x <- gene_info %>% filter(ensembl_id == rv$high_selected_gene)
    g <- cycle_plot(dat=dat, title=sprintf("%s (%s)", x$symbol, x$ensembl_id),
                    subtitle=ifelse(is.na(x$gene_name), "", x$gene_name))

    if (rv$b[2] > rv$a[1]) {
      g +
        annotate("rect", xmin = rv$a[1], xmax = rv$b[2], ymin = Inf, ymax = -Inf,
                 fill="red", alpha = .2)
    } else {
      g +
        annotate("rect", xmin = 0, xmax = rv$b[2], ymin = Inf, ymax = -Inf,
                 fill="red", alpha = .2) +
        annotate("rect", xmin = rv$a[1], xmax = 100, ymin = Inf, ymax = -Inf,
                 fill="red", alpha = .2)
    }
  })
  output$gg_inspect <- renderPlot({
    dat <- sample_df %>%
      mutate(exprs=exprs[rv$inspect_selected_gene, sample_id])
    x <- gene_info %>% filter(ensembl_id == rv$inspect_selected_gene)
    cycle_plot(dat=dat, title=sprintf("%s (%s)", x$symbol, x$ensembl_id),
               subtitle=ifelse(is.na(x$gene_name), "", x$gene_name))
  })


  output$dt_rapid_change <- renderDataTable({
    datatable(rv$rapid_df %>% select(ensembl_id, diff, entrez_id, symbol, gene_name),
              options=list(pageLength=5),
              selection="single")
  })

  output$dt_high <- renderDataTable({
    datatable(rv$high_df %>% select(ensembl_id, diff, entrez_id, symbol, gene_name),
              options=list(pageLength=5),
              selection="single")
  })

  output$dt_inspect <- renderDataTable({
    datatable(rv$inspect_df %>% select(ensembl_id, R2, dev_exp, edf, entrez_id, symbol, gene_name),
              options=list(pageLength=5),
              selection="single")
  })

  output$rapid_text <- renderText({
    x <- c(rv$a[1], rv$a[2], rv$b[1], rv$b[2])
    if (all(x == round(x))) {
      pretty_range <- sprintf("[%d, %d] vs [%d, %d]", rv$a[1], rv$a[2], rv$b[1], rv$b[2])
    } else {
      pretty_range <- sprintf("[%0.1f, %0.1f] vs [%0.1f, %0.1f]", rv$a[1], rv$a[2], rv$b[1], rv$b[2])
    }
    sprintf("Rapidly changing genes between time ranges %s.", pretty_range)
  })

  output$high_text <- renderText({
    x <- c(rv$a[1], rv$a[2], rv$b[1], rv$b[2])
    if (all(x == round(x))) {
      pretty_range <- sprintf("[%d, %d]", rv$a[1], rv$b[2])
    } else {
      pretty_range <- sprintf("[%0.1f, %0.1f]", rv$a[1], rv$b[2])
    }
    sprintf("Highest expressed genes in time range %s compared to everywhere else.",
            pretty_range)
  })


} # End server

shinyApp(ui, server)
