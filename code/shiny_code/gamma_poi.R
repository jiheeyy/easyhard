library(shiny)
library(ggplot2)
library(gridExtra)

# Define UI
ui <- fluidPage(
  titlePanel("Bayesian Inference: Poisson(sλ)"),

  sidebarLayout(
    sidebarPanel(
      h4("Parameters"),
      helpText("Model: x ~ Poisson(sλ)"),

      # 1. New Slider for Exposure 's'
      sliderInput("s",
                  "Exposure  (s):",
                  min = 0.5, max = 10, value = 1, step = 0.1),

      sliderInput("x",
                  "Observed Count (x):",
                  min = 0, max = 50, value = 5, step = 0.1),

      hr(),
      h4("Prior Belief (Gamma)"),
      sliderInput("alpha",
                  "Shape (alpha):",
                  min = 0.1, max = 10, value = 2, step = 0.1),

      sliderInput("beta",
                  "Rate (beta):",
                  min = 0.1, max = 10, value = 1, step = 0.1),

      hr(),
      h5("Calculated Constants"),
      tableOutput("statsTable")
    ),

    mainPanel(
      plotOutput("bayesPlot", height = "700px")
    )
  )
)

# Define Server
server <- function(input, output) {

  # Reactive computations for the distributions
  dist_data <- reactive({
    obs_x <- input$x
    obs_s <- input$s
    p_alpha <- input$alpha
    p_beta <- input$beta

    # 1. Update Rules
    # Posterior Shape = alpha + x
    post_alpha <- p_alpha + obs_x
    # Posterior Rate = beta + s (The exposure adds to the rate)
    post_beta <- p_beta + obs_s

    # 2. Define Lambda Range for plotting
    # We find a range that covers both Prior and Posterior comfortably
    lambda_max <- max(
      qgamma(0.999, shape = p_alpha, rate = p_beta),
      qgamma(0.999, shape = post_alpha, rate = post_beta)
    )
    lambda_seq <- seq(0, lambda_max, length.out = 500)

    # 3. Compute Densities
    prior_y <- dgamma(lambda_seq, shape = p_alpha, rate = p_beta)
    post_y  <- dgamma(lambda_seq, shape = post_alpha, rate = post_beta)

    # Likelihood: Poisson(x | s * lambda)
    # Note: dpois expects the mean parameter. Mean = s * lambda
    lik_raw <- dpois(obs_x, lambda = obs_s * lambda_seq)

    # Scaling Likelihood for visualization (matching Posterior peak)
    scale_factor <- if(max(lik_raw) > 0) max(post_y) / max(lik_raw) else 1
    lik_y <- lik_raw * scale_factor

    list(
      lambda_seq = lambda_seq,
      prior = prior_y,
      posterior = post_y,
      likelihood = lik_y,
      post_alpha = post_alpha,
      post_beta = post_beta
    )
  })

  # Reactive computation for Marginal (Negative Binomial)
  marginal_data <- reactive({
    obs_x <- input$x
    obs_s <- input$s
    p_alpha <- input$alpha
    p_beta <- input$beta

    # Negative Binomial parameters with scaling s
    # Probability of success p = beta / (beta + s)
    nb_prob <- p_beta / (p_beta + obs_s)

    # X Range for plotting
    x_max <- qnbinom(0.995, size = p_alpha, prob = nb_prob)
    # Ensure range is large enough to show the observation
    x_limit <- max(30, x_max, obs_x + 5)
    x_seq <- 0:x_limit

    marginal_prob <- dnbinom(x_seq, size = p_alpha, prob = nb_prob)

    list(
      x_seq = x_seq,
      probs = marginal_prob,
      obs_prob = dnbinom(obs_x, size = p_alpha, prob = nb_prob),
      nb_prob_param = nb_prob
    )
  })

  output$bayesPlot <- renderPlot({
    d <- dist_data()
    m <- marginal_data()

    # --- PLOT 1: Lambda Space ---
    df_lambda <- data.frame(
      Lambda = rep(d$lambda_seq, 3),
      Density = c(d$prior, d$likelihood, d$posterior),
      Type = rep(c("Prior", "Likelihood (Scaled)", "Posterior"), each = length(d$lambda_seq))
    )
    df_lambda$Type <- factor(df_lambda$Type, levels = c("Prior", "Likelihood (Scaled)", "Posterior"))

    p1 <- ggplot(df_lambda, aes(x = Lambda, y = Density, color = Type, fill = Type)) +
      geom_line(size = 1.2) +
      geom_area(alpha = 0.2, position = "identity") +
      scale_color_manual(values = c("Prior"="#56B4E9", "Likelihood (Scaled)"="#E69F00", "Posterior"="#009E73")) +
      scale_fill_manual(values = c("Prior"="#56B4E9", "Likelihood (Scaled)"="#E69F00", "Posterior"="#009E73")) +
      labs(title = "1. Inference Space (Lambda)",
           subtitle = paste0("Prior: Gamma(", input$alpha, ", ", input$beta, ") -> ",
                             "Posterior: Gamma(", d$post_alpha, ", ", d$post_beta, ")"),
           y = "Density", x = expression(lambda)) +
      theme_minimal() +
      theme(legend.position = "top", text = element_text(size=14))

    # --- PLOT 2: X Space ---
    df_x <- data.frame(x = m$x_seq, prob = m$probs)

    p2 <- ggplot(df_x, aes(x = x, y = prob)) +
      geom_col(fill = "gray70", color = "white") +
      geom_col(data = subset(df_x, x == input$x), fill = "#D55E00", width = 0.8) +
      labs(title = "2. Data Space (Marginal P(x))",
           subtitle = paste0("NB(r=", input$alpha, ", p=", round(m$nb_prob_param, 2),
                             ") | Observed P(x=", input$x, ") = ", round(m$obs_prob, 4)),
           y = "Probability", x = "Count (x)") +
      theme_minimal() +
      theme(text = element_text(size=14))

    grid.arrange(p1, p2, nrow = 2, heights = c(1.8, 1))
  })

  output$statsTable <- renderTable({
    d <- dist_data()
    m <- marginal_data()
    obs_s <- input$s

    data.frame(
      Metric = c("Prior Mean E[λ]",
                 "Expected Data E[x] = E[sλ]",
                 "Posterior Mean E[λ | x]"),
      Value = c(input$alpha / input$beta,
                obs_s * (input$alpha / input$beta),
                d$post_alpha / d$post_beta)
    )
  })
}

shinyApp(ui = ui, server = server)
