#' @title Backtest multiple portfolios over multiple datasets of stock prices in a rolling-window basis
#'
#' @description Automated backtesting of multiple portfolios over multiple
#' datasets of stock prices in a rolling-window fashion.
#' Each portfolio design is easily defined as a
#' function that takes as input a window of the stock prices and outputs the
#' portfolio weights. Multiple portfolios can be easily specified as a list
#' of functions or as files in a folder. Multiple datasets can be conveniently
#' obtained with the function \code{\link{financialDataResample}} that resamples
#' the data downloaded with the function \code{\link{stockDataDownload}}.
#' The results can be later assessed and arranged with tables and plots.
#' The backtesting can be highly time-consuming depending on the number of
#' portfolios and datasets can be performed with parallel computation over
#' multiple cores. Errors in functions are properly catched and handled so
#' that the execution of the overal backtesting is not stopped (error messages
#' are stored for debugging purposes). See
#' \href{https://CRAN.R-project.org/package=portfolioBacktest/vignettes/PortfolioBacktest.html}{vignette}
#' for a detailed explanation.
#'
#' @param portfolio_funs List of functions (can also be a single function), each of them taking as input
#'                       a dataset containing a list of \code{xts} objects (following the format of each
#'                       element of the argument \code{dataset_list}) properly windowed (following the
#'                       rolling-window approach) and returning the portfolio as a vector of normalized
#'                       weights. See
#'                       \href{https://CRAN.R-project.org/package=portfolioBacktest/vignettes/PortfolioBacktest.html}{vignette}
#'                       for details.
#' @param dataset_list List of datasets, each containing a list of \code{xts} objects, as generated
#'                     by the function \code{\link{financialDataResample}}.
#' @param folder_path If \code{portfolio_funs} is not defined, this should contain the path to a folder
#'                    containing the portfolio functions saved in files. See
#'                    \href{https://CRAN.R-project.org/package=portfolioBacktest/vignettes/PortfolioBacktest.html}{vignette}
#'                    for details.
#' @param source_to_local Logical value indicating whether to source files to local environment (default is \code{TRUE}).
#'                        It might be dangerous to set it to \code{FALSE} as in such case the global environment may be changed.
#'                        We suggest only to allow \code{FALSE} when the code in the source files does not work when locally
#'                        sourced, e.g., with some versions of package \code{CVXR}. In that case, we further recommend to set
#'                        \code{paral_portfolios > 1} to avoid changing the global environment.
#' @param price_name Name of the \code{xts} object in each dataset that contains the prices to be used in the portfolio return
#'                   computation (default is the name of the first \code{xts} object).
#' @param paral_portfolios Interger indicating number of portfolios to be evaluated in parallel (default is \code{1}),
#'                         see \href{https://CRAN.R-project.org/package=portfolioBacktest/vignettes/PortfolioBacktest.html#parallel-backtesting}{vignette-paralle-mode} for details.
#' @param paral_datasets Interger indicating number of datasets to be evaluated in parallel (default is \code{1}),
#'                        see \href{https://CRAN.R-project.org/package=portfolioBacktest/vignettes/PortfolioBacktest.html#parallel-backtesting}{vignette-paralle-mode} for details.
#' @param show_progress_bar Logical value indicating whether to show progress bar (default is \code{FALSE}).
#' @param benchmarks String vector indicating the benchmark portfolios to be incorporated, currently supports:
#' \itemize{\item{\code{1/N} - the 1/N portfolio, \eqn{w = [1/N, ..., 1/N]} with \eqn{N} be number of stocks;}
#'          \item{\code{IVP} - the inverse-volatility portfolio, with weights be inversely proportional the standard deviation of returns;}
#'          \item{\code{index} - the market index, requires an \code{xts} named `index` in the datasets.}}
#' @param shortselling Logical value indicating whether shortselling is allowed or not
#'                     (default is \code{TRUE}, so no control for shorselling in the backtesting).
#' @param leverage Amount of leverage as in \eqn{||w||_1 <= leverage}
#'                 (default is \code{Inf}, so no control for leverage in the backtesting).
#' @param lookback Length of the lookback rolling window in periods (default is \code{252}).
#' @param T_rolling_window Deprecated: use \code{lookback} instead.
#' @param optimize_every How often the portfolio is to be optimized in periods (default is \code{20}).
#' @param rebalance_every How often the portfolio is to be rebalanced in periods (default is \code{1}).
#' @param bars_per_year Number of bars/periods per year. By default it will be calculated automatically
#'                      (e.g., for daily data there are 252 bars/periods per year).
#' @param execution String that can be either \code{"same period"} (default) or \code{"next period"}.
#'                  At the rebalancing period \code{t}, the portfolio has used information up to (and including)
#'                  period \code{t}. Same period execution means one can get into the position at that period \code{t},
#'                  whereas the next period execution means that one can only get into the position the following period.
#' @param cost List containing four different types of transaction costs (common for all assets)
#'             for buying, selling, shorting, and long leveraging. The default is
#'             \code{cost = list(buy = 0e-4, sell = 0e-4, short = 0e-4, long_leverage = 0e-4)}.
#'             If some elements are not specified then they will be automatically set to zero.
#' @param cpu_time_limit Time limit for executing each portfolio function over a single data set
#'                       (default is \code{Inf}, so no time limit).
#' @param return_portfolio Logical value indicating whether to return the portfolios (default is \code{TRUE}).
#'                         Three portfolio series are returned:
#'                         \code{w_optimized} is the optimized portfolio at each given optimization period
#'                         (using all the information up to and including that period, which can be executed
#'                         either on the same period or the following period),
#'                         \code{w_rebalanced} is the rebalanced portfolio at each given rebalancing period,
#'                         and \code{w_bop} is the "beginning-of-period" portfolio (i.e., at each period it contains
#'                         the weights held in the market in the previous period so that the portfolio return at
#'                         that period is just the product of the asset returns and \code{w_bop} at that period.)
#' @param return_returns Logical value indicating whether to return the portfolio returns (default is \code{TRUE}).
#'                       Three series are returned:
#'                       \code{return} with the portfolio returns,
#'                       \code{wealth} with the portfolio wealth (aka cumulative P&L), and
#'                       \code{X_lin} with the returns of the assets in the universe (note that the portfolio returns
#'                       can also be obtained as \code{rowSums(X_lin * w_bop)} in the absence of transaction costs).
#'
#' @return List with the portfolio backtest results, see
#'         \href{https://CRAN.R-project.org/package=portfolioBacktest/vignettes/PortfolioBacktest.html#result-format}{vignette-result-format}
#'         for details. It can be accessed directly, but we highly recommend the use of the package specific functions to extract
#'         any required information, namely, \code{\link{backtestSelector}}, \code{\link{backtestTable}},
#'         \code{\link{backtestBoxPlot}}, \code{\link{backtestLeaderboard}},
#'         \code{\link{backtestSummary}}, \code{\link{summaryTable}}, \code{\link{summaryBarPlot}}.
#'
#' @author Daniel P. Palomar and Rui Zhou
#'
#' @seealso \code{\link{stockDataDownload}}, \code{\link{financialDataResample}},
#'          \code{\link{backtestSelector}}, \code{\link{backtestTable}},
#'          \code{\link{backtestBoxPlot}}, \code{\link{backtestLeaderboard}},
#'          \code{\link{backtestSummary}}, \code{\link{summaryTable}}, \code{\link{summaryBarPlot}}.
#'
#' @examples
#' \donttest{
#' library(portfolioBacktest)
#' data(dataset10) # load dataset
#'
#' # define your own portfolio function
#' ewp_portfolio <- function(dataset, ...) {
#'   N <- ncol(dataset$adjusted)
#'   return(rep(1 / N, N))
#' }
#'
#' # do backtest
#' bt <- portfolioBacktest(list("EWP" = ewp_portfolio), dataset10)
#'
#' # check your result
#' names(bt)
#' backtestSelector(bt, portfolio_name = "EWP", measures = c("Sharpe ratio", "max drawdown"))
#' backtestTable(bt, measures = c("Sharpe ratio", "max drawdown"))
#' bt_summary <- backtestSummary(bt)
#' summaryTable(bt_summary)
#' }
#'
#' @import xts
#' @importFrom zoo index
#' @importFrom pbapply pblapply pboptions
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom utils object.size
#' @export
portfolioBacktest <- function(portfolio_funs = NULL, dataset_list, folder_path = NULL, source_to_local = TRUE, price_name = NULL,
                              paral_portfolios = 1, paral_datasets = 1, aum = 1,
                              show_progress_bar = FALSE, benchmarks = NULL,
                              shortselling = TRUE, leverage = Inf,
                              lookback = 252, T_rolling_window = NULL, optimize_every = 20, rebalance_every = 1, bars_per_year = 252,
                              execution = c("same period", "next period"),
                              cost = list(buy = 0e-4, sell = 0e-4, short = 0e-4, long_leverage = 0e-4),
                              cpu_time_limit = Inf,
                              return_portfolio = TRUE, return_returns = TRUE) {
  ####### error control ########
  paral_portfolios <- round(paral_portfolios)
  paral_datasets <- round(paral_datasets)
  if (paral_portfolios < 1 || paral_datasets < 1) stop("Parallel number must be a positive interger.")
  if (is.null(folder_path) && is.null(portfolio_funs)) stop("The \"folder_path\" and \"portfolio_fun_list\" cannot be both NULL.")
  if (!is.null(portfolio_funs) && !is.list(portfolio_funs)) portfolio_funs <- list(portfolio_funs)
  cost <- modifyList(list(buy = 0, sell = 0, short = 0, long_leverage = 0), cost)
  if (length(cost) != 4) stop("Problem in specifying the cost: the elements can only be buy, sell, short, and long_leverage.")
  if (is.xts(dataset_list[[1]])) stop("Each element of \"dataset_list\" must be a list of xts objects. Try to surround your passed \"dataset_list\" with list().")
  if (is.null(price_name)) {
    price_name <- names(dataset_list[[1]])[1]
  }
  if (is.null(price_name)) {
    stop("price_name incorrectly specified (or first xts object in each dataset does not have a name).")
  }
  if (!(price_name %in% names(dataset_list[[1]]))) stop("Price data xts element \"", price_name, "\" does not exist in dataset_list.")
  if (!is.xts(dataset_list[[1]][[1]])) stop("prices have to be xts.")
  if (!is.null(T_rolling_window)) {
    stop("Argument ", dQuote("T_rolling_window"), " is deprecated. Use instead ", dQuote("lookback"))
  }
  ##############################

  if (is.null(names(dataset_list))) {
    names(dataset_list) <- paste0("data", 1:length(dataset_list))
  }

  if (!show_progress_bar) {
    pbo <- pboptions(type = "none")
    on.exit(pboptions(pbo), add = TRUE)
  }


  # parse arguments for subsequent calls
  args <- list(
    dataset_list = dataset_list, price_name = price_name,
    paral_datasets = paral_datasets, aum = aum,
    shortselling = shortselling, leverage = leverage,
    lookback = lookback, optimize_every = optimize_every, rebalance_every = rebalance_every,
    bars_per_year = bars_per_year, execution = execution, cost = cost,
    cpu_time_limit = cpu_time_limit, return_portfolio = return_portfolio, return_returns = return_returns
  )

  # when portfolio_funs is passed
  if (!is.null(portfolio_funs)) {
    if (is.null(names(portfolio_funs))) {
      portfolio_names <- paste0("fun", 1:length(portfolio_funs))
    } else {
      portfolio_names <- names(portfolio_funs)
    }
    if ("" %in% portfolio_names) stop("Each element of \"portfolio_funs\" must has a unique name.")
    if (length(portfolio_names) != length(unique(portfolio_names))) stop("\"portfolio_funs\" contains repeated names.")

    message(sprintf(
      "Backtesting %d portfolios over %d datasets (periodicity = %s data)...",
      length(portfolio_names), length(dataset_list), periodicity(dataset_list[[1]][[1]])$scale
    ))

    singlePortfolioBacktest_safe <- function(...) {
      initial_packages <- search() # snapshot of initial packages
      res <- do.call(singlePortfolioBacktest, args = list(...))
      detachPackages(setdiff(search(), initial_packages)) # detach the newly loaded packages
      return(res)
    }

    if (paral_portfolios == 1) {
      result <- list()
      for (i in 1:length(portfolio_funs)) {
        if (show_progress_bar) {
          message("  Backtesting function \"", format(portfolio_names[i], width = 15), "\" (", i, "/", length(portfolio_names), ")")
        }
        result[[i]] <- do.call(singlePortfolioBacktest_safe, args = c(list(portfolio_fun = portfolio_funs[[i]]), args))
      }
    } else {
      cl <- makeCluster(paral_portfolios)
      # export global variables
      exports <- setdiff(ls(envir = .GlobalEnv), c("portfolio_fun", "dataset_list", "show_progress_bar"))
      exports <- exports[sapply(exports, function(x) object.size(get(x, envir = .GlobalEnv))) < 10 * 1024^2] # remove large objects of more than 10 Mb
      clusterExport(cl = cl, varlist = exports, envir = .GlobalEnv)
      # export attached packages
      attached_packages <- .packages()
      clusterExport(cl = cl, varlist = "attached_packages", envir = environment())
      invisible(clusterEvalQ(cl = cl, expr = sapply(attached_packages, function(pkg) require(pkg, character.only = TRUE))))
      # apply!
      result <- pblapply(portfolio_funs,
        function(portfolio_fun) {
          do.call(singlePortfolioBacktest_safe, args = c(list(portfolio_fun = portfolio_fun), args))
        },
        cl = cl
      )
      stopCluster(cl)
    }
  } else { # when folder_path is passed
    files <- list.files(folder_path)
    portfolio_names <- gsub(".R", "", files)
    if (length(portfolio_names) == 0) {
      stop(sprintf("Could not find any .R files in folder: %s", folder_path))
    }

    message(sprintf(
      "Backtesting %d portfolios over %d datasets (periodicity = %s data)...",
      length(portfolio_names), length(dataset_list), periodicity(dataset_list[[1]][[1]])$scale
    ))

    singlePortfolioBacktest_folder_safe <- function(...) {
      args <- list(...)
      initial_packages <- search() # snapshot of initial packages

      error_message <- tryCatch(
        {
          # load the file with portfolio_fun()
          suppressMessages(source(file.path(args$folder_path, args$file), local = args$source_to_local))
          # fix the arguments
          if (args$source_to_local && exists("portfolio_fun", envir = environment())) {
            # args <- c(list(portfolio_fun = environment()$portfolio_fun), args)
            args <- c(list(portfolio_fun = get(x = "portfolio_fun", envir = environment())), args)
          }
          if (!args$source_to_local && exists("portfolio_fun", envir = parent.env(environment()))) {
            # args <- c(list(portfolio_fun = parent.env(environment())$portfolio_fun), args)
            args <- c(list(portfolio_fun = get(x = "portfolio_fun", envir = parent.env(environment()))), args)
          }
          args$folder_path <- args$file <- args$source_to_local <- NULL # remove these arguments
          # do backtest
          res <- do.call(singlePortfolioBacktest, args)
          NULL
        },
        # warning = function(w) return(ifelse(!is.null(w$message), w$message, "")),
        error = function(e) {
          return(ifelse(!is.null(e$message), e$message, ""))
        }
      )
      detachPackages(setdiff(search(), initial_packages)) # detach the newly loaded packages
      if (!is.null(error_message)) {
        return(list(source_error_message = error_message))
      } else {
        return(res)
      }
    }

    if (paral_portfolios == 1) {
      result <- list()
      for (i in 1:length(files)) {
        if (show_progress_bar) {
          message("\n Backtesting file \"", format(portfolio_names[i], width = 15), "\" (", i, "/", length(portfolio_names), ")")
        }
        result[[i]] <- do.call(singlePortfolioBacktest_folder_safe, args = c(list(
          folder_path = folder_path,
          file = files[i],
          source_to_local = source_to_local
        ), args))
      }
    } else {
      cl <- makeCluster(paral_portfolios)

      result <- pblapply(files,
        function(file) {
          do.call(singlePortfolioBacktest_folder_safe, args = c(list(
            folder_path = folder_path,
            file = file,
            source_to_local = source_to_local
          ), args))
        },
        cl = cl
      )
      stopCluster(cl)
    }
  }
  names(result) <- portfolio_names
  attr(result, "portfolio_index") <- 1:length(portfolio_names)


  #
  # add benchmarks if necessary
  #
  if (!is.null(benchmarks)) {
    message("Backtesting benchmarks...")
    benchmark_portfolios <- benchmark_library[names(benchmark_library) %in% benchmarks]
    res_benchmarks <- list()
    for (i in seq_along(benchmark_portfolios)) {
      if (show_progress_bar) {
        message("  Backtesting benchmark \"", format(benchmarks[i], width = 15), "\" (", i, "/", length(benchmarks), ")")
      }
      # change order of xts elements
      new_order_names <- c(price_name, setdiff(names(args$dataset_list[[1]]), price_name))
      args$dataset_list <- lapply(args$dataset_list, function(xts_list) xts_list[new_order_names])
      res_benchmarks[[i]] <- do.call(singlePortfolioBacktest, args = c(list(portfolio_fun = benchmark_portfolios[[i]]), args))
    }
    names(res_benchmarks) <- names(benchmark_portfolios)

    # add index if needed
    if ("index" %in% benchmarks) {
      if (show_progress_bar) {
        message("  Backtesting \"index     ", format("", width = 15), "\" (", i + 1, "/", length(benchmarks), ")")
      }
      args$price_name <- "index"
      args$optimize_every <- 1e10
      args$rebalance_every <- 1e10
      res_benchmarks$index <- do.call(singlePortfolioBacktest, args = c(list(portfolio_fun = function(...) {
        return(1)
      }), args))
    }
    result <- c(result, res_benchmarks)
    attr(result, "contain_benchmark") <- TRUE
    attr(result, "benchmark_index") <- (1:length(result))[-c(1:length(portfolio_names))]
  } else {
    attr(result, "contain_benchmark") <- FALSE
  }

  return(result)
}


singlePortfolioBacktest <- function(...) {
  # arrange arguments
  args <- list(...)
  paral_datasets <- args$paral_datasets
  dataset_list <- args$dataset_list
  args$paral_datasets <- args$dataset_list <- NULL

  # do backtests
  if (paral_datasets == 1) {
    result <- pblapply(dataset_list, function(data) do.call(singlePortfolioSingleXTSBacktest, args = c(list(data = data), args)))
  } else {
    cl <- makeCluster(paral_datasets)
    # export global variables
    exports <- setdiff(ls(envir = .GlobalEnv), c("portfolio_fun", "data"))
    exports <- exports[sapply(exports, function(x) object.size(get(x, envir = .GlobalEnv))) < 10 * 1024^2] # remove large objects of more than 10 Mb
    clusterExport(cl = cl, varlist = exports, envir = .GlobalEnv)
    # export attached packages
    attached_packages <- .packages()
    clusterExport(cl = cl, varlist = "attached_packages", envir = environment())
    invisible(clusterEvalQ(cl = cl, expr = sapply(attached_packages, function(pkg) require(pkg, character.only = TRUE))))
    # apply!
    result <- pblapply(dataset_list,
      function(data) do.call(singlePortfolioSingleXTSBacktest, args = c(list(data = data), args)),
      cl = cl
    )
    stopCluster(cl)
  }

  attr(result, "params") <- attr(args$portfolio_fun, "param") # in case portfolio_fun had params as attribute (for parameter tuning)
  return(result)
}

# Backtesting of one portfolio function on one single xts
# Reoptimizing is how often the portfolio gets to be recomputed (your portfolio function will be called).
# But then you could possibly rebalance more frequently than you reoptimize your portfolio
# this is because naturally your initial portfolio will change over time as the prices change,
# so you may want to reset the portfolio to the one you initially designed again
#
#' @import xts
singlePortfolioSingleXTSBacktest <- function(portfolio_fun, data, price_name,
                                             shortselling, leverage, aum,
                                             lookback, optimize_every, rebalance_every, bars_per_year,
                                             execution = c("same period", "next period"),
                                             cost = list(buy = 0, sell = 0, short = 0, long_leverage = 0),
                                             cpu_time_limit,
                                             return_portfolio, return_returns) {
  ######## error control  #########
  if (!price_name %in% names(data)) {
    stop("Failed to find price data with name \"", price_name, "\"", " in given dataset_list.")
  }
  if (!"volumes" %in% names(data)) {
    stop("Failed to find volumes data in given dataset_list.")
  }
  if (!"slippage" %in% names(data)) {
    stop("Failed to find slippage data in given dataset_list.")
  }
  if (!"vols" %in% names(data)) {
    stop("Failed to find vols data in given dataset_list.")
  }
  prices <- data[[price_name]]
  volumes <- data[["volumes"]]
  slippage <- data[["slippage"]]
  vols <- data[["vols"]]

  if (is.list(prices)) stop("prices have to be xts, not a list, make sure you index the list with double brackets [[.]].")
  if (!is.xts(prices)) stop("prices have to be xts.")
  N <- ncol(prices)
  T <- nrow(prices)
  if (lookback >= T) stop("T is not large enough for the given lookback window length.")
  if (optimize_every %% rebalance_every != 0) stop("The reoptimization period has to be a multiple of the rebalancing period.")
  if (!is.function(portfolio_fun)) stop("portfolio_fun is not a function.")
  if (tzone(prices) != "UTC") tzone(prices) <- "UTC"
  #################################

  delay <- switch(match.arg(execution), # w[t] used info up to (including) price[t]
    "same period" = 1, # w[t] is (ideally) executed at price[t], so will multiply return[t+1]
    "next period" = 2, # w[t] is executed one period later at price[t+1], so will multiply return[t+2]
    stop("Execution method unknown")
  )

  # indices
  optimize_indices <- seq(from = lookback, to = T - delay, by = optimize_every)
  rebalance_indices <- seq(from = lookback, to = T - delay, by = rebalance_every)
  if (any(!(optimize_indices %in% rebalance_indices))) {
    stop("The reoptimization indices have to be a subset of the rebalancing indices.")
  }

  # create variables (w and cash are always normalized wrt NAV_bop)
  compute_tc <- any(cost != 0)
  tc <- 0
  cpu_time <- c()

  X <- prices # meillä prices on jo tuotot

  w_designed <- empty_xts_from(prices)
  w_bop <- empty_xts_from(prices)
  delta_bop <- empty_xts_from(prices)
  w_eop <- empty_xts_from(prices)
  NAV_eop <- empty_xts_from(prices[, 1], colnames = "NAV")

  # initial values
  initial_budget <- aum
  w_eop[lookback, ] <- 0 # w_eop but normalized wrt NAV_bop
  cash_eop <- 1 # all in cash
  NAV_eop[lookback, ] <- (sum(w_eop[lookback, ]) + cash_eop) * initial_budget
  num_rebalances <- 0

  # loop over time
  for (t in lookback:T) {
    # update of w_bop as w_eop
    if (t > lookback) {
      NAV_bop <- as.numeric(NAV_eop[t - 1])
      if (all(is.na(w_bop[t, ]))) {
        w_bop[t, ] <- w_eop[t - 1, ]
      }
      # include tc in cash
      w_prev <- w_eop[t - 1, ]
      w_cur <-  w_bop[t, ]
      
      w_prev[is.na(w_prev)] <- 0 #jos aikaisempi paino on NA laitetaan nollaksi
      w_cur[is.na(w_cur)] <- 0 #jos nykyinen paino on NA laitetaan nollaksi
      
      delta_bop[t, ] <- as.numeric(w_cur) - as.numeric(w_prev)
      
      if (t == lookback + delay) { # this is to avoid the initial huge turnover
        delta_bop[t, ] <- 0
      }
      if (compute_tc) {
        
        #pitaa käsitellä tilanteet
        
        #1. osakkeella ON paino edellisellä salkussa, nyt EI signaalia, osake myydään/ostetaan pois sillä NYKYINEN paino on nolla
        #2. osakkeella EI painoa edellisellä salkussa, nyt ON signaali, osake myydään/ostetaan sillä EDELLINEN paino on nolla
        #   jos ei ole volyymiä tai volaa tai slipagea niin cross section percentiili 0.75 kaikkiin?
    
        adv_pct <- (abs(delta_bop[t, ]) * NAV_bop) / volumes[t, ]
        
        mkt_impact_cost <- mkt_impact(vols[t, ], adv_pct, c = 1)
        slippage_cost <- slippage[t, ]
        
        mkt_impact_cost[is.na(mkt_impact_cost)] <- quantile(mkt_impact_cost, 0.75, na.rm = TRUE)
        slippage_cost[is.na(slippage_cost)] <- quantile(slippage_cost, 0.75, na.rm = TRUE)
        
        mkt_impact_cost <- sum(abs(delta_bop[t, ]) * mkt_impact_cost, na.rm = TRUE)
        
        slippage_cost <- sum(abs(delta_bop[t, ]) * slippage_cost, na.rm = TRUE)

        financing_costs <- cost$long_leverage * max(0, sum(pos(w_bop[t, ]), na.rm = TRUE) - 1) + cost$short * sum(pos(-w_bop[t, ]), na.rm = TRUE)

        tc <- slippage_cost + mkt_impact_cost + financing_costs + 0.0001*sum(abs(delta_bop[t, ]), na.rm = TRUE)
        
      } # borrowing cost
      cash_eop_ <- 1 - sum(w_bop[t, ], na.rm = TRUE) - tc
      
      # include period return in w
      # muutetaan NA tuotot nolliksi
      p_ret <- X[t, ]
      p_ret[is.na(p_ret)] <- 0
      
      w_eop_ <- as.numeric(1 + p_ret) * w_bop[t, ]
      # normalize
      NAV_eop[t] <- (sum(w_eop_, na.rm = TRUE) + cash_eop_) * NAV_bop
      if (as.numeric(NAV_eop[t]) <= 0) { # if bankruptcy
        break
      }
      w_eop[t, ] <- w_eop_ / (sum(w_eop_, na.rm = TRUE) + cash_eop_)
    }

    # rebalance and possibly design of w_designed
    if (t %in% rebalance_indices && t + delay <= T) {
      w_current <- as.matrix(w_eop[t, ])[1, ]
      if (t %in% optimize_indices) {
        # design portfolio
        data_window <- lapply(data, function(x) x[(t - lookback + 1):t, ])
        start_time <- proc.time()[3]
        error_capture <- R.utils::withTimeout(
          expr = evaluate::try_capture_stack(
            last_w_optimized <- do.call(portfolio_fun, list(data_window, w_current = w_current)),
            env = environment()
          ),
          timeout = cpu_time_limit, onTimeout = "silent"
        )
        cpu_time <- c(cpu_time, as.numeric(proc.time()[3] - start_time))
        
        #ne stokit joille signaali on NA, tulee 0 paino
        last_w_optimized[is.na(last_w_optimized)] <- 0
        
        portf <- check_portfolio_errors(error_capture, last_w_optimized, shortselling, leverage)
        if (portf$error) {
          break
        } # return immediately if error happens
      }
     
      w_designed[t, ] <- last_w_optimized
      w_bop[t + delay, ] <- last_w_optimized
      if (!isTRUE(all.equal(last_w_optimized, w_current))) {
        num_rebalances <- num_rebalances + 1
      }
    }
  }

  # in case of error return now
  if (portf$error) {
    return(list(
      performance = portfolioPerformance(rets = NA, bars_per_year = NA),
      cpu_time = NA,
      error = TRUE,
      error_message = portf$error_message
    ))
  }

  # in case of no error, continue normally
  w_optimized <- w_designed[optimize_indices, ]
  w_rebalanced <- w_designed[rebalance_indices, ]
  w_bop <- w_bop[(lookback + delay):T, ]
  delta_bop <- delta_bop[(lookback + delay):T, ]
  w_eop <- w_eop[(lookback + delay - 1):T, ]
  NAV_eop <- NAV_eop[(lookback + delay - 1):T, ]
  returns_eop <- PerformanceAnalytics::CalculateReturns(NAV_eop)[-1]
  colnames(returns_eop) <- "return"

  # compute ROT based on normalized dollars
  sum_turnover_norm <- sum(abs(delta_bop), na.rm = TRUE)
  sum_PnL_norm <- sum(returns_eop[min(1 + rebalance_every, nrow(returns_eop)):nrow(returns_eop)]) # this subsetting is because the initial huge turnover was removed
  # sum_PnL_norm <- sum(returns_eop[1:(nrow(returns_eop) - rebalance_every)])
  ROT_bps <- ifelse(sum_turnover_norm != 0, 1e4 * sum_PnL_norm / sum_turnover_norm, NA)

  # return
  res <- list(
    performance = portfolioPerformance(
      rets = returns_eop, bars_per_year,
      rebalances_per_period = num_rebalances / (T - lookback),
      turnover_per_period = sum_turnover_norm / (T - lookback),
      ROT_bps = ROT_bps
    ),
    cpu_time = mean(cpu_time),
    error = FALSE,
    error_message = NA
  )
  if (return_portfolio) {
    res$w_optimized <- w_optimized
    res$w_rebalanced <- w_rebalanced
    res$w_bop <- w_bop
  }
  if (return_returns) {
    res$return <- returns_eop
    res$wealth <- NAV_eop
    res$X_lin <- X
  }
  return(res)
}

empty_xts_from <- function(orig_xts, fill = NA, colnames = NULL) {
  empty_xts <- xts(matrix(fill, nrow(orig_xts), ncol(orig_xts)),
    order.by = index(orig_xts)
  )
  colnames(empty_xts) <- if (is.null(colnames)) {
    colnames(orig_xts)
  } else {
    colnames
  }
  return(empty_xts)
}

mkt_impact <- function(vol_pct, adv_pct, c = 1) {
  pmax(pmin(c * vol_pct * sqrt(adv_pct), 0.1),0)
}

# tibble(vol_pct = runif(5000, min = 0.005, max = 0.03),
#        adv_pct = runif(5000, min = 0.005, max = 0.1)) %>% 
#   mutate(mkt_impact = mkt_impact(vol_pct = vol_pct, adv_pct = adv_pct)) %>% 
#   ggplot(aes(adv_pct, mkt_impact, color=vol_pct))+
#   geom_point()+
#   geom_smooth()


pos <- function(x) {
  pmax(0, as.numeric(x))
}

# cost <- list(buy = 17*10^(-4), sell = 17*10^(-4), long = 1*10^(-4), short = 1*10^(-4))  # for HK
# cost <- list(buy = 1*10^(-4), sell = 1*10^(-4), long = 1*10^(-4), short = 1*10^(-4))  # for US
# cost <- list(buy = 15*10^(-4), sell = 15*10^(-4), long = 15*10^(-4), short = 15*10^(-4))  # for CH

check_portfolio_errors <- function(error_capture, w, shortselling, leverage) {
  res <- list(error = FALSE, error_message = NA)

  # 1) check code execution errors
  if (is.list(error_capture)) {
    res$error <- TRUE
    res$error_message <- error_capture$message
    attr(res$error_message, "error_stack") <- list(
      "at" = deparse(error_capture$call),
      "stack" = paste(sapply(error_capture$calls[-1], deparse), collapse = "\n")
    )
  } else {
    # 2) check NA portfolio
    if (anyNA(w)) {
      res$error <- TRUE
      res$error_message <- "Returned portfolio contains NA."
    } else {
      # 3) check constraints
      if (!shortselling && any(w < -1e-6)) {
        res$error <- TRUE
        res$error_message <- "No-shortselling constraint not satisfied."
      }
      if (sum(abs(w)) > leverage + 1e-6) {
        res$error <- TRUE
        res$error_message <- c(res$error_message, "Leverage constraint not satisfied.")
      }
    }
  }
  return(res)
}
