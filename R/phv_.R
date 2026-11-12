#' Pharmacovigilance Analysis
#'
#' Pharmacovigilance, also known as drug safety. In the context of
#' pharmacovigilance studies, disproportionality analysis primarily served as a
#' tool to evaluate possible association between a specific adverse event and a
#' particular drug which can then be investigated through clinical assessment of
#' individual case reports.
#'
#' @param a also referred to as `n11` as this is the count of event of interest
#' under exposure of interest.
#' @param b also referred to as `n10` as this is the count of \emph{not} event
#' of interest under exposure of interest.
#' @param c also referred to as `n01` as this is the count of event of interest
#' under \emph{not} exposure of interest.
#' @param d also referred to as `n00` as this is the count of \emph{not} event
#' of interest under \emph{not} exposure of interest.
#' @param methods An atomic character, specifies the method used to signal
#' mining. Currently, only "ror", "prr", "chisq", "bcpnn_norm", "bcpnn_mcmc",
#' "obsexp_shrink", "fisher", and "ebgm" are supported. If `NULL`, all supported
#' methods will be used.
#' @param alpha Level of significance, for construction of the confidence
#'  intervals.
#' @inheritParams BiocParallel::bpmapply
#' @details
#' Note that the `a`, `b`, `c`, `d` inputs can be an atomic vectors of equal
#' length, for which the function will perform the calculations for each
#' individual `(a,b,c,d)`-tuple moving across the vectors.
#'
#' It is assumed that the contingency table under consideration has
#' drugs/exposures in the rows and outcomes/events in the columns. See
#' contingency table section.
#'
#' We use the `distinct` patient count method to obtain the frequency counts of
#' patients exposed to each interested drug, those reporting interested event.
#' As illustrated in the Contingency table, `n` equals the total number of
#' patients in the database, `n11` is the number of patients with exposure to
#' the interested drug during the model period and reporting interested events,
#' `n10` is the number of patients that have used the interested drug but did
#' not experience interested event during any of the model periods associated
#' with the drug, `n01` is the number of patients that did not use the
#' interested drug but experienced interested event, and `n00` is the number of
#' patients that were not exposed to the interested drug and did not report
#' interested condition.
#'
#' @section Contingency table:
#' |                   | ADR of interest | Other ADRs |  Total  |
#' |-------------------|-----------------|------------|---------|
#' | Drug of interest  |      a=n11      |    b=n10   | a+b=n1. |
#' | Other drugs       |      c=n01      |    d=n00   | c+d=n0. |
#' | Total             |     a+c=n.1     |   b+d=n.0  |a+b+c+d=n|
#' @return A [data.table][data.table::data.table] with columns of estimated
#' value and it's confidence interval (`ci_low` and `ci_high`). Estimated column
#' are as follows:
#' - `phv_ror`: reporting odds ratio (`ror`).
#' - `phv_prr`: proportional reporting ratio (`prr`). Signal defined as a `prr`
#'   of at least 2, chi-squared with Yates's correction of at least 4 and `a >=
#'   3`. An equivalent alternative to chi-squared is to calculate a confidence
#'   interval around the `prr`.
#' - `phv_bcpnn_norm`: information component (`ic`).
#' - `phv_bcpnn_mcmc`: information component (`ic`).
#' - `phv_obsexp_shrink`: observed to expected ratio (`oe_ratio`).
#' - `phv_ebgm`: Empirical Bayes Geometric Mean (`ebgm`).
#' @examples
#' phv_signal(122, 1320, 381, 31341)
#' phv_signal(122, 1320, 381, 31341, "ror")
#' phv_ror(122, 1320, 381, 31341)
#' phv_signal(122, 1320, 381, 31341, "prr")
#' phv_prr(122, 1320, 381, 31341)
#' phv_signal(122, 1320, 381, 31341, "chisq")
#' phv_chisq(122, 1320, 381, 31341)
#' phv_signal(122, 1320, 381, 31341, "bcpnn_norm")
#' phv_bcpnn_norm(122, 1320, 381, 31341)
#' phv_signal(122, 1320, 381, 31341, "bcpnn_mcmc")
#' phv_bcpnn_mcmc(122, 1320, 381, 31341)
#' phv_signal(122, 1320, 381, 31341, "obsexp_shrink")
#' phv_obsexp_shrink(122, 1320, 381, 31341)
#' phv_signal(122, 1320, 381, 31341, "fisher")
#' phv_fisher(122, 1320, 381, 31341)
#' phv_signal(122, 1320, 381, 31341, "ebgm")
#' phv_ebgm(122, 1320, 381, 31341)
#' @references
#' - David Olaleye, SAS Institute Inc. (2019), Real-World Evidence and
#'   Population Health Analytics: Intersection and Application,
#'   <https://support.sas.com/resources/papers/proceedings19/3361-2019.pdf>
#' @export
#' @name phv_signal
#' @importFrom BiocParallel SerialParam
phv_signal <- function(a, b, c, d, methods = NULL, alpha = 0.05, correct = TRUE, n_mcmc = 1e5L, alpha1 = 0.5, alpha2 = 0.5, theta_init = NULL, squashing = TRUE, BPPARAM = SerialParam()) {
    allowed_methods <- c(
        "ror", "prr", "chisq", "bcpnn_norm", "bcpnn_mcmc",
        "obsexp_shrink", "fisher", "ebgm"
    )
    assert_inclusive(methods, allowed_methods, allow_null = TRUE)
    methods <- unique(methods %||% allowed_methods)
    out <- data.table(expected = (a + b) / (a + b + c + d) * (a + c))
    args <- list(a = a, b = b, c = c, d = d, alpha = alpha)
    for (method in methods) {
        phv_fn <- sprintf("phv_%s", method)
        cli::cli_alert_info("Running {.fn {phv_fn}}")
        signal_out <- switch(method,
            bcpnn_mcmc = do.call(
                phv_fn,
                c(args, list(n_mcmc = n_mcmc, BPPARAM = BPPARAM))
            ),
            obsexp_shrink = do.call(
                phv_fn,
                c(args, list(
                    alpha1 = alpha1, alpha2 = alpha2, n_mcmc = n_mcmc,
                    BPPARAM = BPPARAM
                ))
            ),
            chisq = do.call(phv_fn, c(
                args[c("a", "b", "c", "d")],
                list(correct = correct, BPPARAM = BPPARAM)
            )),
            ebgm = do.call(phv_fn, c(
                args[c("a", "b", "c", "d")],
                list(theta_init = theta_init, squashing = squashing)
            )),
            fisher = do.call(phv_fn, c(args, list(BPPARAM = BPPARAM))),
            do.call(phv_fn, args)
        )
        added_names <- names(signal_out)
        added_names <- data.table::fcase(
            added_names %chin% c("ci_low", "ci_high"),
            paste(added_names[1L], added_names, sep = "_"),
            added_names == "pvalue",
            paste(method, added_names, sep = "_"),
            rep_len(TRUE, length(added_names)), added_names
        )
        if (startsWith(method, "bcpnn_")) {
            added_names <- paste(method, added_names, sep = "_")
        }
        out[, (added_names) := signal_out]
    }
    out[]
}

#' @export
#' @rdname phv_signal
phv_ror <- function(a, b, c, d, alpha = 0.05) {
    # https://github.com/tystan/pharmsignal/blob/master/R/ror_signal.R
    assert_phv_table(a, b, c, d)

    # ugly way to do ad/bc but this way avoids integer overflow
    ror <- (a / b) * (d / c)
    log_ror <- log(ror)
    var_log_ror <- sqrt(1L / a + 1L / b + 1L / c + 1L / d)
    half <- alpha / 2L
    ci_low <- exp(stats::qnorm(half, log_ror, var_log_ror))
    ci_high <- exp(stats::qnorm(1L - half, log_ror, var_log_ror))

    data.table(ror = ror, ci_low = ci_low, ci_high = ci_high)
}

#' @export
#' @rdname phv_signal
#' @references
#' - Evans, S.J.W., Waller, P.C. and Davis, S. (2001), Use of proportional
#'   reporting ratios (PRRs) for signal generation from spontaneous adverse drug
#'   reaction reports. Pharmacoepidem. Drug Safe., 10: 483-486.
#'   <https://doi.org/10.1002/pds.677>
phv_prr <- function(a, b, c, d, alpha = 0.05) {
    # https://github.com/tystan/pharmsignal/blob/master/R/prr_signal.R
    assert_phv_table(a, b, c, d)

    # run prr analysis
    prr <- (a / (a + b)) / (c / (c + d))
    log_prr <- log(prr)
    sd_log_prr <- sqrt(1L / a - 1L / (a + b) + 1L / c - 1L / (c + d))
    half <- alpha / 2L
    ci_low <- exp(stats::qnorm(half, log_prr, sd_log_prr))
    ci_high <- exp(stats::qnorm(1 - half, log_prr, sd_log_prr))
    data.table(prr = prr, ci_low = ci_low, ci_high = ci_high)
}

#' @param correct A bool indicating whether to apply Yates's continuity
#' correction when computing the chi-squared statistic.
#' @export
#' @rdname phv_signal
phv_chisq <- function(a, b, c, d, correct = TRUE, BPPARAM = SerialParam()) {
    assert_phv_table(a, b, c, d)
    # Zero-length input has no test to run. Returning the full schema keeps the
    # result bindable with a non-empty call (e.g. `phv_signal()` on an empty
    # group), instead of `setnames()` failing on a zero-column data.table.
    if (!length(a)) {
        return(data.table(chisq = numeric(), pvalue = numeric()))
    }
    out <- BiocParallel::bpmapply(
        function(n11, n10, n01, n00) {
            out <- stats::chisq.test(
                matrix(c(n11, n10, n01, n00), nrow = 2L),
                correct = correct
            )
            c(out$statistic, out$p.value)
        },
        n11 = a, n10 = b, n01 = c, n00 = d,
        SIMPLIFY = FALSE, USE.NAMES = FALSE,
        BPPARAM = BPPARAM
    )
    out <- data.table::transpose(out)
    data.table::setDT(out)
    data.table::setnames(out, c("chisq", "pvalue"))
    out[]
}

#' @export
#' @rdname phv_signal
phv_fisher <- function(a, b, c, d, alpha = 0.05, BPPARAM = SerialParam()) {
    assert_phv_table(a, b, c, d)
    # See `phv_chisq()`: keep the schema so a zero-length call stays bindable.
    if (!length(a)) {
        return(data.table(
            odds_ratio = numeric(), ci_low = numeric(),
            ci_high = numeric(), pvalue = numeric()
        ))
    }
    out <- BiocParallel::bpmapply(
        function(n11, n10, n01, n00) {
            out <- stats::fisher.test(
                matrix(c(n11, n10, n01, n00), nrow = 2L),
                conf.int = TRUE, conf.level = 1 - alpha
            )
            c(out$estimate, out$conf.int, out$p.value)
        },
        n11 = a, n10 = b, n01 = c, n00 = d,
        SIMPLIFY = FALSE, USE.NAMES = FALSE,
        BPPARAM = BPPARAM
    )
    out <- data.table::as.data.table(do.call("rbind", out))
    data.table::setnames(out, c("odds_ratio", "ci_low", "ci_high", "pvalue"))
    out[]
}

#' @export
#' @rdname phv_signal
phv_bcpnn_norm <- function(a, b, c, d, alpha = 0.05) {
    assert_phv_table(a, b, c, d)
    # https://github.com/tystan/pharmsignal/blob/master/R/bcpnn_norm_signal.R
    # run bcpnn analysis
    n1. <- a + b
    n.1 <- a + c
    n <- a + b + c + d

    p1 <- 1L + n1.
    p2 <- 1L + n - n1.
    q1 <- 1L + n.1
    q2 <- 1L + n - n.1
    r1 <- 1L + a
    # fix q1 * p1 : NAs produced by integer overflow
    r2b <- n - a - 1L + (2L + n)^2L / q1 / p1

    ic <- (digamma(r1) - digamma(r1 + r2b) -
        (digamma(p1) - digamma(p1 + p2) + digamma(q1) - digamma(q1 + q2))) /
        log(2L)
    var_ic <- (trigamma(r1) - trigamma(r1 + r2b) +
        (trigamma(p1) - trigamma(p1 + p2) + trigamma(q1) -
            trigamma(q1 + q2))) / log(2L)^2L
    half <- alpha / 2L
    ci_low <- stats::qnorm(half, ic, sqrt(var_ic))
    ci_high <- stats::qnorm(1 - half, ic, sqrt(var_ic))
    data.table(ic = ic, ci_low = ci_low, ci_high = ci_high)
}

#' @param n_mcmc Number of MCMC simulations per `(a,b,c,d)`-tuple to calculate
#' confidence intervals.
#' @export
#' @rdname phv_signal
phv_bcpnn_mcmc <- function(a, b, c, d, alpha = 0.05, n_mcmc = 1e5L, BPPARAM = SerialParam()) {
    assert_phv_table(a, b, c, d)
    # run bcpnn analysis
    # NOTE: we could speed up the code by combining some of the expressions
    # here. We decided to keep it like this for readability's sake.
    # https://github.com/bips-hb/pvm/blob/master/R/BCPNN.R
    # determine the marginals and the total
    n1. <- a + b
    n.1 <- a + c
    n <- a + b + c + d

    # determine the q-values (see eq. (5) in the paper by Noren et al.)
    q1. <- (n1. + 0.5) / (n + 1)
    q.1 <- (n.1 + 0.5) / (n + 1)
    q.0 <- (n - n.1 + 0.5) / (n + 1)
    q0. <- (n - n1. + 0.5) / (n + 1)

    # determine the prior parameters (see eq. (3) & (4) in the paper by Noren et al.)
    alpha.. <- 0.5 / (q1. * q.1)
    alpha11 <- q1. * q.1 * alpha..
    alpha10 <- q1. * q.0 * alpha..
    alpha01 <- q0. * q.1 * alpha..
    alpha00 <- q0. * q.0 * alpha..

    # determine the gammas on the basis of the alphas and the observed counts
    # see eq. (2) in the paper by Noren et al.
    gamma11 <- alpha11 + a
    gamma10 <- alpha10 + b
    gamma01 <- alpha01 + c
    gamma00 <- alpha00 + d

    ic <- log2(
        gamma11 * (gamma11 + gamma10 + gamma01 + gamma00) /
            ((gamma11 + gamma10) * (gamma11 + gamma01))
    )
    if (!length(a)) {
        # `ic` is computed above and is already zero-length; the MCMC loop has
        # nothing to sample, so return the schema directly (see `phv_chisq()`).
        return(data.table(ic = ic, ci_low = numeric(), ci_high = numeric()))
    }
    out <- BiocParallel::bpmapply(
        function(g11, g10, g01, g00) {
            # sample from the posterior distribution
            p <- MCMCpack::rdirichlet(n_mcmc, c(g11, g10, g01, g00))

            # compute the corresponding ICs for the samples from the posterior distribution
            p11 <- p[, 1]
            p1. <- p11 + p[, 2]
            p.1 <- p11 + p[, 3]

            ic_monte <- log2(p11 / (p1. * p.1))
            # posterior distribution quantiles:
            # (0.025, 0.975) for alpha = 0.05
            stats::quantile(ic_monte, c(alpha / 2, 1 - alpha / 2))
        },
        g11 = gamma11, g10 = gamma10, g01 = gamma01, g00 = gamma00,
        SIMPLIFY = FALSE, USE.NAMES = FALSE,
        BPPARAM = BPPARAM
    )
    out <- data.table::as.data.table(do.call("rbind", out))
    data.table::setnames(out, c("ci_low", "ci_high"))
    out[, ic := ic]
    data.table::setcolorder(out, "ic")
    out[]
}

#' @param alpha1 Numerator shrinkage parameter `>=0`, default `0.5`.
#' @param alpha2 Denominator shrinkage parameter `>=0`, default `0.5`.
#' @section phv_obsexp_shrink:
#'
#' The observed to expected (OE) ratio with approximate confidence intervals are
#' constructed on the log2 scale as outlined in Norén et al. (2013).
#'
#' Expected value was estimated by `(a + b) / (a + b + c + d) * (a + c)`.
#'
#' The OE ratio with shrinkage estimates is calculated as `(O + alpha1) / (E +
#' alpha2)`.
#'
#' If `(O + alpha1) < 1`, then the exact uncertainty limits should be used. That
#' is the confidence intervals as implemented in `phv_bcpnn_mcmc` (Norén et al.,
#' 2013).
#'
#' `log2(OE)` approximates the Bayesian confidence propagation neural
#' network information component (IC) with reasonable accuracy when `alpha1
#' = alpha2 = 0.5` (Norén et al., 2013).
#' @references
#' - Norén GN, Hopstadius J, Bate A. Shrinkage observed-to-expected ratios for
#'   robust and transparent large-scale pattern discovery.  Statistical methods
#'   in medical research. 2013 Feb;22(1):57-69.
#' @export
#' @rdname phv_signal
phv_obsexp_shrink <- function(a, b, c, d, alpha = 0.05, alpha1 = 0.5, alpha2 = 0.5, n_mcmc = 1e5L, BPPARAM = SerialParam()) {
    # https://github.com/tystan/pharmsignal/blob/master/R/obsexp_shrink_signal.R
    assert_phv_table(a, b, c, d)
    # run bcpnn analysis
    n <- a + b + c + d
    n1. <- a + b
    n.1 <- a + c
    expected <- (n1. / n) * n.1

    oe_ratio <- log2((a + alpha1) / (expected + alpha2))
    ci_low <- oe_ratio - 3.3 * (a + alpha1)^(-1L / 2L) -
        2L * (a + alpha1)^(-3L / 2L)
    ci_high <- oe_ratio + 2.4 * (a + alpha1)^(-1L / 2L) -
        0.5 * (a + alpha1)^(-3L / 2L)

    need_exact_lims <- (a + alpha1) < 1L
    if (any(need_exact_lims)) {
        cli::cli_warn(c(
            `!` = "{sum(need_exact_lims)} case{?s} of (O + alpha1) < 1",
            i = "will use the information components from the BCPNN (using MCMC estimation) to estimate the confidence interval" # nolint
        ))

        # use default mcmc draws
        a_star <- a[need_exact_lims]
        b_star <- b[need_exact_lims]
        c_star <- c[need_exact_lims]
        d_star <- d[need_exact_lims]
        ic <- phv_bcpnn_mcmc(a_star, b_star, c_star, d_star,
            alpha = alpha, n_mcmc = n_mcmc, BPPARAM = BPPARAM
        )
        # replace the CIs for those with (O + alpha1) < 1
        ci_low[need_exact_lims] <- ic$ci_low
        ci_high[need_exact_lims] <- ic$ci_high
    }
    data.table(
        oe_ratio = oe_ratio,
        ci_low = ci_low, ci_high = ci_high
    )
}

#' @param theta_init A data frame of initial hyperparameter guesses with columns
#' ordered as:` alpha1, beta1, alpha2, beta2, P`. See [openEBGM::autoHyper]
#' @param squashing A bool, whether do automated data squashing. If any zeros
#' found in `a`, will always be `TRUE`.
#' @section phv_ebgm:
#' An implementation of the Gamma-Poisson Shrinker (GPS) model for identifying
#' unexpected counts in large contingency tables using an empirical Bayes
#' approach. The Empirical Bayes Geometric Mean (EBGM) and quantile scores are
#' obtained from the GPS model estimates.
#' The GPS was proposed by DuMouchel as a signal detection tool for large
#' frequency tables with both observed (O) and expected (E) counts for each
#' drug-outcome pair. It assumes the observed count of any drug-outcome pair
#' follows the Poisson distribution.
#'
#' For each drug-outcome pair, the primary parameter of interest was the risk
#' ratio. Rather than using the observed over expected (O/E), GPS uses the
#' empirical Bayesian geometric mean (EBGM) posterior distribution of the risk
#' ratio and the surrounding confidence interval for each drug-outcome pair to
#' identify statistical signals of excess risk. To prevent spurious false
#' positives due to implausibly high risk ratios, GPS implements a Bayesian
#' framework that “shrinks” O/E estimates towards a value which is close to the
#' average O/E values for all drug-event pairs at each level of granularity.
#' @references
#' - <https://journal.r-project.org/archive/2017/RJ-2017-063/RJ-2017-063.pdf>
#' @export
#' @rdname phv_signal
phv_ebgm <- function(a, b, c, d, alpha = 0.05, theta_init = NULL, squashing = TRUE) {
    assert_phv_table(a, b, c, d)
    assert_bool(squashing)
    assert_s3_class(
        theta_init,
        is_class = is.data.frame,
        what = "a data frame",
        allow_null = TRUE
    )
    proc <- data.table(N = a, E = (a + b) / (a + b + c + d) * (a + c))
    # For the purpose of reducing computational burden, zero counts should not
    # typically be included for hyperparameter estimation; however, they may
    # help when convergence issues are encountered. If included, the points
    # should be squashed (discussed in later sections) in order to reduce the
    # computation time. These zero counts should not typically be used for EB
    # scores since they are meaningless for studying larger-than-expected counts
    # (which is usually the goal).
    # In some situations, 0 will always cause error, so we just remove zeros
    zeros <- any(a == 0L)
    if (zeros && !squashing) {
        cli::cli_warn("Using {.code squashing = TRUE} since zeros exist")
        squashing <- TRUE
    }
    if (squashing) {
        squashed <- openEBGM::autoSquash(proc)
    } else {
        squashed <- proc
    }
    theta_init <- theta_init %||% data.frame(
        alpha1 = c(0.2, 0.1, 0.3),
        beta1 = c(0.1, 0.2, 0.5),
        alpha2 = c(2, 10, 6),
        beta2 = c(4, 10, 6),
        p = c(1 / 3, 0.2, 0.5)
    )
    # The N_star argument allows the user to choose the smallest value of N used
    # for hyperparameter estimation. The user must be careful to match N_star
    # with the actual minimum count in a. If the user wishes to use a larger
    # value for N_star, the vectors supplied to arguments `ni`, `ei`, and `wi`
    # must be filtered. Here, we are using all counts except zeroes, so no
    # filtering is needed. In general, N_star = 1 should be used whenever
    # practical.
    theta_hats <- tryCatch(
        openEBGM::autoHyper(
            data = squashed,
            theta_init = theta_init,
            squashed = squashing,
            zeroes = zeros,
            N_star = if (zeros) NULL else min(squashed$N)
        )$estimates,
        error = function(cnd) {
            cli::cli_warn(c(
                x = "Cannot estimate hyperparameter using {.pkg openEBGM}",
                i = "Using {.fn nlminb} directly to estimate hyperparameter"
            ))
            # stats::nlminb(
            #     start = c(0.2, 0.1, 2, 4, 1 / 3),
            #     objective = objective,
            #     ni = squashed$N, ei = squashed$E,
            #     wi = squashed$weight, lower = 0,
            #     upper = Inf,
            #     N_star = if (zeros) NULL else min(squashed$N)
            # )$par
            # https://github.com/bips-hb/pvm/blob/master/R/fitPriorParametersGPS.R
            init <- c(0.2, 0.1, 2, 4, 1 / 3)
            if (zeros) {
                out <- stats::nlminb(
                    start = init, objective = openEBGM::negLLzero,
                    N = proc$N, E = proc$E,
                    lower = 0, upper = c(Inf, Inf, Inf, Inf, 1.0),
                    control = list(iter.max = 500L)
                )
            } else {
                out <- stats::nlminb(
                    start = init, objective = openEBGM::negLL,
                    N = proc$N, E = proc$E,
                    lower = 0, upper = c(Inf, Inf, Inf, Inf, 1.0),
                    N_star = if (zeros) NULL else min(proc$N),
                    control = list(iter.max = 500L)
                )
            }
            out$par
        }
    )
    qn <- phv_ebgm_qn(theta_hats, N = proc$N, E = proc$E)
    ebgm <- phv_ebgm_score(theta_hats, N = proc$N, E = proc$E, qn = qn)
    half <- alpha / 2L
    ci_low <- phv_ebgm_quant_bisect(half,
        theta_hat = theta_hats,
        N = proc$N, E = proc$E, qn = qn
    )
    ci_high <- phv_ebgm_quant_bisect(
        1 - half,
        theta_hat = theta_hats,
        N = proc$N, E = proc$E, qn = qn
    )
    data.table(ebgm = ebgm, ci_low = ci_low, ci_high = ci_high)
}

#############################################################
assert_phv_table <- function(a, b, c, d, call = rlang::caller_env()) {
    lst <- list(a, b, c, d)
    right_cls <- vapply(lst, function(x) {
        !anyNA(x) && is.numeric(x) && all(x >= 0L)
    }, logical(1L))
    if (!all(right_cls)) {
        fail_cls <- c("a", "b", "c", "d")[!right_cls] # nolint
        cli::cli_abort(
            "{.arg {fail_cls}} must be non-negative numeric vector",
            call = call
        )
    }
    lens <- lengths(lst)
    if (!all(lens == lens[1L])) {
        cli::cli_abort(
            c("all {.arg a}, {.arg b}, {.arg c}, {.arg d} must be numeric vector and all have the same length",
                i = "length of {.arg a}, {.arg b}, {.arg c}, {.arg d}: {lens}"
            ),
            call = call
        )
    }
}
