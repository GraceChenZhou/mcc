
# mcc: Mean Cumulative Count for Left-Truncated Competing Risks Data

[![R-CMD-check](https://github.com/GraceChenZhou/mcc/actions/workflows/r.yml/badge.svg)](https://github.com/GraceChenZhou/mcc/actions/workflows/r.yml)

The **mcc** R package provides a robust and efficient tool to estimate
the total burden of recurrent events within a population using the Mean
Cumulative Count (MCC) method.

It is specifically designed to handle complex longitudinal survival
data, accurately applying time-dependent censoring and truncation
weights for left-truncated data in the presence of competing risks. The
methodology implements an optimized [Geskus
(2011)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1541-0420.2010.01420.x)
weighting approach, making it highly suitable for complex, long-term
cohort studies.

# Installation

### From GitHub

You can install the development version of the package directly from
GitHub:

``` r
# Install devtools if you haven't already
# install.packages("devtools")

devtools::install_github("GraceChenZhou/mcc")
```

### From Local Source Tarball

If you received the compiled source package (`.tar.gz`), install it
locally using:

``` r
install.packages("path/to/mcc_0.x.0.tar.gz", repos = NULL, type = "source")
```

# Dependencies

To maintain mathematical accuracy when computing weights for
left‑truncated data, this package requires **R version 4.4.0 or later**
and **survival version 3.7-0 or later**.

# Basic Usage

The primary function is `mcc()`. It evaluates unstratified cohorts or
handles group stratification natively while providing convenient S3
methods for reporting and extracting your data. The package also exports
`simulate_subject()` for generating one synthetic subject and
`simulate_mcc_data()` for reproducible multi-group examples.

``` r
library(mcc)

## 1. Simulate recurrent-event data for 300 subjects (150 per group)
# simulate_mcc_data() is included in the mcc package.
# It generates zero or more recurrent events per subject and one terminal
# censoring/competing-event row, with optional delayed entry.
mydata <- simulate_mcc_data(
  n_per_group = 150,
  event_rates = c(A = 0.22, B = 0.38),
  max_followup = 10,
  entry_range = c(0, 2),
  death_rate = 0.03,
  censor_rate = 0.05,
  seed = 2026
)

# 0 = censored, 1 = recurrent event, 2 = competing risk
# Each subject has one final status 0/2 record; subjects can have multiple
# status 1 records. Tstart is constant for each subject.
table(mydata$group[!duplicated(mydata$id)])
#> 
#>   A   B 
#> 150 150
length(unique(mydata$id))  # 300 subjects
#> [1] 300
head(mydata)
#>   id      time status    Tstart group
#> 1  1  1.593865      0 1.3973469     A
#> 2  2 10.000000      0 1.1130610     A
#> 3  3  7.510427      1 0.2802799     A
#> 4  3  8.291473      1 0.2802799     A
#> 5  3  8.769644      1 0.2802799     A
#> 6  3 10.000000      0 0.2802799     A

## 2. Run Stratified MCC with Bootstrap Confidence Intervals
# Set ci = TRUE to compute bootstrap intervals
# For demonstration purposes, only 100 bootstrap iterations (niter = 100) are used to reduce computation time. For actual statistical analyses, a larger number of bootstrap iterations (e.g., niter = 1000 or more) is recommended to obtain more stable confidence interval estimates.

result <- mcc(id = mydata$id, 
              time = mydata$time, 
              status = mydata$status, 
              Tstart = mydata$tstart,
              group = mydata$group,
              ci = TRUE, 
              niter = 100, 
              seed = 123)
                 
## 3. Extract numeric estimates
estimates_df <- as.data.frame(result)
head(estimates_df)
#> # A tibble: 6 × 5
#>   group  time     MCC     lci    uci
#>   <chr> <dbl>   <dbl>   <dbl>  <dbl>
#> 1 A     0.178 0.00667 0       0.0233
#> 2 A     0.193 0.0133  0       0.0333
#> 3 A     0.422 0.02    0       0.0433
#> 4 A     0.480 0.02    0       0.0433
#> 5 A     0.497 0.0267  0.00667 0.0501
#> 6 A     0.632 0.0334  0.00667 0.0601
```

To display a concise numerical summary at selected time points:

``` r
summary(result, times = c(2, 4, 6, 8))
#> Mean Cumulative Count (MCC) Summary
#> =============================================
#>  time       MCC_A (95% CI)       MCC_B (95% CI)
#>     2 0.258 (0.180, 0.342) 0.328 (0.235, 0.428)
#>     4 0.639 (0.521, 0.807) 1.077 (0.917, 1.250)
#>     6 0.994 (0.827, 1.159) 1.849 (1.685, 2.029)
#>     8 1.391 (1.216, 1.568) 2.376 (2.127, 2.639)
```

# Plotting MCC Curves and Number-at-Risk Tables

The package provides two plotting interfaces built with **ggplot2** and
**patchwork**. Group-specific MCC estimates are calculated directly by
`mcc(group = ...)`; there is no need to fit groups separately and
combine their results:

- `plot(result)` is the convenient S3 method for an object returned by
  `mcc()`. It displays the combined MCC curve and number-at-risk table.
- `plot_mcc_with_risk_table(result)` provides the same display and
  returns the individual plot objects and underlying plotting data for
  further customization.

The examples below use the simulated 300-subject `result` from **Basic
Usage**, fitted with `ci = TRUE` and stratified by groups A and B. The
original analysis data are stored in the fitted `mcc` object, so you do
not need to pass them again.

## 1. Basic plot

``` r
# Display MCC step curves with bootstrap confidence ribbons and a risk table
plot(result, risk_times = c(2, 4, 6, 8, 10))
```

![](man/figures/README-basic-plot-1.png)<!-- -->

Confidence ribbons require confidence limits from an `mcc()` fit with
`ci = TRUE`. For a fit with `ci = FALSE`, the function can still display
MCC curves and a risk table without ribbons. The example above uses
`niter = 100` for illustration.

## 2. Highlight a control group and customize colors

``` r
plot(
  result,
  control_group = "A",       # Control is dashed by default
  color_values = c(A = "#2166AC", B = "#B2182B"),
  risk_times = c(2, 4, 6, 8, 10),
  title = "Mean cumulative count by group",
  xlab = "Follow-up time",
  ylab = "Mean cumulative count",
  legend.position = "right"
)
```

![](man/figures/README-control-colors-1.png)<!-- -->

`color_values` and `line_values` are **named vectors indexed by the
original group names**, not by the displayed legend labels. By default,
the control group is represented by a dashed line, but this can be
overridden by specifying a custom line type through line_values.
Additionally, users can customize the risk-table time points using the
risk_times argument.

``` r
plot(
  result,
  control_group = "A",
  line_values = c(A = "dotdash", B = "solid"),
  palette = "Dark2",
  risk_times = c(2, 6, 10)
)
```

![](man/figures/README-control-line-types-1.png)<!-- -->

The risk table counts **distinct subjects**, not recurrent-event
records. For left-truncated analyses, a subject is counted at time `t`
when their entry time is strictly before `t` and their observed exit
time is at or after `t` (`Tstart < t <= time`). Without left truncation,
the criterion is `time >= t`. These counts describe the observed risk
sets, not inverse-probability-weighted effective sample sizes.
`risk_times` sets the displayed risk-table time points and the
corresponding x-axis ticks. When omitted, the function chooses default
time points based on the observed MCC time range.

## 3. Choose the analysis time scale

The `time` and `Tstart` inputs must use the **same time scale**. For an
age-scale analysis, `time` is age at each event or competing risk or
censoring and `Tstart` is age at study entry. The delayed-entry example
below uses the simulated dataset on its generic follow-up time scale.

``` r

fit_delayed <- mcc(
  id = mydata$id, time = mydata$time, status = mydata$status,
  Tstart = mydata$tstart, group = mydata$group
)
plot(fit_delayed, conf_int = FALSE, xlab = "Follow-up time")
```

![](man/figures/README-delayed-entry-1.png)<!-- -->

## 4. Customize themes, labels, and layout

``` r
plot(
  result,
  conf_alpha = 0.20,
  ggtheme = ggplot2::theme_bw(),
  tables.theme = ggplot2::theme_classic(),
  title = "Recurrent-event burden over follow-up",
  title_hjust = 0.5,
  title_size = 15,
  axis_label_size = 12,
  axis_text_size = 10,
  legend.title = "Study group",
  legend.position = "right",
  risk_table_title = "Number at risk",
  risk_table_text_size = 3.5,
  plot_heights = c(4, 1)
)
```

![](man/figures/README-custom-layout-1.png)<!-- -->

## 5. Retrieve and reuse the plot components

Unlike `plot(result)`, the standalone function
`plot_mcc_with_risk_table` returns a list of plot objects and data:

``` r
p <- plot_mcc_with_risk_table(
  result,
  control_group = "A",
  risk_times = c(2, 4, 6, 8, 10)
)

names(p)
#> [1] "plot"       "p_mcc"      "p_risk"     "plot_data"  "risk_table"
#> [6] "risk_times"
head(p$risk_table)
#>   group time n_risk plot_group legend_group
#> 1     A    2    140          A  A (Control)
#> 2     B    2    131          B            B
#> 3     A    4    120          A  A (Control)
#> 4     B    4    115          B            B
#> 5     A    6    105          A  A (Control)
#> 6     B    6     90          B            B

# Show the combined figure
p$plot
```

![](man/figures/README-plot-components-1.png)<!-- -->

``` r

# Save the combined figure
# ggplot2::ggsave("mcc_plot.png", plot = p$plot, width = 8, height = 6, dpi = 300)
```

## 6. Plot an unstratified analysis

``` r
result_overall <- mcc(
  id = mydata$id,
  time = mydata$time,
  status = mydata$status,
  Tstart = mydata$tstart,
  ci = FALSE
)

plot(result_overall, title = "Overall MCC")
```

![](man/figures/README-unstratified-plot-1.png)<!-- -->

An unstratified fit displays one MCC curve and its corresponding
risk-table row. If you already have a data frame of estimates rather
than an `mcc` object, supply the original data explicitly:

``` r
p_df <- plot_mcc_with_risk_table(
  as.data.frame(result),
  data = result$data,
  conf_int = TRUE
)

p_df$plot
```

![](man/figures/README-data-frame-plot-1.png)<!-- -->

For all plotting options, see `?plot_mcc_with_risk_table` and
`?plot.mcc`.

## Simulation assumptions and reproducibility

`simulate_subject()` models recurrent events as a homogeneous Poisson
process after entry, with independent exponential competing-event and
censoring times. A zero hazard disables that terminal mechanism.
`simulate_mcc_data()` assigns unique subject IDs across groups, uses
group-specific recurrent-event rates, and samples entry times from
`entry_range`. These functions create **synthetic demonstration data**,
not validated models of a particular clinical cohort. The `seed`
argument makes the cohort reproducible without permanently changing the
caller’s random-number-generator state.

``` r
# Simulate one individual directly
set.seed(2026)

one_subject <- simulate_subject(
  id = 1,
  group = "A",
  entry_time = 0.5,
  max_followup = 10,
  event_rate = 0.22,
  death_rate = 0.03,
  censor_rate = 0.05
)

one_subject
#>   id    time status Tstart group
#> 1  1 2.76122      0    0.5     A
```

``` r
# Simulate two groups with different recurrent-event rates
mydata <- simulate_mcc_data(
  n_per_group = 100,
  event_rates = c(A = 0.22, B = 0.38),
  max_followup = 10,
  entry_range = c(0, 2),
  seed = 2026
)

# Preview the first five subject
first_5_ids <- unique(mydata$id)[1:5]
mydata_5 <- mydata[mydata$id %in% first_5_ids, ]
print(mydata_5)
#>   id      time status    Tstart group
#> 1  1  7.079862      0 1.3973469     A
#> 2  2 10.000000      0 1.1130610     A
#> 3  3  1.959921      0 0.2802799     A
#> 4  4 10.000000      0 0.5714466     A
#> 5  5  2.239416      1 1.1107380     A
#> 6  5  3.159859      1 1.1107380     A
#> 7  5  7.148878      1 1.1107380     A
#> 8  5  9.700798      2 1.1107380     A

# Number of unique subjects in each group
table(
  mydata$group[
    !duplicated(mydata$id)
  ]
)
#> 
#>   A   B 
#> 100 100

# Distribution of recurrent events, censoring, and competing risks
table(mydata$group, mydata$status)
#>    
#>       0   1   2
#>   A  77 143  23
#>   B  84 266  16
```

For function arguments and return values, see `?simulate_subject` and
`?simulate_mcc_data`.
