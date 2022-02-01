library(broom.mixed)
library(dplyr)
library(ggplot2)
library(lme4)
library(parallel)
library(readxl)
library(tidyr)
library(writexl)

options(mc.cores = detectCores() - 1)

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Puberté)")

# Output directory
outdir <-
  paste0("results/longitudinal_analyses_", format(Sys.Date(), "%Y%m%d"))
outdir <- file.path(outdir, "nonlinear")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Load data
dta <- read_xlsx("data-raw/Evolution-december_2021_AI_updated.xlsx",
                 na = c("", "NA", "-"))
names(dta) <- trimws(gsub("\\(.+\\)|\\r\\n", "", names(dta)))

# Add DX
dta <- dta %>%
  select(-DX) %>%
  left_join(unique(na.omit(dta[c("ID", "DX")])), by = "ID") %>%
  mutate(
    DX = sub(", Kallmann", "", DX),
    DX = factor(DX, c("CDGP", "Partial CHH", "Complete CHH")),
    PCHH = as.numeric(DX == "Partial CHH"),
    CCHH = as.numeric(DX == "Complete CHH")
  )

# Outcomes
Y <- c("TV", "T", "LH", "FSH", "IGF1")

# Nonlinear model
Model <- function(Age, Asym, Asym2, Asym3, xmid, scal, PCHH, CCHH) {
  (Asym + Asym2 * PCHH + Asym3 * CCHH) /
    (1 + exp((xmid - Age) / scal))
}
ModelGradient <- deriv(
  body(Model)[[2]],
  namevec = c("Asym", "Asym2", "Asym3", "xmid", "scal"),
  function.arg = Model
)
start_values <- list(
  TV = c(Asym = 20, Asym2 = -10, Asym3 = -18, xmid = 16, scal = 2),
  T = c(Asym = 20, Asym2 = -15, Asym3 = -18, xmid = 16, scal = 2),
  LH = c(Asym = 5, Asym2 = -2.5, Asym3 = -4, xmid = 16, scal = 2),
  FSH = c(Asym = 5, Asym2 = -2.5, Asym3 = -4, xmid = 16, scal = 2),
  IGF1 = c(Asym = 600, Asym2 = -300, Asym3 = 0, xmid = 16, scal = 2)#,
  #AMH = c(Asym = 500, Asym2 = 0, Asym3 = 0, xmid = 16, scal = 2),
  #INB = c(Asym = 400, Asym2 = -250, Asym3 = -350, xmid = 16, scal = 2)
)

# Analyses
R <- mclapply(setNames(Y, Y), function(y) {
  sdta <- dta %>%
    select(ID, DX, PCHH, CCHH, Age, one_of(y)) %>%
    drop_na()
  fml <- ModelGradient(Age, Asym, Asym2, Asym3, xmid, scal, PCHH, CCHH) ~ 
    Asym | ID
  fml <- update(fml, as.formula(paste(y, "~ . ~ .")))
  fit <- do.call("nlmer", list(formula = fml, data = quote(sdta),
                               start = start_values[[y]]))
  tbl <- tidy(fit, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high) %>%
    mutate(
      term = sub(":", " x ", term),
      term = sub("sd__", "SD ", term)
    )
  figs <- list()
  figs$pred1 <- augment(fit) %>%
    mutate(
      DX = c("CDGP", "Partial CHH", "Complete CHH")[1 + PCHH + 2 * CCHH]
    ) %>%
    group_by(ID) %>%
    filter(n() > 1) %>%
    ggplot(aes(Age, !!sym(y))) +
    geom_point() +
    geom_line(aes(y = .fitted, colour = DX)) +
    facet_wrap(~ID) +
    labs(title = "Individual predictions") +
    theme(legend.position = "bottom", legend.title = element_blank())
  ndta <- expand_grid(
    Age = seq(14, 20, 0.1),
    DX = c("CDGP", "Partial CHH", "Complete CHH")
  ) %>% mutate(
    PCHH = as.numeric(DX == "Partial CHH"),
    CCHH = as.numeric(DX == "Complete CHH")
  ) %>% cbind(as.data.frame(t(setNames(tbl[["estimate"]], tbl[["term"]])))) %>%
  rowwise() %>%
  mutate(predicted = Model(Age, Asym, Asym2, Asym3, xmid, scal, PCHH, CCHH))
  figs$pred2 <- ggplot(ndta, aes(x = Age, y = predicted)) +
    geom_line(aes(colour = DX)) +
    labs(x = "Age", y = y, title = "Fixed effects predictions") +
    theme(legend.position = "bottom", legend.title = element_blank())
  tbl <- tbl %>%
    mutate(
      term = sub("^Asym2", "AsymPartCHH", term),
      term = sub("^Asym3", "AsymCompCHH", term)
    )
  list(fit = fit, tbl = tbl, figs = figs)
})

# Any singular fit ?
b <- sapply(R, function(r) isSingular(r$fit))
if (any(unlist(b))) warning("Singular fit")
rm(b)

# Export results - Tables
for (y in names(R)) {
  f <- paste(y, "regression_coefficients.xlsx", sep = "_")
  f <- file.path(outdir, f)
  write_xlsx(R[[y]]$tbl, f)
}
rm(f, y)

# Export results - Figures
for (y in names(R)) {
  tiff(filename = file.path(outdir, paste0(y, "_individual_predictions.tiff")),
       height = 3600, width = 5400, res = 512, compression = "zip")
  print(R[[y]]$figs$pred1)
  dev.off()
  tiff(filename = file.path(outdir,
                            paste0(y, "_fixed_effects_predictions.tiff")),
       height = 3600, width = 5400, res = 1024, compression = "zip")
  print(R[[y]]$figs$pred2)
  dev.off()
}
rm(y)

# Session Info
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
