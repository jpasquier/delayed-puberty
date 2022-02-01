library(broom.mixed)
library(dplyr)
library(ggeffects)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lmerTest)
library(parallel)
library(readxl)
library(sjPlot)
library(tidyr)
library(writexl)

options(mc.cores = detectCores() - 1)

# Define subset
# k = 1  -> All observations
# k = 2  -> Inculude only observations until 16.5y
k <- 2

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Puberté)")

# Output directory
outdir <-
  paste0("results/longitudinal_analyses_", format(Sys.Date(), "%Y%m%d"))
outdir <- file.path(outdir, c("all_observations", "until_16.5y")[k])
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Load data
dta <- read_xlsx("data-raw/Evolution-december_2021_AI_updated.xlsx",
                 na = c("", "NA", "-"))
names(dta) <- trimws(gsub("\\(.+\\)|\\r\\n", "", names(dta)))

# Subset
dta <- dta[dta$Age >= 14, ]
if (k == 2) {
  dta <- dta[dta$Age <= 16.5, ]
}

# Add DX
dx <- unique(na.omit(dta[c("ID", "DX")]))
dx$DX <- factor(sub(", Kallmann", "", dx$DX),
                c("CDGP", "Partial CHH", "Complete CHH"))
names(dta)[names(dta) == "DX"] <- "DX_complete"
dta <- merge(dta, dx, by = "ID", all.x = TRUE)
rm(dx)

# Outcomes
Y <- c("TV", "T", "LH", "FSH", "IGF1", "AMH", "INB")

# Analyses
options(warn = 1)
R <- mclapply(setNames(Y, Y), function(y) {
  sdta <- dta %>%
    select(ID, DX, Age, one_of(y)) %>%
    drop_na() %>%
    mutate(Age2 = Age^2, Age3 = Age^3)
  #fml <- as.formula(paste(y, "~ DX * poly(Age, 3) + (1 | ID)"))
  fml <- as.formula(paste(y, "~ DX * (Age + I(Age^2) + I(Age^3)) + (1 | ID)"))
  ctrl <- lmerControl(optimizer ="Nelder_Mead")
  fit <- do.call("lmer", list(formula = fml, data = quote(sdta),
                 control = quote(ctrl)))
  tbl <- tidy(fit, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(
      term = sub("^DX", "", term),
      term = sub(":", " x ", term),
      term = sub("sd__", "SD ", term),
    )
  figs <- list()
  dgp <- plot_model(fit, type = "diag")
  dgp[[2]] <- dgp[[2]]$ID
  figs$diag <- ggarrange(plotlist = dgp, nrow = 2, ncol = 2) %>%
    annotate_figure(top = text_grob("Diagnostic plots", 
                    face = "bold", size = 16)) %>%
    suppressMessages()
  figs$pred1 <- augment(fit) %>%
    #mutate(
    #  p = `poly(Age, 3)`,
    #  Age = with(attr(p, "coefs"), alpha[1] + sqrt(norm2)[3] * p[, 1])
    #) %>%
    group_by(ID) %>%
    filter(n() > 1) %>%
    ggplot(aes(Age, !!sym(y))) +
    geom_point() +
    geom_line(aes(y = .fitted, colour = DX)) +
    facet_wrap(~ID) +
    labs(title = "Individual predictions") +
    theme(legend.position = "bottom", legend.title = element_blank())
  ndta <- ggpredict(fit, c("Age [all]", "DX")) %>%
    as_tibble() %>%
    select(x, predicted, conf.low, conf.high, group)
  figs$pred2 <- ggplot(ndta, aes(x = x, y = predicted)) +
    geom_line(aes(colour = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
                alpha = 0.3, show.legend = FALSE) +
    labs(x = "Age", y = y, title = "Fixed effects predictions") +
    theme(legend.position = "bottom", legend.title = element_blank())
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
  tiff(filename = file.path(outdir, paste0(y, "_diagnostic_plots.tiff")),
       height = 3600, width = 5400, res = 384, compression = "zip")
  print(R[[y]]$figs$diag)
  dev.off()
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
