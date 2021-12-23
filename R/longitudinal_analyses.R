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
# k = 1  -> All observation
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
dta <- read_xlsx("data-raw/Evolution-november_AI_updated.xlsx",
                 na = c("", "NA", "-"))
names(dta) <- trimws(gsub("\\(.+\\)|\\r\\n", "", names(dta)))

# Subset
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
Y <- c("TV", "T", "LH", "FSH", "AMH", "INB")

# Analyses
R <- mclapply(setNames(Y, Y), function(y) {
  sdta <- dta %>%
    select(ID, DX, Age, one_of(y)) %>%
    drop_na()
  fml <- as.formula(paste(y, "~ DX * Age + (1 | ID)"))
  ctrl <- lmerControl(optimizer ="Nelder_Mead")
  DX_list <- c("CDGP", "Partial CHH", "Complete CHH")
  DX_list <- setNames(DX_list, sub(" ", "_", DX_list))
  ref_ages <- c(0, seq(14, c(20, 16.5)[k], .5)) %>% setNames(., .)
  fits <- lapply(ref_ages, function(a) {
    lapply(DX_list, function(dx) {
      sdta <- sdta %>%
        mutate(
          DX = relevel(DX, dx),
          Age = Age - a
        )
      fit <- do.call("lmer", list(formula = fml, data = quote(sdta),
                                  control = quote(ctrl)))
      tbl <- tidy(fit, conf.int = TRUE) %>%
        select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
        mutate(
          term = sub("^DX", "", term),
          term = sub(":", " x ", term),
          term = sub("sd__", "SD ", term),
        )
      if (a != 0) {
        tbl <- tbl %>%
          mutate(term = sub("Age", paste0("AgeC", a), term))
      }
      names(tbl)[1] <- paste0("term (ref: ", dx, ")")
      list(fit = fit, tbl = tbl)
    })
  })
  tbls <- lapply(fits, function(z) lapply(z, function(w) w$tbl))
  fits <- lapply(fits, function(z) lapply(z, function(w) w$fit))
  figs <- list()
  dgp <- plot_model(fits[[1]][[1]], type = "diag")
  dgp[[2]] <- dgp[[2]]$ID
  figs$diag <- ggarrange(plotlist = dgp, nrow = 2, ncol = 2) %>%
    annotate_figure(top = text_grob("Diagnostic plots", 
                    face = "bold", size = 16)) %>%
    suppressMessages()
  figs$pred1 <- augment(fits[[1]][[1]]) %>%
    group_by(ID) %>%
    filter(n() > 1) %>%
    ggplot(aes(Age, !!sym(y))) +
    geom_point() +
    geom_line(aes(y = .fitted, colour = DX)) +
    facet_wrap(~ID) +
    labs(title = "Individual predictions") +
    theme(legend.position = "bottom", legend.title = element_blank())
  ndta <- do.call(rbind, lapply(names(fits[[1]]), function(dx) {
    ggpredict(fits[[1]][[dx]], "Age") %>%
      as_tibble() %>%
      select(x, predicted, conf.low, conf.high) %>%
      mutate(dx = dx)
  })) %>%
    mutate(dx = factor(dx, names(fits[[1]])))
  figs$pred2 <- ggplot(ndta, aes(x = x, y = predicted)) +
    geom_line(aes(colour = dx)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = dx), alpha = 0.3,
                show.legend = FALSE) +
    labs(x = "Age", y = y, title = "Fixed effects predictions") +
    theme(legend.position = "bottom", legend.title = element_blank())
  list(fits = fits, tbls = tbls, figs = figs)
})

# Any singular fit ?
b <- lapply(R, function(r) lapply(r$fits, function(z) lapply(z, isSingular)))
if (any(unlist(b))) warning("Singular fit")
rm(b)

# Export results - Tables
Z <- lapply(R, function(r) lapply(r$tbls, function(z) {
  tbl <- lapply(z, cbind, NA) %>%
    append(list(.name_repair = "minimal")) %>%
    do.call(bind_cols, .)
  names(tbl)[names(tbl) == "NA"] <- ""
  return(tbl)
}))
for (y in names(Z)) {
  f <- paste(y, "regression_coefficients.xlsx", sep = "_")
  f <- file.path(outdir, f)
  write_xlsx(Z[[y]], f)
}
rm(Z, f, y)

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
