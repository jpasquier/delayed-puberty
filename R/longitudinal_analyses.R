library(broom.mixed)
library(dplyr)
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

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Puberté)")

# Output directory
outdir <-
  paste0("results/longitudinal_analyses_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(outdir)) dir.create(outdir)

# Load data
dta <- read_xlsx("data-raw/Evolution-november_AI_updated.xlsx",
                 na = c("", "NA", "-"))
names(dta) <- trimws(gsub("\\(.+\\)|\\r\\n", "", names(dta)))

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
    drop_na() %>%
    mutate(AgeC14 = Age - 14)
  fml <- as.formula(paste(y, "~ DX * AgeC14 + (1 | ID)"))
  ctrl <- lmerControl(optimizer ="Nelder_Mead")
  fit <- do.call("lmer", list(formula = fml, data = quote(sdta),
                              control = quote(ctrl)))
  tbl <- tidy(fit, conf.int = TRUE) %>%
    bind_rows(
      do.call(rbind, lapply(c("Partial CHH", "Complete CHH"), function(u) {
        mutate(sdta, DX = relevel(DX, u)) %>%
          lmer(fml, data = ., control = ctrl) %>%
          tidy(conf.int = TRUE) %>%
          filter(term == "AgeC14") %>%
          mutate(term = paste0(term, " (", u, ")"))
      }))
    ) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
    mutate(term = sub("sd__", "SD ", sub(":", " x ", sub("^DX", "", term))))
  figs <- list()
  dgp <- plot_model(fit, type = "diag")
  dgp[[2]] <- dgp[[2]]$ID
  figs$diag <- ggarrange(plotlist = dgp, nrow = 2, ncol = 2) %>%
    annotate_figure(top = text_grob("Diagnostic plots", 
                    face = "bold", size = 16)) %>%
    suppressMessages()
  figs$pred1 <- augment(fit) %>%
    group_by(ID) %>%
    filter(n() > 1) %>%
    mutate(Age = AgeC14 + 14) %>%
    ggplot(aes(Age, !!sym(y))) +
    geom_point() +
    geom_line(aes(y = .fitted, colour = DX)) +
    facet_wrap(~ID) +
    labs(title = "Individual predictions") +
    theme(legend.position = "bottom", legend.title = element_blank())
  ndta <- expand.grid(AgeC14 = range(sdta$AgeC14), DX = levels(sdta$DX))
  ndta[[y]] <- predict(fit, re.form = NA, newdata = ndta)
  figs$pred2 <- ggplot(ndta, aes(x = AgeC14 + 14, y = !!sym(y), colour = DX)) +
    geom_line() +
    labs(title = "Fixed effects predictions", x = "Age") +
    theme(legend.position = "bottom", legend.title = element_blank())
  list(fit = fit, tbl = tbl, figs = figs)
})
if (any(sapply(R, function(r) isSingular(r$fit)))) warning("Singular fit")

# Export results
write_xlsx(lapply(R, function(r) r$tbl),
           file.path(outdir, "regression_coefficients.xlsx"))
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
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
