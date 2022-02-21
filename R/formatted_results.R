library(broom.mixed)
library(ggeffects)
library(ggplot2)
library(lme4)
library(lmerTest)
library(pROC)
library(readxl)
library(svglite)
library(writexl)

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Puberté)")

# Output directory
outdir <- paste0("results/formatted_results_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(outdir)) dir.create(outdir)

# Load data
f <- "data-raw/Database_DECEMBER_21_updated_withoutnames_final.xlsx"
s <- "HEL_BL_CLIN_CX_RAW"
dta <- as.data.frame(read_xlsx(f, sheet = s, na = c("", "NA")))
rownames(dta) <- dta$Pt
lg <- read_xlsx("data-raw/Evolution-december_2021_AI_updated.xlsx",
                na = c("", "NA", "-"))
names(lg) <- trimws(gsub("\\(.+\\)|\\r\\n", "", names(lg)))
rm(f, s)

# Rename variables
names(dta) <- sub("^1st", "First", names(dta))
names(dta) <- sub("/", "_", names(dta))

# New variables
dta$Partial_CHH <- ifelse(grepl("^Partial CHH", dta$Dx), 1, ifelse(
  grepl("^CDGP", dta$Dx), 0, NA))
dta$Complete_CHH <- ifelse(grepl("^Complete CHH", dta$Dx), 1, ifelse(
  grepl("^CDGP", dta$Dx), 0, NA))

# Add DX to longitudinal data
dx <- unique(na.omit(lg[c("ID", "DX")]))
dx$DX <- factor(sub(", Kallmann", "", dx$DX),
                c("CDGP", "Partial CHH", "Complete CHH"))
names(lg)[names(lg) == "DX"] <- "DX_complete"
lg <- merge(lg, dx, by = "ID", all.x = TRUE)
lg$PCHH <- as.numeric(lg$DX == "Partial CHH")
lg$CCHH <- as.numeric(lg$DX == "Complete CHH")
rm(dx)

# ROC curves
for (y in c("Partial_CHH", "Complete_CHH")) {
  X <- c("TV", "T", "LH", "FSH", "INB", "AMH", "IGF1")
  X <- setNames(as.list(paste0("First_", X)), X)
  names(X)[names(X) == "T"] <- "Testosterone"
  names(X)[names(X) == "LH"] <- "Basal LH"
  names(X)[names(X) == "FSH"] <- "Basal FSH"
  names(X)[names(X) == "INB"] <- "Inhibin B"
  if (y == "Complete_CHH") {
    X <- append(X, list(
      `Basal LH + Testosterone` = c("First_LH", "First_T"),
      `Basal LH + Crypt` = c("First_LH", "Crypt"),
      `TV + Crypt` = c("Crypt", "First_TV"),
      `Testosterone + Crypt` = c("First_T", "Crypt"),
      `TV + Micropenis` = c("Micro", "First_TV")
    ))
  }
  for (u in names(X)) {
    x <- X[[u]]
    sdta <- na.omit(dta[c(y, x)])
    if (length(x) == 1) {
      r <- suppressMessages(roc(sdta[[y]], sdta[[x]]))
      if (any(r$sensitivities >= 0.75 & r$specificities >= 0.75)) {
        outdir2 <- file.path(outdir, "ROC_curves_OK")
      } else {
        outdir2 <- file.path(outdir, "ROC_curves_NotOK")
      }
    } else {
      fml <- as.formula(paste(y, "~", paste(x, collapse = " + ")))
      fit <- suppressWarnings(glm(fml, family = binomial, data = sdta))
      r <- suppressMessages(
        roc(fit$model[[y]], predict(fit, type = "response")))
      outdir2 <- file.path(outdir, "ROC_curves_Multi")
    }
    if (!dir.exists(outdir2)) dir.create(outdir2)
    tbl <- as.data.frame(r[c("thresholds", "sensitivities", "specificities")])
    if (length(x) == 1) {
      names(tbl) <- c(x, "sensitivity", "specificity")
    } else {
      names(tbl) <- c("linear_combination", "sensitivity", "specificity")
    }
    tbl$youden_index <- tbl$sensitivity + tbl$specificity - 1
    tbl$youden_max <- tbl$youden_index == max(tbl$youden_index)
    a <- suppressWarnings(ci.auc(r))
    a <- sprintf("AUC %2$1.2f (%1$1.2f;%3$1.2f)", a[1], a[2], a[3])
    fig <- ggplot(tbl, aes(x = 1 - specificity, y = sensitivity)) +
      geom_point(aes(colour = youden_max)) +
      geom_path() +
      geom_abline(slope = 1, intercept = 0, colour = "grey70") +
      annotate("text", x = .5, y = .25, hjust =0, label = a) +
      scale_colour_manual(values = c("black", "red")) +
      lims(x = 0:1, y = 0:1) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position="none",
            plot.title = element_text(hjust = 0.5)) +
      labs(title = u)
    yxx <- paste0(y, "__", paste(x, collapse = "__"))
    yxx2 <- sub("^(P|C)(artial_|omplete_)(.+)$", "\\1\\3", yxx)
    f <- file.path(outdir2, paste0("ROC_curve__", yxx, ".xlsx"))
    write_xlsx(setNames(list(tbl), paste0(yxx2, " (n=", nrow(sdta), ")")), f)
    wth <- 3600
    hgh <- 3600
    res <- 1024
    tiff(filename = sub("\\.xlsx$", ".tiff", f), height = hgh, width = wth,
         res = res, compression = "zip")
    print(fig)
    dev.off()
    svglite(sub("\\.xlsx$", ".svg", f), width = wth / res, height = hgh / res)
    print(fig)
    dev.off()
    if (length(x) > 1) {
      m <- paste(sprintf("%#1.2f", coef(fit)), names(coef(fit)), sep = " * ")
      m <- paste(sub(" \\* \\(Intercept\\)", "", m), collapse = " + ")
      m <- paste("Linear_combination =", gsub("\\+ -", "- ", m))
      cat(m, file = sub("\\.xlsx$", ".txt", f))
    }
  }
}
rm(y, X, u, x, sdta, r, outdir2, fml, fit, tbl, a, fig, yxx, yxx2, f, m)

# Longitudinal analyses
Y <- c(TV = "TV (ml)", T = "Testosterone (nmol/l)", LH = "Basal LH (U/L)",
       FSH = "Basal FSH (U/L)", IGF1 = "IGF1", AMH = "AMH (pmol/L)",
       INB = "Inhibin B (pg/ml)")
for (y in names(Y)) {
  for (k in 1:3) {
    #print(paste(y, k))
    if (k == 3 & y %in% c("AMH", "INB")) next()
    sdta <- na.omit(lg[c("ID", "DX", "PCHH", "CCHH", "Age", y)])
    if (k %in% 1:2) {
      sdta <- subset(sdta, Age >= 14)
      fml <- list(~ DX * Age, ~ DX * (Age + I(Age^2) + I(Age^3)))[[k]]
      fml <- update(fml, as.formula(paste(y, "~ . + (1 | ID)")))
      ctrl <- lmerControl(optimizer ="Nelder_Mead")
      fit <- suppressWarnings(
        do.call("lmer", list(formula = fml, data = quote(sdta),
                             control = quote(ctrl))))
    } else {
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
      )
      fml <- ModelGradient(Age, Asym, Asym2, Asym3, xmid, scal, PCHH, CCHH) ~
        Asym | ID
      fml <- update(fml, as.formula(paste(y, "~ . ~ .")))
      fit <- suppressWarnings(
        do.call("nlmer", list(formula = fml, data = quote(sdta),
                              start = start_values[[y]])))
    }
    v <- c("term", "estimate", "std.error", "conf.low", "conf.high")
    if (k %in% 1:2) v <- c(v, "p.value")
    tbl <- suppressWarnings(tidy(fit, conf.int = TRUE)[v])
    tbl$term = sub("^DX", "", tbl$term)
    tbl$term = sub(":", " x ", tbl$term)
    tbl$term = sub("sd__", "SD ", tbl$term)
    if (k %in% 1:2) {
      ndta <- as.data.frame(ggpredict(fit, c("Age [all]", "DX")))
      ndta <- ndta[c("x", "predicted", "conf.low", "conf.high", "group")]
      names(ndta)[names(ndta) == "x"] <- "Age"
      names(ndta)[names(ndta) == "group"] <- "DX"
    } else {
      ndta <- expand.grid(
        Age = seq(14, 20, 0.1),
        DX = c("CDGP", "Partial CHH", "Complete CHH")
      )
      ndta$PCHH = as.numeric(ndta$DX == "Partial CHH")
      ndta$CCHH = as.numeric(ndta$DX == "Complete CHH")
      ndta <- cbind(ndta, as.data.frame(t(setNames(tbl[["estimate"]],
                                                   tbl[["term"]]))))
      ndta$predicted <- with(ndta, Model(Age, Asym, Asym2, Asym3, xmid,
                                         scal, PCHH, CCHH))
    }
    ndta <- subset(ndta, Age >= 14 & Age <= 19)
    fig <- ggplot(ndta, aes(x = Age, y = predicted)) +
      geom_line(aes(colour = DX)) +
      lims(x = c(14, 19)) +
      labs(x = "Age", y = Y[y], title = "Fixed effects predictions") +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "bottom",
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5))
    if (k == 3) {
      tbl$term = sub("^Asym2", "AsymPartCHH", tbl$term)
      tbl$term = sub("^Asym3", "AsymCompCHH", tbl$term)
    }
    outdir2 <- file.path(outdir, "longitudinal_analyses")
    if (!dir.exists(outdir2)) dir.create(outdir2)
    f <- paste0(y, "_", c("linear", "polynomial", "sigmoidal")[k], ".xlsx")
    f <- file.path(outdir2, f)
    write_xlsx(setNames(list(tbl), paste0(y, " (n=", nrow(sdta), ")")), f)
    f <- sub("\\.xlsx$", "_without_CI.tiff", f)
    wth <- 5400
    hgh <- 3600
    res <- 1024
    tiff(filename = f, height = hgh, width = wth, res = res,
         compression = "zip")
    print(fig)
    dev.off()
    svglite(sub("\\.tiff$", ".svg", f), width = wth / res, height = hgh / res)
    print(fig)
    dev.off()
    if (k %in% 1:2) {
      fig <- fig +
        geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = DX),
                    alpha = 0.3, show.legend = FALSE)
      f <- sub("_without", "_with", f)
      tiff(filename = f, height = hgh, width = wth, res = res,
         compression = "zip")
      print(fig)
      dev.off()
      svglite(sub("\\.tiff$", ".svg", f), width = wth / res,
              height = hgh / res)
      print(fig)
      dev.off()
    }
  }
}
rm(Y, y, k, sdta, fml, ctrl, fit, Model, ModelGradient, start_values, v, tbl,
   ndta, fig, outdir2, f, wth, hgh, res)

# Session Info
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
