library(ggplot2)
library(magrittr)
library(pander)
library(parallel)
library(pROC)
library(readxl)
library(writexl)

options(mc.cores = detectCores() - 1)

# Set working directory
setwd("~/Projects/Consultations/Pitteloud Nelly (Puberté)")

# Output directory
outdir <- paste0("results/analyses_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(outdir)) dir.create(outdir)

# Load data
f <- "data-raw/HEL_Database_october2021_updated_withoutnames_stats.xlsx"
dta <- as.data.frame(read_xlsx(f, na = c("", "NA")))
rownames(dta) <- dta$Pt
rm(f)

# Rename variables
names(dta) <- sub("^1st", "First", names(dta))
names(dta) <- sub("/", "_", names(dta))

# Predictors
X <- names(dta)[!(grepl("^(Pt|Dx)$", names(dta)))]

# Recoding
for (j in which(sapply(dta, class) == "character")) {
  if (all(is.na(dta[[j]]) | dta[[j]] %in% c("N", "Y"))) {
    dta[[j]] <- c(N = FALSE, Y = TRUE)[dta[[j]]]
  }
}
dta$Olfaction <- factor(dta$Olfaction, c("Normosmic", "Anosmic", "Hyposmia"))
dta$Olfaction2 <- dta$Olfaction
levels(dta$Olfaction2) <- c("Normosmic", rep("Anosmic/Hyposmia", 2))
#X[X == "Olfaction"] <- "Olfaction2"
rm(j)

# New variables
dta$Partial_CHH <- ifelse(grepl("^Partial CHH", dta$Dx), 1, ifelse(
  grepl("^CDGP", dta$Dx), 0, NA))
dta$Complete_CHH <- ifelse(grepl("^Complete CHH", dta$Dx), 1, ifelse(
  grepl("^CDGP", dta$Dx), 0, NA))
dta$CHH <- ifelse(grepl("CHH", dta$Dx), 1,
                  ifelse(grepl("^CDGP", dta$Dx), 0, NA))

# Responses
Y <- c("CHH", "Partial_CHH", "Complete_CHH")

# AUC - List of models
f <- file.path(outdir, "AUC_tables_model_list.rds")
if (file.exists(f)) {
  M <- readRDS(f)
} else {
  M <- unlist(recursive = FALSE, lapply(Y, function(y) {
    print(y)
    clps <- function(z) paste(z, collapse = " / ")

    x_list <- unlist(lapply(1:3, function(k) combn(X, k, simplify = FALSE)),
                     recursive = FALSE)
    x_list <- Reduce(append, lapply(x_list, function(x) {
      if (any(x == "Olfaction")) {
        x <- list(x, sub("^Olfaction$", "Olfaction2", x))
      } else {
        x <- list(x)
      }
      return(x)
    }))
    mclapply(x_list, function(x) {
      #print(paste(c(y, x), collapse = " "))
      sdta <- na.omit(dta[c(y, x)])
      fml <- paste(y, "~", paste(x, collapse = " + "))
      n <- nrow(sdta)
      npos <- sum(sdta[[y]] == 1)
      nneg <- sum(sdta[[y]] == 0)
      resp <- sdta[[y]]
      e <- paste0("glm(", fml, ", family = binomial, data = sdta)")
      r <- evals(e, env = environment(), graph.dir = "/tmp")[[1]]
      fit.wrn <- unique(r$msg$warnings)
      fit.wrn <- if (is.null(fit.wrn)) NA else clps(fit.wrn)
      fit.err <- unique(r$msg$errors)
      fit.err <- if (is.null(fit.err)) NA else clps(fit.err)
      if (is.na(fit.err)) {
        fit <- r$result
        e <- "predict(fit, type = 'response')"
        r <- evals(e, env = environment(), graph.dir = "/tmp")[[1]]
        pred.wrn <- unique(r$msg$warnings)
        pred.wrn <- if (is.null(pred.wrn)) NA else clps(pred.wrn)
        pred.err <- unique(r$msg$errors)
        pred.err <- if (is.null(pred.err)) NA else clps(pred.err)
        if (is.na(pred.err)) {
          prob <- r$result
        } else {
          prob <- NULL
        }
        z <- suppressMessages(roc(resp, prob, direction = "<"))
        auc <- as.numeric(z$auc)
        r <- evals("ci.auc(z)", env = environment(), graph.dir = "/tmp")[[1]]
        auc.ci <- r$result
        auc.ci.wrn <- unique(r$msg$warnings)
        auc.ci.wrn <- if (is.null(auc.ci.wrn)) NA else clps(auc.ci.wrn)
        loo <- lapply(1:nrow(sdta), function(i) {
          e <- "glm(fml, family = binomial, data = sdta[-i, ])"
          r <- evals(e, env = environment(), graph.dir = "/tmp")[[1]]
          fit.wrn <- unique(r$msg$warnings)
          fit.wrn <- if (is.null(fit.wrn)) NA else clps(fit.wrn)
          fit.err <- unique(r$msg$errors)
          fit.err <- if (is.null(fit.err)) NA else clps(fit.err)
          if (is.na(fit.err)) {
            fit <- r$result
            e <- paste("predict(fit, newdata = sdta[i, , drop = FALSE],",
                       "type = 'response')")
            r <- evals(e, env = environment(), graph.dir = "/tmp")[[1]]
            pred.wrn <- unique(r$msg$warnings)
            pred.wrn <- if (is.null(pred.wrn)) NA else clps(pred.wrn)
            pred.err <- unique(r$msg$errors)
            pred.err <- if (is.null(pred.err)) NA else clps(pred.err)
            if (is.na(pred.err)) {
              prob <- r$result
            } else {
              prob <- NA
            }
          } else {
            prob <- NA
          }
          list(prob = prob, fit.wrn = fit.wrn, fit.err = fit.err,
               pred.wrn = pred.wrn, pred.err = pred.err, fit = fit)
        })
        loo.prob <- sapply(loo, function(z) z$prob)
        loo.fit.wrn <- na.omit(unique(sapply(loo, function(z) z$fit.wrn)))
        loo.fit.wrn <- if (length(loo.fit.wrn) == 0) NA else clps(loo.fit.wrn)
        loo.fit.err <- na.omit(unique(sapply(loo, function(z) z$fit.err)))
        loo.fit.err <- if (length(loo.fit.err) == 0) NA else clps(loo.fit.err)
        loo.pred.wrn <- na.omit(unique(sapply(loo, function(z) z$pred.wrn)))
        loo.pred.wrn <- if (length(loo.pred.wrn) == 0) NA else clps(loo.pred.wrn)
        loo.pred.err <- na.omit(unique(sapply(loo, function(z) z$pred.err)))
        loo.pred.err <- if (length(loo.pred.err) == 0) NA else clps(loo.pred.err)
        if (is.na(loo.fit.err) & is.na(loo.pred.err)) {
          z <- suppressMessages(roc(resp, loo.prob, direction = "<"))
          loo.auc <- as.numeric(z$auc)
          r <- evals("ci.auc(z)", env = environment(), graph.dir = "/tmp")[[1]]
          loo.auc.ci <- r$result
          w <- unique(r$msg$warnings)
          loo.auc.ci.wrn <- if (is.null(w)) NA else clps(w)
        } else {
          loo.auc <- NA
          loo.auc.ci <- rep(NA, 3)
          loo.auc.ci.wrn <- NA
        }
      } else {
        prob <- NULL
        auc <- NA
        auc.ci <- rep(NA, 3)
        auc.ci.wrn <- NA
        loo.prob <- NULL
        loo.fit.wrn <- NA
        loo.fit.err <- NA
        loo.pred.wrn <- NA
        loo.pred.err <- NA
        loo.auc <- NA
        loo.auc.ci <- rep(NA, 3)
        loo.auc.ci.wrn <- NA
      }
      list(
        row = data.frame(
          Response = y,
          Predictor.1 = x[1],
          Type.1 = class(sdta[[x[1]]]),
          Predictor.2 = x[2],
          Type.2 = class(sdta[[x[2]]]),
          Predictor.3 = x[3],
          Type.3 = class(sdta[[x[3]]]),
          N = n,
          Npos = npos,
          Nneg = nneg,
          AUC = auc,
          AUC.lwr = auc.ci[1],
          AUC.upr = auc.ci[3],
          LOO.AUC = loo.auc,
          LOO.AUC.lwr = loo.auc.ci[1],
          LOO.AUC.upr = loo.auc.ci[3],
          fit.wrn = fit.wrn,
          fit.err = fit.err,
          pred.wrn = pred.wrn,
          pred.err = pred.err,
          auc.ci.wrn = auc.ci.wrn,
          loo.fit.wrn = loo.fit.wrn,
          loo.fit.err = loo.fit.err,
          loo.pred.wrn = loo.pred.wrn,
          loo.pred.err = loo.pred.err,
          loo.auc.ci.wrn = loo.auc.ci.wrn
        ),
        resp = resp,
        prob = prob,
        loo.prob = loo.prob
      )
    })
  }))
  for (i in 1:length(M)) {
    M[[i]]$model <- i
    M[[i]]$row <- cbind(Model = i, M[[i]]$row)
  }
  rm(i)
  saveRDS(M, f)
}
rm(f)

# AUC - Table
tbl <- do.call(rbind, lapply(M, function(z) z$row))
lapply(setNames(Y, Y), function(y) {
  tbl <- tbl[tbl$Response == y, ]
  tbl[order(tbl$LOO.AUC, decreasing = TRUE), ]
}) %>%
  write_xlsx(file.path(outdir, "AUC_tables.xlsx"))

# Help function
get_fit <- function(m) {
  y <- m$row$Response
  x <- m$row[paste0("Predictor.", 1:3)]
  x <- x[!is.na(x)]
  dta <- na.omit(dta[c(y, x)])
  fml <- as.formula(paste(y, "~", paste(x, collapse = " + ")))
  do.call("glm", list(formula = fml, family = quote(binomial),
                      data = quote(dta)))
}

# Save image
save.image(file = file.path(outdir, "AUC_tables_workspace.rda"))

###############################################################################
# resp <- sdta[[y]]
# prob <- predict(fit, type = "response")
#
# prob_loocv <- sapply(1:nrow(dta), function(i) {
#   fit <- glm(fml, family = "binomial", data = dta[-i, ])
#   predict(fit, newdata = dta[i, , drop = FALSE], type = "response")
# })
#
# fml <- CHH ~ Age_1st + First_Ht + First_BMI + First_TV
# fit <- glm(fml, family = "binomial", data = dta)
# x <- predict(fit, type = "response")
# r <- roc(response = dta$CHH, predictor = x)
# x <- sapply(1:nrow(dta), function(i) {
#   fit <- glm(fml, family = "binomial", data = dta[-i, ])
#   predict(fit, newdata = dta[i, , drop = FALSE], type = "response")
# })
# roc(response = dta$CHH, predictor = x)
#
# library(caret)
# dta$y <- factor(paste0("X", dta$CHH))
# fml_2 <- update(fml, y ~ .)
# ctrl <- trainControl(method="LOOCV",
#                      summaryFunction=twoClassSummary,
#                      classProbs=T,
#                      savePredictions = T)
# fit <- train(fml_2, data=dta,
#                method="glm", family = "binomial",
#                trControl=ctrl)
# fit

