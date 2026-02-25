#packages ####
library(tidyverse)
library(readxl)
library(vegan)
library(lme4)
library(ggstatsplot)
library(MASS)
library(car)
library(gridExtra)
library(ordinal)
library(lmerTest)
library(lmeresampler)
library(ggeffects)
library(Hmsc)

setwd("C:/Users/Sam/OneDrive/Documents")
#Interactions new attempt. this is computationally intensive####

#initially we used a less robust, but less computationally intensive,
# approach for approximating species interactions (Veech's co-occurrence).
# Unfortunately, 

data <- read_excel("Data.for.EcoApp25.xlsx", sheet="Build")

# fill NA with 0s for community
data[, 53:ncol(data)][is.na(data[, 53:ncol(data)])] <- 0

#make community dataset (already wide)
comm <- data[,53:ncol(data)]

#also random effects
rand <- data[,c("L3 Ecoregion" , "Year" )]

# env covariates
env <- data[,c("Temp (°C)" , "Turbidity (NTU)",   "DO (mg/L)" ,
               "pH (SU)" , "Conductivity" , "TKN (mg/L)" , "Total Phosphorus (mg/L)" ,
               "Ammonia (mg/L)" ,  "Nitrate (mg/L)"  , "Total Points" )]

#impute a couple values, there's not a lot of missing values so we 
# can just impute them with mean values
env_imp <- apply(env, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

#run pca
pca <- prcomp(env_imp, scale. = TRUE, center = TRUE)

#get principle components as covariates
env.sum <- data.frame(pca$x[,1:4])

#remove rare taxa, they'll automatically get assigned as low
# interacting taxa later
prev <- colSums(comm > 0) / nrow(comm)
comm.final <- comm[, prev >= 0.01, drop = FALSE]  # e.g., taxa present in >=5% samples


#start prepping for the run

#random effects
rand <- data.frame(rand)
rand$Region <- as.factor(rand$L3.Ecoregion)
rand$Year <- as.factor(rand$Year)
rand <- data.frame(rand)
rL.year  <- HmscRandomLevel(units = rand$Year)
rL.basin <- HmscRandomLevel(units = rand$Region)


#create object
m <- Hmsc(
  Y = as.matrix(comm.final),
  XData = env.sum,
  XFormula = ~ PC1+PC2+PC3+PC4,
  studyDesign = rand[,c("Year","Region")],
  ranLevels = list(Year = rL.year, Region = rL.basin),
  distr = "lognormal poisson"
)


# Fit (tune these for your dataset size)
m2 <- sampleMcmc(
  m,
  thin = 50,
  samples = 1000,
  transient = 500,
  nChains = 3,
  verbose = 100
)

saveRDS(m2, "hmscfit.rds")


#diagnostics

m2 <- readRDS("hmscfit.rds")
plot(m2[[2]], parameters = c("Beta", "Omega"))

preds <- computePredictedValues(m2)
MF <- evaluateModelFit(hM = m2, predY = preds)
summary(MF$RMSE)
summary(MF$SR2)

vp <- computeVariancePartitioning(m2)
plotVariancePartitioning(m2, vp)

mcmc <- convertToCodaObject(m2)
effective <- effectiveSize(mcmc$V)
quantile(effective, probs = c(0.05, 0.25, 0.5))

# Gelman-Rubin R-hat (PSRF) for fixed effects
gr_beta <- gelman.diag(mcmc$Beta, autoburnin = FALSE, multivariate = FALSE)

# Vector of R-hat values (one per parameter)
rhat_beta <- gr_beta$psrf[, "Point est."]

# Summaries
summary(rhat_beta)

m3 <- computeAssociations(m2)

saveRDS(m3, "associations.rds")

View(m3[[1]]$mean)


#calculate # of associations #####

m3 <- readRDS("associations.rds")
count_links_per_taxon <- function(S, hi = 0.95, lo = 0.05) {
  # S is the support matrix: P(association > 0)
  diag(S) <- NA  # don't count self
  
  pos <- S >= hi
  neg <- S <= lo
  any <- pos | neg
  
  data.frame(
    taxon   = rownames(S),
    n_pos   = rowSums(pos, na.rm = TRUE),
    n_neg   = rowSums(neg, na.rm = TRUE),
    n_total = rowSums(any, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

links_RE1 <- count_links_per_taxon(m3[[1]]$support)  # e.g., year
links_RE2 <- count_links_per_taxon(m3[[2]]$support)  # e.g., region


associates <- links_RE2

rownames(associates) <- NULL

#Summarize across years/regions
hist(m3[[2]]$mean)

#Read in data, substitute pathways if needed #####

# #older files. files have been combined into one excel file
# genera <- read_excel("ecoBI.xlsx", sheet="Genus2")
# family <- read_excel("ecoBI.xlsx", sheet="Family2")
# allmech <- read_excel("ArkansasValleyOP.xlsx", sheet = "All")
# testdata <- read_excel("indextestdata.xlsx", sheet="Sheet1")

genera <- read_excel("Data.for.EcoApp25.xlsx", sheet="Genus2")
family <- read_excel("Data.for.EcoApp25.xlsx", sheet="Family2")
allmech <- read_excel("Data.for.EcoApp25.xlsx", sheet = "AllMech")
testdata <- read_excel("Data.for.EcoApp25.xlsx", sheet="testdata")


hist(testdata$TN)
#Assign 
generaAll <- genera %>% filter(Ecoregion == "All")
generaAll$scoreTN<- as.integer(ntile(generaAll$TKN,10))
generaAll$scoreTP<- as.integer(ntile(generaAll$TP,10))

familyAll <- family %>% filter(Ecoregion == "All")
familyAll$scoreTN<- as.integer(ntile(familyAll$TKN,10))
familyAll$scoreTP<- as.integer(ntile(familyAll$TP,10))

#Combine data
allindexscores <- rbind(familyAll,generaAll)
allindexscores <- allindexscores[,c(1,2,17:ncol(allindexscores))]

#transform environmental variables
testdata$TP <- log10(testdata$TP)
testdata$TN <- log10(testdata$TN)


#Combine the test data and the scores
allindexscores <- left_join(testdata,allindexscores, by="Species", na_matches = "never")

#also good to add the mechanisms to the dataset, but first we create a new
#column, dispersal, which factors in whether they disperse aquatically
#or terrestrially
allmech$Dispersal <- allmech$SFP #*ifelse(allmech$DispersalType=="Aquatic",0.5,1)

allmech$Dispersal <- ifelse(allmech$Dispersal == 0, mean(allmech$Dispersal),allmech$Dispersal)

#allmech$Dispersal <- allmech$Dispersal *ifelse(allmech$DispersalType=="Aquatic",0.5,1)

allindexscores <- left_join(allindexscores, allmech, by=join_by(Species==Taxa),na_matches = "never")

#Begin index calculation: Score x log-transformed Count
allindexscores <- allindexscores %>% mutate(TNscore = scoreTN*log(1+Count))
allindexscores <- allindexscores %>% mutate(TPscore =scoreTP*log(1+Count))
allindexscores <- allindexscores %>% mutate(logcount = log(1+Count))


#add new associations
allindexscores <- left_join(allindexscores, associates, by=join_by(Species == taxon))

allindexscores$n_neg <- ifelse(is.na(allindexscores$n_neg),0,allindexscores$n_neg)
allindexscores$n_total <- ifelse(is.na(allindexscores$n_total),0,allindexscores$n_total)


#continuous penalty with HMSC ####
rev_minmax <- function(x, trim = 0.01) {
  # trim = winsorization proportion on each tail (set 0 for none)
  lo <- as.numeric(quantile(x, trim, na.rm = TRUE))
  hi <- as.numeric(quantile(x, 1 - trim, na.rm = TRUE))
  xw <- pmin(pmax(x, lo), hi)  # winsorize
  if (hi - lo == 0) return(rep(1, length(x)))
  (hi - xw) / (hi - lo)        # high values -> 0, low values -> 1
}

rev_minmax_half <- function(x, trim = 0.0) {
  w <- rev_minmax(x, trim)
  0.5 + 0.5 * w
}

mid_minmax_half <- function(x, trim = 0.0, center = 0.5) {
  lo <- as.numeric(quantile(x, trim, na.rm = TRUE))
  hi <- as.numeric(quantile(x, 1 - trim, na.rm = TRUE))
  xw <- pmin(pmax(x, lo), hi)  # winsorize
  
  if (hi - lo == 0) return(rep(1, length(x)))
  
  # normalize to [0,1]
  z <- (xw - lo) / (hi - lo)
  
  # distance from center (default 0.5 = midpoint of range)
  d <- abs(z - center) / max(center, 1 - center)  # scaled so max distance = 1
  
  # convert distance to weight: center -> 1, extremes -> 0.5
  1 - 0.5 * d
}
rev_minmax_01 <- function(x, min_w = 0.3, max_w = 1, trim = 0.05) {
  lo <- as.numeric(quantile(x, trim, na.rm = TRUE))
  hi <- as.numeric(quantile(x, 1 - trim, na.rm = TRUE))
  
  xw <- pmin(pmax(x, lo), hi)  # winsorize
  
  if (hi - lo == 0) return(rep(max_w, length(x)))
  
  w <- (hi - xw) / (hi - lo)  # reversed
  min_w + w * (max_w - min_w)
}

allindexscores <- allindexscores %>%
  mutate(
    w_I = rev_minmax_half(n_total),         # associations (I)
    w_D = rev_minmax_half(Dispersal),     # dispersal (D)
    w_N = rev_minmax_half(NicheBreadth)   # niche breadth (N)
  )

allindexscores <- allindexscores %>%
  mutate(
    logcount = log1p(Count),
    TN_base  = scoreTN * logcount,
    TP_base  = scoreTP * logcount
  )

combos <- list(
  base = c(),
  D    = c("w_D"),
  I    = c("w_I"),
  N    = c("w_N"),
  DI   = c("w_D","w_I"),
  DN   = c("w_D","w_N"),
  IN   = c("w_I","w_N"),
  DNI  = c("w_D","w_I","w_N")
)

mul_weights <- function(df, vars) {
  if (length(vars) == 0) return(rep(1, nrow(df)))
  Reduce(`*`, df[vars])
}

for (nm in names(combos)) {
  w <- mul_weights(allindexscores, combos[[nm]])
  allindexscores[[paste0("LC_", nm)]] <- allindexscores$logcount * w
  allindexscores[[paste0("TN_", nm)]] <- allindexscores$TN_base  * w
  allindexscores[[paste0("TP_", nm)]] <- allindexscores$TP_base  * w
}


library(dplyr)

# (A) mean site covariates (keep your columns; just showing your pattern)
testdata1 <- allindexscores[, c(1, 6, 7, 12:21)] %>%
  group_by(UID) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")

# (B) sum index components (all the LC_/TP_/TN_ columns)


testdata2 <- allindexscores %>%
  dplyr::select(UID, starts_with("LC_"), starts_with("TP_"), starts_with("TN_")) %>%
  group_by(UID) %>%
  summarise(across(everything(), ~sum(.x, na.rm = TRUE)), .groups = "drop")

testfinal <- left_join(testdata1, testdata2, by = "UID")


variants <- c("base","D","I","N","DI","DN","IN","DNI")

for (v in variants) {
  # names
  tp <- paste0("TP_", v)
  tn <- paste0("TN_", v)
  lc <- paste0("LC_", v)
  
  # avoid NaN/Inf when lc == 0
  testfinal[[paste0("TPscore_", v)]] <- with(testfinal, ifelse(get(lc) > 0, get(tp) / get(lc), NA_real_))
  testfinal[[paste0("TNscore_", v)]] <- with(testfinal, ifelse(get(lc) > 0, get(tn) / get(lc), NA_real_))
}


test <- testfinal %>%
  filter(!is.na(TPscore_base)&
           !is.na(TP_DNI))


tp_models <- lapply(variants, function(v) {
  pred <- paste0("TPscore_", v)
  dat  <- test[, c("TP", pred)]
  dat  <- na.omit(dat)
  lm(reformulate(pred, response = "TP"), data = dat)
})
names(tp_models) <- paste0("TP_", variants)

tp_models$TP_null <- lm(TP ~ 1, data = na.omit(test[, "TP", drop = FALSE]))

TPAICBIC <- data.frame(
  model_id = names(tp_models),
  AIC = sapply(tp_models, AIC),
  BIC = sapply(tp_models, BIC),
  row.names = NULL
)


test <- testfinal %>%
  filter(!is.na(TNscore_base)&!is.na(TN_DNI))


tn_models <- lapply(variants, function(v) {
  pred <- paste0("TNscore_", v)
  dat  <- test[, c("TN", pred)]
  dat  <- na.omit(dat)
  lm(reformulate(pred, response = "TN"), data = dat)
})

names(tn_models) <- paste0("TN_", variants)

tn_models$TN_null <- lm(TN ~ 1, data = na.omit(test[, "TN", drop = FALSE]))

TNAICBIC <- data.frame(
  model_id = names(tn_models),
  AIC = sapply(tn_models, AIC),
  BIC = sapply(tn_models, BIC),
  row.names = NULL
)

mae  <- function(actual, pred) mean(abs(actual - pred), na.rm = TRUE)
rmse <- function(actual, pred) sqrt(mean((actual - pred)^2, na.rm = TRUE))

get_metrics <- function(mod) {
  mf   <- model.frame(mod)
  y    <- model.response(mf)
  yhat <- fitted(mod)
  s    <- summary(mod)
  
  data.frame(
    model_name = deparse1(formula(mod)),
    r2         = unname(s$r.squared),
    adj_r2     = unname(s$adj.r.squared),
    MAE        = mae(y, yhat),
    RMSE       = rmse(y, yhat),
    nobs       = nrow(mf),
    stringsAsFactors = FALSE
  )
}

get_variable_metrics <- function(model) {
  summ <- summary(model)
  coefs <- summ$coefficients
  
  if (nrow(coefs) < 2) {
    data.frame(
      estimate = NA_real_,
      std_error = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_
    )
  } else {
    data.frame(
      estimate = coefs[2, 1],
      std_error = coefs[2, 2],
      t_value = coefs[2, 3],
      p_value = coefs[2, 4]
    )
  }
}






TPmetrics <- do.call(rbind, lapply(names(tp_models), function(id) {
  mod <- tp_models[[id]]
  out  <- get_metrics(mod)
  out2 <- get_variable_metrics(mod)
  out$model_id <- id
  cbind(out, out2)
}))



TPall <- dplyr::left_join(TPAICBIC, TPmetrics, by = "model_id")


TNmetrics <- do.call(rbind, lapply(names(tn_models), function(id) {
  mod <- tn_models[[id]]
  out  <- get_metrics(mod)
  out2 <- get_variable_metrics(mod)
  out$model_id <- id
  cbind(out, out2)
}))

TNall <- dplyr::left_join(TNAICBIC, TNmetrics, by = "model_id")


View(TNall)
View(TPall)

#continuous plots ####
library(ggeffects)

TPmodel.initial <- ggpredict(tp_models$TP_base,terms = "TPscore_base [0:10]")
TPmodel.corrected <-ggpredict(tp_models$TP_DI,terms = "TPscore_DI [0:10]")
TPmodel.universal <-ggpredict(tp_models$TP_DNI,terms = "TPscore_DNI [0:10]")
initialTP <- ggplot()+  geom_point(data=tp_models$TP_base$model,
                                   aes(x=TPscore_base,y=TP, color="red3"),alpha = 0.6)+
  geom_ribbon(data=TPmodel.initial, aes(x=x,ymin = conf.low,ymax=conf.high),fill = "red3",alpha=0.3)+
  stat_smooth(data=TPmodel.initial, aes(x=x,y=predicted),
              se=F,method = "glm",color = "red3", fill="red3")+ 
  xlab("Initial index")+ylab(expression("Log"["10"]*" transformed total phosphorus ( "*mu*"g/L)"))+
  
  # geom_smooth(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),method="glm",
  #             se=F,aes(x=TPscore,y=TP),size=0.1,linetype="dashed")+
  # geom_smooth(data=TPmodel.initial, se=F,aes(x=x,y=predicted),color="red3")+
  ylim(-0.5,4)+
  labs(title="(a)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()


correctedTP <- ggplot()+  geom_point(data=tp_models$TP_DI$model,
                                     aes(x=TPscore_DI,y=TP), color="green3",alpha = 0.6)+
  geom_ribbon(data=TPmodel.corrected ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "green3",alpha=0.3)+
  stat_smooth(data=TPmodel.corrected, aes(x=x,y=predicted),
              se=F,method = "glm",color = "green3", fill="green3")+ 
  xlab("Top modified index")+ylab('')+
  ylim(-0.5,4)+
  xlim(0,10)+labs(title="(b)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()


universalTP <- ggplot()+  geom_point(data=tp_models$TP_DNI$model,
                                     aes(x=TPscore_DNI,y=TP), color="black",alpha = 0.6)+
  geom_ribbon(data=TPmodel.universal ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "black",alpha=0.3)+
  stat_smooth(data=TPmodel.universal, aes(x=x,y=predicted),
              se=F,method = "glm",color = "black", fill="black")+ 
  xlab("Suggested modified index")+ylab('')+
  ylim(-0.5,4)+
  xlim(0,10)+labs(title="(c)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()




gridExtra::grid.arrange(initialTP,correctedTP,universalTP,ncol=3)





TNmodel.initial <- ggpredict(TNModel,terms = "TNscore [0:10]")
TNmodel.corrected <-ggpredict(TNModel.D,terms = "TNscore.D [0:10]")
TNmodel.universal <-ggpredict(TNModel.D.5N.5I.5,terms = "TNscore.D.5N.5I.5 [0:10]")
initialTN <- ggplot()+  geom_point(data=na.omit(test[,c("TNscore","TNscore.D.5","TN")]),
                                   aes(x=TNscore,y=TN, color="red3"),alpha = 0.6)+
  geom_ribbon(data=TNmodel.initial, aes(x=x,ymin = conf.low,ymax=conf.high),fill = "red3",alpha=0.3)+
  stat_smooth(data=TNmodel.initial, aes(x=x,y=predicted),
              se=F,method = "glm",color = "red3")+ 
  xlab("Initial index")+ylab(expression("Log"["10"]*" transformed total nitrogen (mg/L)"))+
  
  # geom_smooth(data=na.omit(test[,c("TNscore","TNscore.D.5","TN")]),method="glm",
  #             se=F,aes(x=TNscore,y=TN),size=0.1,linetype="dashed")+
  # geom_smooth(data=TNmodel.initial, se=F,aes(x=x,y=predicted),color="red3")+
  labs(title="(d)")+  ylim(-1.6,1.6)+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()

correctedTN <- ggplot()+  geom_point(data=na.omit(test[,c("TNscore","TNscore.D.5","TN")]),
                                     aes(x=TNscore.D.5,y=TN), color="green3",alpha = 0.6)+
  geom_ribbon(data=TNmodel.corrected ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "green3",alpha=0.3)+
  stat_smooth(data=TNmodel.corrected, aes(x=x,y=predicted),
              se=F,method = "glm",color = "green3")+ 
  xlab("Top modified index")+ylab('')+
  ylim(-1.6,1.6)+
  xlim(0,10)+labs(title="(e)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()


universalTN <- ggplot()+  geom_point(data=na.omit(test[,c("TNscore","TNscore.D.5","TNscore.D.5N.5I.5","TN")]),
                                     aes(x=TNscore.D.5N.5I.5,y=TN), color="black",alpha = 0.6)+
  geom_ribbon(data=TNmodel.universal ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "black",alpha=0.3)+
  stat_smooth(data=TNmodel.universal, aes(x=x,y=predicted),
              se=F,method = "glm",color = "black", fill="black")+ 
  xlab("Suggested modified index")+ylab('')+
  ylim(-1.6,1.6)+
  xlim(0,10)+labs(title="(f)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()




gridExtra::grid.arrange(initialTP,correctedTP,universalTP,
                        initialTN,correctedTN,universalTN,ncol=3)


test$TNactual <- 10^test$TN
test$TPactual <- 10^test$TP

test <- test %>%
  mutate(
    TN_class = case_when(
      TNactual <= 0.66 ~ "Low",
      TNactual <= 3.17 ~ "Moderate",
      TNactual > 3.17       ~ "High"
    ),
    TP_class = case_when(
      TPactual <= 50   ~ "Low",
      TPactual <= 280  ~ "Moderate",
      TPactual > 280     ~ "High"
    ),
    TN_class = factor(TN_class, levels = c("Low", "Moderate", "High"), ordered = TRUE),
    TP_class = factor(TP_class, levels = c("Low", "Moderate", "High"), ordered = TRUE)
  ) %>%
  filter(!is.na(TN_class) & !is.na(TP_class))


#effect size 
library(dplyr)
library(effsize)

# Select TP score columns
tp_score_cols <- names(test)[grepl("^TPscore_", names(test))]

rho_tp <- map_dfr(tp_score_cols, function(col) {
  
  df <- test %>%
    select(TPactual, all_of(col)) %>%
    filter(!is.na(TPactual), !is.na(.data[[col]]))
  
  if(nrow(df) < 3) {
    return(tibble(
      score = col,
      n = nrow(df),
      rho = NA_real_,
      p_value = NA_real_
    ))
  }
  
  ct <- cor.test(df$TPactual,
                 df[[col]],
                 method = "spearman",
                 exact = FALSE)
  
  tibble(
    score = col,
    n = nrow(df),
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}) %>%
  arrange(desc(abs(rho)))

rho_tp


tp_score_cols <- names(test)[grepl("^TPscore_", names(test))]


pairwise_auc <- function(data, class_col, score_prefix, pairs = list(c("Low","Moderate"),
                                                                     c("Low","High"),
                                                                     c("Moderate","High"))) {
  class_col <- rlang::ensym(class_col)
  score_cols <- names(data)[grepl(paste0("^", score_prefix), names(data))]
  
  purrr::map_dfr(score_cols, function(sc) {
    purrr::map_dfr(pairs, function(pr) {
      a <- pr[1]; b <- pr[2]
      
      df <- data %>%
        select(!!class_col, all_of(sc)) %>%
        filter(!is.na(!!class_col), !is.na(.data[[sc]])) %>%
        filter((!!class_col) %in% c(a, b)) %>%
        mutate(y = as.integer((!!class_col) == b))  # b is the "positive" class
      
      # need both classes present
      if (nrow(df) < 3 || length(unique(df$y)) < 2) {
        return(tibble(score = sc, contrast = paste(a, "vs", b), n = nrow(df),
                      auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
      }
      
      r <- pROC::roc(response = df$y, predictor = df[[sc]], quiet = TRUE)
      ci <- as.numeric(pROC::ci.auc(r))
      
      tibble(
        score = sc,
        contrast = paste(a, "vs", b),
        n = nrow(df),
        auc = as.numeric(pROC::auc(r)),
        ci_low = ci[1],
        ci_high = ci[3]
      )
    })
  }) %>%
    arrange(desc(auc))
}


auc_tp_pairs <- pairwise_auc(test, TP_class, score_prefix = "TPscore_")
auc_tp_pairs
View(auc_tp_pairs)
tn_score_cols <- names(test)[grepl("^TNscore_", names(test))]

rho_tn <- map_dfr(tn_score_cols, function(col) {
  
  df <- test %>%
    select(TNactual, all_of(col)) %>%
    filter(!is.na(TNactual), !is.na(.data[[col]]))
  
  if(nrow(df) < 3) {
    return(tibble(
      score = col,
      n = nrow(df),
      rho = NA_real_,
      p_value = NA_real_
    ))
  }
  
  ct <- cor.test(df$TNactual,
                 df[[col]],
                 method = "spearman",
                 exact = FALSE)
  
  tibble(
    score = col,
    n = nrow(df),
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}) %>%
  arrange(desc(abs(rho)))

rho_tn

auc_tn_pairs <- pairwise_auc(test, TN_class, score_prefix = "TNscore_")
auc_tn_pairs


TNmetrics


#combine it into a nice table
rho_tp.comb <- rho_tp %>%
  mutate(model_id = str_remove(score, "score"))

rho_tn.comb <- rho_tn %>%
  mutate(model_id = str_remove(score, "score"))

auc_tn_pairs.comb <- auc_tn_pairs %>%   
  mutate(model_id = str_remove(score, "score")) %>%
  filter(contrast == "Low vs High")

auc_tp_pairs.comb <- auc_tp_pairs %>%   
  mutate(model_id = str_remove(score, "score")) %>%
  filter(contrast == "Low vs High")


tp.allmetrics <- TPmetrics %>% 
  left_join(rho_tp.comb, by="model_id") %>% 
  left_join(auc_tp_pairs.comb, by="model_id")

tn.allmetrics <- TNmetrics %>% 
  left_join(rho_tn.comb, by="model_id") %>% 
  left_join(auc_tn_pairs.comb, by="model_id")
  

combined.nutrient.metrics <- bind_rows(tp.allmetrics,tn.allmetrics)


?join_by

summary(test)

library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)

TNinitial_gh_results <- test %>%
  games_howell_test(TNscore_base ~ TN_class) %>%
 # adjust_pvalue(method = "holm") %>%
  #filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TN_class", step.increase = 0.5) # auto-staggers brackets

oneway.test(TNscore_base ~ TN_class, data=test)

TNinitialbox <- ggplot(test, aes(x = TN_class, y = TNscore_base, color = TN_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Initial index", title = "(d)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TNinitial_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )

TNtop_gh_results <- test %>%
  games_howell_test(TNscore_DI ~ TN_class) %>%
#  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TN_class", step.increase = 0.8) # auto-staggers brackets

oneway.test(TNscore_DI ~ TN_class, data=test)

TNtopbox <- ggplot(test, aes(x = TN_class, y = TNscore_DI, color = TN_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "Nitrogen levels", y = "Top modified index", title = "(e)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TNtop_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )


TNuniversal_gh_results <- test %>%
  games_howell_test(TNscore_DNI ~ TN_class) %>%
#  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TN_class", step.increase = 0.8) # auto-staggers brackets

oneway.test(TNscore_DNI ~ TN_class, data=test)


TNuniversalbox <- ggplot(test, aes(x = TN_class, y = TNscore_DNI, color = TN_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Suggested modified index", title = "(f)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TNuniversal_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )



TPinitial_gh_results <- test %>%
  games_howell_test(TPscore_base ~ TP_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TP_class",step.increase = 0.5) # auto-staggers brackets

oneway.test(TPscore_base~ TP_class, data=test)

TPinitialbox <- ggplot(test, aes(x = TP_class, y = TPscore_base, color = TP_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Initial index", title = "(a)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TPinitial_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )

TPtop_gh_results <- test %>%
  games_howell_test(TPscore_DI ~ TP_class) %>%
 # adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TP_class", step.increase = 0.5) # auto-staggers brackets

oneway.test(TPscore_DI~ TP_class, data=test)


TPtopbox <- ggplot(test, aes(x = TP_class, y = TPscore_DI, color = TP_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "Phosphorus levels", y = "Top modified index", title = "(b)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TPtop_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )


TPuniversal_gh_results <- test %>%
  games_howell_test(TPscore_DNI ~ TP_class) %>%
 # adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TP_class", step.increase = 0.5) # auto-staggers brackets

oneway.test(TPscore_DNI ~ TP_class, data=test)

TPuniversalbox <- ggplot(test, aes(x = TP_class, y = TPscore_DNI, color = TP_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Suggested modified index", title = "(c)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TPuniversal_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )

grid.arrange(TPinitialbox,TPtopbox,TPuniversalbox,
             TNinitialbox,TNtopbox,TNuniversalbox,ncol=3)

allTNgames <- rbind(TNuniversal_gh_results,
                    TNtop_gh_results,
                    TNinitial_gh_results)

allTPgames<- rbind(TPuniversal_gh_results,
                   TPtop_gh_results,
                   TPinitial_gh_results)
#with HMSC ####
#Correct counts based on assembly mechanisms
{allindexscores$logcount.D <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)
allindexscores$logcount.I <- allindexscores$logcount*ifelse(allindexscores$n_neg>=CooccurCut,0,1)
allindexscores$logcount.N<- allindexscores$logcount*ifelse(allindexscores$NicheBreadth>=10,0,1)
allindexscores$logcount.DI <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)
allindexscores$logcount.DN <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
allindexscores$logcount.IN <- allindexscores$logcount*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
allindexscores$logcount.DNI <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
allindexscores$logcount.D.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)
allindexscores$logcount.I.5 <- allindexscores$logcount*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)
allindexscores$logcount.N.5<- allindexscores$logcount*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
allindexscores$logcount.D.5I <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)
allindexscores$logcount.DI.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)
allindexscores$logcount.D.5N <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
allindexscores$logcount.DN.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
allindexscores$logcount.I.5N <- allindexscores$logcount*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
allindexscores$logcount.IN.5 <- allindexscores$logcount*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
allindexscores$logcount.D.5NI <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
allindexscores$logcount.DN.5I <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
allindexscores$logcount.DNI.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
allindexscores$logcount.D.5N.5I <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
allindexscores$logcount.DN.5I.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
allindexscores$logcount.D.5NI.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
allindexscores$logcount.D.5N.5I.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
allindexscores$logcount.D.5I.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)
allindexscores$logcount.D.5N.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
allindexscores$logcount.I.5N.5 <- allindexscores$logcount*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)}

#TP site scores
{allindexscores$TPscore.D <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)
  allindexscores$TPscore.I <- allindexscores$TPscore*ifelse(allindexscores$n_neg>=CooccurCut,0,1)
  allindexscores$TPscore.N<- allindexscores$TPscore*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DI <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)
  allindexscores$TPscore.DN <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.IN <- allindexscores$TPscore*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DNI <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.D.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)
  allindexscores$TPscore.I.5 <- allindexscores$TPscore*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)
  allindexscores$TPscore.N.5<- allindexscores$TPscore*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5I <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)
  allindexscores$TPscore.DI.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)
  allindexscores$TPscore.D.5N <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DN.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.I.5N <- allindexscores$TPscore*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.IN.5 <- allindexscores$TPscore*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5NI <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DN.5I <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DNI.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5N.5I <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DN.5I.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5NI.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5N.5I.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5I.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)
  allindexscores$TPscore.D.5N.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.I.5N.5 <- allindexscores$TPscore*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
}

#TN site scores
{allindexscores$TNscore.D <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)
  allindexscores$TNscore.I <- allindexscores$TNscore*ifelse(allindexscores$n_neg>=CooccurCut,0,1)
  allindexscores$TNscore.N<- allindexscores$TNscore*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DI <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)
  allindexscores$TNscore.DN <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.IN <- allindexscores$TNscore*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DNI <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.D.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)
  allindexscores$TNscore.I.5 <- allindexscores$TNscore*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)
  allindexscores$TNscore.N.5<- allindexscores$TNscore*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5I <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)
  allindexscores$TNscore.DI.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)
  allindexscores$TNscore.D.5N <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DN.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.I.5N <- allindexscores$TNscore*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.IN.5 <- allindexscores$TNscore*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5NI <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DN.5I <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DNI.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5N.5I <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DN.5I.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5NI.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5N.5I.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5I.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)
  allindexscores$TNscore.D.5N.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.I.5N.5 <- allindexscores$TNscore*ifelse(allindexscores$n_neg>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
}

#Summarize measured values, dates, etc. 

testdata1 <- allindexscores[,c(1,6,7,12:21)] %>% group_by(UID) %>% summarise(across(.fns=mean, na.rm=T))
#unique_data <- distinct(allindexscores[,1:2])
#testdata1 <- left_join(testdata1,unique_data, by="UID")
testdata2 <- allindexscores[,c(1,34:ncol(allindexscores))] %>% group_by(UID) %>% summarise(across(.fns=sum, na.rm=T))

#View(testfinal)
testfinal <- left_join(testdata1,testdata2,by="UID")

#TP final score
{testfinal$TPscore <- testfinal$TPscore/testfinal$logcount
  testfinal$TPscore.D <- testfinal$TPscore.D/testfinal$logcount.D
  testfinal$TPscore.I <- testfinal$TPscore.I/testfinal$logcount.I
  testfinal$TPscore.N<- testfinal$TPscore.N/testfinal$logcount.N
  testfinal$TPscore.DI <- testfinal$TPscore.DI/testfinal$logcount.DI
  testfinal$TPscore.DN <- testfinal$TPscore.DN/testfinal$logcount.DN
  testfinal$TPscore.IN <- testfinal$TPscore.IN/testfinal$logcount.IN
  testfinal$TPscore.DNI <- testfinal$TPscore.DNI/testfinal$logcount.DNI
  testfinal$TPscore.D.5 <- testfinal$TPscore.D.5/testfinal$logcount.D.5
  testfinal$TPscore.I.5 <- testfinal$TPscore.I.5/testfinal$logcount.I.5
  testfinal$TPscore.N.5<- testfinal$TPscore.N.5/testfinal$logcount.N.5
  testfinal$TPscore.D.5I <- testfinal$TPscore.D.5I/testfinal$logcount.D.5I
  testfinal$TPscore.DI.5 <- testfinal$TPscore.DI.5/testfinal$logcount.DI.5
  testfinal$TPscore.D.5N <- testfinal$TPscore.D.5N/testfinal$logcount.D.5N
  testfinal$TPscore.DN.5 <- testfinal$TPscore.DN.5/testfinal$logcount.DN.5
  testfinal$TPscore.I.5N <- testfinal$TPscore.I.5N/testfinal$logcount.I.5N
  testfinal$TPscore.IN.5 <- testfinal$TPscore.IN.5/testfinal$logcount.IN.5
  testfinal$TPscore.D.5NI <- testfinal$TPscore.D.5NI/testfinal$logcount.D.5NI
  testfinal$TPscore.DN.5I <- testfinal$TPscore.DN.5I/testfinal$logcount.DN.5I
  testfinal$TPscore.DNI.5 <- testfinal$TPscore.DNI.5/testfinal$logcount.DNI.5
  testfinal$TPscore.D.5N.5I <- testfinal$TPscore.D.5N.5I/testfinal$logcount.D.5N.5I
  testfinal$TPscore.DN.5I.5 <- testfinal$TPscore.DN.5I.5/testfinal$logcount.DN.5I.5
  testfinal$TPscore.D.5NI.5 <- testfinal$TPscore.D.5NI.5/testfinal$logcount.D.5NI.5
  testfinal$TPscore.D.5N.5I.5 <- testfinal$TPscore.D.5N.5I.5/testfinal$logcount.D.5N.5I.5
  testfinal$TPscore.D.5I.5 <- testfinal$TPscore.D.5I.5/testfinal$logcount.D.5I.5
  testfinal$TPscore.D.5N.5 <- testfinal$TPscore.D.5N.5/testfinal$logcount.D.5N.5
  testfinal$TPscore.I.5N.5 <- testfinal$TPscore.I.5N.5/testfinal$logcount.I.5N.5
}

#TN final score
{
  testfinal$TNscore <- testfinal$TNscore/testfinal$logcount
  testfinal$TNscore.D <- testfinal$TNscore.D/testfinal$logcount.D
  testfinal$TNscore.I <- testfinal$TNscore.I/testfinal$logcount.I
  testfinal$TNscore.N<- testfinal$TNscore.N/testfinal$logcount.N
  testfinal$TNscore.DI <- testfinal$TNscore.DI/testfinal$logcount.DI
  testfinal$TNscore.DN <- testfinal$TNscore.DN/testfinal$logcount.DN
  testfinal$TNscore.IN <- testfinal$TNscore.IN/testfinal$logcount.IN
  testfinal$TNscore.DNI <- testfinal$TNscore.DNI/testfinal$logcount.DNI
  testfinal$TNscore.D.5 <- testfinal$TNscore.D.5/testfinal$logcount.D.5
  testfinal$TNscore.I.5 <- testfinal$TNscore.I.5/testfinal$logcount.I.5
  testfinal$TNscore.N.5<- testfinal$TNscore.N.5/testfinal$logcount.N.5
  testfinal$TNscore.D.5I <- testfinal$TNscore.D.5I/testfinal$logcount.D.5I
  testfinal$TNscore.DI.5 <- testfinal$TNscore.DI.5/testfinal$logcount.DI.5
  testfinal$TNscore.D.5N <- testfinal$TNscore.D.5N/testfinal$logcount.D.5N
  testfinal$TNscore.DN.5 <- testfinal$TNscore.DN.5/testfinal$logcount.DN.5
  testfinal$TNscore.I.5N <- testfinal$TNscore.I.5N/testfinal$logcount.I.5N
  testfinal$TNscore.IN.5 <- testfinal$TNscore.IN.5/testfinal$logcount.IN.5
  testfinal$TNscore.D.5NI <- testfinal$TNscore.D.5NI/testfinal$logcount.D.5NI
  testfinal$TNscore.DN.5I <- testfinal$TNscore.DN.5I/testfinal$logcount.DN.5I
  testfinal$TNscore.DNI.5 <- testfinal$TNscore.DNI.5/testfinal$logcount.DNI.5
  testfinal$TNscore.D.5N.5I <- testfinal$TNscore.D.5N.5I/testfinal$logcount.D.5N.5I
  testfinal$TNscore.DN.5I.5 <- testfinal$TNscore.DN.5I.5/testfinal$logcount.DN.5I.5
  testfinal$TNscore.D.5NI.5 <- testfinal$TNscore.D.5NI.5/testfinal$logcount.D.5NI.5
  testfinal$TNscore.D.5N.5I.5 <- testfinal$TNscore.D.5N.5I.5/testfinal$logcount.D.5N.5I.5
  testfinal$TNscore.D.5I.5 <- testfinal$TNscore.D.5I.5/testfinal$logcount.D.5I.5
  testfinal$TNscore.D.5N.5 <- testfinal$TNscore.D.5N.5/testfinal$logcount.D.5N.5
  testfinal$TNscore.I.5N.5 <- testfinal$TNscore.I.5N.5/testfinal$logcount.I.5N.5
}

#last step: remove sites without a score
test <- testfinal %>% filter(TNscore != "NaN")

#there is an issue
#test$Ecoregion.x <- stringr::str_to_title(test$Ecoregion.x)

# test <- test %>% mutate(region = case_when(
#   Ecoregion.x %in% c("Arkansas Valley", "Boston Mountains", "Ouachita Mountains","Ozark Highlands") ~ "Temperate Forest",
#   Ecoregion.x  %in% c("Central Great Plains", "Central Irregular Plains", "Cross Timbers") ~ "Central Plains",
#   Ecoregion.x  %in% c("South Central Plains", "Southwestern Tablelands") ~ "Southern Plains"))
# test <- test %>% mutate(season =case_when(Month %in% c(4,5) ~ "Spring",Month %in% c(6,7,8) ~"Summer",Month %in% c(9,10,11)~ "Fall"))
# 
# 
# spring <- test %>% filter(season=="Spring")
# summer <- test %>% filter(season=="Summer")
# fall <- test %>% filter(season=="Fall")
# =
# south <- test %>% filter(region=="Southern Plains")
# central <- test %>% filter(region=="Central Plains")
# forest <- test %>% filter(region=="Temperate Forest")
# 
# cor(spring$TN,spring$TNscore, use = "na.or.complete")

#TP models
{
  TPModel.null <- lm(TP~1,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore','TPscore.I')]))
  TPModel <-lm(TP~1+TPscore ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore','TPscore.I')]))
  TPModel.D <-lm(TP~1+TPscore.D ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D','TPscore.I')]))
  TPModel.I <-lm(TP~1+TPscore.I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I')]))
  TPModel.N<-lm(TP~1+TPscore.N,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.N','TPscore.I')]))
  TPModel.DI <-lm(TP~1+TPscore.DI ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DI')]))
  TPModel.DN <-lm(TP~1+TPscore.DN ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN','TPscore.I')]))
  TPModel.IN <-lm(TP~1+TPscore.IN ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.IN')]))
  TPModel.DNI <-lm(TP~1+TPscore.DNI ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DNI')]))
  TPModel.D.5 <-lm(TP~1+TPscore.D.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5','TPscore.I')]))
  TPModel.I.5 <-lm(TP~1+TPscore.I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I.5','TPscore.I')]))
  TPModel.N.5<-lm(TP~1+TPscore.N.5,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.N.5','TPscore.I')]))
  TPModel.D.5I <-lm(TP~1+TPscore.D.5I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5I')]))
  TPModel.DI.5 <-lm(TP~1+TPscore.DI.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DI.5')]))
  TPModel.D.5N <-lm(TP~1+TPscore.D.5N ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N','TPscore.I')]))
  TPModel.DN.5 <-lm(TP~1+TPscore.DN.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN.5','TPscore.I')]))
  TPModel.I.5N <-lm(TP~1+TPscore.I.5N ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I.5N')]))
  TPModel.IN.5 <-lm(TP~1+TPscore.IN.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.IN.5')]))
  TPModel.D.5NI <-lm(TP~1+TPscore.D.5NI ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5NI')]))
  TPModel.DN.5I <-lm(TP~1+TPscore.DN.5I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN.5I')]))
  TPModel.DNI.5 <-lm(TP~1+TPscore.DNI.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DNI.5')]))
  TPModel.D.5N.5I <-lm(TP~1+TPscore.D.5N.5I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N.5I')]))
  TPModel.DN.5I.5 <-lm(TP~1+TPscore.DN.5I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN.5I.5')]))
  TPModel.D.5NI.5 <-lm(TP~1+TPscore.D.5NI.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5NI.5')]))
  TPModel.D.5N.5I.5 <-lm(TP~1+TPscore.D.5N.5I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N.5I.5')]))
  TPModel.D.5I.5 <-lm(TP~1+TPscore.D.5I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5I.5')]))
  TPModel.D.5N.5 <-lm(TP~1+TPscore.D.5N.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N.5','TPscore.I')]))
  TPModel.I.5N.5 <-lm(TP~1+TPscore.I.5N.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I.5N.5')]))
  TPAICBIC <- cbind(AIC(TPModel,TPModel.D ,TPModel.I ,TPModel.N,TPModel.DI ,TPModel.DN ,TPModel.IN ,
                        TPModel.DNI ,TPModel.D.5 ,TPModel.I.5 ,TPModel.N.5,TPModel.D.5I ,TPModel.DI.5 ,
                        TPModel.D.5N ,TPModel.DN.5 ,TPModel.I.5N ,TPModel.IN.5 ,TPModel.D.5NI ,
                        TPModel.DN.5I ,TPModel.DNI.5 ,TPModel.D.5N.5I ,TPModel.DN.5I.5 ,TPModel.D.5NI.5 ,
                        TPModel.D.5N.5I.5 ,TPModel.D.5I.5 ,TPModel.D.5N.5 ,TPModel.I.5N.5 ,TPModel.null),
                    BIC(TPModel,TPModel.D ,TPModel.I ,TPModel.N,TPModel.DI ,TPModel.DN ,TPModel.IN ,
                        TPModel.DNI ,TPModel.D.5 ,TPModel.I.5 ,TPModel.N.5,TPModel.D.5I ,TPModel.DI.5 ,
                        TPModel.D.5N ,TPModel.DN.5 ,TPModel.I.5N ,TPModel.IN.5 ,TPModel.D.5NI ,
                        TPModel.DN.5I ,TPModel.DNI.5 ,TPModel.D.5N.5I ,TPModel.DN.5I.5 ,TPModel.D.5NI.5 ,
                        TPModel.D.5N.5I.5 ,TPModel.D.5I.5 ,TPModel.D.5N.5 ,TPModel.I.5N.5 ,TPModel.null))}

#TN models
{TNModel.null <- lm(TN~1,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore','TNscore.I')]))
  TNModel <-lm(TN~1+TNscore ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore','TNscore.I')]))
  TNModel.D <-lm(TN~1+TNscore.D ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D','TNscore.I')]))
  TNModel.I <-lm(TN~1+TNscore.I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I')]))
  TNModel.N<-lm(TN~1+TNscore.N,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.N','TNscore.I')]))
  TNModel.DI <-lm(TN~1+TNscore.DI ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DI')]))
  TNModel.DN <-lm(TN~1+TNscore.DN ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN','TNscore.I')]))
  TNModel.IN <-lm(TN~1+TNscore.IN ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.IN')]))
  TNModel.DNI <-lm(TN~1+TNscore.DNI ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DNI')]))
  TNModel.D.5 <-lm(TN~1+TNscore.D.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5','TNscore.I')]))
  TNModel.I.5 <-lm(TN~1+TNscore.I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I.5','TNscore.I')]))
  TNModel.N.5<-lm(TN~1+TNscore.N.5,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.N.5','TNscore.I')]))
  TNModel.D.5I <-lm(TN~1+TNscore.D.5I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5I')]))
  TNModel.DI.5 <-lm(TN~1+TNscore.DI.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DI.5')]))
  TNModel.D.5N <-lm(TN~1+TNscore.D.5N ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N','TNscore.I')]))
  TNModel.DN.5 <-lm(TN~1+TNscore.DN.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN.5','TNscore.I')]))
  TNModel.I.5N <-lm(TN~1+TNscore.I.5N ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I.5N')]))
  TNModel.IN.5 <-lm(TN~1+TNscore.IN.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.IN.5')]))
  TNModel.D.5NI <-lm(TN~1+TNscore.D.5NI ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5NI')]))
  TNModel.DN.5I <-lm(TN~1+TNscore.DN.5I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN.5I')]))
  TNModel.DNI.5 <-lm(TN~1+TNscore.DNI.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DNI.5')]))
  TNModel.D.5N.5I <-lm(TN~1+TNscore.D.5N.5I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N.5I')]))
  TNModel.DN.5I.5 <-lm(TN~1+TNscore.DN.5I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN.5I.5')]))
  TNModel.D.5NI.5 <-lm(TN~1+TNscore.D.5NI.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5NI.5')]))
  TNModel.D.5N.5I.5 <-lm(TN~1+TNscore.D.5N.5I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N.5I.5')]))
  TNModel.D.5I.5 <-lm(TN~1+TNscore.D.5I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5I.5')]))
  TNModel.D.5N.5 <-lm(TN~1+TNscore.D.5N.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N.5','TNscore.I')]))
  TNModel.I.5N.5 <-lm(TN~1+TNscore.I.5N.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I.5N.5')]))
  TNAICBIC <- cbind(AIC(TNModel,TNModel.D ,TNModel.I ,TNModel.N,TNModel.DI ,TNModel.DN ,TNModel.IN ,
                        TNModel.DNI ,TNModel.D.5 ,TNModel.I.5 ,TNModel.N.5,TNModel.D.5I ,TNModel.DI.5 ,
                        TNModel.D.5N ,TNModel.DN.5 ,TNModel.I.5N ,TNModel.IN.5 ,TNModel.D.5NI ,
                        TNModel.DN.5I ,TNModel.DNI.5 ,TNModel.D.5N.5I ,TNModel.DN.5I.5 ,TNModel.D.5NI.5 ,
                        TNModel.D.5N.5I.5 ,TNModel.D.5I.5 ,TNModel.D.5N.5 ,TNModel.I.5N.5 ,TNModel.null),
                    BIC(TNModel,TNModel.D ,TNModel.I ,TNModel.N,TNModel.DI ,TNModel.DN ,TNModel.IN ,
                        TNModel.DNI ,TNModel.D.5 ,TNModel.I.5 ,TNModel.N.5,TNModel.D.5I ,TNModel.DI.5 ,
                        TNModel.D.5N ,TNModel.DN.5 ,TNModel.I.5N ,TNModel.IN.5 ,TNModel.D.5NI ,
                        TNModel.DN.5I ,TNModel.DNI.5 ,TNModel.D.5N.5I ,TNModel.DN.5I.5 ,TNModel.D.5NI.5 ,
                        TNModel.D.5N.5I.5 ,TNModel.D.5I.5 ,TNModel.D.5N.5 ,TNModel.I.5N.5 ,TNModel.null))
  
  
}


View(TNAICBIC)
View(TPAICBIC)

TPAICBIC$model_id <- rownames(TPAICBIC)
TNAICBIC$model_id <- rownames(TNAICBIC)

TPAICBIC <- TPAICBIC[,c("AIC","BIC","model_id")]
TNAICBIC <- TNAICBIC[,c("AIC","BIC","model_id")]

## helper functions
mae  <- function(actual, pred) mean(abs(actual - pred), na.rm = TRUE)
rmse <- function(actual, pred) sqrt(mean((actual - pred)^2, na.rm = TRUE))

get_metrics <- function(mod) {
  mf   <- model.frame(mod)
  y    <- model.response(mf)
  yhat <- fitted(mod)
  s    <- summary(mod)
  
  data.frame(
    model_name = deparse1(formula(mod)),
    r2         = unname(s$r.squared),
    adj_r2     = unname(s$adj.r.squared),
    MAE        = mae(y, yhat),
    RMSE       = rmse(y, yhat)
  )
}

get_variable_metrics <- function(model) {
  summ <- summary(model)
  coefs <- summ$coefficients
  if (nrow(coefs) < 2) {
    return(data.frame(
      estimate = NA,
      std_error = NA,
      t_value = NA,
      p_value = NA
    ))
  } else {
    return(data.frame(
      estimate = coefs[2, 1],
      std_error = coefs[2, 2],
      t_value = coefs[2, 3],
      p_value = coefs[2, 4]
    ))
  }
}

## list of model names (strings only)
model_names <- c(
  "TPModel.null","TPModel","TPModel.D","TPModel.I","TPModel.N","TPModel.DI",
  "TPModel.DN","TPModel.IN","TPModel.DNI",
  "TPModel.D.5","TPModel.I.5","TPModel.N.5","TPModel.D.5I","TPModel.DI.5",
  "TPModel.D.5N","TPModel.DN.5","TPModel.I.5N","TPModel.IN.5",
  "TPModel.D.5NI","TPModel.DN.5I","TPModel.DNI.5",
  "TPModel.D.5N.5I","TPModel.DN.5I.5","TPModel.D.5NI.5","TPModel.D.5N.5I.5",
  "TPModel.D.5I.5","TPModel.D.5N.5","TPModel.I.5N.5"
)

## loop through them
TPmetrics <- do.call(rbind, lapply(model_names, function(nm) {
  mod <- get(nm)     # fetch object by name
  out <- get_metrics(mod)
  out2 <- get_variable_metrics(mod)
  out$model_id <- nm
  out3 <- cbind(out,out2)
  out3
}))

TPall <-  left_join(TPAICBIC,TPmetrics) 


## list of model names (strings only)
model_names <- c(
  "TNModel.null","TNModel","TNModel.D","TNModel.I","TNModel.N","TNModel.DI",
  "TNModel.DN","TNModel.IN","TNModel.DNI",
  "TNModel.D.5","TNModel.I.5","TNModel.N.5","TNModel.D.5I","TNModel.DI.5",
  "TNModel.D.5N","TNModel.DN.5","TNModel.I.5N","TNModel.IN.5",
  "TNModel.D.5NI","TNModel.DN.5I","TNModel.DNI.5",
  "TNModel.D.5N.5I","TNModel.DN.5I.5","TNModel.D.5NI.5","TNModel.D.5N.5I.5",
  "TNModel.D.5I.5","TNModel.D.5N.5","TNModel.I.5N.5"
)

## loop through them
TNmetrics <- do.call(rbind, lapply(model_names, function(nm) {
  mod <- get(nm)     # fetch object by name
  out <- get_metrics(mod)
  out2 <- get_variable_metrics(mod)
  out$model_id <- nm
  out3 <- cbind(out,out2)
  out3
}))

TNall <-  left_join(TNAICBIC,TNmetrics) 


##### Plots #####

library(ggeffects)

TPmodel.initial <- ggpredict(tp_models$TP_base,terms = "TPscore_base [0:10]")
TPmodel.corrected <-ggpredict(tp_models$TP_DI,terms = "TPscore_DI [0:10]")
TPmodel.universal <-ggpredict(tp_models$TP_DNI,terms = "TPscore_DNI [0:10]")
initialTP <- ggplot()+  geom_point(data=test,
                                   aes(x=TPscore_base,y=TP, color="red3"),alpha = 0.6)+
  geom_ribbon(data=TPmodel.initial, aes(x=x,ymin = conf.low,ymax=conf.high),fill = "red3",alpha=0.3)+
  stat_smooth(data=TPmodel.initial, aes(x=x,y=predicted),
              se=F,method = "glm",color = "red3", fill="red3")+ 
  xlab("Initial index")+ylab(expression("Log"["10"]*" transformed total phosphorus ( "*mu*"g/L)"))+
  
  # geom_smooth(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),method="glm",
  #             se=F,aes(x=TPscore,y=TP),size=0.1,linetype="dashed")+
  # geom_smooth(data=TPmodel.initial, se=F,aes(x=x,y=predicted),color="red3")+
  ylim(-0.5,4)+
  labs(title="(a)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()


correctedTP <- ggplot()+  geom_point(data=test,
                                     aes(x=TPscore_DI,y=TP), color="green3",alpha = 0.6)+
  geom_ribbon(data=TPmodel.corrected ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "green3",alpha=0.3)+
  stat_smooth(data=TPmodel.corrected, aes(x=x,y=predicted),
              se=F,method = "glm",color = "green3", fill="green3")+ 
  xlab("Top modified index")+ylab('')+
  ylim(-0.5,4)+
  xlim(0,10)+labs(title="(b)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()


universalTP <- ggplot()+  geom_point(data=test,
                                     aes(x=TPscore_DNI,y=TP), color="black",alpha = 0.6)+
  geom_ribbon(data=TPmodel.universal ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "black",alpha=0.3)+
  stat_smooth(data=TPmodel.universal, aes(x=x,y=predicted),
              se=F,method = "glm",color = "black", fill="black")+ 
  xlab("Suggested modified index")+ylab('')+
  ylim(-0.5,4)+
  xlim(0,10)+labs(title="(c)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()




gridExtra::grid.arrange(initialTP,correctedTP,universalTP,ncol=3)





TNmodel.initial <- ggpredict(tn_models$TN_base,terms = "TNscore_base [0:10]")
TNmodel.corrected <-ggpredict(tn_models$TN_DI,terms = "TNscore_DI [0:10]")
TNmodel.universal <-ggpredict(tn_models$TN_DNI,terms = "TNscore_DNI [0:10]")
initialTN <- ggplot()+  geom_point(data=test,
                                   aes(x=TNscore_base,y=TN, color="red3"),alpha = 0.6)+
  geom_ribbon(data=TNmodel.initial, aes(x=x,ymin = conf.low,ymax=conf.high),fill = "red3",alpha=0.3)+
  stat_smooth(data=TNmodel.initial, aes(x=x,y=predicted),
              se=F,method = "glm",color = "red3")+ 
  xlab("Initial index")+ylab(expression("Log"["10"]*" transformed total nitrogen (mg/L)"))+
  
  # geom_smooth(data=na.omit(test[,c("TNscore","TNscore.D.5","TN")]),method="glm",
  #             se=F,aes(x=TNscore,y=TN),size=0.1,linetype="dashed")+
  # geom_smooth(data=TNmodel.initial, se=F,aes(x=x,y=predicted),color="red3")+
  labs(title="(d)")+  ylim(-1.6,1.6)+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()

correctedTN <- ggplot()+  geom_point(data=test,
                                     aes(x=TNscore_DI,y=TN), color="green3",alpha = 0.6)+
  geom_ribbon(data=TNmodel.corrected ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "green3",alpha=0.3)+
  stat_smooth(data=TNmodel.corrected, aes(x=x,y=predicted),
              se=F,method = "glm",color = "green3")+ 
  xlab("Top modified index")+ylab('')+
  ylim(-1.6,1.6)+
  xlim(0,10)+labs(title="(e)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()


universalTN <- ggplot()+  geom_point(data=test,
                                     aes(x=TNscore_DNI,y=TN), color="black",alpha = 0.6)+
  geom_ribbon(data=TNmodel.universal ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "black",alpha=0.3)+
  stat_smooth(data=TNmodel.universal, aes(x=x,y=predicted),
              se=F,method = "glm",color = "black", fill="black")+ 
  xlab("Suggested modified index")+ylab('')+
  ylim(-1.6,1.6)+
  xlim(0,10)+labs(title="(f)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()




gridExtra::grid.arrange(initialTP,correctedTP,universalTP,
                        initialTN,correctedTN,universalTN,ncol=3)


test$TNactual <- 10^test$TN
test$TPactual <- 10^test$TP

test <- test %>%
  mutate(
    TN_class = case_when(
      TNactual <= 0.66 ~ "Low",
      TNactual <= 3.17 ~ "Medium",
      TNactual > 3.17       ~ "High"
    ),
    TP_class = case_when(
      TPactual <= 50   ~ "Low",
      TPactual <= 280  ~ "Medium",
      TPactual > 280     ~ "High"
    ),
    TN_class = factor(TN_class, levels = c("Low", "Medium", "High"), ordered = TRUE),
    TP_class = factor(TP_class, levels = c("Low", "Medium", "High"), ordered = TRUE)
  ) %>%
  filter(!is.na(TN_class) & !is.na(TP_class))


library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)

TNinitial_gh_results <- test %>%
  games_howell_test(TNscore ~ TN_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TN_class", step.increase = 0.25) # auto-staggers brackets

oneway.test(TNscore ~ TN_class, data=test)

TNinitialbox <- ggplot(test, aes(x = TN_class, y = TNscore, color = TN_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Initial index", title = "(d)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TNinitial_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )

TNtop_gh_results <- test %>%
  games_howell_test(TNscore.D ~ TN_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TN_class", step.increase = 0.25) # auto-staggers brackets

oneway.test(TNscore.D ~ TN_class, data=test)

TNtopbox <- ggplot(test, aes(x = TN_class, y = TNscore.D, color = TN_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Top modified index", title = "(e)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TNtop_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )


TNuniversal_gh_results <- test %>%
  games_howell_test(TNscore.D.5N.5I.5 ~ TN_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TN_class", step.increase = 0.25) # auto-staggers brackets

oneway.test(TNscore.D.5N.5I.5 ~ TN_class, data=test)


TNuniversalbox <- ggplot(test, aes(x = TN_class, y = TNscore.D.5N.5I.5, color = TN_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Suggested modified index", title = "(f)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TNuniversal_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )
TPinitial_gh_results <- test %>%
  games_howell_test(TPscore ~ TP_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TP_class",step.increase = 0.25) # auto-staggers brackets

oneway.test(TPscore~ TP_class, data=test)

TPinitialbox <- ggplot(test, aes(x = TP_class, y = TPscore, color = TP_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Initial index", title = "(a)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TPinitial_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )

TPtop_gh_results <- test %>%
  games_howell_test(TPscore.D.5 ~ TP_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TP_class", step.increase = 0.25) # auto-staggers brackets

oneway.test(TPscore.D.5~ TP_class, data=test)


TPtopbox <- ggplot(test, aes(x = TP_class, y = TPscore.D.5, color = TP_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Top modified index", title = "(b)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TPtop_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )


TPuniversal_gh_results <- test %>%
  games_howell_test(TPscore.D.5N.5I.5 ~ TP_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TP_class", step.increase = 0.25) # auto-staggers brackets

oneway.test(TPscore.D.5N.5I.5 ~ TP_class, data=test)

TPuniversalbox <- ggplot(test, aes(x = TP_class, y = TPscore.D.5N.5I.5, color = TP_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Suggested modified index", title = "(c)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TPuniversal_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )

grid.arrange(TPinitialbox,TPtopbox,TPuniversalbox,
             TNinitialbox,TNtopbox,TNuniversalbox,ncol=3)

allTNgames <- rbind(TNuniversal_gh_results,
                    TNtop_gh_results,
                    TNinitial_gh_results)

allTPgames<- rbind(TPuniversal_gh_results,
                   TPtop_gh_results,
                   TPinitial_gh_results)

#original cooccur ####

allindexscores$Interactions <- ifelse(is.na(allindexscores$Interactions),0,allindexscores$Interactions)
#Correct counts based on assembly mechanisms
{allindexscores$logcount.D <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)
  allindexscores$logcount.I <- allindexscores$logcount*ifelse(allindexscores$Interactions>=CooccurCut,0,1)
  allindexscores$logcount.N<- allindexscores$logcount*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DI <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)
  allindexscores$logcount.DN <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.IN <- allindexscores$logcount*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DNI <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.D.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)
  allindexscores$logcount.I.5 <- allindexscores$logcount*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)
  allindexscores$logcount.N.5<- allindexscores$logcount*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5I <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)
  allindexscores$logcount.DI.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)
  allindexscores$logcount.D.5N <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DN.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.I.5N <- allindexscores$logcount*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.IN.5 <- allindexscores$logcount*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5NI <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DN.5I <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DNI.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5N.5I <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$logcount.DN.5I.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5NI.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5N.5I.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.D.5I.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)
  allindexscores$logcount.D.5N.5 <- allindexscores$logcount*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$logcount.I.5N.5 <- allindexscores$logcount*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)}

#TP site scores
{allindexscores$TPscore.D <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)
  allindexscores$TPscore.I <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=CooccurCut,0,1)
  allindexscores$TPscore.N<- allindexscores$TPscore*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DI <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)
  allindexscores$TPscore.DN <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.IN <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DNI <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.D.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)
  allindexscores$TPscore.I.5 <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)
  allindexscores$TPscore.N.5<- allindexscores$TPscore*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5I <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)
  allindexscores$TPscore.DI.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)
  allindexscores$TPscore.D.5N <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DN.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.I.5N <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.IN.5 <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5NI <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DN.5I <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DNI.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5N.5I <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TPscore.DN.5I.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5NI.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5N.5I.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.D.5I.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)
  allindexscores$TPscore.D.5N.5 <- allindexscores$TPscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TPscore.I.5N.5 <- allindexscores$TPscore*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
}

#TN site scores
{allindexscores$TNscore.D <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)
  allindexscores$TNscore.I <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=CooccurCut,0,1)
  allindexscores$TNscore.N<- allindexscores$TNscore*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DI <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)
  allindexscores$TNscore.DN <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.IN <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DNI <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.D.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)
  allindexscores$TNscore.I.5 <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)
  allindexscores$TNscore.N.5<- allindexscores$TNscore*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5I <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)
  allindexscores$TNscore.DI.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)
  allindexscores$TNscore.D.5N <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DN.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.I.5N <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.IN.5 <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5NI <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DN.5I <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DNI.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5N.5I <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0,1)
  allindexscores$TNscore.DN.5I.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5NI.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5N.5I.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.D.5I.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)
  allindexscores$TNscore.D.5N.5 <- allindexscores$TNscore*ifelse(allindexscores$Dispersal>=DispersalCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
  allindexscores$TNscore.I.5N.5 <- allindexscores$TNscore*ifelse(allindexscores$Interactions>=CooccurCut,0.5,1)*ifelse(allindexscores$NicheBreadth>=10,0.5,1)
}

#Summarize measured values, dates, etc. 

testdata1 <- allindexscores[,c(1,6,7,12:21)] %>% group_by(UID) %>% summarise(across(.fns=mean, na.rm=T))
#unique_data <- distinct(allindexscores[,1:2])
#testdata1 <- left_join(testdata1,unique_data, by="UID")
testdata2 <- allindexscores[,c(1,34:ncol(allindexscores))] %>% group_by(UID) %>% summarise(across(.fns=sum, na.rm=T))

#View(testfinal)
testfinal <- left_join(testdata1,testdata2,by="UID")

#TP final score
{testfinal$TPscore <- testfinal$TPscore/testfinal$logcount
testfinal$TPscore.D <- testfinal$TPscore.D/testfinal$logcount.D
testfinal$TPscore.I <- testfinal$TPscore.I/testfinal$logcount.I
testfinal$TPscore.N<- testfinal$TPscore.N/testfinal$logcount.N
testfinal$TPscore.DI <- testfinal$TPscore.DI/testfinal$logcount.DI
testfinal$TPscore.DN <- testfinal$TPscore.DN/testfinal$logcount.DN
testfinal$TPscore.IN <- testfinal$TPscore.IN/testfinal$logcount.IN
testfinal$TPscore.DNI <- testfinal$TPscore.DNI/testfinal$logcount.DNI
testfinal$TPscore.D.5 <- testfinal$TPscore.D.5/testfinal$logcount.D.5
testfinal$TPscore.I.5 <- testfinal$TPscore.I.5/testfinal$logcount.I.5
testfinal$TPscore.N.5<- testfinal$TPscore.N.5/testfinal$logcount.N.5
testfinal$TPscore.D.5I <- testfinal$TPscore.D.5I/testfinal$logcount.D.5I
testfinal$TPscore.DI.5 <- testfinal$TPscore.DI.5/testfinal$logcount.DI.5
testfinal$TPscore.D.5N <- testfinal$TPscore.D.5N/testfinal$logcount.D.5N
testfinal$TPscore.DN.5 <- testfinal$TPscore.DN.5/testfinal$logcount.DN.5
testfinal$TPscore.I.5N <- testfinal$TPscore.I.5N/testfinal$logcount.I.5N
testfinal$TPscore.IN.5 <- testfinal$TPscore.IN.5/testfinal$logcount.IN.5
testfinal$TPscore.D.5NI <- testfinal$TPscore.D.5NI/testfinal$logcount.D.5NI
testfinal$TPscore.DN.5I <- testfinal$TPscore.DN.5I/testfinal$logcount.DN.5I
testfinal$TPscore.DNI.5 <- testfinal$TPscore.DNI.5/testfinal$logcount.DNI.5
testfinal$TPscore.D.5N.5I <- testfinal$TPscore.D.5N.5I/testfinal$logcount.D.5N.5I
testfinal$TPscore.DN.5I.5 <- testfinal$TPscore.DN.5I.5/testfinal$logcount.DN.5I.5
testfinal$TPscore.D.5NI.5 <- testfinal$TPscore.D.5NI.5/testfinal$logcount.D.5NI.5
testfinal$TPscore.D.5N.5I.5 <- testfinal$TPscore.D.5N.5I.5/testfinal$logcount.D.5N.5I.5
testfinal$TPscore.D.5I.5 <- testfinal$TPscore.D.5I.5/testfinal$logcount.D.5I.5
testfinal$TPscore.D.5N.5 <- testfinal$TPscore.D.5N.5/testfinal$logcount.D.5N.5
testfinal$TPscore.I.5N.5 <- testfinal$TPscore.I.5N.5/testfinal$logcount.I.5N.5
}

#TN final score
{
testfinal$TNscore <- testfinal$TNscore/testfinal$logcount
testfinal$TNscore.D <- testfinal$TNscore.D/testfinal$logcount.D
testfinal$TNscore.I <- testfinal$TNscore.I/testfinal$logcount.I
testfinal$TNscore.N<- testfinal$TNscore.N/testfinal$logcount.N
testfinal$TNscore.DI <- testfinal$TNscore.DI/testfinal$logcount.DI
testfinal$TNscore.DN <- testfinal$TNscore.DN/testfinal$logcount.DN
testfinal$TNscore.IN <- testfinal$TNscore.IN/testfinal$logcount.IN
testfinal$TNscore.DNI <- testfinal$TNscore.DNI/testfinal$logcount.DNI
testfinal$TNscore.D.5 <- testfinal$TNscore.D.5/testfinal$logcount.D.5
testfinal$TNscore.I.5 <- testfinal$TNscore.I.5/testfinal$logcount.I.5
testfinal$TNscore.N.5<- testfinal$TNscore.N.5/testfinal$logcount.N.5
testfinal$TNscore.D.5I <- testfinal$TNscore.D.5I/testfinal$logcount.D.5I
testfinal$TNscore.DI.5 <- testfinal$TNscore.DI.5/testfinal$logcount.DI.5
testfinal$TNscore.D.5N <- testfinal$TNscore.D.5N/testfinal$logcount.D.5N
testfinal$TNscore.DN.5 <- testfinal$TNscore.DN.5/testfinal$logcount.DN.5
testfinal$TNscore.I.5N <- testfinal$TNscore.I.5N/testfinal$logcount.I.5N
testfinal$TNscore.IN.5 <- testfinal$TNscore.IN.5/testfinal$logcount.IN.5
testfinal$TNscore.D.5NI <- testfinal$TNscore.D.5NI/testfinal$logcount.D.5NI
testfinal$TNscore.DN.5I <- testfinal$TNscore.DN.5I/testfinal$logcount.DN.5I
testfinal$TNscore.DNI.5 <- testfinal$TNscore.DNI.5/testfinal$logcount.DNI.5
testfinal$TNscore.D.5N.5I <- testfinal$TNscore.D.5N.5I/testfinal$logcount.D.5N.5I
testfinal$TNscore.DN.5I.5 <- testfinal$TNscore.DN.5I.5/testfinal$logcount.DN.5I.5
testfinal$TNscore.D.5NI.5 <- testfinal$TNscore.D.5NI.5/testfinal$logcount.D.5NI.5
testfinal$TNscore.D.5N.5I.5 <- testfinal$TNscore.D.5N.5I.5/testfinal$logcount.D.5N.5I.5
testfinal$TNscore.D.5I.5 <- testfinal$TNscore.D.5I.5/testfinal$logcount.D.5I.5
testfinal$TNscore.D.5N.5 <- testfinal$TNscore.D.5N.5/testfinal$logcount.D.5N.5
testfinal$TNscore.I.5N.5 <- testfinal$TNscore.I.5N.5/testfinal$logcount.I.5N.5
}

#last step: remove sites without a score
test <- testfinal %>% filter(TNscore != "NaN")

test_old <- na.omit(test[,c('UID','TPscore.DNI','Year','TP','TPscore','TPscore.I')])

these_UIDs <- unique(test_old$UID)

#there is an issue
#test$Ecoregion.x <- stringr::str_to_title(test$Ecoregion.x)

# test <- test %>% mutate(region = case_when(
#   Ecoregion.x %in% c("Arkansas Valley", "Boston Mountains", "Ouachita Mountains","Ozark Highlands") ~ "Temperate Forest",
#   Ecoregion.x  %in% c("Central Great Plains", "Central Irregular Plains", "Cross Timbers") ~ "Central Plains",
#   Ecoregion.x  %in% c("South Central Plains", "Southwestern Tablelands") ~ "Southern Plains"))
# test <- test %>% mutate(season =case_when(Month %in% c(4,5) ~ "Spring",Month %in% c(6,7,8) ~"Summer",Month %in% c(9,10,11)~ "Fall"))
# 
# 
# spring <- test %>% filter(season=="Spring")
# summer <- test %>% filter(season=="Summer")
# fall <- test %>% filter(season=="Fall")
# =
# south <- test %>% filter(region=="Southern Plains")
# central <- test %>% filter(region=="Central Plains")
# forest <- test %>% filter(region=="Temperate Forest")
# 
# cor(spring$TN,spring$TNscore, use = "na.or.complete")

#TP models
{
  TPModel.null <- lm(TP~1,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore','TPscore.I')]))
  TPModel <-lm(TP~1+TPscore ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore','TPscore.I')]))
  TPModel.D <-lm(TP~1+TPscore.D ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D','TPscore.I')]))
  TPModel.I <-lm(TP~1+TPscore.I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I')]))
  TPModel.N<-lm(TP~1+TPscore.N,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.N','TPscore.I')]))
  TPModel.DI <-lm(TP~1+TPscore.DI ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DI')]))
  TPModel.DN <-lm(TP~1+TPscore.DN ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN','TPscore.I')]))
  TPModel.IN <-lm(TP~1+TPscore.IN ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.IN')]))
  TPModel.DNI <-lm(TP~1+TPscore.DNI ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DNI')]))
  TPModel.D.5 <-lm(TP~1+TPscore.D.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5','TPscore.I')]))
  TPModel.I.5 <-lm(TP~1+TPscore.I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I.5','TPscore.I')]))
  TPModel.N.5<-lm(TP~1+TPscore.N.5,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.N.5','TPscore.I')]))
  TPModel.D.5I <-lm(TP~1+TPscore.D.5I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5I')]))
  TPModel.DI.5 <-lm(TP~1+TPscore.DI.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DI.5')]))
  TPModel.D.5N <-lm(TP~1+TPscore.D.5N ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N','TPscore.I')]))
  TPModel.DN.5 <-lm(TP~1+TPscore.DN.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN.5','TPscore.I')]))
  TPModel.I.5N <-lm(TP~1+TPscore.I.5N ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I.5N')]))
  TPModel.IN.5 <-lm(TP~1+TPscore.IN.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.IN.5')]))
  TPModel.D.5NI <-lm(TP~1+TPscore.D.5NI ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5NI')]))
  TPModel.DN.5I <-lm(TP~1+TPscore.DN.5I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN.5I')]))
  TPModel.DNI.5 <-lm(TP~1+TPscore.DNI.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DNI.5')]))
  TPModel.D.5N.5I <-lm(TP~1+TPscore.D.5N.5I ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N.5I')]))
  TPModel.DN.5I.5 <-lm(TP~1+TPscore.DN.5I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.DN.5I.5')]))
  TPModel.D.5NI.5 <-lm(TP~1+TPscore.D.5NI.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5NI.5')]))
  TPModel.D.5N.5I.5 <-lm(TP~1+TPscore.D.5N.5I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N.5I.5')]))
  TPModel.D.5I.5 <-lm(TP~1+TPscore.D.5I.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5I.5')]))
  TPModel.D.5N.5 <-lm(TP~1+TPscore.D.5N.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5N.5','TPscore.I')]))
  TPModel.I.5N.5 <-lm(TP~1+TPscore.I.5N.5 ,data=na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.I.5N.5')]))
  TPAICBIC <- cbind(AIC(TPModel,TPModel.D ,TPModel.I ,TPModel.N,TPModel.DI ,TPModel.DN ,TPModel.IN ,
                        TPModel.DNI ,TPModel.D.5 ,TPModel.I.5 ,TPModel.N.5,TPModel.D.5I ,TPModel.DI.5 ,
                        TPModel.D.5N ,TPModel.DN.5 ,TPModel.I.5N ,TPModel.IN.5 ,TPModel.D.5NI ,
                        TPModel.DN.5I ,TPModel.DNI.5 ,TPModel.D.5N.5I ,TPModel.DN.5I.5 ,TPModel.D.5NI.5 ,
                        TPModel.D.5N.5I.5 ,TPModel.D.5I.5 ,TPModel.D.5N.5 ,TPModel.I.5N.5 ,TPModel.null),
                    BIC(TPModel,TPModel.D ,TPModel.I ,TPModel.N,TPModel.DI ,TPModel.DN ,TPModel.IN ,
                        TPModel.DNI ,TPModel.D.5 ,TPModel.I.5 ,TPModel.N.5,TPModel.D.5I ,TPModel.DI.5 ,
                        TPModel.D.5N ,TPModel.DN.5 ,TPModel.I.5N ,TPModel.IN.5 ,TPModel.D.5NI ,
                        TPModel.DN.5I ,TPModel.DNI.5 ,TPModel.D.5N.5I ,TPModel.DN.5I.5 ,TPModel.D.5NI.5 ,
                        TPModel.D.5N.5I.5 ,TPModel.D.5I.5 ,TPModel.D.5N.5 ,TPModel.I.5N.5 ,TPModel.null))}

#TN models
{TNModel.null <- lm(TN~1,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore','TNscore.I')]))
  TNModel <-lm(TN~1+TNscore ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore','TNscore.I')]))
  TNModel.D <-lm(TN~1+TNscore.D ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D','TNscore.I')]))
  TNModel.I <-lm(TN~1+TNscore.I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I')]))
  TNModel.N<-lm(TN~1+TNscore.N,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.N','TNscore.I')]))
  TNModel.DI <-lm(TN~1+TNscore.DI ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DI')]))
  TNModel.DN <-lm(TN~1+TNscore.DN ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN','TNscore.I')]))
  TNModel.IN <-lm(TN~1+TNscore.IN ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.IN')]))
  TNModel.DNI <-lm(TN~1+TNscore.DNI ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DNI')]))
  TNModel.D.5 <-lm(TN~1+TNscore.D.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5','TNscore.I')]))
  TNModel.I.5 <-lm(TN~1+TNscore.I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I.5','TNscore.I')]))
  TNModel.N.5<-lm(TN~1+TNscore.N.5,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.N.5','TNscore.I')]))
  TNModel.D.5I <-lm(TN~1+TNscore.D.5I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5I')]))
  TNModel.DI.5 <-lm(TN~1+TNscore.DI.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DI.5')]))
  TNModel.D.5N <-lm(TN~1+TNscore.D.5N ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N','TNscore.I')]))
  TNModel.DN.5 <-lm(TN~1+TNscore.DN.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN.5','TNscore.I')]))
  TNModel.I.5N <-lm(TN~1+TNscore.I.5N ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I.5N')]))
  TNModel.IN.5 <-lm(TN~1+TNscore.IN.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.IN.5')]))
  TNModel.D.5NI <-lm(TN~1+TNscore.D.5NI ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5NI')]))
  TNModel.DN.5I <-lm(TN~1+TNscore.DN.5I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN.5I')]))
  TNModel.DNI.5 <-lm(TN~1+TNscore.DNI.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DNI.5')]))
  TNModel.D.5N.5I <-lm(TN~1+TNscore.D.5N.5I ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N.5I')]))
  TNModel.DN.5I.5 <-lm(TN~1+TNscore.DN.5I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.DN.5I.5')]))
  TNModel.D.5NI.5 <-lm(TN~1+TNscore.D.5NI.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5NI.5')]))
  TNModel.D.5N.5I.5 <-lm(TN~1+TNscore.D.5N.5I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N.5I.5')]))
  TNModel.D.5I.5 <-lm(TN~1+TNscore.D.5I.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5I.5')]))
  TNModel.D.5N.5 <-lm(TN~1+TNscore.D.5N.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.D.5N.5','TNscore.I')]))
  TNModel.I.5N.5 <-lm(TN~1+TNscore.I.5N.5 ,data=na.omit(test[,c('TNscore.DNI','Year','TN','TNscore.I.5N.5')]))
  TNAICBIC <- cbind(AIC(TNModel,TNModel.D ,TNModel.I ,TNModel.N,TNModel.DI ,TNModel.DN ,TNModel.IN ,
                        TNModel.DNI ,TNModel.D.5 ,TNModel.I.5 ,TNModel.N.5,TNModel.D.5I ,TNModel.DI.5 ,
                        TNModel.D.5N ,TNModel.DN.5 ,TNModel.I.5N ,TNModel.IN.5 ,TNModel.D.5NI ,
                        TNModel.DN.5I ,TNModel.DNI.5 ,TNModel.D.5N.5I ,TNModel.DN.5I.5 ,TNModel.D.5NI.5 ,
                        TNModel.D.5N.5I.5 ,TNModel.D.5I.5 ,TNModel.D.5N.5 ,TNModel.I.5N.5 ,TNModel.null),
                    BIC(TNModel,TNModel.D ,TNModel.I ,TNModel.N,TNModel.DI ,TNModel.DN ,TNModel.IN ,
                        TNModel.DNI ,TNModel.D.5 ,TNModel.I.5 ,TNModel.N.5,TNModel.D.5I ,TNModel.DI.5 ,
                        TNModel.D.5N ,TNModel.DN.5 ,TNModel.I.5N ,TNModel.IN.5 ,TNModel.D.5NI ,
                        TNModel.DN.5I ,TNModel.DNI.5 ,TNModel.D.5N.5I ,TNModel.DN.5I.5 ,TNModel.D.5NI.5 ,
                        TNModel.D.5N.5I.5 ,TNModel.D.5I.5 ,TNModel.D.5N.5 ,TNModel.I.5N.5 ,TNModel.null))
  
  
  }


View(TNAICBIC)
View(TPAICBIC)

TPAICBIC$model_id <- rownames(TPAICBIC)
TNAICBIC$model_id <- rownames(TNAICBIC)

TPAICBIC <- TPAICBIC[,c("AIC","BIC","model_id")]
TNAICBIC <- TNAICBIC[,c("AIC","BIC","model_id")]

## helper functions
mae  <- function(actual, pred) mean(abs(actual - pred), na.rm = TRUE)
rmse <- function(actual, pred) sqrt(mean((actual - pred)^2, na.rm = TRUE))

get_metrics <- function(mod) {
  mf   <- model.frame(mod)
  y    <- model.response(mf)
  yhat <- fitted(mod)
  s    <- summary(mod)
  
  data.frame(
    model_name = deparse1(formula(mod)),
    r2         = unname(s$r.squared),
    adj_r2     = unname(s$adj.r.squared),
    MAE        = mae(y, yhat),
    RMSE       = rmse(y, yhat)
  )
}

get_variable_metrics <- function(model) {
  summ <- summary(model)
  coefs <- summ$coefficients
  if (nrow(coefs) < 2) {
    return(data.frame(
      estimate = NA,
      std_error = NA,
      t_value = NA,
      p_value = NA
    ))
  } else {
    return(data.frame(
      estimate = coefs[2, 1],
      std_error = coefs[2, 2],
      t_value = coefs[2, 3],
      p_value = coefs[2, 4]
    ))
  }
}

## list of model names (strings only)
model_names <- c(
  "TPModel.null","TPModel","TPModel.D","TPModel.I","TPModel.N","TPModel.DI",
  "TPModel.DN","TPModel.IN","TPModel.DNI",
  "TPModel.D.5","TPModel.I.5","TPModel.N.5","TPModel.D.5I","TPModel.DI.5",
  "TPModel.D.5N","TPModel.DN.5","TPModel.I.5N","TPModel.IN.5",
  "TPModel.D.5NI","TPModel.DN.5I","TPModel.DNI.5",
  "TPModel.D.5N.5I","TPModel.DN.5I.5","TPModel.D.5NI.5","TPModel.D.5N.5I.5",
  "TPModel.D.5I.5","TPModel.D.5N.5","TPModel.I.5N.5"
)

## loop through them
TPmetrics <- do.call(rbind, lapply(model_names, function(nm) {
  mod <- get(nm)     # fetch object by name
  out <- get_metrics(mod)
  out2 <- get_variable_metrics(mod)
  out$model_id <- nm
  out3 <- cbind(out,out2)
  out3
}))

TPall <-  left_join(TPAICBIC,TPmetrics) 


## list of model names (strings only)
model_names <- c(
  "TNModel.null","TNModel","TNModel.D","TNModel.I","TNModel.N","TNModel.DI",
  "TNModel.DN","TNModel.IN","TNModel.DNI",
  "TNModel.D.5","TNModel.I.5","TNModel.N.5","TNModel.D.5I","TNModel.DI.5",
  "TNModel.D.5N","TNModel.DN.5","TNModel.I.5N","TNModel.IN.5",
  "TNModel.D.5NI","TNModel.DN.5I","TNModel.DNI.5",
  "TNModel.D.5N.5I","TNModel.DN.5I.5","TNModel.D.5NI.5","TNModel.D.5N.5I.5",
  "TNModel.D.5I.5","TNModel.D.5N.5","TNModel.I.5N.5"
)

## loop through them
TNmetrics <- do.call(rbind, lapply(model_names, function(nm) {
  mod <- get(nm)     # fetch object by name
  out <- get_metrics(mod)
  out2 <- get_variable_metrics(mod)
  out$model_id <- nm
  out3 <- cbind(out,out2)
  out3
}))

TNall <-  left_join(TNAICBIC,TNmetrics) 


# #plot colors
# 
# region_colors <- c("Central Plains" = "orange", "Temperate Forest" = "forestgreen", "Southern Plains" = "red4")
# season_colors <- c("Fall" = "brown4", "Spring" = "forestgreen", "Summer" = "yellow3")
# 
# ??ggeffects::plot
# #Nitrogen Model Plots
# ?geom_smooth

library(ggeffects)

TPmodel.initial <- ggpredict(TPModel,terms = "TPscore [0:10]")
TPmodel.corrected <-ggpredict(TPModel.D.5,terms = "TPscore.D.5 [0:10]")
TPmodel.universal <-ggpredict(TPModel.D.5N.5I.5,terms = "TPscore.D.5N.5I.5 [0:10]")
initialTP <- ggplot()+  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),
                                 aes(x=TPscore,y=TP, color="red3"),alpha = 0.6)+
  geom_ribbon(data=TPmodel.initial, aes(x=x,ymin = conf.low,ymax=conf.high),fill = "red3",alpha=0.3)+
  stat_smooth(data=TPmodel.initial, aes(x=x,y=predicted),
                                se=F,method = "glm",color = "red3", fill="red3")+ 
  xlab("Initial index")+ylab(expression("Log"["10"]*" transformed total phosphorus ( "*mu*"g/L)"))+

  # geom_smooth(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),method="glm",
  #             se=F,aes(x=TPscore,y=TP),size=0.1,linetype="dashed")+
  # geom_smooth(data=TPmodel.initial, se=F,aes(x=x,y=predicted),color="red3")+
  ylim(-0.5,4)+
  labs(title="(a)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()


correctedTP <- ggplot()+  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TPscore.D.5N.5I.5","TP")]),
                                  aes(x=TPscore.D.5,y=TP), color="green3",alpha = 0.6)+
  geom_ribbon(data=TPmodel.corrected ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "green3",alpha=0.3)+
  stat_smooth(data=TPmodel.corrected, aes(x=x,y=predicted),
              se=F,method = "glm",color = "green3", fill="green3")+ 
  xlab("Top modified index")+ylab('')+
  ylim(-0.5,4)+
  xlim(0,10)+labs(title="(b)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()


universalTP <- ggplot()+  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TPscore.D.5N.5I.5","TP")]),
                                  aes(x=TPscore.D.5N.5I.5,y=TP), color="black",alpha = 0.6)+
  geom_ribbon(data=TPmodel.universal ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "black",alpha=0.3)+
  stat_smooth(data=TPmodel.universal, aes(x=x,y=predicted),
              se=F,method = "glm",color = "black", fill="black")+ 
  xlab("Suggested modified index")+ylab('')+
  ylim(-0.5,4)+
  xlim(0,10)+labs(title="(c)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()




#gridExtra::grid.arrange(initial,corrected,universal,ncol=3)





TNmodel.initial <- ggpredict(TNModel,terms = "TNscore [0:10]")
TNmodel.corrected <-ggpredict(TNModel.D,terms = "TNscore.D [0:10]")
TNmodel.universal <-ggpredict(TNModel.D.5N.5I.5,terms = "TNscore.D.5N.5I.5 [0:10]")
initialTN <- ggplot()+  geom_point(data=na.omit(test[,c("TNscore","TNscore.D.5","TN")]),
                                 aes(x=TNscore,y=TN, color="red3"),alpha = 0.6)+
  geom_ribbon(data=TNmodel.initial, aes(x=x,ymin = conf.low,ymax=conf.high),fill = "red3",alpha=0.3)+
  stat_smooth(data=TNmodel.initial, aes(x=x,y=predicted),
              se=F,method = "glm",color = "red3")+ 
  xlab("Initial index")+ylab(expression("Log"["10"]*" transformed total nitrogen (mg/L)"))+

  # geom_smooth(data=na.omit(test[,c("TNscore","TNscore.D.5","TN")]),method="glm",
  #             se=F,aes(x=TNscore,y=TN),size=0.1,linetype="dashed")+
  # geom_smooth(data=TNmodel.initial, se=F,aes(x=x,y=predicted),color="red3")+
  labs(title="(d)")+  ylim(-1.6,1.6)+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()

correctedTN <- ggplot()+  geom_point(data=na.omit(test[,c("TNscore","TNscore.D.5","TN")]),
                                  aes(x=TNscore.D.5,y=TN), color="green3",alpha = 0.6)+
  geom_ribbon(data=TNmodel.corrected ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "green3",alpha=0.3)+
  stat_smooth(data=TNmodel.corrected, aes(x=x,y=predicted),
              se=F,method = "glm",color = "green3")+ 
  xlab("Top modified index")+ylab('')+
  ylim(-1.6,1.6)+
  xlim(0,10)+labs(title="(e)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()


universalTN <- ggplot()+  geom_point(data=na.omit(test[,c("TNscore","TNscore.D.5","TNscore.D.5N.5I.5","TN")]),
                                   aes(x=TNscore.D.5N.5I.5,y=TN), color="black",alpha = 0.6)+
  geom_ribbon(data=TNmodel.universal ,aes(x=x,ymin = conf.low,ymax=conf.high),fill = "black",alpha=0.3)+
  stat_smooth(data=TNmodel.universal, aes(x=x,y=predicted),
              se=F,method = "glm",color = "black", fill="black")+ 
  xlab("Suggested modified index")+ylab('')+
  ylim(-1.6,1.6)+
  xlim(0,10)+labs(title="(f)")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_classic()




gridExtra::grid.arrange(initialTP,correctedTP,universalTP,
                        initialTN,correctedTN,universalTN,ncol=3)


test$TNactual <- 10^test$TN
test$TPactual <- 10^test$TP

test <- test %>%
  mutate(
    TN_class = case_when(
      TNactual <= 0.66 ~ "Low",
      TNactual <= 3.17 ~ "Medium",
      TNactual > 3.17       ~ "High"
    ),
    TP_class = case_when(
      TPactual <= 50   ~ "Low",
      TPactual <= 280  ~ "Medium",
      TPactual > 280     ~ "High"
    ),
    TN_class = factor(TN_class, levels = c("Low", "Medium", "High"), ordered = TRUE),
    TP_class = factor(TP_class, levels = c("Low", "Medium", "High"), ordered = TRUE)
  ) %>%
  filter(!is.na(TN_class) & !is.na(TP_class))


library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)

TNinitial_gh_results <- test %>%
  games_howell_test(TNscore ~ TN_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TN_class", step.increase = 0.25) # auto-staggers brackets

oneway.test(TNscore ~ TN_class, data=test)

TNinitialbox <- ggplot(test, aes(x = TN_class, y = TNscore, color = TN_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Initial index", title = "(d)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TNinitial_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )

TNtop_gh_results <- test %>%
  games_howell_test(TNscore.D ~ TN_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TN_class", step.increase = 0.25) # auto-staggers brackets

oneway.test(TNscore.D ~ TN_class, data=test)

TNtopbox <- ggplot(test, aes(x = TN_class, y = TNscore.D, color = TN_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Top modified index", title = "(e)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TNtop_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )


TNuniversal_gh_results <- test %>%
  games_howell_test(TNscore.D.5N.5I.5 ~ TN_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "TN_class", step.increase = 0.25) # auto-staggers brackets

oneway.test(TNscore.D.5N.5I.5 ~ TN_class, data=test)


TNuniversalbox <- ggplot(test, aes(x = TN_class, y = TNscore.D.5N.5I.5, color = TN_class)) +
  geom_jitter() +
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  ylim(0, 10) +
  labs(x = "", y = "Suggested modified index", title = "(f)") +
  scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    TNuniversal_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )
 TPinitial_gh_results <- test %>%
   games_howell_test(TPscore ~ TP_class) %>%
   adjust_pvalue(method = "holm") %>%
   filter(p.adj < 0.05) %>%
   mutate(p.adj.signif = "") %>%   # make empty labels
   add_xy_position(x = "TP_class",step.increase = 0.25) # auto-staggers brackets
 
 oneway.test(TPscore~ TP_class, data=test)
 
 TPinitialbox <- ggplot(test, aes(x = TP_class, y = TPscore, color = TP_class)) +
   geom_jitter() +
   geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
   theme_classic() +
   ylim(0, 10) +
   labs(x = "", y = "Initial index", title = "(a)") +
   scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
   theme(legend.position = "none", text = element_text(size = 12)) +
   stat_pvalue_manual(
     TPinitial_gh_results,
     label = "p.adj.signif",   # "" → brackets only
     hide.ns = TRUE,
     tip.length = 0.01,
     size = 0.6,               # bracket line thickness
     bracket.size = 0.8        # widen/thicken brackets
   )
 
 TPtop_gh_results <- test %>%
   games_howell_test(TPscore.D.5 ~ TP_class) %>%
   adjust_pvalue(method = "holm") %>%
   filter(p.adj < 0.05) %>%
   mutate(p.adj.signif = "") %>%   # make empty labels
   add_xy_position(x = "TP_class", step.increase = 0.25) # auto-staggers brackets
 
 oneway.test(TPscore.D.5~ TP_class, data=test)
 

 TPtopbox <- ggplot(test, aes(x = TP_class, y = TPscore.D.5, color = TP_class)) +
   geom_jitter() +
   geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
   theme_classic() +
   ylim(0, 10) +
   labs(x = "", y = "Top modified index", title = "(b)") +
   scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
   theme(legend.position = "none", text = element_text(size = 12)) +
   stat_pvalue_manual(
     TPtop_gh_results,
     label = "p.adj.signif",   # "" → brackets only
     hide.ns = TRUE,
     tip.length = 0.01,
     size = 0.6,               # bracket line thickness
     bracket.size = 0.8        # widen/thicken brackets
   )
 
 
 TPuniversal_gh_results <- test %>%
   games_howell_test(TPscore.D.5N.5I.5 ~ TP_class) %>%
   adjust_pvalue(method = "holm") %>%
   filter(p.adj < 0.05) %>%
   mutate(p.adj.signif = "") %>%   # make empty labels
   add_xy_position(x = "TP_class", step.increase = 0.25) # auto-staggers brackets
 
 oneway.test(TPscore.D.5N.5I.5 ~ TP_class, data=test)
 
 TPuniversalbox <- ggplot(test, aes(x = TP_class, y = TPscore.D.5N.5I.5, color = TP_class)) +
   geom_jitter() +
   geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
   theme_classic() +
   ylim(0, 10) +
   labs(x = "", y = "Suggested modified index", title = "(c)") +
   scale_color_manual(values = c("lightgreen", "green", "darkgreen")) +
   theme(legend.position = "none", text = element_text(size = 12)) +
   stat_pvalue_manual(
     TPuniversal_gh_results,
     label = "p.adj.signif",   # "" → brackets only
     hide.ns = TRUE,
     tip.length = 0.01,
     size = 0.6,               # bracket line thickness
     bracket.size = 0.8        # widen/thicken brackets
   )
 
 grid.arrange(TPinitialbox,TPtopbox,TPuniversalbox,
              TNinitialbox,TNtopbox,TNuniversalbox,ncol=3)
 
 allTNgames <- rbind(TNuniversal_gh_results,
                     TNtop_gh_results,
                     TNinitial_gh_results)
 
 allTPgames<- rbind(TPuniversal_gh_results,
                    TPtop_gh_results,
                    TPinitial_gh_results)
#####
initial2 <- plot(tnmodel.initial,line_size = 1,limit_range = F)+ 
  xlab("Initial index")+ylab("Log10-transformed total nitrogen (mg/L)")+
  geom_point(data=na.omit(test[,c("TNscore","TNscore.D","TN")]),
             aes(x=TNscore,y=TN,group=region,color=season),alpha = 0.7)+
  scale_color_manual(values = season_colors)+xlim(0,10)+labs(title="",color="Season")+scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
corrected2<- plot(tnmodel.corrected,line_size = 1,limit_range = F)+ 
  xlab("Top performing index")+ylab("")+
  geom_point(data=na.omit(test[,c("TNscore","TNscore.D","TN","region","season")]),
             aes(x=TNscore.D,y=TN,group=region,color=season),alpha = 0.7)+
  scale_color_manual(values = season_colors)+xlim(0,10)+labs(title="")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
gridExtra::grid.arrange(initial,corrected,ncol=2)
gridExtra::grid.arrange(initial2,corrected2,ncol=2)


TPmodel.initial <- ggpredict(TPModel,terms = "TPscore")
TPmodel.corrected <-ggpredict(TPModel.D.5,terms = "TPscore.D.5")
initial <- plot(TPmodel.initial,line_size = 1,limit_range = T,colors = "red3")+ 
  xlab("Initial index score")+ylab("Log10-transformed total phosphorus (ug/L)")+
  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),
             aes(x=TPscore,y=TP, color="red3"),alpha = 0.7)+
  # geom_smooth(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),method="glm",
  #             se=F,aes(x=TPscore,y=TP),size=0.1,linetype="dashed")+
  # geom_smooth(data=TPmodel.initial, se=F,aes(x=x,y=predicted),color="red3")+
labs(title="")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
corrected<- plot(TPmodel.corrected,line_size = 1,limit_range = F, color="green3")+ 
  xlab("Corrected index score")+ylab('')+
  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP")]),
             aes(x=TPscore.D.5,y=TP, color="green3"),alpha = 0.7)+
  xlim(0,10)+labs(title="")+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
initial2 <- plot(TPmodel.initial,line_size = 1,limit_range = F)+ 
  xlab("Initial index score")+ylab("Log10-transformed total phosphorus (ug/L)")+
  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP","region","season")]),
             aes(x=TPscore,y=TP,group=region,color=season),alpha = 0.7)+
  scale_color_manual(values = season_colors)+xlim(0,10)+labs(title="",color="Season")+scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')+theme_bw()
corrected2<- plot(TPmodel.corrected,line_size = 1,limit_range = F)+ 
  xlab("Corrected index score")+ylab("")+
  geom_point(data=na.omit(test[,c("TPscore","TPscore.D.5","TP","region","season")]),
             aes(x=TPscore.D.5,y=TP,group=region,color=season),alpha = 0.7)+
  scale_color_manual(values = season_colors)+xlim(0,10)+labs(title="",color="Season")+scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))+guides(color='none')
grid.arrange(initial,corrected,ncol=2)
grid.arrange(initial2,corrected2,ncol=2)

dataTN <- na.omit(test[,c('TNscore.DNI','Year','TN','TNscore','TNscore.D','TNscore.I',"TNscore.DN.5")])
dataTP <- na.omit(test[,c('TPscore.DNI','Year','TP','TPscore.D.5','TPscore.I',"TPscore")])

initialse <- NULL
correctse <- NULL
universalse <- NULL
initialr <- NULL
correctr <- NULL
universalr <- NULL
for (i in 1:1000) {
  #Creating a resampled dataset from the sample data
  sample_d = test[sample(1:nrow(test), size = 50, replace = TRUE), ]
  
  #Running the regression on these data
  initialTN <- summary(lm(TN~1+TNscore ,data=sample_d))
  universalTN <- summary(lm(TN~1+TNscore.D.5N.5I.5 ,data=sample_d))
  correctTN <-summary(lm(TN~1+TNscore.D ,data=sample_d))
  
  
  universalse <- c(universalse, universalTN$coefficients[2,2])
  initialse <- c(initialse, initialTN$coefficients[2,2])
  correctse <- c(correctse,correctTN$coefficients[2,2])
  
  initialr <- c(initialr,initialTN$r.squared)
  universalr <- c(universalr,universalTN$r.squared)
  correctr <- c(correctr, correctTN$r.squared)

}



tncoefse <- cbind(initialse,correctse,universalse)
tncoefse2 <- pivot_longer(data.frame(tncoefse),cols=1:3, names_to = "Index",values_to = "Coefficent SE")

tnr <- cbind("Initial index"=initialr, "Top index"=correctr, "Suggested index" = universalr)
tnr2 <- pivot_longer(data.frame(tnr),cols=1:3, names_to = "Index",
                     values_to = "Rsq")


tnr2 <- tnr2 %>% mutate(Index = case_when(Index == "Initial.index" ~ "Initial index",
                                          Index == "Top.index" ~ "Top index",
                                          Index == "Suggested.index" ~ "Suggested index"))


tnr2$Index <- factor(tnr2$Index, levels = c("Initial index","Top index", "Suggested index"))




ggplot(data=tnr2)+
  geom_jitter(aes(x=Index,y=Rsq,color=Index))+
  geom_boxplot(aes(x=Index,y=Rsq,color=Index),color="black", fill=NA, size=1)+
  theme_classic()




ggstatsplot::ggbetweenstats(data=tnr2 ,x=Index,y=Rsq,plot.type = "box", ylab = "Adjusted r-squared", 
                            results.subtitle = F, bf.message = F, 
                            centrality.plotting = F,xlab=c(" "),ggsignif.args = list(textsize =3, tip_length=0.001),
                            point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), 
                                              alpha = 0.4, size = 2, stroke = 0),
                            violin.args = list(fill = NA, color = "transparent"), boxplot.args= list(alpha=0),
                            ggplot.component = list( scale_color_manual(values=c("red3","green3","black"))))
t.test(Rsq~Index,data=tnr2)
ggstatsplot::ggbetweenstats(data=tncoefse2 ,x=Index,y=`Coefficent SE`,plot.type = "box", ylab = "SE of coefficient", 
                            results.subtitle = F, bf.message = F, 
                            centrality.plotting = F,xlab=c(" "),ggsignif.args = list(textsize =3, tip_length=0.001),
                            point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), 
                                              alpha = 0.4, size = 2, stroke = 0),
                            violin.args = list(fill = NA, color = "transparent"), boxplot.args= list(alpha=0),
)
t.test(`Coefficent SE`~Index,data=tncoefse2)
