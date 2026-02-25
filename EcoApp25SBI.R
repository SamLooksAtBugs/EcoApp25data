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
library(missMDA)

sample_file <- "frag_comm_enviro_trait.xlsx"
traits <- read_excel(sample_file, sheet = "Traits")

traits2 <- traits %>% mutate(Dispersal.weight = 
                               case_when(TraitGroup == "Long-lived & passive" ~ 0.5,
                                         TraitGroup == "Big & winged"   ~ 0.5,
                                         TRUE  ~ 1))

setwd("C:/Users/Sam/OneDrive/Documents") #change to whatever directory you need

data.occ <- read_excel("EcoApp25Data.xlsx", sheet = "OCC")
data.epa <- read_excel("EcoApp25Data.xlsx", sheet = "EPA", guess_max = 10000) #recommend keeping the guess max
sbi <- read_excel("EcoApp25Data.xlsx", sheet = "SBI")
taxa.info <- read_excel("EcoApp25Data.xlsx", sheet = "Groups_Mechanisms")
variables <- read_excel("EcoApp25Data.xlsx", sheet = "Variables")
sbitest <- read_excel("EcoApp25Data.xlsx", sheet="SBItest")
#rs.test <- read_excel("roadsaltop.xlsx", sheet= "Test")

#bring in the new associations
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




columns_to_include <- variables$Columns
bugs <- na.omit(variables$Bugs)



#need to pivot occ to long format
occ.long <- pivot_longer(data.occ, cols = intersect(colnames(data.occ),bugs),
                         names_to = "Taxa", values_to = "Count")


occ <- occ.long[, intersect(colnames(occ.long),columns_to_include)]
epa <- data.epa[, intersect(colnames(data.epa),columns_to_include)]

occ$Dataset <- "OCC"

epa$Dataset <- "EPA"

all <- rbind(occ,epa)





all <- sbitest
#add mechanisms and orders
all <- left_join(all, taxa.info, by="Taxa")

#all <- all %>% filter(Family!="Chironomidae"&Region != "CPL")

all <- all %>% filter(Region != "CPL")

all$Dispersal <- all$SFP #*ifelse(all$DispersalType=="Aquatic",0.5,1)

all$Dispersal <- ifelse(all$Dispersal == 0, median(all$Dispersal),all$Dispersal)

all$Dispersal <- all$Dispersal #*ifelse(all$DispersalType=="Aquatic",0.5,1)

all <- left_join(all, associates, by=join_by(Taxa == taxon))

all$n_neg <- ifelse(is.na(all$n_neg),0,all$n_neg)
all$n_total <- ifelse(is.na(all$n_total),0,all$n_total)

## --- weighting helpers (0.1 to 1) ------------------------------------------
rev_minmax <- function(x, trim = 0.0) {
  lo <- as.numeric(quantile(x, trim, na.rm = TRUE))
  hi <- as.numeric(quantile(x, 1 - trim, na.rm = TRUE))
  xw <- pmin(pmax(x, lo), hi)
  if (hi - lo == 0) return(rep(1, length(x)))
  (hi - xw) / (hi - lo)        # high -> 0, low -> 1
}

rev_minmax_half <- function(x, trim = 0.05) {
  w <- rev_minmax(x, trim)
  0.5 + 0.5 * w
}

minmax_half <- function(x, trim = 0.0) {
  w <- 1- rev_minmax(x, trim)
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

mid_minmax <- function(x, trim = 0.0, center = 0.5) {
  lo <- as.numeric(quantile(x, trim, na.rm = TRUE))
  hi <- as.numeric(quantile(x, 1 - trim, na.rm = TRUE))
  xw <- pmin(pmax(x, lo), hi)  # winsorize
  
  if (hi - lo == 0) return(rep(1, length(x)))
  
  # normalize to [0,1]
  z <- (xw - lo) / (hi - lo)
  
  # distance from center (default 0.5 = midpoint of range)
  d <- abs(z - center) / max(center, 1 - center)  # scaled so max distance = 1
  
  # convert distance to weight: center -> 1, extremes -> 0.5
  1 -  d
}




rev_minmax_01 <- function(x, min_w = 0.1, max_w = 1, trim = 0.01) {
  lo <- as.numeric(quantile(x, trim, na.rm = TRUE))
  hi <- as.numeric(quantile(x, 1 - trim, na.rm = TRUE))
  
  xw <- pmin(pmax(x, lo), hi)  # winsorize
  
  if (hi - lo == 0) return(rep(max_w, length(x)))
  
  w <- (hi - xw) / (hi - lo)  # reversed
  min_w + w * (max_w - min_w)
}


mul_weights <- function(df, vars) {
  if (length(vars) == 0) return(rep(1, nrow(df)))
  Reduce(`*`, df[vars])
}

combos <- list(
  base = character(0),
  D    = c("w_D"),
  I    = c("w_I"),
  N    = c("w_N"),
  DI   = c("w_D","w_I"),
  DN   = c("w_D","w_N"),
  IN   = c("w_I","w_N"),
  DNI  = c("w_D","w_I","w_N")
)

## --- 1) assign Score by region (you already did this) -----------------------
## all currently has correct region-specific Score, plus Dispersal, n_neg, etc.

## --- 2) compute SBI per UID within-region (just to build the score), then bind

all$Chloride <- all$Chloride...3

#all <- left_join(all, traits2[c("Taxa","Dispersal.weight")])

taxon_traits <- all %>%
  group_by(Taxa) %>%
  summarise(
    NicheBreadth_tax = median(NicheBreadth, na.rm = TRUE),
    n_total_tax      = median(n_total, na.rm = TRUE),
    Dispersal_tax    = median(Dispersal, na.rm = TRUE),
   # Dispersal.weight    = median(Dispersal.weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    w_N = rev_minmax_half(NicheBreadth_tax),
    w_I = rev_minmax_half(n_total_tax),
    w_D = rev_minmax_half(Dispersal_tax)
    #w_D = Dispersal.weight
  ) %>%
  dplyr::select(Taxa, w_N, w_I, w_D)

all <- all %>%
 # select(-w_N, -w_I, -w_D) %>%   # drop any old ones
  left_join(taxon_traits, by = "Taxa") %>%
  mutate(
    logcount = log1p(Count),
    logcl    = log10(Chloride)
  )

all.split <- split(all, all$Region)

sbi_long <- sbi %>%
  pivot_longer(
    cols = -c(Phylum, Class, Order, Family, Taxa.original, Taxa),
    names_to = "Region",
    values_to = "Score"
  ) %>%
  filter(!is.na(Score))



scores_by_uid <- imap_dfr(all.split, function(sub, reg) {
  
  # region-specific score table
  sreg_taxa <- sbi_long %>%
    filter(Region == reg) %>%
    dplyr::select(Taxa, Score)
  
  sreg_fam <- sbi_long %>%
    filter(Region == reg) %>%
    dplyr::select(Family, Score)
  
  # join score by Taxa, then fill missing by Family
  sub <- sub %>%
    left_join(sreg_taxa, by = "Taxa") %>%
    left_join(sreg_fam, by = "Family", suffix = c(".taxa", ".fam")) %>%
    mutate(
      Score = coalesce(Score.taxa, Score.fam)
    #  Score =Score.taxa
    ) %>%
    dplyr::select(-Score.taxa, -Score.fam)
  
  # NOTE: weights can be region-specific (as written) OR global.
  # # This version makes weights region-specific because it’s computed inside this split.
  # sub <- sub %>%
  #   mutate(
  #     w_I = rev_minmax_half(n_total),
  #     w_N = rev_minmax_half(NicheBreadth),
  #     w_D = mid_minmax_half(Dispersal),
  #     logcount   = log1p(Count),
  #     Score_base = Score * logcount,
  #     logcl      = log10(Chloride)
  #   )

  sub <- sub %>%
    mutate(
      Score_base = Score * logcount
    )
  
  map_dfr(names(combos), function(nm) {
    w <- mul_weights(sub, combos[[nm]])
    
    sub %>%
      mutate(
        LC_w = logcount   * w,
        TN_w = Score_base * w
      ) %>%
      group_by(UID) %>%
      summarise(
        Region   = first(Region),
       # Dataset  = first(Dataset),
        Chloride = first(Chloride),
        logcl    = first(logcl),
        SBI      = sum(TN_w, na.rm = TRUE) / sum(LC_w, na.rm = TRUE),
        combo    = nm,
        .groups  = "drop"
      )
  })
})

## --- 3) (optional but recommended) enforce SAME UIDs across combos ----------
common_uids <- scores_by_uid %>%
  filter(!is.na(SBI), !is.na(logcl)) %>%
  group_by(UID) %>%
  summarise(n_combo = n_distinct(combo), .groups = "drop") %>%
  filter(n_combo == length(combos)) %>%
  pull(UID)

scores_pooled <- scores_by_uid %>%
  filter(UID %in% common_uids)

## --- 4) run pooled models (NOT by region) -----------------------------------
mae  <- function(actual, pred) mean(abs(actual - pred), na.rm = TRUE)
rmse <- function(actual, pred) sqrt(mean((actual - pred)^2, na.rm = TRUE))

fit_table <- scores_pooled %>%
  group_by(combo) %>%
  group_modify(~{
    df <- .x
    mod <- lm(logcl ~ SBI, data = df)
    s <- summary(mod)
    coefs <- s$coefficients
    
    y <- df$logcl
    yhat <- fitted(mod)
  #  beta_std <- coefs["SBI","Estimate"] * sd(x, na.rm=TRUE) / sd(y, na.rm=TRUE)
    
    tibble(
      n        = nrow(df),
      estimate = coefs["SBI","Estimate"],
      std_error= coefs["SBI","Std. Error"],
      t_value  = coefs["SBI","t value"],
      p_value  = coefs["SBI","Pr(>|t|)"],
     # beta_std = beta_std,
      R2       = unname(s$r.squared),
      R2adj    = unname(s$adj.r.squared),
      RMSE     = rmse(y, yhat),
      MAE      = mae(y, yhat),
      AIC      = AIC(mod),
      BIC      = BIC(mod)
    )
  }) %>%
  ungroup() %>%
  arrange(AIC)



# Split data by combo
scores_split <- split(scores_pooled, scores_pooled$combo)

# Fit models and store in named list
sbi_models <- map(scores_split, ~ lm(logcl ~ SBI, data = .x))


# View(scores_pooled)
# View(fit_table)
fit_table


#boxplots new ####

filtered_datasets <- scores_by_uid %>%
  split(.$combo)
# ---- thresholds (same as you have) ----
thresholds <- tribble(
  ~Region, ~low_cut, ~high_cut,
  "SAP",   5,   30,
  "NAP",   4,   75,
  "UMW",   2,   40,
  "TPL",   40,  120,
  "SPL",   50,  230,
  "XER",   2,   50,
  "WMT",   1,   10,
  "NPL",   3,   50
)

# ---- 1) classify chloride into Low/Moderate/High ----
add_value_class <- function(df, thresholds_df = thresholds) {
  df %>%
    left_join(thresholds_df, by = "Region") %>%
    mutate(
      Value_class = case_when(
        Chloride < low_cut ~ "Low",
        Chloride >= low_cut & Chloride <= high_cut ~ "Moderate",
        Chloride > high_cut ~ "High",
        TRUE ~ NA_character_
      ),
      Value_class = factor(Value_class,
                           levels = c("Low", "Moderate", "High"),
                           ordered = TRUE
      )
    )
}

u.need.classes <- add_value_class(scores_pooled)

library(dplyr)
library(effsize)

library(dplyr)

#spearman rank ordinal
spearman_by_region_combo <- u.need.classes %>%
  group_by(combo, Region) %>%
  group_modify(~{
    df <- .x %>%
      filter(!is.na(SBI), !is.na(Value_class))
    
    # need at least 3 observations and >1 class level to compute correlation
    if (nrow(df) < 3 || n_distinct(df$Value_class) < 2) {
      return(tibble(n = nrow(df), rho = NA_real_, p_value = NA_real_))
    }
    
    # Value_class is ordered factor (Low < Moderate < High)
    vt <- as.numeric(df$Value_class)
    
    ct <- suppressWarnings(cor.test(df$SBI, vt, method = "spearman", exact = FALSE))
    
    tibble(
      n       = nrow(df),
      rho     = unname(ct$estimate),
      p_value = ct$p.value
    )
  }) %>%
  ungroup() %>%
  arrange(Region, combo)

spearman_by_region_combo


spearman_by_combo <- u.need.classes %>%
  group_by(combo) %>%
  group_modify(~{
    df <- .x %>%
      filter(!is.na(SBI), !is.na(Value_class))
    
    # need at least 3 observations and >1 class level to compute correlation
    if (nrow(df) < 3 || n_distinct(df$Value_class) < 2) {
      return(tibble(n = nrow(df), rho = NA_real_, p_value = NA_real_))
    }
    
    # Value_class is ordered factor (Low < Moderate < High)
    vt <- as.numeric(df$Value_class)
    
    ct <- suppressWarnings(cor.test(df$SBI, vt, method = "spearman", exact = FALSE))
    
    tibble(
      n       = nrow(df),
      rho     = unname(ct$estimate),
      p_value = ct$p.value
    )
  }) %>%
  ungroup() %>%
  arrange(combo)

spearman_by_combo


library(dplyr)
library(pROC)
library(tibble)

# helper: compute AUC for a given pair within a data frame
auc_pair <- function(df, g1, g2, min_n = 5) {
  d <- df %>%
    filter(Value_class %in% c(g1, g2)) %>%
    mutate(Value_class = factor(Value_class, levels = c(g1, g2))) %>%
    filter(!is.na(SBI), !is.na(Value_class))
  
  n1 <- sum(d$Value_class == g1, na.rm = TRUE)
  n2 <- sum(d$Value_class == g2, na.rm = TRUE)
  
  if ((n1 < 2) || (n2 < 2) || (n1 + n2) < min_n) {
    return(tibble(
      comparison = paste0(g1, " vs ", g2),
      n = n1 + n2,
      n1 = n1, n2 = n2,
      AUC = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_
    ))
  }
  
  roc_obj <- roc(
    response  = d$Value_class,
    predictor = d$SBI,
    levels    = c(g1, g2),
    direction = "<",   # higher SBI -> g2 (the "positive" class)
    quiet     = TRUE
  )
  ci_vals <- ci.auc(roc_obj)
  
  tibble(
    comparison = paste0(g1, " vs ", g2),
    n = n1 + n2,
    n1 = n1, n2 = n2,
    AUC = as.numeric(auc(roc_obj)),
    CI_low = ci_vals[1],
    CI_high = ci_vals[3]
  )
}

pairs <- list(
  c("Low","Moderate"),
  c("Low","High"),
  c("Moderate","High")
)

auc_pairwise_by_region_combo <- u.need.classes %>%
  group_by(combo, Region) %>%
  group_modify(~{
    df <- .x %>% filter(!is.na(SBI), !is.na(Value_class))
    bind_rows(lapply(pairs, \(p) auc_pair(df, p[1], p[2])))
  }) %>%
  ungroup() %>%
  arrange(Region, combo, comparison)

auc_pairwise_by_region_combo


auc_pairwise_by_combo <- u.need.classes %>%
  group_by(combo) %>%
  group_modify(~{
    df <- .x %>% filter(!is.na(SBI), !is.na(Value_class))
    bind_rows(lapply(pairs, \(p) auc_pair(df, p[1], p[2])))
  }) %>%
  ungroup() %>%
  arrange( combo, comparison)

auc_pairwise_by_combo


#all together now, final boxplots####
SBIinitial_df   <- add_value_class(filtered_datasets$base) %>% mutate(IndexType = "Initial")
SBItop_df     <- add_value_class(filtered_datasets$IN)   %>% mutate(IndexType = "Corrected")
SBIuniversal_df <- add_value_class(filtered_datasets$DNI)  %>% mutate(IndexType = "Universal")

sbi_all <- dplyr::bind_rows(SBIinitial_df, SBItop_df, SBIuniversal_df)

library(dplyr)
library(rstatix)

gh_all <- sbi_all %>%
  filter(!is.na(SBI), !is.na(Value_class)) %>%
  group_by(Region, IndexType) %>%
  games_howell_test(SBI ~ Value_class) %>%
  ungroup()

gh_all <- sbi_all %>%
  filter(!is.na(SBI), !is.na(Value_class)) %>%
  group_by(Region, IndexType) %>%
  games_howell_test(SBI ~ Value_class) %>%
  ungroup()

gh_sig <- gh_all %>%
  filter(p.adj < 0.05)

gh_sig <- gh_sig %>%
  add_xy_position(
    x = "IndexType",
    group = "Value_class",
    dodge = 0.75,
    data = sbi_all
  )%>%
  mutate(y.position = y.position - 1.5)

gh_sig$p.adj.signif <- ""

pd <- position_dodge2(width = 0.75, preserve = "single")

sbi_all <- sbi_all %>% 
  mutate(
    IndexType = case_when(
      IndexType == "Universal" ~ "Suggested", 
      IndexType == "Corrected" ~ "Top",
      IndexType == "Initial"   ~ "Initial"
    ),
    IndexType = factor(IndexType, levels = c("Initial", "Top", "Suggested"))
  )

p <- ggplot(sbi_all, aes(x = IndexType, y = SBI)) +
  geom_point(aes( color = Value_class),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
             alpha = 0.8
  ) +
  geom_boxplot(
    aes(group = interaction(IndexType, Value_class)),color="black",alpha=0,
    position = pd,
    outlier.shape = NA
  ) +

  facet_wrap(~Region, nrow=2) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0,10,2)) +
  scale_color_manual(values = c("lightblue", "skyblue2", "darkblue")) +
  coord_cartesian(ylim = c(0,12)) +
  labs(x="Index",y="SBI score", color="Chloride levels")+
  theme(axis.text.x = element_text(hjust=1,angle=45),text = element_text(size = 16)) +
  stat_pvalue_manual(
    gh_sig,
    label = "p.adj.signif",
    hide.ns = TRUE,      tip.length = 0.01,
    size = 0.6,
    bracket.size = 0.8
  )

p




initial.predict <- ggeffects::ggpredict(sbi_models$base,terms = "SBI [0:10]")

initial <- ggplot()+
  geom_point(data=SBIinitial_df, 
             aes(x=SBI,y=logcl),color="red3",alpha=0.6)+
  geom_ribbon(data=initial.predict, 
              aes(x=x,ymin = conf.low,ymax = conf.high),
              alpha=0.3,fill="red3")+
  geom_line(data=initial.predict,
            aes(x=x,y=predicted),color="red3")+
  theme_classic()+
  labs(title = "(a)", x="Initial SBI", y=expression("Log"["10"]*" transformed chloride (mg/L)"))+
  scale_x_continuous(breaks =c(0,2,4,6,8,10),limits = c(0,10))


correct.predict <- ggeffects::ggpredict(sbi_models$IN,terms = "SBI [0:10]")


final <- ggplot()+
  geom_point(data=SBItop_df , 
             aes(x=SBI,y=logcl),color="green3", alpha =0.6)+  
  geom_ribbon(data=correct.predict, 
              aes(x=x,ymin = conf.low,ymax = conf.high),
              alpha=0.3,fill="green3")+
  geom_line(data=correct.predict,
            aes(x=x,y=predicted),color="green3")+
  theme_classic()+
  labs(title = "(b)", x="Top modified SBI", y="")+
  scale_x_continuous(breaks =c(0,2,4,6,8,10),limits = c(0,10))


universal.predict <- ggeffects::ggpredict(sbi_models$DNI,terms = "SBI [0:10]")


universal <- ggplot()+
  geom_point(data=SBIuniversal_df, 
             aes(x=SBI,y=logcl),color="black", alpha=0.6)+
  geom_ribbon(data=universal.predict, 
              aes(x=x,ymin = conf.low,ymax = conf.high),
              alpha=0.3,fill="black")+
  geom_line(data=universal.predict,
            aes(x=x,y=predicted),color="black")+
  theme_classic()+
  labs(title = "(c)", x="Suggested modified SBI", y="")+
  scale_x_continuous(breaks =c(0,2,4,6,8,10),limits = c(0,10))



gridExtra::grid.arrange(initial,final, universal, ncol=3)


# ---- 2) Games-Howell results (brackets only, no labels) ----
get_gh_brackets_by_region <- function(df, yvar = SBI) {
  yname <- rlang::as_name(rlang::ensym(yvar))
  
  df %>%
    dplyr::group_by(Region) %>%
    dplyr::group_modify(~{
      dsub <- .x %>%
        dplyr::filter(!is.na(.data[[yname]]), !is.na(Value_class)) %>%
        dplyr::mutate(Value_class = droplevels(Value_class))  # ✅ key for SPL
      
      # Need at least 2 classes present (SPL is fine with 2)
      if (dplyr::n_distinct(dsub$Value_class) < 2) {
        return(tibble::tibble())
      }
      
      gh <- rstatix::games_howell_test(dsub, reformulate("Value_class", yname)) %>%
        dplyr::filter(p.adj < 0.05) %>%          # keep sig only
        dplyr::mutate(p.adj.signif = "")         # brackets only
      
      if (nrow(gh) == 0) return(tibble::tibble())
      
      gh <- gh %>%
        rstatix::add_xy_position(x = "Value_class", step.increase = 0.5) %>%
        dplyr::mutate(
          y.position = max(dsub[[yname]], na.rm = TRUE) + 0.5 + dplyr::row_number() * 0.25
        )
      
      gh
    }) %>%
    dplyr::ungroup()
}


# ---- 3) boxplot builder ----
make_sbi_boxplot_region <- function(df, title_tag = "(a)", ylab = "Index",
                                    ybreaks = c(0,2,4,6,8,10), ylims = c(0,12)) {
  
  gh <- get_gh_brackets_by_region(df, SBI)
  
  p <- ggplot(df, aes(x = Value_class, y = SBI, color = Value_class)) +
    geom_jitter(width = 0.15, alpha = 0.8) +
    geom_boxplot(color = "black", fill = "white", alpha = 0.2,
                 outlier.shape = NA, size = 1) +
    theme_classic() +
    facet_wrap(~Region) +
    labs(x = "", y = ylab, title = title_tag) +
    scale_color_manual(values = c("lightblue", "skyblue2", "darkblue")) +
    theme(legend.position = "none", text = element_text(size = 12)) +
    scale_y_continuous(breaks = ybreaks, limits = ylims) +
    ggpubr::stat_pvalue_manual(
      gh,
      label = "p.adj.signif",
      hide.ns = TRUE,
      tip.length = 0.01,
      size = 0.6,
      bracket.size = 0.8
    ) +
    expand_limits(y = ylims[2])
  
  list(plot = p, gh = gh)
}
# ---- 4) apply to each dataset (no repeating) ----
SBIinitial_df   <- add_value_class(filtered_datasets$base)
SBItop_df       <- add_value_class(filtered_datasets$IN)
SBIuniversal_df <- add_value_class(filtered_datasets$DNI)

# optional: quick oneway tests (same as before)
oneway.test(SBI ~ Value_class, data = SBIinitial_df)
oneway.test(SBI ~ Value_class, data = SBItop_df)
oneway.test(SBI ~ Value_class, data = SBIuniversal_df)

# build plots + GH tables
initial_out   <- make_sbi_boxplot_region(SBIinitial_df,   title_tag="(a)", ylab="Initial index")
top_out       <- make_sbi_boxplot_region(SBItop_df,       title_tag="(b)", ylab="Top modified index")
universal_out <- make_sbi_boxplot_region(SBIuniversal_df, title_tag="(c)", ylab="Suggested modified index")

SBIinitialbox   <- initial_out$plot
SBItopbox       <- top_out$plot
SBIuniversalbox <- universal_out$plot

ggpubr::ggarrange(SBIinitialbox, SBItopbox, SBIuniversalbox, ncol = 3)

# combine GH results (now guaranteed consistent)
allSBIgames <- dplyr::bind_rows(
  initial_out$gh  %>% dplyr::mutate(which = "initial"),
  top_out$gh      %>% dplyr::mutate(which = "top"),
  universal_out$gh%>% dplyr::mutate(which = "universal")
)


eta_by_combo <- u.need.classes %>%
  group_by(combo) %>%
  group_modify(~{
    m <- aov(SBI ~ Value_class, data = .x)
    es <- effectsize::cohens_f(m)
    tibble(
      n = nrow(.x),
     epsilon = es
    )
  }) %>%
  ungroup() %>%
  arrange(desc(epsilon))

eta_by_combo



#whatever##### 
# 
# # Apply each correction combination
# datasets <- pmap(
#   correction_grid,
#   function(interactions_corr, niche_corr, dispersal_corr) {
#     all %>%
#       mutate(
#         correction_code = paste0("I", interactions_corr, "_N", niche_corr, "_D", dispersal_corr),
#         Corrected.count = log(1+(Count *
#                                    ifelse(Interactions >= CooccurCut, interactions_corr, 1) *
#                                    ifelse(NicheBreadth >= NicheCut, niche_corr, 1) *
#                                    ifelse(Dispersal >= DispersalCut, dispersal_corr, 1)
#         )))
#   }
# )
# 
# # Name each dataset by its correction code
# names(datasets) <- apply(correction_grid, 1, function(x) {
#   paste0("I", x[1], "_N", x[2], "_D", x[3])
# })


#metric calculation

finaldatas <- list()

for (dataset in names(datasets)){
  
  df <- datasets[[dataset]]
  
  ept_orders <- c("Ephemeroptera", "Plecoptera", "Trichoptera")  # define if not already
  
  scores <- df %>%
    group_by(UID) %>%
    summarise(
      SBI = sum(Corrected.count * Score, na.rm = TRUE) / 
        sum(Corrected.count, na.rm = TRUE),
      
      Region = first(Region),
      Chloride = first(Chloride),
      
      .groups = "drop"
    )
  

  #PCA
  #environmental variables
  enviro <- c("DO", "Temp", "Turb", "Cond", "pH", "Habitat", "TP", "TN", "NH4", "NO3")
  enviro2 <- c("DO", "Temp", "Turb", "Cond", "Habitat")
  #scale
  # score.split <- split(scores, scores$Dataset)
  # 
  # #occ
  # score.split$OCC[enviro] <- scale(score.split$OCC[enviro])
  # 
  # #epa
  # score.split$EPA[enviro] <- scale(score.split$EPA[enviro])
  # 
  # #rejoin
  # newscores <- rbind(score.split$OCC,score.split$EPA)
  # 
  # newscores <- newscores %>% filter(Dataset == "OCC")
  
  newscores <- scores
  
  newscores$logcl <- log10(newscores$Chloride)
  

  
  finaldatas[[dataset]] <- newscores
  
  print(dataset)
  print(summary(newscores[c("SBI")]))
  
}




# Step 1: Get valid SAMPLEs from each dataset (non-NA for both OKSum and MOSum)
valid_samples_list <- lapply(finaldatas, function(df) {
  df %>%
    filter(!is.na(SBI)) %>%
    pull(UID) %>%
    unique()
})

# Step 2: Find the intersection across all datasets
common_samples <- Reduce(intersect, valid_samples_list)

# Step 3: Filter all datasets to include only those common SAMPLE values
filtered_datasets <- lapply(finaldatas, function(df) {
  df %>% filter(UID %in% common_samples)
})



# check to ensure everything looks right
for (dataset in names(datasets)){
  print(summary(filtered_datasets[[dataset]][c("SBI")]))
}


fits <- list()
sbi.models <- list()
corrs <- list()



for (data in names(filtered_datasets)) {
  
  df <- filtered_datasets[[data]]
  
  #df <- df %>% filter(SBI > 0)
  
  #df <- df2 %>% filter(Cond < 2.5)
  
  model.sbi <- lm(logcl~SBI,data=df)
  sbi.sum <- summary(model.sbi)
  coefs <- coefficients(sbi.sum)
  pred.sbi <- predict(model.sbi, type = "response")
  aic.sbi <- AIC(model.sbi)
  r2.sbi <- sbi.sum$r.squared
  r2adj.sbi <- sbi.sum$adj.r.squared
  
  
  obs <- na.omit(df$logcl)
  
  rmse.sbi <- sqrt(mean((obs - pred.sbi)^2))
  mae.sbi <- mean(abs(obs-pred.sbi))
  
  
  fits[[data]] <- data.frame(Correction = data,
                             estimate = coefs[2, 1],
                             std_error = coefs[2, 2],
                             t_value = coefs[2, 3],
                             p_value = coefs[2, 4],
                             R2 = r2.sbi,
                             R2adj =  r2adj.sbi,
                             RMSE = rmse.sbi,
                             MAE = mae.sbi,
                             AIC = aic.sbi
                             )
  
  sbi.models[[data]] <- model.sbi
  #mo.models[[data]] <- model.mo
  
  corrs[[data]] <- cor(df[,c("SBI","logcl","Chloride")],use="complete") 
}

all.fits <- do.call(rbind, fits)
View(all.fits)


#plots ####
initial.predict <- ggeffects::ggpredict(sbi.models$I1_N1_D1,terms = "SBI [0:10]")

initial <- ggplot()+
  geom_point(data=filtered_datasets$I1_N1_D1, 
             aes(x=SBI,y=logcl),color="red3",alpha=0.6)+
  geom_ribbon(data=initial.predict, 
              aes(x=x,ymin = conf.low,ymax = conf.high),
              alpha=0.3,fill="red3")+
  geom_line(data=initial.predict,
            aes(x=x,y=predicted),color="red3")+
  theme_classic()+
  labs(title = "(a)", x="Initial SBI", y=expression("Log"["10"]*" transformed chloride (mg/L)"))+
  scale_x_continuous(breaks =c(0,2,4,6,8,10),limits = c(0,10))


correct.predict <- ggeffects::ggpredict(sbi.models$I1_N1_D1,terms = "SBI [0:10]")


final <- ggplot()+
  geom_point(data=filtered_datasets$I0_N0_D1, 
             aes(x=SBI,y=logcl),color="green3", alpha =0.6)+  
  geom_ribbon(data=initial.predict, 
            aes(x=x,ymin = conf.low,ymax = conf.high),
                alpha=0.3,fill="green3")+
  geom_line(data=correct.predict,
            aes(x=x,y=predicted),color="green3")+
  theme_classic()+
  labs(title = "(b)", x="Top modified SBI", y="")+
  scale_x_continuous(breaks =c(0,2,4,6,8,10),limits = c(0,10))


universal.predict <- ggeffects::ggpredict(sbi.models$I0.5_N0.5_D0.5,terms = "SBI [0:10]")


universal <- ggplot()+
  geom_point(data=filtered_datasets$I0.5_N0.5_D0.5, 
             aes(x=SBI,y=logcl),color="black", alpha=0.6)+
  geom_ribbon(data=universal.predict, 
              aes(x=x,ymin = conf.low,ymax = conf.high),
              alpha=0.3,fill="black")+
  geom_line(data=initial.predict,
            aes(x=x,y=predicted),color="black")+
  theme_classic()+
  labs(title = "(c)", x="Suggested modified SBI", y="")+
  scale_x_continuous(breaks =c(0,2,4,6,8,10),limits = c(0,10))
  


gridExtra::grid.arrange(initial,final, universal, ncol=3)


#boxplots #####




thresholds <- tribble(
  ~Region, ~low_cut, ~high_cut,
  "SAP",   5,   30,
  "NAP",   4,   75,
  "UMW",   2,   40,
  "TPL",   40,  120,
  "SPL",   50,  230,
  "XER",   2,   50,
  "WMT",   1,   10,
  "NPL",   3,   50
)

# Join cutoffs and classify
SBIinitial <- filtered_datasets$I1_N1_D1 %>%
  left_join(thresholds, by = "Region") %>%
  mutate(
    Value_class = case_when(
      Chloride < low_cut              ~ "Low",
      Chloride >= low_cut & Chloride <= high_cut ~ "Moderate",
      Chloride > high_cut             ~ "High",
      TRUE                         ~ NA_character_
    ),
    Value_class = factor(Value_class, 
                         levels = c("Low", "Moderate", "High"), 
                         ordered = TRUE)
  )

SBItop <- filtered_datasets$I0_N0_D1 %>%
  left_join(thresholds, by = "Region") %>%
  mutate(
    Value_class = case_when(
      Chloride < low_cut              ~ "Low",
      Chloride >= low_cut & Chloride <= high_cut ~ "Moderate",
      Chloride > high_cut             ~ "High",
      TRUE                         ~ NA_character_
    ),
    Value_class = factor(Value_class, 
                         levels = c("Low", "Moderate", "High"), 
                         ordered = TRUE)
  )


SBIuniversal <- filtered_datasets$I0.5_N0.5_D0.5 %>%
  left_join(thresholds, by = "Region") %>%
  mutate(
    Value_class = case_when(
      Chloride < low_cut              ~ "Low",
      Chloride >= low_cut & Chloride <= high_cut ~ "Moderate",
      Chloride > high_cut             ~ "High",
      TRUE                         ~ NA_character_
    ),
    Value_class = factor(Value_class, 
                         levels = c("Low", "Moderate", "High"), 
                         ordered = TRUE)
  )


SBIinitial_gh_results <- SBIinitial %>%
  games_howell_test(SBI ~ Value_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%   # make empty labels
  add_xy_position(x = "Value_class", step.increase = 0.5) %>%
  # put the FIRST bracket at (max data + 0.5), then stagger by 0.25
  mutate(y.position = max(SBIinitial$SBI, na.rm = TRUE) + 0.5 + row_number() * 0.25)

oneway.test(SBI ~ Value_class, data=SBIinitial)

SBIinitialbox <- ggplot(data=SBIinitial, aes(x=Value_class,y=SBI, color=Value_class))+
  geom_jitter()+
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0,12))+
  
  labs(x = "", y = "Initial index", title = "(a)") +
  scale_color_manual(values = c("lightblue", "skyblue2", "darkblue")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
stat_pvalue_manual(
  SBIinitial_gh_results,
  label = "p.adj.signif",   # "" → brackets only
  hide.ns = TRUE,
  tip.length = 0.01,
  size = 0.6,               # bracket line thickness
  bracket.size = 0.8        # widen/thicken brackets
)

SBItop_gh_results <- SBItop %>%
  games_howell_test(SBI ~ Value_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%
  add_xy_position(x = "Value_class", step.increase = 0.5) %>%
  # put the FIRST bracket at (max data + 0.5), then stagger by 0.25
  mutate(y.position = max(SBItop$SBI, na.rm = TRUE) + 0.5 + row_number() * 0.25)


oneway.test(SBI ~ Value_class, data=SBItop)


SBItopbox <- ggplot(data=SBItop, aes(x=Value_class,y=SBI, color=Value_class))+
  geom_jitter()+
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  labs(x = "", y = "Top modified index", title = "(b)") +
  scale_color_manual(values = c("lightblue", "skyblue2", "darkblue")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0,12))+
  stat_pvalue_manual(
    SBItop_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )+expand_limits(y = 12) 


SBIuniversal_gh_results <- SBIuniversal %>%
  games_howell_test(SBI ~ Value_class) %>%
  adjust_pvalue(method = "holm") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.adj.signif = "") %>%
  add_xy_position(x = "Value_class", step.increase = 0.5) %>%
  # put the FIRST bracket at (max data + 0.5), then stagger by 0.25
  mutate(y.position = max(SBIuniversal$SBI, na.rm = TRUE) + 0.5 + row_number() * 0.25)

oneway.test(SBI ~ Value_class, data=SBIuniversal)


SBIuniversalbox <- ggplot(data=SBIuniversal, aes(x=Value_class,y=SBI, color=Value_class))+
  geom_jitter()+
  geom_boxplot(color = "black", fill = "white", alpha = 0.2, outlier.shape = NA, size=1) +
  theme_classic() +
  labs(x = "", y = "Suggested modified index", title = "(c)") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0,12))+
  scale_color_manual(values = c("lightblue", "skyblue2", "darkblue")) +
  theme(legend.position = "none", text = element_text(size = 12)) +
  stat_pvalue_manual(
    SBIuniversal_gh_results,
    label = "p.adj.signif",   # "" → brackets only
    hide.ns = TRUE,
    tip.length = 0.01,
    size = 0.6,               # bracket line thickness
    bracket.size = 0.8        # widen/thicken brackets
  )+expand_limits(y = 12) 

ggarrange(SBIinitialbox,SBItopbox,SBIuniversalbox,ncol=3)

allSBIgames<- rbind(SBIuniversal_gh_results,
                   SBItop_gh_results,
                   SBIinitial_gh_results)
