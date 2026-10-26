# STAT201 project: hourly day-ahead prices in NO1, 2016 to early 2022
# Gaussian GLM (identity link) and GAM.
#
# Run from the project folder (or from R/):   Rscript R/project.R
# Input : data/raw/dataset_timeseries.csv, data/raw/nordpool_2016-01-01_2022-11-30.h5
# Output: data/processed/no1_clean.csv, figures/*.png, tables/*.tex, tables/*.csv,
#         tables/numbers.tex (all numbers used in the report, slides and talk script)

rm(list = ls())
if (!dir.exists("data/raw") && dir.exists("../data/raw")) setwd("..")
stopifnot(dir.exists("data/raw"))

suppressPackageStartupMessages({
  library(data.table)  # reading the large csv
  library(hdf5r)       # reading the raw Nord Pool file
  library(car)         # vif()
  library(lmtest)      # dwtest(), coeftest()
  library(sandwich)    # NeweyWest() HAC standard errors
  library(mgcv)        # gam()
})

for (d in c("data/processed", "figures", "tables")) dir.create(d, showWarnings = FALSE)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
num <- list()   # every number used in the text goes in here -> tables/numbers.tex
# numbers for LaTeX: "{,}" as thousands separator and \ensuremath{-} for minus,
# so the same macro works in text and in math mode
f <- function(x, d = 0) {
  x <- as.numeric(x)
  x[abs(x) < 0.5 * 10^-d] <- 0                       # no "-0.00"
  out <- formatC(abs(x), format = "f", digits = d, big.mark = "{,}")
  ifelse(x < 0, paste0("\\ensuremath{-}", out), out)
}
fp <- function(p) {             # p-values
  ifelse(p < 0.001, "\\ensuremath{<0.001}", formatC(p, format = "f", digits = 3))
}
write_tab <- function(df, file, align = NULL) {
  # booktabs tabular (no float), caption and label are set in the .tex files
  if (is.null(align)) align <- paste0("l", strrep("r", ncol(df) - 1))
  out <- c(sprintf("\\begin{tabular}{%s}", align), "\\toprule",
           paste(paste(names(df), collapse = " & "), "\\\\"), "\\midrule",
           apply(df, 1, function(r) paste(paste(r, collapse = " & "), "\\\\")),
           "\\bottomrule", "\\end{tabular}")
  writeLines(out, file.path("tables", file))
}
png_open <- function(file, w = 9, h = 5.5) {
  png(file.path("figures", file), width = w, height = h, units = "in", res = 200)
}

# ------------------------------------------------------------------
# 1. Read raw data
# ------------------------------------------------------------------
raw <- fread("data/raw/dataset_timeseries.csv",
             select = c("Date", "NO1_da_price", "NO1_cons", "NO1_buy",
                        "NO1_Fyllingsgrad_NVE"))
raw[, t_utc := as.POSIXct(Date, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")]
stopifnot(!anyNA(raw$t_utc), !anyDuplicated(raw$t_utc))
setorder(raw, t_utc)
num$nRaw <- f(nrow(raw))
num$rawStart <- format(min(raw$t_utc), "%d %B %Y %H:%M")
num$rawEnd <- format(max(raw$t_utc), "%d %B %Y %H:%M")
cat("Raw rows:", nrow(raw), "\n")

# Missing hours in the UTC time axis (UTC has no daylight saving, so every
# hour should be there exactly once)
full_axis <- seq(min(raw$t_utc), max(raw$t_utc), by = "hour")
missing_hours <- full_axis[!full_axis %in% raw$t_utc]
num$nMissingHours <- length(missing_hours)
num$missingHours <- paste(format(missing_hours, "%d %B %Y %H:%M"), collapse = " and ")
cat("Missing hours:", format(missing_hours), "\n")

# Raw Nord Pool file: tells which NO1 prices were missing in the download.
# In the csv these hours were forward filled with the previous hour's price.
h5 <- H5File$new("data/raw/nordpool_2016-01-01_2022-11-30.h5", mode = "r")
grp <- h5[["nordpool"]]
items <- grp[["block0_items"]]$read()
vals <- grp[["block0_values"]]$read()
idx <- suppressWarnings(as.numeric(grp[["axis1"]]$read()))  # ns since 1970, UTC
h5$close_all()
np <- data.table(t_utc = as.POSIXct(idx / 1e9, origin = "1970-01-01", tz = "UTC"),
                 np_price = vals[match("day_ahead_prices_NO1", items), ],
                 np_cons  = vals[match("consumption_NO1", items), ])
raw <- merge(raw, np, by = "t_utc", all.x = TRUE)
stopifnot(!anyNA(raw$np_cons))

# Cross-checks between the two files
ok <- !is.na(raw$np_price)
num$priceMatch <- f(100 * mean(abs(raw$NO1_da_price[ok] - raw$np_price[ok]) < 1e-6), 1)
cons_ratio <- raw$np_cons / raw$NO1_cons
num$consRatioMin <- f(min(cons_ratio), 3); num$consRatioMax <- f(max(cons_ratio), 3)
num$corConsBuy <- f(cor(raw$NO1_cons, raw$NO1_buy), 3)
num$meanCons <- f(mean(raw$NO1_cons)); num$meanBuy <- f(mean(raw$NO1_buy))
num$meanNpCons <- f(mean(raw$np_cons))
cat("Price identical where raw price present (%):", num$priceMatch, "\n")
cat("Raw-file consumption / csv consumption, range:", range(cons_ratio), "\n")

# ------------------------------------------------------------------
# 2. Local time and factors
# ------------------------------------------------------------------
# Prices and demand follow the Norwegian clock, so the factors are built in
# local time (Europe/Oslo). The data itself stays on the UTC axis, which has no
# duplicated or skipped hours.
raw[, t_local := as.POSIXct(format(t_utc, tz = "Europe/Oslo"), tz = "Europe/Oslo")]
raw[, `:=`(
  hour    = as.integer(format(t_utc, "%H", tz = "Europe/Oslo")),
  month   = as.integer(format(t_utc, "%m", tz = "Europe/Oslo")),
  wday    = as.integer(format(t_utc, "%u", tz = "Europe/Oslo")),   # 1 = Monday
  year    = as.integer(format(t_utc, "%Y", tz = "Europe/Oslo")),
  isoweek = format(t_utc, "%G-W%V", tz = "Europe/Oslo"),
  date    = as.Date(format(t_utc, "%Y-%m-%d", tz = "Europe/Oslo"))
)]

# ------------------------------------------------------------------
# 3. Reservoir: weekly NVE values -> ISO week
# ------------------------------------------------------------------
# The column NO1_Fyllingsgrad_NVE is a weekly series that was linearly
# interpolated to hours in the thesis data. The weekly observations are the
# points where the hourly slope changes ("knots"). I take the value at the knot
# as the observed weekly fill level and give it to every hour of the same ISO
# week (Monday to Sunday, local time).
# The thesis interpolation was done row by row (not by clock time), so the
# slope is taken per row. This matters only around the two missing hours.
r <- raw$NO1_Fyllingsgrad_NVE
slope <- diff(r)
knot <- which(abs(diff(slope)) > 1e-8) + 1
knots <- raw[knot, .(isoweek, t_utc, wday, reservoir_obs = 100 * NO1_Fyllingsgrad_NVE)]
stopifnot(!anyDuplicated(knots$isoweek))
num$nKnots <- f(nrow(knots))
num$knotWeekday <- paste(unique(format(knots$t_utc, "%A", tz = "Europe/Oslo")), collapse = ", ")
cat("Weekly reservoir observations:", nrow(knots), " weekdays:",
    num$knotWeekday, "\n")
raw <- merge(raw, knots[, .(isoweek, reservoir = reservoir_obs)], by = "isoweek", all.x = TRUE)
setorder(raw, t_utc)
# hours after the last weekly value were only copied forward in the csv
num$nAfterLastKnot <- f(sum(raw$t_utc > max(knots$t_utc)))

# ------------------------------------------------------------------
# 4. Remove rows (counted step by step)
# ------------------------------------------------------------------
steps <- data.table(Step = "Hourly rows in the thesis file", Removed = NA_integer_,
                    Left = nrow(raw))
drop_step <- function(dt, keep, label) {
  removed <- sum(!keep)
  steps <<- rbind(steps, data.table(Step = label, Removed = removed,
                                    Left = nrow(dt) - removed))
  dt[keep]
}
d <- drop_step(raw, !is.na(raw$reservoir),
               "No weekly reservoir value in the ISO week")
num$nNoReservoir <- f(steps[.N, Removed])
d <- drop_step(d, !is.na(d$np_price),
               "Price missing in raw Nord Pool file (forward filled in csv)")
num$nMissingPrice <- f(steps[.N, Removed])
num$nMissingPriceAll <- f(sum(is.na(raw$np_price)))
mp_year <- raw[is.na(np_price), .N, by = year][order(-N)]
num$missPriceTopYear <- mp_year$year[1]; num$missPriceTopN <- f(mp_year$N[1])
mp_run <- rle(is.na(raw$np_price)); mp_run <- mp_run$lengths[mp_run$values]
num$missPriceRuns <- f(length(mp_run)); num$missPriceMaxRun <- f(max(mp_run))
num$missPriceRunsShort <- f(100 * mean(mp_run <= 2), 0)
d <- drop_step(d, d$NO1_cons > 0 & d$NO1_da_price > -500 & d$NO1_da_price < 3000,
               "Impossible values (consumption $\\le 0$, price outside market limits)")

d[, `:=`(price = NO1_da_price, cons = NO1_cons / 1000)]   # consumption in GWh per hour
month_lab <- month.abb
wday_lab <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
d[, `:=`(
  fhour  = factor(hour, levels = 0:23),
  fmonth = factor(month_lab[month], levels = month_lab),
  fwday  = factor(wday_lab[wday], levels = wday_lab),
  fyear  = factor(year)
)]
n <- nrow(d)
num$nClean <- f(n)
num$pctRemoved <- f(100 * (nrow(raw) - n) / nrow(raw), 1)
num$cleanStart <- format(min(d$t_local), "%d %B %Y")
num$cleanEnd <- format(max(d$t_local), "%d %B %Y")
num$nNeg <- f(sum(d$price < 0)); num$nZero <- f(sum(d$price == 0))
num$nNonPos <- f(sum(d$price <= 0))
num$minPrice <- f(min(d$price), 2); num$maxPrice <- f(max(d$price), 2)
num$maxPriceTime <- format(d$t_local[which.max(d$price)], "%d %B %Y at %H:00")
num$nAboveTwoHundred <- f(sum(d$price > 200))
num$nAboveTwoHundredPct <- f(100 * mean(d$price > 200), 1)
num$nYearTwentyTwo <- f(sum(d$year == 2022))
cat("Clean rows:", n, "\n")

steps[, Removed := ifelse(is.na(Removed), "", f(Removed))]
steps[, Left := f(Left)]
write_tab(as.data.frame(steps), "tab_cleaning.tex", "lrr")

fwrite(d[, .(time_utc = format(t_utc, "%Y-%m-%d %H:%M"),
             time_local = format(t_local, "%Y-%m-%d %H:%M"),
             isoweek, price_eur_mwh = price, consumption_mwh = NO1_cons,
             reservoir_pct = reservoir, hour, month = fmonth, weekday = fwday,
             year)],
       "data/processed/no1_clean.csv")

# ------------------------------------------------------------------
# 5. Descriptive statistics
# ------------------------------------------------------------------
desc <- function(x, d = 1) c(f(mean(x), d), f(sd(x), d), f(min(x), d),
                             f(quantile(x, .25), d), f(median(x), d),
                             f(quantile(x, .75), d), f(max(x), d))
tab_sum <- data.frame(
  Variable = c("Price (EUR/MWh)", "Consumption (MWh)", "Reservoir fill (\\%)"),
  rbind(desc(d$price, 1), desc(d$NO1_cons, 0), desc(d$reservoir, 1)),
  check.names = FALSE)
names(tab_sum)[-1] <- c("Mean", "SD", "Min", "Q1", "Median", "Q3", "Max")
write_tab(tab_sum, "tab_summary.tex")
num$meanPrice <- f(mean(d$price), 1); num$medianPrice <- f(median(d$price), 1)
num$sdPrice <- f(sd(d$price), 1)
num$meanConsClean <- f(mean(d$NO1_cons)); num$sdCons <- f(sd(d$NO1_cons))
num$minConsClean <- f(min(d$NO1_cons)); num$maxConsClean <- f(max(d$NO1_cons))
num$minRes <- f(min(d$reservoir), 1); num$maxRes <- f(max(d$reservoir), 1)

by_year <- d[, .(Hours = f(.N), Price = f(mean(price), 1), SDp = f(sd(price), 1),
                 Cons = f(mean(NO1_cons)), Res = f(mean(reservoir), 1)), by = year]
setnames(by_year, c("Year", "Hours", "Mean price", "SD price",
                    "Mean consumption", "Mean reservoir (\\%)"))
write_tab(as.data.frame(by_year), "tab_by_year.tex")
yr_mean <- d[, mean(price), by = year]
for (i in seq_len(nrow(yr_mean))) num[[paste0("meanPriceY", as.roman(yr_mean$year[i] - 2015))]] <- f(yr_mean$V1[i], 1)
yr_res <- d[, mean(reservoir), by = year]
for (i in seq_len(nrow(yr_res))) num[[paste0("meanResY", as.roman(yr_res$year[i] - 2015))]] <- f(yr_res$V1[i], 1)

h_mean <- d[, mean(price), by = hour][order(hour)]
num$hourMinMean <- f(min(h_mean$V1), 1); num$hourMinAt <- h_mean$hour[which.min(h_mean$V1)]
num$hourMaxMean <- f(max(h_mean$V1), 1); num$hourMaxAt <- h_mean$hour[which.max(h_mean$V1)]
m_mean <- d[, mean(price), by = month][order(month)]
num$monthMinMean <- f(min(m_mean$V1), 1); num$monthMinAt <- month.name[which.min(m_mean$V1)]
num$monthMaxMean <- f(max(m_mean$V1), 1); num$monthMaxAt <- month.name[which.max(m_mean$V1)]

cors <- cor(d[, .(price, cons, reservoir, year)])
num$corPriceCons <- f(cors["price", "cons"], 2)
num$corPriceRes <- f(cors["price", "reservoir"], 2)
num$corConsRes <- f(cors["cons", "reservoir"], 2)
num$corYearRes <- f(cors["year", "reservoir"], 2)
num$corYearCons <- f(cors["year", "cons"], 2)
num$corYearPrice <- f(cors["year", "price"], 2)
tab_cor <- data.frame(c("Price", "Consumption", "Reservoir", "Year"),
                      apply(round(cors, 2), 2, function(x) f(x, 2)))
names(tab_cor) <- c("", "Price", "Consumption", "Reservoir", "Year")
write_tab(tab_cor, "tab_cor.tex")

# ------------------------------------------------------------------
# 6. Exploratory figures
# ------------------------------------------------------------------
daily <- d[, .(price = mean(price)), by = date]
png_open("fig_price_time.png", 9, 4)
par(mar = c(4, 4.5, 1, 1))
plot(d$t_local, d$price, type = "l", col = "grey70", lwd = 0.4,
     xlab = "Time", ylab = "Day-ahead price (EUR/MWh)")
lines(as.POSIXct(daily$date) + 12 * 3600, daily$price, col = "black", lwd = 0.8)
legend("topleft", c("Hourly price", "Daily mean"), col = c("grey70", "black"),
       lwd = c(1, 1.5), bty = "n")
dev.off()

png_open("fig_price_factors.png", 9, 7)
par(mfrow = c(2, 2), mar = c(4, 4.5, 2, 1))
boxplot(price ~ hour, d, outline = FALSE, col = "grey90", xlab = "Hour of day (local time)",
        ylab = "Price (EUR/MWh)", main = "By hour")
points(1:24, h_mean$V1, pch = 19, col = "red", cex = 0.7)
boxplot(price ~ fmonth, d, outline = FALSE, col = "grey90", xlab = "Month",
        ylab = "Price (EUR/MWh)", main = "By month", las = 2)
points(1:12, m_mean$V1, pch = 19, col = "red", cex = 0.7)
boxplot(price ~ fwday, d, outline = FALSE, col = "grey90", xlab = "Weekday",
        ylab = "Price (EUR/MWh)", main = "By weekday")
points(1:7, d[, mean(price), by = fwday][order(fwday)]$V1, pch = 19, col = "red", cex = 0.7)
boxplot(price ~ fyear, d, outline = FALSE, col = "grey90", xlab = "Year",
        ylab = "Price (EUR/MWh)", main = "By year")
points(seq_along(levels(d$fyear)), d[, mean(price), by = fyear][order(fyear)]$V1,
       pch = 19, col = "red", cex = 0.7)
dev.off()

binned <- function(x, y, nb = 40) {
  b <- cut(x, quantile(x, seq(0, 1, length.out = nb + 1)), include.lowest = TRUE)
  data.table(x = tapply(x, b, mean), y = tapply(y, b, mean))
}
png_open("fig_price_vs_x.png", 9, 4.2)
par(mfrow = c(1, 2), mar = c(4, 4.5, 1, 1))
set.seed(201)
s <- sample(n, 10000)
plot(d$cons[s], d$price[s], pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.15),
     xlab = "Consumption (GWh per hour)", ylab = "Price (EUR/MWh)", ylim = c(0, 250))
b <- binned(d$cons, d$price); lines(b$x, b$y, col = "red", lwd = 2)
plot(d$reservoir[s], d$price[s], pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.15),
     xlab = "Reservoir fill (%)", ylab = "Price (EUR/MWh)", ylim = c(0, 250))
b <- binned(d$reservoir, d$price); lines(b$x, b$y, col = "red", lwd = 2)
dev.off()
num$nScatterSample <- f(length(s))
num$nAboveYlim <- f(sum(d$price > 250))


# ------------------------------------------------------------------
# 7. Gaussian GLM (identity link)
# ------------------------------------------------------------------
# y_i = x_i^T beta + e_i,  e_i ~ N(0, sigma^2)
form_full <- price ~ fhour + fmonth + fwday + fyear + cons + reservoir
m_full <- glm(form_full, family = gaussian(link = "identity"), data = d)
s_full <- summary(m_full)
r2 <- function(m) 1 - deviance(m) / m$null.deviance
rse <- function(m) sqrt(deviance(m) / df.residual(m))
num$nPar <- f(length(coef(m_full)))
num$rTwo <- f(r2(m_full), 3); num$rTwoPct <- f(100 * r2(m_full), 1)
num$rse <- f(rse(m_full), 2)
num$dfRes <- f(df.residual(m_full))
num$nullDev <- f(m_full$null.deviance / 1e6, 2); num$resDev <- f(deviance(m_full) / 1e6, 2)
num$aicFull <- f(AIC(m_full)); num$bicFull <- f(BIC(m_full))

sink("tables/model_output.txt")
cat("=== Full Gaussian GLM ===\n"); print(s_full)
cat("R^2 =", r2(m_full), " RSE =", rse(m_full), "\n")
sink()

b <- coef(m_full); se <- sqrt(diag(vcov(m_full)))
ci <- confint.default(m_full)                # Wald intervals b +/- 1.96 SE
pv <- s_full$coefficients[, 4]
fwrite(data.table(term = names(b), estimate = b, se = se, ci_low = ci[, 1],
                  ci_high = ci[, 2], p = pv), "tables/glm_all_coefficients.csv")

hc <- b[paste0("fhour", 0:23)[-1]]
h_morning <- (1:11)[which.max(hc[1:11])]        # largest hour effect before noon
h_evening <- (12:23)[which.max(hc[12:23])]      # largest hour effect after noon
h_night <- (1:23)[which.min(hc)]
key <- c("(Intercept)", paste0("fhour", c(h_night, h_morning, h_evening)),
         "fmonthApr", "fmonthJul", "fmonthDec", "fwdaySat", "fwdaySun",
         "fyear2018", "fyear2020", "fyear2021", "fyear2022", "cons", "reservoir")
lab <- c("(Intercept)", sprintf("Hour %02d (vs 00)", c(h_night, h_morning, h_evening)),
         "April (vs January)", "July (vs January)", "December (vs January)",
         "Saturday (vs Monday)", "Sunday (vs Monday)",
         "2018 (vs 2016)", "2020 (vs 2016)", "2021 (vs 2016)", "2022 (vs 2016)",
         "Consumption (per GWh/h)", "Reservoir (per \\%-point)")
tab_coef <- data.frame(Term = lab, Estimate = f(b[key], 2), SE = f(se[key], 2),
                       CI = paste0("[", f(ci[key, 1], 2), ", ", f(ci[key, 2], 2), "]"),
                       p = fp(pv[key]), check.names = FALSE)
names(tab_coef) <- c("Term", "Estimate", "SE", "95\\% CI", "$p$-value")
write_tab(tab_coef, "tab_coef.tex", "lrrrr")

nm <- function(term) gsub("[^A-Za-z]", "", term)
cn <- c("Intercept", "HourNight", "HourMorning", "HourEvening", "Apr", "Jul", "Dec",
        "Sat", "Sun", "Yeighteen", "Ytwenty", "Ytwentyone", "Ytwentytwo", "Cons", "Res")
for (i in seq_along(key)) {
  num[[paste0("b", cn[i])]] <- f(b[key[i]], 2)
  num[[paste0("se", cn[i])]] <- f(se[key[i]], 2)
  num[[paste0("ciL", cn[i])]] <- f(ci[key[i], 1], 2)
  num[[paste0("ciU", cn[i])]] <- f(ci[key[i], 2], 2)
}
num$hNight <- sprintf("%02d", h_night); num$hMorning <- sprintf("%02d", h_morning)
num$hEvening <- sprintf("%02d", h_evening)
num$bResTen <- f(10 * b["reservoir"], 1)
num$bResTenAbs <- f(abs(10 * b["reservoir"]), 1)
num$ciLResTen <- f(10 * ci["reservoir", 1], 1); num$ciUResTen <- f(10 * ci["reservoir", 2], 1)
num$bConsOne <- f(b["cons"], 1)
num$bConsSD <- f(b["cons"] * sd(d$cons), 1)            # effect of one SD
num$bResSD <- f(b["reservoir"] * sd(d$reservoir), 1)
num$bResSDAbs <- f(abs(b["reservoir"] * sd(d$reservoir)), 1)
num$sdConsGWh <- f(sd(d$cons), 2); num$sdResPct <- f(sd(d$reservoir), 1)
num$peakDiff <- f(b[paste0("fhour", h_morning)] - b[paste0("fhour", h_night)], 1)
num$wdayPMax <- fp(max(pv[grep("fwday", names(pv))]))
num$nHourNotSig <- f(sum(pv[grep("fhour", names(pv))] > 0.05))

# Figure: hour and month effects with 95% intervals
png_open("fig_effects.png", 9, 4.2)
par(mfrow = c(1, 2), mar = c(4, 4.5, 2, 1))
hb <- c(0, b[grep("fhour", names(b))]); hl <- c(0, ci[grep("fhour", names(b)), 1])
hu <- c(0, ci[grep("fhour", names(b)), 2])
plot(0:23, hb, pch = 19, ylim = range(hl, hu), xlab = "Hour of day (local time)",
     ylab = "Effect vs hour 00 (EUR/MWh)", main = "Hour of day")
suppressWarnings(arrows(0:23, hl, 0:23, hu, angle = 90, code = 3, length = 0.02))
abline(h = 0, lty = 2, col = "grey50")
mb <- c(0, b[grep("fmonth", names(b))]); ml <- c(0, ci[grep("fmonth", names(b)), 1])
mu <- c(0, ci[grep("fmonth", names(b)), 2])
plot(1:12, mb, pch = 19, ylim = range(ml, mu), xaxt = "n", xlab = "",
     ylab = "Effect vs January (EUR/MWh)", main = "Month (adjusted)")
axis(1, 1:12, month.abb, las = 2)
suppressWarnings(arrows(1:12, ml, 1:12, mu, angle = 90, code = 3, length = 0.02))
abline(h = 0, lty = 2, col = "grey50")
dev.off()

png_open("fig_year_effects.png", 5, 4)
par(mar = c(4, 4.5, 1, 1))
yb <- c(0, b[grep("fyear", names(b))]); yl <- c(0, ci[grep("fyear", names(b)), 1])
yu <- c(0, ci[grep("fyear", names(b)), 2]); yy <- as.integer(levels(d$fyear))
plot(yy, yb, pch = 19, ylim = range(yl, yu), xlab = "Year", ylab = "Effect vs 2016 (EUR/MWh)")
suppressWarnings(arrows(yy, yl, yy, yu, angle = 90, code = 3, length = 0.03))
abline(h = 0, lty = 2, col = "grey50")
dev.off()

# ------------------------------------------------------------------
# 8. Which predictors matter: drop-one comparison
# ------------------------------------------------------------------
terms_all <- c("fhour", "fmonth", "fwday", "fyear", "cons", "reservoir")
term_lab <- c(fhour = "Hour", fmonth = "Month", fwday = "Weekday", fyear = "Year",
              cons = "Consumption", reservoir = "Reservoir")
d1 <- drop1(m_full, test = "F")
reduced <- lapply(terms_all, function(tm) update(m_full, as.formula(paste(". ~ . -", tm))))
names(reduced) <- terms_all
tab_d1 <- data.frame(
  Term = term_lab[terms_all],
  df = d1[terms_all, "Df"],
  dDev = f((d1[terms_all, "Deviance"] - deviance(m_full)) / 1e6, 2),
  F = f(d1[terms_all, "F value"], 1),
  dAIC = f(sapply(reduced, AIC) - AIC(m_full)),
  dBIC = f(sapply(reduced, BIC) - BIC(m_full)),
  dR2 = f(100 * (r2(m_full) - sapply(reduced, r2)), 2),
  check.names = FALSE)
names(tab_d1) <- c("Term removed", "df", "Deviance increase ($10^6$)", "$F$",
                   "$\\Delta$AIC", "$\\Delta$BIC", "$\\Delta R^2$ (\\%-points)")
write_tab(tab_d1, "tab_drop1.tex", "lrrrrrr")
dr2 <- 100 * (r2(m_full) - sapply(reduced, r2))
for (tm in terms_all) num[[paste0("dRtwo", nm(term_lab[tm]))]] <- f(dr2[tm], 2)
for (tm in terms_all) num[[paste0("dAIC", nm(term_lab[tm]))]] <- f(AIC(reduced[[tm]]) - AIC(m_full))
for (tm in terms_all) num[[paste0("F", nm(term_lab[tm]))]] <- f(d1[tm, "F value"], 1)
num$minF <- f(min(d1[terms_all, "F value"]), 1)

# calendar-only model and model without both physical variables
m_cal <- update(m_full, . ~ . - cons - reservoir)
num$rTwoCal <- f(r2(m_cal), 3); num$rTwoCalPct <- f(100 * r2(m_cal), 1)
m_year <- glm(price ~ fyear, family = gaussian, data = d)
num$rTwoYearOnly <- f(100 * r2(m_year), 1)

# ------------------------------------------------------------------
# 9. Collinearity
# ------------------------------------------------------------------
v <- vif(m_full)
tab_vif <- data.frame(Term = term_lab[rownames(v)], GVIF = f(v[, 1], 2), df = v[, 2],
                      adj = f(v[, 3], 2), check.names = FALSE)
names(tab_vif) <- c("Term", "GVIF", "df", "GVIF$^{1/(2\\,\\mathrm{df})}$")
write_tab(tab_vif, "tab_vif.tex", "lrrr")
for (tm in rownames(v)) num[[paste0("vif", nm(term_lab[tm]))]] <- f(v[tm, 1], 1)
for (tm in rownames(v)) num[[paste0("vifadj", nm(term_lab[tm]))]] <- f(v[tm, 3], 2)

# auxiliary regressions: how much of reservoir / consumption is explained by
# month, by year, and by all other predictors (VIF = 1/(1 - R^2_all))
aux <- function(y, rhs) summary(lm(as.formula(paste(y, "~", rhs)), data = d))$r.squared
aux_tab <- data.frame(
  Variable = c("Consumption", "Reservoir"),
  Month = f(100 * c(aux("cons", "fmonth"), aux("reservoir", "fmonth")), 1),
  Year = f(100 * c(aux("cons", "fyear"), aux("reservoir", "fyear")), 1),
  Hour = f(100 * c(aux("cons", "fhour"), aux("reservoir", "fhour")), 1),
  All = f(100 * c(aux("cons", "fhour + fmonth + fwday + fyear + reservoir"),
                  aux("reservoir", "fhour + fmonth + fwday + fyear + cons")), 1),
  check.names = FALSE)
names(aux_tab) <- c("Variable", "Month", "Year", "Hour", "All other predictors")
write_tab(aux_tab, "tab_aux.tex", "lrrrr")
num$auxResMonth <- aux_tab$Month[2]; num$auxResYear <- aux_tab$Year[2]
num$auxResAll <- aux_tab[2, 5]; num$auxConsMonth <- aux_tab$Month[1]
num$auxConsYear <- aux_tab$Year[1]; num$auxConsAll <- aux_tab[1, 5]

# ------------------------------------------------------------------
# 10. Residual diagnostics
# ------------------------------------------------------------------
e <- residuals(m_full, type = "response"); fit <- fitted(m_full)
rs <- rstandard(m_full)
num$pctOutThree <- f(100 * mean(abs(rs) > 3), 2)
num$pctOutThreeNormal <- f(100 * 2 * pnorm(-3), 2)
num$nOutThree <- f(sum(abs(rs) > 3))
num$maxResid <- f(max(e), 1)
num$skewResid <- f(mean(((e - mean(e)) / sd(e))^3), 2)
num$kurtResid <- f(mean(((e - mean(e)) / sd(e))^4), 1)
res_sd_year <- d[, .(sd = sd(e[.I])), by = year]
for (i in seq_len(nrow(res_sd_year)))
  num[[paste0("resSDY", as.roman(res_sd_year$year[i] - 2015))]] <- f(res_sd_year$sd[i], 1)

png_open("fig_diagnostics.png", 9, 7.5)
par(mfrow = c(2, 2), mar = c(4, 4.5, 2, 1))
plot(fit, e, pch = 16, cex = 0.25, col = rgb(0, 0, 0, 0.1), xlab = "Fitted values",
     ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "grey50")
lines(lowess(fit[s], e[s]), col = "red", lwd = 2)
qq <- qqnorm(rs, plot.it = FALSE)
plot(qq$x, qq$y, pch = 16, cex = 0.3, xlab = "Theoretical quantiles",
     ylab = "Standardized residuals", main = "Normal Q-Q")
abline(0, 1, lty = 2, col = "red")
plot(fit, sqrt(abs(rs)), pch = 16, cex = 0.25, col = rgb(0, 0, 0, 0.1),
     xlab = "Fitted values", ylab = expression(sqrt("|Standardized residuals|")),
     main = "Scale-Location")
lines(lowess(fit[s], sqrt(abs(rs))[s]), col = "red", lwd = 2)
plot(d$t_local, e, type = "l", lwd = 0.3, xlab = "Time", ylab = "Residuals",
     main = "Residuals over time")
abline(h = 0, lty = 2, col = "red")
dev.off()

# Autocorrelation: residuals in time order. The removed hours are skipped,
# so lag k means k rows back (almost always k hours, see report).
ac <- acf(e, lag.max = 24 * 8, plot = FALSE)
pac <- pacf(e, lag.max = 24 * 8, plot = FALSE)
num$pctLagOneIsHour <- f(100 * mean(diff(as.numeric(d$t_utc)) == 3600), 1)
num$acfOne <- f(ac$acf[2], 3); num$acfTwentyfour <- f(ac$acf[25], 3)
num$acfWeek <- f(ac$acf[169], 3); num$pacfOne <- f(pac$acf[1], 3)
num$pacfTwo <- f(pac$acf[2], 3); num$pacfTwentyfour <- f(pac$acf[24], 3)
png_open("fig_acf.png", 9, 4)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
plot(ac, main = "ACF of residuals", xlab = "Lag (hours)")
abline(v = c(24, 168), lty = 3, col = "red")
plot(pac, main = "PACF of residuals", xlab = "Lag (hours)")
abline(v = c(24, 168), lty = 3, col = "red")
dev.off()

dw <- dwtest(lm(form_full, data = d))
num$dw <- f(dw$statistic, 3)
num$dwP <- fp(dw$p.value)
num$dwRho <- f(1 - dw$statistic / 2, 3)

# HAC (Newey-West) standard errors, Bartlett kernel, 168 lags (one week)
V_hac <- NeweyWest(m_full, lag = 168, prewhite = FALSE, adjust = TRUE)
se_hac <- sqrt(diag(V_hac))
key_hac <- c(paste0("fhour", c(h_morning, h_evening)), "fmonthJul", "fwdaySun",
             "fyear2020", "fyear2021", "cons", "reservoir")
lab_hac <- c(sprintf("Hour %02d (vs 00)", c(h_morning, h_evening)), "July (vs January)",
             "Sunday (vs Monday)", "2020 (vs 2016)", "2021 (vs 2016)",
             "Consumption (per GWh/h)", "Reservoir (per \\%-point)")
p_hac <- 2 * pnorm(-abs(b / se_hac))
tab_hac <- data.frame(Term = lab_hac, Est = f(b[key_hac], 2), SE = f(se[key_hac], 2),
                      SEh = f(se_hac[key_hac], 2), ratio = f(se_hac[key_hac] / se[key_hac], 1),
                      ph = fp(p_hac[key_hac]),
                      CIh = paste0("[", f(b[key_hac] - 1.96 * se_hac[key_hac], 2), ", ",
                                   f(b[key_hac] + 1.96 * se_hac[key_hac], 2), "]"),
                      check.names = FALSE)
names(tab_hac) <- c("Term", "Estimate", "SE (usual)", "SE (HAC)", "Ratio",
                    "$p$ (HAC)", "95\\% CI (HAC)")
write_tab(tab_hac, "tab_hac.tex", "lrrrrrr")
ratio <- se_hac / se
num$hacRatioMin <- f(min(ratio[key_hac]), 1); num$hacRatioMax <- f(max(ratio[key_hac]), 1)
num$seHacCons <- f(se_hac["cons"], 2); num$pHacCons <- fp(p_hac["cons"])
num$ciLHacCons <- f(b["cons"] - 1.96 * se_hac["cons"], 2)
num$ciUHacCons <- f(b["cons"] + 1.96 * se_hac["cons"], 2)
num$seHacRes <- f(se_hac["reservoir"], 3); num$pHacRes <- fp(p_hac["reservoir"])
num$ciLHacResTen <- f(10 * (b["reservoir"] - 1.96 * se_hac["reservoir"]), 1)
num$ciUHacResTen <- f(10 * (b["reservoir"] + 1.96 * se_hac["reservoir"]), 1)
num$seRes <- f(se["reservoir"], 3)
num$seHacMorning <- f(se_hac[paste0("fhour", h_morning)], 2)
num$pHacMorning <- fp(p_hac[paste0("fhour", h_morning)])
num$seHacJul <- f(se_hac["fmonthJul"], 2); num$seHacYtwentyone <- f(se_hac["fyear2021"], 2)
num$pHacSun <- fp(p_hac["fwdaySun"])

# ------------------------------------------------------------------
# 11. GAM: smooth effects of consumption and reservoir
# ------------------------------------------------------------------
# Motivated by the curved patterns in fig_price_vs_x and by the GLM residuals
# against the two variables (binned means below).
binres <- function(x) {
  bb <- cut(x, quantile(x, 0:10 / 10), include.lowest = TRUE)
  tapply(e, bb, mean)
}
br_res <- binres(d$reservoir); br_cons <- binres(d$cons)
num$binResMax <- f(max(br_res), 1); num$binResMin <- f(min(br_res), 1)
num$binConsMax <- f(max(br_cons), 1); num$binConsMin <- f(min(br_cons), 1)
png_open("fig_resid_vs_x.png", 9, 3.8)
par(mfrow = c(1, 2), mar = c(4, 4.5, 1, 1))
plot(d$cons[s], e[s], pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.15), ylim = c(-60, 60),
     xlab = "Consumption (GWh per hour)", ylab = "GLM residual (EUR/MWh)")
bb <- binned(d$cons, e); lines(bb$x, bb$y, col = "red", lwd = 2); abline(h = 0, lty = 2)
plot(d$reservoir[s], e[s], pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.15), ylim = c(-60, 60),
     xlab = "Reservoir fill (%)", ylab = "GLM residual (EUR/MWh)")
bb <- binned(d$reservoir, e); lines(bb$x, bb$y, col = "red", lwd = 2); abline(h = 0, lty = 2)
dev.off()

g_lin <- gam(form_full, data = d, method = "REML")       # same as the GLM, fitted by gam()
g <- gam(price ~ fhour + fmonth + fwday + fyear + s(cons) + s(reservoir),
         data = d, method = "REML")
sg <- summary(g)
sink("tables/model_output.txt", append = TRUE)
cat("\n=== GAM ===\n"); print(sg)
sink()
an <- anova(g_lin, g, test = "F")
set.seed(201)   # k.check uses random sampling
kc <- k.check(g)
sink("tables/model_output.txt", append = TRUE); cat("\n=== k.check ===\n"); print(kc); sink()
num$kIndexCons <- f(kc["s(cons)", "k-index"], 2); num$kIndexRes <- f(kc["s(reservoir)", "k-index"], 2)
bres <- cut(d$reservoir, quantile(d$reservoir, 0:10 / 10), include.lowest = TRUE)
top_bin <- levels(bres)[which.max(tapply(e, bres, mean))]
in_bin <- bres == top_bin
dec21b <- d$year == 2021 & d$month == 12
num$binTopLo <- f(min(d$reservoir[in_bin]), 0); num$binTopHi <- f(max(d$reservoir[in_bin]), 0)
num$binTopPctDec <- f(100 * mean(dec21b[in_bin]), 0)
num$binTopResWithDec <- f(mean(e[in_bin]), 1)
num$binTopResNoDec <- f(mean(e[in_bin & !dec21b]), 1)
print(an)
num$edfCons <- f(sg$s.table["s(cons)", "edf"], 1)
num$edfRes <- f(sg$s.table["s(reservoir)", "edf"], 1)
num$gamDevDiff <- f((deviance(g_lin) - deviance(g)) / 1e6, 2)
num$gamF <- f(an$F[2], 1); num$gamFp <- fp(an$`Pr(>F)`[2])
num$gamDf <- f(an$Df[2], 1)
r2g <- 1 - deviance(g) / m_full$null.deviance
num$rTwoGam <- f(r2g, 3); num$rTwoGamPct <- f(100 * r2g, 1)
num$rseGam <- f(sqrt(sg$scale), 2)
num$aicGam <- f(AIC(g)); num$bicGam <- f(BIC(g))
num$dAICGam <- f(AIC(m_full) - AIC(g)); num$dBICGam <- f(BIC(m_full) - BIC(g))
num$kMax <- f(9)

# predicted smooth effects at a few points (centred smooths)
pg <- function(var, at) {
  nd <- d[rep(1, length(at))]; nd[[var]] <- at
  as.numeric(predict(g, nd, type = "terms")[, paste0("s(", var, ")")])
}
bump_grid <- seq(40, 70, by = 1)
bump_at <- bump_grid[which.max(pg("reservoir", bump_grid))]
num$gamResBumpAt <- f(bump_at)
num$gamResTenToForty <- f(pg("reservoir", 10) - pg("reservoir", 40), 1)
num$gamResBumpToNinetyFive <- f(pg("reservoir", bump_at) - pg("reservoir", 95), 1)
num$gamResFortyToBump <- f(pg("reservoir", bump_at) - pg("reservoir", 40), 1)
num$gamConsTwoToSeven <- f(pg("cons", 7) - pg("cons", 2), 1)
num$gamConsSevenToEight <- f(pg("cons", 8) - pg("cons", 7), 1)
hi <- d$cons > 7
num$nConsAboveSeven <- f(sum(hi))
num$pctConsAboveSeven <- f(100 * mean(hi), 1)
num$pctConsAboveSevenCrisis <- f(100 * mean(d$t_local[hi] >= as.POSIXct("2021-09-01", tz = "Europe/Oslo")), 0)
num$meanPriceConsAboveSeven <- f(mean(d$price[hi]), 1)
crisis <- d$t_local >= as.POSIXct("2021-09-01", tz = "Europe/Oslo")
num$nCrisis <- f(sum(crisis)); num$meanPriceCrisis <- f(mean(d$price[crisis]), 1)
num$meanResCrisis <- f(mean(d$reservoir[crisis]), 1)
dec21 <- d$year == 2021 & d$month == 12
num$resDecTwentyOne <- f(mean(e[dec21]), 1)
num$resvDecTwentyOne <- f(mean(d$reservoir[dec21]), 1)
num$priceDecTwentyOne <- f(mean(d$price[dec21]), 1)
num$pctOutThreeCrisis <- f(100 * mean(crisis[abs(rs) > 3]), 0)

png_open("fig_gam.png", 9, 4)
par(mfrow = c(1, 2), mar = c(4, 4.5, 1, 1))
plot(g, select = 1, shade = TRUE, rug = FALSE, xlab = "Consumption (GWh per hour)",
     ylab = "Effect on price (EUR/MWh)", seWithMean = TRUE)
xx <- seq(min(d$cons), max(d$cons), length.out = 50)
lines(xx, b["cons"] * (xx - mean(d$cons)), col = "red", lty = 2, lwd = 2)
legend("topleft", c("GAM smooth", "GLM straight line"), col = c("black", "red"),
       lty = c(1, 2), lwd = 2, bty = "n")
plot(g, select = 2, shade = TRUE, rug = FALSE, xlab = "Reservoir fill (%)",
     ylab = "Effect on price (EUR/MWh)", seWithMean = TRUE)
xx <- seq(min(d$reservoir), max(d$reservoir), length.out = 50)
lines(xx, b["reservoir"] * (xx - mean(d$reservoir)), col = "red", lty = 2, lwd = 2)
dev.off()

# ------------------------------------------------------------------
# 12. Model comparison table
# ------------------------------------------------------------------
mods <- list("Calendar factors only" = m_cal,
             "Factors + consumption" = reduced$reservoir,
             "Factors + reservoir" = reduced$cons,
             "Full GLM" = m_full)
tab_mod <- data.frame(
  Model = c(names(mods), "GAM (smooth consumption and reservoir)"),
  p = c(sapply(mods, function(m) f(length(coef(m)))), f(sum(g$edf), 1)),
  AIC = c(sapply(mods, function(m) f(AIC(m))), f(AIC(g))),
  BIC = c(sapply(mods, function(m) f(BIC(m))), f(BIC(g))),
  R2 = c(sapply(mods, function(m) f(r2(m), 3)), f(r2g, 3)),
  RSE = c(sapply(mods, function(m) f(rse(m), 2)), f(sqrt(sg$scale), 2)),
  check.names = FALSE)
names(tab_mod) <- c("Model", "Parameters (edf)", "AIC", "BIC", "$R^2$", "RSE")
write_tab(tab_mod, "tab_models.tex", "lrrrrr")
num$edfGamTotal <- f(sum(g$edf), 1)
num$aicCal <- f(AIC(m_cal)); num$bicCal <- f(BIC(m_cal))

# ------------------------------------------------------------------
# 13. Write numbers and R citations
# ------------------------------------------------------------------
bad <- names(num)[!grepl("^[A-Za-z]+$", names(num))]
if (length(bad)) stop("bad macro names: ", paste(bad, collapse = " "))
stopifnot(!anyDuplicated(names(num)))
writeLines(sprintf("\\newcommand{\\%s}{%s}", names(num), unlist(num)), "tables/numbers.tex")

pk <- c("data.table", "hdf5r", "car", "lmtest", "sandwich", "mgcv")
cites <- list(R = citation(), datatable = citation("data.table")[1],
              hdfr = citation("hdf5r")[1], car = citation("car")[1],
              lmtest = citation("lmtest")[1], sandwich = citation("sandwich")[2],
              mgcvREML = citation("mgcv")[1], mgcvBook = citation("mgcv")[4])
bib <- unlist(lapply(names(cites), function(k) {
  x <- as.character(toBibtex(cites[[k]]))
  c(sub("\\{,$", paste0("{", k, ","), x[1]), x[-1], "")
}))
writeLines(bib, "tables/r_citations.bib")
writeLines(c(R.version.string, sapply(pk, function(p) paste(p, packageVersion(p)))),
           "tables/versions.txt")
sink("tables/session_info.txt"); print(sessionInfo()); sink()
cat("Done. Numbers written:", length(num), "\n")
