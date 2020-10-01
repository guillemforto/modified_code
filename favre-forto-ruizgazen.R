# Packages
install.packages("~/Dropbox/Guillem/Cours/Mémoire/revivals_0.0.0.9000.tar",
                 repos = NULL, 
                 type="source")
library(revivals)

##### Data #####
rm(list = ls())
rec99htegne <- read.csv("rec99htegne.csv")
head(rec99htegne, 5)


##### LOGVAC #####
sum(rec99htegne$LOGVAC)
summary(rec99htegne$LOGVAC); var(rec99htegne$LOGVAC)
skewness(rec99htegne$LOGVAC); kurtosis(rec99htegne$LOGVAC)
plot(density(rec99htegne$LOGVAC), xlab="LOGVAC", main="")


##### Scatter plot LOG / LOGVAC #####
library(ggplot2)
ggplot(rec99htegne, aes(x=LOG, y=LOGVAC)) + geom_point() + geom_smooth(method="lm")
cor(rec99htegne$LOGVAC, rec99htegne$LOG)

##### Illustration subsection #####
library(sampling); set.seed(1906)
N = nrow(rec99htegne); n = 80 # population size and sample size

  # SRSWOR
pii <- rep(n/N, n); s <- srswor(n, N)
ech <- rec99htegne[s==1,]; ech$piks <- pii

tb11 <- wrapper(data = ech,
                varname = c("LOGVAC"),
                gn = N,
                est_type = c("BHR", "standard", "DT"),
                method = "si",
                pii = ech$piks,
                id = "CODE_N")

  # poisson
pii <- inclusionprobabilities(rec99htegne$LOG, n); s <- UPpoisson(pii); pii <- pii[s==1]
ech <- rec99htegne[s==1,]; ech$piks <- pii

tb12 <- wrapper(data = ech,
                varname = c("LOGVAC"),
                gn = N,
                est_type = c("BHR", "standard", "DT"),
                method = "poisson",
                pii = ech$piks,
                id = "CODE_N")

  # rejective
pii <- inclusionprobabilities(rec99htegne$LOG, n); s <- UPmaxentropy(pii); pii <- pii[s==1]
ech <- rec99htegne[s==1,]; ech$piks <- pii

tb13 <- wrapper(data = ech,
                varname = c("LOGVAC"),
                gn = N,
                est_type = c("BHR", "standard", "DT"),
                method = "rejective",
                pii = ech$piks,
                id = "CODE_N")


# Application on a Property Values Data set (PVD)
set.seed(1906)
  # Define the population size and the sample size in each stratum (prop. allocation)
n <- 1000
N <- nrow(jeu_donnees)
N_h <- as.vector(table(jeu_donnees$strata))
nh_prop <- vector()
for (i in 1:length(N_h)){ nh_prop[i] = round(n * N_h[i] / N) }

  # STSRSWOR
pii_strata <- as.data.frame(cbind(c(1:length(nh_prop)), nh_prop / N_h))
names(pii_strata) <- c("strata", "piks")
jeu_donnees <- merge(jeu_donnees, pii_strata, by="strata")
st <- strata(jeu_donnees, stratanames=c("strata"), size=nh_prop, method="srswor")

  # Get the sample data
ech <- getdata(jeu_donnees, st)

  # Conditional bias estimation
bi <- strata_HTcondbiasest(data = ech,
                           strataname = "strata",
                           varname = c("Valeur.fonciere_num", "Surface.reelle.bati"),
                           gnh = N_h,
                           method = c("si"),
                           pii = ech$piks,
                           remerge = T)
View(bi)

# Robust estimation using conditional bias
strata_robustest(data = ech,
                 strataname = "strata",
                 varname = c("Valeur.fonciere_num", "Surface.reelle.bati"),
                 gnh = N_h,
                 method = c("si"),
                 pii = ech$piks)

# Associated tuning constant
tuningconst(bi = c(bi$Valeur.fonciere_num, bi$Surface.reelle.bati))

# Application to winsorised estimators
stww <- strata_robustweights.r(data = ech,
                               strataname = "strata", 
                               varname = c("Valeur.fonciere_num", "Surface.reelle.bati"),
                               gnh = N_h,
                               method = "si",
                               pii = ech$piks,
                               typewin = "BHR",
                               remerge = F)
stww

# Winsorisation constant


