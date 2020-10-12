# Packages
install.packages("~/Dropbox/Guillem/Cours/Mémoire/revivals",
                 repos = NULL,
                 type="source")
library(revivals)

##### Data #####
rm(list = ls())
mucensus99 <- read.csv("rec99htegne.csv")
head(mucensus99, 5)


##### LOGVAC #####
sum(mucensus99$LOGVAC)
summary(mucensus99$LOGVAC); var(mucensus99$LOGVAC)
skewness(mucensus99$LOGVAC); kurtosis(mucensus99$LOGVAC)
plot(density(mucensus99$LOGVAC), xlab="LOGVAC", main="")


##### Scatter plot LOG / LOGVAC #####
library(ggplot2)
ggplot(mucensus99, aes(x=LOG, y=LOGVAC)) + geom_point() + geom_smooth(method="lm")
cor(mucensus99$LOGVAC, mucensus99$LOG)

##### Illustration subsection #####
library(sampling); set.seed(1906)
N = nrow(mucensus99); n = 80 # population size and sample size

  # SRSWOR
pii <- rep(n/N, n); s <- srswor(n, N)
ech <- mucensus99[s==1,]; ech$piks <- pii

tb11 <- wrapper(data = ech,
                varname = c("LOGVAC"),
                gn = N,
                est_type = c("BHR", "standard", "DT"),
                method = "si",
                pii = ech$piks,
                id = "CODE_N")

  # poisson
pii <- inclusionprobabilities(mucensus99$LOG, n); s <- UPpoisson(pii); pii <- pii[s==1]
ech <- mucensus99[s==1,]; ech$piks <- pii

tb12 <- wrapper(data = ech,
                varname = c("LOGVAC"),
                gn = N,
                est_type = c("BHR", "standard", "DT"),
                method = "poisson",
                pii = ech$piks,
                id = "CODE_N")

  # rejective
pii <- inclusionprobabilities(mucensus99$LOG, n); s <- UPmaxentropy(pii); pii <- pii[s==1]
ech <- mucensus99[s==1,]; ech$piks <- pii

tb13 <- wrapper(data = ech,
                varname = c("LOGVAC"),
                gn = N,
                est_type = c("BHR", "standard", "DT"),
                method = "rejective",
                pii = ech$piks,
                id = "CODE_N")


# Individual functions
  # Conditional bias estimation
bi <- HTcondbiasest(data = ech,
                    varname = c("LOGVAC"),
                    gn = N,
                    method = c("si"),
                    id = "CODE_N",
                    pii = ech$piks,
                    di = 1 / ech$piks,
                    remerge = T)
View(bi)

# Robust estimation using conditional bias
robustest(data = ech,
          varname = c("LOGVAC"),
          gn = N,
          method = c("rejective"),
          pii = ech$piks)

# Associated tuning constant
tuningconst(bi = bi$condbiasLOGVAC)

# Application to winsorised estimators
ww <- robustweights.r(data = ech,
                      varname = c("LOGVAC"),
                      gn = N,
                      method = c("si"),
                      pii = ech$piks,
                      typewin = "BHR",
                      remerge = F)

# Winsorisation constant
determinconstws(pii = ech$piks, 
                x = ech$LOGVAC,
                bi = bi$condbiasLOGVAC,
                tailleseq = 1000000)

determinconstwDT(pii = ech$piks, 
                 x = ech$LOGVAC, 
                 bi = bi$condbiasLOGVAC, 
                 tailleseq = 1000000)






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

# Histograms


tb <- wrapper(data = ech,
              varname = c("Valeur.fonciere_num", "Surface.reelle.bati"),
              strataname = "strata",
              gnh = N_h,
              est_type = c("BHR", "standard", "DT"),
              method = "si",
              pii = ech$piks)

