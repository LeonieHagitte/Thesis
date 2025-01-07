install.packages('semPower')
library(semPower)

ap <- semPower.aPriori(effect = .05, effect.measure = 'RMSEA', 
                       alpha = .05, power = .80, df = 150)
summary(ap)
#################

# obtain the required N to reject the hypothesis of metric invariance
# in comparison to the configural invariance model
# with a power of 80% on alpha = 5%
# for amodel involving a five factors (= two measurements) which
# is measured by 4 indicators
# loading by .6 each at the first measurement occasion
# loading by .7 each in the second measurement occasion and so on
# and assuming autocorrelated residuals

powerLI <- semPower.powerLI(
  type = 'a-priori', alpha = .05, power = .80,
  comparison = 'configural',
  nullEffect = 'metric',
  nIndicator = c(4, 4, 4, 4, 4),
  loadM = c(.6, .7, .6, .7, .7),
  autocorResiduals = TRUE,
  simulatedPower = TRUE
)
# show summary
summary(powerLI)

