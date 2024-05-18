#Trennschärfen- und Reliabilitätsanalyse

library(psych)
library(MBESS)
library(car)

# Itemanalyse
alpha(x = TAP)

# Trennschärfe wird in der ausgabe als raw.r bezeichnet. Es ist die Korrelation des Items mit dem Mittelwert der restlichen
# items. Bei dieser Trennschärfe wird das item, das mit dem Mittelwert der restlichen items korreliert wird, 
# nicht aus dem Mittelwert entfernt oder dafür kontrolliert. Hier handelt es sich um eine partielle Eigenkorrelation,
# da das item, das korreliert wird, auch im Mittelwert enthalten ist.
# Es überschätzt je nach Konstellation (Höhe der Trennschärfe und Anzahl der Items) die wahre Trennschärfe.
# 
# Die Trennschärfe r.drop bezeichnet und gibt die Korrelation eines items mit dem Mittelwert der restlichen items
# einer latenten variablen an und entspricht der trennschärfe, die in spss angegeben wird.
# Diese Trennschärfe hat den Nachteil, je nachdem welches item weggelassen wurde, unterscheidet sich die Reliabilität
# der restlichen items.

# Reliabilität
MBESS::ci.reliability(TAP, type = "omega", interval.type = "bca", B = 10000, conf.level = 0.95)

# Itemmittelwert

# korrelieren des item MWs mit dem item. 
cor.test(TAP$item1, item1mw)

# Überprüfung der Linearität

#item1 tilde item1mw
# Datendatei = TAP
# regressionslinie reg.line=lm



