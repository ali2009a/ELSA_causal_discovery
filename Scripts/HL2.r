df <- read.csv("/home/ali/Documents/SFU/Research/dementia/refinedDF_heactb.csv", header=TRUE)


m1 <- lm(target ~ 1, data=df)
m2 <- lm(target ~ 1+indager , data=df)
m3 <- lm(target ~ 1+indager+totwq10_bu_s, data=df)
m4 <- lm(target ~ 1+indager+totwq10_bu_s+scfrdm , data=df)
m5 <- lm(target ~ 1+indager+totwq10_bu_s+scfrdm+dhsex , data=df)
m6 <- lm(target ~ 1+indager+totwq10_bu_s+scfrdm+dhsex+baseMemIndex_3, data=df)
m7 <- lm(target ~ 1+indager+totwq10_bu_s+scfrdm+dhsex+baseMemIndex_3+scorg05, data=df)
m8 <- lm(target ~ 1+indager+totwq10_bu_s+scfrdm+dhsex+baseMemIndex_3+scorg05+hehelf, data=df)


m1 <- lm(target ~ 1, data=df)
m2 <- lm(target ~ 1+indager , data=df)
m2 <- lm(target ~ 1+indager+hehelf , data=df)
m3 <- lm(target ~ 1+indager+hehelf+totwq10_bu_s, data=df)
m4 <- lm(target ~ 1+indager+hehelf+totwq10_bu_s+scfrdm , data=df)
m5 <- lm(target ~ 1+indager+hehelf+totwq10_bu_s+scfrdm+dhsex , data=df)
m6 <- lm(target ~ 1+indager+hehelf+totwq10_bu_s+scfrdm+dhsex+baseMemIndex_3, data=df)
m7 <- lm(target ~ 1+indager+hehelf+totwq10_bu_s+scfrdm+dhsex+baseMemIndex_3+scorg05, data=df)


anova(m1,m2,m3,m4,m5,m6,m7,m8)


df <- read.csv("/home/ali/Documents/SFU/Research/dementia/refinedDF_scorg05.csv", header=TRUE)
m1 <- lm(target ~ 1 , data=df)
m2 <- lm(target ~ 1+indager , data=df)
m3 <- lm(target ~ 1+indager+totwq10_bu_s, data=df)
m4 <- lm(target ~ 1+indager+totwq10_bu_s+scfrdm , data=df)
m5 <- lm(target ~ 1+indager+totwq10_bu_s+scfrdm+dhsex , data=df)
m6 <- lm(target ~ 1+indager+totwq10_bu_s+scfrdm+dhsex+baseMemIndex_3 , data=df)
m7 <- lm(target ~ 1+indager+totwq10_bu_s+scfrdm+dhsex+baseMemIndex_3+heactb , data=df)
m8 <- lm(target ~ 1+indager+totwq10_bu_s+scfrdm+dhsex+baseMemIndex_3+heactb+hehelf , data=df)

anova(m1,m2,m3,m4,m5,m6,m7,m8)
