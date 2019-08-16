df <- read.csv("/home/ali/Documents/SFU/Research/dementia/refinedDF_heactb.csv", header=TRUE)


m1 <- lm(target ~ 1, data=df)
m2 <- lm(target ~ 1+hehelf , data=df)
m3 <- lm(target ~ 1+hehelf+indager , data=df)
m4 <- lm(target ~ 1+hehelf+indager+totwq10_bu_s, data=df)
m5 <- lm(target ~ 1+hehelf+indager+totwq10_bu_s+scfrdm , data=df)
m6 <- lm(target ~ 1+hehelf+indager+totwq10_bu_s+scfrdm+dhsex , data=df)
m7 <- lm(target ~ 1+hehelf+indager+totwq10_bu_s+scfrdm+dhsex+memIndex, data=df)
m8 <- lm(target ~ 1+hehelf+indager+totwq10_bu_s+scfrdm+dhsex+memIndex+scorg05, data=df)

anova(m1,m2,m3,m4,m5,m6,m7,m8)


df <- read.csv("/home/ali/Documents/SFU/Research/dementia/refinedDF_scorg05.csv", header=TRUE)
m1 <- lm(target ~ 1 , data=df)
m2 <- lm(target ~ 1+hehelf , data=df)
m3 <- lm(target ~ 1+hehelf+indager , data=df)
m4 <- lm(target ~ 1+hehelf+indager+totwq10_bu_s, data=df)
m5 <- lm(target ~ 1+hehelf+indager+totwq10_bu_s+scfrdm , data=df)
m6 <- lm(target ~ 1+hehelf+indager+totwq10_bu_s+scfrdm+dhsex , data=df)
m7 <- lm(target ~ 1+hehelf+indager+totwq10_bu_s+scfrdm+dhsex+memIndex , data=df)
m8 <- lm(target ~ 1+hehelf+indager+totwq10_bu_s+scfrdm+dhsex+memIndex+heactb , data=df)

anova(m1,m2,m3,m4,m5,m6,m7,m8)
