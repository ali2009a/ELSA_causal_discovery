df <- read.csv(file="C:\\Users\\aarab\\Rep\\dementia\\refinedDF.csv", header=TRUE, sep=",")


m1 <- lm(target ~ scorg05_1+scorg05_2+scorg05_3 , data=df)
m2 <- lm(target ~ scorg05_1+scorg05_2+scorg05_3+hehelf , data=df)
m3 <- lm(target ~ scorg05_1+scorg05_2+scorg05_3+hehelf+indager , data=df)
m4 <- lm(target ~ scorg05_1+scorg05_2+scorg05_3+hehelf+indager+totwq10_bu_s, data=df)
m5 <- lm(target ~ scorg05_1+scorg05_2+scorg05_3+hehelf+indager+totwq10_bu_s+scfrdm , data=df)
m6 <- lm(target ~ scorg05_1+scorg05_2+scorg05_3+hehelf+indager+totwq10_bu_s+scfrdm+dhsex , data=df)
m7 <- lm(target ~ scorg05_1+scorg05_2+scorg05_3+hehelf+indager+totwq10_bu_s+scfrdm+dhsex+memIndex , data=df)

anova(m1,m2,m3,m4,m5,m6,m7)





df <- read.csv(file="C:\\Users\\aarab\\Rep\\dementia\\refinedDF.csv", header=TRUE, sep=",")


m1 <- lm(target ~ heactb_1+heactb_2+heactb_3 , data=df)
m2 <- lm(target ~ heactb_1+heactb_2+heactb_3+hehelf , data=df)
m3 <- lm(target ~ heactb_1+heactb_2+heactb_3+hehelf+indager , data=df)
m4 <- lm(target ~ heactb_1+heactb_2+heactb_3+hehelf+indager+totwq10_bu_s, data=df)
m5 <- lm(target ~ heactb_1+heactb_2+heactb_3+hehelf+indager+totwq10_bu_s+scfrdm , data=df)
m6 <- lm(target ~ heactb_1+heactb_2+heactb_3+hehelf+indager+totwq10_bu_s+scfrdm+dhsex , data=df)
m7 <- lm(target ~ heactb_1+heactb_2+heactb_3+hehelf+indager+totwq10_bu_s+scfrdm+dhsex+memIndex , data=df)

anova(m1,m2,m3,m4,m5,m6,m7)
