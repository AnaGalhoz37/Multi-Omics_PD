# bio statistics 
# example for validation cohort

# nominal variables
table.sex <- data.frame(t(data.frame(CTR=c(0,2),PD=c(4,2))))
fisher.test(table.sex)

# continuous variables 
age.PD <- c(76,
            82,
            85,
            87,
            73,
            93)
age.CTR <- c(60,
             89)
t.test(age.PD,age.CTR)
PMI.PD <- c(15,
            18,
            17,
            21,
            10,
            24)
PMI.CTR <- c(13,
             13)
t.test(PMI.PD,PMI.CTR)
