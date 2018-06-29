##########################################################
# Pairs with DAIC > 20
##########################################################

coev = read.csv('coev_conserv_v2.csv', sep=',', header = TRUE)
coev$X = NULL
par(mfrow=c(1,1))
plot(coev$conservation~log10(coev$pairs), xlab='log10(Amount of coevolving pairs)', ylab='% Conservation')
model = lm(coev$conservation~log10(coev$pairs))
abline(model)
par(mfrow=c(2,2))
plot(model)
modsum = summary(model)

r2 = modsum$adj.r.squared
p_value = modsum$coefficients[2,4]
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(p_value, digits = 2)))[2]

legend('topright', legend = rp, bty = 'n')

##########################################################
# Histogram
##########################################################

sec_struct = read.csv('motifs_patterns.csv', sep=',', header = TRUE)
sec_struct$X = NULL
par(mfrow=c(1,1))
plot(sec_struct, ylim=c(0,60000), ylab='Number of coevolving pairs', xlab='Secondary structures combinations')

ss = sum(sec_struct$X0 == "S-S")
hh = sum(sec_struct$X0 == "H-H")
hs = sum(sec_struct$X0 == "H-S")
ht = sum(sec_struct$X0 == "H-T")
st = sum(sec_struct$X0 == "S-T")
tt = sum(sec_struct$X0 == "T-T")
total = ss+hh+hs+ht+st+tt
per_ss = (ss/total)
per_hh = (hh/total)
per_hs = (hs/total)
per_ht = (ht/total)
per_st = (st/total)
per_tt = (tt/total)

perc_pairs = c(per_ss, per_hs, per_hh, per_st, per_ht, per_tt)
barplot(perc_pairs, ylim=c(0,0.4), ylab='Percentage of coevolving pairs', xlab='Secondary structures combinations', names=c("S-S", "H-S", "H-H", "S-T", "H-T", "T-T"))

##########################################################
# Predicted vs. validation
##########################################################

pred_valid = read.csv("pred_vs_valid.csv", sep=",")
pred_valid$X = NULL
pred = pred_valid$predictions
valid = pred_valid$DAIC

plot(pred, valid, main="Expected against predicted scores", xlab="Predicted scores", ylab="Expected scores", pch=".")
