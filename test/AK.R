rm(list = ls())

library("sketching")

y <- AK$LWKLYWGE
intercept <- AK$CNST
X_end <- AK$EDUC
X_exg <- AK[,3:11]
#X_exg <- AK[,3:5]
x <- cbind(X_exg, X_end)
Z_inst <- AK[,12:(ncol(AK)-1)]
z <- cbind(X_exg, Z_inst)
fullsample <- cbind(y,intercept,x)
n <- nrow(fullsample)
d <- ncol(x)
x <- as.matrix(x)
z <- as.matrix(z)
n0  <- ceiling(n*0.1)

# scaling

sd_x = sd(x)
sd_y = sd(y)
x = scale(x)
y = scale(y)
z = scale(z)

# scaling factor
cns = sd_y/sd_x

# TSLS
tsls = ivreg::ivreg(y ~ x | z) 
# use heteroskedasticity-robust asymptotic variance
#tsls_test = lmtest::coeftest(tsls, df = Inf, vcov = sandwich::vcovHC, type = "HC0")
tsls_ci = lmtest:: coefci(tsls, parm = (d+1), level = 0.95, df = Inf, vcov = sandwich::vcovHC, type = "HC0")

res_tsls = cbind(tsls$coefficients[d+1], tsls_ci)


# sample split 
x0 = x[c(1:n0),]
z0 = z[c(1:n0),]
y0 = y[c(1:n0)]

x1 = x[c((n0+1):n),]
z1 = z[c((n0+1):n),]
y1 = y[c((n0+1):n)]

# SIVE: starting values
Phi_start = (t(z0)%*%x0)/n0
bt_start = solve( (t(x0)%*%z0) %*% solve(t(z0)%*%z0) %*% (t(z0)%*%x0) ) %*% ( (t(x0)%*%z0) %*% solve(t(z0)%*%z0) %*% (t(z0) %*% y0) )
w_start = solve(t(z0)%*%z0/n0)


# random permutation
nper = 200
nper = 1
m = 1
n_stop = m*(n-n0)

res = {}

for (i_p in 1:nper){

ind = sample.int((n-n0),n_stop, replace=TRUE)
x1 = x1[ind,]
z1 = z1[ind,]
y1 = y1[ind]

# SIVE
out = sgmm::sgmm(x=x1, y=y1, z=z1, gamma_0=1, alpha=0.501, bt_start = bt_start,
           inference = "rs", n0 = n0,
           Phi_start = Phi_start, w_start = w_start)

#out_mb = sive_mb(x=x1, y=y1, z=z1, gamma_0=1, alpha=0.501, bt_start = bt_start,
#                 inference = "rs", n0 = n0)

cv_rs = 6.747
est_sive = out$coefficient[d]
est_ci_lb = est_sive - cv_rs * sqrt(out$V_hat[d,d]/(n-n0))
est_ci_ub = est_sive + cv_rs * sqrt(out$V_hat[d,d]/(n-n0))

res_sive = cbind(est_sive, est_ci_lb, est_ci_ub)

#est_sive = out_mb$coefficient[d]
#est_ci_lb = est_sive - cv_rs * sqrt(out_mb$V_hat[d,d]/(n-n0))
#est_ci_ub = est_sive + cv_rs * sqrt(out_mb$V_hat[d,d]/(n-n0))

#res_sive_mb = cbind(est_sive, est_ci_lb, est_ci_ub)

# storing results

out.result = cns * cbind(res_sive)
print(i_p)
print(out.result)

res = rbind(res, out.result)
}

final_median = apply(res, 2, median)
final_mean = apply(res, 2, mean)

print(final_median)
print(final_mean)

print(cns*res_tsls)

