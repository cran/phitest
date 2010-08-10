##############
# This function computes the phi-divergence test to determine whether a 
#     data set come from the uniform(0,1) distribution.
# User inputs are X (the data) and S.
##############

phi.test = function(x, s, fx="punif", ...) {

#To test with known quantiles
#S = 2
#N = 2
#STATISTIC = 100.11/N

X = x
S = s
N = length(X)
if(N < 1) {stop("data vector is empty")}
if(N==1 & S < 1) {stop("for s < 1, n must be at least 2")}
if(N > 10000) {stop("calculations are not supported for n>10,000")}

fx = get(fx, mode="function")
if (mode(fx) != "function") {stop("'fx' must be a string naming a valid function")}
Z = fx(X, ...)
if (sum(Z == 0 | Z == 1) > 0) {stop("data includes values outside the support of 'fx'")}
#print(X)
#print(Z)

STATISTIC = phi.statistic(Z, s)
names(STATISTIC) = paste("S_n(", S, ")", sep="")
PVAL = phi.p(STATISTIC, S, N)
nm_alternative = "two-sided"
METHOD = "One-sample Phi-divergence test"
DNAME = deparse(substitute(x))

RVAL = list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative, method = METHOD, data.name = DNAME)
class(RVAL) = "htest"
return(RVAL)

}


############
# K_s(u, v)
# Inputs are u, v, and the parameter S
############

Ks = function(S, u, v) {

if(S != 0 && S != 1) {
	ks = 1/(S*(1-S))*( 1 - u^S*v^(1-S) - (1-u)^S*(1-v)^(1-S) )
	return(ks)
}
if(S == 1) {
	k0 = u*log(u/v) + (1-u)*log((1-u)/(1-v))
	return(k0)
}
if(S == 0) {
	k1 = v*log(v/u) + (1-v)*log((1-v)/(1-u))
	return(k1)
}

}


############
# This function calculates the statistic for phi.test()
# Inputs are data and S.
# Need to fix this for S not in the interval (0,1)
############

phi.statistic = function(data, S) {

if(S < 1) {
n = length(data)
x = sort(data)

if(n > 2) {
i = 2:(n-1)
STATISTIC = max( c(Ks(S, 1/n, x[1]), Ks(S, (i-1)/n, x[i]), Ks(S, i/n, x[i]), Ks(S, (n-1)/n, x[n]) ) )
}
if(n == 2) {
STATISTIC = max( c(Ks(S, 1/n, x[1]), Ks(S, (n-1)/n, x[n]) ) )
}
}

if(S > 1) {
n = length(data)
x = sort(data)
i = 1:n
STATISTIC = max(c(Ks(S, (i-1)/n, x[i]), Ks(S, i/n, x[i])) )
}

if(S == 1) {
n = length(data)
x = sort(data)

if(n > 2) {
i = 2:(n-1)
STATISTIC = max( c(log(1/(1-x[1])), Ks(S, 1/n, x[1]), Ks(S, (i-1)/n, x[i]), Ks(S, i/n, x[i]), Ks(S, (n-1)/n, x[n]), log(1/x[n]) ) )
}
if(n == 2) {
STATISTIC = max( c(log(1/(1-x[1])), Ks(S, 1/n, x[1]), Ks(S, (n-1)/n, x[n]), log(1/x[n]) ) )
}
if(n < 2) {
STATISTIC = max( c(log(1/(1-x[1])), log(1/x[n]) ) )
}
}

return(STATISTIC)


}


############
# This function calculates the p-value for phi.test()
# Inputs are STATISTIC, S, and N.
# Need to fix this for S not in the interval (0,1)
############

phi.p = function(STATISTIC, S, N) {


#    if(!is.loaded("Pval_wrapper"))
#      dyn.load(paste("phitest", .Platform$dynlib.ext, sep = ""))

    scr = .C("Pval_wrapper", 
      lambda = as.double(STATISTIC),
      s = as.double(S),
      n = as.integer(N),
      scr = double(1))$scr
    
    return(scr)

}



############
# This function finds the (1-ALPHA)th confidence bands for the 
# true underlying distribution function of the data. 
# Inputs are data, S, and ALPHA.
############

phi.bands = function(x, s, alpha=.05) {

N = length(x)
S = s
ALPHA = alpha

if(N < 1) {stop("data vector is empty")}
if(N==1 & S < 1) {stop("for s < 1, n must be at least 2")}
if(N > 10000) {stop("calculations are not supported for n>10,000")}

QUANT = phi.quant(N, S, ALPHA)
b = phi.b(QUANT, S, N)
b = b[2:(N+1)]
a = rep(NA, length(b))
for(i in 1:N) {
    a[i]=1-b[N-i+1]
}
L = c(0,a)
H = c(b,1)
DATA = x
DNAME = deparse(substitute(x))

RVAL = list(l = L, h = H, data = DATA, data.name = DNAME)
class(RVAL) = "phibands"
return(RVAL)
}

#########
# This function finds the (1-ALPHA)th quantile of the null distribution of 
# S_n(s)
# 
#########

phi.quant = function(N, S, ALPHA) {

#    if(!is.loaded("Quant_wrapper"))
#      dyn.load(paste("phitest", .Platform$dynlib.ext, sep = ""))

    scr = .C("Quant_wrapper", 
      alpha = as.double(ALPHA),
      s = as.double(S),
      n = as.integer(N),
      scr = double(1))$scr
    
    return(scr)
}


#########
# This function finds the bi for calculating the confidence bands,
# given a quantile value.
# Inputs are QUANT, S, and N.
#########

phi.b = function(QUANT, S, N) {

#    if(!is.loaded("Bands_wrapper"))
#      dyn.load(paste("phitest", .Platform$dynlib.ext, sep = ""))

    scr = .C("Bands_wrapper", 
      lambda = as.double(QUANT),
      s = as.double(S),
      n = as.integer(N),
      scr = double(N+2))$scr
    
    return(scr)

}


######################
# This function plots the confidence bands with the data.
# Input is the object from phi.bands
######################

plot.phibands = function(x,...) {

    if (!inherits(x, "phibands"))
        stop("use only with \"phibands\" objects")


n = length(x$data)
pm = (range(x$data)[2] - range(x$data)[1])/10 
data.x = c(min(x$data) - pm, sort(x$data), max(x$data) + pm)

plot(data.x, c(0, (1:n)/n, 1), type="s", xlab="x", ylab="Fn", ...)
lines(data.x, c(x$l,x$l[n+1]), type="s", col="gray60")
lines(data.x, c(x$h,x$h[n+1]), type="s", col="gray60")
lines(data.x, c(0, (1:n)/n, 1), type="s")

}


