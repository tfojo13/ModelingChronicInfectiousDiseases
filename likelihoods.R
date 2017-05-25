library(mvtnorm)
source('data.R')

#independend gaussians
independent.gaussian.incidence.likelihood <- function(results, years)
{
    real.data = get.india.data(years)
    model.prediction = get.model.incidence(results, years)
    
    prod(sapply(1:length(years), function(y){
        dnorm(real.data[y,'Incidence'], model.prediction[y], real.data[y,'Incidence.SD'])
    }))
}

independent.gaussian.mort.inc.likelihood <- function(results, years)
{
    real.inc.data = get.india.data(years, incidence = T)
    model.inc = get.model.incidence(results, years)
    
    real.mort.data = get.india.data(years, incidence = F)
    model.mort = get.model.mortality(results, years)
    
    prod(sapply(1:length(years), function(y){
        dnorm(real.inc.data[y,'Incidence'], model.inc[y], real.inc.data[y,'Incidence.SD'])
    })) *
        prod(sapply(1:length(years), function(y){
            dnorm(real.mort.data[y,'Mortality'], model.mort[y], real.mort.data[y,'Mortality.SD'])
        }))
}

mvg.ar.mort.inc.likelihood <- function(results, years, rho)
{
    real.inc.data = get.india.data(years, incidence = T)
    model.inc = get.model.incidence(results, years)
    
    real.mort.data = get.india.data(years, incidence = F)
    model.mort = get.model.mortality(results, years)
    
    observed.inc = real.inc.data[,'Incidence']
    observed.mort = real.mort.data[, 'Mortality']
    
    sds.inc = real.inc.data[,'Incidence.SD']
    sds.mort = real.mort.data[,'Mortality.SD']
    
    sigma.inc = create.autogregressive.covariance.matrix(sds.inc, rho)
    sigma.mort = create.autogregressive.covariance.matrix(sds.mort, rho)
    
    dmvnorm(observed.inc, mean=model.inc, sigma=sigma.inc) *
        dmvnorm(observed.mort, mean=model.mort, sigma=sigma.mort)
}

mvg.ar.incidence.likelihood <- function(results, years, rho)
{
    real.data = get.india.data(years)
    model.prediction = get.model.incidence(results, years)

    observed.inc = real.data[,'Incidence']
    sds = real.data[,'Incidence.SD']
    sigma = create.autogregressive.covariance.matrix(sds, rho)
    
    dmvnorm(observed.inc, mean=model.prediction, sigma=sigma)
}



mvg.cs.mort.inc.likelihood <- function(results, years, rho)
{
    real.inc.data = get.india.data(years, incidence = T)
    model.inc = get.model.incidence(results, years)
    
    real.mort.data = get.india.data(years, incidence = F)
    model.mort = get.model.mortality(results, years)
    
    observed.inc = real.inc.data[,'Incidence']
    observed.mort = real.mort.data[, 'Mortality']
    
    sds.inc = real.inc.data[,'Incidence.SD']
    sds.mort = real.mort.data[,'Mortality.SD']
    
    sigma.inc = create.compound.symmetry.covariance.matrix(sds.inc, rho)
    sigma.mort = create.compound.symmetry.covariance.matrix(sds.mort, rho)
    
    dmvnorm(observed.inc, mean=model.inc, sigma=sigma.inc) *
        dmvnorm(observed.mort, mean=model.mort, sigma=sigma.mort)
}

mvg.cs.incidence.likelihood <- function(results, years, rho)
{
    real.data = get.india.data(years)
    model.prediction = get.model.incidence(results, years)
    
    observed.inc = real.data[,'Incidence']
    sds = real.data[,'Incidence.SD']
    sigma = create.compound.symmetry.covariance.matrix(sds, rho)
    
    dmvnorm(observed.inc, mean=model.prediction, sigma=sigma)
    
}

uniform.incidence.likelihood <- function(results, years)
{
    real.data = get.india.data(years)
    model.prediction = get.model.incidence(results, years)
 
    observed.lower = real.data[,'Incidence.CI.Lower']
    observed.upper = real.data[,'Incidence.CI.Upper']
    
    all(model.prediction >= observed.lower & model.prediction <= observed.upper)
}

uniform.mort.inc.likelihood <- function(results, years)
{
    real.inc.data = get.india.data(years, incidence = T)
    model.inc = get.model.incidence(results, years)
    
    real.mort.data = get.india.data(years, incidence = F)
    model.mort = get.model.mortality(results, years)
    
    observed.inc.lower = real.inc.data[,'Incidence.CI.Lower']
    observed.inc.upper = real.inc.data[,'Incidence.CI.Upper']

    observed.mort.lower = real.mort.data[,'Mortality.CI.Lower']
    observed.mort.upper = real.mort.data[,'Mortality.CI.Upper']
    
    all(model.inc >= observed.inc.lower & model.inc <= observed.inc.upper & 
            model.mort >= observed.mort.lower & model.mort <= observed.mort.upper)
}


independent.binomial.incidence.likelihood <- function(results, years)
{
    real.data = get.india.data(years, per.population=1, just.estimates = T)
    model.prediction = get.model.incidence(results, years, per.population = 1)
    pop = get.india.population(years)
    real.counts = round(real.data*pop)
    
    prod(sapply(1:length(years), function(y){
        dbinom(real.counts[y], pop[y], model.prediction[y])
    }))
}

##-- MULTIVARIATE GAUSSIANS --##


create.autogregressive.covariance.matrix <- function(sds, rho)
{
    n = length(sds)
    exponents = c((n-1):1,0:(n-1))
    coeffs = rho^exponents
    
    corr.mat = sapply(1:n, function(i){
        index = n - i + 1
        coeffs[index:(index+n-1)]
    })
    
    cov.mat = corr.mat * rep(sds, each=n) * rep(sds, n)
    
    cov.mat
}

create.compound.symmetry.covariance.matrix <- function(sds, rho)
{
    n = length(sds)
    
    corr.mat = matrix(rho, ncol=n, nrow=n)
    diag(corr.mat) = 1
    
    cov.mat = corr.mat * rep(sds, each=n) * rep(sds, n)
    
    cov.mat
}

#for testing
make.mvplot <- function(sigma, means=rep(0,dim(sigma)[1]), n=100)
{
    rands = rmvnorm(n, means, sigma)
    dimnames(rands)[[2]] = paste0('x', 1:length(means))
    
    ggplot(as.data.frame(rands), aes(x=x1, y=x2)) + geom_point(size=1)
}

##-- OPTIMIZATION --##
optimize.incidence <- function(model.no=1)
{
    fn = function(trates)
    {
        params = parameter.estimates
        params$trate.2000 = trates[1]
        params$trate.2030 = trates[2]
        params = prepare.parameters(params)
        
        rr = run.model(model.no, params)
        
        lik = gaussian.incidence.likelihood(c(2000,2015), rr)
        
        print(paste0('trate.2000 = ', trates[1], ', trate.2030 = ', trates[2], ' -> likelihood = ', lik))
        
        lik
    }
    
    lower = c(0,0)
    upper = c(15,10)
    
    precision=c(.1,.1)
    
#    maximize.brute.force(lower, upper, precision, fn)
    maximize.no.local.max(c(8,6), lower, upper, precision, fn)
}