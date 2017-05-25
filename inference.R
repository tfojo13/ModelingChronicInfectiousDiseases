source('models.R')

get.incidence.output <- function(rr, params, years, in.intervention)
{
    rv = get.model.incidence(rr, years, in.intervention=in.intervention)
    names(rv) = as.character(years)
    rv
}


get.mortality.output <- function(rr, params, years, in.intervention)
{
    rv = get.model.mortality(rr, years, in.intervention=in.intervention)
    names(rv) = as.character(years)
    rv
}

get.param.output <- function(rr, params)
{
    params
}

get.change.output <- function(rr, params, inc=T, end.year=END.YEAR)
{
    prefix = if (inc) 'inc-change-' else 'mort-change-'
    
    rates.noint = get.model.rates(rr, c(INTERVENTION.YEAR-1, end.year), in.intervention = F, incidence=inc)
    rates.int = get.model.rates(rr, c(INTERVENTION.YEAR-1, end.year), in.intervention = T, incidence=inc)
    
    noint.change = -(rates.noint[2] - rates.noint[1]) / rates.noint[1]
    int.change = -(rates.int[2] - rates.int[1]) / rates.int[1]
    
    c(noint.change=noint.change, 
      int.change=int.change, 
      int.change.over.noint.change=int.change/noint.change,
      abs.reduction=(rates.noint[2]-rates.int[2]),
      rel.reduction=((rates.noint[2]-rates.int[2])/rates.noint[2]))
}

get.main.outcome <- function(rr, params, inc=T, end.year=END.YEAR)
{
    get.change.metrics(rr, params, inc, end.year)['difference.in.2030.rate.relative']
}

get.change.metrics <- function(rr, params, inc=T, end.year=END.YEAR)
{
    rate.pre = get.model.rates(rr, INTERVENTION.YEAR-1, in.intervention = F, incidence=inc)
    rate.noint.post = get.model.rates(rr, end.year, in.intervention=F, incidence=inc)
    rate.int.post = get.model.rates(rr, end.year, in.intervention=T, incidence=inc)
    num.years = end.year - INTERVENTION.YEAR + 1
    
    names(rate.pre) = names(rate.noint.post) = names(rate.int.post) = NULL
    
    cum.rate.noint = sum(get.model.rates(rr, INTERVENTION.YEAR:end.year, in.intervention=F, incidence=inc))
    cum.rate.int = sum(get.model.rates(rr, INTERVENTION.YEAR:end.year, in.intervention=T, incidence=inc))
    
    rv = c(change.per.year.no.intervention.absolute = (rate.noint.post - rate.pre) / num.years,
           change.per.year.yes.intervention.absolute = (rate.int.post - rate.pre) / num.years,
           change.per.year.no.intervention.relative = (rate.noint.post - rate.pre) / num.years / rate.pre,
           change.per.year.yes.intervention.relative = (rate.int.post - rate.pre) / num.years / rate.pre,
           difference.in.2030.rate.absolute = rate.noint.post - rate.int.post,
           difference.in.2030.rate.relative = (rate.noint.post - rate.int.post) / rate.noint.post,
           rate.in.2030.no.intervention = rate.noint.post,
           rate.in.2030.yes.intervention = rate.int.post,
           cases.2018.to.2030.no.intervention = cum.rate.noint,
           cases.2018.to.2030.yes.intervention = cum.rate.int,
           difference.in.cases.to.2030.absolute = cum.rate.noint - cum.rate.int,
           difference.in.cases.to.2030.relative = (cum.rate.noint - cum.rate.int) / cum.rate.noint
      )
    
    
    rv
}

get.paper.metrics <- function(rr, params)
{
    inc.2030.noint = get.model.incidence(rr, 2030, in.intervention=F)
    inc.2030.int = get.model.incidence(rr, 2030, in.intervention=T)

    inc.2025.noint = get.model.incidence(rr, 2025, in.intervention=F)
    inc.2025.int = get.model.incidence(rr, 2025, in.intervention=T)
    
    mort.2030.noint = get.model.mortality(rr, 2030, in.intervention=F)
    mort.2030.int = get.model.mortality(rr, 2030, in.intervention=T)
    
    mort.2025.noint = get.model.mortality(rr, 2025, in.intervention=F)
    mort.2025.int = get.model.mortality(rr, 2025, in.intervention=T)
    
    rv = c(rel.inc.diff.2030 = (inc.2030.int-inc.2030.noint)/inc.2030.noint,
           rel.inc.diff.2025 = (inc.2025.int-inc.2025.noint)/inc.2025.noint,
           rel.mort.diff.2030 = (mort.2030.int-mort.2030.noint)/mort.2030.noint,
           rel.mort.diff.2025 = (mort.2025.int-mort.2025.noint)/mort.2025.noint,
           inc.2030.noint=inc.2030.noint,
           inc.2030.int=inc.2030.int,
           inc.2025.noint=inc.2025.noint,
           inc.2025.int=inc.2025.int,
           mort.2030.noint=mort.2030.noint,
           mort.2030.int=mort.2030.int,
           mort.2025.noint=mort.2025.noint,
           mort.2025.int=mort.2025.int)
}

get.linreg.output <- function(rr, params, inc=T, years=2000:2015)
{
    do.get.linreg(years, get.model.rates(rr, years, incidence=inc))
}

do.get.linreg <- function(years, estimates)
{
    years = years-years[1]
    fit = lm(estimates ~ years)
    rv = fit$coefficients
    names(rv) = c('intercept', 'slope')
    rv
}

cov.to.corr <- function(cov)
{
    n = dim(cov)[1]
    cov / rep(sqrt(diag(cov)), n) / rep(sqrt(diag(cov)), each=n)
}

linreg.matrix <- function(years, relative.to.first.year=T)
{
    if (relative.to.first.year)
        years = years-years[1]
    
    X = cbind(rep(1,length(years)), years)
    
    M=solve(t(X)%*%X)%*%t(X)
    dimnames(M)[[1]] = c('intercept', 'slope')
    M
}

PRE.YEARS = 2000:(INTERVENTION.YEAR-1)
POST.YEARS = (INTERVENTION.YEAR):2030
get.outputs <- function(resamples)
{
    resamples = load.resamples(resamples)

    rv = list()

    rv$resamples = resamples

    rv$pre.inc.dist = perform.inference.on.resample(resamples, get.incidence.output, years=PRE.YEARS, in.intervention=F)
    rv$pre.mort.dist = perform.inference.on.resample(resamples, get.mortality.output, years=PRE.YEARS, in.intervention=F)

    rv$noint.inc.dist = perform.inference.on.resample(resamples, get.incidence.output, years=POST.YEARS, in.intervention=F)
    rv$noint.mort.dist = perform.inference.on.resample(resamples, get.mortality.output, years=POST.YEARS, in.intervention=F)

    rv$int.inc.dist = perform.inference.on.resample(resamples, get.incidence.output, years=POST.YEARS, in.intervention=T)
    rv$int.mort.dist = perform.inference.on.resample(resamples, get.mortality.output, years=POST.YEARS, in.intervention=T)
    
    rv$inc.change.dist = perform.inference.on.resample(resamples, get.change.output, inc=T)
    rv$mort.change.dist = perform.inference.on.resample(resamples, get.change.output, inc=F)
    
    rv$inc.linreg = perform.inference.on.resample(resamples, get.linreg.output, inc=T)
    rv$mort.linreg = perform.inference.on.resample(resamples, get.linreg.output, inc=F)
    
    rv$param.dist = perform.inference.on.resample(resamples, get.param.output)
    
    rv$inc.change.metrics = perform.inference.on.resample(resamples, get.change.metrics, inc=T)
    rv$mort.change.metrics = perform.inference.on.resample(resamples, get.change.metrics, inc=T)
    
    rv$paper.metrics = perform.inference.on.resample(resamples, get.paper.metrics)
    
    rv$pre.dists = list(rv$pre.inc.dist, rv$pre.mort.dist)
    rv$noint.dists = list(rv$noint.inc.dist, rv$noint.mort.dist)
    rv$int.dists = list(rv$int.inc.dist, rv$int.mort.dist)
    rv$change.dists = list(rv$inc.change.dist, rv$mort.change.dist)
    rv$linreg.dists = list(rv$inc.linreg, rv$mort.linreg)
    
    rv
}

resample.and.get.output <- function(directory, likelihood, ...)
{
    res = resample.distributed.samples(directory, NUM.RESAMPLES, likelihood, ...)
    get.outputs(res)
}

get.quantile.matrix <- function(resamples)
{
    
}