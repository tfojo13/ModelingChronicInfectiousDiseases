

#for fitting distributions to the data

get.sd.for.norm.interval <- function(mean, interval.lower, interval.upper, desired.coverage=0.95, tolerance=0.0001)
{
    get.sd.for.interval(cdf=pnorm, mean=mean, interval.lower=interval.lower, interval.upper=interval.upper, desired.coverage=desired.coverage, tolerance=tolerance)
}

library(logitnorm)
get.logitnorm.parameters.for.interval <- function(interval.lower, interval.upper, mean=NA, desired.coverage=0.95, tolerance=0.001)
{
    if (is.na(mean))
        mean = (logit(interval.upper) + logit(interval.lower)) / 2
    sd = get.sd.for.interval(plogitnorm, mean, interval.lower, interval.upper, desired.coverage, tolerance)

    list(mean=mean, sd=sd)
}

get.lognorm.parameters.for.interval <- function(interval.lower, interval.upper, mean=NA, desired.coverage=0.95, tolerance=0.001)
{
    if (is.na(mean))
        mean = (log(interval.upper) + log(max(interval.lower,0.000001))) / 2
    sd = get.sd.for.interval(plnorm, mean, interval.lower, interval.upper, desired.coverage, tolerance)

    list(mean=mean, sd=sd)
}

#calculates the sd such that the given distribution, with the mean given, has ~coverage probability
# over the interval interval.lower to interval.upper
get.sd.for.interval <- function(cdf, mean, interval.lower, interval.upper, desired.coverage=0.95, tolerance=0.001, sd.start=1)
{
    sd = max(sd.start, 0.00001)
    sd.lower.bound = NA
    sd.upper.bound = NA
    iter = 0
    max.iter = 1000
    not.found=T

    while (not.found)
    {
        iter = iter + 1
        if (iter > max.iter)
        {
            stop(paste0('Unable to find solution for mean ', mean, ' and interval [', interval.lower, ', ', interval.upper, '] after ', max.iter, ' iterations. Consider choosing different parameters'))
        }

        actual.coverage = cdf(interval.upper, mean, sd) - cdf(interval.lower, mean, sd)

        if (abs(desired.coverage-actual.coverage) <= tolerance)
            not.found = F
        else if (actual.coverage > desired.coverage) #need to increase sd
        {
            sd.lower.bound = sd
            if (is.na(sd.upper.bound))
                sd = 2 * sd
            else
                sd = (sd.lower.bound + sd.upper.bound) / 2
        }
        else #need to decrease sd
        {
            sd.upper.bound = sd
            if (is.na(sd.lower.bound))
                sd = sd / 2
            else
                sd = (sd.lower.bound + sd.upper.bound) / 2
        }
    }

    sd
}




##---------------------------------------##
##-- ACTUALLY READ AND SET UP THE DATA --##
##---------------------------------------##

INDIA.DATA = read.csv(file='WHO_India.csv')
names(INDIA.DATA) = c('Year', 'Incidence', 'Incidence.CI.Lower', 'Incidence.CI.Upper', 'Mortality', 'Mortality.CI.Lower', 'Mortality.CI.Upper', 'Population')
DATA.DENOMINATOR = 10000
INDIA.DATA$Incidence.SD = apply(INDIA.DATA, 1, function(data){
    get.sd.for.norm.interval(mean=data['Incidence'], interval.lower=data['Incidence.CI.Lower'], interval.upper=data['Incidence.CI.Upper'])
})
INDIA.DATA$Mortality.SD = apply(INDIA.DATA, 1, function(data){
    get.sd.for.norm.interval(mean=data['Mortality'], interval.lower=data['Mortality.CI.Lower'], interval.upper=data['Mortality.CI.Upper'])
})

get.india.data <- function(years, incidence=T, per.population=100000, just.estimates=F)
{
    rv = INDIA.DATA[sapply(INDIA.DATA$Year, function(y){any(y==years)}),]
    if (just.estimates)
    {
        if (incidence)
            rv[,2] * per.population / 100000
        else
            rv[,5] * per.population / 100000
    }
    else
    {
        if (incidence)
            rv[,c(1:4,9)] * per.population / 100000
        else
            rv[,c(1,5:7,10)] * per.population / 100000
    }
}

get.india.population <- function(years)
{
    INDIA.DATA[sapply(INDIA.DATA$Year, function(y){any(y==years)}),'Population']
}

