library(bayesian.modeling)
source('data.R')

parameter.names = c(trate.2000='Transmission Rate in 2002',
                        trate.2030='Transmission Rate by 2030',
                        primary.progression.risk='Risk of Primary Progression',
                        reactivation.rate='Reactivation Rate',
                        reactivation.rate.1='Reactivation Rate 1st 5 years',
                        reactivation.rate.2='Reactivation Rate after 5 years',
                        latent.protection.factor='Latent Protection Factor',
                        self.cure.rate='Self Cure Rate',
                        tb.mortality='Excess TB Mortality Rate',
                        treatment.initiation.rate='Rate of Treatment Initiation',
                        treatment.success.proportion='Treatment Success Proportion',
                        treatment.naive.success.proportion='Initial Treatment Success Proportion',
                        repeat.treatment.success.proportion='Retreatment Success Proportion',
                        relapse.risk='Risk of Relapse after Treatment',
                        relapse.rate='Rate of Relapse'
)

parameter.estimates = c(trate.2000=12,
                        trate.2030=9,
                        primary.progression.risk=0.14,
                        reactivation.rate=0.001,
                        reactivation.rate.1=0.00485,
                        reactivation.rate.2=0.0005,
                        latent.protection.factor=0.5,
                        self.cure.rate=0.13,
                        tb.mortality=0.15,
                        treatment.initiation.rate=0.56,
                        treatment.success.proportion=0.736,
                        treatment.naive.success.proportion=0.74,
                        repeat.treatment.success.proportion=0.65,
                        relapse.risk=0.04,
                        relapse.rate=1.5
                        )

parameters.are.proportions = rep(F, length(parameter.estimates))
names(parameters.are.proportions) = names(parameter.estimates)
parameters.are.proportions[c('primary.progression.risk', 'latent.protection.factor',
                             'treatment.success.proportion', 'treatment.naive.success.proportion',
                             'repeat.treatment.success.proportion', 'relapse.risk')] = T

parameter.lowers = c(trate.2000=0,
                     trate.2030=0,
                     primary.progression.risk=NA,
                     reactivation.rate=0.0005,
                     reactivation.rate.1=NA,
                     reactivation.rate.2=0.00024,
                     latent.protection.factor=0.1,
                     self.cure.rate=0.09,
                     tb.mortality=0.08,
                     treatment.initiation.rate=0.29,
                     treatment.success.proportion=NA,
                     treatment.naive.success.proportion=NA,
                     repeat.treatment.success.proportion=NA,
                     relapse.risk=0.026,
                     relapse.rate=0.9
)

parameter.uppers = c(trate.2000=25,
                     trate.2030=20,
                     primary.progression.risk=NA,
                     reactivation.rate=0.002,
                     reactivation.rate.1=NA,
                     reactivation.rate.2=0.00089,
                     latent.protection.factor=0.9,
                     self.cure.rate=0.2,
                     tb.mortality=0.3,
                     treatment.initiation.rate=1.08,
                     treatment.success.proportion=NA,
                     treatment.naive.success.proportion=NA,
                     repeat.treatment.success.proportion=NA,
                     relapse.risk=0.06,
                     relapse.rate=2.5
)

parameter.bounds.were.unspecified = is.na(parameter.lowers) | is.na(parameter.uppers)

parameter.lowers[is.na(parameter.lowers)] = 0.75 * parameter.estimates[is.na(parameter.lowers)]
parameter.lowers[parameters.are.proportions] = sapply(parameter.lowers[parameters.are.proportions], max, 0)
#proportions can't be less than zero

parameter.uppers[is.na(parameter.uppers)] = 1.25 * parameter.estimates[is.na(parameter.uppers)]
parameter.uppers[parameters.are.proportions] = sapply(parameter.uppers[parameters.are.proportions], min, 1)
    #proportions can't be greater than 1


model.1.parameter.names = c(
    'trate.2000',
    'trate.2030',
    'primary.progression.risk',
    'reactivation.rate',
    'latent.protection.factor',
    'self.cure.rate',
    'tb.mortality',
    'treatment.initiation.rate',
    'treatment.success.proportion',
    'relapse.risk'
)
NUM.PARAMS.1 = length(model.1.parameter.names)

model.2.parameter.names = c(
    'trate.2000',
    'trate.2030',
    'primary.progression.risk',
    'reactivation.rate.1',
    'reactivation.rate.2',
    'latent.protection.factor',
    'self.cure.rate',
    'tb.mortality',
    'treatment.initiation.rate',
    'treatment.success.proportion',
    'relapse.risk'
)
NUM.PARAMS.2 = length(model.2.parameter.names)

model.3.parameter.names = c(
    'trate.2000',
    'trate.2030',
    'primary.progression.risk',
    'reactivation.rate.1',
    'reactivation.rate.2',
    'latent.protection.factor',
    'self.cure.rate',
    'tb.mortality',
    'treatment.initiation.rate',
    'treatment.naive.success.proportion',
    'repeat.treatment.success.proportion',
    'relapse.risk',
    'relapse.rate'
)
NUM.PARAMS.3 = length(model.3.parameter.names)

all.model.param.names = list(model.1.parameter.names, model.2.parameter.names, model.3.parameter.names)

prepare.parameters <- function(parameters)
{
    parameters = c(parameters,
                   mortality.rate=0.0146,
                   ltbi.1.length=5
                   )
    parameters = as.list(parameters)

    trate.2000 = parameters$trate.2000
    trate.2030 = parameters$trate.2030

    START.TRATE.YEAR = 2001
    END.TRATE.YEAR = 2030
    parameters$transmission.rate = function(year)
    {
        if (year < START.TRATE.YEAR)
            trate.2000
        else if (year > END.TRATE.YEAR)
            trate.2030
        else
            trate.2000 + (trate.2030 - trate.2000) / (END.TRATE.YEAR - START.TRATE.YEAR) * (year - START.TRATE.YEAR)
    }

    parameters
}

INTERVENTION.INIT.RATE.MULT = 1.33
modify.parameters.for.intervention <- function(parameters)
{
    parameters$treatment.initiation.rate = parameters$treatment.initiation.rate * INTERVENTION.INIT.RATE.MULT

    parameters
}

NARROW_PRIOR_FACTOR = 0.25 #amount to reduce sampling interval on each end
generate.uniform.params <- function(qq, model.no, narrow=F)
{
    model.names = all.model.param.names[[model.no]]
    names(qq) = model.names

    lowers = parameter.lowers[model.names]
    uppers = parameter.uppers[model.names]

    if (narrow)
    {
        intervals = uppers - lowers
        lowers = lowers + intervals * NARROW_PRIOR_FACTOR
        uppers = uppers - intervals * NARROW_PRIOR_FACTOR
    }
    
    params = qunif(qq, lowers, uppers)
    names(params) = model.names

    params['trate.2030'] = qunif(qq['trate.2030'], lowers['trate.2030'], params['trate.2000'])

    if (any(model.names=='repeat.treatment.success.proportion'))
        params['repeat.treatment.success.proportion'] = qunif(qq['repeat.treatment.success.proportion'], lowers['repeat.treatment.success.proportion'], params['treatment.naive.success.proportion'])

    params
}

generate.uniform.params.narrow <- function(qq, model.no)
{
    generate.uniform.params(qq, model.no, narrow=T)
}


log.and.logit.norm.parameters = lapply(names(parameter.estimates), function(param.name){
    if (parameters.are.proportions[param.name])
    {
        if (parameter.bounds.were.unspecified[param.name])
            get.logitnorm.parameters.for.interval(parameter.lowers[param.name], parameter.uppers[param.name], mean=logit(parameter.estimates[param.name]))
        else
            get.logitnorm.parameters.for.interval(parameter.lowers[param.name], parameter.uppers[param.name], mean=NA)
    }
    else
    {
        if (parameter.bounds.were.unspecified[param.name] || param.name=='trate.2000' || param.name=='trate.2030')
            get.lognorm.parameters.for.interval(parameter.lowers[param.name], parameter.uppers[param.name], mean=log(parameter.estimates[param.name]))
        else
            get.lognorm.parameters.for.interval(parameter.lowers[param.name], parameter.uppers[param.name], mean=NA)
    }
})
names(log.and.logit.norm.parameters) = names(parameter.estimates)

log.and.logit.norm.parameters.narrow = lapply(names(parameter.estimates), function(param.name){
    lowers = parameter.lowers
    uppers = parameter.uppers
    
    intervals = uppers - lowers
    lowers = lowers + intervals * NARROW_PRIOR_FACTOR
    uppers = uppers - intervals * NARROW_PRIOR_FACTOR

    if (parameters.are.proportions[param.name])
    {
        if (parameter.bounds.were.unspecified[param.name])
            get.logitnorm.parameters.for.interval(lowers[param.name], uppers[param.name], mean=logit(parameter.estimates[param.name]))
        else
            get.logitnorm.parameters.for.interval(lowers[param.name], uppers[param.name], mean=NA)
    }
    else
    {
        if (parameter.bounds.were.unspecified[param.name] || param.name=='trate.2000' || param.name=='trate.2030')
            get.lognorm.parameters.for.interval(lowers[param.name], uppers[param.name], mean=log(parameter.estimates[param.name]))
        else
            get.lognorm.parameters.for.interval(lowers[param.name], uppers[param.name], mean=NA)
    }
})
names(log.and.logit.norm.parameters.narrow) = names(parameter.estimates)

generate.log.and.logit.norm.parameters <- function(qq, model.no, narrow=F)
{
    model.names = all.model.param.names[[model.no]]
    names(qq) = model.names
    
    all.hyperparameters = if (narrow) log.and.logit.norm.parameters.narrow else log.and.logit.norm.parameters

    params = sapply(model.names, function(param.name){
        hyperparams = all.hyperparameters[[param.name]]
        if (parameters.are.proportions[param.name])
            qlogitnorm(qq[param.name], hyperparams$mean, hyperparams$sd)
        else
            qlnorm(qq[param.name], hyperparams$mean, hyperparams$sd)
    })
    names(params) = model.names

    hyperparams.trate.2030 = all.hyperparameters[['trate.2030']]
    params['trate.2030'] = qtrunclnorm(qq['trate.2030'], -Inf, log(params['trate.2000']), hyperparams.trate.2030$mean, hyperparams.trate.2030$sd)

    hyperparams.rx = all.hyperparameters[['repeat.treatment.success.proportion']]
    if (any(model.names=='repeat.treatment.success.proportion'))
        params['repeat.treatment.success.proportion'] = qtrunclogitnorm(qq['repeat.treatment.success.proportion'], -Inf, logit(params['treatment.naive.success.proportion']), hyperparams.rx$mean, hyperparams.rx$sd)

    params
}

generate.log.and.logit.norm.parameters.narrow <- function(qq, model.no)
{
    generate.log.and.logit.norm.parameters(qq, model.no, narrow=T)
}

library(truncnorm)
qtrunclnorm <- function(p, alog, blog, meanlog, sdlog)
{
    q.norm = qtruncnorm(p, alog, blog, meanlog, sdlog)
    exp(q.norm)
}

qtrunclogitnorm <- function(p, a, b, mean, sd)
{
    q.norm = qtruncnorm(p, a, b, mean, sd)
    invlogit(q.norm)
}

#testing
if (1==2)
{
    qq=runif(length(model.1.parameter.names))
    generate.log.and.logit.norm.parameters(qq, 3)
}
