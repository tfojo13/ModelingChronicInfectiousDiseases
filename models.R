library(deSolve)

##------------##
##-- STATES --##
##------------##

#For all models

CUMULATIVE.INCIDENCE = 1
CUMULATIVE.MORTALITY = 2

#For models 1 and 2

UNINFECTED = 3
ACTIVE.TB = 4

#For model 1 only

LTBI = 5
NUM.STATES.1 = 5
STATENAMES.1 = c('Incidence', 'Mortality', 'Uninfected', 'Active TB', 'LTBI')

#For model 2 only

LTBI.1 = 5
LTBI.2 = 6
NUM.STATES.2 = 6
STATENAMES.2 = c('Incidence', 'Mortality', 'Uninfected', 'Active TB', 'LTBI.1', 'LTBI.2')

#For model 3
# *suffix .N = never treated, .P = previously treated

UNINFECTED.N = 3
ACTIVE.TB.N = 4
LTBI.1.N = 5
LTBI.2.N = 6
WILL.RELAPSE = 7
UNINFECTED.P = 8
ACTIVE.TB.P = 9
LTBI.1.P = 10
LTBI.2.P = 11

NUM.STATES.3 = 11

LTBI.1.LENGTH = 5


##--------------------------##
##-- DX Functions for ODE --##
##--------------------------##

dx.model.1 <- function(t, y, parameters)
{
    total.pop = sum(y[-c(CUMULATIVE.INCIDENCE, CUMULATIVE.MORTALITY)])
    frac.active = y[ACTIVE.TB]/total.pop
    force.of.infection = frac.active * parameters$transmission.rate(t)

    uninfected.to.ltbi = y[UNINFECTED] * force.of.infection * (1-parameters$primary.progression.risk)
    uninfected.to.active = y[UNINFECTED] * force.of.infection * parameters$primary.progression.risk

    ltbi.to.active = y[LTBI] * (parameters$reactivation.rate +
                                parameters$latent.protection.factor * force.of.infection * parameters$primary.progression.risk)

    active.to.uninfected = y[ACTIVE.TB] * parameters$treatment.initiation.rate * parameters$treatment.success.proportion * (1 - parameters$relapse.risk)
    active.to.ltbi = y[ACTIVE.TB] * parameters$self.cure.rate

    mortality.in.tb = y[ACTIVE.TB] * (parameters$mortality.rate + parameters$tb.mortality)
    total.mortality = total.pop * parameters$mortality.rate + y[ACTIVE.TB] * parameters$tb.mortality


    rv = numeric(NUM.STATES.1)
    rv[CUMULATIVE.INCIDENCE] = ltbi.to.active + uninfected.to.active +
                                y[ACTIVE.TB] * parameters$treatment.initiation.rate * parameters$treatment.success.proportion * parameters$relapse.risk

    rv[CUMULATIVE.MORTALITY] = mortality.in.tb
    rv[UNINFECTED] = active.to.uninfected - uninfected.to.ltbi - uninfected.to.active -
        y[UNINFECTED] * parameters$mortality.rate + total.mortality
    rv[ACTIVE.TB] = uninfected.to.active + ltbi.to.active - active.to.uninfected - active.to.ltbi - mortality.in.tb
    rv[LTBI] = uninfected.to.ltbi + active.to.ltbi - ltbi.to.active - y[LTBI] * parameters$mortality.rate

    list(rv)
}

dx.model.2 <- function(t, y, parameters)
{
    total.pop = sum(y[-c(CUMULATIVE.INCIDENCE, CUMULATIVE.MORTALITY)])
    frac.active = y[ACTIVE.TB]/total.pop
    force.of.infection = frac.active * parameters$transmission.rate(t)

    uninfected.to.ltbi.1 = y[UNINFECTED] * force.of.infection * (1-parameters$primary.progression.risk)
    uninfected.to.active = y[UNINFECTED] * force.of.infection * parameters$primary.progression.risk

    ltbi.1.to.active = y[LTBI.1] * (parameters$reactivation.rate.1 +
                                    parameters$latent.protection.factor * force.of.infection * parameters$primary.progression.risk)
    ltbi.2.to.active = y[LTBI.2] * (parameters$reactivation.rate.2 +
                                        parameters$latent.protection.factor * force.of.infection * parameters$primary.progression.risk)
    ltbi.1.to.2 = y[LTBI.1] / parameters$ltbi.1.length
    ltbi.2.to.1 = y[LTBI.2] * parameters$latent.protection.factor * force.of.infection * (1 - parameters$primary.progression.risk)

    active.to.uninfected = y[ACTIVE.TB] * parameters$treatment.initiation.rate * parameters$treatment.success.proportion * (1-parameters$relapse.risk)
    active.to.ltbi.1 = y[ACTIVE.TB] * parameters$self.cure.rate

    mortality.in.tb = y[ACTIVE.TB] * (parameters$mortality.rate + parameters$tb.mortality)
    total.mortality = total.pop * parameters$mortality.rate + y[ACTIVE.TB] * parameters$tb.mortality

    rv = numeric(NUM.STATES.2)
    rv[CUMULATIVE.INCIDENCE] = ltbi.1.to.active + ltbi.2.to.active + uninfected.to.active +
                                y[ACTIVE.TB] * parameters$treatment.initiation.rate * parameters$treatment.success.proportion * parameters$relapse.risk
    rv[CUMULATIVE.MORTALITY] = mortality.in.tb
    rv[UNINFECTED] = active.to.uninfected - uninfected.to.ltbi.1 - uninfected.to.active -
        y[UNINFECTED] * parameters$mortality.rate + total.mortality
    rv[ACTIVE.TB] = uninfected.to.active + ltbi.1.to.active + ltbi.2.to.active - active.to.uninfected - active.to.ltbi.1 -
        mortality.in.tb
    rv[LTBI.1] = uninfected.to.ltbi.1 + ltbi.2.to.1 + active.to.ltbi.1 - ltbi.1.to.2 - ltbi.1.to.active - y[LTBI.1] * parameters$mortality.rate
    rv[LTBI.2] = ltbi.1.to.2 - ltbi.2.to.1 - ltbi.2.to.active - y[LTBI.2] * parameters$mortality.rate

    list(rv)
}

dx.model.3 <- function(t, y, parameters)
{
    total.pop = sum(y[-c(CUMULATIVE.INCIDENCE, CUMULATIVE.MORTALITY)])
    frac.active = y[ACTIVE.TB]/total.pop
    force.of.infection = frac.active * parameters$transmission.rate(t)

    #FROM UNINFECTED.N
    uninfected.n.to.ltbi.1.n = y[UNINFECTED.N] * force.of.infection * (1-parameters$primary.progression.risk)
    uninfected.n.to.active.n = y[UNINFECTED.N] * force.of.infection * parameters$primary.progression.risk

    #FROM LTBI.N
    ltbi.1.n.to.active.n = y[LTBI.1.N] * (parameters$reactivation.rate.1 +
                                        parameters$latent.protection.factor * force.of.infection * parameters$primary.progression.risk)
    ltbi.2.n.to.active.n = y[LTBI.2.N] * (parameters$reactivation.rate.2 +
                                        parameters$latent.protection.factor * force.of.infection * parameters$primary.progression.risk)
    ltbi.1.n.to.2.n = y[LTBI.1.N] / parameters$ltbi.1.length
    ltbi.2.n.to.1.n = y[LTBI.2.N] * parameters$latent.protection.factor * force.of.infection * (1 - parameters$primary.progression.risk)

    #FROM ACTIVE.N
    active.n.to.ltbi.1.n = y[ACTIVE.TB.N] * parameters$self.cure.rate
    active.n.to.relapse = y[ACTIVE.TB.N] * parameters$treatment.initiation.rate * parameters$treatment.naive.success.proportion * parameters$relapse.risk
    active.n.to.active.p = y[ACTIVE.TB.N] * parameters$treatment.initiation.rate * (1 - parameters$treatment.naive.success.proportion)
    active.n.to.uninfected.p = y[ACTIVE.TB.N] * parameters$treatment.initiation.rate * parameters$treatment.naive.success.proportion * (1 - parameters$relapse.risk)


    #FROM UNINFECTED.P
    uninfected.p.to.ltbi.1.p = y[UNINFECTED.P] * parameters$latent.protection.factor * force.of.infection * (1-parameters$primary.progression.risk)
    uninfected.p.to.active.p = y[UNINFECTED.P] * parameters$latent.protection.factor * force.of.infection * parameters$primary.progression.risk

    #FROM LTBI.P
    ltbi.1.p.to.active.p = y[LTBI.1.P] * (parameters$reactivation.rate.1 +
                                              parameters$latent.protection.factor * force.of.infection * parameters$primary.progression.risk)
    ltbi.2.p.to.active.p = y[LTBI.2.P] * (parameters$reactivation.rate.2 +
                                              parameters$latent.protection.factor * force.of.infection * parameters$primary.progression.risk)
    ltbi.1.p.to.2.p = y[LTBI.1.P] / parameters$ltbi.1.length
    ltbi.2.p.to.1.p = y[LTBI.2.P] * parameters$latent.protection.factor * force.of.infection * (1 - parameters$primary.progression.risk)

    #FROM ACTIVE.P
    active.p.to.ltbi.1.p = y[ACTIVE.TB.P] * parameters$self.cure.rate
    active.p.to.relapse = y[ACTIVE.TB.P] * parameters$treatment.initiation.rate * parameters$repeat.treatment.success.proportion * parameters$relapse.risk
    active.p.to.uninfected.p = y[ACTIVE.TB.P] * parameters$treatment.initiation.rate * parameters$repeat.treatment.success.proportion * (1 - parameters$relapse.risk)


    #FROM RELAPSE
    relapse.to.active.p = y[WILL.RELAPSE] * (parameters$relapse.rate + parameters$latent.protection.factor * force.of.infection * parameters$primary.progression.risk)


    #MORTALITY
    mortality.in.tb.n = y[ACTIVE.TB.N] * (parameters$mortality.rate + parameters$tb.mortality)
    mortality.in.tb.p = y[ACTIVE.TB.P] * (parameters$mortality.rate + parameters$tb.mortality)
    total.mortality = total.pop * parameters$mortality.rate + (y[ACTIVE.TB.N] + y[ACTIVE.TB.P]) * parameters$tb.mortality

    #DDx's
    rv = numeric(NUM.STATES.2)
    rv[CUMULATIVE.INCIDENCE] = ltbi.1.n.to.active.n + ltbi.2.n.to.active.n + uninfected.n.to.active.n +
                                ltbi.1.p.to.active.p + ltbi.2.p.to.active.p + uninfected.p.to.active.p +
                                relapse.to.active.p
    rv[CUMULATIVE.MORTALITY] = mortality.in.tb.n + mortality.in.tb.p

    rv[UNINFECTED.N] = total.mortality - uninfected.n.to.ltbi.1.n - uninfected.n.to.active.n -
        y[UNINFECTED.N] * parameters$mortality.rate
    rv[ACTIVE.TB.N] = uninfected.n.to.active.n + ltbi.1.n.to.active.n + ltbi.2.n.to.active.n -
                    active.n.to.uninfected.p - active.n.to.active.p - active.n.to.relapse - active.n.to.ltbi.1.n -
                    mortality.in.tb.n
    rv[LTBI.1.N] = uninfected.n.to.ltbi.1.n + ltbi.2.n.to.1.n + active.n.to.ltbi.1.n -
                    ltbi.1.n.to.2.n - ltbi.1.n.to.active.n - y[LTBI.1.N] * parameters$mortality.rate
    rv[LTBI.2.N] = ltbi.1.n.to.2.n - ltbi.2.n.to.1.n - ltbi.2.n.to.active.n - y[LTBI.2.N] * parameters$mortality.rate

    rv[UNINFECTED.P] = active.n.to.uninfected.p + active.p.to.uninfected.p -
                        uninfected.p.to.ltbi.1.p - uninfected.p.to.active.p -
                        y[UNINFECTED.P] * parameters$mortality.rate
    rv[ACTIVE.TB.P] = uninfected.p.to.active.p + ltbi.1.p.to.active.p + ltbi.2.p.to.active.p +
                        active.n.to.active.p + relapse.to.active.p -
                        active.p.to.uninfected.p - active.p.to.relapse - active.p.to.ltbi.1.p -
                        mortality.in.tb.p
    rv[LTBI.1.P] = uninfected.p.to.ltbi.1.p + ltbi.2.p.to.1.p + active.p.to.ltbi.1.p -
                    ltbi.1.p.to.2.p - ltbi.1.p.to.active.p - y[LTBI.1.P] * parameters$mortality.rate
    rv[LTBI.2.P] = ltbi.1.p.to.2.p - ltbi.2.p.to.1.p - ltbi.2.p.to.active.p - y[LTBI.2.P] * parameters$mortality.rate

    rv[WILL.RELAPSE] = active.n.to.relapse + active.p.to.relapse - relapse.to.active.p - y[WILL.RELAPSE] * parameters$mortality.rate

    list(rv)
}


##-----------------------##
##-- RUNNING THE MODEL --##
##-----------------------##

END.YEAR = 2030
INTERVENTION.YEAR = 2018

run.model <- function(model.no, parameters)
{
    parameters = prepare.parameters(parameters)

    times = 1500:END.YEAR
    start.year=times[1]
    end.year=times[length(times)]
    intervention.times = (INTERVENTION.YEAR - 1):end.year

    if (model.no == 1)
        ode.output = ode(create.init.state(NUM.STATES.1), times, dx.model.1, parameters)
    else if (model.no == 2)
        ode.output = ode(create.init.state(NUM.STATES.2), times, dx.model.2, parameters)
    else
        ode.output = ode(create.init.state(NUM.STATES.3), times, dx.model.3, parameters)

    before.intervention.index = INTERVENTION.YEAR - start.year
    pre.intervention.state = ode.output[before.intervention.index,-1]
    parameters.intervention = modify.parameters.for.intervention(parameters)

    if (model.no == 1)
        ode.output2 = ode(pre.intervention.state, intervention.times, dx.model.1, parameters.intervention)
    else if (model.no == 2)
        ode.output2 = ode(pre.intervention.state, intervention.times, dx.model.2, parameters.intervention)
    else
        ode.output2 = ode(pre.intervention.state, intervention.times, dx.model.3, parameters.intervention)

    results = list(base.data=process.results(ode.output),
                   intervention.data=process.results(ode.output2),
                   years=times,
                   start.year=start.year,
                   intervention.data.start.year = intervention.times[1],
                   intervention.year = INTERVENTION.YEAR,
                   end.year=end.year)

    results$intervention.data[1, 1+c(CUMULATIVE.INCIDENCE, CUMULATIVE.MORTALITY)] = results$base.data[before.intervention.index, 1+c(CUMULATIVE.INCIDENCE, CUMULATIVE.MORTALITY)]
    
    results
}



process.results <- function(results)
{
    n = dim(results)[1]

    population.size = rowSums(results[,-c(1,1+CUMULATIVE.INCIDENCE,1+CUMULATIVE.MORTALITY)])

    split.mean.population.size = (population.size[-1] + population.size[-n]) / 2
    time.diff = results[-1,1] - results[-n,1]

    results[-1,1+CUMULATIVE.INCIDENCE] = (results[-1,1+CUMULATIVE.INCIDENCE] - results[-n,1+CUMULATIVE.INCIDENCE]) / split.mean.population.size / time.diff
    results[1, 1+CUMULATIVE.INCIDENCE] = NA

    results[-1,1+CUMULATIVE.MORTALITY] = (results[-1,1+CUMULATIVE.MORTALITY] - results[-n,1+CUMULATIVE.MORTALITY]) / split.mean.population.size / time.diff
    results[1, 1+CUMULATIVE.MORTALITY] = NA

    results[,-c(1,1+CUMULATIVE.INCIDENCE,1+CUMULATIVE.MORTALITY)] = results[,-c(1,1+CUMULATIVE.INCIDENCE,1+CUMULATIVE.MORTALITY)] / population.size

    results
}

#creates a state with 1 active case, all else uninfected
create.init.state <- function(num.states, n=100000)
{
    rv = numeric(num.states)
    rv[ACTIVE.TB] = 1
    rv[UNINFECTED] = n-1

    rv
}

get.model.incidence <- function(results, years, per.population=100000, in.intervention=F)
{
    data = if (in.intervention) results$intervention.data else results$base.data
    start.year = if (in.intervention) results$intervention.data.start.year else results$start.year

    rv = data[years - start.year + 1, CUMULATIVE.INCIDENCE+1] * per.population
    names(rv) = as.character(years)
    rv
}


get.model.mortality <- function(results, years, per.population=100000, in.intervention=F)
{
    data = if (in.intervention) results$intervention.data else results$base.data
    start.year = if (in.intervention) results$intervention.data.start.year else results$start.year

    rv = data[years - start.year + 1, CUMULATIVE.MORTALITY+1] * per.population
    names(rv) = as.character(years)
    rv
}

get.model.rates <- function(results, years, per.population=100000, in.intervention=F, incidence=T)
{
    if (incidence)
        get.model.incidence(results, years, per.population, in.intervention)
    else
        get.model.mortality(results, years, per.population, in.intervention)
}
