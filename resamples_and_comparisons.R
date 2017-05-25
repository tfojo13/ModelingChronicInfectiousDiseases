source('priors.R')
source('models.R')
source('likelihoods.R')
source('inference.R')
source('plots.R')
NUM.SAMPLES = 10000
NUM.RESAMPLES = NUM.SAMPLES * 5
ONE.OUTCOME.TO.PLOT='difference.in.2030.rate.relative'
SOME.OUTCOMES.TO.PLOT1 = c('difference.in.2030.rate.relative', 'difference.in.2030.rate.absolute', 'rate.in.2030.no.intervention', 'change.per.year.no.intervention.relative')
OUTCOMES.TO.PLOT=NULL

TEXT_THEME = theme(text=element_text(size=20))
TEXT_THEME_SMALLER = theme(text=element_text(size=16))

priors.u = get.parameters.prior('sampling/model3')
priors.unw = get.parameters.prior('sampling/model3_unarrow')
priors.log = get.parameters.prior('sampling/model3_log')
priors.lognw = get.parameters.prior('sampling/model3_lognarrow')

plot.distributions(priors.u, priors.unw, priors.log, priors.lognw, names=c('Uniforms', 'Uniforms.Narrow', 'Logs', 'Logs.Narrow'), color=c('red','blue','orange','black'), smoothed=F)
plot.distributions(priors.u, priors.unw, priors.log, priors.lognw, names=c('Uniform', 'Uniform.Narrow', 'LogNormal', 'LogNormal.Narrow'), color=c('red','blue','orange','#008000'), smoothed=F, vars='reactivation.rate.2', geom.to.use = function(...){geom_line(...,size=2)}) + TEXT_THEME + theme(legend.position='bottom')
plot.distributions(priors.u, priors.unw, priors.log, priors.lognw, names=c('Uniform', 'Uniform.Narrow', 'LogitNormal', 'LogitNormal.Narrow'), color=c('red','blue','orange','#008000'), smoothed=F, vars='latent.protection.factor', geom.to.use = function(...){geom_line(...,size=2)}) + TEXT_THEME + theme(legend.position='bottom')

plot.four.priors <- function(u1.low, u1.high,
                             u2.low, u2.high,
                             lnorm1.mean, lnorm1.sd,
                             lnorm2.mean, lnorm2.sd,
                             l.is.logit=F,
                             n.x=1000,
                             return.df=F)
{
    d1 = function(x){dunif(x, u1.low, u1.high)}
    d2 = function(x){dunif(x, u2.low, u2.high)}
    
    d.for.l = if (l.is.logit) dlogitnorm else dlnorm
    q.for.l = if (l.is.logit) qlogitnorm else qlnorm
    name.for.l = if (l.is.logit) 'LogitNormal' else 'LogNormal'
    
    d3 = function(x){d.for.l(x,lnorm1.mean,lnorm1.sd)}
    d4 = function(x){d.for.l(x,lnorm2.mean,lnorm2.sd)}
    
    q3 = function(x){q.for.l(x,lnorm1.mean,lnorm1.sd)}
    q4 = function(x){q.for.l(x,lnorm2.mean,lnorm2.sd)}
    
    alpha=0.001
    min.x = min(u1.low, u2.low, q3(alpha), q4(alpha))
    max.x = max(u1.high, u2.high, q3(1-alpha), q4(1-alpha))
    
    x = seq(from=min.x, to=max.x, length.out=n.x)
    
    gen.df = function(d, name){
        data.frame(x=x, d=d(x), Distribution=name)
    }
    
    df = rbind(gen.df(d1, 'Uniform'),
               gen.df(d2, 'Uniform.Narrow'),
               gen.df(d3, name.for.l),
               gen.df(d4, paste0(name.for.l, '.Narrow')))
    
    if (return.df)
        df
    else
        ggplot(df, aes(x,d,color=Distribution)) + geom_line(size=2) + 
            TEXT_THEME + theme(legend.position='bottom')
}



intervals = parameter.uppers-parameter.lowers
narrow.lowers = parameter.lowers + intervals/4
narrow.uppers = parameter.uppers - intervals/4
unif.prior = Joint.Canonical.Distribution('unif', parameter.lowers['reactivation.rate.2'], parameter.uppers['reactivation.rate.2'])
unif.narrow.prior = Joint.Canonical.Distribution('unif', narrow.lowers['reactivation.rate.2'], narrow.uppers['reactivation.rate.2'])
lnorm.prior = Joint.Canonical.Distribution('lnorm', log.and.logit.norm.parameters$reactivation.rate.2$mean, log.and.logit.norm.parameters$reactivation.rate.2$sd)
lnorm.narrow.prior = Joint.Canonical.Distribution('lnorm', log.and.logit.norm.parameters.narrow$reactivation.rate.2$mean, log.and.logit.norm.parameters.narrow$reactivation.rate.2$sd)
plot.distributions(unif.prior, unif.narrow.prior, lnorm.prior, lnorm.narrow.prior, names=c('Uniform', 'Uniform.Narrow', 'LogNormal', 'LogNormal.Narrow'), color=c('red','blue','orange','#008000'), smoothed=F, geom.to.use = function(...){geom_line(...,size=2)}) + TEXT_THEME + theme(legend.position='bottom')

param='reactivation.rate.2'
plot.four.priors(parameter.lowers[param], parameter.uppers[param],
                 narrow.lowers[param], narrow.uppers[param],
                 log.and.logit.norm.parameters[[param]]$mean, log.and.logit.norm.parameters[[param]]$sd,
                 log.and.logit.norm.parameters.narrow[[param]]$mean, log.and.logit.norm.parameters.narrow[[param]]$sd) + 
    ggtitle('Reactivation Rate')
param='latent.protection.factor'
plot.four.priors(parameter.lowers[param], parameter.uppers[param],
                 narrow.lowers[param], narrow.uppers[param],
                 log.and.logit.norm.parameters[[param]]$mean, log.and.logit.norm.parameters[[param]]$sd,
                 log.and.logit.norm.parameters.narrow[[param]]$mean, log.and.logit.norm.parameters.narrow[[param]]$sd,
                 l.is.logit = T) +
    ggtitle('Latent Protection Factor')


param='latent.protection.factor'
df.lpf=plot.four.priors(parameter.lowers[param], parameter.uppers[param],
                 narrow.lowers[param], narrow.uppers[param],
                 log.and.logit.norm.parameters[[param]]$mean, log.and.logit.norm.parameters[[param]]$sd,
                 log.and.logit.norm.parameters.narrow[[param]]$mean, log.and.logit.norm.parameters.narrow[[param]]$sd,
                 l.is.logit = T, return.df=T)

df.lpf.u = df.lpf[df.lpf$Distribution=='Uniform',]
df.lpf.u.logit = rbind(df.lpf.u, df.lpf[df.lpf$Distribution=='LogitNormal',])

ggplot(df.lpf.u, aes(x,d,color=Distribution)) + geom_line(size=2) + 
    TEXT_THEME + theme(legend.position='bottom')+
    ggtitle('Latent Protection Factor') + ylim(NA,1.5)

ggplot(df.lpf.u.logit, aes(x,d,color=Distribution)) + geom_line(size=2) + 
    TEXT_THEME + theme(legend.position='bottom')+
    ggtitle('Latent Protection Factor') + ylim(NA,1.5)


plot.dist.comparison(Model.Run=r3.u.cs5, 
                     dist.name='inc.change.metrics', ncol=2,
                     vars.to.plot=SOME.OUTCOMES.TO.PLOT1) + TEXT_THEME_SMALER + theme(legend.position='bottom')

plot.multiple.results(Uniform.Priors=r3.u.cs5) + TEXT_THEME


#assume distributed sampling has been run on
#-model1
#-model2
#-model3
#-model3_log
#-model3_unarrow
#-model3_lognarrow
#as in run_distributed_2.R


#for comparison of priors
r3.u.cs5 = resample.and.get.output('sampling/model3', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
r3.unw.cs5 = resample.and.get.output('sampling/model3_unarrow', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
r3.log.cs5 = resample.and.get.output('sampling/model3_log', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
r3.lognw.cs5 = resample.and.get.output('sampling/model3_lognarrow', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)


plot.results(r3.u.cs5, incidence=T)
plot.results(r3.unw.cs5, incidence=T)
plot.results(r3.log.cs5, incidence=T)
plot.results(r3.lognw.cs5, incidence=T)

plot.dist.comparison(Uniform.Priors=r3.u.cs5, 
                     Uniform.Priors.Narrow=r3.unw.cs5,
                     Log.Priors=r3.log.cs5,
                     Log.Priors.Narrow=r3.lognw.cs5,
                     dist.name='inc.change.metrics', ncol=3,
                     vars.to.plot=OUTCOMES.TO.PLOT)  + TEXT_THEME + theme(legend.position='bottom')

plot.dist.comparison(Uniform.Priors=r3.u.cs5, 
                     Uniform.Priors.Narrow=r3.unw.cs5,
                     Log.Priors=r3.log.cs5,
                     Log.Priors.Narrow=r3.lognw.cs5,
                     dist.name='inc.change.metrics', ncol=2,
                    vars.to.plot=SOME.OUTCOMES.TO.PLOT1) + TEXT_THEME_SMALLER + theme(legend.position='bottom')


plot.multiple.results(Uniform.Priors=r3.u.cs5, 
                     Uniform.Priors.Narrow=r3.unw.cs5,
                     Log.Priors=r3.log.cs5,
                     Log.Priors.Narrow=r3.lognw.cs5) 

plot.posteriors(Uniform.Priors=r3.u.cs5, 
                Uniform.Priors.Narrow=r3.unw.cs5,
                Log.Priors=r3.log.cs5,
                Log.Priors.Narrow=r3.lognw.cs5,
                color=c('red','blue','orange','black'))


#for comparison of likelihoods
r3.u.ind2 = resample.and.get.output('sampling/model3', independent.gaussian.incidence.likelihood, years=c(2000,2015))
r3.u.indall = resample.and.get.output('sampling/model3', independent.gaussian.incidence.likelihood, years=2000:2015)
r3.u.uall = resample.and.get.output('sampling/model3', uniform.incidence.likelihood, years=2000:2015)
r3.u.cs9 = resample.and.get.output('sampling/model3', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.9)


plot.results(r3.u.cs5, incidence=T)
plot.results(r3.u.ind2, incidence=T)
plot.results(r3.u.indall, incidence=T)
plot.results(r3.u.uall, incidence=T)

plot.dist.comparison(Independent.Gaussians.16=r3.u.indall,
                     Compound.Symmetry.0.5=r3.u.cs5,
#                     Compound.Symmetry.0.9=r3.u.cs9,
                     Independent.Gaussians.2=r3.u.ind2,
                     Uniform=r3.u.uall,
                     dist.name='inc.change.metrics', ncol=3,
                     vars.to.plot=OUTCOMES.TO.PLOT)  + TEXT_THEME

plot.dist.comparison(Independent.Gaussians.16=r3.u.indall,
                     Compound.Symmetry.0.5=r3.u.cs5,
                     #                     Compound.Symmetry.0.9=r3.u.cs9,
                     Independent.Gaussians.2=r3.u.ind2,
                     Uniform=r3.u.uall,
                     dist.name='inc.change.metrics', ncol=2,
vars.to.plot=SOME.OUTCOMES.TO.PLOT1) + TEXT_THEME_SMALLER + theme(legend.position='bottom')


plot.posteriors(Independent.Gaussians.16=r3.u.indall,
                Compound.Symmetry.0.5=r3.u.cs5,
                #                     Compound.Symmetry.0.9=r3.u.cs9,
                Independent.Gaussians.2=r3.u.ind2,
                Uniform=r3.u.uall,
                color=c('red','blue','orange','black'))


plot.multiple.results(Independent.Gaussians.16=r3.u.indall,
                      Compound.Symmetry.0.5=r3.u.cs5,
#                      Compound.Symmetry.0.9=r3.u.cs9,
                      Independent.Gaussians.2=r3.u.ind2,
                      Uniform=r3.u.uall) + TEXT_THEME_SMALLER

plot.posteriors(Independent.Gaussians.16=r3.u.indall,
                Compound.Symmetry.0.5=r3.u.cs5,
                Compound.Symmetry.0.9=r3.u.cs9,
                Independent.Gaussians.2=r3.u.ind2,
                Uniform=r3.u.uall)


r3.u.cs5.m = resample.and.get.output('sampling/model3', mvg.cs.mort.inc.likelihood, years=2000:2015, rho=0.5)
r3.u.indall.m = resample.and.get.output('sampling/model3', independent.gaussian.mort.inc.likelihood, years=2000:2015)
r3.u.ind2.m = resample.and.get.output('sampling/model3', independent.gaussian.mort.inc.likelihood, years=c(2000,2015))
r3.u.uall.m = resample.and.get.output('sampling/model3', uniform.mort.inc.likelihood, years=2000:2015)

plot.multiple.results(Independent.Gaussians.Incidence.Only=r3.u.indall,
                      Compound.Symmetry.Incidence.Only=r3.u.cs5,
                      Independent.Gaussians.With.Mortality=r3.u.indall.m,
                      Compound.Symmetry.With.Mortality=r3.u.cs5.m) + TEXT_THEME_SMALLER


plot.posteriors(Independent.Gaussians.Incidence.Only=r3.u.indall,
                Compound.Symmetry.Incidence.Only=r3.u.cs5,
                Independent.Gaussians.With.Mortality=r3.u.indall.m,
                Compound.Symmetry.With.Mortality=r3.u.cs5.m,
                color=c('red','blue','orange','black'))

plot.multiple.results(Independent.Gaussians.Incidence.Only=r3.u.indall,
                      Compound.Symmetry.Incidence.Only=r3.u.cs5,
                      Independent.Gaussians.With.Mortality=r3.u.indall.m,
                      Compound.Symmetry.With.Mortality=r3.u.cs5.m,
                      incidence = F) + TEXT_THEME_SMALLER

plot.dist.comparison(Compound.Symmetry.Incidence.Only=r3.u.cs5,
                     Compound.Symmetry.With.Mortality=r3.u.cs5.m,
                     Independent.Gaussians.Incidence.Only=r3.u.indall,
                     Independent.Gaussians.With.Mortality=r3.u.indall.m,
                     dist.name='inc.change.metrics', ncol=2,
                     vars.to.plot=SOME.OUTCOMES.TO.PLOT1) + TEXT_THEME_SMALLER




#different compound symmetry coefficients
r3.u.cs7 = resample.and.get.output('sampling/model3', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.7)
r3.u.cs3 = resample.and.get.output('sampling/model3', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.3)
r3.u.cs1 = resample.and.get.output('sampling/model3', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.1)

plot.dist.comparison(rho.0.0.Indep=r3.u.indall,
                     rho.0.1=r3.u.cs1,
                     rho.0.3=r3.u.cs3,
                     rho.0.5=r3.u.cs5,
                     rho.0.7=r3.u.cs7,
                     rho.0.9=r3.u.cs9,
                     dist.name='inc.change.metrics', ncol=3,
                     vars.to.plot=OUTCOMES.TO.PLOT) + ggtitle('Comparison of Compound Symmetry coefficients')

plot.dist.comparison(rho.0.0.Indep=r3.u.indall,
                     rho.0.1=r3.u.cs1,
                     rho.0.3=r3.u.cs3,
                     rho.0.5=r3.u.cs5,
                     rho.0.7=r3.u.cs7,
                     rho.0.9=r3.u.cs9,
                     dist.name='inc.change.metrics', ncol=2,
vars.to.plot=SOME.OUTCOMES.TO.PLOT1) + TEXT_THEME_SMALLER + theme(legend.position='bottom')


plot.multiple.results(rho.0.0.Indep=r3.u.indall,
                       rho.0.1=r3.u.cs1,
                       rho.0.3=r3.u.cs3,
                       rho.0.5=r3.u.cs5,
                       rho.0.7=r3.u.cs7,
                       rho.0.9=r3.u.cs9) + TEXT_THEME 


#model structure comparisons

r.model1 = resample.and.get.output('sampling/model1', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
r.model2 = resample.and.get.output('sampling/model2', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
#r.model3 = resample.and.get.output('sampling/model3', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
r.model3 = r3.u.cs5

r1.u.indall = resample.and.get.output('sampling/model1', independent.gaussian.incidence.likelihood, years=2000:2015)
r2.u.indall = resample.and.get.output('sampling/model2', independent.gaussian.incidence.likelihood, years=2000:2015)


OUTCOMES.TO.PLOT=NULL
OUTCOMES.TO.PLOT=ONE.OUTCOME.TO.PLOT
plot.dist.comparison(Model_1=r.model1,
                     Model_2=r.model2,
                     Model_3=r.model3,
                     dist.name='inc.change.metrics', ncol=3,
                     vars.to.plot=OUTCOMES.TO.PLOT) + TEXT_THEME + theme(legend.position='bottom')

plot.dist.comparison(Model_1=r.model1,
                     Model_2=r.model2,
                     Model_3=r.model3,
                     dist.name='inc.change.metrics', ncol=2,
                     vars.to.plot=SOME.OUTCOMES.TO.PLOT1) + TEXT_THEME_SMALLER + theme(legend.position='bottom')


plot.multiple.results(Model_1=r.model1,
                      Model_2=r.model2,
                      Model_3=r.model3) + TEXT_THEME  + theme(legend.position='bottom')

plot.dist.comparison(Model_1=r1.u.indall,
                     Model_2=r2.u.indall,
                     Model_3=r3.u.indall,
                     dist.name='inc.change.metrics', ncol=3,
                     vars.to.plot=OUTCOMES.TO.PLOT) + ggtitle('Comparison of Compound Symmetry coefficients')
plot.multiple.results(Model_1=r1.u.indall,
                      Model_2=r2.u.indall,
                      Model_3=r3.u.indall) + TEXT_THEME  + theme(legend.position='bottom')


#The single model run graphic
df = get.plot.df(r3.u.indall)
df = df[df$id=='model.pre' | df$id=='model.noint',]
df$ci.lower = df$estimate - 8 - rnorm(31,2,.5)
df$ci.upper = df$estimate + 8 + rnorm(31,2,.5)


ggplot(df, aes(year, estimate, color=id, fill=id, shape=id)) +
    geom_point(size=4) +
    geom_line(size=1) +
    geom_ribbon(aes(ymin=ci.lower, ymax=ci.upper), alpha=0.3) +
    xlab('Year') + ylab('TB Incidence per 100,000 (95% CI)') +
    ylim(0,310) +
    TEXT_THEME


df = get.plot.df(r3.u.indall)
df = df[df$id=='real',]
    
ggplot(df, aes(year, estimate)) +
    geom_point(size=4) +
    geom_line(size=1) +
    xlab('Year') + ylab('TB Incidence per 100,000') +
    ylim(0,310) +
    TEXT_THEME


df = get.plot.df(r3.u.cs5)

df = df[(df$id=='real' | df$id=='model.pre') & df$year<=2015,]

df2 = get.plot.df(r3.u.uall)
df2 = df2[df2$id=='model.pre' & df2$year<=2015,]
df2$id='model2'
df2$estimate = df2$estimate + 15 - (df2$year-2000)/15 * (15+15)

df = rbind(df,df2)

df$ci.lower[df$id=='real'] = INDIA.DATA$Incidence.CI.Lower
df$ci.upper[df$id=='real'] = INDIA.DATA$Incidence.CI.Upper
df$ci.lower[df$id!='real'] = NA
df$ci.upper[df$id!='real'] = NA


df = get.plot.df(r3.u.cs5, incidence = F)
df$ci.lower[df$id=='real'] = INDIA.DATA$Mortality.CI.Lower
df$ci.upper[df$id=='real'] = INDIA.DATA$Mortality.CI.Upper
df$ci.lower[df$id!='real'] = NA
df$ci.upper[df$id!='real'] = NA

df = df[df$id=='real',]

ggplot(df, aes(year, estimate, color=id, fill=id, shape=id)) +
    geom_ribbon(aes(ymin=ci.lower, ymax=ci.upper), alpha=0.3) +
    geom_point(size=4) +
    geom_line(size=1) +
    scale_color_manual(name='Estimate Source:',
                       labels=c(model.pre='Model Run 1', model2='Model Run 2', real='WHO Estimate'),
                       values=c(model.pre=MODEL_COLOR, model2='blue', real=REAL_COLOR)) +
    scale_fill_manual(name='Estimate Source:',
                      labels=c(model.pre='Model Run 1', model2='Model Run 2', real='WHO Estimate'),
                      values=c(model.pre=MODEL_COLOR, model2='blue', real=REAL_COLOR)) +
    scale_shape_manual(name='Estimate Source:',
                       labels=c(model.pre='Model Run 1', model2='Model Run 2', real='WHO Estimate'),
                       values=c(model.pre=15, model2=17, real=19)) +
    xlab('Year') + ylab('TB Incidence per 100,000') +
    ylim(0,NA) +
    TEXT_THEME

sd.2000 = INDIA.DATA$Incidence.SD[1]
real.2000=df$estimate[df$id=='real'&df$year==2000]
lower.2000=INDIA.DATA$Incidence.CI.Lower[1]
upper.2000=INDIA.DATA$Incidence.CI.Upper[1]
run1.2000=df$estimate[df$id=='model.pre'&df$year==2000]
run2.2000=df$estimate[df$id=='model2'&df$year==2000]

plot.dist.with.intercept <- function(param1, param2, dist, intercept, 
                                     color.dist, color.intercept,
                                     n.x=1000)
{
    alpha=0.0001
    
    d.fn = get(paste0('d',dist))
    
    if (dist=='unif')
    {
        x.min = param1 - (param2-param1)/20
        x.max = param2 + (param2-param1)/20
    }
    else
    {
        q.fn = get(paste0('q',dist))
        x.min = q.fn(alpha, param1, param2)
        x.max = q.fn(1-alpha, param1, param2)
    }

    x = seq(from=x.min, to=x.max, length.out=n.x)
    
    df = data.frame(x=x, d=d.fn(x, param1, param2))
print(d.fn(intercept, param1, param2))    
    ggplot(df, aes(x, d)) +
        geom_ribbon(aes(ymin=0, ymax=d), alpha=0.3, fill=color.dist) +
        geom_line(size=1, color=color.dist) +
        geom_vline(xintercept = intercept, color=color.intercept, size=2) +
        geom_point(aes(x=intercept, y=d.fn(intercept, param1, param2)), color=color.intercept, size=6)
}

x.limits = xlim(-50,650)
y.limits = ylim(0,0.0052)
plot.dist.with.intercept(run1.2000, sd.2000, 'norm', real.2000, MODEL_COLOR, REAL_COLOR) + x.limits + y.limits + TEXT_THEME_SMALLER
plot.dist.with.intercept(run2.2000, sd.2000, 'norm', real.2000, 'blue', REAL_COLOR) + x.limits + y.limits + TEXT_THEME_SMALLER

plot.dist.with.intercept(run1.2000 + (lower.2000-real.2000), run1.2000 + (upper.2000-real.2000), 'unif', real.2000, MODEL_COLOR, REAL_COLOR) + x.limits + y.limits + TEXT_THEME_SMALLER
plot.dist.with.intercept(run2.2000 + (lower.2000-real.2000), run2.2000 + (upper.2000-real.2000), 'unif', real.2000, 'blue', REAL_COLOR) + x.limits + y.limits + TEXT_THEME_SMALLER


plot.multiple.results(Independent.Gaussians=r3.u.indall) + TEXT_THEME
plot.with.indiv.resamples(r3.u.indall)

plot.multiple.results(Uniforms=r3.u.uall) + TEXT_THEME

#PRCCs
prccs= calculate.pccs(r3.u.cs5$resamples, get.main.outcome)
plot.prccs(prccs, parameter.names=parameter.names) + TEXT_THEME

#quantiles
quantile.df=get.inference.on.resampled.quantiles.matrix(r3.u.cs5$resamples, get.main.outcome, quantile.lower.bounds=c(low=0,high=0.8), quantile.upper.bounds=c(0.2,1), return.as.data.frame = T, statistic='median')
boxplot.quantile.pairs(quantile.df, parameter.names=parameter.names) + TEXT_THEME
