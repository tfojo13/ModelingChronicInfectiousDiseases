##-------------##
##-- IMPORTS --##
##-------------##

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

MAIN.TEXT.SIZE = 26
AXIS.LABEL.SIZE = 24

#FOUR.COLOR.SCHEME = c('#FF2A00', '#FF9900', '#0FAD00', '#0064B5')
FOUR.COLOR.SCHEME = c('#379F9E', '#FFF87E', '#E1F7A6', '#FFD78D')

TEXT_THEME = theme(text = element_text(size=MAIN.TEXT.SIZE), axis.text.x = element_text(size=AXIS.LABEL.SIZE))

##-------------------------##
##--RESAMPLING FUNCTIONS --##
##-------------------------##

resample.and.output.for.paper <- function(directory, likelihood, groups=NULL, ...)
{
    res = resample.distributed.samples(directory, NUM.RESAMPLES, likelihood, groups = groups...)
    
    rv = list()
    rv$resamples = res
    
    
    rv
}

##-------------------------------------##
##-- DO RESAMPLING OR LOAD RESAMPLES --##
##-------------------------------------##

if (1==2) {
    r1.u.cs5 = resample.and.get.output('sampling/model1', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
    r2.u.cs5 = resample.and.get.output('sampling/model2', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
    r3.u.cs5 = resample.and.get.output('sampling/model3', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5, groups=1:2)
    
    r3.unw.cs5 = resample.and.get.output('sampling/model3_unarrow', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
    r3.log.cs5 = resample.and.get.output('sampling/model3_log', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
    r3.lognw.cs5 = resample.and.get.output('sampling/model3_lognarrow', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
    
    r3.u.ind2 = resample.and.get.output('sampling/model3', independent.gaussian.incidence.likelihood, years=c(2000,2015), groups=1:2)
    r3.u.indall = resample.and.get.output('sampling/model3', independent.gaussian.incidence.likelihood, years=2000:2015, groups=1:2)
    r3.u.uall = resample.and.get.output('sampling/model3', uniform.incidence.likelihood, years=2000:2015, groups=1:2)
    
    r3.u.cs5.m = resample.and.get.output('sampling/model3', mvg.cs.mort.inc.likelihood, years=2000:2015, rho=0.5)
    r3.u.indall.m = resample.and.get.output('sampling/model3', independent.gaussian.mort.inc.likelihood, years=2000:2015)
    
    save(r1.u.cs5, r2.u.cs5, r3.u.cs5,
         r3.unw.cs5, r3.log.cs5, r3.lognw.cs5,
         r3.u.ind2, r3.u.indall, r3.u.uall,
         r3.u.cs5.m, r3.u.indall.m,
         file='saved/saved_for_paper_plots.Rdata')
} else {
    load('saved/saved_for_paper_plots.Rdata')
}

## For comparison of models
if (1==2)
{
    #export 720x540
    plot.results(r1.u.cs5, ylim.upper=420, labels=T) + TEXT_THEME + theme(legend.position="none")
    plot.results(r2.u.cs5, ylim.upper=420, labels=T) + TEXT_THEME + theme(legend.position="none")
    plot.results(r3.u.cs5, ylim.upper=420, labels=T) + TEXT_THEME + theme(legend.position="none")
    
}

## For comparison of priors
if (1==2)
{
    #export all 720x540
#    prior.model.names = c('Uniform Priors\n(Reference Model)', 'Log- and\nLogit-Normal\nPriors', 'Narrow Log-\n and Logit\n-Normal Priors')
    prior.model.names = c('Uniform Priors\n(Reference Model)', 'Shaped\nPriors', 'Narrow\nShaped Priors')
    
    boxplot.dist.comparison(r3.u.cs5, r3.log.cs5, r3.lognw.cs5, vars.to.plot='rel.inc.diff.2030.2030', model.names=prior.model.names, colors=FOUR.COLOR.SCHEME[1:3]) + TEXT_THEME + theme(legend.position='none', axis.title.y=element_blank())

    plot.results(r3.u.cs5, ylim.upper=420) + TEXT_THEME + theme(legend.position="none")
    plot.results(r3.log.cs5, ylim.upper=420) + TEXT_THEME + theme(legend.position="none")
    plot.results(r3.lognw.cs5, ylim.upper=420) + TEXT_THEME + theme(legend.position="none")
    

    round(100*t(sapply(list(r3.u.cs5, r3.log.cs5, r3.lognw.cs5), function(r){
        ci = get.credible.intervals(r$paper.metrics, 'rel.inc.diff.2030.2030', smoothed=F)
        c(mean=as.numeric(get.means(r$paper.metrics, 'rel.inc.diff.2030.2030', smoothed=F)),
          ci.lower = ci[1], ci.upper=ci[2], ci.width=ci[2]-ci[1])
    })))
}

## for comparison of likelihoods
if (1==2)
{
    #export 720 x 540
    plot.results(r3.u.cs5, ylim.upper=475) + TEXT_THEME + theme(legend.position="none")
    plot.results(r3.u.indall, ylim.upper=475) + TEXT_THEME + theme(legend.position="none")
    plot.results(r3.u.ind2, ylim.upper=475) + TEXT_THEME + theme(legend.position="none")
    plot.results(r3.u.uall, ylim.upper=475) + TEXT_THEME + theme(legend.position="none")
    
#    likelihood.model.names = c('Compound\nSymmetry\n(Reference\nModel)', 'Independent\nNormals\n16 points', 'Independent\nNormals\n2 points', 'Uniform\n16 points')
    likelihood.model.names = c('Independent\nNormals,\n2 Points', 'Uniforms', 'Independent\nNormals,\n16 Points', 'Correlated\nNormals\n(Reference\nModel)')
        
    boxplot.dist.comparison(r3.u.ind2, r3.u.uall, r3.u.indall, r3.u.cs5, vars.to.plot='rel.inc.diff.2030.2030', model.names=likelihood.model.names, colors=FOUR.COLOR.SCHEME[c(2,3,4,1)]) + TEXT_THEME + theme(axis.text.x = element_text(size=21)) + theme(legend.position='none', axis.title.y=element_blank()) + ylab('Change in 2030 Incidence (%)')
#    boxplot.dist.comparison(r3.u.cs5, r3.u.indall, r3.u.ind2, r3.u.uall, vars.to.plot='inc.2030.int.2030', model.names=likelihood.model.names, box.and.error=F, percents=F) + TEXT_THEME + theme(legend.position='none', axis.title.y=element_blank()) + theme(axis.title.y=element_text(angle = 90)) + ylab('Estimated 2030 Incidence\n(per 100,000)')
#    boxplot.dist.comparison(r3.u.cs5, r3.u.indall, r3.u.ind2, r3.u.uall, vars.to.plot='inc.2030.noint.2030', model.names=likelihood.model.names, box.and.error=F, percents=F) + TEXT_THEME + theme(legend.position='none', axis.title.y=element_blank()) + theme(axis.title.y=element_text(angle = 90)) + ylab('Estimated 2030 Incidence\n(per 100,000)')
    
    
    plot.results(r3.u.cs5, ylim.upper=420) + TEXT_THEME + theme(legend.position="none")
    plot.results(r3.u.cs5.m, ylim.upper=420) + TEXT_THEME + theme(legend.position="none")
    
    plot.results(r3.u.cs5, ylim.upper=150, incidence=F) + TEXT_THEME + theme(legend.position="none")
    plot.results(r3.u.cs5.m, ylim.upper=150, incidence=F) + TEXT_THEME + theme(legend.position="none")
    
    likelihood.model.names.2 = c('Calibrated to   \nIncidence   \n(Reference   \nModel)   ', '   Calibrated\n   to Incidence\n   and Mortality')
    
    #export 470x540
    boxplot.dist.comparison(r3.u.cs5, r3.u.cs5.m, vars.to.plot='rel.inc.diff.2030.2030', model.names=likelihood.model.names.2, colors=FOUR.COLOR.SCHEME[1:2]) + TEXT_THEME + theme(legend.position='none', axis.title.y=element_blank(), axis.text.x = element_text(size=AXIS.LABEL.SIZE-1))
    boxplot.dist.comparison(r3.u.cs5, r3.u.cs5.m, vars.to.plot='rel.mort.diff.2030.2030', model.names=likelihood.model.names.2, colors=FOUR.COLOR.SCHEME[1:2]) + TEXT_THEME + theme(legend.position='none', axis.title.y=element_blank(), axis.text.x = element_text(size=AXIS.LABEL.SIZE-1))
    
}

## For sensitivity analyses
if (1==2)
{
    #PRCCs - 1050 x 475
    prccs= calculate.pccs(r3.u.cs5$resamples, get.main.outcome)
    plot.prccs(prccs, parameter.names=parameter.names, max.to.plot=6) + TEXT_THEME + theme(axis.title.y=element_blank()) + ylab('Partial Rank Correlation Coefficient')
    
    #quantiles - 1200 x 475
    quantile.df=get.inference.on.resampled.quantiles.matrix(r3.u.cs5$resamples, get.main.outcome, quantile.lower.bounds=c(low=0,high=0.8), quantile.upper.bounds=c(0.2,1), return.as.data.frame = T, statistic='median')
    boxplot.quantile.pairs(quantile.df, parameter.names=parameter.names, max.to.plot = 6) + TEXT_THEME + theme(axis.title.y=element_blank(), legend.title=element_text(size=24), legend.text=element_text(size=21))
    
}

#for distribution graphic
if (1==2)
{
   who = get.india.data(2000:2015, just.estimates = T) 
   sample = get.incidence.output(r3.u.indall$resamples@results[[300]], NULL, 2000:2015, F)
   
   df = data.frame(Year=rep(2000:2015,2),
                   estimate=c(who, sample),
                   Source=rep(c('WHO Estimate', 'Simulation'), each=16))
   
   colors = c(REAL_COLOR, MODEL_COLOR)
   shapes = c(15, 19)
   names(colors) = names(shapes) = unique(df$Source)
   
   FADING=0.4
   df$alpha = sapply(df$Year, function(y){
       if (y==2000 || y==2015)
           'full'
       else
           'faded'
   })
   
   #export 500x680
   ggplot(df, aes(Year, estimate, color=Source, shape=Source)) + 
       geom_point(size=6, aes(alpha=alpha)) + 
       geom_line(size=1, alpha=FADING) + ylim(0, NA) + 
       scale_color_manual(name='', values=colors) + scale_shape_manual(name='', values=shapes) + ylim(0, 350) + ylab('TB Incidence per 100,000') +
       theme(text=element_text(size=34), legend.position = 'bottom') + xlim(NA, 2015.5) + geom_hline(yintercept=0) +
       guides(alpha=F) + scale_alpha_manual(values=c(full=1,faded=FADING))
   

   plot.likelihood.density <- function(d.fn, param1, param2, who, x.min=0, x.max=600)
   {
       x = seq(x.min, x.max, length=1000)
       d = d.fn(x, param1, param2)

       ggplot() + geom_line(aes(x=x, y=d), size=3, color=MODEL_COLOR) + geom_area(aes(x=x, y=d), fill=MODEL_COLOR, alpha=0.4) + 
           geom_vline(xintercept=who, color=REAL_COLOR, size=4) + 
           geom_point(aes(x=who, y=d.fn(who, param1, param2)), color=REAL_COLOR, size=15) +
           theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title=element_blank())

   }
   
   #export 300x400
   plot.likelihood.density(dnorm, sample[16], INDIA.DATA$Incidence.SD[16], who[16]) + ylim(0, 0.007)
   plot.likelihood.density(dnorm, sample[1], INDIA.DATA$Incidence.SD[1], who[1]) + ylim(0, 0.007)
   
   
   unif.lower = sample + who - INDIA.DATA$Incidence.CI.Upper
   unif.upper = sample + who - INDIA.DATA$Incidence.CI.Lower
   
   plot.likelihood.density(dunif, unif.lower[16], unif.upper[16], who[16]) + ylim(0, 0.007)
   plot.likelihood.density(dunif, unif.lower[1], unif.upper[1], who[1]) + ylim(0, 0.007)
}