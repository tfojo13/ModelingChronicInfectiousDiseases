source('data.R')
source('models.R')
source('inference.R')
library(bayesian.modeling)
library(gridExtra)

MODEL_COLOR = "orange"
REAL_COLOR = "#458B00"
MODEL_COLOR2 = '#3278FA'

get.plot.df <- function(dists, incidence=T, smoothed=F)

{
    index = if (incidence) 1 else 2

    real.data = get.india.data(PRE.YEARS, incidence, just.estimates=T)
    
    pre.estimates = get.medians(dists$pre.dists[[index]], smoothed=smoothed)
    pre.cis = get.credible.intervals(dists$pre.dists[[index]], smoothed=smoothed)
    
    noint.estimates = get.medians(dists$noint.dists[[index]], smoothed=smoothed)
    noint.cis = get.credible.intervals(dists$noint.dists[[index]], smoothed=smoothed)
    
    int.estimates = get.medians(dists$int.dists[[index]], smoothed=smoothed)
    int.cis = get.credible.intervals(dists$int.dists[[index]], smoothed=smoothed)

    df.pre = data.frame(year=PRE.YEARS,
                        estimate = pre.estimates,
                        ci.lower = pre.cis[,1],
                        ci.upper = pre.cis[,2],
                        source = 'model',
                        intervention = 'pre',
                        id = 'model.pre')
    
    df.noint = data.frame(year=POST.YEARS,
                          estimate = noint.estimates,
                          ci.lower = noint.cis[,1],
                          ci.upper = noint.cis[,2],
                          source = 'model',
                          intervention = 'noint',
                          id = 'model.noint')

    df.int = data.frame(year=POST.YEARS,
                        estimate = int.estimates,
                        ci.lower = int.cis[,1],
                        ci.upper = int.cis[,2],
                        source = 'model',
                        intervention = 'int',
                        id = 'model.int')
    
    df.real = data.frame(year=2000:2015,
                         estimate = real.data,
                         ci.lower = NA,
                         ci.upper = NA,
                         source = 'WHO',
                         intervention = 'pre',
                         id='real')

    df = rbind(df.pre, df.noint, df.int, df.real)

    df
}
 
plot.results <- function(dists, incidence=T, smoothed=F, labels=F, ylim.upper=NA, point.size=4, line.size=1, label.size=8.5)
{
    df = get.plot.df(dists, incidence, smoothed)
    
    index = if (incidence) 1 else 2
    pre.estimates = get.medians(dists$pre.dists[[index]], smoothed=smoothed)
    pre.cis = get.credible.intervals(dists$pre.dists[[index]], smoothed=smoothed)
    
    noint.estimates = get.medians(dists$noint.dists[[index]], smoothed=smoothed)
    noint.cis = get.credible.intervals(dists$noint.dists[[index]], smoothed=smoothed)
    
    int.estimates = get.medians(dists$int.dists[[index]], smoothed=smoothed)
    int.cis = get.credible.intervals(dists$int.dists[[index]], smoothed=smoothed)
    
    #for label
    max.y = max(df$ci.upper, na.rm = T)
    max.x = max(df$year)
    change.dist = dists$change.dists[[index]]
    change.means = get.means(change.dist, smoothed=smoothed)
    change.cis = get.credible.intervals(change.dist, smoothed=smoothed)
    noint.label = paste0('No Intervention: ', round(100*change.means['noint.change.2030'],1), '% [',round(100*change.cis['noint.change.2030',1],1),'% to ',round(100*change.cis['noint.change.2030',2],1),'%]')
    int.label = paste0('With Intervention: ', round(100*change.means['int.change.2030'],1), '% [',round(100*change.cis['int.change.2030',1],1),'% to ',round(100*change.cis['int.change.2030',2],1),'%]')
 #   delta.label = paste0('Inc. Reduction From Intervention: ', round(100*change.means['rel.reduction.2030'],1), '% [',round(100*change.cis['rel.reduction.2030',1],1),'% to ',round(100*change.cis['rel.reduction.2030',2],1),'%]')
    delta.label = paste0('Change in 2030 Incidence:\n',-round(100*change.means['rel.reduction.2030'],0), '% [',-round(100*change.cis['rel.reduction.2030',2],0),'% to ',-round(100*change.cis['rel.reduction.2030',1],0),'%]')
 
    end = length(noint.estimates)
    inc.noint.label = paste0('Inc. 2030 w/o Intervention: ', round(noint.estimates[end],0), ' [', round(noint.cis[end,1],0),' to ', round(noint.cis[end,2]),']')
    inc.int.label = paste0('Inc. 2030 w/ Intervention: ', round(int.estimates[end],0), ' [', round(int.cis[end,1],0),' to ', round(int.cis[end,2]),']')
    
    rv = ggplot(df, aes(year, estimate, color=id, fill=id, shape=id)) +
        geom_point(size=point.size) +
        geom_line(size=line.size) +
        geom_ribbon(aes(ymin=ci.lower, ymax=ci.upper), alpha=0.3) +
        scale_color_manual(name='Estimate Source:',
                           labels=c(model.pre='Model Estimate Before Intervention', model.noint='Model Estimate Without Intervention', model.int='Model Estimate With Intervention', real='WHO Estimate'),
                           values=c(model.pre=MODEL_COLOR, model.noint=MODEL_COLOR, model.int=MODEL_COLOR2, real=REAL_COLOR)) +
        scale_fill_manual(name='Estimate Source:',
                          labels=c(model.pre='Model Estimate Before Intervention', model.noint='Model Estimate Without Intervention', model.int='Model Estimate With Intervention', real='WHO Estimate'),
                          values=c(model.pre=MODEL_COLOR, model.noint=MODEL_COLOR, model.int=MODEL_COLOR2, real=REAL_COLOR)) +
        scale_shape_manual(name='Estimate Source:',
                           labels=c(model.pre='Model Estimate Before Intervention', model.noint='Model Estimate Without Intervention', model.int='Model Estimate With Intervention', real='WHO Estimate'),
                           values=c(model.pre=19, model.noint=19, model.int=19, real=15)) +
        xlab('Year') + ylab(paste0(if (incidence) 'TB Incidence' else 'TB Mortality', ' per 100,000')) +
        ylim(0,ylim.upper) + geom_hline(yintercept=0)
    
    if (labels)
    {
        label.df = data.frame(x=max.x, y=max.y, label=delta.label, id='real')
        
        rv = rv +
        geom_label(aes(x=max.x, y=max.y, vjust='center', hjust='center', label=delta.label), 
                   fill='white', color='black', 
                   size=label.size, label.padding = unit(0.5, "lines"),
                   nudge_x=-9.5, nudge_y=-10)
    }
    rv
}

plot.multiple.results <- function(..., incidence=T, smoothed=F)
{
    models = list(...)
    model.names = names(models)
    
    df = NULL
    for (model.name in model.names)
    {
        dists = models[[model.name]]
        one.df = get.plot.df(dists, incidence, smoothed)
        one.df$model = model.name
        
        df = rbind(df, one.df)
    }
    ggplot(df, aes(year, estimate, color=id, fill=id, shape=id)) +
        geom_point() +
        geom_line() +
        geom_ribbon(aes(ymin=ci.lower, ymax=ci.upper), alpha=0.3) +
        scale_color_manual(name='Estimate Source:',
                           labels=c(model.pre='Model Estimate Before Intervention', model.noint='Model Estimate Without Intervention', model.int='Model Estimate With Intervention', real='WHO Estimate'),
                           values=c(model.pre=MODEL_COLOR, model.noint=MODEL_COLOR, model.int=MODEL_COLOR2, real=REAL_COLOR)) +
        scale_fill_manual(name='Estimate Source:',
                          labels=c(model.pre='Model Estimate Before Intervention', model.noint='Model Estimate Without Intervention', model.int='Model Estimate With Intervention', real='WHO Estimate'),
                          values=c(model.pre=MODEL_COLOR, model.noint=MODEL_COLOR, model.int=MODEL_COLOR2, real=REAL_COLOR)) +
        scale_shape_manual(name='Estimate Source:',
                           labels=c(model.pre='Model Estimate Before Intervention', model.noint='Model Estimate Without Intervention', model.int='Model Estimate With Intervention', real='WHO Estimate'),
                           values=c(model.pre=19, model.noint=15, model.int=17, real=19)) +
        xlab('Year') + ylab(paste0(if (incidence) 'TB Incidence' else 'TB Mortality', ' per 100,000 (95% CI)')) +
        ylim(0,NA) +
        facet_wrap(~model)
}

plot.with.indiv.resamples <- function(dists, incidence=T, smoothed=F, res=1:5)
{
    BIGGEST.LINE=2
    
    index = if (incidence) 1 else 2
    
    years = 2000:2015
    
    real.data = get.india.data(years, incidence, just.estimates=T)
    
    pre.estimates = get.means(dists$pre.dists[[index]], smoothed=smoothed)[1:length(years)]
    pre.cis = get.credible.intervals(dists$pre.dists[[index]], smoothed=smoothed)[1:length(years),]
    
    df.pre = data.frame(year=years,
                        estimate = pre.estimates,
                        ci.lower = pre.cis[,1],
                        ci.upper = pre.cis[,2],
                        source = 'Posterior.Distribution',
                        type = 'mean',
                        weight=BIGGEST.LINE)
    
    df.real = data.frame(year=years,
                         estimate = real.data,
                         ci.lower = NA,
                         ci.upper = NA,
                         source = 'WHO.Estimates',
                         type = 'real',
                         weight=BIGGEST.LINE)
    
    df = rbind(df.pre, df.real)
    
    max.lik = max(dists$resamples@likelihoods)
    for (i in res)
    {
        rr = dists$resamples@results[[i]]
        lik = dists$resamples@likelihoods[i]
        
        df = rbind(df,
                   data.frame(year=years,
                              estimate=get.model.rates(rr, years, incidence=incidence),
                              ci.lower=NA,
                              ci.upper=NA,
                              source=paste0('Resample ', i, ' (', round(lik/max.lik,2), ')'),
                              type='resample',
                              weight=BIGGEST.LINE * lik / max.lik / 2))
    }
    
    min.y = min(df$ci.lower, na.rm = T)
    min.x = min(df$year)
    linreg.who = do.get.linreg(years, real.data)
    linreg.dist = dists$linreg.dists[[index]]
    linreg.means = get.means(linreg.dist, smoothed=smoothed)
    linreg.cis = get.credible.intervals(linreg.dist, smoothed=smoothed)
    int.label = paste0('Model Intercept: ', round(linreg.means['intercept'],0), ' [',round(linreg.cis['intercept',1],0),' to ',round(linreg.cis['intercept',2],0),']')
    slope.label = paste0('Model Slope: ', round(linreg.means['slope'],1), ' [',round(linreg.cis['slope',1],1),' to ',round(linreg.cis['slope',2],1),']')

    
    ggplot(df, aes(year, estimate, color=source, fill=source, shape=type)) +
        geom_point(aes(size=weight)) +
        geom_line() +
        geom_ribbon(aes(ymin=ci.lower, ymax=ci.upper), alpha=0.3) +
        scale_size_continuous(guide=F) +
        scale_shape_manual(name='estimate type',
                           labels=c(mean='Mean of Posterior Dist', real='WHO Estimate', resample='Individual Model Run'),
                           values=c(mean=19, real=15, resample=17)) +
        xlab('Year') + ylab(paste0(if (incidence) 'TB Incidence' else 'TB Mortality', ' per 100,000 (95% CI)')) +
        geom_label(aes(x=min.x, y=min.y, vjust='bottom', hjust='left',
                       label=paste0("'True' Intercept: ", round(linreg.who['intercept'],0), '\n',
                                    int.label, '\n',
                                    "'True' Slope: ", round(linreg.who['slope'],1), '\n',
                                    slope.label
                       )), fill='white', color='black')
}

plot.with.indiv.resamples.simple <- function(dists, incidence=T, smoothed=F, res=1:5)
{
    BIGGEST.LINE=2
    
    index = if (incidence) 1 else 2
    
    years = 2000:2015
    
    real.data = get.india.data(years, incidence, just.estimates=T)
    
    pre.estimates = get.means(dists$pre.dists[[index]], smoothed=smoothed)[1:length(years)]
    pre.cis = get.credible.intervals(dists$pre.dists[[index]], smoothed=smoothed)[1:length(years),]

    
    df = data.frame(year=years,
                         estimate = real.data,
                         source = 'WHO.Estimates',
                         type = 'real',
                         weight=BIGGEST.LINE)
    
    max.lik = max(dists$resamples@likelihoods)
    for (i in res)
    {
        rr = dists$resamples@results[[i]]
        lik = dists$resamples@likelihoods[i]
        
        df = rbind(df,
                   data.frame(year=years,
                              estimate=get.model.rates(rr, years, incidence=incidence),
                              source=paste0('Resample ', i, ' (', round(lik/max.lik,2), ')'),
                              type='resample',
                              weight=BIGGEST.LINE * lik / max.lik ))
    }


    
    ggplot(df, aes(year, estimate, color=source, fill=source, shape=type)) +
        geom_point(aes(size=weight)) +
        geom_line() +
        scale_size_continuous(guide=F) +
        scale_shape_manual(name='estimate type',
                           labels=c(mean='Mean of Posterior Dist', real='WHO Estimate', resample='Individual Model Run'),
                           values=c(mean=19, real=15, resample=17)) +
        xlab('Year') + ylab(paste0(if (incidence) 'TB Incidence' else 'TB Mortality', ' per 100,000 (95% CI)'))
}

get.dist.comparison.df <- function(..., dist.name, smoothed=F, vars.to.plot=NULL, model.names=NULL)
{
    models = list(...)
    if (is.null(model.names))
        model.names = names(models)
    
    df = NULL
    
    if (is.null(vars.to.plot))
        vars.to.plot = models[[1]][[dist.name]]@var.names
    
    for (i in 1:length(models))
    {
        dist = models[[i]][[dist.name]]
        
        means = get.means(dist, smoothed=smoothed, vars=vars.to.plot)
        medians = get.medians(dist, smoothed=smoothed, vars=vars.to.plot)
        cis = get.credible.intervals(dist, smoothed=smoothed, vars=vars.to.plot)
        iqs = get.credible.intervals(dist, smoothed=smoothed, vars=vars.to.plot, coverage=0.5, type='equal-tailed')
        
        df = rbind(df, data.frame(
                   estimate=medians,
                   ci.lower = cis[,1],
                   ci.upper = cis[,2],
                   iq.lower = iqs[,1],
                   iq.upper = iqs[,2],
                   metric = vars.to.plot,
                   model = model.names[i]))
    }
    
    df$metric = factor(df$metric, levels=vars.to.plot)
    
    df
}

plot.dist.comparison <- function(..., dist.name, smoothed=F, ncol=NA, vars.to.plot=NULL)
{
    df = get.dist.comparison.df(..., dist.name=dist.name, smoothed=smoothed, vars.to.plot=vars.to.plot)


    ggplot(df, aes(x=model, y=estimate, color=model, fill=model)) + 
        geom_bar(stat='identity', position='dodge', alpha=0.5) + 
        geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), size=1) +
        geom_hline(yintercept=0, size=1) + 
        facet_wrap(~metric, scales='free', ncol=ncol) +
        theme(axis.text.x=element_blank()) +
        ylab('Estimate of Metric (95% Credible Interval)') + xlab('Model') +
        scale_fill_discrete(name='Model') + scale_color_discrete('Model')
        
    
}

boxplot.dist.comparison <- function(..., dist.name='paper.metrics', 
                                 smoothed=F, ncol=NA, vars.to.plot=NULL, 
                                 model.names=NULL, box.and.error=F, percents=T,
                                 line.size=2, dodge.width=1, colors=NULL)
{
    df = get.dist.comparison.df(..., dist.name=dist.name, smoothed=smoothed, vars.to.plot=vars.to.plot, model.names=model.names)
   
    rv = ggplot(df, aes(x=model, y=estimate, lower=iq.lower, middle=estimate, upper=iq.upper, ymin=ci.lower, ymax=ci.upper, fill=model)) +
        theme(legend.position='bottom', axis.title.x=element_blank()) +
        geom_hline(yintercept=0, size=1)
   
    if (percents)
        rv = rv + scale_y_continuous(labels=function(y){paste0(round(100*y,0),"%")})
    
    if (box.and.error)
        rv = rv + geom_errorbar(width=0.15, size=1) + geom_point(size=15, shape=22)
    else
        rv = rv + geom_errorbar(width=0.15, size=2, position=position_dodge(width=dodge.width)) + 
                geom_boxplot(aes(width=0.6), position=position_dodge(width=dodge.width), size=1, stat='identity') +
        
    if (length(unique(df$metric)) > 1)
        rv = rv + facet_wrap(~metric, scales='free', ncol=ncol)
    
    if (!is.null(colors))
    {
#        names(colors) = unique(df$model)
        rv = rv + scale_fill_manual(name='Model', values=colors) + scale_color_manual('Model', values=colors)
    }
    else
    {
        rv = rv + scale_fill_discrete(name='Model') + scale_color_discrete('Model')
    }
    
    rv
}



plot.posteriors <- function(..., smoothed=T, ncol=NA, param.names=NULL, color=NULL)
{
    models = list(...)
    model.names = names(models)
    param.dists = lapply(models, function(rr){
        rr$param.dist
    })
    if (is.null(param.names))
        param.names=param.dists[[1]]@var.names
    
    plot.distributions(param.dists, smoothed=smoothed, names=model.names, ncol=ncol, vars=param.names, color=color)
}


plot.prccs <- function(prccs, color='#2F4F4F', parameter.names=NULL, max.to.plot = NA)
{
    o=order(abs(prccs$PRCC))
    names = dimnames(prccs$PRCC)[[1]][o]
    values = prccs$PRCC[o,]
    df = data.frame(Parameter=names, PRCC=values)
    
    df$Parameter = as.character(df$Parameter)
    if (!is.null(parameter.names))
    {
        df$Parameter = parameter.names[df$Parameter]
    }
 
    num.to.plot = min(max.to.plot, length(names))
    params.to.keep = parameter.names[names[length(names) + 1 - 1:num.to.plot]]
    df = df[sapply(df$Parameter, function(pp){any(pp==params.to.keep)}),]
    
    df$Parameter = factor(df$Parameter, levels=factor(df$Parameter))
    
    ggplot(data=df, aes(Parameter, PRCC, fill=Scenario)) + 
        geom_bar(stat='identity', fill=color, width=0.8, alpha=0.85) + coord_flip() +
        ylim(-1,1) + 
        geom_hline(yintercept=1, linetype='solid') + geom_hline(yintercept=-1, linetype='solid') + geom_hline(yintercept=0, linetype='solid') +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
}

boxplot.quantile.pairs <- function(quantile.df, color1=MODEL_COLOR, color2=MODEL_COLOR2, parameter.names=NULL, max.to.plot=NA, percents=T, dodge.width=0.6)
{
    quantile.df$Parameter = as.character(quantile.df$Parameter)
  #  quantile.df$Parameter.Quantile = as.character(quantile.df$Parameter.Quantile)
    
    diffs = quantile.df[quantile.df$Parameter.Quantile=='high', 'median'] - quantile.df[quantile.df$Parameter.Quantile=='low', 'median']
    o = order(abs(diffs), decreasing = F)
    
 #   quantile.df = rbind(quantile.df[quantile.df$Parameter.Quantile=='high',][o,],
  #                      quantile.df[quantile.df$Parameter.Quantile=='low',][o,])
  
    if (!is.null(parameter.names))
    {
        quantile.df$Parameter = parameter.names[quantile.df$Parameter]
    }
    param.names = unique(quantile.df$Parameter)[o]
    
    num.to.plot = min(length(param.names), max.to.plot, na.rm=T)
    param.names = param.names[length(param.names) - num.to.plot + 1:num.to.plot]
    quantile.df = quantile.df[sapply(quantile.df$Parameter, function(pp){any(param.names==pp)}),]
    
    quantile.df$Parameter = factor(quantile.df$Parameter, levels=param.names)
#    quantile.df$Parameter.Quantile = factor(quantile.df$Parameter.Quantile, levels=c('high','low'))
    
#    quantile.df$Parameter.Quantile = factor(quantile.df$Parameter.Quantile, levels=c('high','low'))
    
    rv = ggplot(quantile.df, aes(Parameter, median, color=Parameter.Quantile)) +
        geom_point(position=position_dodge(width=dodge.width), shape=15, size=9) +
        geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), position=position_dodge(width=dodge.width), size=2) + 
        coord_flip() +
        ylab('Estimate of Primary Outcome') +
        scale_color_manual(name='Quantile of\nParameter', guide = guide_legend(reverse=T),
                           labels=c(low='Bottom 20%', high='Top 20%'),
                           values=c(high=color1, low=color2)) + 
        theme(panel.background = element_blank())

    if (percents)
        rv = rv + scale_y_continuous(labels=function(y){paste0(round(100*y,0),"%")})
    
    rv
}