#'@title Class representing a multi-dimensional empiric distribution
#'
#'@description An S4 class that represents the empiric distribution of multivariate data (ie, a distribution defined by a set of values with associated weights)
#'
#'@slot num.vars The number of variables (dimensions) in the distribution
#'@slot var.names The names of the variables (dimensions) in this distribution
#'
#'@details The class can be accessed with a number of methods:\itemize{
#'\item For descriptive statistics: \code{\link{get.means}}, \code{\link{get.medians}}, \code{\link{get.modes}}, and \code{\link{get.quantiles}}
#'\item For inference/uncertainty ranges: \code{\link{get.credible.intervals}}, \code{\link{get.variances}}, and \code{\link{get.sds}}
#'\item For plotting: \code{\link{plot.distributions}} and \code{\link{get.data.frame.for.plot}}
#'}
#'
#'@seealso The constructor function: \code{\link{Empiric.Distribution}}\cr
#'\code{\link{Joint.Canonical.Distribution}}
#'@export
setClass("Empiric.Distribution",
         representation(
             num.vars = 'numeric',
             var.names = 'character',
             samples = 'matrix',
             weights = 'numeric',
             smoothed.marginal.pmfs = 'list',
             raw.marginal.pmfs = 'list'
         ))

#'@title Constructor for an instance of an \code{\link{Empiric.Distribution-class}}
#'
#'@param samples A matrix, where each row corresponds to one set realization values of the variables of this distribution (ie, each column represents a variable, and each row represents one set of values for those variables)
#'@param weights The weights given to each row of samples
#'
#'@return an instance of \code{\link{Empiric.Distribution-class}}
#'@export
setGeneric('Empiric.Distribution', def=function(samples, weights=rep(1, dim(samples)[1]), allow.smoothing=T, bandwidth.adjust=1)
{
    if (class(samples) != 'matrix')
        samples = matrix(samples, ncol=1)
    
    
    num.vars = dim(samples)[2]
#    keep.rows = rep(T, num.vars)
    
    #remove duplicate rows
#    for (r1 in 1:(num.vars-1))
#    {
#        if (keep.rows[r1])
#        {
#            for (r2 in (r1+1):num.vars)
#            {
#                if (all(samples[r1,]==samples[r2,]))
#                {
#                    keep.rows[r2] = F
#                    weights[r1] = weights[r1] + weights[r2]
#                }
#            }
#        }
#    }
    
#    samples = samples[keep.rows,]
#    weights = weights[keep.rows]
    
    num.vars = dim(samples)[2]
    weights = weights/sum(weights)
    
    smoothed.pmfs = list()
    raw.pmfs = list()
    for (i in 1:num.vars)
    {
        raw.pmfs[[i]] = create.marginal.pmf(samples[,i], weights, FALSE)
        if (allow.smoothing)
            smoothed.pmfs[[i]] = create.marginal.pmf(samples[,i], weights, TRUE, bandwidth.adjust)
        else
            smoothed.pmfs[[i]] = raw.pmfs[[i]]
    }
    var.names = dimnames(samples)[[2]]
    if (is.null(var.names))
        var.names = as.character(1:num.vars)
    
    names(smoothed.pmfs) = names(raw.pmfs) = var.names
    names(var.names) = var.names
    
    new('Empiric.Distribution', num.vars=num.vars,
                             var.names=var.names,
                             samples=samples,
                             weights=weights,
                             smoothed.marginal.pmfs=smoothed.pmfs,
                             raw.marginal.pmfs=raw.pmfs)
})


create.marginal.pmf <- function(values, density, smoothed=F, bandwidth.adjust)
{
    if (smoothed)
    {
        dobj = density(values, weights = density/sum(density), adjust=bandwidth.adjust)
        values = dobj$x
        density = dobj$y/sum(dobj$y)
    }
    else
    {
#        o = order(values)
#        values = values[o]
#        density = density[o]
        
#        exclude.indices = integer()
#        i=1
#        while (i < length(values))
#        {
#            j = i+1
#            while (j <= length(values) && values[j]==values[i])
#            {
#                exclude.indices = c(exclude.indices, j)
#                density[i] = density[i] + density[j]
#                j=j+1
#            }
#            i=j
#        }
        
#        values = values[-(exclude.indices)]
#        density = density[-(exclude.indices)]

        new.values = sort(unique(values))
        new.density = numeric(length(new.values))
        for (i in 1:length(new.values))
            new.density[i] = sum(density[values==new.values[i]])
        
        new.density = new.density / sum(new.density)
        
        values = new.values
        density = new.density
    }
    
    list(x=values, d=density)
}

setMethod('show', 'Empiric.Distribution', def=function(object)
{
    variables_text = if(object@num.vars > 1) 'variables:' else 'variable:'
    cat('An Empiric.Distribution with', object@num.vars, variables_text, paste(paste0("'",object@var.names,"'"), collapse=', '))
})

#'@title Marginalize variables out of an Empiric.Distribution
#'
#'@description Marginalizes out some of the variables of an empiric distribution, returning a distribution over the remaining subset of variables
#'
#'@param dist An instance of an \code{\link{Empiric.Distribution-class}}
#'@param vars.to.keep The indices or names of the variables which are NOT to be marginalized
#'
#'@return An instance of an \code{\link{Empiric.Distribution-class}}
#'
#'@export
setGeneric('marginalize', def=function(dist, vars.to.keep){standardGeneric('marginalize')})
setMethod('marginalize', 'Empiric.Distribution', def=function(dist, vars.to.keep)
{
    Empiric.Distribution(dist@samples[,vars.to.keep], dist@weights)
})

#'@title Extract a statistic (mean, median, or mode) plus a credible interval
#'
#'@description A convenience method that combines a call to either \code{\link{get.means}}, \code{\link{get.medians}}, or \code{\link{get.modes}} with a call to \code{\link{get.credible.intervals}}
#'
#'@param dist An instance of an \code{\link{Empiric.Distribution-class}}
#'@param vars The indices or names of the variables to get information for
#'@param statistic The summary statistic to be returned. Can be either 'mean', 'median', or 'mode'
#'@param coverage The probability (0 to 1) which the interval should cover. The returned interval will cover at least this coverage, but may cover more if an exact interval is not possible
#'@param interval.type The type of interval to be returned. Options are 'highest-density' or 'equal-tailed'
#'@param smoothed A logical indicating whether a smoothed distribution over the variables of interest should be used (as opposed to using the exact samples given in constructing the distribution)
#'
#'@return A matrix where each row corresponds to a variable in the vars argument, with three columns: the first column representing either the mean, median, or mode as specified, followed by columns 'ci.lower' and 'ci.upper' denoting the bounds of the interval
#'
#'@export
setGeneric('get.statistic.and.intervals', def=function(dist, vars=dist@var.names, statistic='mean', coverage=0.95, interval.type='highest-density', smoothed=T){standardGeneric('get.statistic.and.intervals')})
setMethod('get.statistic.and.intervals', 'Empiric.Distribution', def=function(dist, vars=dist@var.names, statistic='mean', coverage=0.95, interval.type='highest-density', smoothed=T)
{
    if (statistic=='mean')
        stats = get.means(dist, vars, smoothed)
    else if (statistic=='median')
        stats = get.medians(dist, vars, smoothed)
    else if (statistic=='mode')
        stats = get.modes(dist, vars, smoothed)
    else
        stop("statistic must be either 'mean', 'median', or 'mode")
    
    intervals = matrix(get.credible.intervals(dist, vars, coverage, interval.type, smoothed)[,1:2], ncol=2)
   
    rv = cbind(stats, intervals)
    dimnames(rv)[[2]] = c(statistic, 'ci.lower', 'ci.upper')
    
    rv
})


#'@title Extract a Credible Interval from an Empiric Distribution
#'
#'@description Returns a (marginal) credible interval for each of a specified subset of the variables covered by this Empiric Distribution
#'
#'@param dist An instance of an \code{\link{Empiric.Distribution-class}}
#'@param vars The indices or names of the variables to get intervals for
#'@param coverage The probability (0 to 1) which this interval should cover. The returned interval will cover at least this coverage, but may cover more if an exact interval is not possible
#'@param type The type of interval to be returned. Options are 'highest-density' or 'equal-tailed'
#'@param smoothed A logical indicating whether a smoothed distribution over the variables of interest should be used (as opposed to using the exact samples given in constructing the distribution)
#'
#'@return A matrix with a row for each variable containing 'ci.lower', 'ci.upper', and 'actual.coverage'
#'
#'@export
setGeneric('get.credible.intervals', def=function(dist, vars=dist@var.names, coverage=0.95, type='highest-density', smoothed=T){standardGeneric('get.credible.intervals')})
setMethod('get.credible.intervals', 'Empiric.Distribution', def=function(dist, vars=dist@var.names, coverage=0.95, type='highest-density', smoothed=T)
{
    rv = matrix(0, nrow=length(vars), ncol=3)
    dimnames(rv)[[1]] = as.list(dist@var.names[vars])
    dimnames(rv)[[2]] = c("ci.lower", "ci.upper", "actual.coverage")
    
    for (var in vars)
    {
        pmf = get.pmf(dist, var, smoothed)
        if (type=='highest-density')
            rv[var, ] = calc.highest.density.ci(pmf, coverage)
        else if (type=='equal-tailed')
            rv[var, ] = calc.equal.tailed.ci(pmf, coverage)
        else
            stop("the 'type' argument must be either 'highest-density' or 'equal-tailed'")
    }
    
    attr(rv, 'type') = type
    rv
})

get.pmf <- function(dist, var, smoothed)
{
    if (smoothed)
        dist@smoothed.marginal.pmfs[[var]]
    else
        dist@raw.marginal.pmfs[[var]]
}

calc.highest.density.ci <- function(pmf, coverage)
{
    old.lower = new.lower = 1
    old.upper = new.upper = length(pmf$x)
    
    while(sum(pmf$d[new.lower:new.upper]) >= coverage && (new.lower < length(pmf$x) || new.upper > 1))
    {
        old.lower = new.lower
        old.upper = new.upper
        
        if (pmf$d[old.upper] < pmf$d[old.lower] && old.upper > 1)
            new.upper=old.upper-1
        else
            new.lower=old.lower+1
    }
    
    rv=c(pmf$x[old.lower], pmf$x[old.upper], sum(pmf$d[old.lower:old.upper]))
    names(rv)=c('ci.lower', 'ci.upper', 'actual.coverage')
    rv
}

calc.equal.tailed.ci <- function(pmf, coverage)
{
    n = length(pmf$x)
    old.lower = new.lower = 1
    old.upper = new.upper = n
    
    while(sum(pmf$d[new.lower:new.upper]) >= coverage)
    {
        old.lower = new.lower
        old.upper = new.upper
        
        if (sum(pmf$d[old.upper:n]) < sum(pmf$d[1:old.lower]))
            new.upper=old.upper-1
        else
            new.lower=old.lower+1
    }
    
    rv=c(pmf$x[old.lower], pmf$x[old.upper], sum(pmf$d[old.lower:old.upper]))
    names(rv)=c('ci.lower', 'ci.upper', 'actual.coverage')
    rv
}

#'@title Extract Means from an Empiric Distribution
#'
#'@description Returns (marginal) means for a specified subset of the variables covered by this Empiric Distribution
#'
#'@param dist An instance of an \code{\link{Empiric.Distribution-class}}
#'@param vars The indices or names of the variables to get means for
#'@param smoothed A logical indicating whether a smoothed distribution over the variables of interest should be used (as opposed to using the exact samples given in constructing the distribution)
#'
#'@return A vector of means of the specified variables
#'
#'@export
setGeneric('get.means', def=function(dist, vars=dist@var.names, smoothed=T){standardGeneric('get.means')})
setMethod('get.means', 'Empiric.Distribution', def=function(dist, vars=dist@var.names, smoothed=T)
{
    rv = numeric(length(vars))
    names(rv) = dist@var.names[vars]
    
    for (var in vars)
        rv[var] = calc.mean(get.pmf(dist, var, smoothed))
    
    rv
})

calc.mean <- function(pmf)
{
    sum(pmf$x * pmf$d) / sum(pmf$d)
}

#'@title Extract Medians from an Empiric Distribution
#'
#'@description Returns (marginal) medians for a specified subset of the variables covered by this Empiric Distribution
#'
#'@param dist An instance of an \code{\link{Empiric.Distribution-class}}
#'@param vars The indices or names of the variables to get medians for
#'@param smoothed A logical indicating whether a smoothed distribution over the variables of interest should be used (as opposed to using the exact samples given in constructing the distribution)
#'
#'@return A vector of medians for the specified variables
#'
#'@export
setGeneric('get.medians', def=function(dist, vars=dist@var.names, smoothed=T){standardGeneric('get.medians')})
setMethod('get.medians', 'Empiric.Distribution', def=function(dist, vars=dist@var.names, smoothed=T)
{
    rv = numeric(length(vars))
    names(rv) = dist@var.names[vars]
    
    for (var in vars)
        rv[var] = calc.median(get.pmf(dist, var, smoothed))
    
    rv
})

calc.median <- function(pmf)
{
    cum.prob = cumsum(pmf$d)
    inv.cum.prob = rev(cumsum(rev(pmf$d)))
    
    mean(pmf$x[cum.prob >= 0.5 & inv.cum.prob >= 0.5])
}

#'@title Extract Modes from an Empiric Distribution
#'
#'@description Returns (marginal) modes for a specified subset of the variables covered by this Empiric Distribution
#'
#'@param dist An instance of an \code{\link{Empiric.Distribution-class}}
#'@param vars The indices or names of the variables to get modes for
#'@param smoothed A logical indicating whether a smoothed distribution over the variables of interest should be used (as opposed to using the exact samples given in constructing the distribution)
#'
#'@return A vector of modes for the specified variables
#'
#'@export
setGeneric('get.modes', def=function(dist, vars=dist@var.names, smoothed=T){standardGeneric('get.modes')})
setMethod('get.modes', 'Empiric.Distribution', def=function(dist, vars=dist@var.names, smoothed=T)
{
    rv = numeric(length(vars))
    names(rv) = dist@var.names[vars]
    
    for (var in vars)
        rv[var] = calc.mode(get.pmf(dist, var, smoothed))
    
    rv
})

calc.mode <- function(pmf)
{
    max.d = max(pmf$d)
    modes = pmf$x[pmf$d == max.d]
    
    modes[ceiling(length(modes)/2)]
}


#'@title Extract Quantiles from an Empiric Distribution or Joint Distributution
#'
#'@description Returns (marginal) quantiles for a specified subset of the variables covered by this Empiric Distribution
#'
#'@param dist An instance of an \code{\link{Empiric.Distribution-class}} or a \code{\link{Joint.Canonical.Distribution-class}}
#'@param vars The indices or names of the variables to get quantiles for
#'@param probs Numeric vector of probabilites with values [0,1]
#'@inheritParams get.means
#'
#'@return \itemize{
#'\item If probs is only a single value, returns a vector where each element represents the quantile for a variable.\cr
#'\item If vars refers to only a single value, returns a vector where each element is a quantile of that variable.\cr
#'\item If both probs and vars have length greater than 1, returns a matrix where each row corresponds to a variable in the distribution and each columns corresponds to a quantile defined by the argument in probs (ie, rv[var,i] represents the probs[i]th quantile for variable var)
#'}
#'
#'@export
setGeneric('get.quantiles', def=function(dist, vars=dist@var.names, smoothed=T, probs=seq(0, 1, 0.25)){standardGeneric('get.quantiles')})
setMethod('get.quantiles', 'Empiric.Distribution', def=function(dist, vars=dist@var.names, smoothed=T, probs=seq(0, 1, 0.25))
{
    rv = matrix(0, nrow=length(vars), ncol=length(probs), dimnames=list(dist@var.names[vars], paste0(probs*100,'%')))

    for (var in vars)
        rv[var,] = calc.quantiles(get.pmf(dist, var, smoothed), probs)
    
    if (length(vars) == 1)
        rv[1,]
    else if (length(probs==1))
        rv[,1]
    else
        rv
})

calc.quantiles <- function(pmf, probs)
{
    rv = NULL
    cum.prob = cumsum(pmf$d)
    for (prob in probs)
    {
        q = pmf$x[cum.prob>=prob][1]
        rv = c(rv, q)
    }
    
    rv
}

#'@title Evaluate the Cumulative Distribution Function of an Empiric Distribution or Joint Distribution
#'
#'@param dist An instance of an \code{\link{Empiric.Distribution-class}} or a \code{\link{Joint.Canonical.Distribution-class}}
#'@param var The index or name of the variable to evaluate CDF for
#'@param values Numeric vector of values at which to evaluate the CDF
#'@inheritParams get.means
#'
#'@return A vector of probabilities corresponding to the CDF evaluated at each element of values
#'
#'@export
setGeneric('evaluate.cdf', def=function(dist, values, var, smoothed=T){standardGeneric('evaluate.cdf')})
setMethod('evaluate.cdf', 'Empiric.Distribution', def=function(dist, values, var, smoothed=T)
{
    pmf = get.pmf(dist, var, smoothed)
    sapply(values, calc.cdf, pmf)
})

calc.cdf <- function(value, pmf)
{
    sum(pmf$d[pmf$x<=value])
}

#'@title Extract Variances from an Empiric Distribution
#'
#'@description Returns (marginal) variances for a specified subset of the variables covered by this Empiric Distribution
#'
#'@param dist An instance of an \code{\link{Empiric.Distribution-class}}
#'@param vars The indices or names of the variables to get variances for. Covariances are not given
#'@param smoothed A logical indicating whether a smoothed distribution over the variables of interest should be used (as opposed to using the exact samples given in constructing the distribution)
#'
#'@return A vector of variances for the specified variables
#'
#'@export

setGeneric('get.variances', def=function(dist, vars=dist@var.names, smoothed=T){standardGeneric('get.variances')})
setMethod('get.variances', 'Empiric.Distribution', def=function(dist, vars=dist@var.names, smoothed=T)
{
    rv = numeric(length(vars))
    names(rv) = dist@var.names[vars]
    
    for (var in vars)
        rv[var] = calc.var(get.pmf(dist, var, smoothed))
    
    rv
})

calc.var <- function(pmf)
{
    ex2 = sum(pmf$x^2 * pmf$d) / sum(pmf$d)
    ex2 - calc.mean(pmf)^2
}


#'@title Extract Standard Deviations from an Empiric Distribution
#'
#'@description Returns (marginal) standard deviations for a specified subset of the variables covered by this Empiric Distribution
#'
#'@param dist An instance of an \code{\link{Empiric.Distribution-class}}
#'@param vars The indices or names of the variables to get standard deviations for
#'@param smoothed A logical indicating whether a smoothed distribution over the variables of interest should be used (as opposed to using the exact samples given in constructing the distribution)
#'
#'@return A vector of standard deviations for the specified variables
#'
#'@export
setGeneric('get.sds', def=function(dist, vars=dist@var.names, smoothed=T){standardGeneric('get.sds')})
setMethod('get.sds', 'Empiric.Distribution', def=function(dist, vars=dist@var.names, smoothed=T)
{
    sqrt(var(dist, vars))
})



#'@title Get a data frame to be used by plotting functions
#'
#'@description Generates a data frame that can be passed to plotting functions in order to plot the density. This function is useful to be able to extend plots; if all you want is to plot the densities without adding anything, see \code{\link{plot.distributions}}
#'
#'@param ... One or more instances of \code{\link{Empiric.Distribution-class}} or \code{\link{Joint.Canonical.Distribution-class}}
#'@param vars The names of the variables to plot distributions for. The variables indicated by this parameter must be contained in every distribution passed to ... If passed NULL, chooses the variables that are present in all given distributions
#'@param smoothed A logical indicating whether a smoothed distribution over the variables of interest should be used. Can be either a scalar value that is applied to all distributions, or a logical vector whose ith element indicates whether the ith distribution passed to ... should be smoothed
#'@param names A vector of names, with one name for each distribution passed to ... If there is only one distribution, this argument is ignored.
#'@param lower.quantile.bounds, upper.quantile.bounds Numerics indicating at what quantiles of Joint.Canonical.Distributions the plots should be truncated. Can be either a scalar value that is applied to all distributions, or a numeric vector whose ith element indicates the bound for the ith distribution passed to ...
#'@return A data frame with three or four columns:\itemize{
#'\item Value (containing the x-values at which densities are specified)
#'\item Density (the density for each value)
#'\item Variable (the variable name for each row in the data frame)
#'\item Distribution (only included if there is more than one distribution; this field contains the name as given in the names argument)
#'}
#'
#'@examples
#'dist #A previously calculated empiric distribution object
#'df = get.plot.data.frame(dist)
#'
#'#Now generate a plot of the density, with one facet for each variable
#'ggplot(df, aes(Value, Density)) + geom_line() + facet_wrap(~Variable)
#'
#'@seealso \code{\link{plot.distributions}}
#'
#'@export
get.data.frame.for.plot <- function(..., vars=NULL, smoothed=T,
                                    names=c('Posterior Distribution', 'Prior Distribution'),
                                    lower.quantile.bounds = 0.001, upper.quantile.bounds = 1-lower.quantile.bounds)
{
    df = data.frame()
    dists = flatten.list(list(...))
    if (length(dists) > 1 && length(names) != length(dists))
        stop(paste0('There must be one name for each distribution. ', length(dists), ' distributions and ', length(names), ' names were provided.'))
    
    smoothed = resize.vector(smoothed, length(dists))
    lower.quantile.bounds = resize.vector(lower.quantile.bounds, length(dists))
    upper.quantile.bounds = resize.vector(upper.quantile.bounds, length(dists))
    
    if (is.null(vars) || length(vars) < 1)
    {
        vars = dists[[1]]@var.names
        for (dist in dists)
            vars = intersect(vars, dist@var.names)
        
        if (length(vars)==0)
            stop('The specified distributions do not have any variables in common')
    }
    else
    {
        for (dist in dists)
        {
            if (length(intersect(dist@var.names, vars)) < length(vars))
                stop('The specified variables are NOT contained in all the given distributions')
        }
    }
    

    i = 0 
    for (dist in dists)
    {
        i = i+1
        one.df = get.one.plot.df(dist, vars=vars, smoothed=smoothed[i], lower.quantile.bounds[i], upper.quantile.bounds[i])
        if (length(dists)>1)
            one.df$Distribution = rep(names[i], dim(one.df)[1])
        
        df = rbind(df, one.df)
    }
        
    if (length(dists) > 1)
        df$Distribution = factor(df$Distribution, levels=names)

    df
}

SCALE.TO.POINTS = 512

get.one.plot.df <- function(dist, vars, smoothed, lower.quantile.bound, upper.quantile.bound, n=512)
{
    df = data.frame()

    for (var in vars)
    {
        if (class(dist)=='Empiric.Distribution')
            pmf = get.pmf(dist, var, smoothed)
        else if (class(dist)=='Joint.Canonical.Distribution')
            pmf = get.plot.points(dist, var, lower.quantile.bound, upper.quantile.bound)
        else
            stop("The distributions given must be either class 'Empiric.Distribution' or 'Joint.Canonical.Distribution'")
        
        if (class(dist) == 'Empiric.Distribution' && !smoothed)
        {
            collapse.to = max(16,min(length(pmf$d)/16, SCALE.TO.POINTS/5))
            pmf = collapse.pmf(pmf, collapse.to)
        }
            
        d = length(pmf$d) * pmf$d / sum(pmf$d) / SCALE.TO.POINTS
            
        df = rbind(df,
                   data.frame(Value=pmf$x, Density=d, Variable=rep(dist@var.names[var], length(pmf$x))))
    }

    df
}

collapse.pmf <- function(pmf, len)
{
    len = round(len)
    bounds = seq(min(pmf$x), max(pmf$x), length=len+1)
    new.x = (bounds[-1]+bounds[1:len])/2
    bounds[1] = bounds[1] - (bounds[2]-bounds[1])/1000

    new.d = sapply(1:len, function(i){
        sum(pmf$d[pmf$x <= bounds[i+1] & pmf$x > bounds[i]])
    })
    
    list(x=new.x, d=new.d/sum(new.d))
}

get.plot.df.one.var <- function(dists, var, smoothed, n.points, names, quantile.bounds)
{
    df=data.frame()
    
    min = Inf
    max = -Inf
    for (i in 1:length(dists))
    {
        dist = dists[[i]]
        bounds = get.quantiles(dist, var, smoothed[i], probs=quantile.bounds)
        
        if (bounds[1] < min)
            min = bounds[1]
        if (bounds[2] > max)
            max = bounds[2]
    }
    
    points = seq(min, max, length=n.points)
    separation = points[2] - points[1]
    interval.bounds = c(min-separation, points) + separation/2
    
    for (i in 1:length(dists))
    {
        dist = dists[[i]]
        cdf.values = evaluate.cdf(dist, interval.bounds, var, smoothed[i])
        d = cdf.values[-1] - cdf.values[-(n.points+1)]
        one.df = data.frame(Value=points, Density=d, Variable=rep(dist@var.names[var], n.points))
        if (length(dists) > 1)
            one.df$Distribution = rep(names[i], n.points)
        
        df = rbind(df, one.df)
    }
    
    df
}


#'@title Plot Densities from two Empiric Distributions
#'
#'@description Makes plots of the densities from this Empiric Distribution against other 'comparison' densities defined by the given arguments
#'
#'@param color,shape,linetype,fill Aesthetics which can vary based on distribution - all passed to the ggplot call. Each of these should be a vector of values (or scalar if there is only one distribution) where the ith value is the aesthetic for the ith distribution in ...
#'@param ncol,nrow The final plot will facet wrap the different variables in the multivariate distribution. Use these to specify the number of rows and columns in the wrap, or leave null to let ggplot decide
#'@param geom.to.use The ggplot geom_ to use in plotting (eg geom_point, geom_line, etc)
#'@param hide.legend Logical indicating whether to hide the legend in the plot (by default, with one distribution the legend is hidden)
#'@inheritParams get.data.frame.for.plot
#'
#'@return a ggplot2 object
#'
#'@seealso \code{\link{get.data.frame.for.plot}}
#'
#'@export
setGeneric('plot.distributions', def=function(..., vars=NULL, smoothed=T,
                                   names=c('Posterior Distribution', 'Prior Distribution'),
                                   lower.quantile.bounds=0.001, upper.quantile.bounds=1-lower.quantile.bounds,
                                   geom.to.use=geom_line, 
                                   color = c('red', 'blue', 'green', 'gray', 'orange', 'purple'), 
                                   fill=color,  
                                   shape = c(19, 17, 15, 18, 4, 3), 
                                   linetype = 'solid',
                                   ncol=NULL, nrow=NULL,
                                   hide.legend=F)
{
    df = get.data.frame.for.plot(..., vars=vars, smoothed=smoothed, names=names, lower.quantile.bounds=lower.quantile.bounds, upper.quantile.bounds=upper.quantile.bounds)
    
#    o = order(names)
#    color = color[o]
#    shape = shape[o]
#    linetype = linetype[o]
#    fill = fill[o]

    num.dists = length(flatten.list(list(...)))
    single.dist = is.null(df$Distribution)
    if (single.dist)
        df$Distribution='1'

    color = resize.vector(color, num.dists)
    fill = resize.vector(fill, num.dists)
    shape = resize.vector(shape, num.dists)
    linetype = resize.vector(linetype, num.dists)
    
         
    plot = ggplot(df, aes(Value, Density, color=Distribution, shape=Distribution, linetype=Distribution, fill=Distribution)) + 
        geom.to.use() + 
        scale_color_manual(values=color) + 
        scale_shape_manual(values=shape) + 
        scale_linetype_manual(values=linetype) + 
        scale_fill_manual(values=fill) + 
        facet_wrap(~Variable, scales='free', ncol=ncol, nrow=nrow) +
        ylim(0,NA)
    
    if (single.dist || hide.legend)
        plot + guides(color=F, fill=F, shape=F, linetype=F)
    else
        plot
})

#returns a vector of size n
# if n = length(v), returns v
# if n < length(v), returns v[1:n]
# if n > length(v), wraps v as many times as necessary to reach vector of size n
resize.vector <- function(v, n)
{
    if (n <= length(v))
        v[1:n]
    else
        rep(v, ceiling(n/length(v)))[1:n]
}


flatten.list <- function(l)
{
    rv = list()
    for (elem in l)
    {
        if (class(elem)=='list')
            rv = c(rv, flatten.list(elem))
        else
            rv = c(rv, elem)
    }
    rv
}
