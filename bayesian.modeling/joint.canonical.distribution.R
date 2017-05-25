
#'@title The Joint.Canonical.Distribution class
#'
#'@description A class representing a set of independent distributions (eg Normal, Uniform, Beta - any distribution for which R has functions like qnorm, rnorm, etc)
#'
#'@details For plotting, use \code{\link{plot.distributions}} and \code{\link{get.data.frame.for.plot}}
#'@seealso The constructor \code{\link{Joint.Canonical.Distribution}}\cr
#'\code{\link{Empiric.Distribution}}
#'
#'@export
setClass('Joint.Canonical.Distribution',
         representation(
             dist.names = 'character',
             quantile.functions = 'list',
             pdfs = 'list',
             cdfs = 'list',
             random.functions = 'list',
             parameters = 'list',
             num.vars = 'integer',
             var.names = 'character',
             use.var.names = 'logical'
         ))

#'@title The constructor for a Joint.Canonical.Distribution object
#'
#'@param ... A set of parameters that defines distributions. Each distribution should be defined first by a 'name' that corresponds to the suffix used in R function (eg 'norm' for qnorm, dnorm, 'unif' for qunif, dunif, etc), followed by any number of parameters that would be passed to those functions. For example, 'norm', 0, 1 represents a normal distribution with mean 0 and sd 1.\cr
#'These parameters can be passed either separated by commas, or all in one list, or in multiple lists (see example below)
#'@param var.names A character vector of names for each variable in the distribution
#'
#'@examples
#'#Say we want to construct a distribution consisting of a normal(2,1), a uniform(-1,1), 
#'# and a beta(2,2) with names 'var1', 'var2', and 'var3. 
#'# The following are three equivalent ways of doing so:
#'
#'Joint.Canonical.Distribution('norm',2,1,'unif',-1,1,'beta',2,2,var.names=c('var1','var2','var3'))
#'
#'distributions = list('norm',2,1,'unif',-1,1,'beta',2,2)
#'Joint.Canonical.Distribution(distributions, var.names=c('var1','var2','var3'))
#'
#'dist1 = list('norm', 2, 1)
#'dist2 = list('unif', -1, 1)
#'dist3 = list('beta', 2, 2)
#'Joint.Canonical.Distribution(dist1, dist2, dist3, var.names=c('var1','var2','var3'))
#'
#'@export
setGeneric('Joint.Canonical.Distribution', def=function(..., var.names=NULL)
{
    dist.names = character()
    parameters = list()
    
    i = 0
    for (arg in flatten.list(list(...)))
    {
        if (class(arg)=='character')
        {
            i = i+1
            dist.names[i] = arg
            parameters[[i]] = numeric()
        }
        else
            parameters[[i]] = c(parameters[[i]], arg)
    }
    
    num.vars = as.integer(length(dist.names))
    
    quantile.functions = list()
    pdfs = list()
    cdfs = list()
    random.functions = list()
    for (i in 1:num.vars)
    {
        dist.name = dist.names[i]
        if (!exists(paste0('d',dist.name)))
            stop(paste0("ERROR: '", dist.name, "' is not a valid distribution name in R"))
        
        params = parameters[[i]]
        
        quantile.functions[[i]] = create.dist.function(dist.name, 'q', params)
        pdfs[[i]] = create.dist.function(dist.name, 'd', params)
        cdfs[[i]] = create.dist.function(dist.name, 'p', params)
        random.functions[[i]] = create.dist.function(dist.name, 'r', params)
    }
    
    use.var.names = !is.null(var.names)
    if (is.null(var.names))
        var.names = as.character(1:num.vars)
    else if (num.vars != length(var.names))
        stop(paste0('var.names has length ', length(var.names), ' but ', num.vars, ' distributions were specified'))
    
    names(var.names) = names(quantile.functions) = names(pdfs) = names(cdfs) = names(random.functions) = var.names
        
    new('Joint.Canonical.Distribution',
        dist.names = dist.names,
        parameters = parameters,
        var.names = var.names,
        num.vars = num.vars,
        quantile.functions = quantile.functions,
        pdfs = pdfs,
        cdfs = cdfs,
        random.functions = random.functions,
        use.var.names = use.var.names)
})

setMethod('show', 'Joint.Canonical.Distribution', def=function(object)
{    
    variables_text = if(object@num.vars > 1) 'variables:' else 'variable:'
    cat('A Joint.Canonical.Distribution with', object@num.vars, variables_text)
    for (i in 1:object@num.vars)
        cat('\n',i,") '", object@var.names[i], "' ~ ", object@dist.names[i], '(', paste0(object@parameters[[i]],collapse=','), ')', sep='')

})

#'@title Evaluate the joint density for this Joint Distribution
#'
#'@param dist An instance of \code{\link{Joint.Canonical.Distribution-class}}
#'@param values Numeric vector of values at which to evaluate the joint density. If the values vector has names, the values are mapped to variables by name; otherwise, the 1st value is taken to correspond to the 1st variable, and so forth
#'
#'@return The joint density (ie the product of the density for each variable in the Joint Distribution evaluated at the value in values)
#'
#'@export
setGeneric('get.density', def=function(dist, values){standardGeneric('get.density')})
setMethod('get.density', 'Joint.Canonical.Distribution', def=function(dist, values)
{
    values = map.values(values, dist)
    prod(sapply(1:dist@num.vars, function(i){
        dist@pdfs[[i]](values[i])
    }))
})

#'@title Get a function that evaluates a joint density for all variables in the Joint Distribution
#'
#'@param dist An instance of \code{\link{Joint.Canonical.Distribution-class}}
#'
#'@return A function that takes as its first argument a numeric vector of values and returns the joint density (ie the product of the density for each variable in the Joint Distribution evaluated at the value in values)
#'
#'@export
setGeneric('density.function', def=function(dist){standardGeneric('density.function')})
setMethod('density.function', 'Joint.Canonical.Distribution', def=function(dist)
{
    function(values, ...)
    {
        get.density(dist, values)
    }
})

#'@title Get a function that evaluates quantiles for all variables in the Joint Distribution
#'
#'@param dist An instance of \code{\link{Joint.Canonical.Distribution-class}}
#'
#'@return A function that takes as its first argument a numeric vector of [0,1] values and returns the corresponding quantile for each variable in the Joint Distribution.
#'
#'@export
setGeneric('quantiles.function', def=function(dist){standardGeneric('quantiles.function')})
setMethod('quantiles.function', def=function(dist)
{
    function(values, ...)
    {
        values = map.values(values, dist)
        rv=sapply(1:dist@num.vars, function(i){
            dist@quantile.functions[[i]](values[i])
        })
        names(rv) = dist@var.names
        rv
    }
})

setMethod('get.quantiles', 'Joint.Canonical.Distribution', def=function(dist, vars=dist@var.names, smoothed=T, probs=seq(0, 1, 0.25))
{
    rv = matrix(0, nrow=length(vars), ncol=length(probs), dimnames=list(dist@var.names[vars], paste0(probs*100,'%')))
    
    for (var in vars)
        rv[var,] = sapply(probs, dist@quantile.functions[[var]])
    
    if (length(vars) == 1)
        rv[1,]
    else if (length(probs==1))
        rv[,1]
    else
        rv
})

setMethod('evaluate.cdf', 'Joint.Canonical.Distribution', def=function(dist, values, var, smoothed=T)
{
    sapply(values, dist@cdfs[[var]])
})

get.plot.points <- function(dist, var, lower.quantile, upper.quantile, n=512)
{
    x = seq(dist@quantile.functions[[var]](lower.quantile),
            dist@quantile.functions[[var]](upper.quantile),
            length=n)
    d=sapply(x, dist@pdfs[[var]])
    list(x=x, d=d/sum(d))
}

map.values <- function(values, dist)
{
    if (!dist@use.var.names|| is.null(names(values)))
    {
        if (length(values) != dist@num.vars)
            stop(paste0("'values' should have one value for each of the ", dist@num.vars, " distributions. Instead, it has ", length(values), " elements."))
        else
            values
    }
    else
        values[dist@var.names]
}

#calls fn with val as the first argument, params as the subsequent args
do.dist.fn.call <- function(val, fn, params)
{
    do.call(fn, c(list(val), as.list(params)))
}

#creates a function that takes one parameter (val)
# the function returned calls the function with name paste0(fn.prefix, dist.name)
# with the first argument to it being val, and the subsequent arguments being params
#
#eg, a call to this function with arguments('norm', 'd', c(0,1))
# will return a function, say f
#calling f(val) will return the result of a call to dnorm(val, 0, 1)
create.dist.function <- function(dist.name, fn.prefix, params)
{
    p=params
    fn = get(paste0(fn.prefix, dist.name))
    function(val)
    {
        do.dist.fn.call(val, fn, p)
    }
}