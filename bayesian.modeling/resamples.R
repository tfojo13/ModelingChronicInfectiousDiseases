##------------------------------------------------------------##
## Contains the code and documentation for the resamples
## class and associated methods.
##------------------------------------------------------------##

##-- THE CLASS DEFINITION --##

#'@title The Resamples Class
#'
#'@description An S4 class which represents a 'resampled set' - the results of a importance-resampling algorithm (ie, the results of a call to \code{\link{resample.distributed.samples}})\cr
#'
#'@slot counts A vector of counts representing the number of times each unique model result was resampled (ie, @counts[i] represents the number of times the ith sample was resampled)
#'@slot results A list of unique results objects in the resampled set, where each element is the results as returned by run.model argument passed to \code{\link{prepare.distributed.sampling}}
#'@slot likelihoods A vector likelihoods, where likelihood[i] corresponds to the likelihood for results[[i]]
#'@slot input.parameters A matrix where the ith row represents the input parameters used to generate the ith sample
#'@slot num.original.samples The original number of samples from which this resample set was drawn
#'@slot num.resamples The number of samples (with replacement) that make up this resample set
#'@slot num.unique.resamples The number of UNIQUE samples that constitute the resample set
#'
#'@export
setClass('Resamples',
         representation = list(
             counts = 'integer',
             results = 'list',
             likelihoods = 'numeric',
             input.parameters = 'matrix',
             num.original.samples = 'integer',
             num.resamples = 'integer',
             num.unique.resamples = 'integer',
             maximum.likelihood = 'numeric',
             marginal.likelihood = 'numeric'
         ))

##-- THE (private) CONSTRUCTOR --##

Resamples <- function(counts, results, likelihoods, input.parameters, num.original.samples, max.likelihood, marginal.likelihood)
{
    if (length(counts)==1)
    {
        new('Resamples',
            counts=as.integer(counts),
            results=results,
            likelihoods=likelihoods,
            input.parameters=input.parameters,
            num.original.samples=as.integer(num.original.samples),
            num.resamples=as.integer(sum(counts)),
            num.unique.resamples=as.integer(length(counts)),
            maximum.likelihood=max.likelihood,
            marginal.likelihood=marginal.likelihood)
    }
    else
    {
        o = order(counts, decreasing=T)
        
        new('Resamples',
            counts=as.integer(counts[o]),
            results=results[o],
            likelihoods=likelihoods[o],
            input.parameters=input.parameters[o,],
            num.original.samples=as.integer(num.original.samples),
            num.resamples=as.integer(sum(counts)),
            num.unique.resamples=as.integer(length(counts)),
            maximum.likelihood=max.likelihood,
            marginal.likelihood=marginal.likelihood)
    }
}

##-- SOME METHODS --##

setMethod('show', 'Resamples', def=function(object)
{
    do.summarize.resamples(object, num.to.show=5)
})

setMethod('length', 'Resamples', def=function(x)
{
    length(x@counts)
})


##-- HELPERS --##

do.summarize.resamples <- function(resamples, contained.in=NULL, num.to.show=5)
{
    if (is.null(contained.in))
        title="A resample set consisting of:\n"
    else
        title=paste0("The resample set conatined in ", contained.in, " consists of:\n")
    
    cat(title)
    cat(paste0("-", resamples@num.resamples, " resampled model runs comprising\n"))
    cat(paste0("-", resamples@num.unique.resamples, " unique samples from the original ", resamples@num.original.samples, " samples (", round(100*resamples@num.unique.resamples/resamples@num.original.samples,1), "%)\n"))
    
    cat("\n")
    
    cum.counts = cumsum(resamples@counts)
    
    threshold = 0.5 * resamples@num.resamples
    threshold.index = (1:resamples@num.unique.resamples)[cum.counts>threshold][1]
    pct.threshold = round(100*threshold.index/resamples@num.unique.resamples)
    if (threshold.index==1)
        start = paste0("The single (",pct.threshold,"%) most resampled sample accounts for")
    else
        start = paste0("The ", threshold.index, " / ", resamples@num.unique.resamples, " (",pct.threshold,"%) most resampled samples account for")
    cat(start, paste0(round(100*cum.counts[threshold.index]/resamples@num.resamples,1), '% (', cum.counts[threshold.index], " / ", resamples@num.resamples, ") of the resamples\n"))
    
    threshold = 0.9 * resamples@num.resamples
    threshold.index = (1:resamples@num.unique.resamples)[cum.counts>threshold][1]
    pct.threshold = round(100*threshold.index/resamples@num.unique.resamples)
    if (threshold.index==1)
        start = paste0("The single (",pct.threshold,"%) most resampled sample accounts for")
    else
        start = paste0("The ", threshold.index, " / ", resamples@num.unique.resamples, " (",pct.threshold,"%) most resampled samples account for")
    cat(start, paste0(round(100*cum.counts[threshold.index]/resamples@num.resamples,1), '% (', cum.counts[threshold.index], " / ", resamples@num.resamples, ") of the resamples\n"))
    
    percent.resamples = round(100*resamples@counts/sum(resamples@counts),1)
    num.to.show = min(num.to.show, resamples@num.unique.resamples)
    if (num.to.show > 0)
    {
        for (i in 1:num.to.show)
        {
            if (i==1)
                start='-The most'
            else
                start=paste0('-The ', get.ordinal(i), ' most')
            cat(start, 'frequently resampled sample accounts for', paste0(percent.resamples[i],'% (',resamples@counts[i],')'), "of the resamples\n")
        }
    }
}

#converts a number to an ordinal
get.ordinal <- function(num)
{
    last.digit = num - 10*floor(num/10)
    if (last.digit==1)
        paste0(num, 'st')
    else if (last.digit==2)
        paste0(num, 'nd')
    else if (last.digit==3)
        paste0(num, 'rd')
    else
        paste0(num, 'th')
}

#'@title Calculate PCCs or PRCCs on a Resample Set
#'
#'@param resamples Can be either (1) An instance of \code{\link{Resamples-class}}, (2) a file where such a Resamples object is saved, or (3) a directory (such as one set up by a call to \code{link{prepare.distributed.sampling}} on which \code{\link{resample.distributed.samples}} has been called)
#'@param calculate.response A function that takes as its first argument a result item and as its second argument a input parameter set and returns the response variable to be used in calculating PCCs or PRCCs
#'@param ... Additional arguments to be passed to calculate.response
#'@param rank Boolean indicating whether to calculate PRCCs (TRUE) or PCCs (FALSE)
#'@param model.variable.names A character vector containing the names of the input.parameters to calculate PCCs for. If NULL, all input.parameters are used
#'
#'@seealso \code{\link{pcc}}
#'
#'@export

calculate.pccs <- function(resamples, calculate.response, ..., rank=T, model.variable.names=NULL)
{
    resamples = load.resamples(resamples)
    
    if (is.null(model.variable.names))
        model.variable.names = dimnames(resamples@input.parameters)[[2]]
    
    responses = sapply(1:resamples@num.unique.resamples, function(i){
        calculate.response(resamples@results[[i]], resamples@input.parameters[i,], ...)
    })
    
    df = as.data.frame(resamples@input.parameters[,model.variable.names])
    pcc(df, responses, rank=rank)
}