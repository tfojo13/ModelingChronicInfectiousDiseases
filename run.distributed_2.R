source('priors.R')
source('models.R')
source('likelihoods.R')
source('inference.R')
source('plots.R')

NUM.SAMPLES = 10000
NUM.RESAMPLES = NUM.SAMPLES * 5



prepare.distributed.sampling('sampling/model1', NUM.SAMPLES, 1, 
                             num.input.parameters = NUM.PARAMS.1,
                             generate.input.parameters=generate.uniform.params, 
                             run.model = run.model, 
                             rscript.call="C:/Program Files/R/R-3.3.1/bin/x64/Rscript.exe",
                             batch.extension='bat',
                             model.no=1)
run.distributed.samples('sampling/model1')
augment.distributed.sampling('sampling/model1', NUM.SAMPLES)
run.distributed.samples('sampling/model1', groups=2)

r1.res.1 = resample.distributed.samples('sampling/model1', 1000, gaussian.incidence.likelihood, years=c(2000,2015))
r1.1 = get.outputs(r1.resamples)
plot.results(r1.1)
plot.results(r1.1, incidence = F)

r1.res.3 = resample.distributed.samples('sampling/model1', 1000, gaussian.incidence.likelihood, years=c(2000, 2008, 2015))
r1.3 = get.outputs(r1.res.3)
plot.results(r1.3)

r1.res.all = resample.distributed.samples('sampling/model1', 1000, gaussian.incidence.likelihood, years=2000:2015)
r1.all = get.outputs(r1.res.all)
plot.results(r1.all)
plot.results(r1.all, incidence = F)


prepare.distributed.sampling('sampling/model2', NUM.SAMPLES, 1, 
                             num.input.parameters = NUM.PARAMS.2,
                             generate.input.parameters=generate.uniform.params, 
                             run.model = run.model, 
                             rscript.call="C:/Program Files/R/R-3.3.1/bin/x64/Rscript.exe",
                             batch.extension='bat',
                             model.no=2)
run.distributed.samples('sampling/model2')
augment.distributed.sampling('sampling/model2', NUM.SAMPLES)
run.distributed.samples('sampling/model2', groups=2)

r2.resamples = resample.distributed.samples('sampling/model2', 1000, gaussian.incidence.likelihood, years=c(2000,2015))
r2 = get.outputs(r2.resamples)
plot.results(r2)

prepare.distributed.sampling('sampling/model3', 10000, 1, 
                             num.input.parameters = NUM.PARAMS.3,
                             generate.input.parameters=generate.uniform.params, 
                             run.model = run.model, 
                             rscript.call="C:/Program Files/R/R-3.3.1/bin/x64/Rscript.exe",
                             batch.extension='bat',
                             model.no=3)
run.distributed.samples('sampling/model3')
augment.distributed.sampling('sampling/model3', NUM.SAMPLES)
run.distributed.samples('sampling/model3', groups=2)
augment.distributed.sampling('sampling/model3', 80000)
run.distributed.samples('sampling/model3', groups=3)

r3.u.ind2 = resample.and.get.output('sampling/model3', independent.gaussian.incidence.likelihood, years=c(2000,2015))
r3.u.indall = resample.and.get.output('sampling/model3', independent.gaussian.incidence.likelihood, years=2000:2015)
r3.u.ar9 = resample.and.get.output('sampling/model3', mvg.ar.incidence.likelihood, years=2000:2015, rho=0.9)
r3.u.cs9 = resample.and.get.output('sampling/model3', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.9)
r3.u.cs5 = resample.and.get.output('sampling/model3', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
r3.u.u2 = resample.and.get.output('sampling/model3', uniform.incidence.likelihood, years=c(2000,2015))
r3.u.uall = resample.and.get.output('sampling/model3', uniform.incidence.likelihood, years=2000:2015)

r3.u.ind2.m = resample.and.get.output('sampling/model3', independent.gaussian.mort.inc.likelihood, years=c(2000,2015))
r3.u.indall.m = resample.and.get.output('sampling/model3', independent.gaussian.mort.inc.likelihood, years=2000:2015)
r3.u.ar9.m = resample.and.get.output('sampling/model3', mvg.ar.mort.inc.likelihood, years=2000:2015, rho=0.9)
r3.u.cs9.m = resample.and.get.output('sampling/model3', mvg.cs.mort.inc.likelihood, years=2000:2015, rho=0.9)
r3.u.cs5.m = resample.and.get.output('sampling/model3', mvg.cs.mort.inc.likelihood, years=2000:2015, rho=0.5)
r3.u.u2.m = resample.and.get.output('sampling/model3', uniform.mort.inc.likelihood, years=c(2000,2015))
r3.u.uall.m = resample.and.get.output('sampling/model3', uniform.mort.inc.likelihood, years=2000:2015)

#r3.u.binom2 = resample.and.get.output('sampling/model3', independent.binomial.incidence.likelihood, years=c(2000,2015))
#r3.u.binomall = resample.and.get.output('sampling/model3', independent.binomial.incidence.likelihood, years=2000:2015)

save(r3.u.ind2, r3.u.indall, r3.u.ar9, r3.u.cs9, r3.u.cs5, r3.u.u2, r3.u.uall,
     r3.u.ind2.m, r3.u.indall.m, r3.u.ar9.m, r3.u.cs9.m, r3.u.cs5.m, r3.u.u2.m, r3.u.uall.m,
     file='saved/uniform_prior_results.Rdata')

prepare.distributed.sampling('sampling/model3_log', 10000, 1, 
                             num.input.parameters = NUM.PARAMS.3,
                             generate.input.parameters=generate.log.and.logit.norm.parameters, 
                             run.model = run.model, 
                             rscript.call="C:/Program Files/R/R-3.3.1/bin/x64/Rscript.exe",
                             batch.extension='bat',
                             model.no=3)
run.distributed.samples('sampling/model3_log')
augment.distributed.sampling('sampling/model3_log', NUM.SAMPLES)
run.distributed.samples('sampling/model3_log', groups=2)

r3.log.ind2 = resample.and.get.output('sampling/model3_log', independent.gaussian.incidence.likelihood, years=c(2000,2015))
r3.log.indall = resample.and.get.output('sampling/model3_log', independent.gaussian.incidence.likelihood, years=2000:2015)
r3.log.ar9 = resample.and.get.output('sampling/model3_log', mvg.ar.incidence.likelihood, years=2000:2015, rho=0.9)
r3.log.cs9 = resample.and.get.output('sampling/model3_log', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.9)
r3.log.cs5 = resample.and.get.output('sampling/model3_log', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
r3.log.u2 = resample.and.get.output('sampling/model3_log', uniform.incidence.likelihood, years=c(2000,2015))
r3.log.uall = resample.and.get.output('sampling/model3_log', uniform.incidence.likelihood, years=2000:2015)

#r3.log.ind2.m = resample.and.get.output('sampling/model3_log', independent.gaussian.mort.inc.likelihood, years=c(2000,2015))
#r3.log.indall.m = resample.and.get.output('sampling/model3_log', independent.gaussian.mort.inc.likelihood, years=2000:2015)
#r3.log.ar9.m = resample.and.get.output('sampling/model3_log', mvg.ar.mort.inc.likelihood, years=2000:2015, rho=0.9)
#r3.log.cs9.m = resample.and.get.output('sampling/model3_log', mvg.cs.mort.inc.likelihood, years=2000:2015, rho=0.9)
#r3.log.cs5.m = resample.and.get.output('sampling/model3_log', mvg.cs.mort.inc.likelihood, years=2000:2015, rho=0.5)
#r3.log.u2.m = resample.and.get.output('sampling/model3_log', uniform.mort.inc.likelihood, years=c(2000,2015))
#r3.log.uall.m = resample.and.get.output('sampling/model3_log', uniform.mort.inc.likelihood, years=2000:2015)

#r3.log.binom2 = resample.and.get.output('sampling/model3_log', independent.binomial.incidence.likelihood, years=c(2000,2015))
#r3.log.binomall = resample.and.get.output('sampling/model3_log', independent.binomial.incidence.likelihood, years=2000:2015)

save(r3.log.ind2, r3.log.indall, r3.log.ar9, r3.log.cs9, r3.log.cs5, r3.log.u2, r3.log.uall,
#     r3.log.ind2.m, r3.log.indall.m, r3.log.ar9.m, r3.log.cs9.m, r3.log.cs5.m, r3.log.u2.m, r3.log.uall.m,
     file='saved/log_prior_results.Rdata')


#lets resample with: ind2, indall, ar9, cs9, cs5, uniform

prepare.distributed.sampling('sampling/model3_unarrow', 10000, 1, 
                             num.input.parameters = NUM.PARAMS.3,
                             generate.input.parameters=generate.uniform.params.narrow, 
                             run.model = run.model, 
                             rscript.call="C:/Program Files/R/R-3.3.1/bin/x64/Rscript.exe",
                             batch.extension='bat',
                             model.no=3)
run.distributed.samples('sampling/model3_unarrow')
augment.distributed.sampling('sampling/model3_unarrow', NUM.SAMPLES)
run.distributed.samples('sampling/model3_unarrow', groups=2)

r3.unw.ind2 = resample.and.get.output('sampling/model3_unarrow', independent.gaussian.incidence.likelihood, years=c(2000,2015))
r3.unw.indall = resample.and.get.output('sampling/model3_unarrow', independent.gaussian.incidence.likelihood, years=2000:2015)
r3.unw.ar9 = resample.and.get.output('sampling/model3_unarrow', mvg.ar.incidence.likelihood, years=2000:2015, rho=0.9)
r3.unw.cs9 = resample.and.get.output('sampling/model3_unarrow', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.9)
r3.unw.cs5 = resample.and.get.output('sampling/model3_unarrow', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
r3.unw.u2 = resample.and.get.output('sampling/model3_unarrow', uniform.incidence.likelihood, years=c(2000,2015))
r3.unw.uall = resample.and.get.output('sampling/model3_unarrow', uniform.incidence.likelihood, years=2000:2015)

save(r3.unw.ind2, r3.unw.indall, r3.unw.ar9, r3.unw.cs9, r3.unw.cs5, r3.unw.u2, r3.unw.uall,
     file='saved/uniform_narrow_prior_results.Rdata')

prepare.distributed.sampling('sampling/model3_lognarrow', 10000, 1, 
                             num.input.parameters = NUM.PARAMS.3,
                             generate.input.parameters=generate.log.and.logit.norm.parameters.narrow, 
                             run.model = run.model,  
                             rscript.call="C:/Program Files/R/R-3.3.1/bin/x64/Rscript.exe",
                             batch.extension='bat',
                             model.no=3)
run.distributed.samples('sampling/model3_lognarrow')
augment.distributed.sampling('sampling/model3_lognarrow', NUM.SAMPLES)
run.distributed.samples('sampling/model3_lognarrow', groups=2)

r3.lognw.ind2 = resample.and.get.output('sampling/model3_lognarrow', independent.gaussian.incidence.likelihood, years=c(2000,2015))
r3.lognw.indall = resample.and.get.output('sampling/model3_lognarrow', independent.gaussian.incidence.likelihood, years=2000:2015)
r3.lognw.ar9 = resample.and.get.output('sampling/model3_lognarrow', mvg.ar.incidence.likelihood, years=2000:2015, rho=0.9)
r3.lognw.cs9 = resample.and.get.output('sampling/model3_lognarrow', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.9)
r3.lognw.cs5 = resample.and.get.output('sampling/model3_lognarrow', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
r3.lognw.u2 = resample.and.get.output('sampling/model3_lognarrow', uniform.incidence.likelihood, years=c(2000,2015))
r3.lognw.uall = resample.and.get.output('sampling/model3_lognarrow', uniform.incidence.likelihood, years=2000:2015)

save(r3.lognw.ind2, r3.lognw.indall, r3.lognw.ar9, r3.lognw.cs9, r3.lognw.cs5, r3.lognw.u2, r3.lognw.uall,
     file='saved/log_narrow_prior_results.Rdata')



prepare.distributed.sampling('sampling/test', 2500, 1, 
                             num.input.parameters = NUM.PARAMS.3,
                             generate.input.parameters=generate.uniform.params, 
                             run.model = run.model, 
                             rscript.call="C:/Program Files/R/R-3.3.1/bin/x64/Rscript.exe",
                             batch.extension='bat',
                             model.no=3)
run.distributed.samples('sampling/test')
augment.distributed.sampling('sampling/test', 2500, 1)
run.distributed.samples('sampling/test', groups=2)

test.cs5 = resample.and.get.output('sampling/test', mvg.cs.incidence.likelihood, years=2000:2015, rho=0.5)
plot.results(test.cs5)
