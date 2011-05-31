setGeneric("get.subnets", function( model, get.names = TRUE, min.size = 2, max.size = Inf, min.responses = 2){ standardGeneric ("get.subnets") })
setGeneric("get.qofz", function( model, subnet.id, log = FALSE){ standardGeneric ("get.qofz") })
setGeneric("get.dat", function( model, subnet.id, sample = NULL){ standardGeneric ("get.dat") })
setGeneric("get.P.r", function(model, subnet.id, log = TRUE){ standardGeneric ("get.P.r") })
setGeneric("get.P.Sr", function(sample, model, subnet.id, log = TRUE ){ standardGeneric ("get.P.Sr") })
setGeneric("get.P.rs.joint", function(sample, model, subnet.id, log = TRUE ){ standardGeneric ("get.P.rs.joint") })
setGeneric("get.P.rS", function( model, subnet.id, log = TRUE){ standardGeneric ("get.P.rS") })
setGeneric("get.P.rs.joint.individual", function(sample, model, subnet.id, log = TRUE ){ standardGeneric ("get.P.rs.joint.individual") })
setGeneric("get.P.s.individual", function(sample, model, subnet.id, log = TRUE ){ standardGeneric ("get.P.s.individual") })
setGeneric("sample.densities", function( sample, model, subnet.id, log = TRUE, summarize = FALSE){ standardGeneric ("sample.densities") })
setGeneric("get.P.s", function( sample, model, subnet.id, log = TRUE){ standardGeneric ("get.P.s") })
setGeneric("get.P.rs", function(model, subnet.id, log = FALSE ){ standardGeneric ("get.P.rs") })

#setGeneric("response.enrichment", function(model, subnet.id, log = FALSE ){ standardGeneric ("response.enrichment") })
#setGeneric("plot.associations", function( x, subnet.id, labels, method = "hypergeometric", ...){ standardGeneric ("plot.associations") })
#setGeneric("plot.pca", function( x, subnet.id, labels = NULL, confidence = 0.95, ...){ standardGeneric ("plot.pca") })
#setGeneric("get.model.parameters", function( model, ... ){ standardGeneric ("get.model.parameters") })                                        


