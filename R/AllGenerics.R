#' @importFrom rmarkdown render

setGeneric("get.subnets", function(model, get.names = TRUE,
    min.size = 2, max.size = Inf, 
    min.responses = 2) {
    standardGeneric("get.subnets")
})
setGeneric("getqofz", function(model, subnet.id, log = FALSE) {
    standardGeneric("getqofz")
})
setGeneric("get.dat", function(model, subnet.id, sample = NULL) {
    standardGeneric("get.dat")
})
setGeneric("get.P.r", function(model, subnet.id, log = TRUE) {
    standardGeneric("get.P.r")
})
setGeneric("get.P.Sr", function(sample, model, subnet.id, log = TRUE) {
    standardGeneric("get.P.Sr")
})
setGeneric("get.P.rs.joint", function(sample, model, subnet.id, log = TRUE) {
    standardGeneric("get.P.rs.joint")
})
setGeneric("get.P.rS", function(model, subnet.id, log = TRUE) {
    standardGeneric("get.P.rS")
})
setGeneric("get.P.rs.joint.individual", function(sample, model,
    subnet.id, log = TRUE) {
    standardGeneric("get.P.rs.joint.individual")
})
setGeneric("get.P.s.individual", function(sample, model,
    subnet.id, log = TRUE) {
    standardGeneric("get.P.s.individual")
})
setGeneric("sample.densities", function(sample, model, subnet.id,
    log = TRUE, summarize = FALSE) {
    standardGeneric("sample.densities")
})
setGeneric("get.P.s", function(sample, model, subnet.id, log = TRUE) {
    standardGeneric("get.P.s")
})
setGeneric("get.P.rs", function(model, subnet.id, log = FALSE) {
    standardGeneric("get.P.rs")
})

