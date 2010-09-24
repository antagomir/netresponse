#setMethod("get.model.parameters", "NetResponseModel",
#          function (model, subnet.idx, level = NULL) {})

setMethod(f = "[[", signature("NetResponseModel"),
   definition = (function(x, i, j = "missing", ..., exact = TRUE) {
      if (typeof(i) == "numeric"){i <- names(x)[[i]]}
      #new("rpa.list", list(d = x$d[i,], sigma2 = x$sigma2[[i]], cind = x$cind, set = x$sets[[i]]))
      get.model.parameters(x, subnet.id = i, level = j)
   })
)

#setReplaceMethod(f="[[",signature("ChromosomeArmModels"),
#                                definition=(function(x,i,j,value) {
#                                        x@models[[i]] <- value
#                                        return(x)
#                                }
#))
 
