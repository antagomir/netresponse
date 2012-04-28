# Copyright (C) 2008-2012 Olli-Pekka Huovilainen and Leo Lahti 
# Contact: Leo Lahti <leo.lahti@iki.fi>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
2# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

setGeneric("get.subnets", function( model, get.names = TRUE, min.size = 2, max.size = Inf, min.responses = 2){ standardGeneric ("get.subnets") })
setGeneric("getqofz", function( model, subnet.id, log = FALSE){ standardGeneric ("getqofz") })
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