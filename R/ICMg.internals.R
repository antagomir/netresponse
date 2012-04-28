# Copyright (C) 2008-2011 Juuso Parkkinen
# Contact: Juuso Parkkinen <juuso.parkkinen@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

`ICMg.combined.iteration` <-
function(L, X, Niter, N, M, Lindices, Nindices,
                        D, C, z, w, n, m, q, E,
                        alpha, beta, pm0, V0, V,
                        convl, convn) {
  V0i <- V0^(-1)
  Vi <- V^(-1)
  
  ## Main iteration loop
  for (s in 1:Niter) {
    cat(".")
    
    ## Link loop
    for (li in 1:N) {
      
      l = Lindices[li]
      i = L[l,1]
      j = L[l,2]
      
      ## Subtract the contribution of the link from the counts
      q[z[l],i] <- q[z[l],i]-1
      q[z[l],j] <- q[z[l],j]-1
      n[z[l]] <- n[z[l]]-1

      ## Loop for computing probabilities for the components to be sampled
      uz <- vector("numeric", C)
      for (p in 1:C) {
        A <- (n[p] +m[p] +alpha)
        B <- (q[p,i] +beta)*(q[p,j] +beta) / ((2*n[p] +m[p] +1 +M*beta)*(2*n[p] +m[p] +M*beta))
        uz[p] <- A *B
      }
            
      ## Draw a new component for the links and update the counts */
      newz <- ICMg.multinom.single(uz)
      convl[l] <- uz[newz]/sum(uz)
      n[newz] <- n[newz]+1
      q[newz,i] <- q[newz,i]+1
      q[newz,j] <- q[newz,j]+1
      z[l] <- newz
    }

    ## Node loop
    for (ki in 1:M) {
      
      k = Nindices[ki]
      
      ## Subtract the contribution of the node from the counts
      q[w[k],k] <- q[w[k],k]-1
      m[w[k]] <- m[w[k]]-1
      E[w[k], ] <- E[w[k], ] - X[k, ]

      ## Loop for computing probabilities for the components to be sampled
      uw <- vector("numeric", C)
      AB <- vector("numeric", C)
      H1 <- vector("numeric", C)
      He <- vector("numeric", C)
      
      for (p in 1:C) {
        
        AB[p] <- (n[p] +m[p] +alpha) * (q[p,k] +beta) / (2*n[p] +m[p] +(M-1)*beta)

        S <- (V0i + (m[p]+1)*Vi)^(-1)
        Sd <-(V0i + m[p]*Vi)^(-1)     
        
        A <- as.vector(S * (V0i * pm0 + Vi * E[p, ] + Vi * X[k, ]))
        Ad <- as.vector(Sd * (V0i * pm0 + Vi * E[p, ]))
      
        H1[p] <- ( S / Sd )^(D/2)  
        
        H31e <- crossprod(A) / S
        H32e <- crossprod(Ad) / Sd
        He[p] <- 1/2*(H31e -H32e)
      }

      for (p in 1:C) 
        uw[p] <- AB[p] *H1[p] * exp(He[p] - max(He)) 
      
      ## Draw a new component for the links and update the counts */
      neww <- ICMg.multinom.single(uw)
      convn[k] <- uw[neww]/sum(uw)
      m[neww] <- m[neww]+1
      q[neww,k] <- q[neww,k]+1
      E[neww, ] <- E[neww, ] + X[k, ]
      w[k] <- neww
    }
  }
  cat("\n")
  
  return(list(z=z, w=w, n=n, m=m, q=q, E=E, convl=convl, convn=convn))
}

`ICMg.combined.wrapper` <-
function(L, X, Niter, N, M, Lindices, Nindices,
                          D, C, z, w, n, m, q, E,
                          alpha, beta, pm0, V0, V,
                          convl, convn, C.boost) {

  if (C.boost) {
    out <- .C("ICMgCombinedIteration",PACKAGE="netresponse",L=as.integer(L),X=as.double(X),
              Niter=as.integer(Niter),N=as.integer(N),
              M=as.integer(M),Lindices=as.integer(Lindices),
              Nindices=as.integer(Nindices),D=as.integer(D),
              C=as.integer(C),z=as.integer(z),w=as.integer(w),
              n=as.integer(n), m=as.integer(m),q=as.integer(q),
              E=as.double(E),alpha=as.double(alpha),
              beta=as.double(beta), pm0=as.double(pm0),
              V0=as.double(V0), V=as.double(V),
              convl=as.double(convl),convn=as.double(convn))
  } else {
    out <- ICMg.combined.iteration(L, X, Niter, N, M, Lindices, Nindices,
                        D, C, z, w, n, m, q, E,
                        alpha, beta, pm0, V0, V,
                        convl, convn)
  }
  return(out)
}

`ICMg.links.iteration` <-
function(L, Niter, N, M, Lindices,
                       C, z, q, n, alpha, beta, conv) {
  ## Main iteration loop
  for (s in 1:Niter) {
    cat(".")
    ## Sample new component for each link
    for (li in 1:N) {

      l = Lindices[li]
      i = L[l,1]
      j = L[l,2]
      
      ## Subtract the contribution of the link from the counts
      q[z[l],i] <- q[z[l],i]-1
      q[z[l],j] <- q[z[l],j]-1
      n[z[l]] <- n[z[l]]-1

      ## Loop for computing probabilities for the components to be sampled
      uz <- vector("numeric", C)
      for (p in 1:C) {
        A <- (n[p] +alpha)
        B <- (q[p,i] +beta)*(q[p,j] +beta) / ((2*n[p] +1 +M*beta)*(2*n[p] +M*beta))
        uz[p] <- A *B
      }
            
      ## Draw a new component for the links and update the counts */
      newz <- ICMg.multinom.single(uz)
      conv[l] <- uz[newz]/sum(uz)
      n[newz] <- n[newz]+1
      q[newz,i] <- q[newz,i]+1
      q[newz,j] <- q[newz,j]+1
      z[l] <- newz
    }
  }
  cat("\n")
  return(list(z=z, q=q, n=n, conv=conv))
}

`ICMg.links.wrapper` <-
function(L, Niter, N, M, Lindices,
                        C, z, q, n,
                        alpha, beta, conv, C.boost) {

  if (C.boost) {
    out <- .C("ICMgLinksIteration",PACKAGE="netresponse", L=as.integer(L),Niter=as.integer(Niter),
              N=as.integer(N),M=as.integer(M),
              Lindices=as.integer(Lindices),C=as.integer(C),
              z=as.integer(z),q=as.integer(q),
              n=as.integer(n),alpha=as.double(alpha),
              beta=as.double(beta), conv=as.double(conv))
  } else {
    out <- ICMg.links.iteration(L, Niter, N, M, Lindices,
                      C, z, q, n, alpha, beta, conv)
  }
  return(out)
}

`ICMg.multinom.single` <-
function(prob) {
  cs <- cumsum(prob)
  which.max(runif(1) <= cs/cs[length(cs)])
}

`ICMg.randominit` <-
function() {
  .C("ICMgRandominit", PACKAGE="netresponse")
}

