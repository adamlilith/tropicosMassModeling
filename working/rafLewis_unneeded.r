#' Raftery-Lewis procedure for determining iterations needed for convergence
#'
#' This function is \emph{exactly} the same as the code in the \code{\link[coda]{raftery.diag}} function in the \pkg{coda} package except that it anticipates some issues related to the \code{thin} function when adding additional iterations to a chain.
#' @param data An object of class \code{mcmc}.
#' @param thin Thinning interval.
#' @param q A quantile to be estimated.
#' @param r The precision (margin of error) with which to estimate \code{q}.
#' @param s Probability of obtaining an estimate in the interval (\code{q - r}, \code{q + r}).
#' @param converge.eps Precision for the estimate of time to converge.
#' @details Please see the help for \code{\link[coda]{raftery.diag}}.
#' @export

rafLewis <- function(
	data,
	thin = 1,
	q = 0.025,
	r = 0.005,
	s = 0.95,
	converge.eps = 0.001
) {
    if (is.mcmc.list(data)) 
        return(lapply(data, raftery.diag, q, r, s, converge.eps))
    data <- as.mcmc(data)
    resmatrix <- matrix(nrow = nvar(data), ncol = 4, dimnames = list(varnames(data, 
        allow.null = TRUE), c("M", "N", "Nmin", 
        "I")))
    phi <- qnorm(0.5 * (1 + s))
    nmin <- as.integer(ceiling((q * (1 - q) * phi^2)/r^2))
    if (nmin > niter(data)) 
        resmatrix <- c("Error", nmin)
    else for (i in 1:nvar(data)) {
        if (is.matrix(data)) {
            # quant <- quantile(data[ , i, drop = TRUE], probs = q)
            quant <- quantile(data[ , i, drop = TRUE], probs = q)
            # dichot <- mcmc(data[, i, drop = TRUE] <= quant, start = start(data), 
                # end = end(data), thin = thin(data))
            dichot <- mcmc(data[ , i, drop = TRUE] <= quant, start = start(data), 
                end = end(data), thin = thin)
        }
        else {
            quant <- quantile(data, probs = q)
            dichot <- mcmc(data <= quant, start = start(data), 
                # end = end(data), thin = thin(mcmc))
                end = end(data), thin = thin)
        }
        kthin <- 0
        bic <- 1
        while (bic >= 0) {
            kthin <- kthin + thin(data)
            testres <- as.vector(window.mcmc(dichot, thin = kthin))
            testres <- factor(testres, levels = c(FALSE, TRUE))
            newdim <- length(testres)
            testtran <- table(testres[1:(newdim - 2)], testres[2:(newdim - 
                1)], testres[3:newdim])
            testtran <- array(as.double(testtran), dim = dim(testtran))
            g2 <- 0
            for (i1 in 1:2) {
                for (i2 in 1:2) {
                  for (i3 in 1:2) {
                    if (testtran[i1, i2, i3] != 0) {
                      fitted <- (sum(testtran[i1, i2, 1:2]) * 
                        sum(testtran[1:2, i2, i3]))/(sum(testtran[1:2, 
                        i2, 1:2]))
                      g2 <- g2 + testtran[i1, i2, i3] * log(testtran[i1, 
                        i2, i3]/fitted) * 2
                    }
                  }
                }
            }
            bic <- g2 - log(newdim - 2) * 2
        }
        finaltran <- table(testres[1:(newdim - 1)], testres[2:newdim])
        alpha <- finaltran[1, 2]/(finaltran[1, 1] + finaltran[1, 
            2])
        beta <- finaltran[2, 1]/(finaltran[2, 1] + finaltran[2, 
            2])
        tempburn <- log((converge.eps * (alpha + beta))/max(alpha, 
            beta))/(log(abs(1 - alpha - beta)))
        nburn <- as.integer(ceiling(tempburn) * kthin)
        tempprec <- ((2 - alpha - beta) * alpha * beta * phi^2)/(((alpha + 
            beta)^3) * r^2)
        nkeep <- as.integer(ceiling(tempprec) * kthin)
        iratio <- (nburn + nkeep)/nmin
        resmatrix[i, 1] <- nburn
        resmatrix[i, 2] <- nkeep + nburn
        resmatrix[i, 3] <- nmin
        resmatrix[i, 4] <- signif(iratio, digits = 3)
    }
    y <- list(params = c(r = r, s = s, q = q), resmatrix = resmatrix)
    class(y) <- "raftery.diag"
    return(y)
}
