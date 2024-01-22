

function bw_tll_nn(udata, deg)


    zdata = qnorm(udata);
    n = size(zdata,1);
    d = size(zdata,2);
%    transform to uncorrelated data
    [pca qrs] = princomp(zdata);
%     B <- unclass(pca$loadings)
    B = B * sign(diag(B));
%     qrs <- unclass(pca$scores)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LATER 


%     ## find optimal alpha in for each principal component
    alphsq <- seq(size(zdata,2)^(-1/5), 1, l = 50)
    
    opt <- function(i) {
        val <- lscvplot(~qrs[, i],
                        alpha  = alphsq,
                        deg    = deg,
                        kern   = "gauss",
                        maxk   = 512)$value
        mean(alphsq[which.min(val)])
    }
    
    alpha.vec <- sapply(1:d, opt)

%     ## adjustments for multivariate estimation and transformation
    kappa = alpha.vec[1]/alpha.vec;
    
    dimnames(B) <- NULL
    
    if (deg == 1) 
        alpha <- n^(1/5 - d/(4 + d)) * alpha.vec[1]
    else 
        alpha <- n^(1/9 - d/(8 + d)) * alpha.vec[1]
    end

%     ## return results
    list(B = B, alpha = alpha, kappa = kappa)
