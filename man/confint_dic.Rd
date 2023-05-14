\name{confint_dic}
\alias{confint_dic}
\title{
Confidence interval
}
\description{
Confidence interval
}
\usage{
confint_dic(out_em, prob = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out_em}{
result of \code{eccdf_dic_em}
}
  \item{prob}{
confidence level
}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (out_em, prob = 0.95) 
{
    ind = with(out_em, estimates$event > 0)
    variance <- diag(solve(out_em$I[ind, ind]))
    z <- qnorm(1 - 0.5 * (1 - prob))
    value = out_em$estimates$value[ind]
    ccdf = out_em$estimates$ccdf[ind]
    b <- sapply(z, function(z0) rev(z0 * sqrt(cumsum(rev(variance)))))
    lower = as.matrix(ifelse(ccdf - b < 0, 0, ccdf - b))
    upper = as.matrix(ifelse(ccdf + b > 1, 1, ccdf + b))
    res <- lapply(1:length(prob), function(i) data.frame(value = value, 
        ccdf = ccdf, se = b[, i], lower = lower[, i], upper = upper[, 
            i], level = prob[i]))
    res <- do.call("rbind", res)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
