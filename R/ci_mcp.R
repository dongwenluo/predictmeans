#' Multiple Comparisons Based on the Confidence Intervals
#'
#' This function produces letter representations for a multiple comparison test
#' by analyzing the confidence intervals associated with the mean values of
#' different treatments. In particular, if the confidence intervals of two
#' treatments overlap, it indicates that there is no significant difference
#' between them. Conversely, if the confidence intervals do not overlap, it
#' indicates that the treatments are significantly different from each other.
#'
#'
#' @param LL Lower limits of treatments' confidence interval.
#' @param UL Upper limits of treatments' confidence interval.
#' @param trt_n Treatments' names.
#' @author Dongwen Luo, Siva Ganesh and John Koolaard
#' @references Vanessa, C. (05 October 2022), \emph{Confidence tricks: the
#' 83.4\% confidence interval for comparing means},
#' https://vsni.co.uk/blogs/confidence_trick.
#' @examples
#'
#' library(predictmeans)
#' ci_mcp(LL=c(68.2566,  87.7566, 103.0899, 112.2566), UL=c(90.5212, 110.0212, 125.3545, 134.5212))
#' #--------------------------------------------------------
#'   data("Oats", package="nlme")
#'   Oats$nitro <- factor(Oats$nitro)
#'   fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
#'   predictmeans(fm, "nitro", adj="BH", plot=FALSE)$mean_table
#'   predictmeans(fm, "nitro", level=0.05, letterCI = TRUE, plot=FALSE)$mean_table
#'   predictmeans(fm, "nitro", level=0.01, letterCI = TRUE, plot=FALSE)$mean_table
#'   predictmeans(fm, "nitro", level=0.1, letterCI = TRUE, plot=FALSE)$mean_table
#'
#' @export

ci_mcp <- function(LL, UL, trt_n = NULL) {
  stopifnot("Check your LL and UL input!" = {
    is.numeric(LL)
    is.numeric(UL)
    length(LL) == length(UL)
    all(LL <= UL)
  })
  trt_len <- length(LL)
  if (is.null(trt_n) || length(unique(trt_n)) != trt_len) {
    trt_n <- as.character(1:trt_len)
  }

  ci_mcp_letters_0 <- rep("A", trt_len)
  names(ci_mcp_letters_0) <- trt_n

  results <- matrix(NA_real_, nrow = trt_len, ncol = trt_len)

  for (i in 1:trt_len) {
    for (j in (i + 1):trt_len) {
      if (j > trt_len) {
        break
      }
      ci1 <- c(LL[i], UL[i])
      ci2 <-  c(LL[j], UL[j])
      if (max(ci1) < min(ci2) || max(ci2) < min(ci1)) {
        results[i, j] <- 0.01
      } else {
        results[i, j] <- 0.08
      }
    }
  }

  if (all(unique(na.omit(as.vector(results))) == 0.08)) {
    ci_mcp_letters <- ci_mcp_letters_0
  } else {
    rownames(results) <- colnames(results) <- trt_n
    results[lower.tri(results)] <- t(results)[lower.tri(results)]
    ci_mcp_letters <- multcompLetters(results, Letters = LETTERS)
  }
  return(ci_mcp_letters)
}
