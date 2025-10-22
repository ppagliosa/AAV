#' @title Create a plot showing the estimated magnitude of an ecological impact
#'
#' @description
#' It provides a customised 'ggplot' for displaying the magnitude of an ecological
#' impact based on the results from \code{\link[AAV]{magnitude}}.
#'
#' @param data a list genereted from \code{\link[AAV]{magnitude}}.
#' @param tail only for \code{aci} post-impact design with only one disturbed location
#'     and multiple control locations.
#'     Perform an \code{upper} tailed, \code{lower} tailed or both join
#'     (\code{upper_lower}) plot for all terms of the \code{aav} analyses.
#'     Default is \code{upper}.
#'
#' @details
#' The magnitude of an ecological impact is estimated using the same data types, designs, factors
#' and parameters as in \code{\link[AAV]{aav}}, but only considers one variable at a time.
#'
#' The estimates are plotted to each relevant term from the asymmetrical analysis of variance
#' in order to determine the magnitude of the impact.
#'
#' The p-values in \code{aav} analyses came from an upper-tail test, except
#' wheater the analyses for \code{aci} post-impact design with only one disturbed location
#' are performed for a two-tailed test (i.e. aci.tailed.test = 2).
#' In this case, a plot showing an upper-tailed, lower-tailed or both types of tail
#' can be produced using the \code{tail} parameter.
#'
#' The terms can be interpreted separated each other.
#' To each term there is a chart showing the observed p-values from the original data,
#' and the expected p-value from a random data (resampled the original data)
#' and from each estimated alteration on the mean and variance of the data.
#' The magnitude of the impact can be discovered when the original p-value from a term
#' is significant (p < 0.05). In this case, the first estimate of the altered mean and variance
#' that is significant indicate the minimum magnitude of the impact that the original dataset can detect.
#'
#' Attention should be paid to whether the random data is significant.
#' In these cases, some spurious significant p-values may occur.
#'
#' Usually, lower significance (e.g. p < 0.05) can be seen flowing from the bottom right to the bottom left side of the chart.
#' Conversely, when visualising lower-tail tests the lower significance can be seen flowing from the top right to the top left of the chart.
#' In the cases of spurious significance, the lower significance can be seen flowing from the bottom left to the bottom right.
#'
#' Wheather a term is significant, it could detect:
#'
#'
#' In baci desing:
#' * Before vs After x Impact vs Controls = B x I:   Long-term local-scale impact
#' * B x Site(I) = B x S(I):   Long-term small-scale impact
#' * Perido(After) x I = P(Aft) x I:   Medium-term local-scale impact
#' * P(Aft) x S(I):   Medium-term small-scale impact
#' * Time(P(Aft)) x I = T(P(Aft)) x I or T(Aft) x I:   Short-term local-scale impact
#' * T(P(Aft)) x S(I) or T(Aft) x S(I):   Short-term small-scale impact
#'
#' In post-impact design with only one disturbed location and multiple control locations:
#' * Impact vs Controls = I:   Long-term local-scale impact
#' * Site(I) = S(I):   Long-term small-scale impact
#' * Period x I = P x I:   Medium-term local-scale impact
#' * P x S(I):   Medium-term small-scale impact
#' * T(P) x I or T x I:   Short-term local-scale impact
#' * T(P) x S(I) or T x S(I):   Short-term and small-scale impact
#'
#' In post-impact design involving more than one place with one disturbed
#' location and multiple control locations:
#' * Treatment = T:   Generalized local-scale impact
#' * Site(T) = S(T):   Generalized site-scale impact
#' * Place x T = P x T:   Local-specific scale impact
#' * P x S(T):   Site-specific scale impact
#'
#' @returns
#' The function returns a xy plot with the original (top right), a random (= 0; bottom left)
#' and the estimated alterations on the mean and variance of the data (center).
#'
#' @author Paulo Pagliosa \email{paulo.pagliosa@ufsc.br}
#'
#' @references
#'
#' Underwood, A.J. 1997. \emph{Experiments in Ecology: Their Logical Design and Interpretation Using
#' Analysis of Variance}. Cambridge University Press, Cambridge. \doi{10.1017/CBO9780511806407}.
#'
#' @seealso
#' \code{\link[AAV]{aav}},
#' \code{\link[GAD]{magnitude}}
#'
#' @examples
#' # Example for BACI design containing 2 temporal factors and 1 spatial factor.
#' # using the vectors:
#' # magnitude_md = c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1) and
#' # magnitude_var = c(0.1, 1, 2, 5, 10, 25, 50, 75, 100)
#' library(ggplot2)
#' data(ex.ggmag)
#' ggmagnitude(data = ex.ggmag)
#'
#' @export

ggmagnitude <- function(data, tail = "upper") {
    data<-data$extract_aav
    ad<- unique(data$aav.design)
    nft<-unique(data$n.ftemp)
    nt<-length(unique(data$term))
    size<-ifelse(nt == 2, 4.5, ifelse(nt == 3, 3.5, ifelse(nt == 4, 3.1, ifelse(nt == 6, 3))))

    ifelse(ad == "aci" & nft != 0 & tail == "lower", data<-data[data$two.tail == "lower",],
           ifelse(ad == "aci" & nft != 0 & tail == "upper_lower", data<-data[data$two.tail == "upper_lower",],
                  data<-data[data$two.tail == "upper",]))

    data$magVAR <- ordered(data$magVAR, levels=c(unique(data$magVAR[!data$magVAR == "obs"]),"obs"))
    data$magMD <- ordered(data$magMD, levels=c(unique(data$magMD[!data$magMD == "obs"]),"obs"))
    if(ad == "baci" & nt == 4) data$term <- ordered(data$term, levels=c("B x I", "T(Aft) x I","B x S(I)","T(Aft) x S(I)"))
    if(ad == "aci" & nt == 4) data$term <- ordered(data$term, levels=c("I","T x I","S(I)","T x S(I)"))
    if(ad == "baci" & nt == 6) data$term <- ordered(data$term, levels=c("B x I", "P(Aft) x I","T(P(Aft)) x I",
                                                                        "B x S(I)","P(Aft) x S(I)",  "T(P(Aft)) x S(I)"))
    if(ad == "aci" & nt == 6) data$term <- ordered(data$term, levels=c("I","P x I","T(P) x I",
                                                                       "S(I)","P x S(I)","T(P) x S(I)"))

    data$threshold[data$p.value <= 0.055]<-"low"
    data$threshold[data$p.value > 0.055]<-"high"
    to_string <- as_labeller(c(`B x I` = "LONG-TERM LOCAL-SCALE IMPACT:  B x I",
                               `B x S(I)` = "LONG-TERM SMALL-SCALE IMPACT:  B x S(I)",
                               `P(Aft) x I` = "MEDIUM-TERM LOCAL-SCALE:  P(Aft) x I",
                               `P(Aft) x S(I)` = "MEDIUM-TERM SMALL-SCALE:  P(Aft) x S(I)",
                               `T(P(Aft)) x I` = "SHORT-TERM LOCAL-SCALE IMPACT:  T(P(Aft)) x I",
                               `T(P(Aft)) x S(I)` = "SHORT-TERM SMALL-SCALE IMPACT:  T(P(Aft)) x S(I)",
                               `T(Aft) x I` = "SHORT-TERM LOCAL-SCALE IMPACT:  T(Aft) x I",
                               `T(Aft) x S(I)` = "SHORT-TERM SMALL-SCALE IMPACT:  T(Aft) x S(I)",
                               `I` = "LONG-TERM LOCAL-SCALE IMPACT:  I",
                               `S(I)` = "LONG-TERM SMALL-SCALE IMPACT:  S(I)",
                               `P x I` = "MEDIUM-TERM LOCAL-SCALE:  P x I",
                               `P x S(I)` = "MEDIUM-TERM SMALL-SCALE:  P x S(I)",
                               `T(P) x I` = "SHORT-TERM LOCAL-SCALE IMPACT:  T(P) x I",
                               `T(P) x S(I)` = "SHORT-TERM SMALL-SCALE IMPACT:  T(P) x S(I)",
                               `T x I` = "SHORT-TERM LOCAL-SCALE IMPACT:  T x I",
                               `T x S(I)` = "SHORT-TERM SMALL-SCALE IMPACT:  T x S(I)",
                               `T` = "GENERALIZED LOCAL-SCALE IMPACT:  T",
                               `P x T` = "LOCAL-SPECIFIC SCALE IMPACT:  P x T",
                               `S(T)` = "GENERALIZED SITE-SCALE IMPACT:  S(T)",
                               `P x S(T)` = "SITE-SPECIFIC SCALE IMPACT:  P x S(T)"))

    p<- ggplot(data, aes(x = magMD, y = magVAR, fill = threshold)) +
      xlab("Magnitude of the impact on the observed data mean") +
      ylab("Magnitude of the impact on the observed data variance") +
      geom_tile() +
      geom_text(aes(label = round(p.value,2)), color = "white", size = size) +
      coord_fixed() +
      scale_fill_manual(values = c("low" = "#1a9850", "high" = "#d73027")) +
      theme_classic() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, face="bold"),
            strip.text = element_text(face = "bold", color = "black"),
            strip.background = element_rect(linetype = "solid",color = "black", linewidth = 0.5)) +
      facet_wrap(~term, labeller = to_string)
    print(p)
}
