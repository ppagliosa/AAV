#' @title Create a plot showing the estimated magnitude of an ecological impact
#'
#' @description
#' It provides a customised 'ggplot' for displaying the magnitude of an ecological
#' impact based on the results from \code{\link[AAV]{magnitude}}.
#'
#' @param data a list genereted from \code{\link[AAV]{magnitude}}.
#' @param aav.design \code{baci} for beyond Before/After-Control/Impact design or \code{aci} for Post-Impact design.
#' @param n.ftemp number of temporal factors.
#' @param n.fspac number of spatial factors.
#'
#' @details
#' The magnitude of an ecological impact is estimated using the same data types, designs, factors
#' and parameters as in \code{\link[AAV]{aav}}, but only considers one variable at a time.
#'
#' The estimates are plotted to each relevant term from the asymmetrical analysis of variance
#' in order to determine the magnitude of the impact.
#' The terms can be interpreted separated each other.
#'
#' To each term there is a chart showing the observed p-values from the original data,
#' and the expected p-value from a random data (resampled the original data)
#' and from each estimated alteration on the mean and variance of the data.
#' The magnitude of the impact can be discovered when the original p-value from a term
#' is significant (p < 0.05). In this case, the first estimate of the altered mean and variance
#' that is significant indicate the minimum magnitude of the impact that the original dataset can detect.
#'
#' Wheather a term is significant, it could detect:
#'
#'
#' In baci desing:
#' * Before vs After x Impact vs Controls = B x I:   Long-term local-scale impact
#' * B x Site(I) = B x S(I):   Long-term small-scale impact
#' * Perido(After) x I = P(Aft) x I:   Medium-term local-scale impact
#' * Time(P(Aft)) x I = T(P(Aft)) x I:   Medium-term small-scale impact
#' * P(Aft) x S(I) or T(Aft) x I:   Short-term local-scale impact
#' * T(P(Aft)) x S(I) or T(Aft) x S(I):   Short-term local-scale impact
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
#' ggmagnitude(data = ex.ggmag, aav.design = "baci", n.ftemp = 2, n.fspac = 1)
#'
#' @export

ggmagnitude <- function(data, aav.design, n.ftemp, n.fspac) {
  te<-data.frame(unique(data$extract_aav$impType[!data$extract_aav$impType == "obs"]), unique(data$extract_aav$term))

  for (i in 1:nrow(te)) {
    dfl<-data$extract_aav[data$extract_aav$impType == te[i,1] | (data$extract_aav$impType == "obs" & data$extract_aav$term == te[i,2]),]
    dfl$magVAR <- ordered(dfl$magVAR, levels=c(unique(dfl$magVAR[!dfl$magVAR == "obs"]),"obs"))
    dfl$magMD <- ordered(dfl$magMD, levels=c(unique(dfl$magMD[!dfl$magMD == "obs"]),"obs"))
    dfl$threshold[dfl$p.value <= 0.055]<-"low"
    dfl$threshold[dfl$p.value > 0.055]<-"high"

    p<- ggplot(dfl, aes(x = magMD, y = magVAR, fill = threshold)) +
      xlab("Magnitude of the impact on the observed data mean") +
      ylab("Magnitude of the impact on the observed data variance") +
      geom_tile() +
      geom_text(aes(label = round(p.value,2)), color = "white", size = 4) +
      coord_fixed() +
      scale_fill_manual(values = c("low" = "#1a9850", "high" = "#d73027")) +
      theme_classic() +
      theme(legend.position = "none",
            panel.spacing = unit(-2, 'mm'),
            plot.title = element_text(hjust = 0.5, size=12, face="bold"),
            strip.text = element_text(face = "bold", color = "black", size = 10),
            strip.background = element_rect(linetype = "solid",color = "black", linewidth = 0.5))
    p<- p + list(
        if(aav.design == "baci" | aav.design == "aci" & n.ftemp != 0){
        list(if(te[i,1] == "L.imp") {labs(title = paste0("LONG-TERM LOCAL-SCALE IMPACT:  ", te[i,2]))},
             if(te[i,1] == "S.imp") {labs(title = paste0("LONG-TERM SMALL-SCALE IMPACT:  ", te[i,2]))},
             if(te[i,1] == "P.imp") {labs(title = paste0("MEDIUM-TERM LOCAL-SCALE:  ", te[i,2]))},
             if(te[i,1] == "SP.imp") {labs(title = paste0("MEDIUM-TERM SMALL-SCALE:  ", te[i,2]))},
             if(te[i,1] == "TP.imp") {labs(title = paste0("SHORT-TERM LOCAL-SCALE IMPACT:  ", te[i,2]))},
             if(te[i,1] == "T.imp") {labs(title = paste0("SHORT-TERM LOCAL-SCALE IMPACT:  ", te[i,2]))},
             if(te[i,1] == "ST.imp") {labs(title = paste0("SHORT-TERM SMALL-SCALE IMPACT:  ", te[i,2]))},
             if(te[i,1] == "STP.imp") {labs(title = paste0("SHORT-TERM SMALL-SCALE IMPACT:  ", te[i,2]))}
        )},
      if(aav.design == "aci" & n.ftemp == 0){
        list(if(te[i,1] == "L.imp") {labs(title = paste0("GENERALIZED LOCAL-SCALE IMPACT:  ", te[i,2]))},
             if(te[i,1] == "LA.imp") {labs(title = paste0("LOCAL-SPECIFIC SCALE IMPACT:  ", te[i,2]))},
             if(te[i,1] == "S.imp") {labs(title = paste0("GENERALIZED SITE-SCALE IMPACT:  ", te[i,2]))},
             if(te[i,1] == "SA.imp") {labs(title = paste0("SITE-SPECIFIC SCALE IMPACT:  ", te[i,2]))}
        )}
    )
    print(p)
  }
}
