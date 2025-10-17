#' @title Estimates the magnitude of an ecological impact
#'
#' @description
#' It provides an estimate of the magnitude of an ecological impact on the mean and variance of the data
#' that can be detected, given the observed P-Value from the asymmetrical analysis of variance.
#'
#' @param data a data frame with temporal factors, spatial factors and variable(s).
#' @param aav.design \code{baci} for beyond Before/After-Control/Impact design or \code{aci} for Post-Impact design.
#' @param n.ftemp number of temporal factors.
#' @param n.fspac number of spatial factors.
#' @param names.impact name of the asymmetrical level of the Location factor (e.g. 'impact').
#' @param names.before only for \code{baci} design. Name of the before level of the Before/After factor (e.g. 'before').
#' @param p.aav level of significance for the asymmetrical ANOVA. Default is \code{0.05}.
#' @param pooling remove a term from the analysis using a post-hoc pooling procedure. Default is \code{TRUE}.
#' @param method only when \code{baci} design and \code{pooling = TRUE}.
#'     Method to be used for the pooling procedure.
#'     \code{pool} to pool a term or a set of terms with zero contribution to the model or
#'     \code{eliminate} to completely exclude a term from the model. Default is \code{pool}.
#' @param p.pooling only when \code{pooling = TRUE}.
#'     Level of significance for the pooling procedure. Default is \code{0.25}.
#' @param aci.tailed.test only for \code{aci} post-impact design with only one disturbed location
#'     and multiple control locations. Perform a 1-tailed or a 2-tailed test for all terms.
#'     Default is \code{1}.
#' @param magnitude_md a vector to estimate alterations to the mean values of the data.
#' @param magnitude_var a vector to estimate alterations to the variance values of the data.
#'
#' @details
#' The magnitude of an ecological impact is estimated using the same data types, designs, factors
#' and parameters as in \code{\link[AAV]{aav}}, but only considers one variable at a time.
#'
#' The estimates are made using either the original data from the control group
#' (post-impact designs) or data from the control group before the impact (BACI design).
#' These estimated values are then applied to the level of the relevant factors in the
#' impact group (for post-impact designs) or to the impact group before the impact for BACI designs.
#'
#' The relevant factors depend on the number of temporal and spatial factors analysed.
#' To estimate long-term or local-scale impacts, we manipulate all the data involved.
#' To estimate medium-term, short-term or small-scale impacts, we modify just one
#' relevant factor (one period, one time, one site, or one combination of these).
#'
#' The magnitude of the ecological impact on the mean and variance of the data
#' is simulated using a linear function:
#' Y = a + bX,
#' where 'a' is the y-intercept (the additive constant, shift or constant), and 'b'
#' is the slope (the multiplicative constant, scale factor or coefficient).
#' The logic behind the calculation is that changing the mean involves modifying
#' the intercept, while changing the variance involves modifying the slope.
#' The random data is a resample from the original control or control-before data.
#'
#' The terms can be interpreted separated each other. The relevant terms are:
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
#' The function returns a list containing the following items:
#' \item{data_magimp}{the original dataset, along with a random data
#' and the simulated  magnitude of the impact on the mean and variance of the data.}
#' \item{aav_test}{the result of the complete asymmetrical ANOVA tables for the original
#' dataset and for each simulation of the magnitude of the impact on the mean and variance of the data.}
#' \item{extract_aav}{a table with the characterization of the design and factors used
#' ('aav.design','n.ftemp','n.fspac'),
#' a code for the type of the impact ('impType'), the magnitude of the variance alteration ('magVAR'),
#' the magnitude of the mean alteration ('magMD'), the corresponding term of the analysis ('term'),
#' and the P-Value from the specific analysis ('p.value').}
#'
#' @author Paulo Pagliosa \email{paulo.pagliosa@ufsc.br}
#'
#' @references
#' Underwood, A.J. 1997. \emph{Experiments in Ecology: Their Logical Design and Interpretation Using
#' Analysis of Variance}. Cambridge University Press, Cambridge. \doi{10.1017/CBO9780511806407}.
#'
#' @seealso
#' \code{\link[AAV]{aav}},
#' \code{\link[GAD]{ggmagnitude}}
#'
#' @examples
#' # Example for BACI design containing 3 temporal factors and 2 spatial factor:
#' data(baci.3t2s)
#' ex.baci.3t2s<-magnitude(data = baci.3t2s, aav.design = "baci", n.ftemp = 3, n.fspac = 2,
#'     names.impact = "impact", names.before = "before",
#'     pooling = TRUE, method = "eliminate", p.aav = 0.05, p.pooling = 0.05,
#'     magnitude_md = c(0.50, 1),
#'     magnitude_var = c(0.1, 1))
#' ex.baci.3t2s
#'
#' @include aav.R
#'
#' @export

magnitude <- function(data, aav.design, n.ftemp, n.fspac, names.impact, names.before = NULL,
                 p.aav = 0.05, pooling = TRUE, method = "pool", p.pooling = 0.25, aci.tailed.test = 1,
                 magnitude_md = c(0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1),
                 magnitude_var = c(0.1,1,2,5,10,25,50,75,100)) {

  if(aav.design %in% c("baci","aci") != T) stop("The 'aav.design' object does not exist. Check your command line.")
  if(names.impact %in% unique(data$location) != T) stop("The 'names.impact' object does not exist. Check your command line.")
  if(ncol(data) != (n.ftemp + n.fspac + 1)) message("The data has more than one variable. Using the first one.")
  options(dplyr.summarise.inform = FALSE) # To suppress dplyr "regrouping output by" message from summarise()

  ###--------------------------------------------------------- Impact Simulation

    data<-data[,c(1:(n.ftemp+n.fspac+1))] # ensure use for only one variable

    if(aav.design == "baci") {
      # Stops
      if(names.before %in% unique(data$before_after) != T) stop("The 'names.before' object does not exist. Check your command line")
      if(n.ftemp <= 1 | n.ftemp >= 4 | n.fspac == 0 | n.fspac >= 3) stop("Select different values for 'n.ftemp' or 'n.fspac'.")

      # Factor names
      if(n.ftemp == 3 & n.fspac == 2) {
        names(data)<-c("before_after", "period", "time", "location", "site", "obs_obs_obs")
        period1<-unique(data$period)[1]
        time1<-unique(data$time)[1]
        site1<-unique(data$site)[1]
      }
      if(n.ftemp == 3 & n.fspac == 1) {
        names(data)<-c("before_after", "period", "time", "location", "obs_obs_obs")
        period1<-unique(data$period)[1]
        time1<-unique(data$time)[1]
      }
      if(n.ftemp == 2 & n.fspac == 2) {
        names(data)<-c("before_after", "time", "location", "site", "obs_obs_obs")
        time1<-unique(data$time)[1]
        site1<-unique(data$site)[1]
      }
      if(n.ftemp == 2 & n.fspac == 1) {
        names(data)<-c("before_after", "time", "location", "obs_obs_obs")
        time1<-unique(data$time)[1]
      }

      ### THE RANDOM CASES
      ### LONG-TERM LOCAL-SCALE IMPACT (press impact)
      data$L.imp_0_0<-data$obs_obs_obs
      size<-length(data$obs_obs_obs[data$before_after != names.before & data$location == names.impact])
      data$L.imp_0_0[data$before_after != names.before & data$location == names.impact]<- sample(
        data$obs_obs_obs[data$before_after == names.before & data$location == names.impact], size = size, replace = FALSE) # resample obs_obs_obs data

      if(n.fspac == 2) {
        ### LONG-TERM SMALL-SCALE IMPACT (press impact)
        data$S.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$before_after != names.before & data$location == names.impact & data$site == site1])
        data$S.imp_0_0[data$before_after != names.before & data$location == names.impact & data$site == site1]<- sample(
          data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1], size = size, replace = FALSE) # resample obs_obs_obs data
      }
      if(n.ftemp == 3) {
        ### MEDIUM-TERM LOCAL-SCALE (pulse impact)
        data$P.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$before_after != names.before & data$location == names.impact & data$period == period1])
        data$P.imp_0_0[data$before_after != names.before & data$location == names.impact & data$period == period1]<- sample(
          data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$period == period1], size = size, replace = FALSE) # resample obs_obs_obs data

        ### SHORT-TERM LOCAL-SCALE IMPACT (pulse impact)
        data$TP.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$before_after != names.before & data$location == names.impact & data$time == time1 & data$period == period1])
        data$TP.imp_0_0[data$before_after != names.before & data$location == names.impact & data$time == time1 & data$period == period1]<- sample(
          data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$time == time1 & data$period == period1], size = size, replace = FALSE) # resample obs_obs_obs data
      }
      if(n.ftemp == 3 & n.fspac == 2) {
        ### MEDIUM-TERM SMALL-SCALE (pulse impact)
        data$SP.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$before_after != names.before & data$location == names.impact & data$site == site1 & data$period == period1])
        data$SP.imp_0_0[data$before_after != names.before & data$location == names.impact & data$site == site1 & data$period == period1]<- sample(
          data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1 & data$period == period1], size = size, replace = FALSE) # resample obs_obs_obs data

        ### SHORT-TERM SMALL-SCALE IMPACT (pulse impact)
        data$STP.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$before_after != names.before & data$location == names.impact & data$site == site1 & data$time == time1 & data$period == period1])
        data$STP.imp_0_0[data$before_after != names.before & data$location == names.impact & data$site == site1 & data$time == time1 & data$period == period1]<- sample(
          data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1 & data$time == time1 & data$period == period1], size = size, replace = FALSE) # resample obs_obs_obs data
      }
      if(n.ftemp == 2) {
        ### SHORT-TERM LOCAL-SCALE IMPACT (pulse impact)
        data$T.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$before_after != names.before & data$location == names.impact & data$time == time1])
        data$T.imp_0_0[data$before_after != names.before & data$location == names.impact & data$time == time1]<- sample(
          data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$time == time1], size = size, replace = FALSE) # resample obs_obs_obs data
      }
      if(n.ftemp == 2 & n.fspac == 2) {
        ### SHORT-TERM SMALL-SCALE IMPACT (pulse impact)
        data$ST.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$before_after != names.before & data$location == names.impact & data$site == site1 & data$time == time1])
        data$ST.imp_0_0[data$before_after != names.before & data$location == names.impact & data$site == site1 & data$time == time1]<- sample(
          data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1 & data$time == time1], size = size, replace = FALSE) # resample obs_obs_obs data
      }

      # Apply the magnitudes of the impact to the real data
      # Apply a linear function: Y = a + bX
      # a = y-intercept (the additive constant, the shift; the constant)
      ## to change the mean, modify the intercept
      # b = slope (the multiplicative constant, the scale factor; the coefficient)
      ## to change variance, modify the slope
      # Define the new variable as a linear function of the old one (the random variable).
      # b_new <- sqrt(1 + magnitude of impact) # to  increase the variance
      # b_new <- sqrt(1 - magnitude of impact) # to decrease the variance
      # To keep the variance the same, you must set the b_new = 1

      # Magnitude over the data variation
      for (magVAR in magnitude_var) {
        b_new <- magVAR

        # Magnitude over the data mean
        for (magMD in magnitude_md) {

        ### LONG-TERM LOCAL-SCALE IMPACT (press impact)
        md_imp<-mean(data$obs_obs_obs[data$before_after == names.before & data$location == names.impact]) # md from impact location before
        a_new<-(md_imp * (1 + magMD)) - b_new * md_imp
        data$new<-data$obs_obs_obs
        data$new[data$before_after != names.before & data$location == names.impact]<- round(a_new + b_new * data$obs_obs_obs[data$before_after == names.before & data$location == names.impact])
        names(data)[length(data)]<-paste0("L.imp_",magVAR,"_",magMD)

        if(n.fspac == 2) {
          ### LONG-TERM SMALL-SCALE IMPACT (press impact)
          md_s1<-mean(data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1]) # md from impact s1 before
          a_new<-(md_s1 * (1 + magMD)) - b_new * md_s1
          data$new<-data$obs_obs_obs
          data$new[data$before_after != names.before & data$location == names.impact & data$site == site1]<- round(a_new + b_new * data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1])
          names(data)[length(data)]<-paste0("S.imp_",magVAR,"_",magMD)
        }
        if(n.ftemp == 3) {
          ### MEDIUM-TERM LOCAL-SCALE (pulse impact)
          md_p1<- mean(data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$period == period1]) # md from impact p1 before
          a_new<-(md_p1 * (1 + magMD)) - b_new * md_p1
          data$new<-data$obs_obs_obs
          data$new[data$before_after != names.before & data$location == names.impact & data$period == period1]<- round(a_new + b_new * data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$period == period1])
          names(data)[length(data)]<-paste0("P.imp_",magVAR,"_",magMD)

          ### SHORT-TERM LOCAL-SCALE IMPACT (pulse impact)
          md_t1p1<- mean(data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$time == time1 & data$period == period1]) # md from impact t1p1 before
          a_new<-(md_t1p1 * (1 + magMD)) - b_new * md_t1p1
          data$new<-data$obs_obs_obs
          data$new[data$before_after != names.before & data$location == names.impact & data$time == time1 & data$period == period1]<- round(a_new + b_new * data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$time == time1 & data$period == period1])
          names(data)[length(data)]<-paste0("TP.imp_",magVAR,"_",magMD)
        }
        if(n.ftemp == 3 & n.fspac == 2) {
          ### MEDIUM-TERM SMALL-SCALE (pulse impact)
          md_s1p1<-mean(data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1 & data$period == period1]) # md from impact s1p1 before
          a_new<-(md_s1p1 * (1 + magMD)) - b_new * md_s1p1
          data$new<-data$obs_obs_obs
          data$new[data$before_after != names.before & data$location == names.impact & data$site == site1 & data$period == period1]<- round(a_new + b_new * data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1 & data$period == period1])
          names(data)[length(data)]<-paste0("SP.imp_",magVAR,"_",magMD)

          ### SHORT-TERM SMALL-SCALE IMPACT (pulse impact)
          md_s1t1p1<-mean(data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1 & data$time == time1 & data$period == period1]) # md from impact s1t1p1 before
          a_new<-(md_s1t1p1 * (1 + magMD)) - b_new * md_s1t1p1
          data$new<-data$obs_obs_obs
          data$new[data$before_after != names.before & data$location == names.impact & data$site == site1 & data$time == time1 & data$period == period1]<- round(a_new + b_new * data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1 & data$time == time1 & data$period == period1])
          names(data)[length(data)]<-paste0("STP.imp_",magVAR,"_",magMD)
        }
        if(n.ftemp == 2) {
          ### SHORT-TERM LOCAL-SCALE IMPACT (pulse impact)
          md_t1<- mean(data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$time == time1]) # md from t1 impact location before
          a_new<-(md_t1 * (1 + magMD)) - b_new * md_t1
          data$new<-data$obs_obs_obs
          data$new[data$before_after != names.before & data$location == names.impact & data$time == time1]<- round(a_new + b_new * data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$time == time1])
          names(data)[length(data)]<-paste0("T.imp_",magVAR,"_",magMD)
        }
        if(n.ftemp == 2 & n.fspac == 2) {
          ### SHORT-TERM SMALL-SCALE IMPACT (pulse impact)
          md_s1t1<-mean(data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1 & data$time == time1])
          a_new<-(md_s1t1 * (1 + magMD)) - b_new * md_s1t1
          data$new<-data$obs_obs_obs
          data$new[data$before_after != names.before & data$location == names.impact & data$site == site1 & data$time == time1]<- round(a_new + b_new * data$obs_obs_obs[data$before_after == names.before & data$location == names.impact & data$site == site1 & data$time == time1])
          names(data)[length(data)]<-paste0("ST.imp_",magVAR,"_",magMD)
        }
      }
      }
    }
    if(aav.design == "aci" & n.ftemp != 0) {
      # Stops
      if(n.ftemp >= 3 | n.fspac == 0 | n.fspac >= 3) stop("Select different values for 'n.ftemp' or 'n.fspac'.")

      # Factors names
      if(n.ftemp == 2 & n.fspac == 2) {
        names(data)<-c("period", "time", "location", "site", "obs_obs_obs")
        period1<-unique(data$period)[1]
        time1<-unique(data$time)[1]
        site1<-unique(data$site)[1]
        ulmc<-unique(data$location[data$location != names.impact])
        up<-unique(data$period)
        ut<-unique(data$time)
        us<-unique(data$site)
      }
      if(n.ftemp == 2 & n.fspac == 1) {
        names(data)<-c("period", "time", "location", "obs_obs_obs")
        period1<-unique(data$period)[1]
        time1<-unique(data$time)[1]
        ulmc<-unique(data$location[data$location != names.impact])
        up<-unique(data$period)
        ut<-unique(data$time)
      }
      if(n.ftemp == 1 & n.fspac == 2) {
        names(data)<-c("time", "location", "site", "obs_obs_obs")
        time1<-unique(data$time)[1]
        site1<-unique(data$site)[1]
        ulmc<-unique(data$location[data$location != names.impact])
        ut<-unique(data$time)
        us<-unique(data$site)
      }
      if(n.ftemp == 1 & n.fspac == 1) {
        names(data)<-c("time", "location", "obs_obs_obs")
        time1<-unique(data$time)[1]
        ulmc<-unique(data$location[data$location != names.impact])
        ut<-unique(data$time)
      }

      ### LONG-TERM LOCAL-SCALE IMPACT (press impact)
      data$L.imp_0_0<-data$obs_obs_obs
      size<-length(data$obs_obs_obs[data$location == names.impact])
      data$L.imp_0_0[data$location == names.impact]<- sample(
        data$obs_obs_obs[data$location != names.impact], size = size, replace = FALSE) # resample obs_obs_obs data

      if(n.fspac == 2) {
        ### LONG-TERM SMALL-SCALE IMPACT (press impact)
        data$S.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$location == names.impact & data$site == site1])
        data$S.imp_0_0[data$location == names.impact & data$site == site1]<- sample(
          data$obs_obs_obs[data$location != names.impact & data$site == site1], size = size, replace = FALSE) # resample obs_obs_obs data
      }
      if(n.ftemp == 2) {
        ### MEDIUM-TERM LOCAL-SCALE (pulse impact)
        data$P.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$location == names.impact & data$period == period1])
        data$P.imp_0_0[data$location == names.impact & data$period == period1]<- sample(
          data$obs_obs_obs[data$location != names.impact & data$period == period1], size = size, replace = FALSE) # resample obs_obs_obs data

        ### SHORT-TERM LOCAL-SCALE IMPACT (pulse impact)
        data$TP.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$location == names.impact & data$time == time1 & data$period == period1])
        data$TP.imp_0_0[data$location == names.impact & data$time == time1 & data$period == period1]<- sample(
          data$obs_obs_obs[data$location != names.impact & data$time == time1 & data$period == period1], size = size, replace = FALSE) # resample obs_obs_obs data
      }
      if(n.ftemp == 2 & n.fspac == 2) {
        ### MEDIUM-TERM SMALL-SCALE (pulse impact)
        data$SP.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$location == names.impact & data$site == site1 & data$period == period1])
        data$SP.imp_0_0[data$location == names.impact & data$site == site1 & data$period == period1]<- sample(
          data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$period == period1], size = size, replace = FALSE) # resample obs_obs_obs data

        ### SHORT-TERM SMALL-SCALE IMPACT (pulse impact)
        data$STP.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$location == names.impact & data$site == site1 & data$time == time1 & data$period == period1])
        data$STP.imp_0_0[data$location == names.impact & data$site == site1 & data$time == time1 & data$period == period1]<- sample(
          data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$time == time1 & data$period == period1], size = size, replace = FALSE) # resample obs_obs_obs data
      }
      if(n.ftemp == 1) {
        ### SHORT-TERM LOCAL-SCALE IMPACT (pulse impact)
        data$T.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$location == names.impact & data$time == time1])
        data$T.imp_0_0[data$location == names.impact & data$time == time1]<- sample(
          data$obs_obs_obs[data$location != names.impact & data$time == time1], size = size, replace = FALSE) # resample obs_obs_obs data
      }
      if(n.ftemp == 1 & n.fspac == 2) {
        ### SHORT-TERM SMALL-SCALE IMPACT (pulse impact)
        data$ST.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$location == names.impact & data$site == site1 & data$time == time1])
        data$ST.imp_0_0[data$location == names.impact & data$site == site1 & data$time == time1]<- sample(
          data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$time == time1], size = size, replace = FALSE) # resample obs_obs_obs data
      }

      for (magVAR in magnitude_var) {
        b_new <- magVAR
        for (magMD in magnitude_md) {

        ### LONG-TERM LOCAL-SCALE IMPACT (press impact)
        md_control<-mean(data$obs_obs_obs[data$location != names.impact]) # md from control locations
        a_new<-(md_control * (1 + magMD)) - b_new * md_control

        ### Criar novo conjunto de dados que preserva média e variância dos conjuntos originais
        l.out<-length(data$obs_obs_obs[data$location != names.impact])/length(ulmc)
        base <- seq(-l.out/2, l.out/2, length.out = l.out)
        tab_stat<-data.frame()
        for (xx in ulmc) {
          y<- data$obs_obs_obs[data$location == xx]
          this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
          tab_stat<-rbind(tab_stat,this.ts)
        }
        values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
        data$new<-data$obs_obs_obs
        data$new[data$location == names.impact]<- round(a_new + b_new * values_adj)
        names(data)[length(data)]<-paste0("L.imp_",magVAR,"_",magMD)

        if(n.fspac == 2) {
          ### LONG-TERM SMALL-SCALE IMPACT (press impact)
          md_s1<-mean(data$obs_obs_obs[data$location != names.impact & data$site == site1]) # md from controls s1
          a_new<-(md_s1 * (1 + magMD)) - b_new * md_s1
          l.out<-length(data$obs_obs_obs[data$location != names.impact & data$site == site1])/length(ulmc)
          base <- seq(-l.out/2, l.out/2, length.out = l.out)
          tab_stat<-data.frame()
          for (xx in ulmc) {
            y<- data$obs_obs_obs[data$location == xx & data$site == site1]
            this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
            tab_stat<-rbind(tab_stat,this.ts)
          }
          values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
          data$new<-data$obs_obs_obs
          data$new[data$location == names.impact & data$site == site1]<- round(a_new + b_new * values_adj)
          names(data)[length(data)]<-paste0("S.imp_",magVAR,"_",magMD)
        }
        if(n.ftemp == 2) {
          ### MEDIUM-TERM LOCAL-SCALE (pulse impact)
          md_p1<- mean(data$obs_obs_obs[data$location != names.impact & data$period == period1]) # md from controls p1
          a_new<-(md_p1 * (1 + magMD)) - b_new * md_p1
          l.out<-length(data$obs_obs_obs[data$location != names.impact & data$period == period1])/length(ulmc)
          base <- seq(-l.out/2, l.out/2, length.out = l.out)
          tab_stat<-data.frame()
          for (xx in ulmc) {
            y<- data$obs_obs_obs[data$location == xx & data$period == period1]
            this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
            tab_stat<-rbind(tab_stat,this.ts)
          }
          values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
          data$new<-data$obs_obs_obs
          data$new[data$location == names.impact & data$period == period1]<- round(a_new + b_new * values_adj)
          names(data)[length(data)]<-paste0("P.imp_",magVAR,"_",magMD)

          ### SHORT-TERM LOCAL-SCALE IMPACT (pulse impact)
          md_t1p1<- mean(data$obs_obs_obs[data$location != names.impact & data$time == time1 & data$period == period1]) # md from controls t1p1
          a_new<-(md_t1p1 * (1 + magMD)) - b_new * md_t1p1
          l.out<-length(data$obs_obs_obs[data$location != names.impact & data$time == time1 & data$period == period1])/length(ulmc)
          base <- seq(-l.out/2, l.out/2, length.out = l.out)
          tab_stat<-data.frame()
          for (xx in ulmc) {
            y<- data$obs_obs_obs[data$location == xx & data$time == time1 & data$period == period1]
            this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
            tab_stat<-rbind(tab_stat,this.ts)
          }
          values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
          data$new<-data$obs_obs_obs
          data$new[data$location == names.impact & data$time == time1 & data$period == period1]<- round(a_new + b_new * values_adj)
          names(data)[length(data)]<-paste0("TP.imp_",magVAR,"_",magMD)
        }
        if(n.ftemp == 2 & n.fspac == 2) {
          ### MEDIUM-TERM SMALL-SCALE (pulse impact)
          md_s1p1<-mean(data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$period == period1]) # md from controls s1p1
          a_new<-(md_s1p1 * (1 + magMD)) - b_new * md_s1p1
          l.out<-length(data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$period == period1])/length(ulmc)
          base <- seq(-l.out/2, l.out/2, length.out = l.out)
          tab_stat<-data.frame()
          for (xx in ulmc) {
            y<- data$obs_obs_obs[data$location == xx & data$site == site1 & data$period == period1]
            this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
            tab_stat<-rbind(tab_stat,this.ts)
          }
          values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
          data$new<-data$obs_obs_obs
          data$new[data$location == names.impact & data$site == site1 & data$period == period1]<- round(a_new + b_new * values_adj)
          names(data)[length(data)]<-paste0("SP.imp_",magVAR,"_",magMD)

          ### SHORT-TERM SMALL-SCALE IMPACT (pulse impact)
          md_s1t1p1<-mean(data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$time == time1 & data$period == period1]) # md from controls s1t1p1
          a_new<-(md_s1t1p1 * (1 + magMD)) - b_new * md_s1t1p1
          l.out<-length(data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$time == time1 & data$period == period1])/length(ulmc)
          base <- seq(-l.out/2, l.out/2, length.out = l.out)
          tab_stat<-data.frame()
          for (xx in ulmc) {
            y<- data$obs_obs_obs[data$location == xx & data$site == site1 & data$time == time1 & data$period == period1]
            this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
            tab_stat<-rbind(tab_stat,this.ts)
          }
          values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
          data$new<-data$obs_obs_obs
          data$new[data$location == names.impact & data$site == site1 & data$time == time1 & data$period == period1]<- round(a_new + b_new * values_adj)
          names(data)[length(data)]<-paste0("STP.imp_",magVAR,"_",magMD)
        }
        if(n.ftemp == 1) {
          ### SHORT-TERM LOCAL-SCALE IMPACT (pulse impact)
          md_t1<- mean(data$obs_obs_obs[data$location != names.impact & data$time == time1]) # md from controls t1
          a_new<-(md_t1 * (1 + magMD)) - b_new * md_t1
          l.out<-length(data$obs_obs_obs[data$location != names.impact & data$time == time1])/length(ulmc)
          base <- seq(-l.out/2, l.out/2, length.out = l.out)
          tab_stat<-data.frame()
          for (xx in ulmc) {
            y<- data$obs_obs_obs[data$location == xx & data$time == time1]
            this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
            tab_stat<-rbind(tab_stat,this.ts)
          }
          values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
          data$new<-data$obs_obs_obs
          data$new[data$location == names.impact & data$time == time1]<- round(a_new + b_new * values_adj)
          names(data)[length(data)]<-paste0("T.imp_",magVAR,"_",magMD)
        }
        if(n.ftemp == 1 & n.fspac == 2) {
          ### SHORT-TERM SMALL-SCALE IMPACT (pulse impact)
          md_s1t1<-mean(data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$time == time1]) # md from controls s1t1
          a_new<-(md_s1t1 * (1 + magMD)) - b_new * md_s1t1
          l.out<-length(data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$time == time1])/length(ulmc)
          base <- seq(-l.out/2, l.out/2, length.out = l.out)
          tab_stat<-data.frame()
          for (xx in ulmc) {
            y<- data$obs_obs_obs[data$location == xx & data$site == site1 & data$time == time1]
            this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
            tab_stat<-rbind(tab_stat,this.ts)
          }
          values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
          data$new<-data$obs_obs_obs
          data$new[data$location == names.impact & data$site == site1 & data$time == time1]<- round(a_new + b_new * values_adj)
          names(data)[length(data)]<-paste0("ST.imp_",magVAR,"_",magMD)
        }
      }
      }
    }
    if(aav.design == "aci" & n.ftemp == 0) {
      # Stops
      if(n.fspac <= 1 | n.fspac >= 4) stop("Select different values for 'n.ftemp' or 'n.fspac'.")

      # Factors names
      if(n.fspac == 2) {
        names(data)<-c("place", "location", "obs_obs_obs")
        place1<-unique(data$place)[1]
        ulmc<-unique(data$location[data$location != names.impact])
        up<-unique(data$place)
      }
      if(n.fspac == 3) {
        names(data)<-c("place", "location", "site", "obs_obs_obs")
        place1<-unique(data$place)[1]
        site1<-unique(data$site)[1]
        ulmc<-unique(data$location[data$location != names.impact])
        up<-unique(data$place)
        us<-unique(data$site)
      }

      ### GENERALIZED LOCAL-SCALE IMPACT
      data$L.imp_0_0<-data$obs_obs_obs
      size<-length(data$obs_obs_obs[data$location == names.impact])
      data$L.imp_0_0[data$location == names.impact]<- sample(
        data$obs_obs_obs[data$location != names.impact], size = size, replace = FALSE) # resample obs_obs_obs data

      ### LOCAL-SPECIFIC SCALE IMPACT
      data$LA.imp_0_0<-data$obs_obs_obs
      size<-length(data$obs_obs_obs[data$location == names.impact & data$place == place1])
      data$LA.imp_0_0[data$location == names.impact & data$place == place1]<- sample(
        data$obs_obs_obs[data$location != names.impact & data$place == place1], size = size, replace = FALSE) # resample obs_obs_obs data

      if(n.fspac == 3) {
        ### GENERALIZED SITE-SCALE IMPACT
        data$S.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$location == names.impact & data$site == site1])
        data$S.imp_0_0[data$location == names.impact & data$site == site1]<- sample(
          data$obs_obs_obs[data$location != names.impact & data$site == site1], size = size, replace = FALSE) # resample obs_obs_obs data

        ### SITE-SPECIFIC SCALE
        data$SA.imp_0_0<-data$obs_obs_obs
        size<-length(data$obs_obs_obs[data$location == names.impact & data$site == site1 & data$place == place1])
        data$SA.imp_0_0[data$location == names.impact & data$site == site1 & data$place == place1]<- sample(
          data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$place == place1], size = size, replace = FALSE) # resample obs_obs_obs data
      }

      for (magVAR in magnitude_var) {
        b_new <- magVAR
        for (magMD in magnitude_md) {

        ### GENERALIZED LOCAL-SCALE IMPACT
        md_control<-mean(data$obs_obs_obs[data$location != names.impact]) # md from control locations
        a_new<-(md_control * (1 + magMD)) - b_new * md_control
        l.out<-length(data$obs_obs_obs[data$location != names.impact])/length(ulmc)
        base <- seq(-l.out/2, l.out/2, length.out = l.out)
        tab_stat<-data.frame()
        for (xx in ulmc) {
          y<- data$obs_obs_obs[data$location == xx]
          this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
          tab_stat<-rbind(tab_stat,this.ts)
        }
        values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
        data$new<-data$obs_obs_obs
        data$new[data$location == names.impact]<- round(a_new + b_new * values_adj)
        names(data)[length(data)]<-paste0("L.imp_",magVAR,"_",magMD)

        ### LOCAL-SPECIFIC SCALE IMPACT
        md_la1<- mean(data$obs_obs_obs[data$location != names.impact & data$place == place1]) # md from controls a1
        a_new<-(md_la1 * (1 + magMD)) - b_new * md_la1
        l.out<-length(data$obs_obs_obs[data$location != names.impact & data$place == place1])/length(ulmc)
        base <- seq(-l.out/2, l.out/2, length.out = l.out)
        tab_stat<-data.frame()
        for (xx in ulmc) {
          y<- data$obs_obs_obs[data$location == xx & data$place == place1]
          this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
          tab_stat<-rbind(tab_stat,this.ts)
        }
        values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
        data$new<-data$obs_obs_obs
        data$new[data$location == names.impact & data$place == place1]<- round(a_new + b_new * values_adj)
        names(data)[length(data)]<-paste0("LA.imp_",magVAR,"_",magMD)

        if(n.fspac == 3) {
          ### GENERALIZED SITE-SCALE IMPACT
          md_s1<-mean(data$obs_obs_obs[data$location != names.impact & data$site == site1]) # md from controls s1
          a_new<-(md_s1 * (1 + magMD)) - b_new * md_s1
          l.out<-length(data$obs_obs_obs[data$location != names.impact & data$site == site1])/length(ulmc)
          base <- seq(-l.out/2, l.out/2, length.out = l.out)
          tab_stat<-data.frame()
          for (xx in ulmc) {
            y<- data$obs_obs_obs[data$location == xx & data$site == site1]
            this.ts<-data.frame(term = xx,media = mean(y),desvio_padrao = sd(y))
            tab_stat<-rbind(tab_stat,this.ts)
          }
          values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
          data$new<-data$obs_obs_obs
          data$new[data$location == names.impact & data$site == site1]<- round(a_new + b_new * values_adj)
          names(data)[length(data)]<-paste0("S.imp_",magVAR,"_",magMD)

          ### SITE-SPECIFIC SCALE IMPACT
          md_s1a1<-mean(data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$place == place1]) # md from controls s1a1
          a_new<-(md_s1a1 * (1 + magMD)) - b_new * md_s1a1
          l.out<-length(data$obs_obs_obs[data$location != names.impact & data$site == site1 & data$place == place1])/length(ulmc)
          base <- seq(-l.out/2, l.out/2, length.out = l.out)
          tab_stat<-data.frame()
          for (xx in ulmc) {
            y<- data$obs_obs_obs[data$location == xx & data$site == site1 & data$place == place1]
            this.ts<-data.frame(term = xx, media = mean(y), desvio_padrao = sd(y))
            tab_stat<-rbind(tab_stat,this.ts)
          }
          values_adj <- mean(tab_stat$media) + base * mean(tab_stat$desvio_padrao) / sd(base)
          data$new<-data$obs_obs_obs
          data$new[data$location == names.impact & data$site == site1 & data$place == place1]<- round(a_new + b_new * values_adj)
          names(data)[length(data)]<-paste0("SA.imp_",magVAR,"_",magMD)
        }
      }
      }
    }

    ###--------------------------------------------------------------- AAV_tests

    aav_tests<-aav(data = data, aav.design = aav.design, n.ftemp = n.ftemp, n.fspac = n.fspac,
                    names.impact = names.impact, names.before = names.before,
                    p.aav = p.aav, pooling = pooling, method = method, p.pooling = p.pooling,
                    aci.tailed.test = aci.tailed.test)


    ###----------------------------------- extract p-values for each impact type

    # real data (obs_obs_obs)
    tab.final<-data.frame(matrix(data=NA, nrow=0, ncol=8))
    names(tab.final)<-c("aav.design","n.ftemp","n.fspac","impType","magVAR","magMD","term","p.value")

  for (fc in 1:length(aav_tests$asym.anova_complete)) {
    tab<-data.frame(matrix(data=NA, nrow=8, ncol=8))
    names(tab)<-c("aav.design","n.ftemp","n.fspac","impType","magVAR","magMD","term","p.value")
    tab$aav.design<-aav.design
    tab$n.ftemp<-n.ftemp
    tab$n.fspac<-n.fspac
    fc.n<-names(aav_tests$asym.anova_complete[fc])
    tab$impType<-strsplit(fc.n,"_")[[1]][1]
    tab$magVAR<-strsplit(fc.n,"_")[[1]][2]
    tab$magMD<-strsplit(fc.n,"_")[[1]][3]
    aav.res<-aav_tests$asym.anova_complete[[fc]]

    if(aav.design == "baci") {
      if(tab$impType[1] == "L.imp" | tab$impType[1] == "obs") {
        tab$term[1]<- "B x I"
        #p.value
        tab$p.value[1]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "B x I" & aav.res$Final.Table == "Yes"]
      }

      if(n.fspac == 2) {
        if(tab$impType[2] == "S.imp" | tab$impType[1] == "obs") {
          tab$term[2]<- "B x S(I)"
          #p.value
          tab$p.value[2]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "B x S(I)" & aav.res$Final.Table == "Yes"]
        }
      }
      if(n.ftemp == 3) {
        if(tab$impType[3] == "P.imp" | tab$impType[1] == "obs") {
          tab$term[3]<- "P(Aft) x I"
          #p.value
          tab$p.value[3]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "P(Aft) x I"  & aav.res$Final.Table == "Yes"]
        }
        if(tab$impType[4] == "TP.imp" | tab$impType[1] == "obs") {
          tab$term[4]<- "T(P(Aft)) x I"
          #p.value
          tab$p.value[4]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "T(P(Aft)) x I"  & aav.res$Final.Table == "Yes"]
        }
      }
      if(n.ftemp == 3 & n.fspac == 2) {
        if(tab$impType[5] == "SP.imp" | tab$impType[1] == "obs") {
          tab$term[5]<- "P(Aft) x S(I)"
          #p.value
          tab$p.value[5]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "P(Aft) x S(I)"  & aav.res$Final.Table == "Yes"]
        }
        if(tab$impType[6] == "STP.imp" | tab$impType[1] == "obs") {
          tab$term[6]<- "T(P(Aft)) x S(I)"
          #p.value
          tab$p.value[6]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "T(P(Aft)) x S(I)"  & aav.res$Final.Table == "Yes"]
        }
      }
      if(n.ftemp == 2) {
        if(tab$impType[7] == "T.imp" | tab$impType[1] == "obs") {
          tab$term[7]<- "T(Aft) x I"
          #p.value
          tab$p.value[7]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "T(Aft) x I"  & aav.res$Final.Table == "Yes"]
        }
      }
      if(n.ftemp == 2 & n.fspac == 2) {
        if(tab$impType[8] == "ST.imp" | tab$impType[1] == "obs") {
          tab$term[8]<- "T(Aft) x S(I)"
          #p.value
          tab$p.value[8]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "T(Aft) x S(I)"  & aav.res$Final.Table == "Yes"]
        }
      }
    }
    if(aav.design == "aci" & n.ftemp != 0) {
      if(tab$impType[1] == "L.imp" | tab$impType[1] == "obs") {
        tab$term[1]<- "I"
        #p.value
        tab$p.value[1]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "Impact vs Controls = I" & aav.res$Final.Table == "Yes"]
      }
      if(n.fspac == 2) {
        if(tab$impType[2] == "S.imp" | tab$impType[1] == "obs") {
          tab$term[2]<- "S(I)"
          #p.value
          tab$p.value[2]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "S(I)" & aav.res$Final.Table == "Yes"]
        }
      }
      if(n.ftemp == 2) {
        if(tab$impType[3] == "P.imp" | tab$impType[1] == "obs") {
          tab$term[3]<- "P x I"
          #p.value
          tab$p.value[3]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "P x I"  & aav.res$Final.Table == "Yes"]
        }
        if(tab$impType[4] == "TP.imp" | tab$impType[1] == "obs") {
          tab$term[4]<- "T(P) x I"
          #p.value
          tab$p.value[4]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "T(P) x I"  & aav.res$Final.Table == "Yes"]
        }
      }
      if(n.ftemp == 2 & n.fspac == 2) {
        if(tab$impType[5] == "SP.imp" | tab$impType[1] == "obs") {
          tab$term[5]<- "P x S(I)"
          #p.value
          tab$p.value[5]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "P x S(I)"  & aav.res$Final.Table == "Yes"]
        }
        if(tab$impType[6] == "STP.imp" | tab$impType[1] == "obs") {
          tab$term[6]<- "T(P) x S(I)"
          #p.value
          tab$p.value[6]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "T(P) x S(I)"  & aav.res$Final.Table == "Yes"]
        }
      }
      if(n.ftemp == 1) {
        if(tab$impType[7] == "T.imp" | tab$impType[1] == "obs") {
          tab$term[7]<- "T x I"
          #p.value
          tab$p.value[7]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "T x I"  & aav.res$Final.Table == "Yes"]
        }
      }
      if(n.ftemp == 1 & n.fspac == 2) {
        if(tab$impType[8] == "ST.imp" | tab$impType[1] == "obs") {
          tab$term[8]<- "T x S(I)"
          #p.value
          tab$p.value[8]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "T x S(I)"  & aav.res$Final.Table == "Yes"]
        }
      }
    }
    if(aav.design == "aci" & n.ftemp == 0) {
      if(tab$impType[1] == "L.imp" | tab$impType[1] == "obs") {
        tab$term[1]<- "T"
        #p.value
        tab$p.value[1]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "Treatment = T" & aav.res$Final.Table == "Yes"]
      }
      if(tab$impType[2] == "LA.imp" | tab$impType[1] == "obs") {
        tab$term[2]<- "P x T"
        #p.value
        tab$p.value[2]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "P x T"  & aav.res$Final.Table == "Yes"]
      }
      if(n.fspac == 3) {
        if(tab$impType[3] == "S.imp" | tab$impType[1] == "obs") {
          tab$term[3]<- "S(T)"
          #p.value
          tab$p.value[3]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "Sites(T) = S(T)" & aav.res$Final.Table == "Yes"]
        }
        if(tab$impType[4] == "SA.imp" | tab$impType[1] == "obs") {
          tab$term[4]<- "P x S(T)"
          #p.value
          tab$p.value[4]<-aav.res$P.Value.Upper.Tail[aav.res$Source.of.Variation == "P x S(T)"  & aav.res$Final.Table == "Yes"]
        }
      }
    }
    tab<-tab[!is.na(tab$term), ] # remove rows not used
    tab.final<-rbind(tab.final,tab)
  } # end extract

  magnitude.list<-list()
  magnitude.list[[1]]<-data
  names(magnitude.list)[[1]]<-"data_magimp"
  magnitude.list[[2]]<-aav_tests
  names(magnitude.list)[[2]]<-"aav_test"
  magnitude.list[[3]]<-tab.final
  names(magnitude.list)[[3]]<-"extract_aav"

  return(magnitude.list)
}
