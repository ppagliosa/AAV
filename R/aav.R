#' @title Asymmetrical ANOVA From Beyond BACI and Post-Impact Designs
#'
#' @description
#' Provides an asymmetrical ANOVA from a beyond BACI (before/after-control/impact)
#' design with multiple control locations, a post-impact design with only one
#' disturbed location and multiple control locations, or a post-impact design
#' involving more than one place with one disturbed location and multiple control locations.
#'
#' Asymmetrical ANOVAs are essentially multi-factorial.
#' The beyond BACI design can be used for models with 2 or 3 temporal factors
#' and 1 or 2 spatial factors. A post-impact design with only one disturbed location
#' can be used for models with 1 or 2 temporal factors
#' and 1 or 2 spatial factors. A post-impact design with more than one place with
#' one disturbed location can be used for models
#' with 2 or 3 spatial factors.
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
#'
#' @details
#' Asymmetrical analyses of variance are constructed by combining
#' the Sum of Squares values from separate (symmetrical) ANOVAs.
#' In beyond BACI designs, four fully orthogonal ANOVAs are required: one on all the data,
#' one on 'control' data, one on 'before' data, and one on 'control' and 'before' data.
#' In post-impact designs with only one disturbed location two fully orthogonal
#' ANOVAs are performed: one on all the data, and one on the 'control' data.
#' In post-impact designs involving more than one place with one disturbed location,
#' two fully orthogonal and one fully nested ANOVA need to be performed: one on all the data
#' (orthogonal) and two on the 'control' data (orthogonal and nested).
#' \code{aav} uses \code{\link[GAD]{GAD}} (Sandrini-Neto et al. 2024) to perform symmetrical ANOVAs.
#'
#' The \code{aav} analysis fits mixed asymmetrical ANOVA designs with a specific
#' combination of orthogonal/nested and fixed/random factors:
#' The  beyond BACI design can be used with five factors (3 temporal and 2 spatial),
#' four factors (3 temporal and 1 spatial or 2 temporal and 2 spatial), and
#' three factors (2 temporal and 1 spatial).
#' The post-impact design with only one disturbed location runs for
#' four factors (2 temporal and 2 spatial), three factors (2 temporal and 1 spatial or
#' 1 temporal and 2 spatial), and two factors (1 temporal and 1 spatial).
#' The post-impact design involving more than one place with one disturbed location runs for
#' three factors (0 temporal and 3 spatial) and two factors (0 temporal and 2 spatial).
#'
#' In all models, the Location factor is asymmetric (e.g., one 'impact' group and
#' two or more 'control' groups), meaning there are different numbers of levels of this
#' nested factor at each  level of its upper-level factor (e.g., Place factor) or
#' contrasting factors (e.g., Before-After, Period, Time).
#'
#' In all models the temporal and spatial factors are orthogonal.
#' Within each factor type (temporal or spatial) the factors are hierarchically nested.
#' In the beyond BACI and post-impact design with only one disturbed location models,
#' the first factor and the Location factor are fixed, while all others are random.
#' In the post-impact design involving than one place with one disturbed location,
#' the first factor is fixed and the others are random.
#'
#' The data frame must follow a specific structure:
#' 1- first the temporal factors, then the spatial factors, and then the variable(s);
#' 2- each factor type must appear in a hierarchical order.
#'
#' For beyond BACI design the order of the temporal factors is Before-after, Period
#' and Time (3 factors) or Before_after and Time (2 factors).
#' For spatial factors, the order is Location and Site (2 factors) or just Location (1 factor).
#' In a post-impact design involving only one disturbed location and multiple control locations,
#' the temporal factor order must be either Period and Time (2 factors) or just Time (1 factor), and
#' the spatial factor order must be either Location and Site (2 factors) or just Location (1 factor).
#' In a post-impact design involving more than one place with one disturbed location and multiple control locations,
#' the order of the spatial factors must be Place, Location and Site (3 factors) or Place and Location (2 factors).
#'
#' If the names of the factors do not match those listed above, \code{AAV} will change them.
#'
#' When the data frame contains more than one variable, separate univariate ANOVAs are run for each one.
#'
#' In the case of a pooling structure, both the \code{pool} and \code{eliminate} methods
#' are available for the beyond BACI design. For post-impact design when
#' \code{pooling = TRUE} the analysis is running using the method \code{pool}.
#'
#' For \code{pool}, the degrees of freedom (Df) and sums of squares (SS) for that term are
#' pooled with the term(s) that have an equivalent mean square expectation
#' after the component of variation of that term has been set to zero. This
#' component of variation is set to zero if it is >= \code{p.pooling} value.
#' When \code{eliminate}, the Df and SS for a term is added to the residual
#' Df and SS and the term is treated as ever having been part of the model.
#'
#' The magnitude of interactions in time and space can be analysed using a two-tailed tests.
#' For beyond BACI design the two-tailed tests are automatically performed to contrast
#' the differences in the space and time interactions from before to after the disturbance occurs.
#' This is crucial for establishing a cause-and-effect relationship in ecological impact detection.
#' for Post-Impact design involving only one impact location,
#' where temporal interactions among locations occur just after the disturbance,
#' a two-tailed test can be performed (by selecting \code{aci.tailed.test = 2}) to determine
#' whether and in which direction the time trends differ between impact and controls locations.
#' For a Post-Impact design involving more than one place with one disturbed location,
#' there is no need for a two-tailed test. In this design, F-tests contrast treatments
#' against an interaction term of spatial scale and treatments and are properly run as one-tailed tests.
#'
#' @returns
#' The function returns a list containing the following items:
#' \item{sym.anova}{the symmetrical ANOVA tables used for running asymmetrical ANOVAs.}
#' \item{asym.anova_final}{the final asymmetrical ANOVA tables.
#' The 'final table' is a typical ANOVA table
#' ('Source.of.Variation','Df', 'Sum.Sq','Mean.Sq','F.value','P.Value.Upper.Tail','P.Value.Lower.Tail')
#' with the partitioning of the variance components between control and impacts
#' groups and between before and after disturbance groups (only for BACI designs).
#' The 'final table' has two additional step-by-step elements:
#' a 'Next.Step' that indicates where to go from here, according to the
#' calculation of the mean square estimates for that term,
#' and an 'Interpretation' that states the results found.
#' There is no 'P.Value.Lower.Tail' for the a post-impact design with more than
#' one place with one disturbed location. For beyond BACI design 'P.Value.Lower.Tail'
#' appears only when there is a valid test for a two-tailed analysis.}
#' \item{asym.anova_complete}{the complete asymmetrical ANOVA tables.
#' The 'complete table' presents all elements used in the asymmetrical analysis,
#' with the 'final table' elements plus:
#' The complement for a two-tailed analysis ('P.Value.Lower.Tail');
#' the degrees of freedom formula ('Df.Formula');
#' from which symmetrical ANOVA derived the SS and Df values ('Calculated.From');
#' who is the denominator to the calculus of the F-ratio for that term ('Denominator');
#' the Df, SS and MS for the pooled term(s) ('Df.Pooled','Sum.Sq.Pooled','Mean.Sq.Pooled');
#' based on the \code{p.aav} and \code{p.pooling} values, a statement if that term is significant,
#' non-significant or can be eliminated ('Significant');
#' based on the \code{p.pooling} value and the mean square estimates, a statement on
#' which of the possible pooled term were used ('Post.Hoc.Pooling');
#' some designs has two-tailed analysis ('P.Value.Lower.Tail').}
#'
#' When analysing more than one variable at once, the results for each variable
#' appear in the sublists.
#'
#' @author Paulo Pagliosa \email{paulo.pagliosa@ufsc.br}
#'
#' @references
#' Glasby, T.M. 1997. Analysing data from post-impact studies using asymmetrical analyses
#' of variance: A case study of epibiota on marinas. Aust. J. Ecol., 22, 448-459.
#' \doi{10.1111/j.1442-9993.1997.tb00696.x}.
#'
#' Sandrini-Neto, L., Gilbert, E., Camargo, M. 2024. \emph{GAD: Analysis of Variance from General Principles}.
#' R package version 2.0, \doi{10.32614/CRAN.package.GAD}.
#'
#' Underwood, A.J. 1993. The mechanics of spatially replicated sampling programmes
#' to detect environmental impacts in a variable world. Aust. J. Ecol. 18,99-116.
#' \doi{10.1111/j.1442-9993.1993.tb00437.x}.
#'
#' Underwood, A.J. 1997. \emph{Experiments in Ecology: Their Logical Design and Interpretation Using
#' Analysis of Variance}. Cambridge University Press, Cambridge. \doi{10.1017/CBO9780511806407}.
#'
#' @seealso
#' \code{\link[AAV]{sym4asym}},
#' \code{\link[AAV]{asym.aci}},
#' \code{\link[AAV]{asym.baci}},
#' \code{\link[GAD]{gad}}
#'
#' @examples
#' # Example for BACI design containing 3 temporal factors and 2 spatial factor:
#' data(baci.3t2s)
#' ex.baci.3t2s<-aav(data = baci.3t2s, aav.design = "baci", n.ftemp = 3, n.fspac = 2,
#'     names.impact = "impact", names.before = "before",
#'     pooling = TRUE, method = "eliminate", p.aav = 0.05, p.pooling = 0.05)
#' ex.baci.3t2s
#'
#' # Example for post-impact design with only one disturbed location
#' # containing 2 temporal factors and 2 spatial factor:
#' data(aci.2t2s)
#' ex.aci.2t2s<-aav(data = aci.2t2s, aav.design = "aci", n.ftemp = 2, n.fspac = 2,
#'    names.impact = "impact", pooling = TRUE, p.aav = 0.05, p.pooling = 0.25)
#' ex.aci.2t2s
#'
#' # Example for post-impact design with more than one place with one
#' # disturbed location containing 0 temporal factors and 3 spatial factor:
#'  data(aci.0t3s)
#' ex.aci.0t3s<-aav(data = aci.0t3s, aav.design = "aci", n.ftemp = 0, n.fspac = 3,
#'     names.impact = "impact", pooling = TRUE, p.aav = 0.05, p.pooling = 0.25)
#' ex.aci.0t3s
#'
#' @import GAD
#'         dplyr
#'
#' @include sym4asym.R
#' @include asym.aci.R
#' @include asym.baci.R
#'
#' @export

aav<- function(data, aav.design, n.ftemp, n.fspac, names.impact, names.before = NULL,
               p.aav = 0.05, pooling = TRUE, method = "pool", p.pooling = 0.25,
               aci.tailed.test = 1) {

  # Stops
  if(aav.design %in% c("baci","aci") != T) stop("The 'aav.design' object does not exist. Check your command line.")
  if(names.impact %in% unique(data$location) != T) stop("The 'names.impact' object does not exist. Check your command line")

  md <- data.frame(data[,(n.ftemp+n.fspac+1):length(data)])
  if(length(md) == 1) names(md)<-colnames(data[length(data)])
  if(0 %in% colSums(md)) stop("At least one of the variable responses has only 0 values. Check your dataset.")

  if(aav.design == "aci") {
    if(n.ftemp >= 3 | n.fspac == 0 | n.fspac >= 4) stop("Select different values for 'n.ftemp' or 'n.fspac'.")

    if(n.ftemp == 2) names(data)[3]<-"location"
    if(n.ftemp < 2) names(data)[2]<-"location"
    md.ca<-data.frame(data[data$location != names.impact,(n.ftemp+n.fspac+1):length(data)])
    if(length(md.ca) == 1) names(md.ca)<-colnames(data[length(data)])
    if(0 %in% colSums(md.ca)) stop("At least one variable response has only 0 values in the 'Controls' location. Check your dataset.")
  }
  if(aav.design == "baci") {
    if(names.before %in% unique(data$before_after) != T) stop("The 'names.before' object does not exist. Check your command line")
    if(n.ftemp <= 1 | n.ftemp >= 4 | n.fspac == 0 | n.fspac >= 3) stop("Select different values for 'n.ftemp' or 'n.fspac'.")

    if(n.ftemp == 3) names(data)[1:4]<-c("before_after", "period", "time", "location")
    if(n.ftemp == 2) names(data)[1:3]<-c("before_after", "time", "location")
    md.ca<-data.frame(data[data$location != names.impact & data$before_after != names.before,(n.ftemp+n.fspac+1):length(data)])
    if(length(md.ca) == 1) names(md.ca)<-colnames(data[length(data)])
    if(0 %in% colSums(md.ca)) stop("At least one variable response has only 0 values in the 'Controls' location after the disturbance has occurred. Check your dataset.")
  }

  # Running the symmetrical ANOVAs
  sym.anova<-sym4asym(data = data, aav.design = aav.design,
           n.ftemp = n.ftemp, n.fspac = n.fspac,
           names.impact = names.impact, names.before = names.before)

  # Running the asymmetrical ANOVAs
  asym.anova_complete<-list()
  asym.anova_final<-list()
  for (i in 1:ncol(md)) {

    if(aav.design == "aci") {
      sym.F<-cbind(Source.of.Variation = rownames(sym.anova$full_design[[i]]), data.frame(sym.anova$full_design[[i]][1:2], row.names=NULL))
      sym.C<-cbind(Source.of.Variation = rownames(sym.anova$control_design[[i]]), data.frame(sym.anova$control_design[[i]][1:2], row.names=NULL))
      sym.NC<-NULL
        if(n.ftemp == 0) {
          sym.NC<-cbind(Source.of.Variation = rownames(sym.anova$nested_control_design[[i]]), data.frame(sym.anova$nested_control_design[[i]][1:2], row.names=NULL))
        }

      asym.anova_complete[[i]]<-asym.aci(data = data, n.ftemp = n.ftemp, n.fspac = n.fspac,
                               names.impact = names.impact,
                               sym.F= sym.F, sym.C = sym.C, sym.NC = sym.NC,
                               p.aav = p.aav, pooling = pooling, p.pooling = p.pooling,
                               aci.tailed.test = aci.tailed.test)
      names(asym.anova_complete)[i]<-colnames(md)[i]

      asym.anova_final[[i]]<-asym.anova_complete[[i]]
      names(asym.anova_final)[i]<-colnames(md)[i]
      asym.anova_final[[i]]<-asym.anova_final[[i]][asym.anova_final[[i]]$Final.Table == "Yes",]
      asym.anova_final[[i]]<-asym.anova_final[[i]][,names(asym.anova_final[[i]]) %in% c(
        "Source.of.Variation","Df","Sum.Sq","Mean.Sq","F.Value",
        "P.Value.Upper.Tail","P.Value.Lower.Tail","Next.Step", "Interpretation")]
      if(all(is.na(asym.anova_final[[i]]$P.Value.Lower.Tail))) {
        asym.anova_final[[i]]<-asym.anova_final[[i]][,!names(asym.anova_final[[i]]) %in% c("P.Value.Lower.Tail")]
      }
    }

    if(aav.design == "baci") {
        sym.F<-cbind(Source.of.Variation = rownames(sym.anova$full_design[[i]]), data.frame(sym.anova$full_design[[i]][1:2], row.names=NULL))
        sym.C<-cbind(Source.of.Variation = rownames(sym.anova$control_design[[i]]), data.frame(sym.anova$control_design[[i]][1:2], row.names=NULL))
        sym.B<-cbind(Source.of.Variation = rownames(sym.anova$before_design[[i]]), data.frame(sym.anova$before_design[[i]][1:2], row.names=NULL))
        sym.CB<-cbind(Source.of.Variation = rownames(sym.anova$control_before_design[[i]]), data.frame(sym.anova$control_before_design[[i]][1:2], row.names=NULL))

      asym.anova_complete[[i]]<-asym.baci(data = data, n.ftemp = n.ftemp, n.fspac = n.fspac,
                                names.impact = names.impact, names.before = names.before,
                                sym.F = sym.F, sym.C = sym.C, sym.B = sym.B, sym.CB = sym.CB,
                                p.aav = p.aav, pooling = pooling, method = method, p.pooling = p.pooling)
      names(asym.anova_complete)[i]<-colnames(md)[i]

      asym.anova_final[[i]]<-asym.anova_complete[[i]]
      names(asym.anova_final)[i]<-colnames(md)[i]
      asym.anova_final[[i]]<-asym.anova_final[[i]][asym.anova_final[[i]]$Final.Table == "Yes",]
      asym.anova_final[[i]]<-asym.anova_final[[i]][,names(asym.anova_final[[i]]) %in% c(
        "Source.of.Variation","Df","Sum.Sq","Mean.Sq","F.Value",
        "P.Value.Upper.Tail","P.Value.Lower.Tail","Next.Step", "Interpretation")]
      if(all(is.na(asym.anova_final[[i]]$P.Value.Lower.Tail))) {
        asym.anova_final[[i]]<-asym.anova_final[[i]][,!names(asym.anova_final[[i]]) %in% c("P.Value.Lower.Tail")]
      }
    }
  }
  aav.list<-list()
  aav.list[[1]]<-sym.anova
  names(aav.list)[[1]]<-"sym.anova"
  aav.list[[2]]<-asym.anova_final
  names(aav.list)[[2]]<-"asym.anova_final"
  aav.list[[3]]<-asym.anova_complete
  names(aav.list)[[2]]<-"asym.anova_complete"

  return(aav.list)
}

