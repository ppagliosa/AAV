#' @title Asymmetrical ANOVA From Post-Impact Design
#'
#' @description
#' Provides an asymmetrical ANOVA from post-impact design with only one
#' disturbed location and multiple control locations or
#' a post-impact design involving more than one place with one
#' disturbed location and multiple control.
#'
#' Post-impact design with only one disturbed location can be used for models
#' with 1 or 2 temporal factors and 1 or 2 spatial factors.
#' Post-impact design involving more than one place with one disturbed location
#' can be used for models with 2 or 3 spatial factors.
#'
#' @param data a data frame with temporal factors, spatial factors and variable(s).
#' @param n.ftemp number of temporal factors.
#' @param n.fspac number of spatial factors.
#' @param names.impact name of the asymmetrical level of the Location factor (e.g. 'impact').
#' @param sym.F an object containing the results returned by \code{sym4asym} function or
#'     an symmetrical analysis of variance table run with all the data,
#'     meaning the asymmetry is deal as a factor level (e.g., Control_1,
#'     Control_2 and Impact are the factor levels).
#' @param sym.C an object containing the results returned by \code{sym4asym} function or
#'     an symmetrical analysis of variance table run omitting the asymmetrical data
#'     (e.g., only the Control_1 and Control_2 are the factor levels).
#' @param sym.NC only for a design with more than one place with one disturbed location.
#'     An object containing the results returned by \code{sym4asym} function or
#'     an symmetrical analysis of variance table run omitting the asymmetrical data
#'     and with all factors nested (e.g., Control location nested in Place).
#' @param p.aav level of significance for the asymmetrical ANOVA. Default is \code{0.05}.
#' @param pooling remove a term from the analysis using a post-hoc pooling procedure. Default is \code{TRUE}.
#' @param p.pooling only when \code{pooling = TRUE}.
#'     Level of significance for the pooling procedure. Default is \code{0.25}.
#' @param aci.tailed.test only for post-impact design with only one disturbed location
#'     and multiple control locations. Perform a 1-tailed or a 2-tailed test for all terms.
#'     Default is \code{1}.
#'
#' @details
#' Asymmetrical analyses of variance are constructed by combining
#' the Sum of Squares values from separate (symmetrical) ANOVAs.
#' In post-impact designs with only one disturbed location two fully orthogonal
#' ANOVAs are performed: one on all the data, and one on the 'control' data.
#' In post-impact designs involving more than one place with one disturbed location,
#' two fully orthogonal and one fully nested ANOVA need to be performed: one on all the data
#' (orthogonal) and two on the 'control' data (orthogonal and nested).
#' \code{aav} uses \code{\link[GAD]{GAD}} (Sandrini-Neto et al. 2024) to perform symmetrical ANOVAs.
#'
#' The \code{asym.aci} analysis fits mixed asymmetrical ANOVA designs with a specific
#' combination of orthogonal/nested and fixed/random factors:
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
#' In the post-impact design with only one disturbed location models,
#' the first factor and the Location factor are fixed and all others are random.
#' In the post-impact design involving than one place with one disturbed location,
#' the first factor is fixed and the others are random.
#'
#' The data frame must follow a specific structure:
#' 1- first the temporal factors, then the spatial factors, and then the variable(s);
#' 2- each factor type must appear in a hierarchical order.
#'
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
#' \code{pooling} is \code{TRUE} the analysis is running using the method \code{pool}.
#'
#' In the case of a pooling structure, when \code{pooling = TRUE},
#' the analysis runs using the \code{pool} method, where the
#' degrees of freedom (Df) and sums of squares (SS) for that term are
#' pooled with the term(s) that have an equivalent mean square expectation
#' after the component of variation of that term has been set to zero. This
#' component of variation is set to zero if it is >= \code{p.pooling} value.
#'
#' The magnitude of interactions in time and space can be analysed using a two-tailed tests.
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
#' Underwood, A.J. 1993. The mechanics of spatially replicated sampling programmes
#' to detect environmental impacts in a variable world. Aust. J. Ecol. 18,99-116.
#' \doi{10.1111/j.1442-9993.1993.tb00437.x}.
#'
#' Underwood, A.J. 1997. \emph{Experiments in Ecology: Their Logical Design and Interpretation Using
#' Analysis of Variance}. Cambridge University Press, Cambridge. \doi{10.1017/CBO9780511806407}.
#'
#' @seealso
#' \code{\link[AAV]{sym4asym}},
#' \code{\link[AAV]{asym.baci}},
#' \code{\link[AAV]{aav}}
#'
#' @examples
#' # Example for post-impact design with only one disturbed location
#' # containing 2 temporal factors and 2 spatial factor:
#' data(aci.2t2s) # 4-factors dataset
#' data(symANOVA_aci.2t2s) # symmetrical ANOVAS from aci.2t2s data
#' ex.symANOVA_aci.2t2s<-asym.aci(data = aci.2t2s, n.ftemp = 2, n.fspac = 2,
#'    names.impact = "impact",
#'    sym.F = symANOVA_aci.2t2s$full_design,
#'    sym.C = symANOVA_aci.2t2s$control_design,
#'    p.aav = 0.05, pooling = TRUE, p.pooling = 0.25, aci.tailed.test = 1)
#' ex.symANOVA_aci.2t2s
#'
#' # Example for post-impact design with more than one place with one
#' # disturbed location containing 0 temporal factors and 3 spatial factor:
#' data(aci.0t3s) # 3-factors dataset
#' data(symANOVA_aci.0t3s) # symmetrical ANOVAS from aci.0t3s data
#' ex.symANOVA_aci.0t3s<-asym.aci(data = aci.0t3s, n.ftemp = 0, n.fspac = 3,
#'     names.impact = "impact",
#'     sym.F = symANOVA_aci.0t3s$full_design,
#'     sym.C = symANOVA_aci.0t3s$control_design,
#'     sym.NC = symANOVA_aci.0t3s$nested_control_design,
#'     p.aav = 0.05, pooling = TRUE, p.pooling = 0.25)
#' ex.symANOVA_aci.0t3s
#'
#' @import dplyr
#'
#' @include sym4asym.R
#'
#' @export

asym.aci<- function(data, n.ftemp, n.fspac, names.impact,
                    sym.F, sym.C, sym.NC = NULL,
                    p.aav = 0.05, pooling = TRUE, p.pooling = 0.25, aci.tailed.test = 1) {

  # Stops
  if(n.ftemp >= 3 | n.fspac == 0 | n.fspac >= 4) stop("Select different values for 'n.ftemp' or 'n.fspac'.")

  # Factors names
  if(n.ftemp == 2 & n.fspac == 2) {
    names(data)[1:4]<-c("period", "time", "location", "site")
    aav.model <-"BEGC"
    aav.model.lm <-c("B","E","BE","BC","EG","BEC","BEG","BEGC")
  }
  if(n.ftemp == 2 & n.fspac == 1) {
    names(data)[1:3]<-c("period", "time", "location")
    aav.model <-"BEC"
    aav.model.lm <-c("B","E","BE","BC","BEC")
  }
  if(n.ftemp == 1 & n.fspac == 2) {
    names(data)[1:3]<-c("time", "location", "site")
    aav.model <-"CEG"
    aav.model.lm <-c("C","E","CE","EG","CEG")
  }
  if(n.ftemp == 1 & n.fspac == 1) {
    names(data)[1:2]<-c("time", "location")
    aav.model <-"CE"
    aav.model.lm <-c("C","E","CE")
  }
  if(n.ftemp == 0 & n.fspac == 2) {
    names(data)[1:2]<-c("place", "location")
    aav.model <-"DE"
    aav.model.lm <-c("D","E","DE")
  }
  if(n.ftemp == 0 & n.fspac == 3) {
    names(data)[1:3]<-c("place", "location", "site")
    aav.model <-"DEG"
    aav.model.lm <-c("D","E","DE","EG","DEG")
  }

  # Additional Stops
  md <- data.frame(data[,(n.ftemp+n.fspac+1):length(data)])
  if(length(md) == 1) names(md)<-colnames(data[length(data)])
  if(0 %in% colSums(md)) stop("At least one of the variable responses has only 0 values. Check your dataset.")
  md.ca<-data.frame(data[data$location != names.impact,(n.ftemp+n.fspac+1):length(data)])
  if(length(md.ca) == 1) names(md.ca)<-colnames(data[length(data)])
  if(0 %in% colSums(md.ca)) stop("At least one variable response has only 0 values in the 'Controls' location. Check your dataset.")

  # prevenir "no visible binding for global variable"
  ID<-NULL

  ### Cria tabela vazia para colocar os dados
  # importante especificar o tipo de dados para mais tarde poder juntar as linhas
  # (em geral, as funcoes so juntam dados de mesmo tipo)
  aci<-NULL
  aci<-data.frame(matrix(data=character(0), nrow = 0, ncol=19))
  names(aci)<-c("ID","Source.of.Variation","Df.Formula",
                "Df", "Sum.Sq", "Calculated.From", "Mean.Sq",
                "Denominator","Df.Pooled","Sum.Sq.Pooled","Mean.Sq.Pooled",
                "F.Value", "P.Value.Upper.Tail","P.Value.Lower.Tail",
                "Significant", "Post.Hoc.Pooling", "Next.Step", "Interpretation", "Final.Table")
  aci$Df<-as.numeric(aci$Df)
  aci$Sum.Sq<-as.numeric(aci$Sum.Sq)
  aci$Mean.Sq<-as.numeric(aci$Mean.Sq)
  aci$Df.Pooled<-as.numeric(aci$Df.Pooled)
  aci$Sum.Sq.Pooled<-as.numeric(aci$Sum.Sq.Pooled)
  aci$Mean.Sq.Pooled<-as.numeric(aci$Mean.Sq.Pooled)
  aci$F.Value<-as.numeric(aci$F.Value)
  aci$P.Value.Upper.Tail<-as.numeric(aci$P.Value.Upper.Tail)
  aci$P.Value.Lower.Tail<-as.numeric(aci$P.Value.Lower.Tail)

  ###----------------------------------------------------------- FROM sym.F
  ##---------------------------------------------- code xx.0000
  colnames(sym.F)[1:3] <- c("Source.of.Variation","Df", "Sum.Sq")
  colnames(sym.C)[1:3] <- c("Source.of.Variation","Df", "Sum.Sq")
  if(is.null(sym.NC) == F) colnames(sym.NC)[1:3] <- c("Source.of.Variation","Df", "Sum.Sq")

  f.an<-aci[F,]
  if("B" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "B",1:3])
    f.an$ID[f.an$Source.of.Variation == "B"]<-"01B.0000"
    f.an$Df.Formula[f.an$ID == "01B.0000"]<-"(p-1)"
    f.an$Calculated.From[f.an$ID == "01B.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "01B.0000"]<-"Periods = P"
  }
  if("C" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "C",1:3])
    f.an$ID[f.an$Source.of.Variation == "C"]<-"01C.0000"
    f.an$Df.Formula[f.an$ID == "01C.0000"]<-"(t-1)"
    f.an$Calculated.From[f.an$ID == "01C.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "01C.0000"]<-"Times = (T)"
  }
  if("E" %in% aav.model.lm == T & "D" %in% aav.model.lm != T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "E",1:3])
    f.an$ID[f.an$Source.of.Variation == "E"]<-"02E.0000"
    f.an$Df.Formula[f.an$ID == "02E.0000"]<-"(l-1)"
    f.an$Calculated.From[f.an$ID == "02E.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "02E.0000"]<-"Locations = L"
  }
  if("BE" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "B:E",1:3])
    f.an$ID[f.an$Source.of.Variation == "B:E"]<-"03BE.0000"
    f.an$Df.Formula[f.an$ID == "03BE.0000"]<-"(p-1)(l-1)"
    f.an$Calculated.From[f.an$ID == "03BE.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "03BE.0000"]<-"P x L"
  }
  if("EG" %in% aav.model.lm == T & "D" %in% aav.model.lm != T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "E:G",1:3])
    f.an$ID[f.an$Source.of.Variation == "E:G"]<-"05EG.0000"
    f.an$Df.Formula[f.an$ID == "05EG.0000"]<-"l(s-1)"
    f.an$Calculated.From[f.an$ID == "05EG.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "05EG.0000"]<-"Sites(L) = S(L)"
  }
  if("BC" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "B:C",1:3])
    f.an$ID[f.an$Source.of.Variation == "B:C"]<-"06BC.0000"
    f.an$Df.Formula[f.an$ID == "06BC.0000"]<-"p(t-1)"
    f.an$Calculated.From[f.an$ID == "06BC.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "06BC.0000"]<-"Times(P) = T(P)"
  }
  if("BEG" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "B:E:G",1:3])
    f.an$ID[f.an$Source.of.Variation == "B:E:G"]<-"09BEG.0000"
    f.an$Df.Formula[f.an$ID == "09BEG.0000"]<-"(t-1)l(s-1)"
    f.an$Calculated.From[f.an$ID == "09BEG.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "09BEG.0000"]<-"P x S(L)"
  }
  if("BEC" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "B:E:C",1:3])
    f.an$ID[f.an$Source.of.Variation == "B:E:C"]<-"10BEC.0000"
    f.an$Df.Formula[f.an$ID == "10BEC.0000"]<-"p(t-1)(l-1)"
    f.an$Calculated.From[f.an$ID == "10BEC.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "10BEC.0000"]<-"T(P) x L"
  }
  if("BEGC" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "B:E:G:C",1:3])
    f.an$ID[f.an$Source.of.Variation == "B:E:G:C"]<-"11BEGC.0000"
    f.an$Df.Formula[f.an$ID == "11BEGC.0000"]<-"p(t-1)l(s-1)"
    f.an$Calculated.From[f.an$ID == "11BEGC.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "11BEGC.0000"]<-"T(P) x S(L)"
  }
  if("CE" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "C:E",1:3])
    f.an$ID[f.an$Source.of.Variation == "C:E"]<-"10CE.0000"
    f.an$Df.Formula[f.an$ID == "10CE.0000"]<-"(t-1)(l-1)"
    f.an$Calculated.From[f.an$ID == "10CE.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "10CE.0000"]<-"T x L"
  }
  if("CEG" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "C:E:G",1:3])
    f.an$ID[f.an$Source.of.Variation == "C:E:G"]<-"11CEG.0000"
    f.an$Df.Formula[f.an$ID == "11CEG.0000"]<-"(t-1)l(s-1)"
    f.an$Calculated.From[f.an$ID == "11CEG.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "11CEG.0000"]<-"T x S(L)"
  }
  # Residuals
  f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "Residuals",1:3])
  f.an$ID[f.an$Source.of.Variation == "Residuals"]<-"12R.0000"
  f.an$Calculated.From[f.an$ID == "12R.0000"]<-"sym.F"
  f.an$Mean.Sq[f.an$ID == "12R.0000"]<-
    f.an$Sum.Sq[f.an$ID == "12R.0000"] / f.an$Df[f.an$ID == "12R.0000"]
  if(aav.model == "BEGC") f.an$Df.Formula[f.an$ID == "12R.0000"]<-"ptls(n-1)"
  if(aav.model == "BEC") f.an$Df.Formula[f.an$ID == "12R.0000"]<-"ptl(n-1)"
  if(aav.model == "CE") f.an$Df.Formula[f.an$ID == "12R.0000"]<-"tl(n-1)"
  if(aav.model == "CEG") f.an$Df.Formula[f.an$ID == "12R.0000"]<-"tls(n-1)"
  if(aav.model == "DE") f.an$Df.Formula[f.an$ID == "12R.0000"]<-"pt(n-1)"
  if(aav.model == "DEG") f.an$Df.Formula[f.an$ID == "12R.0000"]<-"pts(n-1)"

  ###----------------------------------------- AMONG Controls = From sym.C
  ##------------------------- code xx.2000
  c.an<-aci[F,]
  ###----------------------------------- IMPACTs = From sym.F - sym.C
  ##------------------------- codigo xx.1000
  i.an<-aci[F,]

  ###------------------ INTERACTIONS with TWO temporal and TWO spatial factors

  if(aav.model == "BEGC") {
    # 	F = (MST(P)xS(CC) / MSRes)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "B:E:G:C",1:3])
    c.an$ID[c.an$Source.of.Variation == "B:E:G:C"]<-"11BEGC.0040"
    c.an$Df.Formula[c.an$ID == "11BEGC.0040"]<-"p(t-1)c(s-1)"
    c.an$Calculated.From[c.an$ID == "11BEGC.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "11BEGC.0040"] <- c.an$Sum.Sq[c.an$ID == "11BEGC.0040"] / c.an$Df[c.an$ID == "11BEGC.0040"]
    c.an$Source.of.Variation[c.an$ID == "11BEGC.0040"]<-"T(P) x S(C)"
    c.an$Denominator[c.an$ID == "11BEGC.0040"]<-"Residuals"
    c.an$F.Value[c.an$ID == "11BEGC.0040"]<-c.an$Mean.Sq[c.an$ID == "11BEGC.0040"] / f.an$Mean.Sq[f.an$ID == "12R.0000"]
    c.an$P.Value.Upper.Tail[c.an$ID == "11BEGC.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "11BEGC.0040"], c.an$Df[c.an$ID == "11BEGC.0040"], f.an$Df[f.an$ID == "12R.0000"], lower.tail = F)

    # 	F = (MST(P)xS(IC) / MST(P)xS(CC)), or
    dii.11BEGC.0030<-cbind(Source.of.Variation = "B:E:G:C",
                           sym.F[sym.F$Source.of.Variation == "B:E:G:C",2:3] - sym.C[sym.C$Source.of.Variation == "B:E:G:C",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.11BEGC.0030)
    i.an$ID[i.an$Source.of.Variation == "B:E:G:C"]<-"11BEGC.0030"
    i.an$Df.Formula[i.an$ID == "11BEGC.0030"]<-"p(t-1)(s-1)"
    i.an$Calculated.From[i.an$ID == "11BEGC.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "11BEGC.0030"]<-i.an$Sum.Sq[i.an$ID == "11BEGC.0030"] / i.an$Df[i.an$ID == "11BEGC.0030"]
    i.an$Source.of.Variation[i.an$ID == "11BEGC.0030"]<-"T(P) x S(I)"
    i.an$Denominator[i.an$ID == "11BEGC.0030"]<-"T(P)xS(C)"
    i.an$F.Value[i.an$ID == "11BEGC.0030"]<-i.an$Mean.Sq[i.an$ID == "11BEGC.0030"] / c.an$Mean.Sq[c.an$ID == "11BEGC.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "11BEGC.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "11BEGC.0030"], i.an$Df[i.an$ID == "11BEGC.0030"], c.an$Df[c.an$ID == "11BEGC.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "11BEGC.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "11BEGC.0030"], i.an$Df[i.an$ID == "11BEGC.0030"], c.an$Df[c.an$ID == "11BEGC.0040"], lower.tail = T)

    # 	F = (MST(P)xCC / MST(P)xS(CC))  10BEC.0040
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "B:E:C",1:3])
    c.an$ID[c.an$Source.of.Variation == "B:E:C"]<-"10BEC.0040"
    c.an$Df.Formula[c.an$ID == "10BEC.0040"]<-"p(t-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "10BEC.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "10BEC.0040"] <- c.an$Sum.Sq[c.an$ID == "10BEC.0040"] / c.an$Df[c.an$ID == "10BEC.0040"]
    c.an$Source.of.Variation[c.an$ID == "10BEC.0040"]<-"T(P) x C"
    c.an$Denominator[c.an$ID == "10BEC.0040"]<-"T(P)xS(C)"
    c.an$F.Value[c.an$ID == "10BEC.0040"]<-c.an$Mean.Sq[c.an$ID == "10BEC.0040"] / c.an$Mean.Sq[c.an$ID == "11BEGC.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "10BEC.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "10BEC.0040"], c.an$Df[c.an$ID == "10BEC.0040"], c.an$Df[c.an$ID == "11BEGC.0040"], lower.tail = F)

    # 	F = (MST(P)xIC) / MST(P)xCC), or
    dii.10BEC.0030<-cbind(Source.of.Variation = "B:E:C",
                          sym.F[sym.F$Source.of.Variation == "B:E:C",2:3] - sym.C[sym.C$Source.of.Variation == "B:E:C",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.10BEC.0030)
    i.an$ID[i.an$Source.of.Variation == "B:E:C"]<-"10BEC.0030"
    i.an$Df.Formula[i.an$ID == "10BEC.0030"]<-"p(t-1)"
    i.an$Calculated.From[i.an$ID == "10BEC.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "10BEC.0030"]<-i.an$Sum.Sq[i.an$ID == "10BEC.0030"] / i.an$Df[i.an$ID == "10BEC.0030"]
    i.an$Source.of.Variation[i.an$ID == "10BEC.0030"]<-"T(P) x I"
    i.an$Denominator[i.an$ID == "10BEC.0030"]<-"T(P)xC"
    i.an$F.Value[i.an$ID == "10BEC.0030"]<-i.an$Mean.Sq[i.an$ID == "10BEC.0030"] / c.an$Mean.Sq[c.an$ID == "10BEC.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "10BEC.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "10BEC.0030"], i.an$Df[i.an$ID == "10BEC.0030"], c.an$Df[c.an$ID == "10BEC.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "10BEC.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "10BEC.0030"], i.an$Df[i.an$ID == "10BEC.0030"], c.an$Df[c.an$ID == "10BEC.0040"], lower.tail = T)

    # 	F = (MSPxS(CC) / MST(P)xS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "B:E:G",1:3])
    c.an$ID[c.an$Source.of.Variation == "B:E:G"]<-"09BEG.0040"
    c.an$Df.Formula[c.an$ID == "09BEG.0040"]<-"(p-1)c(s-1)"
    c.an$Calculated.From[c.an$ID == "09BEG.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "09BEG.0040"] <- c.an$Sum.Sq[c.an$ID == "09BEG.0040"] / c.an$Df[c.an$ID == "09BEG.0040"]
    c.an$Source.of.Variation[c.an$ID == "09BEG.0040"]<-"P x S(C)"
    c.an$Denominator[c.an$ID == "09BEG.0040"]<-"T(P)xS(C)"
    c.an$F.Value[c.an$ID == "09BEG.0040"]<-c.an$Mean.Sq[c.an$ID == "09BEG.0040"] / c.an$Mean.Sq[c.an$ID == "11BEGC.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "09BEG.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "09BEG.0040"], c.an$Df[c.an$ID == "09BEG.0040"], c.an$Df[c.an$ID == "11BEGC.0040"], lower.tail = F)

    # 	F = (MSPxS(IC) / MSPxS(CC)), or
    dii.09BEG.0030<-cbind(Source.of.Variation = "B:E:G",
                          sym.F[sym.F$Source.of.Variation == "B:E:G",2:3] - sym.C[sym.C$Source.of.Variation == "B:E:G",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.09BEG.0030)
    i.an$ID[i.an$Source.of.Variation == "B:E:G"]<-"09BEG.0030"
    i.an$Df.Formula[i.an$ID == "09BEG.0030"]<-"(p-1)(s-1)"
    i.an$Calculated.From[i.an$ID == "09BEG.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "09BEG.0030"]<-i.an$Sum.Sq[i.an$ID == "09BEG.0030"] / i.an$Df[i.an$ID == "09BEG.0030"]
    i.an$Source.of.Variation[i.an$ID == "09BEG.0030"]<-"P x S(I)"
    i.an$Denominator[i.an$ID == "09BEG.0030"]<-"PxS(C)"
    i.an$F.Value[i.an$ID == "09BEG.0030"]<-i.an$Mean.Sq[i.an$ID == "09BEG.0030"] / c.an$Mean.Sq[c.an$ID == "09BEG.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "09BEG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "09BEG.0030"], i.an$Df[i.an$ID == "09BEG.0030"], c.an$Df[c.an$ID == "09BEG.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "09BEG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "09BEG.0030"], i.an$Df[i.an$ID == "09BEG.0030"], c.an$Df[c.an$ID == "09BEG.0040"], lower.tail = T)

    # 	F = (MSS(CC) / MSPxS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E:G",1:3])
    c.an$ID[c.an$Source.of.Variation == "E:G"]<-"05EG.0040"
    c.an$Df.Formula[c.an$ID == "05EG.0040"]<-"c(s-1)"
    c.an$Calculated.From[c.an$ID == "05EG.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "05EG.0040"] <- c.an$Sum.Sq[c.an$ID == "05EG.0040"] / c.an$Df[c.an$ID == "05EG.0040"]
    c.an$Source.of.Variation[c.an$ID == "05EG.0040"]<-"S(C)"
    c.an$Denominator[c.an$ID == "05EG.0040"]<-"PxS(C)"
    c.an$F.Value[c.an$ID == "05EG.0040"]<-c.an$Mean.Sq[c.an$ID == "05EG.0040"] / c.an$Mean.Sq[c.an$ID == "09BEG.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "05EG.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "05EG.0040"], c.an$Df[c.an$ID == "05EG.0040"], c.an$Df[c.an$ID == "09BEG.0040"], lower.tail = F)

    # 	F = (MSS(IC) / MSS(CC)), or
    dii.05EG.0030<-cbind(Source.of.Variation = "E:G",
                         sym.F[sym.F$Source.of.Variation == "E:G",2:3] - sym.C[sym.C$Source.of.Variation == "E:G",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.05EG.0030)
    i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0030"
    i.an$Df.Formula[i.an$ID == "05EG.0030"]<-"(s-1)"
    i.an$Calculated.From[i.an$ID == "05EG.0030"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$ID == "05EG.0030"]<-"S(I)"
    i.an$Denominator[i.an$ID == "05EG.0030"]<-"S(C)"
    i.an$Mean.Sq[i.an$ID == "05EG.0030"]<-i.an$Sum.Sq[i.an$ID == "05EG.0030"] / i.an$Df[i.an$ID == "05EG.0030"]
    i.an$F.Value[i.an$ID == "05EG.0030"]<-i.an$Mean.Sq[i.an$ID == "05EG.0030"] / c.an$Mean.Sq[c.an$ID == "05EG.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "05EG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0030"], i.an$Df[i.an$ID == "05EG.0030"], c.an$Df[c.an$ID == "05EG.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "05EG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0030"], i.an$Df[i.an$ID == "05EG.0030"], c.an$Df[c.an$ID == "05EG.0040"], lower.tail = T)

    # 	If varPxS(CC) <= 0 and varT(P)xCC > 0), then F = (MSPxCC / MST(P)xCC)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "B:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "B:E"]<-"03BE.0040"
    c.an$Df.Formula[c.an$ID == "03BE.0040"]<-"(p-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03BE.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03BE.0040"] <- c.an$Sum.Sq[c.an$ID == "03BE.0040"] / c.an$Df[c.an$ID == "03BE.0040"]
    c.an$Source.of.Variation[c.an$ID == "03BE.0040"]<-"P x C"
    c.an$Denominator[c.an$ID == "03BE.0040"]<-"T(P)xC"
    c.an$F.Value[c.an$ID == "03BE.0040"]<-c.an$Mean.Sq[c.an$ID == "03BE.0040"] / c.an$Mean.Sq[c.an$ID == "10BEC.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03BE.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "03BE.0040"], c.an$Df[c.an$ID == "03BE.0040"], c.an$Df[c.an$ID == "10BEC.0040"], lower.tail = F)

    # 	If varT(P)xCC <= 0 and varPxS(CC) > 0, then F = (MSPxCC / MSPxS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "B:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "B:E"]<-"03BE.0041"
    c.an$Df.Formula[c.an$ID == "03BE.0041"]<-"(p-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03BE.0041"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03BE.0041"] <- c.an$Sum.Sq[c.an$ID == "03BE.0041"] / c.an$Df[c.an$ID == "03BE.0041"]
    c.an$Source.of.Variation[c.an$ID == "03BE.0041"]<-"P x C"
    c.an$Denominator[c.an$ID == "03BE.0041"]<-"PxS(C)"
    c.an$F.Value[c.an$ID == "03BE.0041"]<-c.an$Mean.Sq[c.an$ID == "03BE.0041"] / c.an$Mean.Sq[c.an$ID == "09BEG.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03BE.0041"]<- stats::pf(c.an$F.Value[c.an$ID == "03BE.0041"], c.an$Df[c.an$ID == "03BE.0041"], c.an$Df[c.an$ID == "09BEG.0040"], lower.tail = F)

    # 	If varT(P)xCC <= 0, varPxS(CC) <= 0, then F = (MSPxCC / MSpooled), where MSpooled = (SST(P)x CC + SSPxS(CC) + SST(P)xS(CC)) / (dfT(P)x CC + dfPxS(CC) + dfT(P)xS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "B:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "B:E"]<-"03BE.0042"
    c.an$Df.Formula[c.an$ID == "03BE.0042"]<-"(p-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03BE.0042"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03BE.0042"] <- c.an$Sum.Sq[c.an$ID == "03BE.0042"] / c.an$Df[c.an$ID == "03BE.0042"]
    c.an$Source.of.Variation[c.an$ID == "03BE.0042"]<-"P x C"
    c.an$Denominator[c.an$ID == "03BE.0042"]<-"T(P)xC+PxS(C)+T(P)xS(C)"
    c.an$Df.Pooled[c.an$ID == "03BE.0042"]<- c.an$Df[c.an$ID == "10BEC.0040"] + c.an$Df[c.an$ID == "09BEG.0040"] + c.an$Df[c.an$ID == "11BEGC.0040"]
    c.an$Sum.Sq.Pooled[c.an$ID == "03BE.0042"]<- c.an$Sum.Sq[c.an$ID == "10BEC.0040"] + c.an$Sum.Sq[c.an$ID == "09BEG.0040"] + c.an$Sum.Sq[c.an$ID == "11BEGC.0040"]
    c.an$Mean.Sq.Pooled[c.an$ID == "03BE.0042"]<- c.an$Sum.Sq.Pooled[c.an$ID == "03BE.0042"] / c.an$Df.Pooled[c.an$ID == "03BE.0042"]
    c.an$F.Value[c.an$ID == "03BE.0042"]<- c.an$Mean.Sq[c.an$ID == "03BE.0042"] / c.an$Mean.Sq.Pooled[c.an$ID == "03BE.0042"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03BE.0042"]<- stats::pf(c.an$F.Value[c.an$ID == "03BE.0042"], c.an$Df[c.an$ID == "03BE.0042"], c.an$Df.Pooled[c.an$ID == "03BE.0042"], lower.tail = F)

    # 	F = (MSPxIC) / MSPxCC), or
    dii.03BE.0030<-cbind(Source.of.Variation = "B:E",
                         sym.F[sym.F$Source.of.Variation == "B:E",2:3] - sym.C[sym.C$Source.of.Variation == "B:E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.03BE.0030)
    i.an$ID[i.an$Source.of.Variation == "B:E"]<-"03BE.0030"
    i.an$Df.Formula[i.an$ID == "03BE.0030"]<-"(p-1)"
    i.an$Calculated.From[i.an$ID == "03BE.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "03BE.0030"]<-i.an$Sum.Sq[i.an$ID == "03BE.0030"] / i.an$Df[i.an$ID == "03BE.0030"]
    i.an$Source.of.Variation[i.an$ID == "03BE.0030"]<-"P x I"
    i.an$Denominator[i.an$ID == "03BE.0030"]<-"PxC"
    i.an$F.Value[i.an$ID == "03BE.0030"]<-i.an$Mean.Sq[i.an$ID == "03BE.0030"] / c.an$Mean.Sq[c.an$ID == "03BE.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "03BE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0030"], i.an$Df[i.an$ID == "03BE.0030"], c.an$Df[c.an$ID == "03BE.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "03BE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0030"], i.an$Df[i.an$ID == "03BE.0030"], c.an$Df[c.an$ID == "03BE.0040"], lower.tail = T)

    # 	If varS(CC) <= 0 and varPxCC > 0, then F = (MSCC / MSPxCC)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0040"
    c.an$Df.Formula[c.an$ID == "02E.0040"]<-"(c-1)"
    c.an$Calculated.From[c.an$ID == "02E.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "02E.0040"] <- c.an$Sum.Sq[c.an$ID == "02E.0040"] / c.an$Df[c.an$ID == "02E.0040"]
    c.an$Source.of.Variation[c.an$ID == "02E.0040"]<-"Among Controls = C"
    c.an$Denominator[c.an$ID == "02E.0040"]<-"PxC"
    c.an$F.Value[c.an$ID == "02E.0040"]<-c.an$Mean.Sq[c.an$ID == "02E.0040"] / c.an$Mean.Sq[c.an$ID == "03BE.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "02E.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "02E.0040"], c.an$Df[c.an$ID == "02E.0040"], c.an$Df[c.an$ID == "03BE.0040"], lower.tail = F)

    # 	If varPxCC <= 0 and varS(CC) > 0, then F = (MSCC / MSS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0041"
    c.an$Df.Formula[c.an$ID == "02E.0041"]<-"(c-1)"
    c.an$Calculated.From[c.an$ID == "02E.0041"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "02E.0041"] <- c.an$Sum.Sq[c.an$ID == "02E.0041"] / c.an$Df[c.an$ID == "02E.0041"]
    c.an$Source.of.Variation[c.an$ID == "02E.0041"]<-"Among Controls = C"
    c.an$Denominator[c.an$ID == "02E.0041"]<-"S(C)"
    c.an$F.Value[c.an$ID == "02E.0041"]<-c.an$Mean.Sq[c.an$ID == "02E.0041"] / c.an$Mean.Sq[c.an$ID == "05EG.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "02E.0041"]<- stats::pf(c.an$F.Value[c.an$ID == "02E.0041"], c.an$Df[c.an$ID == "02E.0041"], c.an$Df[c.an$ID == "05EG.0040"], lower.tail = F)

    # 	If varPxCC <= 0, varS(CC) <=, varPxS(CC) <= 0 and varT(P)xCC > 0, then F = (MSCC / MSpooled), Where MSpooled = (SSPxCC + SSS(CC) + SST(P)xCC) / (dfPxCC + dfS(CC) + dfT(P)xCC)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0042"
    c.an$Df.Formula[c.an$ID == "02E.0042"]<-"(c-1)"
    c.an$Calculated.From[c.an$ID == "02E.0042"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "02E.0042"] <- c.an$Sum.Sq[c.an$ID == "02E.0042"] / c.an$Df[c.an$ID == "02E.0042"]
    c.an$Source.of.Variation[c.an$ID == "02E.0042"]<-"Among Controls = C"
    c.an$Denominator[c.an$ID == "02E.0042"]<-"PxC+S(C)+T(P)xC"
    c.an$Df.Pooled[c.an$ID == "02E.0042"]<- c.an$Df[c.an$ID == "03BE.0040"] + c.an$Df[c.an$ID == "05EG.0040"] + c.an$Df[c.an$ID == "10BEC.0040"]
    c.an$Sum.Sq.Pooled[c.an$ID == "02E.0042"]<- c.an$Sum.Sq[c.an$ID == "03BE.0040"] + c.an$Sum.Sq[c.an$ID == "05EG.0040"] + c.an$Sum.Sq[c.an$ID == "10BEC.0040"]
    c.an$Mean.Sq.Pooled[c.an$ID == "02E.0042"]<- c.an$Sum.Sq.Pooled[c.an$ID == "02E.0042"] / c.an$Df.Pooled[c.an$ID == "02E.0042"]
    c.an$F.Value[c.an$ID == "02E.0042"]<- c.an$Mean.Sq[c.an$ID == "02E.0042"] / c.an$Mean.Sq.Pooled[c.an$ID == "02E.0042"]
    c.an$P.Value.Upper.Tail[c.an$ID == "02E.0042"]<- stats::pf(c.an$F.Value[c.an$ID == "02E.0042"], c.an$Df[c.an$ID == "02E.0042"], c.an$Df.Pooled[c.an$ID == "02E.0042"], lower.tail = F)

    # 	If varPxCC <= 0, varS(CC) <=, varT(P)xCC <= 0 and varPxS(CC) > 0, then F = (MSCC / MSpooled), Where MSpooled = (SSPxCC + SSS(CC) + SSPxS(CC)) / (dfPxCC + dfS(CC) + dfPxS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0043"
    c.an$Df.Formula[c.an$ID == "02E.0043"]<-"(c-1)"
    c.an$Calculated.From[c.an$ID == "02E.0043"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "02E.0043"] <- c.an$Sum.Sq[c.an$ID == "02E.0043"] / c.an$Df[c.an$ID == "02E.0043"]
    c.an$Source.of.Variation[c.an$ID == "02E.0043"]<-"Among Controls = C"
    c.an$Denominator[c.an$ID == "02E.0043"]<-"PxC+S(C)+PxS(C)"
    c.an$Df.Pooled[c.an$ID == "02E.0043"]<- c.an$Df[c.an$ID == "03BE.0040"] + c.an$Df[c.an$ID == "05EG.0040"] + c.an$Df[c.an$ID == "09BEG.0040"]
    c.an$Sum.Sq.Pooled[c.an$ID == "02E.0043"]<- c.an$Sum.Sq[c.an$ID == "03BE.0040"] + c.an$Sum.Sq[c.an$ID == "05EG.0040"] + c.an$Sum.Sq[c.an$ID == "09BEG.0040"]
    c.an$Mean.Sq.Pooled[c.an$ID == "02E.0043"]<- c.an$Sum.Sq.Pooled[c.an$ID == "02E.0043"] / c.an$Df.Pooled[c.an$ID == "02E.0043"]
    c.an$F.Value[c.an$ID == "02E.0043"]<- c.an$Mean.Sq[c.an$ID == "02E.0043"] / c.an$Mean.Sq.Pooled[c.an$ID == "02E.0043"]
    c.an$P.Value.Upper.Tail[c.an$ID == "02E.0043"]<- stats::pf(c.an$F.Value[c.an$ID == "02E.0043"], c.an$Df[c.an$ID == "02E.0043"], c.an$Df.Pooled[c.an$ID == "02E.0043"], lower.tail = F)

    # 	If varPxCC <= 0, varS(CC) <=, varT(P)xCC <= 0, varPxS(CC) <= 0, then F = (MSCC / MSpooled), Where MSpooled = (SSPxCC + SSS(CC) + SST(P)xCC + SSPxS(CC) + SST(P)xS(CC)) / (dfPxCC + dfS(CC) + dfT(P)xCC + dfPxS(CC) + dfT(P)xS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0044"
    c.an$Df.Formula[c.an$ID == "02E.0044"]<-"(c-1)"
    c.an$Calculated.From[c.an$ID == "02E.0044"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "02E.0044"] <- c.an$Sum.Sq[c.an$ID == "02E.0044"] / c.an$Df[c.an$ID == "02E.0044"]
    c.an$Source.of.Variation[c.an$ID == "02E.0044"]<-"Among Controls = C"
    c.an$Denominator[c.an$ID == "02E.0044"]<-"PxC+S(C)+T(P)xC+PxS(C)+T(P)xS(C)"
    c.an$Df.Pooled[c.an$ID == "02E.0044"]<- c.an$Df[c.an$ID == "03BE.0040"] + c.an$Df[c.an$ID == "05EG.0040"]  + c.an$Df[c.an$ID == "10BEC.0040"]  + c.an$Df[c.an$ID == "09BEG.0040"] + c.an$Df[c.an$ID == "11BEGC.0040"]
    c.an$Sum.Sq.Pooled[c.an$ID == "02E.0044"]<- c.an$Sum.Sq[c.an$ID == "03BE.0040"] + c.an$Sum.Sq[c.an$ID == "05EG.0040"]  + c.an$Sum.Sq[c.an$ID == "10BEC.0040"]  + c.an$Sum.Sq[c.an$ID == "09BEG.0040"] + c.an$Sum.Sq[c.an$ID == "11BEGC.0040"]
    c.an$Mean.Sq.Pooled[c.an$ID == "02E.0044"]<- c.an$Sum.Sq.Pooled[c.an$ID == "02E.0044"] / c.an$Df.Pooled[c.an$ID == "02E.0044"]
    c.an$F.Value[c.an$ID == "02E.0044"]<- c.an$Mean.Sq[c.an$ID == "02E.0044"] / c.an$Mean.Sq.Pooled[c.an$ID == "02E.0044"]
    c.an$P.Value.Upper.Tail[c.an$ID == "02E.0044"]<- stats::pf(c.an$F.Value[c.an$ID == "02E.0044"], c.an$Df[c.an$ID == "02E.0044"], c.an$Df.Pooled[c.an$ID == "02E.0044"], lower.tail = F)

    # 	F = (MSIC /MSCC), or
    dii.02E.0030<-cbind(Source.of.Variation = "E",
                        sym.F[sym.F$Source.of.Variation == "E",2:3] - sym.C[sym.C$Source.of.Variation == "E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0030"
    i.an$Df.Formula[i.an$ID == "02E.0030"]<-"1"
    i.an$Calculated.From[i.an$ID == "02E.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "02E.0030"]<-i.an$Sum.Sq[i.an$ID == "02E.0030"] / i.an$Df[i.an$ID == "02E.0030"]
    i.an$Source.of.Variation[i.an$ID == "02E.0030"]<-"Impact vs Controls = I"
    i.an$Denominator[i.an$ID == "02E.0030"]<-"C"
    i.an$F.Value[i.an$ID == "02E.0030"]<-i.an$Mean.Sq[i.an$ID == "02E.0030"] / c.an$Mean.Sq[c.an$ID == "02E.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], c.an$Df[c.an$ID == "02E.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], c.an$Df[c.an$ID == "02E.0040"], lower.tail = T)

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	If varT(P)xS(CC) <= 0, then F = (MST(P)xS(IC) / MSpooled), where MSpooled = (SST(P)xS(CC) + SSRes) / (dfT(P)xS(CC) + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.11BEGC.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E:G:C"]<-"11BEGC.0031"
      i.an$Df.Formula[i.an$ID == "11BEGC.0031"]<-"p(t-1)(s-1)"
      i.an$Calculated.From[i.an$ID == "11BEGC.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "11BEGC.0031"]<-i.an$Sum.Sq[i.an$ID == "11BEGC.0031"] / i.an$Df[i.an$ID == "11BEGC.0031"]
      i.an$Source.of.Variation[i.an$ID == "11BEGC.0031"]<-"T(P) x S(I)"
      i.an$Denominator[i.an$ID == "11BEGC.0031"]<-"T(P)xS(C)+Residuals"
      i.an$Df.Pooled[i.an$ID == "11BEGC.0031"]<-c.an$Df[c.an$ID == "11BEGC.0040"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "11BEGC.0031"]<-c.an$Sum.Sq[c.an$ID == "11BEGC.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "11BEGC.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "11BEGC.0031"] / i.an$Df.Pooled[i.an$ID == "11BEGC.0031"]
      i.an$F.Value[i.an$ID == "11BEGC.0031"]<- i.an$Mean.Sq[i.an$ID == "11BEGC.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "11BEGC.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "11BEGC.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "11BEGC.0031"], i.an$Df[i.an$ID == "11BEGC.0031"], i.an$Df.Pooled[i.an$ID == "11BEGC.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "11BEGC.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "11BEGC.0031"], i.an$Df[i.an$ID == "11BEGC.0031"], i.an$Df.Pooled[i.an$ID == "11BEGC.0031"], lower.tail = T)

      # 	If varT(P)xCC <= 0 and varT(P)xS(IC) > 0, then F = (MST(P)xIC / MSpooled), where MSpooled = (SST(P)xCC + SST(P)xS(IC)) / (dfT(P)xCC) + dfT(P)xS(IC)) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.10BEC.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E:C"]<-"10BEC.0031"
      i.an$Df.Formula[i.an$ID == "10BEC.0031"]<-"p(t-1)"
      i.an$Calculated.From[i.an$ID == "10BEC.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "10BEC.0031"]<-i.an$Sum.Sq[i.an$ID == "10BEC.0031"] / i.an$Df[i.an$ID == "10BEC.0031"]
      i.an$Source.of.Variation[i.an$ID == "10BEC.0031"]<-"T(P) x I"
      i.an$Denominator[i.an$ID == "10BEC.0031"]<-"T(P)xC+T(P)xS(I)"
      i.an$Df.Pooled[i.an$ID == "10BEC.0031"]<-c.an$Df[c.an$ID == "10BEC.0040"] + i.an$Df[i.an$ID == "11BEGC.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "10BEC.0031"]<-c.an$Sum.Sq[c.an$ID == "10BEC.0040"] + i.an$Sum.Sq[i.an$ID == "11BEGC.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "10BEC.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "10BEC.0031"] / i.an$Df.Pooled[i.an$ID == "10BEC.0031"]
      i.an$F.Value[i.an$ID == "10BEC.0031"]<-i.an$Mean.Sq[i.an$ID == "10BEC.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "10BEC.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "10BEC.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "10BEC.0031"], i.an$Df[i.an$ID == "10BEC.0031"], i.an$Df.Pooled[i.an$ID == "10BEC.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "10BEC.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "10BEC.0031"], i.an$Df[i.an$ID == "10BEC.0031"], i.an$Df.Pooled[i.an$ID == "10BEC.0031"], lower.tail = T)

      # 	If varT(P)xCC <= 0 and varT(P)xS(IC) <= 0, then F = (MST(P)xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.10BEC.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E:C"]<-"10BEC.0032"
      i.an$Df.Formula[i.an$ID == "10BEC.0032"]<-"p(t-1)"
      i.an$Calculated.From[i.an$ID == "10BEC.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "10BEC.0032"]<-i.an$Sum.Sq[i.an$ID == "10BEC.0032"] / i.an$Df[i.an$ID == "10BEC.0032"]
      i.an$Source.of.Variation[i.an$ID == "10BEC.0032"]<-"T(P) x I"
      i.an$Denominator[i.an$ID == "10BEC.0032"]<-"T(P)xC+T(P)xS(I)+Residuals"
      i.an$Df.Pooled[i.an$ID == "10BEC.0032"]<- i.an$Df.Pooled[i.an$ID == "10BEC.0031"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "10BEC.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "10BEC.0031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "10BEC.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "10BEC.0032"] / i.an$Df.Pooled[i.an$ID == "10BEC.0032"]
      i.an$F.Value[i.an$ID == "10BEC.0032"]<-i.an$Mean.Sq[i.an$ID == "10BEC.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "10BEC.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "10BEC.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "10BEC.0032"], i.an$Df[i.an$ID == "10BEC.0032"], i.an$Df.Pooled[i.an$ID == "10BEC.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "10BEC.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "10BEC.0032"], i.an$Df[i.an$ID == "10BEC.0032"], i.an$Df.Pooled[i.an$ID == "10BEC.0032"], lower.tail = T)

      # 	If varPxS(CC) <= 0 and varT(P)xS(IC) > 0, then F = (MSPxS(IC) / MSpooled), where MSpooled = (SSPxS(CC) + SST(P)xS(IC)) / (dfPxS(CC) + dfT(P)xS(IC)) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.09BEG.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E:G"]<-"09BEG.0031"
      i.an$Df.Formula[i.an$ID == "09BEG.0031"]<-"(p-1)(s-1)"
      i.an$Calculated.From[i.an$ID == "09BEG.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "09BEG.0031"]<-i.an$Sum.Sq[i.an$ID == "09BEG.0031"] / i.an$Df[i.an$ID == "09BEG.0031"]
      i.an$Source.of.Variation[i.an$ID == "09BEG.0031"]<-"P x S(I)"
      i.an$Denominator[i.an$ID == "09BEG.0031"]<-"PxS(C)+T(P)xS(I)"
      i.an$Df.Pooled[i.an$ID == "09BEG.0031"]<-c.an$Df[c.an$ID == "09BEG.0040"] + i.an$Df[i.an$ID == "11BEGC.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "09BEG.0031"]<-c.an$Sum.Sq[c.an$ID == "09BEG.0040"] + i.an$Sum.Sq[i.an$ID == "11BEGC.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "09BEG.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "09BEG.0031"] / i.an$Df.Pooled[i.an$ID == "09BEG.0031"]
      i.an$F.Value[i.an$ID == "09BEG.0031"]<-i.an$Mean.Sq[i.an$ID == "09BEG.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "09BEG.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "09BEG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "09BEG.0031"], i.an$Df[i.an$ID == "09BEG.0031"], i.an$Df.Pooled[i.an$ID == "09BEG.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "09BEG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "09BEG.0031"], i.an$Df[i.an$ID == "09BEG.0031"], i.an$Df.Pooled[i.an$ID == "09BEG.0031"], lower.tail = T)

      # 	If varPxS(CC) <= 0 and varT(P)xS(IC) <= 0, then F = (MSPxS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.09BEG.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E:G"]<-"09BEG.0032"
      i.an$Df.Formula[i.an$ID == "09BEG.0032"]<-"(p-1)(s-1)"
      i.an$Calculated.From[i.an$ID == "09BEG.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "09BEG.0032"]<-i.an$Sum.Sq[i.an$ID == "09BEG.0032"] / i.an$Df[i.an$ID == "09BEG.0032"]
      i.an$Source.of.Variation[i.an$ID == "09BEG.0032"]<-"P x S(I)"
      i.an$Denominator[i.an$ID == "09BEG.0032"]<-"PxS(C)+T(P)xS(I)+Residuals"
      i.an$Df.Pooled[i.an$ID == "09BEG.0032"]<- i.an$Df.Pooled[i.an$ID == "09BEG.0031"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "09BEG.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "09BEG.0031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "09BEG.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "09BEG.0032"] / i.an$Df.Pooled[i.an$ID == "09BEG.0032"]
      i.an$F.Value[i.an$ID == "09BEG.0032"]<-i.an$Mean.Sq[i.an$ID == "09BEG.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "09BEG.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "09BEG.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "09BEG.0032"], i.an$Df[i.an$ID == "09BEG.0032"], i.an$Df.Pooled[i.an$ID == "09BEG.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "09BEG.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "09BEG.0032"], i.an$Df[i.an$ID == "09BEG.0032"], i.an$Df.Pooled[i.an$ID == "09BEG.0032"], lower.tail = T)

      # 	If varS(CC) <= 0 and varPxS(IC) > 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SSS(CC) + SSPxS(IC)) / (dfS(CC) + dfPxS(IC)) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.05EG.0030)
      i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0031"
      i.an$Df.Formula[i.an$ID == "05EG.0031"]<-"(s-1)"
      i.an$Calculated.From[i.an$ID == "05EG.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "05EG.0031"]<-i.an$Sum.Sq[i.an$ID == "05EG.0031"] / i.an$Df[i.an$ID == "05EG.0031"]
      i.an$Source.of.Variation[i.an$ID == "05EG.0031"]<-"S(I)"
      i.an$Denominator[i.an$ID == "05EG.0031"]<-"S(C)+PxS(I)"
      i.an$Df.Pooled[i.an$ID == "05EG.0031"]<- c.an$Df[c.an$ID == "05EG.0040"] + i.an$Df[i.an$ID == "09BEG.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0031"]<- c.an$Sum.Sq[c.an$ID == "05EG.0040"] + i.an$Sum.Sq[i.an$ID == "09BEG.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0031"] / i.an$Df.Pooled[i.an$ID == "05EG.0031"]
      i.an$F.Value[i.an$ID == "05EG.0031"]<-i.an$Mean.Sq[i.an$ID == "05EG.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "05EG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0031"], i.an$Df[i.an$ID == "05EG.0031"], i.an$Df.Pooled[i.an$ID == "05EG.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "05EG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0031"], i.an$Df[i.an$ID == "05EG.0031"], i.an$Df.Pooled[i.an$ID == "05EG.0031"], lower.tail = T)

      # 	If varS(CC) <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) > 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS1 + SST(P)xS(IC)) / (df1 + dfT(P)xS(IC)) = SS2/df2
      i.an<-dplyr::bind_rows(i.an,dii.05EG.0030)
      i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0032"
      i.an$Df.Formula[i.an$ID == "05EG.0032"]<-"(s-1)"
      i.an$Calculated.From[i.an$ID == "05EG.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "05EG.0032"]<-i.an$Sum.Sq[i.an$ID == "05EG.0032"] / i.an$Df[i.an$ID == "05EG.0032"]
      i.an$Source.of.Variation[i.an$ID == "05EG.0032"]<-"S(I)"
      i.an$Denominator[i.an$ID == "05EG.0032"]<-"S(C)+PxS(I)+T(P)xS(I)"
      i.an$Df.Pooled[i.an$ID == "05EG.0032"]<- i.an$Df.Pooled[i.an$ID == "05EG.0031"] + i.an$Df[i.an$ID == "11BEGC.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0031"] + i.an$Sum.Sq[i.an$ID == "11BEGC.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0032"] / i.an$Df.Pooled[i.an$ID == "05EG.0032"]
      i.an$F.Value[i.an$ID == "05EG.0032"]<-i.an$Mean.Sq[i.an$ID == "05EG.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "05EG.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0032"], i.an$Df[i.an$ID == "05EG.0032"], i.an$Df.Pooled[i.an$ID == "05EG.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "05EG.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0032"], i.an$Df[i.an$ID == "05EG.0032"], i.an$Df.Pooled[i.an$ID == "05EG.0032"], lower.tail = T)

      # 	If varS(CC) <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS2 + SSRes) / (df2 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.05EG.0030)
      i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0033"
      i.an$Df.Formula[i.an$ID == "05EG.0033"]<-"(s-1)"
      i.an$Calculated.From[i.an$ID == "05EG.0033"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "05EG.0033"]<-i.an$Sum.Sq[i.an$ID == "05EG.0033"] / i.an$Df[i.an$ID == "05EG.0033"]
      i.an$Source.of.Variation[i.an$ID == "05EG.0033"]<-"S(I)"
      i.an$Denominator[i.an$ID == "05EG.0033"]<-"S(C)+PxS(I)+T(P)xS(I)+Residuals"
      i.an$Df.Pooled[i.an$ID == "05EG.0033"]<- i.an$Df.Pooled[i.an$ID == "05EG.0032"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0033"]<- i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0032"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0033"] <- i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0033"] / i.an$Df.Pooled[i.an$ID == "05EG.0033"]
      i.an$F.Value[i.an$ID == "05EG.0033"]<-i.an$Mean.Sq[i.an$ID == "05EG.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0033"]
      i.an$P.Value.Upper.Tail[i.an$ID == "05EG.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0033"], i.an$Df[i.an$ID == "05EG.0033"], i.an$Df.Pooled[i.an$ID == "05EG.0033"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "05EG.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0033"], i.an$Df[i.an$ID == "05EG.0033"], i.an$Df.Pooled[i.an$ID == "05EG.0033"], lower.tail = T)

      # 	If varPxCC <= 0, varPxS(IC) <= 0 and varT(P)xIC > 0), then F = (MSPxIC / MSpooled), where MSpooled = (SSPxCC + SST(P)xIC) / (dfPxCC) + dfT(P)xIC)
      i.an<-dplyr::bind_rows(i.an,dii.03BE.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E"]<-"03BE.0031"
      i.an$Df.Formula[i.an$ID == "03BE.0031"]<-"(p-1)"
      i.an$Calculated.From[i.an$ID == "03BE.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "03BE.0031"]<-i.an$Sum.Sq[i.an$ID == "03BE.0031"] / i.an$Df[i.an$ID == "03BE.0031"]
      i.an$Source.of.Variation[i.an$ID == "03BE.0031"]<-"P x I"
      i.an$Denominator[i.an$ID == "03BE.0031"]<-"PxC+T(P)xI"
      i.an$Df.Pooled[i.an$ID == "03BE.0031"]<-c.an$Df[c.an$ID == "03BE.0040"] + i.an$Df[i.an$ID == "10BEC.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0031"]<-c.an$Sum.Sq[c.an$ID == "03BE.0040"] + i.an$Sum.Sq[i.an$ID == "10BEC.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0031"] / i.an$Df.Pooled[i.an$ID == "03BE.0031"]
      i.an$F.Value[i.an$ID == "03BE.0031"]<-i.an$Mean.Sq[i.an$ID == "03BE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "03BE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0031"], i.an$Df[i.an$ID == "03BE.0031"], i.an$Df.Pooled[i.an$ID == "03BE.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "03BE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0031"], i.an$Df[i.an$ID == "03BE.0031"], i.an$Df.Pooled[i.an$ID == "03BE.0031"], lower.tail = T)

      # 	If varPxCC <= 0, varT(P)xIC <= 0 and varPxS(IC) > 0, then F = (MSPxIC / MSpooled), where MSpooled = (SSPxCC + SSPxS(IC)) / (dfPxCC) + dfPxS(IC))
      i.an<-dplyr::bind_rows(i.an,dii.03BE.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E"]<-"03BE.0032"
      i.an$Df.Formula[i.an$ID == "03BE.0032"]<-"(p-1)"
      i.an$Calculated.From[i.an$ID == "03BE.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "03BE.0032"]<-i.an$Sum.Sq[i.an$ID == "03BE.0032"] / i.an$Df[i.an$ID == "03BE.0032"]
      i.an$Source.of.Variation[i.an$ID == "03BE.0032"]<-"P x I"
      i.an$Denominator[i.an$ID == "03BE.0032"]<-"PxC+PxS(I)"
      i.an$Df.Pooled[i.an$ID == "03BE.0032"]<-c.an$Df[c.an$ID == "03BE.0040"] + i.an$Df[i.an$ID == "09BEG.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0032"]<-c.an$Sum.Sq[c.an$ID == "03BE.0040"] + i.an$Sum.Sq[i.an$ID == "09BEG.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0032"] / i.an$Df.Pooled[i.an$ID == "03BE.0032"]
      i.an$F.Value[i.an$ID == "03BE.0032"]<-i.an$Mean.Sq[i.an$ID == "03BE.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "03BE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0032"], i.an$Df[i.an$ID == "03BE.0032"], i.an$Df.Pooled[i.an$ID == "03BE.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "03BE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0032"], i.an$Df[i.an$ID == "03BE.0032"], i.an$Df.Pooled[i.an$ID == "03BE.0032"], lower.tail = T)

      # 	If varPxCC <= 0, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) > 0, then F = (MSPxIC / MSpooled), where MSpooled = (SSPxCC + SST(P)xIC + SSPxS(IC) + SST(P)xS(IC)) / (dfTxCC) + dfT(P)xIC + dfPxS(IC) + dfT(P)xS(IC)) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.03BE.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E"]<-"03BE.0033"
      i.an$Df.Formula[i.an$ID == "03BE.0033"]<-"(p-1)"
      i.an$Calculated.From[i.an$ID == "03BE.0033"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "03BE.0033"]<-i.an$Sum.Sq[i.an$ID == "03BE.0033"] / i.an$Df[i.an$ID == "03BE.0033"]
      i.an$Source.of.Variation[i.an$ID == "03BE.0033"]<-"P x I"
      i.an$Denominator[i.an$ID == "03BE.0033"]<-"PxC+T(P)xI+PxS(I)+T(P)xS(I)"
      i.an$Df.Pooled[i.an$ID == "03BE.0033"]<-c.an$Df[c.an$ID == "03BE.0040"] + i.an$Df[i.an$ID == "10BEC.0030"] + i.an$Df[i.an$ID == "09BEG.0030"] + i.an$Df[i.an$ID == "11BEGC.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0033"]<-c.an$Sum.Sq[c.an$ID == "03BE.0040"] + i.an$Sum.Sq[i.an$ID == "09BEG.0030"] + i.an$Sum.Sq[i.an$ID == "09BEG.0030"] + i.an$Sum.Sq[i.an$ID == "11BEGC.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0033"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0033"] / i.an$Df.Pooled[i.an$ID == "03BE.0033"]
      i.an$F.Value[i.an$ID == "03BE.0033"]<-i.an$Mean.Sq[i.an$ID == "03BE.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0033"]
      i.an$P.Value.Upper.Tail[i.an$ID == "03BE.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0033"], i.an$Df[i.an$ID == "03BE.0033"], i.an$Df.Pooled[i.an$ID == "03BE.0033"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "03BE.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0033"], i.an$Df[i.an$ID == "03BE.0033"], i.an$Df.Pooled[i.an$ID == "03BE.0033"], lower.tail = T)

      # 	If varPxCC <= 0, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) <= 0, then F = (MSPxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.03BE.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E"]<-"03BE.0034"
      i.an$Df.Formula[i.an$ID == "03BE.0034"]<-"(p-1)"
      i.an$Calculated.From[i.an$ID == "03BE.0034"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "03BE.0034"]<-i.an$Sum.Sq[i.an$ID == "03BE.0031"] / i.an$Df[i.an$ID == "03BE.0031"]
      i.an$Source.of.Variation[i.an$ID == "03BE.0034"]<-"P x I"
      i.an$Denominator[i.an$ID == "03BE.0034"]<-"PxC+T(P)xI+PxS(I)+T(P)xS(I)+Residuals)"
      i.an$Df.Pooled[i.an$ID == "03BE.0034"]<- i.an$Df.Pooled[i.an$ID == "03BE.0033"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0034"]<- i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0033"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0034"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0034"] / i.an$Df.Pooled[i.an$ID == "03BE.0034"]
      i.an$F.Value[i.an$ID == "03BE.0034"]<-i.an$Mean.Sq[i.an$ID == "03BE.0034"] / i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0034"]
      i.an$P.Value.Upper.Tail[i.an$ID == "03BE.0034"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0034"], i.an$Df[i.an$ID == "03BE.0034"], i.an$Df.Pooled[i.an$ID == "03BE.0034"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "03BE.0034"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0034"], i.an$Df[i.an$ID == "03BE.0034"], i.an$Df.Pooled[i.an$ID == "03BE.0034"], lower.tail = T)

      # 	If varCC <= 0, varS(IC) <= 0 and varPxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSPxIC) / (dfCC + dfPxIC)
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0031"
      i.an$Df.Formula[i.an$ID == "02E.0031"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0031"]<-i.an$Sum.Sq[i.an$ID == "02E.0031"] / i.an$Df[i.an$ID == "02E.0031"]
      i.an$Source.of.Variation[i.an$ID == "02E.0031"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0031"]<-"C+PxI"
      i.an$Df.Pooled[i.an$ID == "02E.0031"]<-c.an$Df[c.an$ID == "02E.0040"] + i.an$Df[i.an$ID == "03BE.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"]<-c.an$Sum.Sq[c.an$ID == "02E.0040"] + i.an$Sum.Sq[i.an$ID == "03BE.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"] / i.an$Df.Pooled[i.an$ID == "02E.0031"]
      i.an$F.Value[i.an$ID == "02E.0031"]<-i.an$Mean.Sq[i.an$ID == "02E.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df.Pooled[i.an$ID == "02E.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df.Pooled[i.an$ID == "02E.0031"], lower.tail = T)

      # 	If varCC <= 0, varPxIC <= 0 and varS(IC) > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSS(IC)) / (dfCC + dfS(IC))
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0032"
      i.an$Df.Formula[i.an$ID == "02E.0032"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0032"]<-i.an$Sum.Sq[i.an$ID == "02E.0032"] / i.an$Df[i.an$ID == "02E.0032"]
      i.an$Source.of.Variation[i.an$ID == "02E.0032"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0032"]<-"C+S(I)"
      i.an$Df.Pooled[i.an$ID == "02E.0032"]<-c.an$Df[c.an$ID == "02E.0040"] + i.an$Df[i.an$ID == "05EG.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"]<-c.an$Sum.Sq[c.an$ID == "02E.0040"] + i.an$Sum.Sq[i.an$ID == "05EG.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"] / i.an$Df.Pooled[i.an$ID == "02E.0032"]
      i.an$F.Value[i.an$ID == "02E.0032"]<-i.an$Mean.Sq[i.an$ID == "02E.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = T)

      # 	If varCC <= 0, varPxIC <= 0, varS(IC) <=, varPxS(IC) <= 0 and varT(P)xIC > 0, then F = (MSIC / MSpooled), Where MSpooled = (SSCC + SSPxIC + SSS(IC) + SST(P)xIC) / (dfCC + dfPxIC + dfS(IC) + dfT(P)xIC)
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0033"
      i.an$Df.Formula[i.an$ID == "02E.0033"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0033"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0033"]<-i.an$Sum.Sq[i.an$ID == "02E.0033"] / i.an$Df[i.an$ID == "02E.0033"]
      i.an$Source.of.Variation[i.an$ID == "02E.0033"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0033"]<-"C+PxI+S(I)+T(P)xI"
      i.an$Df.Pooled[i.an$ID == "02E.0033"]<-c.an$Df[c.an$ID == "02E.0040"] + i.an$Df[i.an$ID == "03BE.0030"] + i.an$Df[i.an$ID == "05EG.0030"] + i.an$Df[i.an$ID == "10BEC.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0033"]<-c.an$Sum.Sq[c.an$ID == "02E.0040"] + i.an$Sum.Sq[i.an$ID == "03BE.0030"] + i.an$Sum.Sq[i.an$ID == "05EG.0030"] + i.an$Sum.Sq[i.an$ID == "10BEC.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0033"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0033"] / i.an$Df.Pooled[i.an$ID == "02E.0033"]
      i.an$F.Value[i.an$ID == "02E.0033"]<-i.an$Mean.Sq[i.an$ID == "02E.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0033"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0033"], i.an$Df[i.an$ID == "02E.0033"], i.an$Df.Pooled[i.an$ID == "02E.0033"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0033"], i.an$Df[i.an$ID == "02E.0033"], i.an$Df.Pooled[i.an$ID == "02E.0033"], lower.tail = T)

      # 	If varCC <= 0, varPxIC <= 0, varS(IC) <=, varT(P)xIC <= 0 and varPxS(IC) > 0, then F = (MSIC / MSpooled), Where MSpooled = (SSCC + SSPxIC + SSS(IC) + SSPxS(IC)) / (dfCC + dfPxIC + dfS(IC) + dfPxS(IC))
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0034"
      i.an$Df.Formula[i.an$ID == "02E.0034"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0034"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0034"]<-i.an$Sum.Sq[i.an$ID == "02E.0034"] / i.an$Df[i.an$ID == "02E.0034"]
      i.an$Source.of.Variation[i.an$ID == "02E.0034"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0034"]<-"C+PxI+S(I)+PxS(I)"
      i.an$Df.Pooled[i.an$ID == "02E.0034"]<-c.an$Df[c.an$ID == "02E.0040"] + i.an$Df[i.an$ID == "03BE.0030"] + i.an$Df[i.an$ID == "05EG.0030"] + i.an$Df[i.an$ID == "09BEG.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0034"]<-c.an$Sum.Sq[c.an$ID == "02E.0040"] + i.an$Sum.Sq[i.an$ID == "03BE.0030"] + i.an$Sum.Sq[i.an$ID == "05EG.0030"] + i.an$Sum.Sq[i.an$ID == "09BEG.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0034"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0034"] / i.an$Df.Pooled[i.an$ID == "02E.0034"]
      i.an$F.Value[i.an$ID == "02E.0034"]<-i.an$Mean.Sq[i.an$ID == "02E.0034"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0034"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0034"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0034"], i.an$Df[i.an$ID == "02E.0034"], i.an$Df.Pooled[i.an$ID == "02E.0034"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0034"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0034"], i.an$Df[i.an$ID == "02E.0034"], i.an$Df.Pooled[i.an$ID == "02E.0034"], lower.tail = T)

      # 	If varCC <= 0, varPxIC <= 0, varS(IC) <=, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) > 0, then F = (MSIC / MSpooled), Where MSpooled = (SSCC + SSPxIC + SSS(IC) + SST(P)xIC + SSPxS(IC) + SST(P)xS(IC)) / (dfCC + dfPxIC + dfS(IC) + dfT(P)xIC + dfPxS(IC) + dfT(P)xS(IC)) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0035"
      i.an$Df.Formula[i.an$ID == "02E.0035"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0035"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0035"]<-i.an$Sum.Sq[i.an$ID == "02E.0035"] / i.an$Df[i.an$ID == "02E.0035"]
      i.an$Source.of.Variation[i.an$ID == "02E.0035"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0035"]<-"C+PxI+S(I)+PxS(I)+T(P)xI+T(P)xS(I)"
      i.an$Df.Pooled[i.an$ID == "02E.0035"]<- c.an$Df[c.an$ID == "02E.0040"] + i.an$Df[i.an$ID == "03BE.0030"] + i.an$Df[i.an$ID == "05EG.0030"] + i.an$Df[i.an$ID == "09BEG.0030"] + i.an$Df[i.an$ID == "10BEC.0030"] + i.an$Df[i.an$ID == "11BEGC.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0035"]<- c.an$Sum.Sq[c.an$ID == "02E.0040"] + i.an$Sum.Sq[i.an$ID == "03BE.0030"] + i.an$Sum.Sq[i.an$ID == "05EG.0030"] + i.an$Sum.Sq[i.an$ID == "09BEG.0030"] + i.an$Sum.Sq[i.an$ID == "10BEC.0030"] + i.an$Df[i.an$ID == "11BEGC.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0035"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0035"] / i.an$Df.Pooled[i.an$ID == "02E.0035"]
      i.an$F.Value[i.an$ID == "02E.0035"]<-i.an$Mean.Sq[i.an$ID == "02E.0035"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0035"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0035"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0035"], i.an$Df[i.an$ID == "02E.0035"], i.an$Df.Pooled[i.an$ID == "02E.0035"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0035"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0035"], i.an$Df[i.an$ID == "02E.0035"], i.an$Df.Pooled[i.an$ID == "02E.0035"], lower.tail = T)

      # 	If varCC <= 0, varPxIC <= 0, varS(IC) <= 0, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) <= 0, then F = (MSIC / MSpooled), Where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0036"
      i.an$Df.Formula[i.an$ID == "02E.0036"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0036"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0036"]<-i.an$Sum.Sq[i.an$ID == "02E.0036"] / i.an$Df[i.an$ID == "02E.0036"]
      i.an$Source.of.Variation[i.an$ID == "02E.0036"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0036"]<-"C+PxI+S(I)+PxS(I)+T(P)xI+T(P)xS(I)+Residuals"
      i.an$Df.Pooled[i.an$ID == "02E.0036"]<-  i.an$Df.Pooled[i.an$ID == "02E.0035"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0036"]<- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0035"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0036"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0036"] / i.an$Df.Pooled[i.an$ID == "02E.0036"]
      i.an$F.Value[i.an$ID == "02E.0036"]<-i.an$Mean.Sq[i.an$ID == "02E.0036"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0036"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0036"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0036"], i.an$Df[i.an$ID == "02E.0036"], i.an$Df.Pooled[i.an$ID == "02E.0036"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0036"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0036"], i.an$Df[i.an$ID == "02E.0036"], i.an$Df.Pooled[i.an$ID == "02E.0036"], lower.tail = T)
    }
  }


  ###------------------ INTERACTIONS with TWO temporal and ONE spatial factors

  if(aav.model == "BEC") {
    # 	F = (MST(P)xCC / MSRes)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "B:E:C",1:3])
    c.an$ID[c.an$Source.of.Variation == "B:E:C"]<-"10BEC.0040"
    c.an$Df.Formula[c.an$ID == "10BEC.0040"]<-"p(t-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "10BEC.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "10BEC.0040"] <- c.an$Sum.Sq[c.an$ID == "10BEC.0040"] / c.an$Df[c.an$ID == "10BEC.0040"]
    c.an$Source.of.Variation[c.an$ID == "10BEC.0040"]<-"T(P) x C"
    c.an$Denominator[c.an$ID == "10BEC.0040"]<-"Residuals"
    c.an$F.Value[c.an$ID == "10BEC.0040"]<-c.an$Mean.Sq[c.an$ID == "10BEC.0040"] / f.an$Mean.Sq[f.an$ID == "12R.0000"]
    c.an$P.Value.Upper.Tail[c.an$ID == "10BEC.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "10BEC.0040"], c.an$Df[c.an$ID == "10BEC.0040"], f.an$Df[f.an$ID == "12R.0000"], lower.tail = F)

    # 	F = (MST(P)xIC) / MST(P)xCC), or
    dii.10BEC.0030<-cbind(Source.of.Variation = "B:E:C",
                          sym.F[sym.F$Source.of.Variation == "B:E:C",2:3] - sym.C[sym.C$Source.of.Variation == "B:E:C",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.10BEC.0030)
    i.an$ID[i.an$Source.of.Variation == "B:E:C"]<-"10BEC.0030"
    i.an$Df.Formula[i.an$ID == "10BEC.0030"]<-"p(t-1)"
    i.an$Calculated.From[i.an$ID == "10BEC.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "10BEC.0030"]<-i.an$Sum.Sq[i.an$ID == "10BEC.0030"] / i.an$Df[i.an$ID == "10BEC.0030"]
    i.an$Source.of.Variation[i.an$ID == "10BEC.0030"]<-"T(P) x I"
    i.an$Denominator[i.an$ID == "10BEC.0030"]<-"T(P)xC"
    i.an$F.Value[i.an$ID == "10BEC.0030"]<-i.an$Mean.Sq[i.an$ID == "10BEC.0030"] / c.an$Mean.Sq[c.an$ID == "10BEC.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "10BEC.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "10BEC.0030"], i.an$Df[i.an$ID == "10BEC.0030"], c.an$Df[c.an$ID == "10BEC.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "10BEC.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "10BEC.0030"], i.an$Df[i.an$ID == "10BEC.0030"], c.an$Df[c.an$ID == "10BEC.0040"], lower.tail = T)

    # 	F = (MS_PxCC / MS_T(P)xCC)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "B:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "B:E"]<-"03BE.0040"
    c.an$Df.Formula[c.an$ID == "03BE.0040"]<-"(p-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03BE.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03BE.0040"] <- c.an$Sum.Sq[c.an$ID == "03BE.0040"] / c.an$Df[c.an$ID == "03BE.0040"]
    c.an$Source.of.Variation[c.an$ID == "03BE.0040"]<-"P x C"
    c.an$Denominator[c.an$ID == "03BE.0040"]<-"Residuals"
    c.an$F.Value[c.an$ID == "03BE.0040"]<-c.an$Mean.Sq[c.an$ID == "03BE.0040"] / c.an$Mean.Sq[c.an$ID == "10BEC.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03BE.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "03BE.0040"], c.an$Df[c.an$ID == "03BE.0040"], c.an$Df[c.an$ID == "10BEC.0040"], lower.tail = F)

    # 	F = (MSPxIC) / MSPxCC), or
    dii.03BE.0030<-cbind(Source.of.Variation = "B:E",
                         sym.F[sym.F$Source.of.Variation == "B:E",2:3] - sym.C[sym.C$Source.of.Variation == "B:E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.03BE.0030)
    i.an$ID[i.an$Source.of.Variation == "B:E"]<-"03BE.0030"
    i.an$Df.Formula[i.an$ID == "03BE.0030"]<-"(p-1)"
    i.an$Calculated.From[i.an$ID == "03BE.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "03BE.0030"]<-i.an$Sum.Sq[i.an$ID == "03BE.0030"] / i.an$Df[i.an$ID == "03BE.0030"]
    i.an$Source.of.Variation[i.an$ID == "03BE.0030"]<-"P x I"
    i.an$Denominator[i.an$ID == "03BE.0030"]<-"PxC"
    i.an$F.Value[i.an$ID == "03BE.0030"]<-i.an$Mean.Sq[i.an$ID == "03BE.0030"] / c.an$Mean.Sq[c.an$ID == "03BE.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "03BE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0030"], i.an$Df[i.an$ID == "03BE.0030"], c.an$Df[c.an$ID == "03BE.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "03BE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0030"], i.an$Df[i.an$ID == "03BE.0030"], c.an$Df[c.an$ID == "03BE.0040"], lower.tail = T)

    # 	F = (MSCC / MSPxCC)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0040"
    c.an$Df.Formula[c.an$ID == "02E.0040"]<-"(c-1)"
    c.an$Calculated.From[c.an$ID == "02E.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "02E.0040"] <- c.an$Sum.Sq[c.an$ID == "02E.0040"] / c.an$Df[c.an$ID == "02E.0040"]
    c.an$Source.of.Variation[c.an$ID == "02E.0040"]<-"Among Controls = C"
    c.an$Denominator[c.an$ID == "02E.0040"]<-"PxC"
    c.an$F.Value[c.an$ID == "02E.0040"]<-c.an$Mean.Sq[c.an$ID == "02E.0040"] / c.an$Mean.Sq[c.an$ID == "03BE.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "02E.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "02E.0040"], c.an$Df[c.an$ID == "02E.0040"], c.an$Df[c.an$ID == "03BE.0040"], lower.tail = F)

    # 	F = (MSIC / MSCC), or
    dii.02E.0030<-cbind(Source.of.Variation = "E",
                        sym.F[sym.F$Source.of.Variation == "E",2:3] - sym.C[sym.C$Source.of.Variation == "E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0030"
    i.an$Df.Formula[i.an$ID == "02E.0030"]<-"1"
    i.an$Calculated.From[i.an$ID == "02E.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "02E.0030"]<-i.an$Sum.Sq[i.an$ID == "02E.0030"] / i.an$Df[i.an$ID == "02E.0030"]
    i.an$Source.of.Variation[i.an$ID == "02E.0030"]<-"Impact vs Controls = I"
    i.an$Denominator[i.an$ID == "02E.0030"]<-"C"
    i.an$F.Value[i.an$ID == "02E.0030"]<-i.an$Mean.Sq[i.an$ID == "02E.0030"] / c.an$Mean.Sq[c.an$ID == "02E.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], c.an$Df[c.an$ID == "02E.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], c.an$Df[c.an$ID == "02E.0040"], lower.tail = T)

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	If varT(P)xCC <= 0, then F = (MS_T(P)xIC / MSpooled), where MSpooled = (SS_T(P)xCC + SSRes) / (dfT(P)xCC) + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.10BEC.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E:C"]<-"10BEC.0031"
      i.an$Df.Formula[i.an$ID == "10BEC.0031"]<-"p(t-1)"
      i.an$Calculated.From[i.an$ID == "10BEC.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "10BEC.0031"]<-i.an$Sum.Sq[i.an$ID == "10BEC.0031"] / i.an$Df[i.an$ID == "10BEC.0031"]
      i.an$Source.of.Variation[i.an$ID == "10BEC.0031"]<-"T(P) x I"
      i.an$Denominator[i.an$ID == "10BEC.0031"]<-"T(P)xC+Residuals"
      i.an$Df.Pooled[i.an$ID == "10BEC.0031"]<-c.an$Df[c.an$ID == "10BEC.0040"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "10BEC.0031"]<-c.an$Sum.Sq[c.an$ID == "10BEC.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "10BEC.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "10BEC.0031"] / i.an$Df.Pooled[i.an$ID == "10BEC.0031"]
      i.an$F.Value[i.an$ID == "10BEC.0031"]<-i.an$Mean.Sq[i.an$ID == "10BEC.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "10BEC.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "10BEC.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "10BEC.0031"], i.an$Df[i.an$ID == "10BEC.0031"], i.an$Df.Pooled[i.an$ID == "10BEC.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "10BEC.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "10BEC.0031"], i.an$Df[i.an$ID == "10BEC.0031"], i.an$Df.Pooled[i.an$ID == "10BEC.0031"], lower.tail = T)

      # 	If varPxCC <= 0 and varT(P)xIC > 0, then F = (MSPxIC / MSpooled), where MSpooled = (SSPxCC + SST(P)xIC) / (dfPxCC) + dfT(P)xIC) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.03BE.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E"]<-"03BE.0031"
      i.an$Df.Formula[i.an$ID == "03BE.0031"]<-"(p-1)"
      i.an$Calculated.From[i.an$ID == "03BE.0031"]<-"sym.F-sym.C"
      i.an$Source.of.Variation[i.an$ID == "03BE.0031"]<-"P x I"
      i.an$Mean.Sq[i.an$ID == "03BE.0031"]<-i.an$Sum.Sq[i.an$ID == "03BE.0031"] / i.an$Df[i.an$ID == "03BE.0031"]
      i.an$Denominator[i.an$ID == "03BE.0031"]<-"PxC+T(P)xI"
      i.an$Df.Pooled[i.an$ID == "03BE.0031"]<-c.an$Df[c.an$ID == "03BE.0040"] + i.an$Df[i.an$ID == "10BEC.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0031"]<-c.an$Sum.Sq[c.an$ID == "03BE.0040"] + i.an$Sum.Sq[i.an$ID == "10BEC.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0031"] / i.an$Df.Pooled[i.an$ID == "03BE.0031"]
      i.an$F.Value[i.an$ID == "03BE.0031"]<-i.an$Mean.Sq[i.an$ID == "03BE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "03BE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0031"], i.an$Df[i.an$ID == "03BE.0031"], i.an$Df.Pooled[i.an$ID == "03BE.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "03BE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0031"], i.an$Df[i.an$ID == "03BE.0031"], i.an$Df.Pooled[i.an$ID == "03BE.0031"], lower.tail = T)

      # 	If varPxCC <= 0 and varT(P)xIC <= 0, then F = (MSPxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.03BE.0030)
      i.an$ID[i.an$Source.of.Variation == "B:E"]<-"03BE.0032"
      i.an$Df.Formula[i.an$ID == "03BE.0032"]<-"(p-1)"
      i.an$Calculated.From[i.an$ID == "03BE.0032"]<-"sym.F-sym.C"
      i.an$Source.of.Variation[i.an$ID == "03BE.0032"]<-"P x I"
      i.an$Mean.Sq[i.an$ID == "03BE.0032"]<-i.an$Sum.Sq[i.an$ID == "03BE.0032"] / i.an$Df[i.an$ID == "03BE.0032"]
      i.an$Denominator[i.an$ID == "03BE.0032"]<-"PxC+T(P)xI+Residuals"
      i.an$Df.Pooled[i.an$ID == "03BE.0032"]<-i.an$Df.Pooled[i.an$ID == "03BE.0031"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0032"]<-i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03BE.0032"] / i.an$Df.Pooled[i.an$ID == "03BE.0032"]
      i.an$F.Value[i.an$ID == "03BE.0032"]<-i.an$Mean.Sq[i.an$ID == "03BE.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "03BE.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "03BE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0032"], i.an$Df[i.an$ID == "03BE.0032"], i.an$Df.Pooled[i.an$ID == "03BE.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "03BE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "03BE.0032"], i.an$Df[i.an$ID == "03BE.0032"], i.an$Df.Pooled[i.an$ID == "03BE.0032"], lower.tail = T)

      # 	If varCC <= 0 and varPxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSPxIC) / (dfCC + dfPxIC) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0031"
      i.an$Df.Formula[i.an$ID == "02E.0031"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0031"]<-i.an$Sum.Sq[i.an$ID == "02E.0031"] / i.an$Df[i.an$ID == "02E.0031"]
      i.an$Source.of.Variation[i.an$ID == "02E.0031"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0031"]<-"C+PxI"
      i.an$Df.Pooled[i.an$ID == "02E.0031"]<- c.an$Df[c.an$ID == "02E.0040"] + i.an$Df[i.an$ID == "03BE.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"]<- c.an$Sum.Sq[c.an$ID == "02E.0040"] + i.an$Sum.Sq[i.an$ID == "03BE.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"] / i.an$Df.Pooled[i.an$ID == "02E.0031"]
      i.an$F.Value[i.an$ID == "02E.0031"]<-i.an$Mean.Sq[i.an$ID == "02E.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df.Pooled[i.an$ID == "02E.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df.Pooled[i.an$ID == "02E.0031"], lower.tail = T)

      # 	If varCC <= 0 and varPxIC <= 0 and varT(P)xIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SS1 + SST(P)xIC) / (df1 + dfT(P)xIC) = SS2/df2
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0032"
      i.an$Df.Formula[i.an$ID == "02E.0032"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0032"]<-i.an$Sum.Sq[i.an$ID == "02E.0032"] / i.an$Df[i.an$ID == "02E.0032"]
      i.an$Source.of.Variation[i.an$ID == "02E.0032"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0032"]<-"C+PxI+T(P)xI"
      i.an$Df.Pooled[i.an$ID == "02E.0032"]<- i.an$Df.Pooled[i.an$ID == "02E.0031"] + i.an$Df[i.an$ID == "10BEC.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"] + i.an$Sum.Sq[i.an$ID == "10BEC.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"] / i.an$Df.Pooled[i.an$ID == "02E.0032"]
      i.an$F.Value[i.an$ID == "02E.0032"]<-i.an$Mean.Sq[i.an$ID == "02E.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = T)

      # 	If varCC <= 0 and varPxIC <= 0 and varT(P)xIC <= 0, then F = (MSIC / MSpooled), Where MSpooled = (SS2 + SSRes) / (df2 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0033"
      i.an$Df.Formula[i.an$ID == "02E.0033"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0033"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0033"]<-i.an$Sum.Sq[i.an$ID == "02E.0033"] / i.an$Df[i.an$ID == "02E.0033"]
      i.an$Source.of.Variation[i.an$ID == "02E.0033"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0033"]<-"C+PxI+T(P)xI+Residuals"
      i.an$Df.Pooled[i.an$ID == "02E.0033"]<- i.an$Df.Pooled[i.an$ID == "02E.0032"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0033"]<- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0033"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0033"] / i.an$Df.Pooled[i.an$ID == "02E.0033"]
      i.an$F.Value[i.an$ID == "02E.0033"]<-i.an$Mean.Sq[i.an$ID == "02E.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0033"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0033"], i.an$Df[i.an$ID == "02E.0033"], i.an$Df.Pooled[i.an$ID == "02E.0033"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0033"], i.an$Df[i.an$ID == "02E.0033"], i.an$Df.Pooled[i.an$ID == "02E.0033"], lower.tail = T)
    }
  }


  ###------------------ INTERACTIONS with ONE temporal and TWO spatial factors

  if(aav.model == "CEG") {
    # 	F = (MSTxS(CC)) / MSRes)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "C:E:G",1:3])
    c.an$ID[c.an$Source.of.Variation == "C:E:G"]<-"11CEG.0040"
    c.an$Df.Formula[c.an$ID == "11CEG.0040"]<-"(t-1)c(s-1)"
    c.an$Calculated.From[c.an$ID == "11CEG.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "11CEG.0040"] <- c.an$Sum.Sq[c.an$ID == "11CEG.0040"] / c.an$Df[c.an$ID == "11CEG.0040"]
    c.an$Source.of.Variation[c.an$ID == "11CEG.0040"]<-"T x S(C)"
    c.an$Denominator[c.an$ID == "11CEG.0040"]<-"Residuals"
    c.an$F.Value[c.an$ID == "11CEG.0040"]<-c.an$Mean.Sq[c.an$ID == "11CEG.0040"] / f.an$Mean.Sq[f.an$ID == "12R.0000"]
    c.an$P.Value.Upper.Tail[c.an$ID == "11CEG.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "11CEG.0040"], c.an$Df[c.an$ID == "11CEG.0040"], f.an$Df[f.an$ID == "12R.0000"], lower.tail = F)

    # 	F = (MSTxS(IC) / MSTxS(CC)), or
    dii.11CEG.0030<-cbind(Source.of.Variation = "C:E:G",
                          sym.F[sym.F$Source.of.Variation == "C:E:G",2:3] - sym.C[sym.C$Source.of.Variation == "C:E:G",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.11CEG.0030)
    i.an$ID[i.an$Source.of.Variation == "C:E:G"]<-"11CEG.0030"
    i.an$Df.Formula[i.an$ID == "11CEG.0030"]<-"(t-1)(s-1)"
    i.an$Mean.Sq[i.an$ID == "11CEG.0030"]<-i.an$Sum.Sq[i.an$ID == "11CEG.0030"] / i.an$Df[i.an$ID == "11CEG.0030"]
    i.an$Calculated.From[i.an$ID == "11CEG.0030"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$ID == "11CEG.0030"]<-"T x S(I)"
    i.an$Denominator[i.an$ID == "11CEG.0030"]<-"TxS(C)"
    i.an$F.Value[i.an$ID == "11CEG.0030"]<-i.an$Mean.Sq[i.an$ID == "11CEG.0030"] / c.an$Mean.Sq[c.an$ID == "11CEG.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "11CEG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "11CEG.0030"], i.an$Df[i.an$ID == "11CEG.0030"], c.an$Df[c.an$ID == "11CEG.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "11CEG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "11CEG.0030"], i.an$Df[i.an$ID == "11CEG.0030"], c.an$Df[c.an$ID == "11CEG.0040"], lower.tail = T)

    # 	F = (MSTxCC) / MSTxS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "C:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "C:E"]<-"10CE.0040"
    c.an$Df.Formula[c.an$ID == "10CE.0040"]<-"(t-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "10CE.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "10CE.0040"] <- c.an$Sum.Sq[c.an$ID == "10CE.0040"] / c.an$Df[c.an$ID == "10CE.0040"]
    c.an$Source.of.Variation[c.an$ID == "10CE.0040"]<-"T x C"
    c.an$Denominator[c.an$ID == "10CE.0040"]<-"TxS(C)"
    c.an$F.Value[c.an$ID == "10CE.0040"]<-c.an$Mean.Sq[c.an$ID == "10CE.0040"] / c.an$Mean.Sq[c.an$ID == "11CEG.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "10CE.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "10CE.0040"], c.an$Df[c.an$ID == "10CE.0040"], c.an$Df[c.an$ID == "11CEG.0040"], lower.tail = F)

    # 	F = (MSTxIC) / MSTxCC), or
    dii.10CE.0030<-cbind(Source.of.Variation = "C:E",
                         sym.F[sym.F$Source.of.Variation == "C:E",2:3] - sym.C[sym.C$Source.of.Variation == "C:E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.10CE.0030)
    i.an$ID[i.an$Source.of.Variation == "C:E"]<-"10CE.0030"
    i.an$Df.Formula[i.an$ID == "10CE.0030"]<-"(t-1)"
    i.an$Calculated.From[i.an$ID == "10CE.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "10CE.0030"]<-i.an$Sum.Sq[i.an$ID == "10CE.0030"] / i.an$Df[i.an$ID == "10CE.0030"]
    i.an$Source.of.Variation[i.an$ID == "10CE.0030"]<-"T x I"
    i.an$Denominator[i.an$ID == "10CE.0030"]<-"TxC"
    i.an$F.Value[i.an$ID == "10CE.0030"]<-i.an$Mean.Sq[i.an$ID == "10CE.0030"] / c.an$Mean.Sq[c.an$ID == "10CE.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "10CE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "10CE.0030"], i.an$Df[i.an$ID == "10CE.0030"], c.an$Df[c.an$ID == "10CE.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "10CE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "10CE.0030"], i.an$Df[i.an$ID == "10CE.0030"], c.an$Df[c.an$ID == "10CE.0040"], lower.tail = T)

    # 	F = (MSS(CC)) / MSTxS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E:G",1:3])
    c.an$ID[c.an$Source.of.Variation == "E:G"]<-"05EG.0040"
    c.an$Df.Formula[c.an$ID == "05EG.0040"]<-"c(s-1)"
    c.an$Calculated.From[c.an$ID == "05EG.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "05EG.0040"] <- c.an$Sum.Sq[c.an$ID == "05EG.0040"] / c.an$Df[c.an$ID == "05EG.0040"]
    c.an$Source.of.Variation[c.an$ID == "05EG.0040"]<-"S(C)"
    c.an$Denominator[c.an$ID == "05EG.0040"]<-"TxS(C)"
    c.an$F.Value[c.an$ID == "05EG.0040"]<-c.an$Mean.Sq[c.an$ID == "05EG.0040"] / c.an$Mean.Sq[c.an$ID == "11CEG.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "05EG.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "05EG.0040"], c.an$Df[c.an$ID == "05EG.0040"], c.an$Df[c.an$ID == "11CEG.0040"], lower.tail = F)

    # 	F = (MSS(IC) / MSS(CC)), or
    dii.05EG.0030<-cbind(Source.of.Variation = "E:G",
                         sym.F[sym.F$Source.of.Variation == "E:G",2:3] - sym.C[sym.C$Source.of.Variation == "E:G",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.05EG.0030)
    i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0030"
    i.an$Df.Formula[i.an$ID == "05EG.0030"]<-"(s-1)"
    i.an$Calculated.From[i.an$ID == "05EG.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "05EG.0030"]<-i.an$Sum.Sq[i.an$ID == "05EG.0030"] / i.an$Df[i.an$ID == "05EG.0030"]
    i.an$Source.of.Variation[i.an$ID == "05EG.0030"]<-"S(I)"
    i.an$Denominator[i.an$ID == "05EG.0030"]<-"S(C)"
    i.an$F.Value[i.an$ID == "05EG.0030"]<-i.an$Mean.Sq[i.an$ID == "05EG.0030"] / c.an$Mean.Sq[c.an$ID == "05EG.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "05EG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0030"], i.an$Df[i.an$ID == "05EG.0030"], c.an$Df[c.an$ID == "05EG.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "05EG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0030"], i.an$Df[i.an$ID == "05EG.0030"], c.an$Df[c.an$ID == "05EG.0040"], lower.tail = T)

    # 	If varS(CC) <= 0 and varTxCC > 0, then F = (MSCC / MSTxCC)  02E.0040
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0040"
    c.an$Df.Formula[c.an$ID == "02E.0040"]<-"(c-1)"
    c.an$Calculated.From[c.an$ID == "02E.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "02E.0040"] <- c.an$Sum.Sq[c.an$ID == "02E.0040"] / c.an$Df[c.an$ID == "02E.0040"]
    c.an$Source.of.Variation[c.an$ID == "02E.0040"]<-"Among Controls = C"
    c.an$Denominator[c.an$ID == "02E.0040"]<-"TxC"
    c.an$F.Value[c.an$ID == "02E.0040"]<-c.an$Mean.Sq[c.an$ID == "02E.0040"] / c.an$Mean.Sq[c.an$ID == "10CE.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "02E.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "02E.0040"], c.an$Df[c.an$ID == "02E.0040"], c.an$Df[c.an$ID == "10CE.0040"], lower.tail = F)

    # 	If varTxCC <= 0, varS(CC) > 0, then F = (MSCC / MSS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0041"
    c.an$Df.Formula[c.an$ID == "02E.0041"]<-"(c-1)"
    c.an$Calculated.From[c.an$ID == "02E.0041"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "02E.0041"] <- c.an$Sum.Sq[c.an$ID == "02E.0041"] / c.an$Df[c.an$ID == "02E.0041"]
    c.an$Source.of.Variation[c.an$ID == "02E.0041"]<-"Among Controls = C"
    c.an$Denominator[c.an$ID == "02E.0041"]<-"S(C)"
    c.an$F.Value[c.an$ID == "02E.0041"]<-c.an$Mean.Sq[c.an$ID == "02E.0041"] / c.an$Mean.Sq[c.an$ID == "05EG.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "02E.0041"]<- stats::pf(c.an$F.Value[c.an$ID == "02E.0041"], c.an$Df[c.an$ID == "02E.0041"], c.an$Df[c.an$ID == "05EG.0040"], lower.tail = F)

    # 	If varTxCC <= 0, varS(CC) <= 0, then F = (MSCC / MSpooled), Where MSpooled = (SSTxCC + SSS(CC) + SSTxS(CC)) / (dfTxCC + dfS(CC) + dfTxS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0042"
    c.an$Df.Formula[c.an$ID == "02E.0042"]<-"(c-1)"
    c.an$Calculated.From[c.an$ID == "02E.0042"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "02E.0042"] <- c.an$Sum.Sq[c.an$ID == "02E.0042"] / c.an$Df[c.an$ID == "02E.0042"]
    c.an$Source.of.Variation[c.an$ID == "02E.0042"]<-"Among Controls = C"
    c.an$Denominator[c.an$ID == "02E.0042"]<-"TxC+S(C)+TxS(C)"
    c.an$Df.Pooled[c.an$ID == "02E.0042"]<- c.an$Df[c.an$ID == "10CE.0040"] + c.an$Df[c.an$ID == "05EG.0040"] + c.an$Df[c.an$ID == "11CEG.0040"]
    c.an$Sum.Sq.Pooled[c.an$ID == "02E.0042"]<- c.an$Sum.Sq[c.an$ID == "10CE.0040"] + c.an$Sum.Sq[c.an$ID == "05EG.0040"] + c.an$Sum.Sq[c.an$ID == "11CEG.0040"]
    c.an$Mean.Sq.Pooled[c.an$ID == "02E.0042"]<- c.an$Sum.Sq.Pooled[c.an$ID == "02E.0042"] / c.an$Df.Pooled[c.an$ID == "02E.0042"]
    c.an$F.Value[c.an$ID == "02E.0042"]<- c.an$Mean.Sq[c.an$ID == "02E.0042"] / c.an$Mean.Sq.Pooled[c.an$ID == "02E.0042"]
    c.an$P.Value.Upper.Tail[c.an$ID == "02E.0042"]<- stats::pf(c.an$F.Value[c.an$ID == "02E.0042"], c.an$Df[c.an$ID == "02E.0042"], c.an$Df.Pooled[c.an$ID == "02E.0042"], lower.tail = F)

    # 	F = (MSIC / MSCC), or  02E.0030
    dii.02E.0030<-cbind(Source.of.Variation = "E",
                        sym.F[sym.F$Source.of.Variation == "E",2:3] - sym.C[sym.C$Source.of.Variation == "E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0030"
    i.an$Df.Formula[i.an$ID == "02E.0030"]<-"1"
    i.an$Calculated.From[i.an$ID == "02E.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "02E.0030"]<-i.an$Sum.Sq[i.an$ID == "02E.0030"] / i.an$Df[i.an$ID == "02E.0030"]
    i.an$Source.of.Variation[i.an$ID == "02E.0030"]<-"Impact vs Controls = I"
    i.an$Denominator[i.an$ID == "02E.0030"]<-"C"
    i.an$F.Value[i.an$ID == "02E.0030"]<-i.an$Mean.Sq[i.an$ID == "02E.0030"] / c.an$Mean.Sq[c.an$ID == "02E.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], c.an$Df[c.an$ID == "02E.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], c.an$Df[c.an$ID == "02E.0040"], lower.tail = T)

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	If varTxS(CC) <= 0, then F = (MSTxS(IC) / MSpooled), where MSpooled = (SSTxS(CC) + SSRes) / (dfTxS(CC) + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.11CEG.0030)
      i.an$ID[i.an$Source.of.Variation == "C:E:G"]<-"11CEG.0031"
      i.an$Df.Formula[i.an$ID == "11CEG.0031"]<-"(t-1)(s-1)"
      i.an$Mean.Sq[i.an$ID == "11CEG.0031"]<-i.an$Sum.Sq[i.an$ID == "11CEG.0031"] / i.an$Df[i.an$ID == "11CEG.0031"]
      i.an$Calculated.From[i.an$ID == "11CEG.0031"]<-"sym.F-sym.C"
      i.an$Source.of.Variation[i.an$ID == "11CEG.0031"]<-"T x S(I)"
      i.an$Denominator[i.an$ID == "11CEG.0031"]<-"TxS(C)+Residuals"
      i.an$Df.Pooled[i.an$ID == "11CEG.0031"]<-c.an$Df[c.an$ID == "11CEG.0040"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "11CEG.0031"]<-c.an$Sum.Sq[c.an$ID == "11CEG.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "11CEG.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "11CEG.0031"] / i.an$Df.Pooled[i.an$ID == "11CEG.0031"]
      i.an$F.Value[i.an$ID == "11CEG.0031"]<-i.an$Mean.Sq[i.an$ID == "11CEG.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "11CEG.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "11CEG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "11CEG.0031"], i.an$Df[i.an$ID == "11CEG.0031"], i.an$Df.Pooled[i.an$ID == "11CEG.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "11CEG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "11CEG.0031"], i.an$Df[i.an$ID == "11CEG.0031"], i.an$Df.Pooled[i.an$ID == "11CEG.0031"], lower.tail = T)

      # 	If varTxCC <= 0 and varTxS(IC) > 0, then F = (MSTxIC / MSpooled), where MSpooled = (SSTxCC + SSTxS(IC)) / (dfTxCC) + dfTxS(IC)) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.10CE.0030)
      i.an$ID[i.an$Source.of.Variation == "C:E"]<-"10CE.0031"
      i.an$Df.Formula[i.an$ID == "10CE.0031"]<-"(t-1)"
      i.an$Calculated.From[i.an$ID == "10CE.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "10CE.0031"]<-i.an$Sum.Sq[i.an$ID == "10CE.0031"] / i.an$Df[i.an$ID == "10CE.0031"]
      i.an$Source.of.Variation[i.an$ID == "10CE.0031"]<-"T x I"
      i.an$Denominator[i.an$ID == "10CE.0031"]<-"TxC+TxS(I)"
      i.an$Df.Pooled[i.an$ID == "10CE.0031"]<-c.an$Df[c.an$ID == "10CE.0040"] + i.an$Df[i.an$ID == "11CEG.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "10CE.0031"]<-c.an$Sum.Sq[c.an$ID == "10CE.0040"] + i.an$Sum.Sq[i.an$ID == "11CEG.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "10CE.0031"]<- i.an$Sum.Sq.Pooled[i.an$ID == "10CE.0031"] / i.an$Df.Pooled[i.an$ID == "10CE.0031"]
      i.an$F.Value[i.an$ID == "10CE.0031"]<-i.an$Mean.Sq[i.an$ID == "10CE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "10CE.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "10CE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "10CE.0031"], i.an$Df[i.an$ID == "10CE.0031"], i.an$Df.Pooled[i.an$ID == "10CE.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "10CE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "10CE.0031"], i.an$Df[i.an$ID == "10CE.0031"], i.an$Df.Pooled[i.an$ID == "10CE.0031"], lower.tail = T)

      # 	If varTxCC <= 0 and varTxS(IC) <= 0, then F = (MSTxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.10CE.0030)
      i.an$ID[i.an$Source.of.Variation == "C:E"]<-"10CE.0032"
      i.an$Df.Formula[i.an$ID == "10CE.0032"]<-"(t-1)"
      i.an$Calculated.From[i.an$ID == "10CE.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "10CE.0032"]<-i.an$Sum.Sq[i.an$ID == "10CE.0032"] / i.an$Df[i.an$ID == "10CE.0032"]
      i.an$Source.of.Variation[i.an$ID == "10CE.0032"]<-"T x I"
      i.an$Denominator[i.an$ID == "10CE.0032"]<-"TxC+TxS(I)+Residuals"
      i.an$Df.Pooled[i.an$ID == "10CE.0032"]<-i.an$Df.Pooled[i.an$ID == "10CE.0031"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "10CE.0032"]<-i.an$Sum.Sq.Pooled[i.an$ID == "10CE.0031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "10CE.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "10CE.0032"] / i.an$Df.Pooled[i.an$ID == "10CE.0032"]
      i.an$F.Value[i.an$ID == "10CE.0032"]<-i.an$Mean.Sq[i.an$ID == "10CE.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "10CE.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "10CE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "10CE.0032"], i.an$Df[i.an$ID == "10CE.0032"], i.an$Df.Pooled[i.an$ID == "10CE.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "10CE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "10CE.0032"], i.an$Df[i.an$ID == "10CE.0032"], i.an$Df.Pooled[i.an$ID == "10CE.0032"], lower.tail = T)

      # 	If varS(CC) <= 0 and varTxS(IC) > 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SSS(CC) + SSTxS(IC)) / (dfS(CC) + dfTxS(IC)) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.05EG.0030)
      i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0031"
      i.an$Df.Formula[i.an$ID == "05EG.0031"]<-"(s-1)"
      i.an$Calculated.From[i.an$ID == "05EG.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "05EG.0031"]<-i.an$Sum.Sq[i.an$ID == "05EG.0031"] / i.an$Df[i.an$ID == "05EG.0031"]
      i.an$Source.of.Variation[i.an$ID == "05EG.0031"]<-"S(I)"
      i.an$Denominator[i.an$ID == "05EG.0031"]<-"S(C)+TxS(I)"
      i.an$Df.Pooled[i.an$ID == "05EG.0031"]<-c.an$Df[c.an$ID == "05EG.0040"] + i.an$Df[i.an$ID == "11CEG.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0031"]<-c.an$Sum.Sq[c.an$ID == "05EG.0040"] + i.an$Sum.Sq[i.an$ID == "11CEG.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0031"]<- i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0031"] / i.an$Df.Pooled[i.an$ID == "05EG.0031"]
      i.an$F.Value[i.an$ID == "05EG.0031"]<-i.an$Mean.Sq[i.an$ID == "05EG.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "05EG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0031"], i.an$Df[i.an$ID == "05EG.0031"], i.an$Df.Pooled[i.an$ID == "05EG.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "05EG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0031"], i.an$Df[i.an$ID == "05EG.0031"], i.an$Df.Pooled[i.an$ID == "05EG.0031"], lower.tail = T)

      # 	If varS(CC) <= 0 and varTxS(IC) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.05EG.0030)
      i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0032"
      i.an$Df.Formula[i.an$ID == "05EG.0032"]<-"(s-1)"
      i.an$Calculated.From[i.an$ID == "05EG.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "05EG.0032"]<-i.an$Sum.Sq[i.an$ID == "05EG.0032"] / i.an$Df[i.an$ID == "05EG.0032"]
      i.an$Source.of.Variation[i.an$ID == "05EG.0032"]<-"S(I)"
      i.an$Denominator[i.an$ID == "05EG.0032"]<-"S(C)+TxS(I)+Residuals"
      i.an$Df.Pooled[i.an$ID == "05EG.0032"]<-i.an$Df.Pooled[i.an$ID == "05EG.0031"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0032"]<-i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0032"] / i.an$Df.Pooled[i.an$ID == "05EG.0032"]
      i.an$F.Value[i.an$ID == "05EG.0032"]<-i.an$Mean.Sq[i.an$ID == "05EG.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "05EG.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0032"], i.an$Df[i.an$ID == "05EG.0032"], i.an$Df.Pooled[i.an$ID == "05EG.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "05EG.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0032"], i.an$Df[i.an$ID == "05EG.0032"], i.an$Df.Pooled[i.an$ID == "05EG.0032"], lower.tail = T)

      # 	If varCC <= 0, varS(IC) <= 0 and varTxIC > 0, then F = (MS_IC / MSpooled), where MSpooled = (SS_CC + SS_TxIC) / (dfCC + dfTxIC)
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0031"
      i.an$Df.Formula[i.an$ID == "02E.0031"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0031"]<-i.an$Sum.Sq[i.an$ID == "02E.0031"] / i.an$Df[i.an$ID == "02E.0031"]
      i.an$Source.of.Variation[i.an$ID == "02E.0031"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0031"]<-"C+TxI"
      i.an$Df.Pooled[i.an$ID == "02E.0031"]<-c.an$Df[c.an$ID == "02E.0040"] + i.an$Df[i.an$ID == "10CE.0031"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"]<-c.an$Sum.Sq[c.an$ID == "02E.0040"] + i.an$Sum.Sq[i.an$ID == "10CE.0031"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0031"]<- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"] / i.an$Df.Pooled[i.an$ID == "02E.0031"]
      i.an$F.Value[i.an$ID == "02E.0031"]<-i.an$Mean.Sq[i.an$ID == "02E.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df.Pooled[i.an$ID == "02E.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df.Pooled[i.an$ID == "02E.0031"], lower.tail = T)

      # 	If varCC <= 0, varTxIC <= 0, varS(IC) > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSS(IC)) / (dfCC + dfS(IC))
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0032"
      i.an$Df.Formula[i.an$ID == "02E.0032"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0032"]<-i.an$Sum.Sq[i.an$ID == "02E.0032"] / i.an$Df[i.an$ID == "02E.0032"]
      i.an$Source.of.Variation[i.an$ID == "02E.0032"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0032"]<-"C+S(I)"
      i.an$Df.Pooled[i.an$ID == "02E.0032"]<-c.an$Df[c.an$ID == "02E.0040"] + i.an$Df[i.an$ID == "05EG.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"]<-c.an$Sum.Sq[c.an$ID == "02E.0040"] + i.an$Sum.Sq[i.an$ID == "05EG.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"] / i.an$Df.Pooled[i.an$ID == "02E.0032"]
      i.an$F.Value[i.an$ID == "02E.0032"]<-i.an$Mean.Sq[i.an$ID == "02E.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = T)

      # 	If varCC <= 0, varTxIC <= 0, varS(IC) <= 0 and varTxS(IC) > 0, then F = (MSIC / MSpooled), Where MSpooled = (SSCC + SSTxIC + SSS(IC) + SSTxS(IC)) / (dfCC + dfTxIC + dfS(IC) + dfTxS(IC)) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0033"
      i.an$Df.Formula[i.an$ID == "02E.0033"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0033"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0033"]<-i.an$Sum.Sq[i.an$ID == "02E.0033"] / i.an$Df[i.an$ID == "02E.0033"]
      i.an$Source.of.Variation[i.an$ID == "02E.0033"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0033"]<-"C+TxI+S(I)+TxS(I)"
      i.an$Df.Pooled[i.an$ID == "02E.0033"]<-c.an$Df[c.an$ID == "02E.0040"] + i.an$Df[i.an$ID == "10CE.0031"] + i.an$Df[i.an$ID == "05EG.0030"] + i.an$Df[i.an$ID == "11CEG.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0033"]<-c.an$Sum.Sq[c.an$ID == "02E.0040"] + i.an$Sum.Sq[i.an$ID == "10CE.0031"] + i.an$Sum.Sq[i.an$ID == "05EG.0030"] + i.an$Sum.Sq[i.an$ID == "11CEG.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0033"]<- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0033"] / i.an$Df.Pooled[i.an$ID == "02E.0033"]
      i.an$F.Value[i.an$ID == "02E.0033"]<-i.an$Mean.Sq[i.an$ID == "02E.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0033"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0033"], i.an$Df[i.an$ID == "02E.0033"], i.an$Df.Pooled[i.an$ID == "02E.0033"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0033"], i.an$Df[i.an$ID == "02E.0033"], i.an$Df.Pooled[i.an$ID == "02E.0033"], lower.tail = T)

      # 	If varCC <= 0, varTxIC <= 0, varS(IC) <= 0 and varTxS(IC) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0034"
      i.an$Df.Formula[i.an$ID == "02E.0034"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0034"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0034"]<-i.an$Sum.Sq[i.an$ID == "02E.0034"] / i.an$Df[i.an$ID == "02E.0034"]
      i.an$Source.of.Variation[i.an$ID == "02E.0034"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0034"]<-"C+TxI+S(I)+TxS(I)+Residuals"
      i.an$Df.Pooled[i.an$ID == "02E.0034"]<-i.an$Df.Pooled[i.an$ID == "02E.0033"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0034"]<-i.an$Mean.Sq.Pooled[i.an$ID == "02E.0033"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0034"]<- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0034"] / i.an$Df.Pooled[i.an$ID == "02E.0034"]
      i.an$F.Value[i.an$ID == "02E.0034"]<-i.an$Mean.Sq[i.an$ID == "02E.0034"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0034"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0034"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0034"], i.an$Df[i.an$ID == "02E.0034"], i.an$Df.Pooled[i.an$ID == "02E.0034"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0034"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0034"], i.an$Df[i.an$ID == "02E.0034"], i.an$Df.Pooled[i.an$ID == "02E.0034"], lower.tail = T)
    }
  }


  ###------------------ INTERACTIONS with ONE temporal and ONE spatial factors

  if(aav.model == "CE") {
    # 	F = (MSTxCC) / MSRes)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "C:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "C:E"]<-"10CE.0040"
    c.an$Df.Formula[c.an$ID == "10CE.0040"]<-"(t-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "10CE.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "10CE.0040"] <- c.an$Sum.Sq[c.an$ID == "10CE.0040"] / c.an$Df[c.an$ID == "10CE.0040"]
    c.an$Source.of.Variation[c.an$ID == "10CE.0040"]<-"T x C"
    c.an$Denominator[c.an$ID == "10CE.0040"]<-"Residuals"
    c.an$F.Value[c.an$ID == "10CE.0040"]<-c.an$Mean.Sq[c.an$ID == "10CE.0040"] / f.an$Mean.Sq[f.an$ID == "12R.0000"]
    c.an$P.Value.Upper.Tail[c.an$ID == "10CE.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "10CE.0040"], c.an$Df[c.an$ID == "10CE.0040"], f.an$Df[f.an$ID == "12R.0000"], lower.tail = F)

    # 	F = (MSTxIC) / MSTxCC), or
    dii.10CE.0030<-cbind(Source.of.Variation = "C:E",
                         sym.F[sym.F$Source.of.Variation == "C:E",2:3] - sym.C[sym.C$Source.of.Variation == "C:E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.10CE.0030)
    i.an$ID[i.an$Source.of.Variation == "C:E"]<-"10CE.0030"
    i.an$Df.Formula[i.an$ID == "10CE.0030"]<-"(t-1)"
    i.an$Calculated.From[i.an$ID == "10CE.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "10CE.0030"]<-i.an$Sum.Sq[i.an$ID == "10CE.0030"] / i.an$Df[i.an$ID == "10CE.0030"]
    i.an$Source.of.Variation[i.an$ID == "10CE.0030"]<-"T x I"
    i.an$Denominator[i.an$ID == "10CE.0030"]<-"TxC"
    i.an$F.Value[i.an$ID == "10CE.0030"]<-i.an$Mean.Sq[i.an$ID == "10CE.0030"] / c.an$Mean.Sq[c.an$ID == "10CE.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "10CE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "10CE.0030"], i.an$Df[i.an$ID == "10CE.0030"], c.an$Df[c.an$ID == "10CE.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "10CE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "10CE.0030"], i.an$Df[i.an$ID == "10CE.0030"], c.an$Df[c.an$ID == "10CE.0040"], lower.tail = T)

    # 	F = (MSCC) / MSTxCC)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0040"
    c.an$Df.Formula[c.an$ID == "02E.0040"]<-"(c-1)"
    c.an$Calculated.From[c.an$ID == "02E.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "02E.0040"] <- c.an$Sum.Sq[c.an$ID == "02E.0040"] / c.an$Df[c.an$ID == "02E.0040"]
    c.an$Source.of.Variation[c.an$ID == "02E.0040"]<-"Among Controls = C"
    c.an$Denominator[c.an$ID == "02E.0040"]<-"TxC"
    c.an$F.Value[c.an$ID == "02E.0040"]<-c.an$Mean.Sq[c.an$ID == "02E.0040"] / c.an$Mean.Sq[c.an$ID == "10CE.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "02E.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "02E.0040"], c.an$Df[c.an$ID == "02E.0040"], c.an$Df[c.an$ID == "10CE.0040"], lower.tail = F)

    # 	F = (MSIC / MSCC), or
    dii.02E.0030<-cbind(Source.of.Variation = "E",
                        sym.F[sym.F$Source.of.Variation == "E",2:3] - sym.C[sym.C$Source.of.Variation == "E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0030"
    i.an$Df.Formula[i.an$ID == "02E.0030"]<-"1"
    i.an$Calculated.From[i.an$ID == "02E.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "02E.0030"]<-i.an$Sum.Sq[i.an$ID == "02E.0030"] / i.an$Df[i.an$ID == "02E.0030"]
    i.an$Source.of.Variation[i.an$ID == "02E.0030"]<-"Impact vs Controls = I"
    i.an$Denominator[i.an$ID == "02E.0030"]<-"C"
    i.an$F.Value[i.an$ID == "02E.0030"]<-i.an$Mean.Sq[i.an$ID == "02E.0030"] / c.an$Mean.Sq[c.an$ID == "02E.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], c.an$Df[c.an$ID == "02E.0040"], lower.tail = F)
    i.an$P.Value.Lower.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], c.an$Df[c.an$ID == "02E.0040"], lower.tail = T)

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	If varTxCC <= 0, then F = (MSTxIC / MSpooled), where MSpooled = (SSTxCC + SSRes) / (dfTxCC) + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.10CE.0030)
      i.an$ID[i.an$Source.of.Variation == "C:E"]<-"10CE.0031"
      i.an$Df.Formula[i.an$ID == "10CE.0031"]<-"(t-1)"
      i.an$Calculated.From[i.an$ID == "10CE.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "10CE.0031"]<-i.an$Sum.Sq[i.an$ID == "10CE.0031"] / i.an$Df[i.an$ID == "10CE.0031"]
      i.an$Source.of.Variation[i.an$ID == "10CE.0031"]<-"T x I"
      i.an$Denominator[i.an$ID == "10CE.0031"]<-"TxC+Residuals"
      i.an$Df.Pooled[i.an$ID == "10CE.0031"]<- c.an$Df[c.an$ID == "10CE.0040"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "10CE.0031"]<- c.an$Sum.Sq[c.an$ID == "10CE.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "10CE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "10CE.0031"] / i.an$Df.Pooled[i.an$ID == "10CE.0031"]
      i.an$F.Value[i.an$ID == "10CE.0031"]<-i.an$Mean.Sq[i.an$ID == "10CE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "10CE.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "10CE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "10CE.0031"], i.an$Df[i.an$ID == "10CE.0031"], i.an$Df.Pooled[i.an$ID == "10CE.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "10CE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "10CE.0031"], i.an$Df[i.an$ID == "10CE.0031"], i.an$Df.Pooled[i.an$ID == "10CE.0031"], lower.tail = T)

      # 	If varCC <= 0 and varTxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSTxIC) / (dfCC + dfTxIC) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0031"
      i.an$Df.Formula[i.an$ID == "02E.0031"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0031"]<-i.an$Sum.Sq[i.an$ID == "02E.0031"] / i.an$Df[i.an$ID == "02E.0031"]
      i.an$Source.of.Variation[i.an$ID == "02E.0031"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0031"]<-"C+TxI"
      i.an$Df.Pooled[i.an$ID == "02E.0031"]<- c.an$Df[c.an$ID == "02E.0040"] + i.an$Df[i.an$ID == "10CE.0030"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"]<- c.an$Sum.Sq[c.an$ID == "02E.0040"] + i.an$Sum.Sq[i.an$ID == "10CE.0030"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"] / i.an$Df.Pooled[i.an$ID == "02E.0031"]
      i.an$F.Value[i.an$ID == "02E.0031"]<-i.an$Mean.Sq[i.an$ID == "02E.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df.Pooled[i.an$ID == "02E.0031"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df.Pooled[i.an$ID == "02E.0031"], lower.tail = T)

      # 	If varCC <= 0 and varTxIC <= 0, then F = (MSIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0032"
      i.an$Df.Formula[i.an$ID == "02E.0032"]<-"1"
      i.an$Calculated.From[i.an$ID == "02E.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0032"]<-i.an$Df.Pooled[i.an$ID == "02E.0031"] / i.an$Df[i.an$ID == "02E.0032"]
      i.an$Source.of.Variation[i.an$ID == "02E.0032"]<-"Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0032"]<-"C+TxI+Residuals"
      i.an$Df.Pooled[i.an$ID == "02E.0032"]<- i.an$Df.Pooled[i.an$ID == "02E.0031"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"] / i.an$Df.Pooled[i.an$ID == "02E.0032"]
      i.an$F.Value[i.an$ID == "02E.0032"]<-i.an$Mean.Sq[i.an$ID == "02E.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = F)
      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = T)
    }
  }


  ###------------------ INTERACTIONS with THREE spatial factors

  if(aav.model == "DEG") {
    # 	F = (MS_S(CC(P)) / MSRes)
    c.an<-dplyr::bind_rows(c.an,sym.NC[sym.NC$Source.of.Variation == "D:E:G",1:3])
    c.an$ID[c.an$Source.of.Variation == "D:E:G"]<-"08DEGn.0040"
    c.an$Df.Formula[c.an$ID == "08DEGn.0040"]<-"pc(s-1)"
    c.an$Calculated.From[c.an$ID == "08DEGn.0040"]<-"sym.NC"
    c.an$Mean.Sq[c.an$ID == "08DEGn.0040"] <- c.an$Sum.Sq[c.an$ID == "08DEGn.0040"] / c.an$Df[c.an$ID == "08DEGn.0040"]
    c.an$Source.of.Variation[c.an$ID == "08DEGn.0040"]<-"Control sites(C(P)) = S(C(P))"
    c.an$Denominator[c.an$ID == "08DEGn.0040"]<-"Residuals"
    c.an$F.Value[c.an$ID == "08DEGn.0040"]<-c.an$Mean.Sq[c.an$ID == "08DEGn.0040"] / f.an$Mean.Sq[f.an$ID == "12R.0000"]
    c.an$P.Value.Upper.Tail[c.an$ID == "08DEGn.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "08DEGn.0040"], c.an$Df[c.an$ID == "08DEGn.0040"], f.an$Df[f.an$ID == "12R.0000"], lower.tail = F)

    # 	F = (MS_CC(P) / MSS(CC(P))), or
    c.an<-dplyr::bind_rows(c.an,sym.NC[sym.NC$Source.of.Variation == "D:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "D:E"]<-"07DEn.0040"
    c.an$Df.Formula[c.an$ID == "07DEn.0040"]<-"p(c-1)"
    c.an$Calculated.From[c.an$ID == "07DEn.0040"]<-"sym.NC"
    c.an$Mean.Sq[c.an$ID == "07DEn.0040"] <- c.an$Sum.Sq[c.an$ID == "07DEn.0040"] / c.an$Df[c.an$ID == "07DEn.0040"]
    c.an$Source.of.Variation[c.an$ID == "07DEn.0040"]<-"Control location(P) = C(P)"
    c.an$Denominator[c.an$ID == "07DEn.0040"]<-"S(C(P))"
    c.an$F.Value[c.an$ID == "07DEn.0040"]<-c.an$Mean.Sq[c.an$ID == "07DEn.0040"] / c.an$Mean.Sq[c.an$ID == "08DEGn.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "07DEn.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "07DEn.0040"], c.an$Df[c.an$ID == "07DEn.0040"], c.an$Df[c.an$ID == "08DEGn.0040"], lower.tail = F)

    # 	F = (MS_PxS(IC) / MSS(CC(P))), or
    dii.06DEG.0030<-cbind(Source.of.Variation = "D:E:G",
                          sym.F[sym.F$Source.of.Variation == "D:E:G",2:3] - sym.C[sym.C$Source.of.Variation == "D:E:G",2:3])
    i.an<-dplyr::bind_rows(i.an, dii.06DEG.0030)
    i.an$ID[i.an$Source.of.Variation == "D:E:G"]<-"06DEG.0030"
    i.an$Df.Formula[i.an$ID == "06DEG.0030"]<-"(p-1)(s-1)"
    i.an$Calculated.From[i.an$ID == "06DEG.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "06DEG.0030"] <- i.an$Sum.Sq[i.an$ID == "06DEG.0030"] / i.an$Df[i.an$ID == "06DEG.0030"]
    i.an$Source.of.Variation[i.an$ID == "06DEG.0030"]<-"P x S(I)"
    i.an$Denominator[i.an$ID == "06DEG.0030"]<-"S(C(P))"
    i.an$F.Value[i.an$ID == "06DEG.0030"]<-i.an$Mean.Sq[i.an$ID == "06DEG.0030"] / c.an$Mean.Sq[c.an$ID == "08DEGn.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "06DEG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "06DEG.0030"], i.an$Df[i.an$ID == "06DEG.0030"], c.an$Df[c.an$ID == "08DEGn.0040"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "06DEG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "06DEG.0030"], i.an$Df[i.an$ID == "06DEG.0030"], c.an$Df[c.an$ID == "08DEGn.0040"], lower.tail = T)

    # 	F = (MS_S(IC) / MSPxS(IC)), or
    dii.05EG.0030<-cbind(Source.of.Variation = "E:G",
                         sym.F[sym.F$Source.of.Variation == "E:G",2:3] - sym.C[sym.C$Source.of.Variation == "E:G",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.05EG.0030)
    i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0030"
    i.an$Df.Formula[i.an$ID == "05EG.0030"]<-"(s-1)"
    i.an$Calculated.From[i.an$ID == "05EG.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "05EG.0030"]<-i.an$Sum.Sq[i.an$ID == "05EG.0030"] / i.an$Df[i.an$ID == "05EG.0030"]
    i.an$Source.of.Variation[i.an$ID == "05EG.0030"]<-"Sites(I) = S(I)"
    i.an$Denominator[i.an$ID == "05EG.0030"]<-"PxS(I)"
    i.an$F.Value[i.an$ID == "05EG.0030"]<-i.an$Mean.Sq[i.an$ID == "05EG.0030"] / i.an$Mean.Sq[i.an$ID == "06DEG.0030"]
    i.an$P.Value.Upper.Tail[i.an$ID == "05EG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0030"], i.an$Df[i.an$ID == "05EG.0030"], i.an$Df[i.an$ID == "06DEG.0030"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "05EG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0030"], i.an$Df[i.an$ID == "05EG.0030"], i.an$Df[i.an$ID == "06DEG.0030"], lower.tail = T)

    # 	If varPxS(IC) <= 0 and varCC(P) > 0, then F = (MS_PxIC / MS_CC(P))
    dii.03DE.0030<-cbind(Source.of.Variation = "D:E",
                         sym.F[sym.F$Source.of.Variation == "D:E",2:3] - sym.C[sym.C$Source.of.Variation == "D:E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.03DE.0030)
    i.an$ID[i.an$Source.of.Variation == "D:E"]<-"03DE.0030"
    i.an$Df.Formula[i.an$ID == "03DE.0030"]<-"(p-1)(t-1)"
    i.an$Calculated.From[i.an$ID == "03DE.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "03DE.0030"]<-i.an$Sum.Sq[i.an$ID == "03DE.0030"] / i.an$Df[i.an$ID == "03DE.0030"]
    i.an$Source.of.Variation[i.an$ID == "03DE.0030"]<-"P x T = P x I"
    i.an$Denominator[i.an$ID == "03DE.0030"]<-"C(P)"
    i.an$F.Value[i.an$ID == "03DE.0030"]<-i.an$Mean.Sq[i.an$ID == "03DE.0030"] / c.an$Mean.Sq[c.an$ID == "07DEn.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "03DE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0030"], i.an$Df[i.an$ID == "03DE.0030"], c.an$Df[c.an$ID == "07DEn.0040"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "03DE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0030"], i.an$Df[i.an$ID == "03DE.0030"], c.an$Df[c.an$ID == "07DEn.0040"], lower.tail = T)

    # 	If varCC(P) <= 0 and varPxS(IC) > 0, then F = (MS_PxIC / MSPxS(IC))
    i.an<-dplyr::bind_rows(i.an,dii.03DE.0030)
    i.an$ID[i.an$Source.of.Variation == "D:E"]<-"03DE.0031"
    i.an$Df.Formula[i.an$ID == "03DE.0031"]<-"(p-1)(t-1)"
    i.an$Calculated.From[i.an$ID == "03DE.0031"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "03DE.0031"]<-i.an$Sum.Sq[i.an$ID == "03DE.0031"] / i.an$Df[i.an$ID == "03DE.0031"]
    i.an$Source.of.Variation[i.an$ID == "03DE.0031"]<-"P x T = P x I"
    i.an$Denominator[i.an$ID == "03DE.0031"]<-"PxS(I)"
    i.an$F.Value[i.an$ID == "03DE.0031"]<-i.an$Mean.Sq[i.an$ID == "03DE.0031"] / i.an$Mean.Sq[i.an$ID == "06DEG.0030"]
    i.an$P.Value.Upper.Tail[i.an$ID == "03DE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0031"], i.an$Df[i.an$ID == "03DE.0031"], i.an$Df[i.an$ID == "06DEG.0030"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "03DE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0031"], i.an$Df[i.an$ID == "03DE.0031"], i.an$Df[i.an$ID == "06DEG.0030"], lower.tail = T)

    # 	If varCC(P) <= 0, varPxS(IC) <= 0 and varS(CC(P)) > 0, then F = (MS_PxIC / MSpooled), where MSpooled = (SS_CC(P) + SS_PxS(IC) + SS_S(CC(P))) / (dfCC(P) + dfPxS(IC) + dfS(CC(P))) = SS1/df1
    i.an<-dplyr::bind_rows(i.an,dii.03DE.0030)
    i.an$ID[i.an$Source.of.Variation == "D:E"]<-"03DE.0032"
    i.an$Df.Formula[i.an$ID == "03DE.0032"]<-"(p-1)(t-1)"
    i.an$Calculated.From[i.an$ID == "03DE.0032"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "03DE.0032"]<-i.an$Sum.Sq[i.an$ID == "03DE.0032"] / i.an$Df[i.an$ID == "03DE.0032"]
    i.an$Source.of.Variation[i.an$ID == "03DE.0032"]<-"P x T = P x I"
    i.an$Denominator[i.an$ID == "03DE.0032"]<-"C(P)+PxS(I)+S(C(P))"
    i.an$Df.Pooled[i.an$ID == "03DE.0032"]<-c.an$Df[c.an$ID == "07DEn.0040"] + i.an$Df[i.an$ID == "06DEG.0030"] + c.an$Df[c.an$ID == "08DEGn.0040"]
    i.an$Sum.Sq.Pooled[i.an$ID == "03DE.0032"]<-c.an$Sum.Sq[c.an$ID == "07DEn.0040"] + i.an$Sum.Sq[i.an$ID == "06DEG.0030"] + c.an$Sum.Sq[c.an$ID == "08DEGn.0040"]
    i.an$Mean.Sq.Pooled[i.an$ID == "03DE.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03DE.0032"] / i.an$Df.Pooled[i.an$ID == "03DE.0032"]
    i.an$F.Value[i.an$ID == "03DE.0032"]<-i.an$Mean.Sq[i.an$ID == "03DE.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "03DE.0032"]
    i.an$P.Value.Upper.Tail[i.an$ID == "03DE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0032"], i.an$Df[i.an$ID == "03DE.0032"], i.an$Df.Pooled[i.an$ID == "03DE.0032"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "03DE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0032"], i.an$Df[i.an$ID == "03DE.0032"], i.an$Df.Pooled[i.an$ID == "03DE.0032"], lower.tail = T)

    # 	If varS(IC) <= 0 and varPxIC > 0, then F = (MS_IC / MS_PxIC)
    dii.02E.0030<-cbind(Source.of.Variation = "E",
                        sym.F[sym.F$Source.of.Variation == "E",2:3] - sym.C[sym.C$Source.of.Variation == "E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0030"
    i.an$Df.Formula[i.an$ID == "02E.0030"]<-"(t-1)"
    i.an$Calculated.From[i.an$ID == "02E.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "02E.0030"]<-i.an$Sum.Sq[i.an$ID == "02E.0030"] / i.an$Df[i.an$ID == "02E.0030"]
    i.an$Source.of.Variation[i.an$ID == "02E.0030"]<-"Treatment = T = Impact vs Controls = I"
    i.an$Denominator[i.an$ID == "02E.0030"]<-"PxI"
    i.an$F.Value[i.an$ID == "02E.0030"]<-i.an$Mean.Sq[i.an$ID == "02E.0030"] / i.an$Mean.Sq[i.an$ID == "03DE.0030"]
    i.an$P.Value.Upper.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "03DE.0030"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "03DE.0030"], lower.tail = T)

    # 	If varPxIC <= 0 and varS(IC) > 0, then F = (MS_IC / MS_S(IC))
    i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0031"
    i.an$Df.Formula[i.an$ID == "02E.0031"]<-"(t-1)"
    i.an$Calculated.From[i.an$ID == "02E.0031"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "02E.0031"]<-i.an$Sum.Sq[i.an$ID == "02E.0031"] / i.an$Df[i.an$ID == "02E.0031"]
    i.an$Source.of.Variation[i.an$ID == "02E.0031"]<-"Treatment = T = Impact vs Controls = I"
    i.an$Denominator[i.an$ID == "02E.0031"]<-"S(I)"
    i.an$F.Value[i.an$ID == "02E.0031"]<-i.an$Mean.Sq[i.an$ID == "02E.0031"] / i.an$Mean.Sq[i.an$ID == "05EG.0030"]
    i.an$P.Value.Upper.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "05EG.0030"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "05EG.0030"], lower.tail = T)

    # 	If varPxIC <= 0, varS(IC) <= 0, varCC(P) <= 0, varPxS(IC) > 0, then F = (MSIC / MSpooled), where MSpooled = (SS_PxIC + SS_S(IC) + SS_PxS(IC)) / (dfPxIC + dfS(IC) + dfPxS(IC))
    i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0032"
    i.an$Df.Formula[i.an$ID == "02E.0032"]<-"(t-1)"
    i.an$Calculated.From[i.an$ID == "02E.0032"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "02E.0032"]<-i.an$Sum.Sq[i.an$ID == "02E.0032"] / i.an$Df[i.an$ID == "02E.0032"]
    i.an$Source.of.Variation[i.an$ID == "02E.0032"]<-"Treatment = T = Impact vs Controls = I"
    i.an$Denominator[i.an$ID == "02E.0032"]<-"PxI+S(I)+PxS(I)"
    i.an$Df.Pooled[i.an$ID == "02E.0032"]<-i.an$Df[i.an$ID == "03DE.0030"] + i.an$Df[i.an$ID == "05EG.0030"] + i.an$Df[i.an$ID == "06DEG.0030"]
    i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"]<-i.an$Sum.Sq[i.an$ID == "03DE.0030"] + i.an$Sum.Sq[i.an$ID == "05EG.0030"] +i.an$Sum.Sq[i.an$ID == "06DEG.0030"]
    i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"] / i.an$Df.Pooled[i.an$ID == "02E.0032"]
    i.an$F.Value[i.an$ID == "02E.0032"]<-i.an$Mean.Sq[i.an$ID == "02E.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"]
    i.an$P.Value.Upper.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = T)

    # 	If varPxIC <= 0, varS(IC) <= 0, varPxS(IC) <= 0 or varCC(P) > 0, then F = (MSIC / MSpooled), where MSpooled = (SSPxIC + SSS(IC) + SS_CC(P)) / (dfPxIC + dfS(IC) + dfCC(P))
    i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0033"
    i.an$Df.Formula[i.an$ID == "02E.0033"]<-"(t-1)"
    i.an$Calculated.From[i.an$ID == "02E.0033"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "02E.0033"]<-i.an$Sum.Sq[i.an$ID == "02E.0033"] / i.an$Df[i.an$ID == "02E.0033"]
    i.an$Source.of.Variation[i.an$ID == "02E.0033"]<-"Treatment = T = Impact vs Controls = I"
    i.an$Denominator[i.an$ID == "02E.0033"]<-"PxI+S(I)+C(P)"
    i.an$Df.Pooled[i.an$ID == "02E.0033"]<-i.an$Df[i.an$ID == "03DE.0030"]  + i.an$Df[i.an$ID == "05EG.0030"] + c.an$Df[c.an$ID == "07DEn.0040"]
    i.an$Sum.Sq.Pooled[i.an$ID == "02E.0033"]<-i.an$Sum.Sq[i.an$ID == "03DE.0030"] + i.an$Sum.Sq[i.an$ID == "05EG.0030"] + c.an$Sum.Sq[c.an$ID == "07DEn.0040"]
    i.an$Mean.Sq.Pooled[i.an$ID == "02E.0033"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0033"] / i.an$Df.Pooled[i.an$ID == "02E.0033"]
    i.an$F.Value[i.an$ID == "02E.0033"]<-i.an$Mean.Sq[i.an$ID == "02E.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0033"]
    i.an$P.Value.Upper.Tail[i.an$ID == "02E.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0033"], i.an$Df[i.an$ID == "02E.0033"], i.an$Df.Pooled[i.an$ID == "02E.0033"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0033"], i.an$Df[i.an$ID == "02E.0033"], i.an$Df.Pooled[i.an$ID == "02E.0033"], lower.tail = T)

    # 	If varPxIC <= 0, varS(IC) <= 0, varPxS(IC) <= 0, varCC(P) <= 0 and varS(CC(P)) > 0, then F = (MSIC / MSpooled), where MSpooled = (SS_PxIC + SS_S(IC) + SS_PxS(IC) + SS_CC(P) + SS_S(CC(P))) / (dfPxIC + dfS(IC) + dfPxS(IC) + dfCC(P) + dfS(CC(P))) = SS1/df1
    i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0034"
    i.an$Df.Formula[i.an$ID == "02E.0034"]<-"(t-1)"
    i.an$Calculated.From[i.an$ID == "02E.0034"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "02E.0034"]<-i.an$Sum.Sq[i.an$ID == "02E.0034"] / i.an$Df[i.an$ID == "02E.0034"]
    i.an$Source.of.Variation[i.an$ID == "02E.0034"]<-"Treatment = T = Impact vs Controls = I"
    i.an$Denominator[i.an$ID == "02E.0034"]<-"PxI+S(I)+PxS(I)+C(P)+S(C(P))"
    i.an$Df.Pooled[i.an$ID == "02E.0034"]<-i.an$Df[i.an$ID == "03DE.0030"]  + i.an$Df[i.an$ID == "05EG.0030"] + i.an$Df[i.an$ID == "06DEG.0030"] + c.an$Df[c.an$ID == "07DEn.0040"] + c.an$Df[c.an$ID == "08DEGn.0040"]
    i.an$Sum.Sq.Pooled[i.an$ID == "02E.0034"]<-i.an$Sum.Sq[i.an$ID == "03DE.0030"] + i.an$Sum.Sq[i.an$ID == "05EG.0030"] + c.an$Sum.Sq[c.an$ID == "07DEn.0040"] + c.an$Sum.Sq[c.an$ID == "08DEGn.0040"]
    i.an$Mean.Sq.Pooled[i.an$ID == "02E.0034"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0034"] / i.an$Df.Pooled[i.an$ID == "02E.0034"]
    i.an$F.Value[i.an$ID == "02E.0034"]<-i.an$Mean.Sq[i.an$ID == "02E.0034"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0034"]
    i.an$P.Value.Upper.Tail[i.an$ID == "02E.0034"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0034"], i.an$Df[i.an$ID == "02E.0034"], i.an$Df.Pooled[i.an$ID == "02E.0034"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0034"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0034"], i.an$Df[i.an$ID == "02E.0034"], i.an$Df.Pooled[i.an$ID == "02E.0034"], lower.tail = T)

    # 	F = (MSP / MS_CC(P))
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "D",1:3])
    f.an$ID[f.an$Source.of.Variation == "D"]<-"01D.0000"
    f.an$Df.Formula[f.an$ID == "01D.0000"]<-"(p-1)"
    f.an$Calculated.From[f.an$ID == "01D.0000"]<-"sym.F"
    f.an$Mean.Sq[f.an$ID == "01D.0000"]<-f.an$Sum.Sq[f.an$ID == "01D.0000"] / f.an$Df[f.an$ID == "01D.0000"]
    f.an$Source.of.Variation[f.an$ID == "01D.0000"]<-"Place = P"
    f.an$Denominator[f.an$ID == "01D.0000"]<-"C(P)"
    f.an$F.Value[f.an$ID == "01D.0000"]<-f.an$Mean.Sq[f.an$ID == "01D.0000"] / c.an$Mean.Sq[c.an$ID == "07DEn.0040"]
    f.an$P.Value.Upper.Tail[f.an$ID == "01D.0000"]<- stats::pf(f.an$F.Value[f.an$ID == "01D.0000"], f.an$Df[f.an$ID == "01D.0000"], c.an$Df[c.an$ID == "07DEn.0040"], lower.tail = F)
    #      f.an$P.Value.Lower.Tail[f.an$ID == "01D.0000"]<- stats::pf(f.an$F.Value[f.an$ID == "01D.0000"], f.an$Df[f.an$ID == "01D.0000"], c.an$Df[c.an$ID == "07DEn.0040"], lower.tail = T)

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	If varS(CC(P)) <= 0, then F = (MSCC(P) / MSpooled), where MSpooled = (SS_S(CC(P)) + SSRes) / (dfS(CC(P)) + dfRes)
      c.an<-dplyr::bind_rows(c.an,sym.NC[sym.NC$Source.of.Variation == "D:E",1:3])
      c.an$ID[c.an$Source.of.Variation == "D:E"]<-"07DEn.0041"
      c.an$Df.Formula[c.an$ID == "07DEn.0041"]<-"p(c-1)"
      c.an$Calculated.From[c.an$ID == "07DEn.0041"]<-"sym.NC"
      c.an$Mean.Sq[c.an$ID == "07DEn.0041"] <- c.an$Sum.Sq[c.an$ID == "07DEn.0041"] / c.an$Df[c.an$ID == "07DEn.0041"]
      c.an$Source.of.Variation[c.an$ID == "07DEn.0041"]<-"Control location in Place = C(P)"
      c.an$Denominator[c.an$ID == "07DEn.0041"]<-"S(C(P))+Residuals"
      c.an$Df.Pooled[c.an$ID == "07DEn.0041"]<- c.an$Df[c.an$ID == "08DEGn.0040"] + f.an$Df[f.an$ID == "12R.0000"]
      c.an$Sum.Sq.Pooled[c.an$ID == "07DEn.0041"]<- c.an$Sum.Sq[c.an$ID == "08DEGn.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      c.an$Mean.Sq.Pooled[c.an$ID == "07DEn.0041"]<- c.an$Sum.Sq.Pooled[c.an$ID == "07DEn.0041"] / c.an$Df.Pooled[c.an$ID == "07DEn.0041"]
      c.an$F.Value[c.an$ID == "07DEn.0041"]<- c.an$Mean.Sq[c.an$ID == "07DEn.0041"] / c.an$Mean.Sq.Pooled[c.an$ID == "07DEn.0041"]
      c.an$P.Value.Upper.Tail[c.an$ID == "07DEn.0041"]<- stats::pf(c.an$F.Value[c.an$ID == "07DEn.0041"], c.an$Df[c.an$ID == "07DEn.0041"], c.an$Df.Pooled[c.an$ID == "07DEn.0041"], lower.tail = F)

      # 	If varS(CC(P)) <= 0, then F = (MSPxS(IC) / MSpooled), where MSpooled = (SS_S(CC(P)) + SSRes) / (dfS(CC(P)) + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.06DEG.0030)
      i.an$ID[i.an$Source.of.Variation == "D:E:G"]<-"06DEG.0031"
      i.an$Df.Formula[i.an$ID == "06DEG.0031"]<-"(p-1)(s-1)"
      i.an$Calculated.From[i.an$ID == "06DEG.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "06DEG.0031"]<-i.an$Sum.Sq[i.an$ID == "06DEG.0031"] / i.an$Df[i.an$ID == "06DEG.0031"]
      i.an$Source.of.Variation[i.an$ID == "06DEG.0031"]<-"P x S(I)"
      i.an$Denominator[i.an$ID == "06DEG.0031"]<-"S(C(P))+Residuals"
      i.an$Df.Pooled[i.an$ID == "06DEG.0031"]<- c.an$Df[c.an$ID == "08DEGn.0040"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "06DEG.0031"]<- c.an$Sum.Sq[c.an$ID == "08DEGn.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "06DEG.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "06DEG.0031"] / i.an$Df.Pooled[i.an$ID == "06DEG.0031"]
      i.an$F.Value[i.an$ID == "06DEG.0031"]<-i.an$Mean.Sq[i.an$ID == "06DEG.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "06DEG.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "06DEG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "06DEG.0031"], i.an$Df[i.an$ID == "06DEG.0031"], i.an$Df.Pooled[i.an$ID == "06DEG.0031"], lower.tail = F)
      #        i.an$P.Value.Lower.Tail[i.an$ID == "06DEG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "06DEG.0031"], i.an$Df[i.an$ID == "06DEG.0031"], i.an$Df.Pooled[i.an$ID == "06DEG.0031"], lower.tail = T)

      # 	If varPxS(IC) <= 0 and varS(CC)P)) > 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SSPxS(IC) + SSS(CC(P))) / (dfPxS(IC) + dfS(CC(P))) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.05EG.0030)
      i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0031"
      i.an$Df.Formula[i.an$ID == "05EG.0031"]<-"(s-1)"
      i.an$Calculated.From[i.an$ID == "05EG.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "05EG.0031"]<-i.an$Sum.Sq[i.an$ID == "05EG.0031"] / i.an$Df[i.an$ID == "05EG.0031"]
      i.an$Source.of.Variation[i.an$ID == "05EG.0031"]<-"Sites(I) = S(I)"
      i.an$Denominator[i.an$ID == "05EG.0031"]<-"PxS(I)+S(C(P))"
      i.an$Df.Pooled[i.an$ID == "05EG.0031"]<- i.an$Df[i.an$ID == "06DEG.0030"] + c.an$Df[c.an$ID == "08DEGn.0040"]
      i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0031"]<- i.an$Sum.Sq[i.an$ID == "06DEG.0030"] + c.an$Sum.Sq[c.an$ID == "08DEGn.0040"]
      i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0031"] / i.an$Df.Pooled[i.an$ID == "05EG.0031"]
      i.an$F.Value[i.an$ID == "05EG.0031"]<-i.an$Mean.Sq[i.an$ID == "05EG.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "05EG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0031"], i.an$Df[i.an$ID == "05EG.0031"], i.an$Df.Pooled[i.an$ID == "05EG.0031"], lower.tail = F)
      #        i.an$P.Value.Lower.Tail[i.an$ID == "05EG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0031"], i.an$Df[i.an$ID == "05EG.0031"], i.an$Df.Pooled[i.an$ID == "05EG.0031"], lower.tail = T)

      # 	If varPxS(IC) <= 0 and varS(CC(P)) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.05EG.0030)
      i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0032"
      i.an$Df.Formula[i.an$ID == "05EG.0032"]<-"(s-1)"
      i.an$Calculated.From[i.an$ID == "05EG.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "05EG.0032"]<-i.an$Sum.Sq[i.an$ID == "05EG.0032"] / i.an$Df[i.an$ID == "05EG.0032"]
      i.an$Source.of.Variation[i.an$ID == "05EG.0032"]<-"Sites(I) = S(I)"
      i.an$Denominator[i.an$ID == "05EG.0032"]<-"PxS(I)+S(C(P))+Residuals"
      i.an$Df.Pooled[i.an$ID == "05EG.0032"]<- i.an$Df.Pooled[i.an$ID == "05EG.0031"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "05EG.0032"] / i.an$Df.Pooled[i.an$ID == "05EG.0032"]
      i.an$F.Value[i.an$ID == "05EG.0032"]<-i.an$Mean.Sq[i.an$ID == "05EG.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "05EG.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "05EG.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0032"], i.an$Df[i.an$ID == "05EG.0032"], i.an$Df.Pooled[i.an$ID == "05EG.0032"], lower.tail = F)
      #        i.an$P.Value.Lower.Tail[i.an$ID == "05EG.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "05EG.0032"], i.an$Df[i.an$ID == "05EG.0032"], i.an$Df.Pooled[i.an$ID == "05EG.0032"], lower.tail = T)

      # 	If varCC(P) <= 0, varPxS(IC) <= 0 and varS(CC(P)) <= 0, then F = (MS_PxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.03DE.0030)
      i.an$ID[i.an$Source.of.Variation == "D:E"]<-"03DE.0033"
      i.an$Df.Formula[i.an$ID == "03DE.0033"]<-"(p-1)(t-1)"
      i.an$Calculated.From[i.an$ID == "03DE.0033"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "03DE.0033"]<-i.an$Sum.Sq[i.an$ID == "03DE.0033"] / i.an$Df[i.an$ID == "03DE.0033"]
      i.an$Source.of.Variation[i.an$ID == "03DE.0033"]<-"P x T = P x I"
      i.an$Denominator[i.an$ID == "03DE.0033"]<-"C(P)+PxS(I)+S(C(P))+Residuals"
      i.an$Df.Pooled[i.an$ID == "03DE.0033"]<-i.an$Df.Pooled[i.an$ID == "03DE.0032"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "03DE.0033"]<-i.an$Sum.Sq.Pooled[i.an$ID == "03DE.0032"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "03DE.0033"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03DE.0033"] / i.an$Df.Pooled[i.an$ID == "02E.0032"]
      i.an$F.Value[i.an$ID == "03DE.0033"]<-i.an$Mean.Sq[i.an$ID == "03DE.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "03DE.0033"]
      i.an$P.Value.Upper.Tail[i.an$ID == "03DE.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0033"], i.an$Df[i.an$ID == "03DE.0033"], i.an$Df.Pooled[i.an$ID == "03DE.0033"], lower.tail = F)
      #        i.an$P.Value.Lower.Tail[i.an$ID == "03DE.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0033"], i.an$Df[i.an$ID == "03DE.0033"], i.an$Df.Pooled[i.an$ID == "03DE.0033"], lower.tail = T)

      # 	If varPxIC <= 0, varS(IC) <= 0, varPxS(IC) <= 0, varCC(P) <= 0 and varS(CC(P)) <= 0, then F = (MSIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0035"
      i.an$Df.Formula[i.an$ID == "02E.0035"]<-"(t-1)"
      i.an$Calculated.From[i.an$ID == "02E.0035"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0035"]<-i.an$Sum.Sq[i.an$ID == "02E.0035"] / i.an$Df[i.an$ID == "02E.0035"]
      i.an$Source.of.Variation[i.an$ID == "02E.0035"]<-"Treatment = T = Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0035"]<-"PxI+S(I)+PxS(I)+C(P)+S(C(P)+Residuals"
      i.an$Df.Pooled[i.an$ID == "02E.0035"]<- i.an$Df.Pooled[i.an$ID == "02E.0034"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0035"]<-i.an$Sum.Sq.Pooled[i.an$ID == "02E.0034"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0035"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0035"] / i.an$Df.Pooled[i.an$ID == "02E.0035"]
      i.an$F.Value[i.an$ID == "02E.0035"]<-i.an$Mean.Sq[i.an$ID == "02E.0035"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0035"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0035"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0035"], i.an$Df[i.an$ID == "02E.0035"], i.an$Df.Pooled[i.an$ID == "02E.0035"], lower.tail = F)
      #        i.an$P.Value.Lower.Tail[i.an$ID == "02E.0035"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0035"], i.an$Df[i.an$ID == "02E.0035"], i.an$Df.Pooled[i.an$ID == "02E.0035"], lower.tail = T)

      # 	If varCC(P) <= 0 and varS(CC(P)) > 0, then F = (MSP / MSpooled), where MSpooled = (SSCC(P) + SS_S(CC(P))) / (dfCC(P) + dfS(CC(P))) = SS1/df1
      f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "D",1:3])
      f.an$ID[f.an$Source.of.Variation == "D"]<-"01D.0001"
      f.an$Df.Formula[f.an$ID == "01D.0001"]<-"(p-1)"
      f.an$Calculated.From[f.an$ID == "01D.0001"]<-"sym.F"
      f.an$Mean.Sq[f.an$ID == "01D.0001"]<-f.an$Sum.Sq[f.an$ID == "01D.0001"] / f.an$Df[f.an$ID == "01D.0001"]
      f.an$Source.of.Variation[f.an$ID == "01D.0001"]<-"Place = P"
      f.an$Denominator[f.an$ID == "01D.0001"]<-"C(P)+S(C(P)"
      f.an$Df.Pooled[f.an$ID == "01D.0001"]<-c.an$Df[c.an$ID == "07DEn.0040"] +  + c.an$Df[c.an$ID == "08DEGn.0040"]
      f.an$Sum.Sq.Pooled[f.an$ID == "01D.0001"]<-c.an$Sum.Sq[c.an$ID == "07DEn.0040"] +  + c.an$Sum.Sq[c.an$ID == "08DEGn.0040"]
      f.an$Mean.Sq.Pooled[f.an$ID == "01D.0001"] <- f.an$Sum.Sq.Pooled[f.an$ID == "01D.0001"] / f.an$Df.Pooled[f.an$ID == "01D.0001"]
      f.an$F.Value[f.an$ID == "01D.0001"]<-f.an$Mean.Sq[f.an$ID == "01D.0001"] / f.an$Mean.Sq.Pooled[f.an$ID == "01D.0001"]
      f.an$P.Value.Upper.Tail[f.an$ID == "01D.0001"]<- stats::pf(f.an$F.Value[f.an$ID == "01D.0001"], f.an$Df[f.an$ID == "01D.0001"], f.an$Df.Pooled[f.an$ID == "01D.0001"], lower.tail = F)
      #        f.an$P.Value.Lower.Tail[f.an$ID == "01D.0001"]<- stats::pf(f.an$F.Value[f.an$ID == "01D.0001"], f.an$Df[f.an$ID == "01D.0001"], f.an$Df.Pooled[f.an$ID == "01D.0001"], lower.tail = T)

      # 	If varCC(P) <= 0 and varS(CC(P)) <= 0, then F = (MSP / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "D",1:3])
      f.an$ID[f.an$Source.of.Variation == "D"]<-"01D.0002"
      f.an$Df.Formula[f.an$ID == "01D.0002"]<-"(p-1)"
      f.an$Calculated.From[f.an$ID == "01D.0002"]<-"sym.F"
      f.an$Mean.Sq[f.an$ID == "01D.0002"]<-f.an$Sum.Sq[f.an$ID == "01D.0002"] / f.an$Df[f.an$ID == "01D.0002"]
      f.an$Source.of.Variation[f.an$ID == "01D.0002"]<-"Place = P"
      f.an$Denominator[f.an$ID == "01D.0002"]<-"C(P)+S(C(P)+Residuals"
      f.an$Df.Pooled[f.an$ID == "01D.0002"]<- f.an$Df.Pooled[f.an$ID == "01D.0001"] + f.an$Df[f.an$ID == "12R.0000"]
      f.an$Sum.Sq.Pooled[f.an$ID == "01D.0002"]<- f.an$Sum.Sq.Pooled[f.an$ID == "01D.0001"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      f.an$Mean.Sq.Pooled[f.an$ID == "01D.0002"] <- f.an$Sum.Sq.Pooled[f.an$ID == "01D.0002"] / f.an$Df.Pooled[f.an$ID == "01D.0002"]
      f.an$F.Value[f.an$ID == "01D.0002"]<-f.an$Mean.Sq[f.an$ID == "01D.0002"] / f.an$Mean.Sq.Pooled[f.an$ID == "01D.0002"]
      f.an$P.Value.Upper.Tail[f.an$ID == "01D.0002"]<- stats::pf(f.an$F.Value[f.an$ID == "01D.0002"], f.an$Df[f.an$ID == "01D.0002"], f.an$Df.Pooled[f.an$ID == "01D.0002"], lower.tail = F)
      #        f.an$P.Value.Lower.Tail[f.an$ID == "01D.0002"]<- stats::pf(f.an$F.Value[f.an$ID == "01D.0002"], f.an$Df[f.an$ID == "01D.0002"], f.an$Df.Pooled[f.an$ID == "01D.0002"], lower.tail = T)
    }
  }



  ###------------------ INTERACTIONS with TWO spatial factors

  if(aav.model == "DE") {
    # 	F = (MSCC(P) / MSRes)
    c.an<-dplyr::bind_rows(c.an,sym.NC[sym.NC$Source.of.Variation == "D:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "D:E"]<- "07DEn.0040"
    c.an$Df.Formula[c.an$ID == "07DEn.0040"]<-"p(c-1)"
    c.an$Calculated.From[c.an$ID == "07DEn.0040"]<-"sym.NC"
    c.an$Mean.Sq[c.an$ID == "07DEn.0040"] <- c.an$Sum.Sq[c.an$ID == "07DEn.0040"] / c.an$Df[c.an$ID == "07DEn.0040"]
    c.an$Source.of.Variation[c.an$ID == "07DEn.0040"]<-"Control location(P) = C(P)"
    c.an$Denominator[c.an$ID == "07DEn.0040"]<-"Residuals"
    c.an$F.Value[c.an$ID == "07DEn.0040"]<-c.an$Mean.Sq[c.an$ID == "07DEn.0040"] / f.an$Mean.Sq[f.an$ID == "12R.0000"]
    c.an$P.Value.Upper.Tail[c.an$ID == "07DEn.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "07DEn.0040"], c.an$Df[c.an$ID == "07DEn.0040"], f.an$Df[f.an$ID == "12R.0000"], lower.tail = F)

    # 	F = (MSPxIC / MSCC(P)) 03DE.0030
    dii.03DE.0030<-cbind(Source.of.Variation = "D:E",
                         sym.F[sym.F$Source.of.Variation == "D:E",2:3] - sym.C[sym.C$Source.of.Variation == "D:E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.03DE.0030)
    i.an$ID[i.an$Source.of.Variation == "D:E"]<-"03DE.0030"
    i.an$Df.Formula[i.an$ID == "03DE.0030"]<-"(p-1)(t-1)"
    i.an$Mean.Sq[i.an$ID == "03DE.0030"]<-i.an$Sum.Sq[i.an$ID == "03DE.0030"] / i.an$Df[i.an$ID == "03DE.0030"]
    i.an$Calculated.From[i.an$ID == "03DE.0030"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$ID == "03DE.0030"]<-"P x T = P x I"
    i.an$Denominator[i.an$ID == "03DE.0030"]<-"C(P)"
    i.an$F.Value[i.an$ID == "03DE.0030"]<-i.an$Mean.Sq[i.an$ID == "03DE.0030"] / c.an$Mean.Sq[c.an$ID == "07DEn.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "03DE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0030"], i.an$Df[i.an$ID == "03DE.0030"], c.an$Df[c.an$ID == "07DEn.0040"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "03DE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0030"], i.an$Df[i.an$ID == "03DE.0030"], c.an$Df[c.an$ID == "07DEn.0040"], lower.tail = T)

    # 	F = (MSIC / MSPxIC) 02E.0030
    dii.02E.0030<-cbind(Source.of.Variation = "E",
                        sym.F[sym.F$Source.of.Variation == "E",2:3] - sym.C[sym.C$Source.of.Variation == "E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0030"
    i.an$Df.Formula[i.an$ID == "02E.0030"]<-"(t-1)"
    i.an$Calculated.From[i.an$ID == "02E.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "02E.0030"]<-i.an$Sum.Sq[i.an$ID == "02E.0030"] / i.an$Df[i.an$ID == "02E.0030"]
    i.an$Source.of.Variation[i.an$ID == "02E.0030"]<-"Treatment = T = Impact vs Controls = I"
    i.an$Denominator[i.an$ID == "02E.0030"]<-"PxI"
    i.an$F.Value[i.an$ID == "02E.0030"]<-i.an$Mean.Sq[i.an$ID == "02E.0030"] / i.an$Mean.Sq[i.an$ID == "03DE.0030"]
    i.an$P.Value.Upper.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "03DE.0030"], lower.tail = F)
    #      i.an$P.Value.Lower.Tail[i.an$ID == "02E.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "02E.0030"], i.an$Df[i.an$ID == "03DE.0030"], lower.tail = T)

    # 	F = (MSP / MSCC(P)) 01D.0000
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "D",1:3])
    f.an$ID[f.an$Source.of.Variation == "D"]<-"01D.0000"
    f.an$Df.Formula[f.an$ID == "01D.0000"]<-"(p-1)"
    f.an$Calculated.From[f.an$ID == "01D.0000"]<-"sym.F"
    f.an$Mean.Sq[f.an$ID == "01D.0000"]<-f.an$Sum.Sq[f.an$ID == "01D.0000"] / f.an$Df[f.an$ID == "01D.0000"]
    f.an$Source.of.Variation[f.an$ID == "01D.0000"]<-"Place = P"
    f.an$Denominator[f.an$ID == "01D.0000"]<-"C(P)"
    f.an$F.Value[f.an$ID == "01D.0000"]<-f.an$Mean.Sq[f.an$ID == "01D.0000"] / c.an$Mean.Sq[c.an$ID == "07DEn.0040"]
    f.an$P.Value.Upper.Tail[f.an$ID == "01D.0000"]<- stats::pf(f.an$F.Value[f.an$ID == "01D.0000"], f.an$Df[f.an$ID == "01D.0000"], c.an$Df[c.an$ID == "07DEn.0040"], lower.tail = F)
    #      f.an$P.Value.Lower.Tail[f.an$ID == "01D.0000"]<- stats::pf(f.an$F.Value[f.an$ID == "01D.0000"], f.an$Df[f.an$ID == "01D.0000"], c.an$Df[c.an$ID == "07DEn.0040"], lower.tail = T)

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	If varCC(P) <= 0, then F = (MSPxIC / MSpooled), where MSpooled = (SSCC(P) + SSRes) / (dfCC(P) + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.03DE.0030)
      i.an$ID[i.an$Source.of.Variation == "D:E"]<-"03DE.0031"
      i.an$Df.Formula[i.an$ID == "03DE.0031"]<-"(p-1)(t-1)"
      i.an$Mean.Sq[i.an$ID == "03DE.0031"]<-i.an$Sum.Sq[i.an$ID == "03DE.0031"] / i.an$Df[i.an$ID == "03DE.0031"]
      i.an$Calculated.From[i.an$ID == "03DE.0031"]<-"sym.F-sym.C"
      i.an$Source.of.Variation[i.an$ID == "03DE.0031"]<-"P x T = P x I"
      i.an$Denominator[i.an$ID == "03DE.0031"]<-"C(P)+Residuals"
      i.an$Df.Pooled[i.an$ID == "03DE.0031"]<-c.an$Df[c.an$ID == "07DEn.0040"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "03DE.0031"]<-c.an$Sum.Sq[c.an$ID == "07DEn.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "03DE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03DE.0031"] / i.an$Df.Pooled[i.an$ID == "03DE.0031"]
      i.an$F.Value[i.an$ID == "03DE.0031"]<-i.an$Mean.Sq[i.an$ID == "03DE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03DE.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "03DE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0031"], i.an$Df[i.an$ID == "03DE.0031"], i.an$Df.Pooled[i.an$ID == "03DE.0031"], lower.tail = F)
      #        i.an$P.Value.Lower.Tail[i.an$ID == "03DE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03DE.0031"], i.an$Df[i.an$ID == "03DE.0031"], i.an$Df.Pooled[i.an$ID == "03DE.0031"], lower.tail = T)

      # 	If varPxIC <= 0 and varCC(P) > 0, then F = (MS_IC / MSpooled), where MSpooled = (SSPxIC + SSCC(P)) / (dfPxIC + dfCC(P)) = SS1/df1
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0031"
      i.an$Df.Formula[i.an$ID == "02E.0031"]<-"(t-1)"
      i.an$Calculated.From[i.an$ID == "02E.0031"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0031"]<-i.an$Sum.Sq[i.an$ID == "02E.0031"] / i.an$Df[i.an$ID == "02E.0031"]
      i.an$Source.of.Variation[i.an$ID == "02E.0031"]<-"Treatment = T = Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0031"]<-"PxI+C(P)"
      i.an$Df.Pooled[i.an$ID == "02E.0031"]<-i.an$Df[i.an$ID == "03DE.0030"] + c.an$Df[c.an$ID == "07DEn.0040"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"]<-i.an$Sum.Sq[i.an$ID == "03DE.0030"] + c.an$Sum.Sq[c.an$ID == "07DEn.0040"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"] / i.an$Df.Pooled[i.an$ID == "02E.0031"]
      i.an$F.Value[i.an$ID == "02E.0031"]<-i.an$Mean.Sq[i.an$ID == "02E.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0031"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df.Pooled[i.an$ID == "02E.0031"], lower.tail = F)
      #        i.an$P.Value.Lower.Tail[i.an$ID == "02E.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0031"], i.an$Df[i.an$ID == "02E.0031"], i.an$Df.Pooled[i.an$ID == "02E.0031"], lower.tail = T)

      # 	If varPxIC <= 0 and varCC(P) <= 0, then F = (MSIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      i.an<-dplyr::bind_rows(i.an,dii.02E.0030)
      i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0032"
      i.an$Df.Formula[i.an$ID == "02E.0032"]<-"(t-1)"
      i.an$Calculated.From[i.an$ID == "02E.0032"]<-"sym.F-sym.C"
      i.an$Mean.Sq[i.an$ID == "02E.0032"]<-i.an$Sum.Sq[i.an$ID == "02E.0032"] / i.an$Df[i.an$ID == "02E.0032"]
      i.an$Source.of.Variation[i.an$ID == "02E.0032"]<-"Treatment = T = Impact vs Controls = I"
      i.an$Denominator[i.an$ID == "02E.0032"]<-"PxI+C(P)+Residuals"
      i.an$Df.Pooled[i.an$ID == "02E.0032"]<- i.an$Df.Pooled[i.an$ID == "02E.0031"] + f.an$Df[f.an$ID == "12R.0000"]
      i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"]<-i.an$Sum.Sq.Pooled[i.an$ID == "02E.0031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "02E.0032"] / i.an$Df.Pooled[i.an$ID == "02E.0032"]
      i.an$F.Value[i.an$ID == "02E.0032"]<-i.an$Mean.Sq[i.an$ID == "02E.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "02E.0032"]
      i.an$P.Value.Upper.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = F)
      #        i.an$P.Value.Lower.Tail[i.an$ID == "02E.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "02E.0032"], i.an$Df[i.an$ID == "02E.0032"], i.an$Df.Pooled[i.an$ID == "02E.0032"], lower.tail = T)

      # 	If varCC(P) <= 0, then F = (MSP / MSpooled), where MSpooled = (SSCC(P) + SSRes) / (dfCC(P) + dfRes)
      f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "D",1:3])
      f.an$ID[f.an$Source.of.Variation == "D"]<-"01D.0001"
      f.an$Df.Formula[f.an$ID == "01D.0001"]<-"(p-1)"
      f.an$Calculated.From[f.an$ID == "01D.0001"]<-"sym.F"
      f.an$Mean.Sq[f.an$ID == "01D.0001"]<-f.an$Sum.Sq[f.an$ID == "01D.0001"] / f.an$Df[f.an$ID == "01D.0001"]
      f.an$Source.of.Variation[f.an$ID == "01D.0001"]<-"Place = P"
      f.an$Denominator[f.an$ID == "01D.0001"]<-"C(P)+Residuals"
      f.an$Df.Pooled[f.an$ID == "01D.0001"]<-c.an$Df[c.an$ID == "07DEn.0040"] + f.an$Df[f.an$ID == "12R.0000"]
      f.an$Sum.Sq.Pooled[f.an$ID == "01D.0001"]<-c.an$Sum.Sq[c.an$ID == "07DEn.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      f.an$Mean.Sq.Pooled[f.an$ID == "01D.0001"] <- f.an$Sum.Sq.Pooled[f.an$ID == "01D.0001"] / f.an$Df.Pooled[f.an$ID == "01D.0001"]
      f.an$F.Value[f.an$ID == "01D.0001"]<-f.an$Mean.Sq[f.an$ID == "01D.0001"] / f.an$Mean.Sq.Pooled[f.an$ID == "01D.0001"]
      f.an$P.Value.Upper.Tail[f.an$ID == "01D.0001"]<- stats::pf(f.an$F.Value[f.an$ID == "01D.0001"], f.an$Df[f.an$ID == "01D.0001"], f.an$Df.Pooled[f.an$ID == "01D.0001"], lower.tail = F)
      #        f.an$P.Value.Lower.Tail[f.an$ID == "01D.0001"]<- stats::pf(f.an$F.Value[f.an$ID == "01D.0001"], f.an$Df[f.an$ID == "01D.0001"], f.an$Df.Pooled[f.an$ID == "01D.0001"], lower.tail = T)
    }
  }


  ###----------------------------------------------------------------------- End

  ###----------------------------------------- #################################
  ###----------------------------------------- ### Table: Asymmetrical ANOVA ###
  ###----------------------------------------- #################################

  aci<-rbind(f.an,c.an,i.an)
  aci<-dplyr::arrange(aci,ID)


  ###----------------------------------------------------------- TWO-TAILED test
  if(aci.tailed.test == 2) {

  ###------------------------------------------------------------ Significant
  ## The F-ratio for some terms of the aav.models BEGC, CEG, DEG
  ## should be made by pooling procedures. Otherwise there is no test!

  # agora considerando P/2 para os casos de 2-tail
  if(pooling != T | (pooling == T & p.aav < p.pooling)) {
    # se um deles for menor que p.aav; se um deles for menor que p.pooling; todos os outros/ todos maior que p.pooling
    aci$Significant<-ifelse(aci$P.Value.Upper.Tail <= (p.aav/2) & !is.na(aci$P.Value.Upper.Tail) |
             aci$P.Value.Lower.Tail <= (p.aav/2) & !is.na(aci$P.Value.Lower.Tail), "*",
           ifelse(aci$P.Value.Upper.Tail <= (p.pooling/2) & !is.na(aci$P.Value.Upper.Tail) |
                    aci$P.Value.Lower.Tail <= (p.pooling/2) & !is.na(aci$P.Value.Lower.Tail), "ns",
                  "ns (Can be Eliminated)"))
  }
  if(pooling == T) {
    if(p.aav > p.pooling)   stop("p.aav must be equal or less than p.pooling", "\n")
    if(p.aav == p.pooling) {
      aci$Significant<-ifelse(aci$P.Value.Upper.Tail <= (p.aav/2) & !is.na(aci$P.Value.Upper.Tail) |
               aci$P.Value.Lower.Tail <= (p.aav/2) & !is.na(aci$P.Value.Lower.Tail), "*",
             "ns (Can be Eliminated)")
    }
  }


  ###------------------------------------------------------------ Post.Hoc.Pooling

  if(aav.model == "CE") {
    # 	F = (MSTxCC) / MSRes)
    aci$Post.Hoc.Pooling[aci$ID == "10CE.0040"]<-"Cannot be pooled"

    # 	F = (MSCC) / MSTxCC)
    aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Cannot be pooled"

    if(pooling != T) {
      # 	F = (MSTxIC) / MSTxCC)
      aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]<-"Cannot be pooled"

      # 	F = (MSIC / MSCC)
      aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
    }

    if(pooling == T) {
      # 	F = (MSTxIC) / MSTxCC)
      if(aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]<-"Cannot be pooled"
      # 	If varTxCC <= 0, then F = (MSTxIC / MSpooled), where MSpooled = (SSTxCC + SSRes) / (dfTxCC) + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10CE.0031"]<-"Can be pooled"
      CE10<- c("10CE.0030","10CE.0031")
      CE10n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0031"])) == F)

      # 	F = (MSIC / MSCC)
      if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      # 	If varCC <= 0 and varTxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSTxIC) / (dfCC + dfTxIC) = SS1/df1
      if((aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == "10CE.0030"] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
      # 	If varCC <= 0 and varTxIC <= 0, then F = (MSIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == "10CE.0030"] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"
    }
  }


  if(aav.model == "BEC") {
    # 	F = (MST(P)xCC / MSRes)
    aci$Post.Hoc.Pooling[aci$ID == "10BEC.0040"]<-"Cannot be pooled"

    # 	F = (MS_PxCC / MS_T(P)xCC)
    aci$Post.Hoc.Pooling[aci$ID == "03BE.0040"]<-"Cannot be pooled"

    # 	F = (MSCC / MSPxCC)
    aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Cannot be pooled"

    if(pooling != T) {
      # 	F = (MST(P)xIC) / MST(P)xCC)
      aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]<-"Cannot be pooled"

      # 	F = (MSPxIC) / MSPxCC)
      aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]<-"Cannot be pooled"

      # 	F = (MSIC / MSCC)
      aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
    }

    if(pooling == T) {
      # 	F = (MST(P)xIC) / MST(P)xCC)
      if(aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]<-"Cannot be pooled"
      # 	If varT(P)xCC <= 0, then F = (MS_T(P)xIC / MSpooled), where MSpooled = (SS_T(P)xCC + SSRes) / (dfT(P)xCC) + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10BEC.0031"]<-"Can be pooled"
      BEC10<- c("10BEC.0030","10BEC.0031")
      BEC10n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0031"])) == F)

      # 	F = (MSPxIC) / MSPxCC)
      if(aci$P.Value.Upper.Tail[aci$ID == "03BE.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]<-"Cannot be pooled"
      # 	If varPxCC <= 0 and varT(P)xIC > 0, then F = (MSPxIC / MSpooled), where MSpooled = (SSPxCC + SST(P)xIC) / (dfPxCC) + dfT(P)xIC) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == "03BE.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0031"]<-"Can be pooled"
      # 	If varPxCC <= 0 and varT(P)xIC <= 0, then F = (MSPxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "03BE.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0032"]<-"Can be pooled"
      BE03<- c("03BE.0030","03BE.0031","03BE.0032")
      BE03n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0031"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0032"])) == F)

      # 	F = (MSIC / MSCC)
      if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      # 	If varCC <= 0 and varPxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSPxIC) / (dfCC + dfPxIC) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BE03[BE03n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BE03[BE03n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
      # 	If varCC <= 0 and varPxIC <= 0 and varT(P)xIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SS1 + SST(P)xIC) / (df1 + dfT(P)xIC) = SS2/df2
      if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BE03[BE03n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BE03[BE03n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"
      # 	If varCC <= 0 and varPxIC <= 0 and varT(P)xIC <= 0, then F = (MSIC / MSpooled), Where MSpooled = (SS2 + SSRes) / (df2 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BE03[BE03n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BE03[BE03n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0033"]<-"Can be pooled"
    }
  }


  if(aav.model == "CEG") {
    # 	F = (MSTxS(CC)) / MSRes)
    aci$Post.Hoc.Pooling[aci$ID == "11CEG.0040"]<-"Cannot be pooled"

    # 	F = (MSTxCC) / MSTxS(CC))
    aci$Post.Hoc.Pooling[aci$ID == "10CE.0040"]<-"Cannot be pooled"

    # 	F = (MSS(CC)) / MSTxS(CC))
    aci$Post.Hoc.Pooling[aci$ID == "05EG.0040"]<-"Cannot be pooled"

    # If the two terms cannot be eliminated (* or ns)
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Cannot be pooled"
    # 	If varS(CC) <= 0 and varTxCC > 0, then F = (MSCC / MSTxCC)
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Can be pooled"
    # 	If varTxCC <= 0, varS(CC) > 0, then F = (MSCC / MSS(CC))
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0041"]<-"Can be pooled"
    # 	If varTxCC <= 0, varS(CC) <= 0, then F = (MSCC / MSpooled),
    #Where MSpooled = (SSTxCC + SSS(CC) + SSTxS(CC)) / (dfTxCC + dfS(CC) + dfTxS(CC))
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0042"]<-"Can be pooled"
    E02<- c("02E.0040","02E.0041","02E.0042")
    E02n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]),
                  is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0041"]),
                  is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0042"])) == F)

    if(pooling != T) {
      # 	F = (MSTxS(IC) / MSTxS(CC))
      aci$Post.Hoc.Pooling[aci$ID == "11CEG.0030"]<-"Cannot be pooled"

      # 	F = (MSTxIC) / MSTxCC), or
      aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]<-"Cannot be pooled"

      # 	F = (MSS(IC) / MSS(CC)), or
      aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]<-"Cannot be pooled"

      # 	F = (MSIC / MSCC)
      aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
    }

    if(pooling == T) {
      # 	F = (MSTxS(IC) / MSTxS(CC))
      if(aci$P.Value.Upper.Tail[aci$ID == "11CEG.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "11CEG.0030"]<-"Cannot be pooled"
      # 	If varTxS(CC) <= 0, then F = (MSTxS(IC) / MSpooled), where MSpooled = (SSTxS(CC) + SSRes) / (dfTxS(CC) + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "11CEG.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "11CEG.0031"]<-"Can be pooled"
      CEG11<- c("11CEG.0030","11CEG.0031")
      CEG11n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "11CEG.0030"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "11CEG.0031"])) == F)

      # 	F = (MSTxIC) / MSTxCC), or
      if(aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]<-"Cannot be pooled"
      # 	If varTxCC <= 0 and varTxS(IC) > 0, then F = (MSTxIC / MSpooled), where MSpooled = (SSTxCC + SSTxS(IC)) / (dfTxCC) + dfTxS(IC)) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == CEG11[CEG11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "10CE.0031"]<-"Can be pooled"
      # 	If varTxCC <= 0 and varTxS(IC) <= 0, then F = (MSTxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == CEG11[CEG11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "10CE.0032"]<-"Can be pooled"
      CE10<- c("10CE.0030","10CE.0031","10CE.0032")
      CE10n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0031"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0032"])) == F)

      # 	F = (MSS(IC) / MSS(CC)), or
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]<-"Cannot be pooled"
      # 	If varS(CC) <= 0 and varTxS(IC) > 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SSS(CC) + SSTxS(IC)) / (dfS(CC) + dfTxS(IC)) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
# 17_01_25        aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == CEG11[CEG11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]<-"Can be pooled"

      # 	If varS(CC) <= 0 and varTxS(IC) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
# 17_01_25        aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == CEG11[CEG11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"]<-"Can be pooled"
      EG05<- c("05EG.0030","05EG.0031","05EG.0032")
      EG05n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"])) == F)

      # 	F = (MSIC / MSCC)
      if(aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      # If the two terms cannot be eliminated (* or ns)
      if((aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] < p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == CE10[CE10n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      # 	If varCC <= 0, varS(IC) <= 0 and varTxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSTxIC) / (dfCC + dfTxIC)
      if((aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == CE10[CE10n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
      # 	If varCC <= 0, varTxIC <= 0, varS(IC) > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSS(IC)) / (dfCC + dfS(IC))
      if((aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] < p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == CE10[CE10n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"
      # 	If varCC <= 0, varTxIC <= 0, varS(IC) <= 0 and varTxS(IC) > 0, then F = (MSIC / MSpooled),
      # Where MSpooled = (SSCC + SSTxIC + SSS(IC) + SSTxS(IC)) / (dfCC + dfTxIC + dfS(IC) + dfTxS(IC)) = SS1/df1
      if((aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == CE10[CE10n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == CEG11[CEG11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0033"]<-"Can be pooled"
      # 	If varCC <= 0, varTxIC <= 0, varS(IC) <= 0 and varTxS(IC) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if((aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == CE10[CE10n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == CEG11[CEG11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0034"]<-"Can be pooled"
    }
  }



  if(aav.model == "BEGC") {
    # 	F = (MS_T(P)xS(CC) / MSRes)
    aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0040"]<-"Cannot be pooled"

    # 	F = (MST(P)xCC / MST(P)xS(CC))
    aci$Post.Hoc.Pooling[aci$ID == "10BEC.0040"]<-"Cannot be pooled"

    # 	F = (MSPxS(CC) / MST(P)xS(CC))
    aci$Post.Hoc.Pooling[aci$ID == "09BEG.0040"]<-"Cannot be pooled"

    # 	F = (MSS(CC) / MSPxS(CC))
    aci$Post.Hoc.Pooling[aci$ID == "05EG.0040"]<-"Cannot be pooled"

    # If the two terms cannot be eliminated (* or ns)
    if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] < p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0040"]<-"Cannot be pooled"
    # 	If varPxS(CC) <= 0 and varT(P)xCC > 0), then F = (MSPxCC / MST(P)xCC)
    if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0040"]<-"Can be pooled"
    # 	If varT(P)xCC <= 0 and varPxS(CC) > 0, then F = (MSPxCC / MSPxS(CC))
    if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] < p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0041"]<-"Can be pooled"
    # 	If varT(P)xCC <= 0, varPxS(CC) <= 0, then F = (MSPxCC / MSpooled),
    # where MSpooled = (SST(P)x CC + SSPxS(CC) + SST(P)xS(CC)) / (dfT(P)x CC + dfPxS(CC) + dfT(P)xS(CC))
    if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0042"]<-"Can be pooled"
    BE03.40<- c("03BE.0040","03BE.0041","03BE.0042")
    BE03.40n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0040"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0041"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0042"])) == F)

    # If the two terms cannot be eliminated (* or ns)
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Cannot be pooled"
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] < p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Cannot be pooled"
    # 	If varS(CC) <= 0 and varPxCC > 0, then F = (MSCC / MSPxCC)
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Can be pooled"
    # 	If varPxCC <= 0 and varS(CC) > 0, then F = (MSCC / MSS(CC))
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0041"]<-"Can be pooled"
    # 	If varPxCC <= 0, varS(CC) <=, varPxS(CC) <= 0 and varT(P)xCC > 0, then F = (MSCC / MSpooled),
    # Where MSpooled = (SSPxCC + SSS(CC) + SST(P)xCC) / (dfPxCC + dfS(CC) + dfT(P)xCC)
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0042"]<-"Can be pooled"
    # 	If varPxCC <= 0, varS(CC) <=, varT(P)xCC <= 0 and varPxS(CC) > 0, then F = (MSCC / MSpooled),
    # Where MSpooled = (SSPxCC + SSS(CC) + SSPxS(CC)) / (dfPxCC + dfS(CC) + dfPxS(CC))
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] < p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0043"]<-"Can be pooled"
    # 	If varPxCC <= 0, varS(CC) <=, varT(P)xCC <= 0, varPxS(CC) <= 0, then F = (MSCC / MSpooled),
    # Where MSpooled = (SSPxCC + SSS(CC) + SST(P)xCC + SSPxS(CC) + SST(P)xS(CC)) / (dfPxCC + dfS(CC) + dfT(P)xCC + dfPxS(CC) + dfT(P)xS(CC))
    if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling &
       aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0044"]<-"Can be pooled"
    E02.40<- c("02E.0040","02E.0041","02E.0042","02E.0043","02E.0044")
    E02.40n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0041"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0042"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0043"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0044"])) == F)

    if(pooling != T) {
      # 	F = (MST(P)xS(IC) / MST(P)xS(CC)), or
      aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0030"]<-"Cannot be pooled"

      # 	F = (MST(P)xIC) / MST(P)xCC), or
      aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]<-"Cannot be pooled"

      # 	F = (MSPxS(IC) / MSPxS(CC)), or
      aci$Post.Hoc.Pooling[aci$ID == "09BEG.0030"]<-"Cannot be pooled"

      # 	F = (MSS(IC) / MSS(CC)), or
      aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]<-"Cannot be pooled"

      # 	F = (MSPxIC) / MSPxCC), or
      aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]<-"Cannot be pooled"

      # 	F = (MSIC /MSCC), or
      aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
    }

    if(pooling == T) {
      # 	F = (MST(P)xS(IC) / MST(P)xS(CC)), or
      if(aci$P.Value.Upper.Tail[aci$ID == "11BEGC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0030"]<-"Cannot be pooled"
      # 	If varT(P)xS(CC) <= 0, then F = (MST(P)xS(IC) / MSpooled), where MSpooled = (SST(P)xS(CC) + SSRes) / (dfT(P)xS(CC) + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "11BEGC.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0031"]<-"Can be pooled"
      BEGC11<- c("11BEGC.0030","11BEGC.0031")
      BEGC11n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0030"]),
                       is.na(aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0031"])) == F)

      # 	F = (MST(P)xIC) / MST(P)xCC), or
      if(aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]<-"Cannot be pooled"
      # 	If varT(P)xCC <= 0 and varT(P)xS(IC) > 0, then F = (MST(P)xIC / MSpooled),
      # where MSpooled = (SST(P)xCC + SST(P)xS(IC)) / (dfT(P)xCC) + dfT(P)xS(IC)) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "10BEC.0031"]<-"Can be pooled"
      # 	If varT(P)xCC <= 0 and varT(P)xS(IC) <= 0, then F = (MST(P)xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "10BEC.0032"]<-"Can be pooled"
      BEC10<- c("10BEC.0030","10BEC.0031","10BEC.0032")
      BEC10n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0031"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0032"])) == F)

      # 	F = (MSPxS(IC) / MSPxS(CC)), or
      if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "09BEG.0030"]<-"Cannot be pooled"
      # 	If varPxS(CC) <= 0 and varT(P)xS(IC) > 0, then F = (MSPxS(IC) / MSpooled),
      # where MSpooled = (SSPxS(CC) + SST(P)xS(IC)) / (dfPxS(CC) + dfT(P)xS(IC)) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "09BEG.0031"]<-"Can be pooled"
      # 	If varPxS(CC) <= 0 and varT(P)xS(IC) <= 0, then F = (MSPxS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "09BEG.0032"]<-"Can be pooled"
      BEG09<- c("09BEG.0030","09BEG.0031","09BEG.0032")
      BEG09n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "09BEG.0030"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "09BEG.0031"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "09BEG.0032"])) == F)

      # 	F = (MSS(IC) / MSS(CC)), or
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]<-"Cannot be pooled"
      # 	If varS(CC) <= 0 and varPxS(IC) > 0, then F = (MSS(IC) / MSpooled),
      # where MSpooled = (SSS(CC) + SSPxS(IC)) / (dfS(CC) + dfPxS(IC)) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]<-"Can be pooled"
      # 	If varS(CC) <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) > 0, then F = (MSS(IC) / MSpooled),
      # where MSpooled = (SS1 + SST(P)xS(IC)) / (df1 + dfT(P)xS(IC)) = SS2/df2
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"]<-"Can be pooled"
      # 	If varS(CC) <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS2 + SSRes) / (df2 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "05EG.0033"]<-"Can be pooled"
      EG05<- c("05EG.0030","05EG.0031","05EG.0032","05EG.0033")
      EG05n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0033"])) == F)

      # 	F = (MSPxIC) / MSPxCC), or
      if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]<-"Cannot be pooled"
      # If the two terms cannot be eliminated (* or ns)
      if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] < p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]<-"Cannot be pooled"
      # 	If varPxCC <= 0, varPxS(IC) <= 0 and varT(P)xIC > 0), then F = (MSPxIC / MSpooled),
      # where MSpooled = (SSPxCC + SST(P)xIC) / (dfPxCC) + dfT(P)xIC)
      if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0031"]<-"Can be pooled"
      # 	If varPxCC <= 0, varT(P)xIC <= 0 and varPxS(IC) > 0, then F = (MSPxIC / MSpooled),
      # where MSpooled = (SSPxCC + SSPxS(IC)) / (dfPxCC) + dfPxS(IC))
      if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] < p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0032"]<-"Can be pooled"
      # 	If varPxCC <= 0, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) > 0,
      # then F = (MSPxIC / MSpooled), where MSpooled = (SSPxCC + SST(P)xIC + SSPxS(IC) + SST(P)xS(IC)) / (dfTxCC) + dfT(P)xIC + dfPxS(IC) + dfT(P)xS(IC)) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0033"]<-"Can be pooled"
      # 	If varPxCC <= 0, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) <= 0,
      # then F = (MSPxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0034"]<-"Can be pooled"
      BE03.30<- c("03BE.0030","03BE.0031","03BE.0032","03BE.0033","03BE.0034")
      BE03.30n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]),
                        is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0031"]),
                        is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0032"]),
                        is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0033"]),
                        is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0034"])) == F)

      # 	F = (MSIC /MSCC), or
      if(aci$P.Value.Upper.Tail[aci$ID == E02.40[E02.40n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      # If the two terms cannot be eliminated (* or ns)
      if(aci$P.Value.Upper.Tail[aci$ID == E02.40[E02.40n]] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] < p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BE03.30[BE03.30n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      if(aci$P.Value.Upper.Tail[aci$ID == E02.40[E02.40n]] >= p.pooling &
         (aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] < p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] < p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      # 	If varCC <= 0, varS(IC) <= 0 and varPxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSPxIC) / (dfCC + dfPxIC)
      if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BE03.30[BE03.30n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
      # 	If varCC <= 0, varPxIC <= 0 and varS(IC) > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSS(IC)) / (dfCC + dfS(IC))
      if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] < p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"
      # 	If varCC <= 0, varPxIC <= 0, varS(IC) <=, varPxS(IC) <= 0 and varT(P)xIC > 0, then F = (MSIC / MSpooled), Where MSpooled = (SSCC + SSPxIC + SSS(IC) + SST(P)xIC) / (dfCC + dfPxIC + dfS(IC) + dfT(P)xIC)
      if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] < p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0033"]<-"Can be pooled"
      # 	If varCC <= 0, varPxIC <= 0, varS(IC) <=, varT(P)xIC <= 0 and varPxS(IC) > 0, then F = (MSIC / MSpooled), Where MSpooled = (SSCC + SSPxIC + SSS(IC) + SSPxS(IC)) / (dfCC + dfPxIC + dfS(IC) + dfPxS(IC))
      if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0034"]<-"Can be pooled"
      # 	If varCC <= 0, varPxIC <= 0, varS(IC) <=, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) > 0, then F = (MSIC / MSpooled), Where MSpooled = (SSCC + SSPxIC + SSS(IC) + SST(P)xIC + SSPxS(IC) + SST(P)xS(IC)) / (dfCC + dfPxIC + dfS(IC) + dfT(P)xIC + dfPxS(IC) + dfT(P)xS(IC)) = SS1/df1
      if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling |
          aci$P.Value.Lower.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0035"]<-"Can be pooled"
      # 	If varCC <= 0, varPxIC <= 0, varS(IC) <= 0, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) <= 0, then F = (MSIC / MSpooled), Where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
         (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling &
          aci$P.Value.Lower.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0036"]<-"Can be pooled"
    }
  }
  }



    ###-------------------------------------------------------- ONE-TAILED test
    if(aci.tailed.test == 1) {

      ###------------------------------------------------------------ Significant
      ## The F-ratio for some terms of the aav.models BEGC, CEG, DEG
      ## should be made by pooling procedures. Otherwise there is no test!

      if(pooling != T | (pooling == T & p.aav < p.pooling)) {
        # se um deles for menor que p.aav; se um deles for menor que p.pooling; todos os outros/ todos maior que p.pooling
        aci$Significant<-ifelse(aci$P.Value.Upper.Tail <= p.aav & !is.na(aci$P.Value.Upper.Tail), "*",
                                ifelse(aci$P.Value.Upper.Tail <= p.pooling & !is.na(aci$P.Value.Upper.Tail), "ns",
                                       "ns (Can be Eliminated)"))
      }
      if(pooling == T) {
        if(p.aav > p.pooling)   stop("p.aav must be equal or less than p.pooling", "\n")
        if(p.aav == p.pooling) {
          aci$Significant<-ifelse(aci$P.Value.Upper.Tail <= p.aav & !is.na(aci$P.Value.Upper.Tail), "*",
                                  "ns (Can be Eliminated)")
        }
      }

      ###------------------------------------------------------------ Post.Hoc.Pooling

      if(aav.model == "CE") {
        # 	F = (MSTxCC) / MSRes)
        aci$Post.Hoc.Pooling[aci$ID == "10CE.0040"]<-"Cannot be pooled"

        # 	F = (MSCC) / MSTxCC)
        aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Cannot be pooled"

        if(pooling != T) {
          # 	F = (MSTxIC) / MSTxCC)
          aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]<-"Cannot be pooled"

          # 	F = (MSIC / MSCC)
          aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
        }

        if(pooling == T) {
          # 	F = (MSTxIC) / MSTxCC)
          if(aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]<-"Cannot be pooled"
          # 	If varTxCC <= 0, then F = (MSTxIC / MSpooled), where MSpooled = (SSTxCC + SSRes) / (dfTxCC) + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10CE.0031"]<-"Can be pooled"
          CE10<- c("10CE.0030","10CE.0031")
          CE10n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0031"])) == F)

          # 	F = (MSIC / MSCC)
          if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
          # 	If varCC <= 0 and varTxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSTxIC) / (dfCC + dfTxIC) = SS1/df1
          if((aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] >= p.pooling) &
             aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
          # 	If varCC <= 0 and varTxIC <= 0, then F = (MSIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"
        }
      }


      if(aav.model == "BEC") {
        # 	F = (MST(P)xCC / MSRes)
        aci$Post.Hoc.Pooling[aci$ID == "10BEC.0040"]<-"Cannot be pooled"

        # 	F = (MS_PxCC / MS_T(P)xCC)
        aci$Post.Hoc.Pooling[aci$ID == "03BE.0040"]<-"Cannot be pooled"

        # 	F = (MSCC / MSPxCC)
        aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Cannot be pooled"

        if(pooling != T) {
          # 	F = (MST(P)xIC) / MST(P)xCC)
          aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]<-"Cannot be pooled"

          # 	F = (MSPxIC) / MSPxCC)
          aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]<-"Cannot be pooled"

          # 	F = (MSIC / MSCC)
          aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
        }

        if(pooling == T) {
          # 	F = (MST(P)xIC) / MST(P)xCC)
          if(aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]<-"Cannot be pooled"
          # 	If varT(P)xCC <= 0, then F = (MS_T(P)xIC / MSpooled), where MSpooled = (SS_T(P)xCC + SSRes) / (dfT(P)xCC) + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10BEC.0031"]<-"Can be pooled"
          BEC10<- c("10BEC.0030","10BEC.0031")
          BEC10n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]),
                          is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0031"])) == F)

          # 	F = (MSPxIC) / MSPxCC)
          if(aci$P.Value.Upper.Tail[aci$ID == "03BE.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]<-"Cannot be pooled"
          # 	If varPxCC <= 0 and varT(P)xIC > 0, then F = (MSPxIC / MSpooled), where MSpooled = (SSPxCC + SST(P)xIC) / (dfPxCC) + dfT(P)xIC) = SS1/df1
          if(aci$P.Value.Upper.Tail[aci$ID == "03BE.0040"] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0031"]<-"Can be pooled"
          # 	If varPxCC <= 0 and varT(P)xIC <= 0, then F = (MSPxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "03BE.0040"] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0032"]<-"Can be pooled"
          BE03<- c("03BE.0030","03BE.0031","03BE.0032")
          BE03n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0031"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0032"])) == F)

          # 	F = (MSIC / MSCC)
          if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
          # 	If varCC <= 0 and varPxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSPxIC) / (dfCC + dfPxIC) = SS1/df1
          if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == BE03[BE03n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
          # 	If varCC <= 0 and varPxIC <= 0 and varT(P)xIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SS1 + SST(P)xIC) / (df1 + dfT(P)xIC) = SS2/df2
          if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == BE03[BE03n]] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"
          # 	If varCC <= 0 and varPxIC <= 0 and varT(P)xIC <= 0, then F = (MSIC / MSpooled), Where MSpooled = (SS2 + SSRes) / (df2 + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "02E.0040"] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == BE03[BE03n]] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0033"]<-"Can be pooled"
        }
      }


      if(aav.model == "CEG") {
        # 	F = (MSTxS(CC)) / MSRes)
        aci$Post.Hoc.Pooling[aci$ID == "11CEG.0040"]<-"Cannot be pooled"

        # 	F = (MSTxCC) / MSTxS(CC))
        aci$Post.Hoc.Pooling[aci$ID == "10CE.0040"]<-"Cannot be pooled"

        # 	F = (MSS(CC)) / MSTxS(CC))
        aci$Post.Hoc.Pooling[aci$ID == "05EG.0040"]<-"Cannot be pooled"

        # If the two terms cannot be eliminated (* or ns)
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Cannot be pooled"
        # 	If varS(CC) <= 0 and varTxCC > 0, then F = (MSCC / MSTxCC)
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Can be pooled"
        # 	If varTxCC <= 0, varS(CC) > 0, then F = (MSCC / MSS(CC))
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0041"]<-"Can be pooled"
        # 	If varTxCC <= 0, varS(CC) <= 0, then F = (MSCC / MSpooled),
        #Where MSpooled = (SSTxCC + SSS(CC) + SSTxS(CC)) / (dfTxCC + dfS(CC) + dfTxS(CC))
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0042"]<-"Can be pooled"
        E02<- c("02E.0040","02E.0041","02E.0042")
        E02n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0041"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0042"])) == F)

        if(pooling != T) {
          # 	F = (MSTxS(IC) / MSTxS(CC))
          aci$Post.Hoc.Pooling[aci$ID == "11CEG.0030"]<-"Cannot be pooled"

          # 	F = (MSTxIC) / MSTxCC), or
          aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]<-"Cannot be pooled"

          # 	F = (MSS(IC) / MSS(CC)), or
          aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]<-"Cannot be pooled"

          # 	F = (MSIC / MSCC)
          aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
        }

        if(pooling == T) {
          # 	F = (MSTxS(IC) / MSTxS(CC))
          if(aci$P.Value.Upper.Tail[aci$ID == "11CEG.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "11CEG.0030"]<-"Cannot be pooled"
          # 	If varTxS(CC) <= 0, then F = (MSTxS(IC) / MSpooled), where MSpooled = (SSTxS(CC) + SSRes) / (dfTxS(CC) + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "11CEG.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "11CEG.0031"]<-"Can be pooled"
          CEG11<- c("11CEG.0030","11CEG.0031")
          CEG11n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "11CEG.0030"]),
                          is.na(aci$Post.Hoc.Pooling[aci$ID == "11CEG.0031"])) == F)

          # 	F = (MSTxIC) / MSTxCC), or
          if(aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]<-"Cannot be pooled"
          # 	If varTxCC <= 0 and varTxS(IC) > 0, then F = (MSTxIC / MSpooled), where MSpooled = (SSTxCC + SSTxS(IC)) / (dfTxCC) + dfTxS(IC)) = SS1/df1
          if(aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10CE.0031"]<-"Can be pooled"
          # 	If varTxCC <= 0 and varTxS(IC) <= 0, then F = (MSTxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10CE.0032"]<-"Can be pooled"
          CE10<- c("10CE.0030","10CE.0031","10CE.0032")
          CE10n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0031"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0032"])) == F)

          # 	F = (MSS(IC) / MSS(CC)), or
          if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]<-"Cannot be pooled"
          # 	If varS(CC) <= 0 and varTxS(IC) > 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SSS(CC) + SSTxS(IC)) / (dfS(CC) + dfTxS(IC)) = SS1/df1
          if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
             # 17_01_25        aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]<-"Can be pooled"

          # 	If varS(CC) <= 0 and varTxS(IC) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
             # 17_01_25        aci$P.Value.Upper.Tail[aci$ID == "10CE.0040"] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"]<-"Can be pooled"
          EG05<- c("05EG.0030","05EG.0031","05EG.0032")
          EG05n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"])) == F)

          # 	F = (MSIC / MSCC)
          if(aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
          # If the two terms cannot be eliminated (* or ns)
          if((aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] >= p.pooling) &
             aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] < p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
          # 	If varCC <= 0, varS(IC) <= 0 and varTxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSTxIC) / (dfCC + dfTxIC)
          if((aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] >= p.pooling) &
             aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
          # 	If varCC <= 0, varTxIC <= 0, varS(IC) > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSS(IC)) / (dfCC + dfS(IC))
          if((aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] >= p.pooling) &
             aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] < p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"
          # 	If varCC <= 0, varTxIC <= 0, varS(IC) <= 0 and varTxS(IC) > 0, then F = (MSIC / MSpooled),
          # Where MSpooled = (SSCC + SSTxIC + SSS(IC) + SSTxS(IC)) / (dfCC + dfTxIC + dfS(IC) + dfTxS(IC)) = SS1/df1
          if((aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] >= p.pooling) &
             aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] >= p.pooling &
             aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0033"]<-"Can be pooled"
          # 	If varCC <= 0, varTxIC <= 0, varS(IC) <= 0 and varTxS(IC) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
          if((aci$P.Value.Upper.Tail[aci$ID == E02[E02n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == CE10[CE10n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == CEG11[CEG11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0034"]<-"Can be pooled"
        }
      }



      if(aav.model == "BEGC") {
        # 	F = (MS_T(P)xS(CC) / MSRes)
        aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0040"]<-"Cannot be pooled"

        # 	F = (MST(P)xCC / MST(P)xS(CC))
        aci$Post.Hoc.Pooling[aci$ID == "10BEC.0040"]<-"Cannot be pooled"

        # 	F = (MSPxS(CC) / MST(P)xS(CC))
        aci$Post.Hoc.Pooling[aci$ID == "09BEG.0040"]<-"Cannot be pooled"

        # 	F = (MSS(CC) / MSPxS(CC))
        aci$Post.Hoc.Pooling[aci$ID == "05EG.0040"]<-"Cannot be pooled"

        # If the two terms cannot be eliminated (* or ns)
        if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] < p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0040"]<-"Cannot be pooled"
        # 	If varPxS(CC) <= 0 and varT(P)xCC > 0), then F = (MSPxCC / MST(P)xCC)
        if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0040"]<-"Can be pooled"
        # 	If varT(P)xCC <= 0 and varPxS(CC) > 0, then F = (MSPxCC / MSPxS(CC))
        if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] < p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0041"]<-"Can be pooled"
        # 	If varT(P)xCC <= 0, varPxS(CC) <= 0, then F = (MSPxCC / MSpooled),
        # where MSpooled = (SST(P)x CC + SSPxS(CC) + SST(P)xS(CC)) / (dfT(P)x CC + dfPxS(CC) + dfT(P)xS(CC))
        if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0042"]<-"Can be pooled"
        BE03.40<- c("03BE.0040","03BE.0041","03BE.0042")
        BE03.40n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0040"]),
                          is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0041"]),
                          is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0042"])) == F)

        # If the two terms cannot be eliminated (* or ns)
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Cannot be pooled"
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] < p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Cannot be pooled"
        # 	If varS(CC) <= 0 and varPxCC > 0, then F = (MSCC / MSPxCC)
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]<-"Can be pooled"
        # 	If varPxCC <= 0 and varS(CC) > 0, then F = (MSCC / MSS(CC))
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0041"]<-"Can be pooled"
        # 	If varPxCC <= 0, varS(CC) <=, varPxS(CC) <= 0 and varT(P)xCC > 0, then F = (MSCC / MSpooled),
        # Where MSpooled = (SSPxCC + SSS(CC) + SST(P)xCC) / (dfPxCC + dfS(CC) + dfT(P)xCC)
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0042"]<-"Can be pooled"
        # 	If varPxCC <= 0, varS(CC) <=, varT(P)xCC <= 0 and varPxS(CC) > 0, then F = (MSCC / MSpooled),
        # Where MSpooled = (SSPxCC + SSS(CC) + SSPxS(CC)) / (dfPxCC + dfS(CC) + dfPxS(CC))
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] < p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0043"]<-"Can be pooled"
        # 	If varPxCC <= 0, varS(CC) <=, varT(P)xCC <= 0, varPxS(CC) <= 0, then F = (MSCC / MSpooled),
        # Where MSpooled = (SSPxCC + SSS(CC) + SST(P)xCC + SSPxS(CC) + SST(P)xS(CC)) / (dfPxCC + dfS(CC) + dfT(P)xCC + dfPxS(CC) + dfT(P)xS(CC))
        if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling &
           aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0044"]<-"Can be pooled"
        E02.40<- c("02E.0040","02E.0041","02E.0042","02E.0043","02E.0044")
        E02.40n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0041"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0042"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0043"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0044"])) == F)

        if(pooling != T) {
          # 	F = (MST(P)xS(IC) / MST(P)xS(CC)), or
          aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0030"]<-"Cannot be pooled"

          # 	F = (MST(P)xIC) / MST(P)xCC), or
          aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]<-"Cannot be pooled"

          # 	F = (MSPxS(IC) / MSPxS(CC)), or
          aci$Post.Hoc.Pooling[aci$ID == "09BEG.0030"]<-"Cannot be pooled"

          # 	F = (MSS(IC) / MSS(CC)), or
          aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]<-"Cannot be pooled"

          # 	F = (MSPxIC) / MSPxCC), or
          aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]<-"Cannot be pooled"

          # 	F = (MSIC /MSCC), or
          aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
        }

        if(pooling == T) {
          # 	F = (MST(P)xS(IC) / MST(P)xS(CC)), or
          if(aci$P.Value.Upper.Tail[aci$ID == "11BEGC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0030"]<-"Cannot be pooled"
          # 	If varT(P)xS(CC) <= 0, then F = (MST(P)xS(IC) / MSpooled), where MSpooled = (SST(P)xS(CC) + SSRes) / (dfT(P)xS(CC) + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "11BEGC.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0031"]<-"Can be pooled"
          BEGC11<- c("11BEGC.0030","11BEGC.0031")
          BEGC11n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0030"]),
                           is.na(aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0031"])) == F)

          # 	F = (MST(P)xIC) / MST(P)xCC), or
          if(aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]<-"Cannot be pooled"
          # 	If varT(P)xCC <= 0 and varT(P)xS(IC) > 0, then F = (MST(P)xIC / MSpooled),
          # where MSpooled = (SST(P)xCC + SST(P)xS(IC)) / (dfT(P)xCC) + dfT(P)xS(IC)) = SS1/df1
          if(aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "10BEC.0031"]<-"Can be pooled"
          # 	If varT(P)xCC <= 0 and varT(P)xS(IC) <= 0, then F = (MST(P)xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "10BEC.0040"] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "10BEC.0032"]<-"Can be pooled"
          BEC10<- c("10BEC.0030","10BEC.0031","10BEC.0032")
          BEC10n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]),
                          is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0031"]),
                          is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0032"])) == F)

          # 	F = (MSPxS(IC) / MSPxS(CC)), or
          if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "09BEG.0030"]<-"Cannot be pooled"
          # 	If varPxS(CC) <= 0 and varT(P)xS(IC) > 0, then F = (MSPxS(IC) / MSpooled),
          # where MSpooled = (SSPxS(CC) + SST(P)xS(IC)) / (dfPxS(CC) + dfT(P)xS(IC)) = SS1/df1
          if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "09BEG.0031"]<-"Can be pooled"
          # 	If varPxS(CC) <= 0 and varT(P)xS(IC) <= 0, then F = (MSPxS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "09BEG.0040"] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "09BEG.0032"]<-"Can be pooled"
          BEG09<- c("09BEG.0030","09BEG.0031","09BEG.0032")
          BEG09n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "09BEG.0030"]),
                          is.na(aci$Post.Hoc.Pooling[aci$ID == "09BEG.0031"]),
                          is.na(aci$Post.Hoc.Pooling[aci$ID == "09BEG.0032"])) == F)

          # 	F = (MSS(IC) / MSS(CC)), or
          if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]<-"Cannot be pooled"
          # 	If varS(CC) <= 0 and varPxS(IC) > 0, then F = (MSS(IC) / MSpooled),
          # where MSpooled = (SSS(CC) + SSPxS(IC)) / (dfS(CC) + dfPxS(IC)) = SS1/df1
          if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]<-"Can be pooled"
          # 	If varS(CC) <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) > 0, then F = (MSS(IC) / MSpooled),
          # where MSpooled = (SS1 + SST(P)xS(IC)) / (df1 + dfT(P)xS(IC)) = SS2/df2
          if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"]<-"Can be pooled"
          # 	If varS(CC) <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS2 + SSRes) / (df2 + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0040"] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "05EG.0033"]<-"Can be pooled"
          EG05<- c("05EG.0030","05EG.0031","05EG.0032","05EG.0033")
          EG05n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"]),
                         is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0033"])) == F)

          # 	F = (MSPxIC) / MSPxCC), or
          if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]<-"Cannot be pooled"
          # If the two terms cannot be eliminated (* or ns)
          if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] < p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]<-"Cannot be pooled"
          # 	If varPxCC <= 0, varPxS(IC) <= 0 and varT(P)xIC > 0), then F = (MSPxIC / MSpooled),
          # where MSpooled = (SSPxCC + SST(P)xIC) / (dfPxCC) + dfT(P)xIC)
          if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0031"]<-"Can be pooled"
          # 	If varPxCC <= 0, varT(P)xIC <= 0 and varPxS(IC) > 0, then F = (MSPxIC / MSpooled),
          # where MSpooled = (SSPxCC + SSPxS(IC)) / (dfPxCC) + dfPxS(IC))
          if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] < p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0032"]<-"Can be pooled"
          # 	If varPxCC <= 0, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) > 0,
          # then F = (MSPxIC / MSpooled), where MSpooled = (SSPxCC + SST(P)xIC + SSPxS(IC) + SST(P)xS(IC)) / (dfTxCC) + dfT(P)xIC + dfPxS(IC) + dfT(P)xS(IC)) = SS1/df1
          if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0033"]<-"Can be pooled"
          # 	If varPxCC <= 0, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) <= 0,
          # then F = (MSPxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
          if(aci$P.Value.Upper.Tail[aci$ID == BE03.40[BE03.40n]] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "03BE.0034"]<-"Can be pooled"
          BE03.30<- c("03BE.0030","03BE.0031","03BE.0032","03BE.0033","03BE.0034")
          BE03.30n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]),
                            is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0031"]),
                            is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0032"]),
                            is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0033"]),
                            is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0034"])) == F)

          # 	F = (MSIC /MSCC), or
          if(aci$P.Value.Upper.Tail[aci$ID == E02.40[E02.40n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
          # If the two terms cannot be eliminated (* or ns)
          if(aci$P.Value.Upper.Tail[aci$ID == E02.40[E02.40n]] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] < p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
          if(aci$P.Value.Upper.Tail[aci$ID == E02.40[E02.40n]] >= p.pooling &
             (aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
          # 	If varCC <= 0, varS(IC) <= 0 and varPxIC > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSPxIC) / (dfCC + dfPxIC)
          if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
          # 	If varCC <= 0, varPxIC <= 0 and varS(IC) > 0, then F = (MSIC / MSpooled), where MSpooled = (SSCC + SSS(IC)) / (dfCC + dfS(IC))
          if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] < p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"
          # 	If varCC <= 0, varPxIC <= 0, varS(IC) <=, varPxS(IC) <= 0 and varT(P)xIC > 0, then F = (MSIC / MSpooled), Where MSpooled = (SSCC + SSPxIC + SSS(IC) + SST(P)xIC) / (dfCC + dfPxIC + dfS(IC) + dfT(P)xIC)
          if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] < p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0033"]<-"Can be pooled"
          # 	If varCC <= 0, varPxIC <= 0, varS(IC) <=, varT(P)xIC <= 0 and varPxS(IC) > 0, then F = (MSIC / MSpooled), Where MSpooled = (SSCC + SSPxIC + SSS(IC) + SSPxS(IC)) / (dfCC + dfPxIC + dfS(IC) + dfPxS(IC))
          if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0034"]<-"Can be pooled"
          # 	If varCC <= 0, varPxIC <= 0, varS(IC) <=, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) > 0, then F = (MSIC / MSpooled), Where MSpooled = (SSCC + SSPxIC + SSS(IC) + SST(P)xIC + SSPxS(IC) + SST(P)xS(IC)) / (dfCC + dfPxIC + dfS(IC) + dfT(P)xIC + dfPxS(IC) + dfT(P)xS(IC)) = SS1/df1
          if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] < p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0035"]<-"Can be pooled"
          # 	If varCC <= 0, varPxIC <= 0, varS(IC) <= 0, varT(P)xIC <= 0, varPxS(IC) <= 0 and varT(P)xS(IC) <= 0, then F = (MSIC / MSpooled), Where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
          if((aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BE03.30[BE03.30n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEC10[BEC10n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEG09[BEG09n]] >= p.pooling) &
             (aci$P.Value.Upper.Tail[aci$ID == BEGC11[BEGC11n]] >= p.pooling)) aci$Post.Hoc.Pooling[aci$ID == "02E.0036"]<-"Can be pooled"
        }
      }
    }


###--------------------------------------------------------------------------
###--------------------------------------- Post-impact (Glasby) nao tem 2-tailed
  if(aav.model == "DE") {
    # 	F = (MSCC(P) / MSRes) 07DEn.0040
    aci$Post.Hoc.Pooling[aci$ID == "07DEn.0040"]<-"Cannot be pooled"

    if(pooling != T) {
      # 	F = (MSPxIC / MSCC(P)) 03DE.0030
      aci$Post.Hoc.Pooling[aci$ID == "03DE.0030"]<-"Cannot be pooled"

      # 	F = (MSIC / MSPxIC) 02E.0030
      aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"

      # 	F = (MSP / MSCC(P)) 01D.0000
      aci$Post.Hoc.Pooling[aci$ID == "01D.0000"]<-"Cannot be pooled"
    }

    if(pooling == T) {
      # 	F = (MSPxIC / MSCC(P)) 03DE.0030
      if(aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0030"]<-"Cannot be pooled"
      # 	If varCC(P) <= 0, then F = (MSPxIC / MSpooled), where MSpooled = (SSCC(P) + SSRes) / (dfCC(P) + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0031"]<-"Can be pooled"
      DE03<- c("03DE.0030","03DE.0031")
      DE03n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0030"]), is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0031"])) == F)

      # 	F = (MSIC / MSPxIC) 02E.0030
      if(aci$P.Value.Upper.Tail[aci$ID == "03DE.0030"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      # 	If varPxIC <= 0 and varCC(P) > 0, then F = (MS_IC / MSpooled), where MSpooled = (SSPxIC + SSCC(P)) / (dfPxIC + dfCC(P)) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == DE03[DE03n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
      # 	If varPxIC <= 0 and varCC(P) <= 0, then F = (MSIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == DE03[DE03n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"

      # 	F = (MSP / MSCC(P)) 01D.0000
      if(aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "01D.0000"]<-"Cannot be pooled"
      # 	If varCC(P) <= 0, then F = (MSP / MSpooled), where MSpooled = (SSCC(P) + SSRes) / (dfCC(P) + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "01D.0001"]<-"Can be pooled"
    }
  }


  if(aav.model == "DEG") {
    # 	F = (MS_S(CC(P)) / MSRes)
    aci$Post.Hoc.Pooling[aci$ID == "08DEGn.0040"]<-"Cannot be pooled"

    if(pooling != T) {
      # 	F = (MS_CC(P) / MSS(CC(P)))
      aci$Post.Hoc.Pooling[aci$ID == "07DEn.0040"]<-"Cannot be pooled"

      # 	F = (MS_PxS(IC) / MSS(CC(P))), or
      aci$Post.Hoc.Pooling[aci$ID == "06DEG.0030"]<-"Cannot be pooled"

      # 	F = (MS_S(IC) / MSPxS(IC)), or
      aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]<-"Cannot be pooled"

      # If the two terms cannot be eliminated (* or ns)
      if(aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "06DEG.0030"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0030"]<-"Cannot be pooled"
      # 	If varPxS(IC) <= 0 and varCC(P) > 0, then F = (MS_PxIC / MS_CC(P))
      if(aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "06DEG.0030"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0030"]<-"Can be pooled"
      # 	If varCC(P) <= 0 and varPxS(IC) > 0, then F = (MS_PxIC / MSPxS(IC))
      if(aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "06DEG.0030"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0031"]<-"Can be pooled"
      # 	If varCC(P) <= 0, varPxS(IC) <= 0, then F = (MSPxIC / MSpooled), where MSpooled = (SSCC(P) + SSPxS(IC) + SSS(CC(P))) / (dfCC(P) + dfPxS(IC) + dfS(CC(P))) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "06DEG.0030"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0032"]<-"Can be pooled"
      DE03<- c("03DE.0030","03DE.0031","03DE.0032")
      DE03n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0030"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0031"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0032"])) == F)

      # If the two terms cannot be eliminated (* or ns)
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0030"] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03[DE03n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0030"] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03[DE03n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "06DEG.0030"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      # 	If varS(IC) <= 0 and varPxIC > 0, then F = (MS_IC / MS_PxIC)
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0030"] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03[DE03n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Can be pooled"
      # 	If varPxIC <= 0 and varS(IC) > 0, then F = (MS_IC / MS_S(IC))
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0030"] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03[DE03n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
      # 	If varPxIC <= 0, varS(IC) <= 0, varCC(P) <= 0, varPxS(IC) > 0, then F = (MSIC / MSpooled), where MSpooled = (SS_PxIC + SS_S(IC) + SS_PxS(IC)) / (dfPxIC + dfS(IC) + dfPxS(IC))
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0030"] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03[DE03n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "06DEG.0030"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"
      # 	If varPxIC <= 0, varS(IC) <= 0, varPxS(IC) <= 0 or varCC(P) > 0, then F = (MSIC / MSpooled), where MSpooled = (SSPxIC + SSS(IC) + SS_CC(P)) / (dfPxIC + dfS(IC) + dfCC(P))
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0030"] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03[DE03n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "06DEG.0030"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0033"]<-"Can be pooled"
      # 	If varPxIC <= 0, varS(IC) <= 0, varPxS(IC) <= 0 and varCC(P) <= 0, then F = (MSIC / MSpooled), where MSpooled = (SSPxIC + SSS(IC) + SSPxS(IC) + SSCC(P) + SSS(CC(P))) / (dfPxIC + dfS(IC) + dfPxS(IC) + dfCC(P) + dfS(CC(P))) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == "05EG.0030"] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03[DE03n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "06DEG.0030"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0034"]<-"Can be pooled"

      # 	F = (MSP / MS_CC(P))
      aci$Post.Hoc.Pooling[aci$ID == "01D.0000"]<-"Cannot be pooled"
    }

    if(pooling == T) {
      # 	F = (MS_CC(P) / MSS(CC(P)))
      if(aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "07DEn.0040"]<-"Cannot be pooled"
      # 	If varS(CC(P)) <= 0, then F = (MSCC(P) / MSpooled), where MSpooled = (SS_S(CC(P)) + SSRes) / (dfS(CC(P)) + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "07DEn.0041"]<-"Can be pooled"
      DEn07<- c("07DEn.0040","07DEn.0041")
      DEn07n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "07DEn.0040"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "07DEn.0041"])) == F)

      # 	F = (MS_PxS(IC) / MSS(CC(P))), or
      if(aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "06DEG.0030"]<-"Cannot be pooled"
      # 	If varS(CC(P)) <= 0, then F = (MSPxS(IC) / MSpooled), where MSpooled = (SS_S(CC(P)) + SSRes) / (dfS(CC(P)) + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "06DEG.0031"]<-"Can be pooled"
      DEG06<- c("06DEG.0030","06DEG.0031")
      DEG06n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "06DEG.0030"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "06DEG.0031"])) == F)

      # 	F = (MS_S(IC) / MSPxS(IC)), or
      if(aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]<-"Cannot be pooled"
      # 	If varPxS(IC) <= 0 and varS(CC)P)) > 0, then F = (MSS(IC) / MSpooled),
      # where MSpooled = (SSPxS(IC) + SSS(CC(P))) / (dfPxS(IC) + dfS(CC(P))) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]<-"Can be pooled"
      # 	If varPxS(IC) <= 0 and varS(CC(P)) <= 0, then F = (MSS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"]<-"Can be pooled"
      EG05<- c("05EG.0030","05EG.0031","05EG.0032")
      EG05n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"])) == F)

      # If the two terms cannot be eliminated (* or ns)
      if(aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0030"]<-"Cannot be pooled"
      # 	If varPxS(IC) <= 0 and varCC(P) > 0, then F = (MS_PxIC / MS_CC(P))
      if(aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0030"]<-"Can be pooled"
      # 	If varCC(P) <= 0 and varPxS(IC) > 0, then F = (MS_PxIC / MSPxS(IC))
      if(aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0031"]<-"Can be pooled"
      # 	If varCC(P) <= 0, varPxS(IC) <= 0 and varS(CC(P)) > 0, then F = (MSPxIC / MSpooled),
      # where MSpooled = (SSCC(P) + SSPxS(IC) + SSS(CC(P))) / (dfCC(P) + dfPxS(IC) + dfS(CC(P))) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0032"]<-"Can be pooled"
      # 	If varCC(P) <= 0, varPxS(IC) <= 0 and varS(CC(P)) <= 0, then F = (MS_PxIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "03DE.0033"]<-"Can be pooled"
      DE03p<- c("03DE.0030","03DE.0031","03DE.0032","03DE.0033")
      DE03pn<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0030"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0031"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0032"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0033"])) == F)

      # If the two terms cannot be eliminated (* or ns)
      if(aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03p[DE03pn]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      if(aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03p[DE03pn]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Cannot be pooled"
      # 	If varS(IC) <= 0 and varPxIC > 0, then F = (MS_IC / MS_PxIC)
      if(aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03p[DE03pn]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]<-"Can be pooled"
      # 	If varPxIC <= 0 and varS(IC) > 0, then F = (MS_IC / MS_S(IC))
      if(aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03p[DE03pn]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]<-"Can be pooled"
      # 	If varPxIC <= 0, varS(IC) <= 0, varCC(P) <= 0, varPxS(IC) > 0, then F = (MSIC / MSpooled), where MSpooled = (SS_PxIC + SS_S(IC) + SS_PxS(IC)) / (dfPxIC + dfS(IC) + dfPxS(IC))
      if(aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03p[DE03pn]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]<-"Can be pooled"
      # 	If varPxIC <= 0, varS(IC) <= 0, varPxS(IC) <= 0 or varCC(P) > 0, then F = (MSIC / MSpooled),
      # where MSpooled = (SSPxIC + SSS(IC) + SS_CC(P)) / (dfPxIC + dfS(IC) + dfCC(P))
      if(aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03p[DE03pn]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] < p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0033"]<-"Can be pooled"
      # 	If varPxIC <= 0, varS(IC) <= 0, varPxS(IC) <= 0, varCC(P) <= 0 and varS(CC(P)) < 0, then F = (MSIC / MSpooled),
      # where MSpooled = (SSPxIC + SSS(IC) + SSPxS(IC) + SSCC(P) + SSS(CC(P))) / (dfPxIC + dfS(IC) + dfPxS(IC) + dfCC(P) + dfS(CC(P))) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03p[DE03pn]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0034"]<-"Can be pooled"
      # 	If varPxIC <= 0, varS(IC) <= 0, varPxS(IC) <= 0, varCC(P) <= 0 and varS(CC(P)) <= 0,
      # then F = (MSIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == EG05[EG05n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DE03p[DE03pn]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == DEG06[DEG06n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "02E.0035"]<-"Can be pooled"

      # 	F = (MSP / MS_CC(P))
      if(aci$P.Value.Upper.Tail[aci$ID == "07DEn.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "01D.0000"]<-"Cannot be pooled"
      # 	If varCC(P) <= 0 and varS(CC(P)) > 0, then F = (MSP / MSpooled), where MSpooled = (SSCC(P) + SS_S(CC(P))) / (dfCC(P) + dfS(CC(P))) = SS1/df1
      if(aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] < p.pooling) aci$Post.Hoc.Pooling[aci$ID == "01D.0001"]<-"Can be pooled"
      # 	If varCC(P) <= 0 and varS(CC(P)) <= 0, then F = (MSP / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
      if(aci$P.Value.Upper.Tail[aci$ID == DEn07[DEn07n]] >= p.pooling &
         aci$P.Value.Upper.Tail[aci$ID == "08DEGn.0040"] >= p.pooling) aci$Post.Hoc.Pooling[aci$ID == "01D.0002"]<-"Can be pooled"
    }
  }

  ###------------------------------------------ Next.Step & Interpretation


  ### Sites vs Times - T(P)xS(L)
  if(aav.model == "BEGC") {
    BEGC11<- c("11BEGC.0030","11BEGC.0031")
    BEGC11n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0030"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "11BEGC.0031"])) == F)
    # T(P)xS(C) = ns; T(P)xS(I) = ns
    if(aci$Significant[aci$ID == "11BEGC.0040"] != "*") {
      if(aci$Significant[aci$ID == BEGC11[BEGC11n]] != "*") {
        aci$Next.Step[aci$ID == BEGC11[BEGC11n]]<-aci$Next.Step[aci$ID == "11BEGC.0040"]<- "Go to PxS(C) and T(P)xC"
        aci$Interpretation[aci$ID == BEGC11[BEGC11n]]<-aci$Interpretation[aci$ID == "11BEGC.0040"]<- "No short-term and small-scale interactions"
      }
      # T(P)xS(C) = ns; T(P)xS(I) = *
      if(aci$Significant[aci$ID == BEGC11[BEGC11n]] == "*") {
        aci$Next.Step[aci$ID == BEGC11[BEGC11n]]<-aci$Next.Step[aci$ID == "11BEGC.0040"]<- "Go to PxS(C) and T(P)xC"
        aci$Interpretation[aci$ID == BEGC11[BEGC11n]]<-aci$Interpretation[aci$ID == "11BEGC.0040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and small-scale: temporal trends in impact sites differs from controls"
      }
    }
    # T(P)xS(C) = *; T(P)xS(I) = ns
    if(aci$Significant[aci$ID == "11BEGC.0040"] == "*") {
      if(aci$Significant[aci$ID == BEGC11[BEGC11n]] != "*") {
        aci$Next.Step[aci$ID == BEGC11[BEGC11n]]<-aci$Next.Step[aci$ID == "11BEGC.0040"]<- "Go to PxS(C) and T(P)xC"
        aci$Interpretation[aci$ID == BEGC11[BEGC11n]]<-aci$Interpretation[aci$ID == "11BEGC.0040"]<- "Short-term and small-scale interactions: temporal trends are not associated with impact sites"
      }
      # T(P)xS(C) = *; T(P)xS(I) = *
      if(aci$Significant[aci$ID == BEGC11[BEGC11n]] == "*") {
        aci$Next.Step[aci$ID == BEGC11[BEGC11n]]<-aci$Next.Step[aci$ID == "11BEGC.0040"]<- "Go to PxS(C) and T(P)xC"
        aci$Interpretation[aci$ID == BEGC11[BEGC11n]]<-aci$Interpretation[aci$ID == "11BEGC.0040"]<- "CAUTION - to many variability to detec a impact from short-term and small-scale: temporal trends in impact sites differs from controls"
      }
    }

    ### Sites vs Periods - PxS(L)
    BEG09<- c("09BEG.0030","09BEG.0031","09BEG.0032")
    BEG09n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "09BEG.0030"]),
                    is.na(aci$Post.Hoc.Pooling[aci$ID == "09BEG.0031"]),
                    is.na(aci$Post.Hoc.Pooling[aci$ID == "09BEG.0032"])) == F)
    # PxS(C) = ns; PxS(I) = ns
    if(aci$Significant[aci$ID == "09BEG.0040"] != "*") {
      if(aci$Significant[aci$ID == BEG09[BEG09n]] != "*")
        aci$Next.Step[aci$ID == BEG09[BEG09n]]<-aci$Next.Step[aci$ID == "09BEG.0040"]<- "Go to S(C) and PxC"
      aci$Interpretation[aci$ID == BEG09[BEG09n]]<-aci$Interpretation[aci$ID == "09BEG.0040"]<-"No Medium-term and small-scale interactions"
      # PxS(C) = ns; PxS(I) = *
      if(aci$Significant[aci$ID == BEG09[BEG09n]] == "*") {
        aci$Next.Step[aci$ID == BEG09[BEG09n]]<-aci$Next.Step[aci$ID == "09BEG.0040"]<- "Go to S(C) and PxC"
        aci$Interpretation[aci$ID == BEG09[BEG09n]]<-aci$Interpretation[aci$ID == "09BEG.0040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from medium-term and small-scale: temporal trends in impact sites differs from controls"
      }
    }
    # PxS(C) = *; PxS(I) = ns
    if(aci$Significant[aci$ID == "09BEG.0040"] == "*") {
      if(aci$Significant[aci$ID == BEG09[BEG09n]] != "*")
        aci$Next.Step[aci$ID == BEG09[BEG09n]]<-aci$Next.Step[aci$ID == "09BEG.0040"]<- "Go to S(C) and PxC"
      aci$Interpretation[aci$ID == BEG09[BEG09n]]<-aci$Interpretation[aci$ID == "09BEG.0040"]<-"Medium-term and small-scale interactions: temporal trends are not associated with impact sites"
      # PxS(C) = *; PxS(I) = *
      if(aci$Significant[aci$ID == BEG09[BEG09n]] == "*") {
        aci$Next.Step[aci$ID == BEG09[BEG09n]]<-aci$Next.Step[aci$ID == "09BEG.0040"]<- "Go to S(C) and PxC"
        aci$Interpretation[aci$ID == BEG09[BEG09n]]<-aci$Interpretation[aci$ID == "09BEG.0040"]<- "CAUTION - to many variability to detec a impact from medium-term and small-scale: temporal trends in impact sites differs from controls"
      }
    }
  }

  ### Sites - S(L)
  if(aav.model == "BEGC" | aav.model == "CEG") {
    EG05<- c("05EG.0030","05EG.0031","05EG.0032","05EG.0033")
    EG05n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0033"])) == F)

    # S(C) = ns; S(I) = ns
    if(aci$Significant[aci$ID == "05EG.0040"] != "*") {
      if(aci$Significant[aci$ID == EG05[EG05n]] != "*")
        aci$Next.Step[aci$ID == EG05[EG05n]]<-aci$Next.Step[aci$ID == "05EG.0040"]<- "Go to C"
      aci$Interpretation[aci$ID == EG05[EG05n]]<-aci$Interpretation[aci$ID == "05EG.0040"]<- "No long-term and small-scale differences among control and impact sites"
      # S(C) = ns; S(I) = *
      if(aci$Significant[aci$ID == EG05[EG05n]] == "*") {
        aci$Next.Step[aci$ID == EG05[EG05n]]<-aci$Next.Step[aci$ID == "05EG.0040"]<- "Go to C"
        aci$Interpretation[aci$ID == EG05[EG05n]]<-aci$Interpretation[aci$ID == "05EG.0040"]<- "PRESS DISTURBANCE - IMPACT DETECTED from long-term and small-scale: impact sites differs from controls"
      }
    }
    # S(C) = *; S(I) = ns
    if(aci$Significant[aci$ID == "05EG.0040"] == "*") {
      if(aci$Significant[aci$ID == EG05[EG05n]] != "*")
        aci$Next.Step[aci$ID == EG05[EG05n]]<-aci$Next.Step[aci$ID == "05EG.0040"]<- "Go to C"
      aci$Interpretation[aci$ID == EG05[EG05n]]<-aci$Interpretation[aci$ID == "05EG.0040"]<- "Long-term and small-scale differences: spatial trends are not associated with impact sites"

      # S(C) = *; S(I) = *
      if(aci$Significant[aci$ID == EG05[EG05n]] == "*") {
        aci$Next.Step[aci$ID == EG05[EG05n]]<-aci$Next.Step[aci$ID == "05EG.0040"]<- "Go to C"
        aci$Interpretation[aci$ID == EG05[EG05n]]<-aci$Interpretation[aci$ID == "05EG.0040"]<- "CAUTION - to many variability to detec a impact from long-term and small-scale: impact sites differs from controls"
      }
    }
  }


  ### Location vs Period in Time - T(P) x L
  if(aav.model == "BEGC" | aav.model == "BEC") {
    BEC10<- c("10BEC.0030","10BEC.0031","10BEC.0032")
    BEC10n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0030"]),
                    is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0031"]),
                    is.na(aci$Post.Hoc.Pooling[aci$ID == "10BEC.0032"])) == F)

    # T(P)xC = ns; T(P)xI = ns
    if(aci$Significant[aci$ID == "10BEC.0040"] != "*") {
      if(aci$Significant[aci$ID == BEC10[BEC10n]] != "*") {
        aci$Next.Step[aci$ID == BEC10[BEC10n]]<-aci$Next.Step[aci$ID == "10BEC.0040"]<- "Go to PxC"
        aci$Interpretation[aci$ID == BEC10[BEC10n]]<-aci$Interpretation[aci$ID == "10BEC.0040"]<- "No short-term and local-scale interactions"
      }
      # T(P)xC = ns; T(P)xI = *
      if(aci$Significant[aci$ID == BEC10[BEC10n]] == "*") {
        aci$Next.Step[aci$ID == BEC10[BEC10n]]<-aci$Next.Step[aci$ID == "10BEC.0040"]<- "Go to PxC"
        aci$Interpretation[aci$ID == BEC10[BEC10n]]<-aci$Interpretation[aci$ID == "10BEC.0040"]<- "PRESS DISTURBANCE - IMPACT DETECTED from short-term and local-scale: temporal trends in impact location differs from controls"
      }
    }
    # T(P)xC = *; T(P)xI = ns
    if(aci$Significant[aci$ID == "10BEC.0040"] == "*") {
      if(aci$Significant[aci$ID == BEC10[BEC10n]] != "*") {
        aci$Next.Step[aci$ID == BEC10[BEC10n]]<-aci$Next.Step[aci$ID == "10BEC.0040"]<- "Go to PxC"
        aci$Interpretation[aci$ID == BEC10[BEC10n]]<-aci$Interpretation[aci$ID == "10BEC.0040"]<- "Short-term and local-scale interactions: temporal trends are not associated with impact location"
      }
      # T(P)xC = *; T(P)xI = *
      if(aci$Significant[aci$ID == BEC10[BEC10n]] == "*") {
        aci$Next.Step[aci$ID == BEC10[BEC10n]]<-aci$Next.Step[aci$ID == "10BEC.0040"]<- "Go to PxC"
        aci$Interpretation[aci$ID == BEC10[BEC10n]]<-aci$Interpretation[aci$ID == "10BEC.0040"]<- "CAUTION - to many variability to detec a impact: temporal trends in impact location differs from controls"
      }
    }

    # Location vs Periods - P x L
    BE03.40<- c("03BE.0040","03BE.0041","03BE.0042")
    BE03.40n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0040"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0041"]),
                      is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0042"])) == F)

    BE03<- c("03BE.0030","03BE.0031","03BE.0032","03BE.0033","03BE.0034")
    BE03n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0030"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0031"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0032"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0033"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "03BE.0034"])) == F)

    # PxC = ns; PxI = ns
    if(aci$Significant[aci$ID == BE03.40[BE03.40n]] != "*") {
      if(aci$Significant[aci$ID == BE03[BE03n]] != "*") {
        aci$Next.Step[aci$ID == BE03[BE03n]]<-aci$Next.Step[aci$ID == BE03.40[BE03.40n]]<- "Go to C"
        aci$Interpretation[aci$ID == BE03[BE03n]]<-aci$Interpretation[aci$ID == BE03.40[BE03.40n]]<- "No medium-term and local-scale interactions"
      }
      # PxC = ns; PxI = *
      if(aci$Significant[aci$ID == BE03[BE03n]] == "*") {
        aci$Next.Step[aci$ID == BE03[BE03n]]<-aci$Next.Step[aci$ID == BE03.40[BE03.40n]]<- "Go to C"
        aci$Interpretation[aci$ID == BE03[BE03n]]<-aci$Interpretation[aci$ID == BE03.40[BE03.40n]]<- "PULSE DISTURBANCE - IMPACT DETECTED from medium-term and local-scale: temporal trends in impact location differs from controls"
      }
    }
    # PxC = *;PxI = ns
    if(aci$Significant[aci$ID == BE03.40[BE03.40n]] == "*") {
      if(aci$Significant[aci$ID == BE03[BE03n]] != "*") {
        aci$Next.Step[aci$ID == BE03[BE03n]]<- aci$Next.Step[aci$ID == BE03.40[BE03.40n]]<- "Go to C"
        aci$Interpretation[aci$ID == BE03[BE03n]]<-aci$Interpretation[aci$ID == BE03.40[BE03.40n]]<- "Medium-term and local-scale interactions: temporal trends are not associated with location impact"
      }
      # PxC = *; PxI = *
      if(aci$Significant[aci$ID == BE03[BE03n]] == "*") {
        aci$Next.Step[aci$ID == BE03[BE03n]]<-aci$Next.Step[aci$ID == BE03.40[BE03.40n]]<- "Go to C"
        aci$Interpretation[aci$ID == BE03[BE03n]]<-aci$Interpretation[aci$ID == BE03.40[BE03.40n]]<- "CAUTION - to many variability to detec a impact: temporal trends in impact location differs from controls"
      }
    }
  }


  # Location - L
  if(aav.model == "BEGC" | aav.model == "BEC" | aav.model == "CEG" | aav.model == "CE") {
    E02.40<- c("02E.0040","02E.0041","02E.0042","02E.0043","02E.0044")
    E02.40n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0040"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0041"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0042"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0043"]),
                     is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0044"])) == F)

    E02<- c("02E.0030","02E.0031","02E.0032","02E.0033","02E.0034", "02E.0035", "02E.0036")
    E02n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]),
                  is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]),
                  is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]),
                  is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0033"]),
                  is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0034"]),
                  is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0035"]),
                  is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0036"])) == F)
    # C = ns; I = ns
    if(aci$Significant[aci$ID == E02.40[E02.40n]] != "*") {
      if(aci$Significant[aci$ID == E02[E02n]] != "*") {
        aci$Next.Step[aci$ID == E02[E02n]]<-aci$Next.Step[aci$ID == E02.40[E02.40n]]<- "END"
        aci$Interpretation[aci$ID == E02[E02n]]<-aci$Interpretation[aci$ID == E02.40[E02.40n]]<- "No long-term and local-scale differences among control and impact locations"
      }
      # C = ns; I = *
      if(aci$Significant[aci$ID == E02[E02n]] == "*") {
        aci$Next.Step[aci$ID == E02[E02n]]<-aci$Next.Step[aci$ID == E02.40[E02.40n]]<- "END"
        aci$Interpretation[aci$ID == E02[E02n]]<-aci$Interpretation[aci$ID == E02.40[E02.40n]]<- "PRESS DISTURBANCE - IMPACT DETECTED from long-term and local-scale: impact location differs from controls"
      }
    }
    # C = *; I = ns
    if(aci$Significant[aci$ID == E02.40[E02.40n]] == "*") {
      if(aci$Significant[aci$ID == E02[E02n]] != "*") {
        aci$Next.Step[aci$ID == E02[E02n]]<-aci$Next.Step[aci$ID == E02.40[E02.40n]]<- "END"
        aci$Interpretation[aci$ID == E02[E02n]]<-aci$Interpretation[aci$ID == E02.40[E02.40n]]<- "Long-term and local-scale differences: spatial trends are not associated with impact location"
      }
      # C = *; I = *
      if(aci$Significant[aci$ID == E02[E02n]] == "*") {
        aci$Next.Step[aci$ID == E02[E02n]]<-aci$Next.Step[aci$ID == E02.40[E02.40n]]<- "END"
        aci$Interpretation[aci$ID == E02[E02n]]<-aci$Interpretation[aci$ID == E02.40[E02.40n]]<- "CAUTION - to many variability to detec a impact: impact location differs from controls"
      }
    }
  }


  ### Sites vs Times (no Periods) - T x S(L)
  if(aav.model == "CEG") {
    CEG11<- c("11CEG.0030","11CEG.0031")
    CEG11n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "11CEG.0030"]), is.na(aci$Post.Hoc.Pooling[aci$ID == "11CEG.0031"])) == F)

    # TxS(C) = ns; TxS(I) = ns
    if(aci$Significant[aci$ID == "11CEG.0040"] != "*") {
      if(aci$Significant[aci$ID == CEG11[CEG11n]] != "*") {
        aci$Next.Step[aci$ID == CEG11[CEG11n]]<-aci$Next.Step[aci$ID == "11CEG.0040"]<- "Go to S(C) and TxC"
        aci$Interpretation[aci$ID == CEG11[CEG11n]]<-aci$Interpretation[aci$ID == "11CEG.0040"]<- "No short-term and small-scale interactions"
      }
      # TxS(C) = ns; TxS(I) = *
      if(aci$Significant[aci$ID == CEG11[CEG11n]] == "*") {
        aci$Next.Step[aci$ID == CEG11[CEG11n]]<-aci$Next.Step[aci$ID == "11CEG.0040"]<- "Go to S(C) and TxC"
        aci$Interpretation[aci$ID == CEG11[CEG11n]]<-aci$Interpretation[aci$ID == "11CEG.0040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and small-scale: temporal trends in impact sites differs from controls"
      }
    }
    # TxS(C) = *; TxS(I) = ns
    if(aci$Significant[aci$ID == "11CEG.0040"] == "*") {
      if(aci$Significant[aci$ID == CEG11[CEG11n]] != "*") {
        aci$Next.Step[aci$ID == CEG11[CEG11n]]<-aci$Next.Step[aci$ID == "11CEG.0040"]<- "Go to S(C) and TxC"
        aci$Interpretation[aci$ID == CEG11[CEG11n]]<-aci$Interpretation[aci$ID == "11CEG.0040"]<- "Short-term and small-scale interactions: temporal trends are not associated with impact sites"
      }
      # TxS(C) = *; TxS(I) = *
      if(aci$Significant[aci$ID == CEG11[CEG11n]] == "*") {
        aci$Next.Step[aci$ID == CEG11[CEG11n]]<-aci$Next.Step[aci$ID == "11CEG.0040"]<- "Go to S(C) and TxC"
        aci$Interpretation[aci$ID == CEG11[CEG11n]]<-aci$Interpretation[aci$ID == "11CEG.0040"]<- "CAUTION - to many variability to detec a impact: temporal trends in impact sites differs from controls"
      }
    }
  }

  # Location vs Times (no Periods) - T x L
  if(aav.model == "CEG" | aav.model == "CE") {
    CE10<- c("10CE.0030","10CE.0031","10CE.0032")
    CE10n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0030"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0031"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "10CE.0032"])) == F)

    # TxC = ns; TxI = ns
    if(aci$Significant[aci$ID == "10CE.0040"] != "*") {
      if(aci$Significant[aci$ID == CE10[CE10n]] != "*") {
        aci$Next.Step[aci$ID == CE10[CE10n]]<-aci$Next.Step[aci$ID == "10CE.0040"]<- "Go to C"
        aci$Interpretation[aci$ID == CE10[CE10n]]<-aci$Interpretation[aci$ID == "10CE.0040"]<- "No short-term and local-scale interactions"
      }
      # TxC = ns; TxI = *
      if(aci$Significant[aci$ID == CE10[CE10n]] == "*") {
        aci$Next.Step[aci$ID == CE10[CE10n]]<-aci$Next.Step[aci$ID == "10CE.0040"]<- "Go to C"
        aci$Interpretation[aci$ID == CE10[CE10n]]<-aci$Interpretation[aci$ID == "10CE.0040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and local-scale: temporal trends in impact location differs from controls"
      }
    }
    # TxC = *; TxI = ns
    if(aci$Significant[aci$ID == "10CE.0040"] == "*") {
      if(aci$Significant[aci$ID == CE10[CE10n]] != "*") {
        aci$Next.Step[aci$ID == CE10[CE10n]]<-aci$Next.Step[aci$ID == "10CE.0040"]<- "Go to C"
        aci$Interpretation[aci$ID == CE10[CE10n]]<-aci$Interpretation[aci$ID == "10CE.0040"]<- "short-term and local-scale interactions: temporal trends are not associate with impact location"
      }
      # TxC) = *; TxI = *
      if(aci$Significant[aci$ID == CE10[CE10n]] == "*") {
        aci$Next.Step[aci$ID == CE10[CE10n]]<-aci$Next.Step[aci$ID == "10CE.0040"]<- "Go to C"
        aci$Interpretation[aci$ID == CE10[CE10n]]<-aci$Interpretation[aci$ID == "10CE.0040"]<- "CAUTION - to many variability to detec a impact: temporal trends in impact location differs from controls"
      }
    }
  }

  #  Sites(in Locations(in Places))
  if(aav.model == "DEG") {
    # S(C(P))
    aci$Next.Step[aci$ID == "08DEGn.0040"]<-"Go to PxS(Tr) and C(P)"
    ifelse(aci$Significant[aci$ID == "08DEGn.0040"] != "*",
           aci$Interpretation[aci$ID == "08DEGn.0040"]<- "Control sites do not differ",
           aci$Interpretation[aci$ID == "08DEGn.0040"]<-aci$Interpretation[aci$ID == "08DEGn.0040"]<- "Control sites differ")

    # PxS(I)
    DEG06<- c("06DEG.0030","06DEG.0031")
    DEG06n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "06DEG.0030"]),
                    is.na(aci$Post.Hoc.Pooling[aci$ID == "06DEG.0031"])) == F)
    aci$Next.Step[aci$ID == DEG06[DEG06n]]<-"Go to S(Tr) and PxTr"
    ifelse(aci$Significant[aci$ID == DEG06[DEG06n]] != "*",
           aci$Interpretation[aci$ID == DEG06[DEG06n]]<- "Differences between impact and control sites do not vary from place to place",
           aci$Interpretation[aci$ID == DEG06[DEG06n]]<-aci$Interpretation[aci$ID == DEG06[DEG06n]]<- "CAUTION - Differences between impact and control sites vary from place to place")

    # S(I)
    EG05<- c("05EG.0030","05EG.0031","05EG.0032")
    EG05n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0030"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0031"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "05EG.0032"])) == F)
    aci$Next.Step[aci$ID == EG05[EG05n]]<-"Go to Tr"
    ifelse(aci$Significant[aci$ID == EG05[EG05n]] != "*",
           aci$Interpretation[aci$ID == EG05[EG05n]]<- "Impact and control sites do not differ",
           aci$Interpretation[aci$ID == EG05[EG05n]]<-aci$Interpretation[aci$ID == EG05[EG05n]]<- "IMPACT DETECTED from small-scale: Impact and control sites differ")
  }


  #  Locations(in Places)
  if(aav.model == "DEG" | aav.model == "DE") {

    # C(P)
    DEn07<- c("07DEn.0040","07DEn.0041")
    DEn07n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "07DEn.0040"]),
                    is.na(aci$Post.Hoc.Pooling[aci$ID == "07DEn.0041"])) == F)
    aci$Next.Step[aci$ID == DEn07[DEn07n]]<-"Go to P and PxTr"
    ifelse(aci$Significant[aci$ID == DEn07[DEn07n]] != "*",
           aci$Interpretation[aci$ID == DEn07[DEn07n]]<- "Control locations do not differ",
           aci$Interpretation[aci$ID == DEn07[DEn07n]]<-aci$Interpretation[aci$ID == DEn07[DEn07n]]<- "Control locations differ")

    # PxI
    DE03p<- c("03DE.0030","03DE.0031","03DE.0032","03DE.0033")
    DE03pn<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0030"]),
                    is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0031"]),
                    is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0032"]),
                    is.na(aci$Post.Hoc.Pooling[aci$ID == "03DE.0033"])) == F)
    aci$Next.Step[aci$ID == DE03p[DE03pn]]<-"Go to Tr"
    ifelse(aci$Significant[aci$ID == DE03p[DE03pn]] != "*",
           aci$Interpretation[aci$ID == DE03p[DE03pn]]<- "Differences between impact and control locations do not vary from place to place",
           aci$Interpretation[aci$ID == DE03p[DE03pn]]<-aci$Interpretation[aci$ID == DE03p[DE03pn]]<- "CAUTION - Differences between impact and control locations vary from place to place")

    # I = Treatment
    DE02<- c("02E.0030","02E.0031","02E.0032","02E.0033","02E.0034","02E.0035")
    DE02n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0030"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0031"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0032"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0033"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0034"]),
                   is.na(aci$Post.Hoc.Pooling[aci$ID == "02E.0035"])) == F)
    aci$Next.Step[aci$ID == DE02[DE02n]]<-"END"
    ifelse(aci$Significant[aci$ID == DE02[DE02n]] != "*",
           aci$Interpretation[aci$ID == DE02[DE02n]]<- "Impact and control locations do not differ",
           aci$Interpretation[aci$ID == DE02[DE02n]]<-aci$Interpretation[aci$ID == DE02[DE02n]]<- "IMPACT DETECTED from local-scale: Impact and control locations differ")

    #P = Place
    D01<- c("01D.0000","01D.0001","01D.0002")
    D01n<-which(c(is.na(aci$Post.Hoc.Pooling[aci$ID == "01D.0000"]),
                  is.na(aci$Post.Hoc.Pooling[aci$ID == "01D.0001"]),
                  is.na(aci$Post.Hoc.Pooling[aci$ID == "01D.0002"])) == F)
    aci$Next.Step[aci$ID == D01[D01n]]<-"END"
    ifelse(aci$Significant[aci$ID == D01[D01n]] != "*",
           aci$Interpretation[aci$ID == D01[D01n]]<- "Places do not differ",
           aci$Interpretation[aci$ID == D01[D01n]]<-aci$Interpretation[aci$ID == D01[D01n]]<- "Places differ")
  }


  # For certain taxa, there were differences between marina and control locations at some of the places.
  # It is important to note that these differences are only correlative;
  #these data do not demonstrate that the marinas themselves are the cause of any differences,
  # only that some activity (either past and/or present) associated with those particular locations is responsible
  # for the differences.

  #If the same patterns had occurred at all three marinas and if the locations were sampled
  # independently at a number of times and similar patterns emerged each time, this would implicate the
  #marinas more convincingly. Nevertheless, these results support the model that marinas are affecting assemblages
  #of epibiota growing on pilings and that this effect varies among locations. There were, however, a few
  #confounding factors in these comparisons.


  ###---------------------------------------------------------------  Final.Table
  if(aav.model == "BEGC" | aav.model == "BEC" | aav.model == "CEG" | aav.model == "CE") {
    aci$Final.Table<-ifelse(is.na(aci$P.Value.Upper.Tail) == T, "Yes",
                            ifelse(!is.na(aci$Post.Hoc.Pooling) == T, "Yes", "No"))
  }

  if(aav.model == "DE" | aav.model == "DEG") {
    aci$Final.Table<-ifelse(!is.na(aci$Post.Hoc.Pooling) == T, "Yes", "No")
  }

  ###------------------------------------------------------------------------- END

  return(aci)
}

###------------------------------------------------------------------------F I M

