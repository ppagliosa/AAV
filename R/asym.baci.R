#' @title Asymmetrical ANOVA From Beyond BACI Design
#'
#' @description
#' Provides an asymmetrical ANOVA from beyond BACI (before/after-control/impact) design
#' for models with 2 or 3 temporal factors and with 1 or 2 spatial factors.
#'
#' @param data a data frame with temporal factors, spatial factors and variable(s).
#' @param n.ftemp number of temporal factors.
#' @param n.fspac number of spatial factors.
#' @param names.impact name of the asymmetrical level of the Location factor (e.g. 'impact').
#' @param names.before Name of the before level of the Before/After factor (e.g. 'before').
#' @param sym.F an object containing the results returned by \code{sym4asym} function or
#'     an symmetrical analysis of variance table run with all the data,
#'     meaning the asymmetry is deal as a factor level (e.g., Control_1,
#'     Control_2 and Impact are the factor levels).
#' @param sym.C an object containing the results returned by \code{sym4asym} function or
#'     an symmetrical analysis of variance table run omitting the asymmetrical data
#'     (e.g., Control_1 and Control_2 are the factor levels).
#' @param sym.B an object containing the results returned by \code{sym4asym} function or
#'     an symmetrical analysis of variance table run omitting the after impact data.
#' @param sym.CB an object containing the results returned by \code{sym4asym} function or
#'     an symmetrical analysis of variance table run omitting the after impact
#'     and the asymmetrical data.
#' @param p.aav level of significance for the asymmetrical ANOVA. Default is \code{0.05}.
#' @param pooling remove a term from the analysis using a post-hoc pooling procedure. Default is \code{TRUE}.
#' @param method only when \code{pooling = TRUE}.
#'     Method to be used for the pooling procedure.
#'     \code{pool} to pool a term or a set of terms with zero contribution to the model or
#'     \code{eliminate} to completely exclude a term from the model. Default is \code{pool}.
#' @param p.pooling only when \code{pooling = TRUE}.
#'     Level of significance for the pooling procedure. Default is \code{0.25}.
#'
#' @details
#' Asymmetrical analyses of variance are constructed by combining
#' the Sum of Squares values from separate (symmetrical) ANOVAs.
#' In beyond BACI designs, four fully orthogonal ANOVAs are required: one on all the data,
#' one on 'control' data, one on 'before' data, and one on 'control' and 'before' data.
#'
#' The \code{asy.baci} analysis fits mixed asymmetrical ANOVA designs with a specific
#' combination of orthogonal/nested and fixed/random factors:
#' The  beyond BACI design can be used with five factors (3 temporal and 2 spatial),
#' four factors (3 temporal and 1 spatial or 2 temporal and 2 spatial), and
#' three factors (2 temporal and 1 spatial).
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
#'
#' The data frame must follow a specific structure:
#' 1- first the temporal factors, then the spatial factors, and then the variable(s);
#' 2- each factor type must appear in a hierarchical order.
#'
#' For beyond BACI design the order of the temporal factors is Before-after, Period
#' and Time (3 factors) or Before_after and Time (2 factors).
#' For spatial factors, the order is Location and Site (2 factors) or just Location (1 factor).
#'
#' If the names of the factors do not match those listed above, \code{AAV} will change them.
#'
#' When the data frame contains more than one variable, separate univariate ANOVAs are run for each one.
#'
#' In the case of a pooling structure, both the \code{pool} and \code{eliminate} methods
#' are available for the beyond BACI design.
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
#' \code{\link[AAV]{aav}}
#'
#' @examples
#' # Example for BACI design containing 3 temporal factors and 2 spatial factor:
#' data(baci.3t2s) # 5-factors and one variable dataset
#' data(symANOVA_baci.3t2s) # symmetrical ANOVAS from baci.3t2s data
#' ex.asym.baci.3t2s<-asym.baci(data = baci.3t2s, n.ftemp = 3, n.fspac = 2,
#'     names.impact = "impact",  names.before = "before",
#'     sym.F = symANOVA_baci.3t2s$full_design,
#'     sym.C = symANOVA_baci.3t2s$control_design,
#'     sym.B = symANOVA_baci.3t2s$before_design,
#'     sym.CB = symANOVA_baci.3t2s$control_before_design,
#'     p.aav = 0.05, pooling = TRUE, method = "pool", p.pooling = 0.25)
#' ex.asym.baci.3t2s
#'
#' @import dplyr
#'
#' @include sym4asym.R
#'
#' @export

asym.baci <- function(data, n.ftemp, n.fspac, names.impact, names.before = NULL,
                      sym.F, sym.C, sym.B, sym.CB,
                      p.aav = 0.05, pooling = TRUE, method = "pool", p.pooling = 0.25) {

  # Stops
  if(n.ftemp <= 1 | n.ftemp >= 4 | n.fspac == 0 | n.fspac >= 3) stop("Select different values for 'n.ftemp' or 'n.fspac'.")

  # Factors names
  if(n.ftemp == 3 & n.fspac == 2) {
    names(data)[1:5]<-c("before_after", "period", "time", "location", "site")
    aav.model <-"AEBGC"
    aav.model.lm <-  c("A","E","AE","AB","EG","ABC","AEB","AEG","AEBC","AEBG","AEBGC")
  }
  if(n.ftemp == 3 & n.fspac == 1) {
    names(data)[1:4]<-c("before_after", "period", "time", "location")
    aav.model <-"AEBC"
    aav.model.lm <-c("A","E","AE","AB","ABC","AEB","AEBC")
  }
  if(n.ftemp == 2 & n.fspac == 2) {
    names(data)[1:4]<-c("before_after", "time", "location", "site")
    aav.model <-"AEGC"
    aav.model.lm <-c("A","E","AE","EG","AC","AEC","AEG","AEGC")
  }
  if(n.ftemp == 2 & n.fspac == 1) {
    names(data)[1:3]<-c("before_after", "time", "location")
    aav.model <-"AEC"
    aav.model.lm <-c("A","E","AE","AC","AEC")
  }

  # Additional Stops
  md <- data.frame(data[,(n.ftemp+n.fspac+1):length(data)])
  if(length(md) == 1) names(md)<-colnames(data[length(data)])
  if(0 %in% colSums(md)) stop("At least one of the variable responses has only 0 values. Check your dataset.")
  md.ca<-data.frame(data[data$location != names.impact & data$before_after != names.before,(n.ftemp+n.fspac+1):length(data)])
  if(length(md.ca) == 1) names(md.ca)<-colnames(data[length(data)])
  if(0 %in% colSums(md.ca)) stop("At least one variable response has only 0 values in the 'Controls' location after the disturbance has occurred. Check your dataset.")

  # prevenir "no visible binding for global variable"
  ID<-NULL

  ### Cria tabela vazia para colocar os dados
  # importante especificar o tipo de dados para mais tarde poder juntar as linhas
  # (em geral, as funcoes so juntam dados de mesmo tipo)
  baci<-NULL
  baci<-data.frame(matrix(data=character(0), nrow = 0, ncol=19))
  names(baci)<-c("ID","Source.of.Variation","Df.Formula",
                 "Df", "Sum.Sq", "Calculated.From", "Mean.Sq",
                 "Denominator","Df.Pooled","Sum.Sq.Pooled","Mean.Sq.Pooled",
                 "F.Value","P.Value.Upper.Tail","P.Value.Lower.Tail",
                 "Significant", "Post.Hoc.Pooling", "Next.Step", "Interpretation", "Final.Table")
  baci$Df<-as.numeric(baci$Df)
  baci$Sum.Sq<-as.numeric(baci$Sum.Sq)
  baci$Mean.Sq<-as.numeric(baci$Mean.Sq)
  baci$Df.Pooled<-as.numeric(baci$Df.Pooled)
  baci$Sum.Sq.Pooled<-as.numeric(baci$Sum.Sq.Pooled)
  baci$Mean.Sq.Pooled<-as.numeric(baci$Mean.Sq.Pooled)
  baci$F.Value<-as.numeric(baci$F.Value)
  baci$P.Value.Upper.Tail<-as.numeric(baci$P.Value.Upper.Tail)
  baci$P.Value.Lower.Tail<-as.numeric(baci$P.Value.Lower.Tail)

  # Sequencia da programacao:
  # Seguindo a ordem das anovas simetricas (full, only control,)
  # Primeiro termos sem interacao
  # Depois, termos com interacao
  ###----------------------------------------------------------- FROM sym.F
  ##---------------------------------------------- code xx.0000

  # garantir q tenha o nome correto das colunas
  colnames(sym.F)[1:3] <- c("Source.of.Variation","Df", "Sum.Sq")
  colnames(sym.C)[1:3] <- c("Source.of.Variation","Df", "Sum.Sq")
  colnames(sym.B)[1:3] <- c("Source.of.Variation","Df", "Sum.Sq")
  colnames(sym.CB)[1:3] <- c("Source.of.Variation","Df", "Sum.Sq")

  f.an<-baci[F,]

  if("A" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A",1:3])
    f.an$ID[f.an$Source.of.Variation == "A"] <- "01A.0000"
    f.an$Df.Formula[f.an$ID == "01A.0000"]<-"(b-1)"
    f.an$Calculated.From[f.an$ID == "01A.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "01A.0000"]<- "Before vs After = B"
  }
  if("E" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "E",1:3])
    f.an$ID[f.an$Source.of.Variation == "E"]<-"02E.0000"
    f.an$Df.Formula[f.an$ID == "02E.0000"]<-"(l-1)"
    f.an$Calculated.From[f.an$ID == "02E.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "02E.0000"]<-"Locations = L"
  }
  if("AE" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:E",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:E"]<-"03AE.0000"
    f.an$Df.Formula[f.an$ID == "03AE.0000"]<-"(b-1)(l-1)"
    f.an$Calculated.From[f.an$ID == "03AE.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "03AE.0000"]<-"B x L"
  }
  if("AB" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:B",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:B"]<-"04AB.0000"
    f.an$Df.Formula[f.an$ID == "04AB.0000"]<-"b(p-1)"
    f.an$Calculated.From[f.an$ID == "04AB.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "04AB.0000"]<-"Period(B) = P(B)"
  }
  if("EG" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "E:G",1:3])
    f.an$ID[f.an$Source.of.Variation == "E:G"]<-"05EG.0000"
    f.an$Df.Formula[f.an$ID == "05EG.0000"]<-"l(s-1)"
    f.an$Calculated.From[f.an$ID == "05EG.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "05EG.0000"]<-"Site(L) = S(L)"
  }
  if("ABC" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:B:C",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:B:C"]<-"06ABC.0000"
    f.an$Df.Formula[f.an$ID == "06ABC.0000"]<-"bp(t-1)"
    f.an$Calculated.From[f.an$ID == "06ABC.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "06ABC.0000"]<-"Times(P(B)) = T(P(B))"
  }
  if("AEG" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:E:G",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:E:G"]<-"07AEG.0000"
    f.an$Df.Formula[f.an$ID == "07AEG.0000"]<-"(b-1)l(s-1)"
    f.an$Calculated.From[f.an$ID == "07AEG.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "07AEG.0000"]<-"B x S(L)"
  }
  if("AEB" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:E:B",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:E:B"]<-"08AEB.0000"
    f.an$Df.Formula[f.an$ID == "08AEB.0000"]<-"b(p-1)(l-1)"
    f.an$Calculated.From[f.an$ID == "08AEB.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "08AEB.0000"]<-"P(B) x L"
  }
  if("AEBG" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:E:B:G",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:E:B:G"]<-"09AEBG.0000"
    f.an$Df.Formula[f.an$ID == "09AEBG.0000"]<-"b(p-1)l(s-1)"
    f.an$Calculated.From[f.an$ID == "09AEBG.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "09AEBG.0000"]<-"P(B) x S(L)"
  }
  if("AEBC" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:E:B:C",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:E:B:C"]<-"10AEBC.0000"
    f.an$Df.Formula[f.an$ID == "10AEBC.0000"]<-"bp(t-1)(l-1)"
    f.an$Calculated.From[f.an$ID == "10AEBC.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "10AEBC.0000"]<-"T(P(B)) x L"
  }
  if("AEBGC" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:E:B:G:C",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:E:B:G:C"]<-"11AEBGC.0000"
    f.an$Df.Formula[f.an$ID == "11AEBGC.0000"]<-"bp(t-1)l(s-1)"
    f.an$Calculated.From[f.an$ID == "11AEBGC.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "11AEBGC.0000"]<-"T(P(B)) x S(L)"
  }
  if("AC" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:C",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:C"]<-"06AC.0000"
    f.an$Df.Formula[f.an$ID == "06AC.0000"]<-"b(t-1)"
    f.an$Calculated.From[f.an$ID == "06AC.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "06AC.0000"]<-"Times(B) = T(B)"
  }
  if("AEC" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:E:C",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:E:C"]<-"10AEC.0000"
    f.an$Df.Formula[f.an$ID == "10AEC.0000"]<-"b(t-1)(l-1)"
    f.an$Calculated.From[f.an$ID == "10AEC.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "10AEC.0000"]<-"T(B) x L"
  }
  if("AEGC" %in% aav.model.lm == T) {
    f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "A:E:G:C",1:3])
    f.an$ID[f.an$Source.of.Variation == "A:E:G:C"]<-"11AEGC.0000"
    f.an$Df.Formula[f.an$ID == "11AEGC.0000"]<-"b(t-1)l(s-1)"
    f.an$Calculated.From[f.an$ID == "11AEGC.0000"]<-"sym.F"
    f.an$Source.of.Variation[f.an$ID == "11AEGC.0000"]<-"T(B) x S(L)"
  }
  # RESIDUALS
  f.an<-dplyr::bind_rows(f.an,sym.F[sym.F$Source.of.Variation == "Residuals",1:3])
  f.an$ID[f.an$Source.of.Variation == "Residuals"]<-"12R.0000"
  f.an$Calculated.From[f.an$ID == "12R.0000"]<-"sym.F"
  f.an$Mean.Sq[f.an$ID == "12R.0000"]<- f.an$Sum.Sq[f.an$ID == "12R.0000"] / f.an$Df[f.an$ID == "12R.0000"]
  if(aav.model == "AEBGC") f.an$Df.Formula[f.an$ID == "12R.0000"]<-"bptls(n-1)"
  if(aav.model == "AEBC") f.an$Df.Formula[f.an$ID == "12R.0000"]<-"bptl(n-1)"
  if(aav.model == "AEGC") f.an$Df.Formula[f.an$ID == "12R.0000"]<-"btls(n-1)"
  if(aav.model == "AEC") f.an$Df.Formula[f.an$ID == "12R.0000"]<-"btl(n-1)"


  ##-------------------------------------------- AMONG Before = From sym.B
  ##----------------------------- code xx.1000
  b.an<-baci[F,]
  ###NO IF CHAMA A COMBINACAO LM ORIGINAL, MAS NA ANOVA ESTA SEM O "A"

  if("AEB" %in% aav.model.lm == T) {
    b.an<-dplyr::bind_rows(b.an,sym.B[sym.B$Source.of.Variation == "B:E",1:3])
    b.an$ID[b.an$Source.of.Variation == "B:E"]<-"08AEB.1000"
    b.an$Df.Formula[b.an$Source.of.Variation == "B:E"]<-"(p-1)(l-1)"
    b.an$Calculated.From[b.an$Source.of.Variation == "B:E"]<-"sym.B"
    b.an$Source.of.Variation[b.an$Source.of.Variation == "B:E"]<-"P(Bef) x L"
  }
  if("AEBG" %in% aav.model.lm == T) {
    b.an<-dplyr::bind_rows(b.an,sym.B[sym.B$Source.of.Variation == "B:E:G",1:3])
    b.an$ID[b.an$Source.of.Variation == "B:E:G"]<-"09AEBG.1000"
    b.an$Df.Formula[b.an$Source.of.Variation == "B:E:G"]<-"(p-1)l(s-1)"
    b.an$Calculated.From[b.an$Source.of.Variation == "B:E:G"]<-"sym.B"
    b.an$Source.of.Variation[b.an$Source.of.Variation == "B:E:G"]<-"P(Bef) x S(L)"
  }
  if("AEBC" %in% aav.model.lm == T) {
    b.an<-dplyr::bind_rows(b.an,sym.B[sym.B$Source.of.Variation == "B:E:C",1:3])
    b.an$ID[b.an$Source.of.Variation == "B:E:C"]<-"10AEBC.1000"
    b.an$Df.Formula[b.an$Source.of.Variation == "B:E:C"]<-"p(t-1)(l-1)"
    b.an$Calculated.From[b.an$Source.of.Variation == "B:E:C"]<-"sym.B"
    b.an$Source.of.Variation[b.an$Source.of.Variation == "B:E:C"]<-"T(P(Bef)) x L"
  }
  if("AEBGC" %in% aav.model.lm == T) {
    b.an<-dplyr::bind_rows(b.an,sym.B[sym.B$Source.of.Variation == "B:E:G:C",1:3])
    b.an$ID[b.an$Source.of.Variation == "B:E:G:C"]<-"11AEBGC.1000"
    b.an$Df.Formula[b.an$Source.of.Variation == "B:E:G:C"]<-"p(t-1)l(s-1)"
    b.an$Calculated.From[b.an$Source.of.Variation == "B:E:G:C"]<-"sym.B"
    b.an$Source.of.Variation[b.an$Source.of.Variation == "B:E:G:C"]<-"T(P(Bef)) x S(L)"
  }
  if("AEC" %in% aav.model.lm == T) {
    b.an<-dplyr::bind_rows(b.an,sym.B[sym.B$Source.of.Variation == "C:E",1:3])
    b.an$ID[b.an$Source.of.Variation == "C:E"]<-"10AEC.1000"
    b.an$Df.Formula[b.an$ID == "10AEC.1000"]<-"(t-1)(l-1)"
    b.an$Calculated.From[b.an$ID == "10AEC.1000"]<-"sym.B"
    b.an$Source.of.Variation[b.an$ID == "10AEC.1000"]<-"T(Bef) x L"
  }
  if("AEGC" %in% aav.model.lm == T) {
    b.an<-dplyr::bind_rows(b.an,sym.B[sym.B$Source.of.Variation == "C:E:G",1:3])
    b.an$ID[b.an$Source.of.Variation == "C:E:G"]<-"11AEGC.1000"
    b.an$Df.Formula[b.an$ID == "11AEGC.1000"]<-"(t-1)l(s-1)"
    b.an$Calculated.From[b.an$ID == "11AEGC.1000"]<-"sym.B"
    b.an$Source.of.Variation[b.an$ID == "11AEGC.1000"]<-"T(Bef) x S(L)"
  }


  ###------------------------------- ONLY AFTER = From sym.F - sym.B
  ##----------------------------- codigo xx.2000

  a.an<-baci[F,]
  ###NO IF CHAMA A COMBINACAO LM ORIGINAL, MAS NA ANOVA ESTA SEM O "A"

  if("AEB" %in% aav.model.lm == T) {
    dia.an<-cbind(Source.of.Variation = "B:E",
                  sym.F[sym.F$Source.of.Variation == "A:E:B",2:3] - sym.B[sym.B$Source.of.Variation == "B:E",2:3])
    a.an<-dplyr::bind_rows(a.an,dia.an)
    a.an$ID[a.an$Source.of.Variation == "B:E"]<-"08AEB.2000"
    a.an$Df.Formula[a.an$Source.of.Variation == "B:E"]<-"(p-1)(l-1)"
    a.an$Calculated.From[a.an$Source.of.Variation == "B:E"]<-"sym.F-sym.B"
    a.an$Source.of.Variation[a.an$Source.of.Variation == "B:E"]<-"P(Aft) x L"
  }
  if("AEBG" %in% aav.model.lm == T) {
    dia.an<-cbind(Source.of.Variation = "B:E:G",
                  sym.F[sym.F$Source.of.Variation == "A:E:B:G",2:3] - sym.B[sym.B$Source.of.Variation == "B:E:G",2:3])
    a.an<-dplyr::bind_rows(a.an,dia.an)
    a.an$ID[a.an$Source.of.Variation == "B:E:G"]<-"09AEBG.2000"
    a.an$Df.Formula[a.an$Source.of.Variation == "B:E:G"]<-"(p-1)l(s-1)"
    a.an$Calculated.From[a.an$Source.of.Variation == "B:E:G"]<-"sym.F-sym.B"
    a.an$Source.of.Variation[a.an$Source.of.Variation == "B:E:G"]<-"P(Aft) x S(L)"
  }
  if("AEBC" %in% aav.model.lm == T) {
    dia.an<-cbind(Source.of.Variation = "B:E:C",
                  sym.F[sym.F$Source.of.Variation == "A:E:B:C",2:3] - sym.B[sym.B$Source.of.Variation == "B:E:C",2:3])
    a.an<-dplyr::bind_rows(a.an,dia.an)
    a.an$ID[a.an$Source.of.Variation == "B:E:C"]<-"10AEBC.2000"
    a.an$Df.Formula[a.an$Source.of.Variation == "B:E:C"]<-"p(t-1)(l-1)"
    a.an$Calculated.From[a.an$Source.of.Variation == "B:E:C"]<-"sym.F-sym.B"
    a.an$Source.of.Variation[a.an$Source.of.Variation == "B:E:C"]<-"T(P(Aft)) x L"
  }
  if("AEBGC" %in% aav.model.lm == T) {
    dia.an<-cbind(Source.of.Variation = "B:E:G:C",
                  sym.F[sym.F$Source.of.Variation == "A:E:B:G:C",2:3] - sym.B[sym.B$Source.of.Variation == "B:E:G:C",2:3])
    a.an<-dplyr::bind_rows(a.an,dia.an)
    a.an$ID[a.an$Source.of.Variation == "B:E:G:C"]<-"11AEBGC.2000"
    a.an$Df.Formula[a.an$Source.of.Variation == "B:E:G:C"]<-"p(t-1)l(s-1)"
    a.an$Calculated.From[a.an$Source.of.Variation == "B:E:G:C"]<-"sym.F-sym.B"
    a.an$Source.of.Variation[a.an$Source.of.Variation == "B:E:G:C"]<-"T(P(Aft)) x S(L)"
  }
  if("AEC" %in% aav.model.lm == T) {
    dia.an<-cbind(Source.of.Variation = "C:E",
                  sym.F[sym.F$Source.of.Variation == "A:E:C",2:3] - sym.B[sym.B$Source.of.Variation == "C:E",2:3])
    a.an<-dplyr::bind_rows(a.an,dia.an)
    a.an$ID[a.an$Source.of.Variation == "C:E"]<-"10AEC.2000"
    a.an$Df.Formula[a.an$Source.of.Variation == "C:E"]<-"(t-1)(l-1)"
    a.an$Calculated.From[a.an$Source.of.Variation == "C:E"]<-"sym.F-sym.B"
    a.an$Source.of.Variation[a.an$Source.of.Variation == "C:E"]<-"T(Aft) x L"
  }
  if("AEGC" %in% aav.model.lm == T) {
    dia.an<-cbind(Source.of.Variation = "C:E:G",
                  sym.F[sym.F$Source.of.Variation == "A:E:G:C",2:3] - sym.B[sym.B$Source.of.Variation == "C:E:G",2:3])
    a.an<-dplyr::bind_rows(a.an,dia.an)
    a.an$ID[a.an$Source.of.Variation == "C:E:G"]<-"11AEGC.2000"
    a.an$Df.Formula[a.an$ID == "11AEGC.2000"]<-"(t-1)l(s-1)"
    a.an$Calculated.From[a.an$ID == "11AEGC.2000"]<-"sym.F-sym.B"
    a.an$Source.of.Variation[a.an$ID == "11AEGC.2000"]<-"T(Aft) x S(L)"
  }



  ###----------------------------------------- AMONG Controls = From sym.C
  ##------------------------- code xx.0040

  c.an<-baci[F,]
  if("AEGC" %in% aav.model.lm == T) {
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E:G:C",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E:G:C"]<-"11AEGC.0040"
    c.an$Df.Formula[c.an$ID == "11AEGC.0040"]<-"b(t-1)c(s-1)"
    c.an$Calculated.From[c.an$ID == "11AEGC.0040"]<-"sym.C"
    c.an$Source.of.Variation[c.an$ID == "11AEGC.0040"]<-"T(B) x S(C)"
  }
  if("AEC" %in% aav.model.lm == T) {
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E:C",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E:C"]<-"10AEC.0040"
    c.an$Df.Formula[c.an$Source.of.Variation == "A:E:C"]<-"b(t-1)(c-1)"
    c.an$Calculated.From[c.an$Source.of.Variation == "A:E:C"]<-"sym.C"
    c.an$Source.of.Variation[c.an$Source.of.Variation == "A:E:C"]<-"T(B) x C"
  }
  if("AEBGC" %in% aav.model.lm == T) {
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E:B:G:C",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E:B:G:C"]<-"11AEBGC.0040"
    c.an$Df.Formula[c.an$Source.of.Variation == "A:E:B:G:C"]<-"bp(t-1)c(s-1)"
    c.an$Calculated.From[c.an$Source.of.Variation == "A:E:B:G:C"]<-"sym.C"
    c.an$Source.of.Variation[c.an$Source.of.Variation == "A:E:B:G:C"]<-"T(P(B)) x S(C)"
  }
  if("AEBC" %in% aav.model.lm == T) {
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E:B:C",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E:B:C"]<-"10AEBC.0040"
    c.an$Df.Formula[c.an$Source.of.Variation == "A:E:B:C"]<-"bp(t-1)(c-1)"
    c.an$Calculated.From[c.an$Source.of.Variation == "A:E:B:C"]<-"sym.C"
    c.an$Source.of.Variation[c.an$Source.of.Variation == "A:E:B:C"]<-"T(P(B)) x C"
  }
  if("AEBG" %in% aav.model.lm == T) {
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E:B:G",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E:B:G"]<-"09AEBG.0040"
    c.an$Df.Formula[c.an$Source.of.Variation == "A:E:B:G"]<-"b(p-1)c(s-1)"
    c.an$Calculated.From[c.an$Source.of.Variation == "A:E:B:G"]<-"sym.C"
    c.an$Source.of.Variation[c.an$Source.of.Variation == "A:E:B:G"]<-"P(B) x S(C)"
  }
  if("AEB" %in% aav.model.lm == T) {
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E:B",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E:B"]<-"08AEB.0040"
    c.an$Df.Formula[c.an$Source.of.Variation == "A:E:B"]<-"b(p-1)(c-1)"
    c.an$Calculated.From[c.an$Source.of.Variation == "A:E:B"]<-"sym.C"
    c.an$Source.of.Variation[c.an$Source.of.Variation == "A:E:B"]<-"P(B) x C"
  }
  if("EG" %in% aav.model.lm == T) {
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E:G",1:3])
    c.an$ID[c.an$Source.of.Variation == "E:G"]<-"05EG.0040"
    c.an$Df.Formula[c.an$Source.of.Variation == "E:G"]<-"c(s-1)"
    c.an$Calculated.From[c.an$Source.of.Variation == "E:G"]<-"sym.C"
    c.an$Source.of.Variation[c.an$Source.of.Variation == "E:G"]<-"S(C)"
  }
  if("E" %in% aav.model.lm == T) {
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "E",1:3])
    c.an$ID[c.an$Source.of.Variation == "E"]<-"02E.0040"
    c.an$Df.Formula[c.an$Source.of.Variation == "E"]<-"(c-1)"
    c.an$Calculated.From[c.an$Source.of.Variation == "E"]<-"sym.C"
    c.an$Source.of.Variation[c.an$Source.of.Variation == "E"]<-"Among Controls = C"
  }



  ###----------------------------------- IMPACTs = From sym.F - sym.C
  ##------------------------- codigo xx.0030
  i.an<-baci[F,]
  if("E" %in% aav.model.lm == T) {
    dii.an<-cbind(Source.of.Variation = "E",
                  sym.F[sym.F$Source.of.Variation == "E",2:3] - sym.C[sym.C$Source.of.Variation == "E",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.an)
    i.an$ID[i.an$Source.of.Variation == "E"]<-"02E.0030"
    i.an$Df.Formula[i.an$Source.of.Variation == "E"]<-"1"
    i.an$Calculated.From[i.an$Source.of.Variation == "E"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$Source.of.Variation == "E"]<-"Impact vs Controls = I"
  }
  if("EG" %in% aav.model.lm == T) {
    dii.an<-cbind(Source.of.Variation = "E:G",
                  sym.F[sym.F$Source.of.Variation == "E:G",2:3] - sym.C[sym.C$Source.of.Variation == "E:G",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.an)
    i.an$ID[i.an$Source.of.Variation == "E:G"]<-"05EG.0030"
    i.an$Df.Formula[i.an$Source.of.Variation == "E:G"]<-"(s-1)"
    i.an$Calculated.From[i.an$Source.of.Variation == "E:G"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$Source.of.Variation == "E:G"]<-"S(I)"
  }
  if("AEB" %in% aav.model.lm == T) {
    dii.an<-cbind(Source.of.Variation = "A:E:B",
                  sym.F[sym.F$Source.of.Variation == "A:E:B",2:3] - sym.C[sym.C$Source.of.Variation == "A:E:B",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.an)
    i.an$ID[i.an$Source.of.Variation == "A:E:B"]<-"08AEB.0030"
    i.an$Df.Formula[i.an$Source.of.Variation == "A:E:B"]<-"b(p-1)"
    i.an$Calculated.From[i.an$Source.of.Variation == "A:E:B"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$Source.of.Variation == "A:E:B"]<-"P(B) x I"
  }
  if("AEBG" %in% aav.model.lm == T) {
    dii.an<-cbind(Source.of.Variation = "A:E:B:G",
                  sym.F[sym.F$Source.of.Variation == "A:E:B:G",2:3] - sym.C[sym.C$Source.of.Variation == "A:E:B:G",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.an)
    i.an$ID[i.an$Source.of.Variation == "A:E:B:G"]<-"09AEBG.0030"
    i.an$Df.Formula[i.an$Source.of.Variation == "A:E:B:G"]<-"b(p-1)(s-1)"
    i.an$Calculated.From[i.an$Source.of.Variation == "A:E:B:G"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$Source.of.Variation == "A:E:B:G"]<-"P(B) x S(I)"
  }
  if("AEBC" %in% aav.model.lm == T) {
    dii.an<-cbind(Source.of.Variation = "A:E:B:C",
                  sym.F[sym.F$Source.of.Variation == "A:E:B:C",2:3] - sym.C[sym.C$Source.of.Variation == "A:E:B:C",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.an)
    i.an$ID[i.an$Source.of.Variation == "A:E:B:C"]<-"10AEBC.0030"
    i.an$Df.Formula[i.an$Source.of.Variation == "A:E:B:C"]<-"bp(t-1)"
    i.an$Calculated.From[i.an$Source.of.Variation == "A:E:B:C"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$Source.of.Variation == "A:E:B:C"]<-"T(P(B)) x I"
  }
  if("AEBGC" %in% aav.model.lm == T) {
    dii.an<-cbind(Source.of.Variation = "A:E:B:G:C",
                  sym.F[sym.F$Source.of.Variation == "A:E:B:G:C",2:3] - sym.C[sym.C$Source.of.Variation == "A:E:B:G:C",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.an)
    i.an$ID[i.an$Source.of.Variation == "A:E:B:G:C"]<-"11AEBGC.0030"
    i.an$Df.Formula[i.an$Source.of.Variation == "A:E:B:G:C"]<-"bp(t-1)(s-1)"
    i.an$Calculated.From[i.an$Source.of.Variation == "A:E:B:G:C"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$Source.of.Variation == "A:E:B:G:C"]<-"T(P(B)) x S(I)"
  }
  if("AEC" %in% aav.model.lm == T) {
    dii.an<-cbind(Source.of.Variation = "A:E:C",
                  sym.F[sym.F$Source.of.Variation == "A:E:C",2:3] - sym.C[sym.C$Source.of.Variation == "A:E:C",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.an)
    i.an$ID[i.an$Source.of.Variation == "A:E:C"]<-"10AEC.0030"
    i.an$Df.Formula[i.an$Source.of.Variation == "A:E:C"]<-"b(t-1)"
    i.an$Calculated.From[i.an$Source.of.Variation == "A:E:C"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$Source.of.Variation == "A:E:C"]<-"T(B) x I"
  }
  if("AEGC" %in% aav.model.lm == T) {
    dii.an<-cbind(Source.of.Variation = "A:E:G:C",
                  sym.F[sym.F$Source.of.Variation == "A:E:G:C",2:3] - sym.C[sym.C$Source.of.Variation == "A:E:G:C",2:3])
    i.an<-dplyr::bind_rows(i.an,dii.an)
    i.an$ID[i.an$Source.of.Variation == "A:E:G:C"]<-"11AEGC.0030"
    i.an$Df.Formula[i.an$ID == "11AEGC.0030"]<-"b(t-1)(s-1)"
    i.an$Calculated.From[i.an$ID == "11AEGC.0030"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$ID == "11AEGC.0030"]<-"T(B) x S(I)"
  }



  ###-------------------------------------------------------------------------
  ###------------------------------------ INTERACTIONS (F-ratio por aav.model)

  ###------ Controls After  = sym.C - sym.CB
  ###------ codigo xx.2040
  ac.an<-baci[F,]
  ###----- IMPACTs (AFTER) = From sym.F - sym.B - sym.C + sym.CB
  ###------ codigo xx.2030
  ai.an<-baci[F,]


  ###------------------ INTERACTIONS with TWO temporal and ONE spatial factors

  if(aav.model == "AEC") {
    # 	F = (MS_T(Aft) xCC / MS_Res)
    di10AEC.2040<-cbind(Source.of.Variation = "C:E",
                        sym.C[sym.C$Source.of.Variation == "A:E:C",2:3] - sym.CB[sym.CB$Source.of.Variation == "C:E",2:3])
    ac.an<-dplyr::bind_rows(ac.an,di10AEC.2040)
    ac.an$ID[ac.an$Source.of.Variation == "C:E"] <- "10AEC.2040"
    ac.an$Df.Formula[ac.an$ID == "10AEC.2040"]<-"(t-1)(c-1)"
    ac.an$Calculated.From[ac.an$ID == "10AEC.2040"]<-"sym.C-sym.CB"
    ac.an$Mean.Sq[ac.an$ID == "10AEC.2040"] <- ac.an$Sum.Sq[ac.an$ID == "10AEC.2040"] / ac.an$Df[ac.an$ID == "10AEC.2040"]
    ac.an$Source.of.Variation[ac.an$ID == "10AEC.2040"] <-"T(Aft) x C"
    ac.an$Denominator[ac.an$ID == "10AEC.2040"]<-"Residuals"
    ac.an$F.Value[ac.an$ID == "10AEC.2040"]<-ac.an$Mean.Sq[ac.an$ID == "10AEC.2040"] / f.an$Mean.Sq[f.an$ID == "12R.0000"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "10AEC.2040"]<- stats::pf(ac.an$F.Value[ac.an$ID == "10AEC.2040"], ac.an$Df[ac.an$ID == "10AEC.2040"], f.an$Df[f.an$ID == "12R.0000"], lower.tail = F)

    # 	F = (MST(Aft) xIC / MST(Aft) xCC), or
    di10AEC.2030<-cbind(Source.of.Variation = "C:E",
                        sym.F[sym.F$Source.of.Variation == "A:E:C",2:3] -
                          sym.B[sym.B$Source.of.Variation == "C:E",2:3] -
                          sym.C[sym.C$Source.of.Variation == "A:E:C",2:3] +
                          sym.CB[sym.CB$Source.of.Variation == "C:E",2:3])
    ai.an<-dplyr::bind_rows(ai.an,di10AEC.2030)
    ai.an$ID[ai.an$Source.of.Variation == "C:E"] <- "10AEC.2030"
    ai.an$Df.Formula[ai.an$ID == "10AEC.2030"]<-"(t-1)"
    ai.an$Calculated.From[ai.an$ID == "10AEC.2030"]<-"sym.F-sym.B-sym.C+sym.BC"
    ai.an$Mean.Sq[ai.an$ID == "10AEC.2030"]<-ai.an$Sum.Sq[ai.an$ID == "10AEC.2030"] / ai.an$Df[ai.an$ID == "10AEC.2030"]
    ai.an$Source.of.Variation[ai.an$ID == "10AEC.2030"] <-"T(Aft) x I"
    ai.an$Denominator[ai.an$ID == "10AEC.2030"]<-"T(Aft)xC"
    ai.an$F.Value[ai.an$ID == "10AEC.2030"]<-ai.an$Mean.Sq[ai.an$ID == "10AEC.2030"] / ac.an$Mean.Sq[ac.an$ID == "10AEC.2040"]
    ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEC.2030"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEC.2030"], ai.an$Df[ai.an$ID == "10AEC.2030"], ac.an$Df[ac.an$ID == "10AEC.2040"], lower.tail = F)

    # 	F = (MSB xCC / MST(Aft) xCC)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0040"
    c.an$Calculated.From[c.an$ID == "03AE.0040"]<-"sym.C"
    c.an$Df.Formula[c.an$ID == "03AE.0040"]<-"(b-1)(c-1)"
    c.an$Mean.Sq[c.an$ID == "03AE.0040"] <- c.an$Sum.Sq[c.an$ID == "03AE.0040"] / c.an$Df[c.an$ID == "03AE.0040"]
    c.an$Source.of.Variation[c.an$ID == "03AE.0040"]<-"B x C"
    c.an$Denominator[c.an$ID == "03AE.0040"]<-"T(Aft)xC"
    c.an$F.Value[c.an$ID == "03AE.0040"]<-c.an$Mean.Sq[c.an$ID == "03AE.0040"] / ac.an$Mean.Sq[ac.an$ID == "10AEC.2040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0040"], c.an$Df[c.an$ID == "03AE.0040"], ac.an$Df[ac.an$ID == "10AEC.2040"], lower.tail = F)

    # 	F = (MSB xIC / MSB xCC), or
    di03AE.0030<-cbind(Source.of.Variation = "A:E",
                       sym.F[sym.F$Source.of.Variation == "A:E",2:3] - sym.C[sym.C$Source.of.Variation == "A:E",2:3])
    i.an<-dplyr::bind_rows(i.an,di03AE.0030)
    i.an$ID[i.an$Source.of.Variation == "A:E"]<-"03AE.0030"
    i.an$Df.Formula[i.an$ID == "03AE.0030"]<-"(b-1)"
    i.an$Calculated.From[i.an$Source.of.Variation == "A:E"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "03AE.0030"]<-i.an$Sum.Sq[i.an$ID == "03AE.0030"] / i.an$Df[i.an$ID == "03AE.0030"]
    i.an$Source.of.Variation[i.an$ID == "03AE.0030"]<-"B x I"
    i.an$Denominator[i.an$ID == "03AE.0030"]<-"BxC"
    i.an$F.Value[i.an$ID == "03AE.0030"]<-i.an$Mean.Sq[i.an$ID == "03AE.0030"] / c.an$Mean.Sq[c.an$ID == "03AE.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0030"], i.an$Df[i.an$ID == "03AE.0030"], c.an$Df[c.an$ID == "03AE.0040"], lower.tail = F)

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	If varT(Aft) xCC <= 0, then F = (MST(Aft) xIC / MSpooled), where MSpooled = (SST(Aft) xCC + SSRes) / (dfT(Aft) xCC + dfRes)
      ai.an<-dplyr::bind_rows(ai.an,di10AEC.2030)
      ai.an$ID[ai.an$Source.of.Variation == "C:E"]<-"10AEC.2031"
      ai.an$Df.Formula[ai.an$ID == "10AEC.2031"]<-"(t-1)"
      ai.an$Calculated.From[ai.an$ID == "10AEC.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
      ai.an$Mean.Sq[ai.an$ID == "10AEC.2031"]<- ai.an$Sum.Sq[ai.an$ID == "10AEC.2031"] / ai.an$Df[ai.an$ID == "10AEC.2031"]
      ai.an$Source.of.Variation[ai.an$ID == "10AEC.2031"] <-"T(Aft) x I"
      ai.an$Denominator[ai.an$ID == "10AEC.2031"]<-"T(Aft)xC+Residuals"
      ai.an$Df.Pooled[ai.an$ID == "10AEC.2031"]<- ac.an$Df[ac.an$ID == "10AEC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
      ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEC.2031"]<- ac.an$Sum.Sq[ac.an$ID == "10AEC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEC.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEC.2031"] / ai.an$Df.Pooled[ai.an$ID == "10AEC.2031"]
      ai.an$F.Value[ai.an$ID == "10AEC.2031"]<-ai.an$Mean.Sq[ai.an$ID == "10AEC.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEC.2031"]
      ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEC.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEC.2031"], ai.an$Df[ai.an$ID == "10AEC.2031"], ai.an$Df.Pooled[ai.an$ID == "10AEC.2031"], lower.tail = F)

      if("eliminate" %in% method == T) {
        # 	F = (MSB xCC / MSpooled), where MSpooled = (SST(Aft) xCC + SSRes) / (dfT(Aft) xCC + dfRes)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0041"
        c.an$Df.Formula[c.an$ID == "03AE.0041"]<-"(b-1)(c-1)"
        c.an$Calculated.From[c.an$ID == "03AE.0041"]<-"sym.C"
        c.an$Mean.Sq[c.an$ID == "03AE.0041"]<-c.an$Sum.Sq[c.an$ID == "03AE.0041"] / c.an$Df[c.an$ID == "03AE.0041"]
        c.an$Source.of.Variation[c.an$ID == "03AE.0041"]<-"B x C"
        c.an$Denominator[c.an$ID == "03AE.0041"]<-"T(Aft)xC+Residual"
        c.an$Df.Pooled[c.an$ID == "03AE.0041"]<- ac.an$Df[ac.an$ID == "10AEC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0041"]<- ac.an$Sum.Sq[ac.an$ID == "10AEC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0041"]<- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0041"] / c.an$Df.Pooled[c.an$ID == "03AE.0041"]
        c.an$F.Value[c.an$ID == "03AE.0041"]<- c.an$Mean.Sq[c.an$ID == "03AE.0041"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0041"]
        c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0041"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0041"], c.an$Df[c.an$ID == "03AE.0041"], c.an$Df.Pooled[c.an$ID == "03AE.0041"], lower.tail = F)

        # 	If varB xCC <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SSRes) / (dfB xCC + dfRes)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"]<-"03AE.0031"
        i.an$Df.Formula[i.an$ID == "03AE.0031"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0031"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0031"]<-i.an$Sum.Sq[i.an$ID == "03AE.0031"] / i.an$Df[i.an$ID == "03AE.0031"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0031"]<-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0031"]<-"BxC+Residuals"
        i.an$Df.Pooled[i.an$ID == "03AE.0031"]<- c.an$Df[c.an$ID == "03AE.0040"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"] / i.an$Df.Pooled[i.an$ID == "03AE.0031"]
        i.an$F.Value[i.an$ID == "03AE.0031"]<-i.an$Mean.Sq[i.an$ID == "03AE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0031"], i.an$Df[i.an$ID == "03AE.0031"], i.an$Df.Pooled[i.an$ID == "03AE.0031"], lower.tail = F)
      }
      if("pool" %in% method == T) {
        # 	If varB xCC <= 0 and varT(Aft) xIC > 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SST(Aft) xIC) / (dfB xCC + dfT(Aft) xIC) = SS1/df1
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"]<-"03AE.0031"
        i.an$Df.Formula[i.an$ID == "03AE.0031"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0031"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0031"]<-i.an$Sum.Sq[i.an$ID == "03AE.0031"] / i.an$Df[i.an$ID == "03AE.0031"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0031"]<-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0031"]<-"BxC+T(Aft)xI"
        i.an$Df.Pooled[i.an$ID == "03AE.0031"]<- c.an$Df[c.an$ID == "03AE.0040"] + ai.an$Df[ai.an$ID == "10AEC.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + ai.an$Sum.Sq[ai.an$ID == "10AEC.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"]<- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"] / i.an$Df.Pooled[i.an$ID == "03AE.0031"]
        i.an$F.Value[i.an$ID == "03AE.0031"]<-i.an$Mean.Sq[i.an$ID == "03AE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0031"], i.an$Df[i.an$ID == "03AE.0031"], i.an$Df.Pooled[i.an$ID == "03AE.0031"], lower.tail = F)

        # 	If varB xCC <= 0 and varT(Aft) xIC <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"]<-"03AE.0032"
        i.an$Df.Formula[i.an$ID == "03AE.0032"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0032"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0032"]<-i.an$Sum.Sq[i.an$ID == "03AE.0032"] / i.an$Df[i.an$ID == "03AE.0032"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0032"]<-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0032"]<-"BxC+T(Aft)xI+Residuals"
        i.an$Df.Pooled[i.an$ID == "03AE.0032"]<- i.an$Df.Pooled[i.an$ID == "03AE.0031"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0032"] / i.an$Df.Pooled[i.an$ID == "03AE.0032"]
        i.an$F.Value[i.an$ID == "03AE.0032"]<-i.an$Mean.Sq[i.an$ID == "03AE.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0032"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0032"], i.an$Df[i.an$ID == "03AE.0032"], i.an$Df.Pooled[i.an$ID == "03AE.0032"], lower.tail = F)
      }
    }
  }



  ###---------------- INTERACTIONS with THREE temporal and ONE spatial factors

  if(aav.model == "AEBC") {
    # 	F = (MS_T(P(Aft)) xCC / MS_Res)
    di10AEBC.2040<-cbind(Source.of.Variation = "B:E:C",
                         sym.C[sym.C$Source.of.Variation == "A:E:B:C",2:3] - sym.CB[sym.CB$Source.of.Variation == "B:E:C",2:3])
    ac.an<-dplyr::bind_rows(ac.an,di10AEBC.2040)
    ac.an$ID[ac.an$Source.of.Variation == "B:E:C"] <- "10AEBC.2040"
    ac.an$Df.Formula[ac.an$ID == "10AEBC.2040"]<-"p(t-1)(c-1)"
    ac.an$Calculated.From[ac.an$ID == "10AEBC.2040"]<-"sym.C-sym.CB" #AMONG Controls After
    ac.an$Mean.Sq[ac.an$ID == "10AEBC.2040"] <- ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] / ac.an$Df[ac.an$ID == "10AEBC.2040"]
    ac.an$Source.of.Variation[ac.an$ID == "10AEBC.2040"] <-"T(P(Aft)) x C"
    ac.an$Denominator[ac.an$ID == "10AEBC.2040"]<-"Residuals"
    ac.an$F.Value[ac.an$ID == "10AEBC.2040"]<-ac.an$Mean.Sq[ac.an$ID == "10AEBC.2040"] / f.an$Mean.Sq[f.an$ID == "12R.0000"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "10AEBC.2040"]<- stats::pf(ac.an$F.Value[ac.an$ID == "10AEBC.2040"], ac.an$Df[ac.an$ID == "10AEBC.2040"], f.an$Df[f.an$ID == "12R.0000"], lower.tail = F)

    # 	F = (MS_T(P(Aft)) xIC / MS_T(P(Aft)) xCC), or
    di10AEBC.2030<-cbind(Source.of.Variation = "B:E:C",
                         sym.F[sym.F$Source.of.Variation == "A:E:B:C",2:3] -
                           sym.B[sym.B$Source.of.Variation == "B:E:C",2:3] -
                           sym.C[sym.C$Source.of.Variation == "A:E:B:C",2:3] +
                           sym.CB[sym.CB$Source.of.Variation == "B:E:C",2:3])
    ai.an<-dplyr::bind_rows(ai.an,di10AEBC.2030)
    ai.an$ID[ai.an$Source.of.Variation == "B:E:C"] <- "10AEBC.2030"
    ai.an$Df.Formula[ai.an$ID == "10AEBC.2030"]<-"p(t-1)"
    ai.an$Calculated.From[ai.an$ID == "10AEBC.2030"]<-"sym.F-sym.B-sym.C+sym.BC" # #AMONG impact After
    ai.an$Mean.Sq[ai.an$ID == "10AEBC.2030"]<-ai.an$Sum.Sq[ai.an$ID == "10AEBC.2030"] / ai.an$Df[ai.an$ID == "10AEBC.2030"]
    ai.an$Source.of.Variation[ai.an$ID == "10AEBC.2030"] <-"T(P(Aft)) x I"
    ai.an$Denominator[ai.an$ID == "10AEBC.2030"]<-"T(P(Aft))xC"
    ai.an$F.Value[ai.an$ID == "10AEBC.2030"]<-ai.an$Mean.Sq[ai.an$ID == "10AEBC.2030"] / ac.an$Mean.Sq[ac.an$ID == "10AEBC.2040"]
    ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEBC.2030"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEBC.2030"], ai.an$Df[ai.an$ID == "10AEBC.2030"], ac.an$Df[ac.an$ID == "10AEBC.2040"], lower.tail = F)

    # 	F = (MS_P(Aft) xCC / MS_T(P(Aft)) xCC), or
    di08AEB.2040<-cbind(Source.of.Variation = "B:E",
                        sym.C[sym.C$Source.of.Variation == "A:E:B",2:3] - sym.CB[sym.CB$Source.of.Variation == "B:E",2:3])
    ac.an<-dplyr::bind_rows(ac.an,di08AEB.2040)
    ac.an$ID[ac.an$Source.of.Variation == "B:E"] <- "08AEB.2040"
    ac.an$Df.Formula[ac.an$ID == "08AEB.2040"]<-"(p-1)(c-1)"
    ac.an$Calculated.From[ac.an$ID == "08AEB.2040"]<-"sym.C-sym.CB"
    ac.an$Mean.Sq[ac.an$ID == "08AEB.2040"] <- ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] / ac.an$Df[ac.an$ID == "08AEB.2040"]
    ac.an$Source.of.Variation[ac.an$ID == "08AEB.2040"]<-"P(Aft) x C"
    ac.an$Denominator[ac.an$ID == "08AEB.2040"]<-"T(P(Aft))xC"
    ac.an$F.Value[ac.an$ID == "08AEB.2040"]<-ac.an$Mean.Sq[ac.an$ID == "08AEB.2040"] / ac.an$Mean.Sq[ac.an$ID == "10AEBC.2040"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "08AEB.2040"]<- stats::pf(ac.an$F.Value[ac.an$ID == "08AEB.2040"], ac.an$Df[ac.an$ID == "08AEB.2040"], ac.an$Df[ac.an$ID == "10AEBC.2040"], lower.tail = F)

    # 	F = (MS_P(Aft) xIC / MS_P(Aft) xCC), or
    di08AEB.2030<-cbind(Source.of.Variation = "B:E",
                        sym.F[sym.F$Source.of.Variation == "A:E:B",2:3] -
                          sym.B[sym.B$Source.of.Variation == "B:E",2:3] -
                          sym.C[sym.C$Source.of.Variation == "A:E:B",2:3] +
                          sym.CB[sym.CB$Source.of.Variation == "B:E",2:3])
    ai.an<-dplyr::bind_rows(ai.an,di08AEB.2030)
    ai.an$ID[ai.an$Source.of.Variation == "B:E"] <- "08AEB.2030"
    ai.an$Df.Formula[ai.an$ID == "08AEB.2030"]<-"(p-1)"
    ai.an$Calculated.From[ai.an$ID == "08AEB.2030"]<-"sym.F-sym.B-sym.C+sym.BC"
    ai.an$Mean.Sq[ai.an$ID == "08AEB.2030"]<-ai.an$Sum.Sq[ai.an$ID == "08AEB.2030"] / ai.an$Df[ai.an$ID == "08AEB.2030"]
    ai.an$Source.of.Variation[ai.an$ID == "08AEB.2030"] <-"P(Aft) x I"
    ai.an$Denominator[ai.an$ID == "08AEB.2030"]<-"P(Aft)xC"
    ai.an$F.Value[ai.an$ID == "08AEB.2030"]<-ai.an$Mean.Sq[ai.an$ID == "08AEB.2030"] / ac.an$Mean.Sq[ac.an$ID == "08AEB.2040"]
    ai.an$P.Value.Upper.Tail[ai.an$ID == "08AEB.2030"]<- stats::pf(ai.an$F.Value[ai.an$ID == "08AEB.2030"], ai.an$Df[ai.an$ID == "08AEB.2030"], ac.an$Df[ac.an$ID == "08AEB.2040"], lower.tail = F)

    # 	F = (MS_B xCC / MS_P(Aft) xCC), or
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0040"
    c.an$Df.Formula[c.an$ID == "03AE.0040"]<-"(b-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03AE.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03AE.0040"] <- c.an$Sum.Sq[c.an$ID == "03AE.0040"] / c.an$Df[c.an$ID == "03AE.0040"]
    c.an$Source.of.Variation[c.an$ID == "03AE.0040"]<-"B x C"
    c.an$Denominator[c.an$ID == "03AE.0040"]<-"P(Aft)xC"
    c.an$F.Value[c.an$ID == "03AE.0040"]<-c.an$Mean.Sq[c.an$ID == "03AE.0040"] / ac.an$Mean.Sq[ac.an$ID == "08AEB.2040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0040"], c.an$Df[c.an$ID == "03AE.0040"], ac.an$Df[ac.an$ID == "08AEB.2040"], lower.tail = F)

    #F = (MSB xIC / MSB xCC), or
    di03AE.0030<-cbind(Source.of.Variation = "A:E",
                       sym.F[sym.F$Source.of.Variation == "A:E",2:3] - sym.C[sym.C$Source.of.Variation == "A:E",2:3])
    i.an<-dplyr::bind_rows(i.an,di03AE.0030)
    i.an$ID[i.an$Source.of.Variation == "A:E"]<-"03AE.0030"
    i.an$Df.Formula[i.an$ID == "03AE.0030"]<-"(b-1)"
    i.an$Calculated.From[i.an$Source.of.Variation == "A:E"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "03AE.0030"]<-i.an$Sum.Sq[i.an$ID == "03AE.0030"] / i.an$Df[i.an$ID == "03AE.0030"]
    i.an$Source.of.Variation[i.an$ID == "03AE.0030"]<-"B x I"
    i.an$Denominator[i.an$ID == "03AE.0030"]<-"BxC"
    i.an$F.Value[i.an$ID == "03AE.0030"]<-i.an$Mean.Sq[i.an$ID == "03AE.0030"] / c.an$Mean.Sq[c.an$ID == "03AE.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0030"], i.an$Df[i.an$ID == "03AE.0030"], c.an$Df[c.an$ID == "03AE.0040"], lower.tail = F)

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	If varT(P(Aft)) xCC <= 0, then F = (MS_T(P(Aft)) xIC / MSpooled), where MSpooled = (SST(P(Aft)) xCC + SSRes) / (dfT(P(Aft)) xCC + dfRes)
      ai.an<-dplyr::bind_rows(ai.an,di10AEBC.2030)
      ai.an$ID[ai.an$Source.of.Variation == "B:E:C"]<-"10AEBC.2031"
      ai.an$Df.Formula[ai.an$ID == "10AEBC.2031"]<-"(t-1)"
      ai.an$Calculated.From[ai.an$ID == "10AEBC.2031"]<-"sym.F-sym.B-sym.C+sym.BC" # impact after
      ai.an$Mean.Sq[ai.an$ID == "10AEBC.2031"]<-ai.an$Sum.Sq[ai.an$ID == "10AEBC.2031"] / ai.an$Df[ai.an$ID == "10AEBC.2031"]
      ai.an$Source.of.Variation[ai.an$ID == "10AEBC.2031"] <-"T(P(Aft)) x I"
      ai.an$Denominator[ai.an$ID == "10AEBC.2031"]<-"T(P(Aft))xC+Residuals"
      ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"]<- ac.an$Df[ac.an$ID == "10AEBC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
      ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEBC.2031"]<- ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEBC.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEBC.2031"] / ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"]
      ai.an$F.Value[ai.an$ID == "10AEBC.2031"]<-ai.an$Mean.Sq[ai.an$ID == "10AEBC.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEBC.2031"]
      ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEBC.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEBC.2031"], ai.an$Df[ai.an$ID == "10AEBC.2031"], ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"], lower.tail = F)

      if("eliminate" %in% method == T) {
        # 	F = (MS_P(Aft) xCC / MS_pooled), where MS_pooled = (SS_T(P(Aft)) xCC + SS_Res) / (df_T(P(Aft)) xCC + df_Res)
        ac.an<-dplyr::bind_rows(ac.an,di08AEB.2040)
        ac.an$ID[ac.an$Source.of.Variation == "B:E"] <- "08AEB.2041"
        ac.an$Df.Formula[ac.an$ID == "08AEB.2041"]<-"(p-1)(c-1)"
        ac.an$Calculated.From[ac.an$ID == "08AEB.2041"]<-"sym.C-sym.CB"
        ac.an$Mean.Sq[ac.an$ID == "08AEB.2041"] <- ac.an$Sum.Sq[ac.an$ID == "08AEB.2041"] / ac.an$Df[ac.an$ID == "08AEB.2041"]
        ac.an$Source.of.Variation[ac.an$ID == "08AEB.2041"]<-"P(Aft) x C"
        ac.an$Denominator[ac.an$ID == "08AEB.2041"]<-"T(P(Aft))xC+Residuals"
        ac.an$Df.Pooled[ac.an$ID == "08AEB.2041"]<- ac.an$Df[ac.an$ID == "10AEBC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ac.an$Sum.Sq.Pooled[ac.an$ID == "08AEB.2041"]<- ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ac.an$Mean.Sq.Pooled[ac.an$ID == "08AEB.2041"] <- ac.an$Sum.Sq.Pooled[ac.an$ID == "08AEB.2041"] / ac.an$Df.Pooled[ac.an$ID == "08AEB.2041"]
        ac.an$F.Value[ac.an$ID == "08AEB.2041"]<-ac.an$Mean.Sq[ac.an$ID == "08AEB.2041"] / ac.an$Mean.Sq.Pooled[ac.an$ID == "08AEB.2041"]
        ac.an$P.Value.Upper.Tail[ac.an$ID == "08AEB.2041"]<- stats::pf(ac.an$F.Value[ac.an$ID == "08AEB.2041"], ac.an$Df[ac.an$ID == "08AEB.2041"], ac.an$Df.Pooled[ac.an$ID == "08AEB.2041"], lower.tail = F)

        # 	If varP(Aft) xCC <= 0, then F = (MS_P(Aft) xIC / MS_pooled), where MS_pooled = (SS_P(Aft) xCC + SS_Res) / (df_P(Aft) xCC + df_Res)
        ai.an<-dplyr::bind_rows(ai.an,di08AEB.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E"] <- "08AEB.2031"
        ai.an$Df.Formula[ai.an$ID == "08AEB.2031"]<-"(p-1)"
        ai.an$Calculated.From[ai.an$ID == "08AEB.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "08AEB.2031"]<- ai.an$Sum.Sq[ai.an$ID == "08AEB.2031"] / ai.an$Df[ai.an$ID == "08AEB.2031"]
        ai.an$Source.of.Variation[ai.an$ID == "08AEB.2031"] <-"P(Aft) x I"
        ai.an$Denominator[ai.an$ID == "08AEB.2031"]<-"P(Aft)xC+Residuals"
        ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"]<- ac.an$Df[ac.an$ID == "08AEB.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2031"]<- ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2031"] / ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"]
        ai.an$F.Value[ai.an$ID == "08AEB.2031"]<-ai.an$Mean.Sq[ai.an$ID == "08AEB.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2031"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "08AEB.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "08AEB.2031"], ai.an$Df[ai.an$ID == "08AEB.2031"], ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"], lower.tail = F)

        # 	F = (MS_B xCC / MS_pooled), where MS_pooled = (SS_P(Aft) xCC + SS_Res) / (df_P(Aft) xCC + df_Res)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0041"
        c.an$Df.Formula[c.an$ID == "03AE.0041"]<-"(b-1)(c-1)"
        c.an$Calculated.From[c.an$ID == "03AE.0041"]<-"sym.C"
        c.an$Mean.Sq[c.an$ID == "03AE.0041"] <- c.an$Sum.Sq[c.an$ID == "03AE.0041"] / c.an$Df[c.an$ID == "03AE.0041"]
        c.an$Source.of.Variation[c.an$ID == "03AE.0041"]<-"B x C"
        c.an$Denominator[c.an$ID == "03AE.0041"]<-"P(Aft)xC+Residuals"
        c.an$Df.Pooled[c.an$ID == "03AE.0041"]<- ac.an$Df[ac.an$ID == "08AEB.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0041"]<- ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0041"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0041"] / c.an$Df.Pooled[c.an$ID == "03AE.0041"]
        c.an$F.Value[c.an$ID == "03AE.0041"]<-c.an$Mean.Sq[c.an$ID == "03AE.0041"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0041"]
        c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0041"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0041"], c.an$Df[c.an$ID == "03AE.0041"], c.an$Df.Pooled[c.an$ID == "03AE.0041"], lower.tail = F)

        # 	If varB xCC <= 0, then F = (MS_B xIC / MSpooled), where MSpooled = (SSB xCC + SSRes) / (dfB xCC + dfRes)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"]<-"03AE.0031"
        i.an$Df.Formula[i.an$ID == "03AE.0031"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0031"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0031"]<-i.an$Sum.Sq[i.an$ID == "03AE.0031"] / i.an$Df[i.an$ID == "03AE.0031"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0031"]<-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0031"]<-"BxC+Residuals"
        i.an$Df.Pooled[i.an$ID == "03AE.0031"]<- c.an$Df[c.an$ID == "03AE.0040"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"] / i.an$Df.Pooled[i.an$ID == "03AE.0031"]
        i.an$F.Value[i.an$ID == "03AE.0031"]<-i.an$Mean.Sq[i.an$ID == "03AE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0031"], i.an$Df[i.an$ID == "03AE.0031"], i.an$Df.Pooled[i.an$ID == "03AE.0031"], lower.tail = F)
      }
      if("pool" %in% method == T) {
        # 	If varP(Aft) xCC <= 0 and varT(P(Aft)) xIC > 0, then F = (MSP(Aft) xIC / MSpooled), where MSpooled = (SSP(Aft) xCC + SST(P(Aft)) xIC) / (dfP(Aft) xCC + dfT(P(Aft)) xIC) = SS1/df1
        ai.an<-dplyr::bind_rows(ai.an,di08AEB.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E"] <- "08AEB.2031"
        ai.an$Df.Formula[ai.an$ID == "08AEB.2031"]<-"(p-1)"
        ai.an$Calculated.From[ai.an$ID == "08AEB.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "08AEB.2031"]<-ai.an$Sum.Sq[ai.an$ID == "08AEB.2031"] / ai.an$Df[ai.an$ID == "08AEB.2031"]
        ai.an$Source.of.Variation[ai.an$ID == "08AEB.2031"] <-"P(Aft) x I"
        ai.an$Denominator[ai.an$ID == "08AEB.2031"]<-"P(Aft)xC+T(P(Aft))xI"
        ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"]<- ac.an$Df[ac.an$ID == "08AEB.2040"] + ai.an$Df[ai.an$ID == "10AEBC.2030"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2031"]<- ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + ai.an$Sum.Sq[ai.an$ID == "10AEBC.2030"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2031"] / ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"]
        ai.an$F.Value[ai.an$ID == "08AEB.2031"]<-ai.an$Mean.Sq[ai.an$ID == "08AEB.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2031"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "08AEB.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "08AEB.2031"], ai.an$Df[ai.an$ID == "08AEB.2031"], ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"], lower.tail = F)

        # 	If varP(Aft) xCC <= 0 and varT(P(Aft)) xIC <= 0, then F = (MSP(Aft) xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        ai.an<-dplyr::bind_rows(ai.an,di08AEB.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E"] <- "08AEB.2032"
        ai.an$Df.Formula[ai.an$ID == "08AEB.2032"]<-"(p-1)"
        ai.an$Calculated.From[ai.an$ID == "08AEB.2032"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "08AEB.2032"]<-ai.an$Sum.Sq[ai.an$ID == "08AEB.2032"] / ai.an$Df[ai.an$ID == "08AEB.2032"]
        ai.an$Source.of.Variation[ai.an$ID == "08AEB.2032"] <-"P(Aft) x I"
        ai.an$Denominator[ai.an$ID == "08AEB.2032"]<-"P(Aft)xC+T(P(Aft))xI+Residuals"
        ai.an$Df.Pooled[ai.an$ID == "08AEB.2032"]<- ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"] + f.an$Df[f.an$ID == "12R.0000"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2032"]<- ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2032"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2032"] / ai.an$Df.Pooled[ai.an$ID == "08AEB.2032"]
        ai.an$F.Value[ai.an$ID == "08AEB.2032"]<-ai.an$Mean.Sq[ai.an$ID == "08AEB.2032"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2032"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "08AEB.2032"]<- stats::pf(ai.an$F.Value[ai.an$ID == "08AEB.2032"], ai.an$Df[ai.an$ID == "08AEB.2032"], ai.an$Df.Pooled[ai.an$ID == "08AEB.2032"], lower.tail = F)

        # 	If varB xCC <= 0 and varP(Aft) xIC > 0, then F = (MS_B xIC / MS_pooled), where MSpooled = (SS_B xCC + SS_P(Aft) xIC) / (df_B xCC + df_P(Aft) xIC) = SS1/df1
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"]<-"03AE.0031"
        i.an$Df.Formula[i.an$ID == "03AE.0031"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0031"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0031"]<-i.an$Sum.Sq[i.an$ID == "03AE.0031"] / i.an$Df[i.an$ID == "03AE.0031"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0031"]<-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0031"]<-"BxC+P(Aft)xI"
        i.an$Df.Pooled[i.an$ID == "03AE.0031"]<- c.an$Df[c.an$ID == "03AE.0040"] + ai.an$Df[ai.an$ID == "08AEB.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + ai.an$Sum.Sq[ai.an$ID == "08AEB.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"] / i.an$Df.Pooled[i.an$ID == "03AE.0031"]
        i.an$F.Value[i.an$ID == "03AE.0031"]<-i.an$Mean.Sq[i.an$ID == "03AE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0031"], i.an$Df[i.an$ID == "03AE.0031"], i.an$Df.Pooled[i.an$ID == "03AE.0031"], lower.tail = F)

        # 	If varB xCC <= 0, varP(Aft) xIC <= 0 and varT(P(Aft)) xIC > 0, then F = (MS_B xIC / MS_pooled), where MS_pooled = (SS1 + SS_T(P(Aft)) xIC) / (df1 + df_T(P(Aft)) xIC) = SS2/df2
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"]<-"03AE.0032"
        i.an$Df.Formula[i.an$ID == "03AE.0032"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0032"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0032"]<-i.an$Sum.Sq[i.an$ID == "03AE.0032"] / i.an$Df[i.an$ID == "03AE.0032"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0032"]<-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0032"]<-"BxC+P(Aft)xI+T(P(Aft))xI"
        i.an$Df.Pooled[i.an$ID == "03AE.0032"]<- i.an$Df.Pooled[i.an$ID == "03AE.0031"] + ai.an$Df[ai.an$ID == "10AEBC.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"] + ai.an$Sum.Sq[ai.an$ID == "10AEBC.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0032"] / i.an$Df.Pooled[i.an$ID == "03AE.0032"]
        i.an$F.Value[i.an$ID == "03AE.0032"]<-i.an$Mean.Sq[i.an$ID == "03AE.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0032"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0032"], i.an$Df[i.an$ID == "03AE.0032"], i.an$Df.Pooled[i.an$ID == "03AE.0032"], lower.tail = F)

        # 	If varB xCC <= 0, varP(Aft) xIC <= 0 and varT(P(Aft)) xIC <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SS2 + SSRes) / (df2 + dfRes)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"]<-"03AE.0033"
        i.an$Df.Formula[i.an$ID == "03AE.0033"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0033"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0033"]<-i.an$Sum.Sq[i.an$ID == "03AE.0033"] / i.an$Df[i.an$ID == "03AE.0033"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0033"]<-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0033"]<-"BxC+P(Aft)xI+T(P(Aft))xI+Residuals"
        i.an$Df.Pooled[i.an$ID == "03AE.0033"]<- i.an$Df.Pooled[i.an$ID == "03AE.0032"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0033"]<- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0032"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0033"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0033"] / i.an$Df.Pooled[i.an$ID == "03AE.0033"]
        i.an$F.Value[i.an$ID == "03AE.0033"]<-i.an$Mean.Sq[i.an$ID == "03AE.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0033"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0033"], i.an$Df[i.an$ID == "03AE.0033"], i.an$Df.Pooled[i.an$ID == "03AE.0033"], lower.tail = F)
      }
    }
  }



  ###------------------ INTERACTIONS with TWO temporal and TWO spatial factors

  if(aav.model == "AEGC") {
    #	F = (MS_T(Aft) xS(CC) / MS_Res)
    di11AEGC.2040<-cbind(Source.of.Variation = "C:E:G",
                         sym.C[sym.C$Source.of.Variation == "A:E:G:C",2:3] - sym.CB[sym.CB$Source.of.Variation == "C:E:G",2:3])
    ac.an<-dplyr::bind_rows(ac.an,di11AEGC.2040)
    ac.an$ID[ac.an$Source.of.Variation == "C:E:G"] <- "11AEGC.2040"
    ac.an$Df.Formula[ac.an$ID == "11AEGC.2040"]<-"(t-1)c(s-1)"
    ac.an$Calculated.From[ac.an$ID == "11AEGC.2040"]<-"sym.C-sym.CB"
    ac.an$Mean.Sq[ac.an$ID == "11AEGC.2040"] <- ac.an$Sum.Sq[ac.an$ID == "11AEGC.2040"] / ac.an$Df[ac.an$ID == "11AEGC.2040"]
    ac.an$Source.of.Variation[ac.an$ID == "11AEGC.2040"] <-"T(Aft) x S(C)"
    ac.an$Denominator[ac.an$ID == "11AEGC.2040"]<-"Residuals"
    ac.an$F.Value[ac.an$ID == "11AEGC.2040"]<-ac.an$Mean.Sq[ac.an$ID == "11AEGC.2040"] / f.an$Mean.Sq[f.an$ID == "12R.0000"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "11AEGC.2040"]<- stats::pf(ac.an$F.Value[ac.an$ID == "11AEGC.2040"], ac.an$Df[ac.an$ID == "11AEGC.2040"], f.an$Df[f.an$ID == "12R.0000"], lower.tail = F)

    #  	F = (MST(Aft) xS(IC) / MST(Aft) xS(CC)), or
    di11AEGC.2030<-cbind(Source.of.Variation = "C:E:G",
                         sym.F[sym.F$Source.of.Variation == "A:E:G:C",2:3] -
                           sym.B[sym.B$Source.of.Variation == "C:E:G",2:3] -
                           sym.C[sym.C$Source.of.Variation == "A:E:G:C",2:3] +
                           sym.CB[sym.CB$Source.of.Variation == "C:E:G",2:3])
    ai.an<-dplyr::bind_rows(ai.an,di11AEGC.2030)
    ai.an$ID[ai.an$Source.of.Variation == "C:E:G"] <- "11AEGC.2030"
    ai.an$Df.Formula[ai.an$ID == "11AEGC.2030"]<-"(t-1)(s-1)"
    ai.an$Calculated.From[ai.an$ID == "11AEGC.2030"]<-"sym.F-sym.B-sym.C+sym.BC"
    ai.an$Mean.Sq[ai.an$ID == "11AEGC.2030"]<-ai.an$Sum.Sq[ai.an$ID == "11AEGC.2030"] / ai.an$Df[ai.an$ID == "11AEGC.2030"]
    ai.an$Source.of.Variation[ai.an$ID == "11AEGC.2030"] <-"T(Aft) x S(I)"
    ai.an$Denominator[ai.an$ID == "11AEGC.2030"]<-"T(Aft)xS(C)"
    ai.an$F.Value[ai.an$ID == "11AEGC.2030"]<-ai.an$Mean.Sq[ai.an$ID == "11AEGC.2030"] / ac.an$Mean.Sq[ac.an$ID == "11AEGC.2040"]
    ai.an$P.Value.Upper.Tail[ai.an$ID == "11AEGC.2030"]<- stats::pf(ai.an$F.Value[ai.an$ID == "11AEGC.2030"], ai.an$Df[ai.an$ID == "11AEGC.2030"], ac.an$Df[ac.an$ID == "11AEGC.2040"], lower.tail = F)

    #  	F = (MST(Aft) xCC / MS T(Aft) xS(CC)), or
    di10AEC.2040<-cbind(Source.of.Variation = "C:E",
                        sym.C[sym.C$Source.of.Variation == "A:E:C",2:3] - sym.CB[sym.CB$Source.of.Variation == "C:E",2:3])
    ac.an<-dplyr::bind_rows(ac.an,di10AEC.2040)
    ac.an$ID[ac.an$Source.of.Variation == "C:E"] <- "10AEC.2040"
    ac.an$Df.Formula[ac.an$ID == "10AEC.2040"]<-"(t-1)(c-1)"
    ac.an$Calculated.From[ac.an$ID == "10AEC.2040"]<-"sym.C-sym.CB"
    ac.an$Mean.Sq[ac.an$ID == "10AEC.2040"] <- ac.an$Sum.Sq[ac.an$ID == "10AEC.2040"] / ac.an$Df[ac.an$ID == "10AEC.2040"]
    ac.an$Source.of.Variation[ac.an$ID == "10AEC.2040"] <-"T(Aft) x C"
    ac.an$Denominator[ac.an$ID == "10AEC.2040"]<-"T(Aft)xS(C)"
    ac.an$F.Value[ac.an$ID == "10AEC.2040"]<-ac.an$Mean.Sq[ac.an$ID == "10AEC.2040"] / ac.an$Mean.Sq[ac.an$ID == "11AEGC.2040"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "10AEC.2040"]<- stats::pf(ac.an$F.Value[ac.an$ID == "10AEC.2040"], ac.an$Df[ac.an$ID == "10AEC.2040"], ac.an$Df[ac.an$ID == "11AEGC.2040"], lower.tail = F)

    #F = (MST(Aft) xIC / MST(Aft) xCC), or
    di10AEC.2030<-cbind(Source.of.Variation = "C:E",
                        sym.F[sym.F$Source.of.Variation == "A:E:C",2:3] -
                          sym.B[sym.B$Source.of.Variation == "C:E",2:3] -
                          sym.C[sym.C$Source.of.Variation == "A:E:C",2:3] +
                          sym.CB[sym.CB$Source.of.Variation == "C:E",2:3])
    ai.an<-dplyr::bind_rows(ai.an,di10AEC.2030)
    ai.an$ID[ai.an$Source.of.Variation == "C:E"] <- "10AEC.2030"
    ai.an$Df.Formula[ai.an$ID == "10AEC.2030"]<-"(t-1)"
    ai.an$Calculated.From[ai.an$ID == "10AEC.2030"]<-"sym.F-sym.B-sym.C+sym.BC"
    ai.an$Mean.Sq[ai.an$ID == "10AEC.2030"]<-ai.an$Sum.Sq[ai.an$ID == "10AEC.2030"] / ai.an$Df[ai.an$ID == "10AEC.2030"]
    ai.an$Source.of.Variation[ai.an$ID == "10AEC.2030"] <-"T(Aft) x I"
    ai.an$Denominator[ai.an$ID == "10AEC.2030"]<-"T(Aft)xC"
    ai.an$F.Value[ai.an$ID == "10AEC.2030"]<-ai.an$Mean.Sq[ai.an$ID == "10AEC.2030"] / ac.an$Mean.Sq[ac.an$ID == "10AEC.2040"]
    ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEC.2030"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEC.2030"], ai.an$Df[ai.an$ID == "10AEC.2030"], ac.an$Df[ac.an$ID == "10AEC.2040"], lower.tail = F)

    #  	F = (MSB xS(CC) / MST(Aft) xS(CC)), or
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E:G",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E:G"] <- "07AEG.0040"
    c.an$Df.Formula[c.an$ID == "07AEG.0040"]<-"(b-1)c(s-1)"
    c.an$Calculated.From[c.an$ID == "07AEG.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "07AEG.0040"] <- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] / c.an$Df[c.an$ID == "07AEG.0040"]
    c.an$Source.of.Variation[c.an$ID == "07AEG.0040"] <-"B x S(C)"
    c.an$Denominator[c.an$ID == "07AEG.0040"]<-"T(Aft)xS(C)"
    c.an$F.Value[c.an$ID == "07AEG.0040"]<-c.an$Mean.Sq[c.an$ID == "07AEG.0040"] / ac.an$Mean.Sq[ac.an$ID == "11AEGC.2040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "07AEG.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "07AEG.0040"], c.an$Df[c.an$ID == "07AEG.0040"], ac.an$Df[ac.an$ID == "11AEGC.2040"], lower.tail = F)

    # 	F = (MSB xS(IC) / MSB xS(CC)), or
    di07AEG.0030<-cbind(Source.of.Variation = "A:E:G",
                        sym.F[sym.F$Source.of.Variation == "A:E:G",2:3] -
                          sym.C[sym.C$Source.of.Variation == "A:E:G",2:3])
    i.an<-dplyr::bind_rows(i.an,di07AEG.0030)
    i.an$ID[i.an$Source.of.Variation == "A:E:G"] <- "07AEG.0030"
    i.an$Df.Formula[i.an$ID == "07AEG.0030"]<-"(b-1)(s-1)"
    i.an$Calculated.From[i.an$ID == "07AEG.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "07AEG.0030"]<-i.an$Sum.Sq[i.an$ID == "07AEG.0030"] / i.an$Df[i.an$ID == "07AEG.0030"]
    i.an$Source.of.Variation[i.an$ID == "07AEG.0030"] <-"B x S(I)"
    i.an$Denominator[i.an$ID == "07AEG.0030"]<-"BxS(C)"
    i.an$F.Value[i.an$ID == "07AEG.0030"]<-i.an$Mean.Sq[i.an$ID == "07AEG.0030"] / c.an$Mean.Sq[c.an$ID == "07AEG.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "07AEG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "07AEG.0030"], i.an$Df[i.an$ID == "07AEG.0030"], c.an$Df[c.an$ID == "07AEG.0040"], lower.tail = F)

    # 	If varB xS(CC) <= 0 and varT(Aft) xCC > 0, then F = (MSB xCC / MST(Aft) xCC)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0040"
    c.an$Df.Formula[c.an$ID == "03AE.0040"]<-"(b-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03AE.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03AE.0040"] <- c.an$Sum.Sq[c.an$ID == "03AE.0040"] / c.an$Df[c.an$ID == "03AE.0040"]
    c.an$Source.of.Variation[c.an$ID == "03AE.0040"]<-"B x C"
    c.an$Denominator[c.an$ID == "03AE.0040"]<-"T(Aft)xC"
    c.an$F.Value[c.an$ID == "03AE.0040"]<-c.an$Mean.Sq[c.an$ID == "03AE.0040"] / ac.an$Mean.Sq[ac.an$ID == "10AEC.2040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0040"], c.an$Df[c.an$ID == "03AE.0040"], ac.an$Df[ac.an$ID == "10AEC.2040"], lower.tail = F)

    # 	If varT(Aft) xCC <= 0 and varB xS(CC) > 0, then F = (MSB xCC / MSB xS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0041"
    c.an$Df.Formula[c.an$ID == "03AE.0041"]<-"(b-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03AE.0041"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03AE.0041"] <- c.an$Sum.Sq[c.an$ID == "03AE.0041"] / c.an$Df[c.an$ID == "03AE.0041"]
    c.an$Source.of.Variation[c.an$ID == "03AE.0041"]<-"B x C"
    c.an$Denominator[c.an$ID == "03AE.0041"]<-"BxS(C)"
    c.an$F.Value[c.an$ID == "03AE.0041"]<-c.an$Mean.Sq[c.an$ID == "03AE.0041"] / c.an$Mean.Sq[c.an$ID == "07AEG.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0041"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0041"], c.an$Df[c.an$ID == "03AE.0041"], c.an$Df[c.an$ID == "07AEG.0040"], lower.tail = F)

    # 	If varT(Aft) xCC <= 0, varB xS(CC) <= 0, then F = (MSB xCC / MSpooled), where MSpooled = (SS_T(Aft) xCC + SS_B xS(CC) + SS_T(Aft) xS(CC)) / (df_T(Aft) xCC + df_B xS(CC) + df_T(Aft) xS(CC)), or
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0042"
    c.an$Df.Formula[c.an$ID == "03AE.0042"]<-"(b-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03AE.0042"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03AE.0042"] <- c.an$Sum.Sq[c.an$ID == "03AE.0042"] / c.an$Df[c.an$ID == "03AE.0042"]
    c.an$Source.of.Variation[c.an$ID == "03AE.0042"]<-"B x C"
    c.an$Denominator[c.an$ID == "03AE.0042"]<-"T(Aft)xC+BxS(C)+T(Aft)xS(C)"
    c.an$Df.Pooled[c.an$ID == "03AE.0042"]<- c.an$Df[c.an$ID == "07AEG.0040"] + ac.an$Df[ac.an$ID == "10AEC.2040"] + ac.an$Df[ac.an$ID == "11AEGC.2040"]
    c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0042"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + ac.an$Sum.Sq[ac.an$ID == "10AEC.2040"] + ac.an$Sum.Sq[ac.an$ID == "11AEGC.2040"]
    c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0042"]<- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0042"] / c.an$Df.Pooled[c.an$ID == "03AE.0042"]
    c.an$F.Value[c.an$ID == "03AE.0042"]<- c.an$Mean.Sq[c.an$ID == "03AE.0042"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0042"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0042"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0042"], c.an$Df[c.an$ID == "03AE.0042"], c.an$Df.Pooled[c.an$ID == "03AE.0042"], lower.tail = F)

    # 	F = (MSB xIC / MSB xCC), or
    di03AE.0030<-cbind(Source.of.Variation = "A:E",
                       sym.F[sym.F$Source.of.Variation == "A:E",2:3] -
                         sym.C[sym.C$Source.of.Variation == "A:E",2:3])
    i.an<-dplyr::bind_rows(i.an,di03AE.0030)
    i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0030"
    i.an$Df.Formula[i.an$ID == "03AE.0030"]<-"(b-1)"
    i.an$Calculated.From[i.an$ID == "03AE.0030"]<-"sym.F-sym.C"
    i.an$Source.of.Variation[i.an$ID == "03AE.0030"] <-"B x I"
    i.an$Denominator[i.an$ID == "03AE.0030"]<-"BxC"
    i.an$Mean.Sq[i.an$ID == "03AE.0030"]<-i.an$Sum.Sq[i.an$ID == "03AE.0030"] / i.an$Df[i.an$ID == "03AE.0030"]
    i.an$F.Value[i.an$ID == "03AE.0030"]<-i.an$Mean.Sq[i.an$ID == "03AE.0030"] / c.an$Mean.Sq[c.an$ID == "03AE.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0030"], i.an$Df[i.an$ID == "03AE.0030"], c.an$Df[c.an$ID == "03AE.0040"], lower.tail = F)

    ### POOLING STRUCTURE
    if(pooling == T) {
      #  	If varT(Aft) xS(CC) <= 0, then F = (MST(Aft) xS(IC) / MSpooled), where MSpooled = (SST(Aft) xS(CC) + SSRes) / (dfT(Aft) xS(CC) + dfRes)
      ai.an<-dplyr::bind_rows(ai.an,di11AEGC.2030)
      ai.an$ID[ai.an$Source.of.Variation == "C:E:G"] <- "11AEGC.2031"
      ai.an$Df.Formula[ai.an$ID == "11AEGC.2031"]<-"(t-1)(s-1)"
      ai.an$Calculated.From[ai.an$ID == "11AEGC.2031"]<-"sym.F-sym.C"
      ai.an$Mean.Sq[ai.an$ID == "11AEGC.2031"]<-ai.an$Sum.Sq[ai.an$ID == "11AEGC.2031"] / ai.an$Df[ai.an$ID == "11AEGC.2031"]
      ai.an$Source.of.Variation[ai.an$ID == "11AEGC.2031"] <-"T(Aft) x S(I)"
      ai.an$Denominator[ai.an$ID == "11AEGC.2031"]<-"T(Aft)xS(C)+Residuals"
      ai.an$Df.Pooled[ai.an$ID == "11AEGC.2031"]<- ac.an$Df[ac.an$ID == "11AEGC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
      ai.an$Sum.Sq.Pooled[ai.an$ID == "11AEGC.2031"]<- ac.an$Sum.Sq[ac.an$ID == "11AEGC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      ai.an$Mean.Sq.Pooled[ai.an$ID == "11AEGC.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "11AEGC.2031"] / ai.an$Df.Pooled[ai.an$ID == "11AEGC.2031"]
      ai.an$F.Value[ai.an$ID == "11AEGC.2031"]<-ai.an$Mean.Sq[ai.an$ID == "11AEGC.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "11AEGC.2031"]
      ai.an$P.Value.Upper.Tail[ai.an$ID == "11AEGC.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "11AEGC.2031"], ai.an$Df[ai.an$ID == "11AEGC.2031"], ai.an$Df.Pooled[ai.an$ID == "11AEGC.2031"], lower.tail = F)

      if("eliminate" %in% method == T) {
        # 	F = (MST(Aft) xCC / MSpooled), where MSpooled = (SST(Aft) xS(CC) + SSRes) / (dfT(Aft) xS(CC) + dfRes)
        ac.an<-dplyr::bind_rows(ac.an,di10AEC.2040)
        ac.an$ID[ac.an$Source.of.Variation == "C:E"] <- "10AEC.2041"
        ac.an$Df.Formula[ac.an$ID == "10AEC.2041"]<-"(t-1)(c-1)"
        ac.an$Calculated.From[ac.an$ID == "10AEC.2041"]<-"sym.C-sym.CB"
        ac.an$Mean.Sq[ac.an$ID == "10AEC.2041"] <- ac.an$Sum.Sq[ac.an$ID == "10AEC.2041"] / ac.an$Df[ac.an$ID == "10AEC.2041"]
        ac.an$Source.of.Variation[ac.an$ID == "10AEC.2041"] <-"T(Aft) x C"
        ac.an$Denominator[ac.an$ID == "10AEC.2041"]<-"T(Aft)xS(C)+Residuals"
        ac.an$Df.Pooled[ac.an$ID == "10AEC.2041"]<- ac.an$Df[ac.an$ID == "11AEGC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ac.an$Sum.Sq.Pooled[ac.an$ID == "10AEC.2041"]<- ac.an$Sum.Sq[ac.an$ID == "11AEGC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ac.an$Mean.Sq.Pooled[ac.an$ID == "10AEC.2041"] <- ac.an$Sum.Sq.Pooled[ac.an$ID == "10AEC.2041"] / ac.an$Df.Pooled[ac.an$ID == "10AEC.2041"]
        ac.an$F.Value[ac.an$ID == "10AEC.2041"]<-ac.an$Mean.Sq[ac.an$ID == "10AEC.2041"] / ac.an$Mean.Sq.Pooled[ac.an$ID == "10AEC.2041"]
        ac.an$P.Value.Upper.Tail[ac.an$ID == "10AEC.2041"]<- stats::pf(ac.an$F.Value[ac.an$ID == "10AEC.2041"], ac.an$Df[ac.an$ID == "10AEC.2041"], ac.an$Df.Pooled[ac.an$ID == "10AEC.2041"], lower.tail = F)

        # 	If varT(Aft) xCC <= 0, then F = (MST(Aft) xIC / MSpooled), where MSpooled = (SST(Aft) xCC + SSRes) / (dfT(Aft) xCC + dfRes)
        ai.an<-dplyr::bind_rows(ai.an,di10AEC.2030)
        ai.an$ID[ai.an$Source.of.Variation == "C:E"] <- "10AEC.2031"
        ai.an$Df.Formula[ai.an$ID == "10AEC.2031"]<-"(t-1)"
        ai.an$Calculated.From[ai.an$ID == "10AEC.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "10AEC.2031"]<-ai.an$Sum.Sq[ai.an$ID == "10AEC.2031"] / ai.an$Df[ai.an$ID == "10AEC.2031"]
        ai.an$Source.of.Variation[ai.an$ID == "10AEC.2031"] <-"T(Aft) x I"
        ai.an$Denominator[ai.an$ID == "10AEC.2031"]<-"T(Aft)xC+Residuals"
        ai.an$Df.Pooled[ai.an$ID == "10AEC.2031"]<- ac.an$Df[ac.an$ID == "10AEC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEC.2031"]<- ac.an$Sum.Sq[ac.an$ID == "10AEC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEC.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEC.2031"] / ai.an$Df.Pooled[ai.an$ID == "10AEC.2031"]
        ai.an$F.Value[ai.an$ID == "10AEC.2031"]<-ai.an$Mean.Sq[ai.an$ID == "10AEC.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEC.2031"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEC.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEC.2031"], ai.an$Df[ai.an$ID == "10AEC.2031"], ai.an$Df.Pooled[ai.an$ID == "10AEC.2031"], lower.tail = F)

        # 	F = (MSB xS(CC) / MSpooled), where MSpooled = (SST(Aft) xS(CC) + SSRes) / (dfT(Aft) xS(CC) + dfRes)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E:G",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E:G"] <- "07AEG.0041"
        c.an$Df.Formula[c.an$ID == "07AEG.0041"]<-"(b-1)c(s-1)"
        c.an$Calculated.From[c.an$ID == "07AEG.0041"]<-"sym.C-sym.CB"
        c.an$Mean.Sq[c.an$ID == "07AEG.0041"] <- c.an$Sum.Sq[c.an$ID == "07AEG.0041"] / c.an$Df[c.an$ID == "07AEG.0041"]
        c.an$Source.of.Variation[c.an$ID == "07AEG.0041"] <-"B x S(C)"
        c.an$Denominator[c.an$ID == "07AEG.0041"]<-"T(Aft)xS(C)+Residuals"
        c.an$Df.Pooled[c.an$ID == "07AEG.0041"]<- ac.an$Df[ac.an$ID == "11AEGC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "07AEG.0041"]<- ac.an$Sum.Sq[ac.an$ID == "11AEGC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "07AEG.0041"] <- c.an$Sum.Sq.Pooled[c.an$ID == "07AEG.0041"] / c.an$Df.Pooled[c.an$ID == "07AEG.0041"]
        c.an$F.Value[c.an$ID == "07AEG.0041"]<-c.an$Mean.Sq[c.an$ID == "07AEG.0041"] / c.an$Mean.Sq.Pooled[c.an$ID == "07AEG.0041"]
        c.an$P.Value.Upper.Tail[c.an$ID == "07AEG.0041"]<- stats::pf(c.an$F.Value[c.an$ID == "07AEG.0041"], c.an$Df[c.an$ID == "07AEG.0041"], c.an$Df.Pooled[c.an$ID == "07AEG.0041"], lower.tail = F)

        # 	If varB xS(CC) <= 0, then F = (MSB xS(IC) / MSpooled), where MSpooled = (SSB xS(CC) + SSRes) / (dfB xS(CC) + dfRes)
        i.an<-dplyr::bind_rows(i.an,di07AEG.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E:G"] <- "07AEG.0031"
        i.an$Df.Formula[i.an$ID == "07AEG.0031"]<-"(b-1)(s-1)"
        i.an$Calculated.From[i.an$ID == "07AEG.0031"]<-"sym.F-sym.B-sym.C+sym.BC"
        i.an$Mean.Sq[i.an$ID == "07AEG.0031"]<-i.an$Sum.Sq[i.an$ID == "07AEG.0031"] / i.an$Df[i.an$ID == "07AEG.0031"]
        i.an$Source.of.Variation[i.an$ID == "07AEG.0031"] <-"B x S(I)"
        i.an$Denominator[i.an$ID == "07AEG.0031"]<-"BxS(C)+Residuals"
        i.an$Df.Pooled[i.an$ID == "07AEG.0031"]<- c.an$Df[c.an$ID == "07AEG.0040"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0031"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0031"] / i.an$Df.Pooled[i.an$ID == "07AEG.0031"]
        i.an$F.Value[i.an$ID == "07AEG.0031"]<-i.an$Mean.Sq[i.an$ID == "07AEG.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "07AEG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "07AEG.0031"], i.an$Df[i.an$ID == "07AEG.0031"], i.an$Df.Pooled[i.an$ID == "07AEG.0031"], lower.tail = F)

        # 	If varB xS(CC) <= 0 and varT(Aft) xCC > 0, then F = (MSB xCC / MSpooled), where MSpooled = (SST(Aft) xCC + SSRes) / (dfT(Aft) xCC  + dfRes)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0043"
        c.an$Df.Formula[c.an$ID == "03AE.0043"]<-"(b-1)(c-1)"
        c.an$Calculated.From[c.an$ID == "03AE.0043"]<-"sym.C"
        c.an$Mean.Sq[c.an$ID == "03AE.0043"] <- c.an$Sum.Sq[c.an$ID == "03AE.0043"] / c.an$Df[c.an$ID == "03AE.0043"]
        c.an$Source.of.Variation[c.an$ID == "03AE.0043"]<-"B x C"
        c.an$Denominator[c.an$ID == "03AE.0043"]<-"T(Aft)xC+Residuals"
        c.an$Df.Pooled[c.an$ID == "03AE.0043"]<- ac.an$Df[ac.an$ID == "10AEC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0043"]<- ac.an$Sum.Sq[ac.an$ID == "10AEC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0043"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0043"] / c.an$Df.Pooled[c.an$ID == "03AE.0043"]
        c.an$F.Value[c.an$ID == "03AE.0043"]<-c.an$Mean.Sq[c.an$ID == "03AE.0043"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0043"]
        c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0043"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0043"], c.an$Df[c.an$ID == "03AE.0043"], c.an$Df.Pooled[c.an$ID == "03AE.0043"], lower.tail = F)

        # 	If varT(Aft) xCC <= 0 and varB xS(CC) > 0, then F = (MSB xCC / MSpooled), where MSpooled = (SSB xS(CC) + SSRes) / (dfB xS(CC) + dfRes)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0044"
        c.an$Df.Formula[c.an$ID == "03AE.0044"]<-"(b-1)(c-1)"
        c.an$Calculated.From[c.an$ID == "03AE.0044"]<-"sym.C"
        c.an$Mean.Sq[c.an$ID == "03AE.0044"] <- c.an$Sum.Sq[c.an$ID == "03AE.0044"] / c.an$Df[c.an$ID == "03AE.0044"]
        c.an$Source.of.Variation[c.an$ID == "03AE.0044"]<-"B x C"
        c.an$Denominator[c.an$ID == "03AE.0044"]<-"BxS(C)+Residuals"
        c.an$Df.Pooled[c.an$ID == "03AE.0044"]<- c.an$Df[c.an$ID == "07AEG.0040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0044"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0044"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0044"] / c.an$Df.Pooled[c.an$ID == "03AE.0044"]
        c.an$F.Value[c.an$ID == "03AE.0044"]<-c.an$Mean.Sq[c.an$ID == "03AE.0044"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0044"]
        c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0044"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0044"], c.an$Df[c.an$ID == "03AE.0044"], c.an$Df.Pooled[c.an$ID == "03AE.0044"], lower.tail = F)

        # 	If varT(Aft) xCC <= 0, varB xS(CC) <= 0, then F = (MSB xCC / MSpooled), where MSpooled = (SST(Aft) xCC + SSB xS(CC) + SST(Aft) xS(CC) + SSRes) / (dfT(Aft) xCC + dfB xS(CC) + dfT(Aft) xS(CC) + dfRes)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0045"
        c.an$Df.Formula[c.an$ID == "03AE.0045"]<-"(b-1)(c-1)"
        c.an$Calculated.From[c.an$ID == "03AE.0045"]<-"sym.C"
        c.an$Mean.Sq[c.an$ID == "03AE.0045"] <- c.an$Sum.Sq[c.an$ID == "03AE.0045"] / c.an$Df[c.an$ID == "03AE.0045"]
        c.an$Source.of.Variation[c.an$ID == "03AE.0045"]<-"B x C"
        c.an$Denominator[c.an$ID == "03AE.0045"]<-"T(Aft)xC+BxS(C)+T(Aft)xS(C)+Residuals"
        c.an$Df.Pooled[c.an$ID == "03AE.0045"]<- ac.an$Df[ac.an$ID == "10AEC.2040"] + c.an$Df[c.an$ID == "07AEG.0040"] + ac.an$Df[ac.an$ID == "11AEGC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0045"]<- ac.an$Sum.Sq[ac.an$ID == "10AEC.2040"] + c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + ac.an$Sum.Sq[ac.an$ID == "11AEGC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0045"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0045"] / c.an$Df.Pooled[c.an$ID == "03AE.0045"]
        c.an$F.Value[c.an$ID == "03AE.0045"]<-c.an$Mean.Sq[c.an$ID == "03AE.0045"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0045"]
        c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0045"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0045"], c.an$Df[c.an$ID == "03AE.0045"], c.an$Df.Pooled[c.an$ID == "03AE.0045"], lower.tail = F)

        # 	If varB xCC <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SSRes) / (dfB xCC + dfRes)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0031"
        i.an$Df.Formula[i.an$ID == "03AE.0031"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0031"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0031"]<-i.an$Sum.Sq[i.an$ID == "03AE.0031"] / i.an$Df[i.an$ID == "03AE.0031"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0031"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0031"]<-"BxC+Residuals"
        i.an$Df.Pooled[i.an$ID == "03AE.0031"]<- c.an$Df[c.an$ID == "03AE.0040"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"] / i.an$Df.Pooled[i.an$ID == "03AE.0031"]
        i.an$F.Value[i.an$ID == "03AE.0031"]<-i.an$Mean.Sq[i.an$ID == "03AE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0031"], i.an$Df[i.an$ID == "03AE.0031"], i.an$Df.Pooled[i.an$ID == "03AE.0031"], lower.tail = F)
      }

      if("pool" %in% method == T) {
        # 	If varT(Aft) xCC <= 0 and varT(Aft) xS(IC) > 0, then F = (MST(Aft) xIC / MSpooled), where MSpooled = (SST(Aft) xCC + SST(Aft) xS(IC) / (dfT(Aft) xCC + dfT(Aft) xS(IC) = SS1/df1
        ai.an<-dplyr::bind_rows(ai.an,di10AEC.2030)
        ai.an$ID[ai.an$Source.of.Variation == "C:E"] <- "10AEC.2031"
        ai.an$Df.Formula[ai.an$ID == "10AEC.2031"]<-"(t-1)"
        ai.an$Calculated.From[ai.an$ID == "10AEC.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "10AEC.2031"]<-ai.an$Sum.Sq[ai.an$ID == "10AEC.2031"] / ai.an$Df[ai.an$ID == "10AEC.2031"]
        ai.an$Source.of.Variation[ai.an$ID == "10AEC.2031"] <-"T(Aft) x I"
        ai.an$Denominator[ai.an$ID == "10AEC.2031"]<-"T(Aft)xC+T(Aft)xS(I)"
        ai.an$Df.Pooled[ai.an$ID == "10AEC.2031"]<- ac.an$Df[ac.an$ID == "10AEC.2040"] + ai.an$Df[ai.an$ID == "11AEGC.2030"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEC.2031"]<- ac.an$Sum.Sq[ac.an$ID == "10AEC.2040"] + ai.an$Sum.Sq[ai.an$ID == "11AEGC.2030"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEC.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEC.2031"] / ai.an$Df.Pooled[ai.an$ID == "10AEC.2031"]
        ai.an$F.Value[ai.an$ID == "10AEC.2031"]<-ai.an$Mean.Sq[ai.an$ID == "10AEC.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEC.2031"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEC.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEC.2031"], ai.an$Df[ai.an$ID == "10AEC.2031"], ai.an$Df.Pooled[ai.an$ID == "10AEC.2031"], lower.tail = F)

        # 	If varT(Aft) xCC <= 0 and varT(Aft) xS(IC) <= 0, then F = (MST(Aft) xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        ai.an<-dplyr::bind_rows(ai.an,di10AEC.2030)
        ai.an$ID[ai.an$Source.of.Variation == "C:E"] <- "10AEC.2032"
        ai.an$Df.Formula[ai.an$ID == "10AEC.2032"]<-"(t-1)"
        ai.an$Calculated.From[ai.an$ID == "10AEC.2032"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "10AEC.2032"]<-ai.an$Sum.Sq[ai.an$ID == "10AEC.2032"] / ai.an$Df[ai.an$ID == "10AEC.2032"]
        ai.an$Source.of.Variation[ai.an$ID == "10AEC.2032"] <-"T(Aft) x I"
        ai.an$Denominator[ai.an$ID == "10AEC.2032"]<-"T(Aft)xC+T(Aft)xS(I)+Residuals"
        ai.an$Df.Pooled[ai.an$ID == "10AEC.2032"]<- ai.an$Df.Pooled[ai.an$ID == "10AEC.2031"] + f.an$Df[f.an$ID == "12R.0000"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEC.2032"]<- ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEC.2031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEC.2032"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEC.2032"] / ai.an$Df.Pooled[ai.an$ID == "10AEC.2032"]
        ai.an$F.Value[ai.an$ID == "10AEC.2032"]<-ai.an$Mean.Sq[ai.an$ID == "10AEC.2032"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEC.2032"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEC.2032"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEC.2032"], ai.an$Df[ai.an$ID == "10AEC.2032"], ai.an$Df.Pooled[ai.an$ID == "10AEC.2032"], lower.tail = F)

        # 	If varB xS(CC) <= 0 and varT(Aft) xS(IC) > 0, then F = (MSB xS(IC) / MSpooled), where MSpooled = (SSB xS(CC) + SST(Aft) xS(IC)) / (dfB xS(CC) + dfT(Aft) xS(IC)) = SS1/df1
        i.an<-dplyr::bind_rows(i.an,di07AEG.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E:G"] <- "07AEG.0031"
        i.an$Df.Formula[i.an$ID == "07AEG.0031"]<-"(b-1)(s-1)"
        i.an$Calculated.From[i.an$ID == "07AEG.0031"]<-"sym.F-sym.B-sym.C+sym.BC"
        i.an$Mean.Sq[i.an$ID == "07AEG.0031"]<-i.an$Sum.Sq[i.an$ID == "07AEG.0031"] / i.an$Df[i.an$ID == "07AEG.0031"]
        i.an$Source.of.Variation[i.an$ID == "07AEG.0031"] <-"B x S(I)"
        i.an$Denominator[i.an$ID == "07AEG.0031"]<-"BxS(C)+T(Aft)xS(I)"
        i.an$Df.Pooled[i.an$ID == "07AEG.0031"]<- c.an$Df[c.an$ID == "07AEG.0040"] + ai.an$Df[ai.an$ID == "11AEGC.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0031"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + ai.an$Sum.Sq[ai.an$ID == "11AEGC.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0031"] / i.an$Df.Pooled[i.an$ID == "07AEG.0031"]
        i.an$F.Value[i.an$ID == "07AEG.0031"]<-i.an$Mean.Sq[i.an$ID == "07AEG.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "07AEG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "07AEG.0031"], i.an$Df[i.an$ID == "07AEG.0031"], i.an$Df.Pooled[i.an$ID == "07AEG.0031"], lower.tail = F)

        # 	If varB xS(CC) <= 0 and varT(Aft) xS(IC) <= 0, then F = (MSB xS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        i.an<-dplyr::bind_rows(i.an,di07AEG.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E:G"] <- "07AEG.0032"
        i.an$Df.Formula[i.an$ID == "07AEG.0032"]<-"(b-1)(s-1)"
        i.an$Calculated.From[i.an$ID == "07AEG.0032"]<-"sym.F-sym.B-sym.C+sym.BC"
        i.an$Mean.Sq[i.an$ID == "07AEG.0032"]<-i.an$Sum.Sq[i.an$ID == "07AEG.0032"] / i.an$Df[i.an$ID == "07AEG.0032"]
        i.an$Source.of.Variation[i.an$ID == "07AEG.0032"] <-"B x S(I)"
        i.an$Denominator[i.an$ID == "07AEG.0032"]<-"BxS(C)+T(Aft)xS(I)+Residuals"
        i.an$Df.Pooled[i.an$ID == "07AEG.0032"]<- i.an$Df.Pooled[i.an$ID == "07AEG.0031"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0032"] / i.an$Df.Pooled[i.an$ID == "07AEG.0032"]
        i.an$F.Value[i.an$ID == "07AEG.0032"]<-i.an$Mean.Sq[i.an$ID == "07AEG.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0032"]
        i.an$P.Value.Upper.Tail[i.an$ID == "07AEG.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "07AEG.0032"], i.an$Df[i.an$ID == "07AEG.0032"], i.an$Df.Pooled[i.an$ID == "07AEG.0032"], lower.tail = F)

        # 	If varB xCC <= 0, varB xS(IC) <= 0 and varT(Aft) xIC > 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SST(Aft) xIC) / (dfB xCC + dfT(Aft) xIC)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0031"
        i.an$Df.Formula[i.an$ID == "03AE.0031"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0031"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0031"]<-i.an$Sum.Sq[i.an$ID == "03AE.0031"] / i.an$Df[i.an$ID == "03AE.0031"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0031"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0031"]<-"BxC+T(Aft)xI"
        i.an$Df.Pooled[i.an$ID == "03AE.0031"]<- c.an$Df[c.an$ID == "03AE.0040"] + ai.an$Df[ai.an$ID == "10AEC.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + ai.an$Sum.Sq[ai.an$ID == "10AEC.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"] / i.an$Df.Pooled[i.an$ID == "03AE.0031"]
        i.an$F.Value[i.an$ID == "03AE.0031"]<-i.an$Mean.Sq[i.an$ID == "03AE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0031"], i.an$Df[i.an$ID == "03AE.0031"], i.an$Df.Pooled[i.an$ID == "03AE.0031"], lower.tail = F)

        # 	If varB xCC <= 0, varT(Aft) xIC <= 0 and varB xS(IC) > 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SSB xS(IC)) / (dfB xCC + dfB xS(IC))
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0032"
        i.an$Df.Formula[i.an$ID == "03AE.0032"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0032"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0032"]<-i.an$Sum.Sq[i.an$ID == "03AE.0032"] / i.an$Df[i.an$ID == "03AE.0032"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0032"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0032"]<-"BxC+BxS(I)"
        i.an$Df.Pooled[i.an$ID == "03AE.0032"]<- c.an$Df[c.an$ID == "03AE.0040"] + i.an$Df[i.an$ID == "07AEG.0030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0032"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + i.an$Sum.Sq[i.an$ID == "07AEG.0030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0032"] / i.an$Df.Pooled[i.an$ID == "03AE.0032"]
        i.an$F.Value[i.an$ID == "03AE.0032"]<-i.an$Mean.Sq[i.an$ID == "03AE.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0032"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0032"], i.an$Df[i.an$ID == "03AE.0032"], i.an$Df.Pooled[i.an$ID == "03AE.0032"], lower.tail = F)

        # 	If varB xCC <= 0, varT(Aft) xIC <= 0, varB xS(IC) <= 0 and varT(Aft) xS(IC) > 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SST(Aft) xIC + SSB xS(IC) + SST(Aft) xS(IC)) / (dfB xCC + dfT(Aft) xIC + dfB xS(IC) + dfT(Aft) xS(IC)) = SS1/df1
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0033"
        i.an$Df.Formula[i.an$ID == "03AE.0033"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0033"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0033"]<-i.an$Sum.Sq[i.an$ID == "03AE.0033"] / i.an$Df[i.an$ID == "03AE.0033"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0033"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0033"]<-"BxC+T(Aft)xI+BxS(I)+T(Aft)xS(I)"
        i.an$Df.Pooled[i.an$ID == "03AE.0033"]<- c.an$Df[c.an$ID == "03AE.0040"] + ai.an$Df[ai.an$ID == "10AEC.2030"] + i.an$Df[i.an$ID == "07AEG.0030"] + ai.an$Df[ai.an$ID == "11AEGC.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0033"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + ai.an$Sum.Sq[ai.an$ID == "10AEC.2030"]  + i.an$Sum.Sq[i.an$ID == "07AEG.0030"] + ai.an$Sum.Sq[ai.an$ID == "11AEGC.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0033"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0033"] / i.an$Df.Pooled[i.an$ID == "03AE.0033"]
        i.an$F.Value[i.an$ID == "03AE.0033"]<-i.an$Mean.Sq[i.an$ID == "03AE.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0033"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0033"], i.an$Df[i.an$ID == "03AE.0033"], i.an$Df.Pooled[i.an$ID == "03AE.0033"], lower.tail = F)

        # 	If varB xCC <= 0, varB xS(IC) <= 0, varT(B) xIC <= 0, varT(B) xS(IC) <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0034"
        i.an$Df.Formula[i.an$ID == "03AE.0034"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0034"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0034"]<-i.an$Sum.Sq[i.an$ID == "03AE.0034"] / i.an$Df[i.an$ID == "03AE.0034"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0034"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0034"]<-"BxC+T(Aft)xI+BxS(I)+T(Aft)xS(I)+Residuals"
        i.an$Df.Pooled[i.an$ID == "03AE.0034"]<- i.an$Df.Pooled[i.an$ID == "03AE.0033"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0034"]<- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0033"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0034"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0034"] / i.an$Df.Pooled[i.an$ID == "03AE.0034"]
        i.an$F.Value[i.an$ID == "03AE.0034"]<-i.an$Mean.Sq[i.an$ID == "03AE.0034"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0034"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0034"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0034"], i.an$Df[i.an$ID == "03AE.0034"], i.an$Df.Pooled[i.an$ID == "03AE.0034"], lower.tail = F)
      }
    }
  }



  ###---------------- INTERACTIONS with THREE temporal and TWO spatial factors

  if(aav.model == "AEBGC") {
    # 	F = (MS_T(P(Aft)) xS(CC) / MS_Res)
    di11AEBGC.2040<-cbind(Source.of.Variation = "B:E:G:C",
                          sym.C[sym.C$Source.of.Variation == "A:E:B:G:C",2:3] - sym.CB[sym.CB$Source.of.Variation == "B:E:G:C",2:3])
    ac.an<-dplyr::bind_rows(ac.an,di11AEBGC.2040)
    ac.an$ID[ac.an$Source.of.Variation == "B:E:G:C"] <- "11AEBGC.2040"
    ac.an$Df.Formula[ac.an$ID == "11AEBGC.2040"]<-"p(t-1)c(s-1)"
    ac.an$Calculated.From[ac.an$ID == "11AEBGC.2040"]<-"sym.C-sym.CB"
    ac.an$Mean.Sq[ac.an$ID == "11AEBGC.2040"] <- ac.an$Sum.Sq[ac.an$ID == "11AEBGC.2040"] / ac.an$Df[ac.an$ID == "11AEBGC.2040"]
    ac.an$Source.of.Variation[ac.an$ID == "11AEBGC.2040"] <-"T(P(Aft)) x S(C)"
    ac.an$Denominator[ac.an$ID == "11AEBGC.2040"]<-"Residuals"
    ac.an$F.Value[ac.an$ID == "11AEBGC.2040"]<-ac.an$Mean.Sq[ac.an$ID == "11AEBGC.2040"] / f.an$Mean.Sq[f.an$ID == "12R.0000"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "11AEBGC.2040"]<- stats::pf(ac.an$F.Value[ac.an$ID == "11AEBGC.2040"], ac.an$Df[ac.an$ID == "11AEBGC.2040"], f.an$Df[f.an$ID == "12R.0000"], lower.tail = F)

    # 	F = (MS_T(P(Aft)) xS(IC) / MS_T(P(Aft)) xS(CC)), or
    di11AEBGC.2030<-cbind(Source.of.Variation = "B:E:G:C",
                          sym.F[sym.F$Source.of.Variation == "A:E:B:G:C",2:3] -
                            sym.B[sym.B$Source.of.Variation == "B:E:G:C",2:3] -
                            sym.C[sym.C$Source.of.Variation == "A:E:B:G:C",2:3] +
                            sym.CB[sym.CB$Source.of.Variation == "B:E:G:C",2:3])
    ai.an<-dplyr::bind_rows(ai.an,di11AEBGC.2030)
    ai.an$ID[ai.an$Source.of.Variation == "B:E:G:C"] <- "11AEBGC.2030"
    ai.an$Df.Formula[ai.an$ID == "11AEBGC.2030"]<-"p(t-1)(s-1)"
    ai.an$Calculated.From[ai.an$ID == "11AEBGC.2030"]<-"sym.F-sym.B-sym.C+sym.BC"
    ai.an$Mean.Sq[ai.an$ID == "11AEBGC.2030"]<-ai.an$Sum.Sq[ai.an$ID == "11AEBGC.2030"] / ai.an$Df[ai.an$ID == "11AEBGC.2030"]
    ai.an$Source.of.Variation[ai.an$ID == "11AEBGC.2030"] <-"T(P(Aft)) x S(I)"
    ai.an$Denominator[ai.an$ID == "11AEBGC.2030"]<-"T(P(Aft))xS(C)"
    ai.an$F.Value[ai.an$ID == "11AEBGC.2030"]<-ai.an$Mean.Sq[ai.an$ID == "11AEBGC.2030"] / ac.an$Mean.Sq[ac.an$ID == "11AEBGC.2040"]
    ai.an$P.Value.Upper.Tail[ai.an$ID == "11AEBGC.2030"]<- stats::pf(ai.an$F.Value[ai.an$ID == "11AEBGC.2030"], ai.an$Df[ai.an$ID == "11AEBGC.2030"], ac.an$Df[ac.an$ID == "11AEBGC.2040"], lower.tail = F)

    # 	F = (MS_T(P(Aft)) xCC / MS_T(P(Aft)) xS(CC)), or
    di10AEBC.2040<-cbind(Source.of.Variation = "B:E:C",
                         sym.C[sym.C$Source.of.Variation == "A:E:B:C",2:3] - sym.CB[sym.CB$Source.of.Variation == "B:E:C",2:3])
    ac.an<-dplyr::bind_rows(ac.an,di10AEBC.2040)
    ac.an$ID[ac.an$Source.of.Variation == "B:E:C"] <- "10AEBC.2040"
    ac.an$Df.Formula[ac.an$ID == "10AEBC.2040"]<-"p(t-1)(c-1)"
    ac.an$Calculated.From[ac.an$ID == "10AEBC.2040"]<-"sym.C-sym.CB"
    ac.an$Mean.Sq[ac.an$ID == "10AEBC.2040"] <- ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] / ac.an$Df[ac.an$ID == "10AEBC.2040"]
    ac.an$Source.of.Variation[ac.an$ID == "10AEBC.2040"] <-"T(P(Aft)) x C"
    ac.an$Denominator[ac.an$ID == "10AEBC.2040"]<-"T(P(Aft))xS(C)"
    ac.an$F.Value[ac.an$ID == "10AEBC.2040"]<-ac.an$Mean.Sq[ac.an$ID == "10AEBC.2040"] / ac.an$Mean.Sq[ac.an$ID == "11AEBGC.2040"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "10AEBC.2040"]<- stats::pf(ac.an$F.Value[ac.an$ID == "10AEBC.2040"], ac.an$Df[ac.an$ID == "10AEBC.2040"], ac.an$Df[ac.an$ID == "11AEBGC.2040"], lower.tail = F)

    # 	F = (MS_T(P(Aft)) xIC / MS_T(P(Aft)) xCC), or
    di10AEBC.2030<-cbind(Source.of.Variation = "B:E:C",
                         sym.F[sym.F$Source.of.Variation == "A:E:B:C",2:3] -
                           sym.B[sym.B$Source.of.Variation == "B:E:C",2:3] -
                           sym.C[sym.C$Source.of.Variation == "A:E:B:C",2:3] +
                           sym.CB[sym.CB$Source.of.Variation == "B:E:C",2:3])
    ai.an<-dplyr::bind_rows(ai.an,di10AEBC.2030)
    ai.an$ID[ai.an$Source.of.Variation == "B:E:C"] <- "10AEBC.2030"
    ai.an$Df.Formula[ai.an$ID == "10AEBC.2030"]<-"p(t-1)"
    ai.an$Calculated.From[ai.an$ID == "10AEBC.2030"]<-"sym.F-sym.B-sym.C+sym.BC"
    ai.an$Mean.Sq[ai.an$ID == "10AEBC.2030"]<-ai.an$Sum.Sq[ai.an$ID == "10AEBC.2030"] / ai.an$Df[ai.an$ID == "10AEBC.2030"]
    ai.an$Source.of.Variation[ai.an$ID == "10AEBC.2030"] <-"T(P(Aft)) x I"
    ai.an$Denominator[ai.an$ID == "10AEBC.2030"]<-"T(P(Aft))xC"
    ai.an$F.Value[ai.an$ID == "10AEBC.2030"]<-ai.an$Mean.Sq[ai.an$ID == "10AEBC.2030"] / ac.an$Mean.Sq[ac.an$ID == "10AEBC.2040"]
    ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEBC.2030"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEBC.2030"], ai.an$Df[ai.an$ID == "10AEBC.2030"], ac.an$Df[ac.an$ID == "10AEBC.2040"], lower.tail = F)

    # 	F = (MS_P(Aft) xS(CC) / MS_T(P(Aft)) xS(CC)), or
    di09AEBG.2040<-cbind(Source.of.Variation = "B:E:G",
                         sym.C[sym.C$Source.of.Variation == "A:E:B:G",2:3] - sym.CB[sym.CB$Source.of.Variation == "B:E:G",2:3])
    ac.an<-dplyr::bind_rows(ac.an,di09AEBG.2040)
    ac.an$ID[ac.an$Source.of.Variation == "B:E:G"] <- "09AEBG.2040"
    ac.an$Df.Formula[ac.an$ID == "09AEBG.2040"]<-"(p-1)c(s-1)"
    ac.an$Calculated.From[ac.an$ID == "09AEBG.2040"]<-"sym.C-sym.CB"
    ac.an$Mean.Sq[ac.an$ID == "09AEBG.2040"] <- ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"] / ac.an$Df[ac.an$ID == "09AEBG.2040"]
    ac.an$Source.of.Variation[ac.an$ID == "09AEBG.2040"] <-"P(Aft) x S(C)"
    ac.an$Denominator[ac.an$ID == "09AEBG.2040"]<-"T(P(Aft))xS(C)"
    ac.an$F.Value[ac.an$ID == "09AEBG.2040"]<-ac.an$Mean.Sq[ac.an$ID == "09AEBG.2040"] / ac.an$Mean.Sq[ac.an$ID == "11AEBGC.2040"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "09AEBG.2040"]<- stats::pf(ac.an$F.Value[ac.an$ID == "09AEBG.2040"], ac.an$Df[ac.an$ID == "09AEBG.2040"], ac.an$Df[ac.an$ID == "11AEBGC.2040"], lower.tail = F)

    # 	F = (MS_P(Aft) xS(IC) / MS_P(Aft) xS(CC)), or
    di09AEBG.2030<-cbind(Source.of.Variation = "B:E:G",
                         sym.F[sym.F$Source.of.Variation == "A:E:B:G",2:3] -
                           sym.B[sym.B$Source.of.Variation == "B:E:G",2:3] -
                           sym.C[sym.C$Source.of.Variation == "A:E:B:G",2:3] +
                           sym.CB[sym.CB$Source.of.Variation == "B:E:G",2:3])
    ai.an<-dplyr::bind_rows(ai.an,di09AEBG.2030)
    ai.an$ID[ai.an$Source.of.Variation == "B:E:G"] <- "09AEBG.2030"
    ai.an$Df.Formula[ai.an$ID == "09AEBG.2030"]<-"(p-1)(s-1)"
    ai.an$Calculated.From[ai.an$ID == "09AEBG.2030"]<-"sym.F-sym.B-sym.C+sym.BC"
    ai.an$Mean.Sq[ai.an$ID == "09AEBG.2030"]<-ai.an$Sum.Sq[ai.an$ID == "09AEBG.2030"] / ai.an$Df[ai.an$ID == "09AEBG.2030"]
    ai.an$Source.of.Variation[ai.an$ID == "09AEBG.2030"] <-"P(Aft) x S(I)"
    ai.an$Denominator[ai.an$ID == "09AEBG.2030"]<-"P(Aft)xS(C)"
    ai.an$F.Value[ai.an$ID == "09AEBG.2030"]<-ai.an$Mean.Sq[ai.an$ID == "09AEBG.2030"] / ac.an$Mean.Sq[ac.an$ID == "09AEBG.2040"]
    ai.an$P.Value.Upper.Tail[ai.an$ID == "09AEBG.2030"]<- stats::pf(ai.an$F.Value[ai.an$ID == "09AEBG.2030"], ai.an$Df[ai.an$ID == "09AEBG.2030"], ac.an$Df[ac.an$ID == "09AEBG.2040"], lower.tail = F)

    # 	If varT(P(Aft)) xCC <= and varP(Aft) xS(CC) > 0, then F = (MS_P(Aft) xCC / MS_P(Aft) xS(CC))
    di08AEB.2040<-cbind(Source.of.Variation = "B:E",
                        sym.C[sym.C$Source.of.Variation == "A:E:B",2:3] - sym.CB[sym.CB$Source.of.Variation == "B:E",2:3])
    ac.an<-dplyr::bind_rows(ac.an,di08AEB.2040)
    ac.an$ID[ac.an$Source.of.Variation == "B:E"] <- "08AEB.2040"
    ac.an$Df.Formula[ac.an$ID == "08AEB.2040"]<-"(p-1)(c-1)"
    ac.an$Calculated.From[ac.an$ID == "08AEB.2040"]<-"sym.C-sym.CB"
    ac.an$Mean.Sq[ac.an$ID == "08AEB.2040"] <- ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] / ac.an$Df[ac.an$ID == "08AEB.2040"]
    ac.an$Source.of.Variation[ac.an$ID == "08AEB.2040"]<-"P(Aft) x C"
    ac.an$Denominator[ac.an$ID == "08AEB.2040"]<-"P(Aft)xS(C)"
    ac.an$F.Value[ac.an$ID == "08AEB.2040"]<-ac.an$Mean.Sq[ac.an$ID == "08AEB.2040"] / ac.an$Mean.Sq[ac.an$ID == "09AEBG.2040"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "08AEB.2040"]<- stats::pf(ac.an$F.Value[ac.an$ID == "08AEB.2040"], ac.an$Df[ac.an$ID == "08AEB.2040"], ac.an$Df[ac.an$ID == "09AEBG.2040"], lower.tail = F)

    # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xCC > 0, then F = (MS_P(Aft) xCC / MS_T(P(Aft)) xCC)
    ac.an<-dplyr::bind_rows(ac.an,di08AEB.2040)
    ac.an$ID[ac.an$Source.of.Variation == "B:E"] <- "08AEB.2041"
    ac.an$Df.Formula[ac.an$ID == "08AEB.2041"]<-"(p-1)(c-1)"
    ac.an$Calculated.From[ac.an$ID == "08AEB.2041"]<-"sym.C-sym.CB"
    ac.an$Mean.Sq[ac.an$ID == "08AEB.2041"] <- ac.an$Sum.Sq[ac.an$ID == "08AEB.2041"] / ac.an$Df[ac.an$ID == "08AEB.2041"]
    ac.an$Source.of.Variation[ac.an$ID == "08AEB.2041"]<-"P(Aft) x C"
    ac.an$Denominator[ac.an$ID == "08AEB.2041"]<-"T(P(Aft))xC"
    ac.an$F.Value[ac.an$ID == "08AEB.2041"]<-ac.an$Mean.Sq[ac.an$ID == "08AEB.2041"] / ac.an$Mean.Sq[ac.an$ID == "10AEBC.2040"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "08AEB.2041"]<- stats::pf(ac.an$F.Value[ac.an$ID == "08AEB.2041"], ac.an$Df[ac.an$ID == "08AEB.2041"], ac.an$Df[ac.an$ID == "10AEBC.2040"], lower.tail = F)

    # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xCC <= 0, then F = (MS_P(Aft) xCC / MS_pooled), where MS_pooled = (SS_P(Aft) xS(CC) + SS_T(P(Aft)) xCC + SS_T(P(Aft)) xS(CC)) / (df_P(Aft) xS(CC) + df_T(P(Aft)) xCC + df_T(P(Aft)) xS(CC))
    ac.an<-dplyr::bind_rows(ac.an,di08AEB.2040)
    ac.an$ID[ac.an$Source.of.Variation == "B:E"] <- "08AEB.2042"
    ac.an$Df.Formula[ac.an$ID == "08AEB.2042"]<-"(p-1)(c-1)"
    ac.an$Calculated.From[ac.an$ID == "08AEB.2042"]<-"sym.C-sym.CB"
    ac.an$Mean.Sq[ac.an$ID == "08AEB.2042"] <- ac.an$Sum.Sq[ac.an$ID == "08AEB.2042"] / ac.an$Df[ac.an$ID == "08AEB.2042"]
    ac.an$Source.of.Variation[ac.an$ID == "08AEB.2042"]<-"P(Aft) x C"
    ac.an$Denominator[ac.an$ID == "08AEB.2042"]<-"P(Aft)xS(C)+T(P(Aft)) xC+T(P(Aft))xS(C)"
    ac.an$Df.Pooled[ac.an$ID == "08AEB.2042"]<- ac.an$Df[ac.an$ID == "09AEBG.2040"] + ac.an$Df[ac.an$ID == "10AEBC.2040"] + ac.an$Df[ac.an$ID == "11AEBGC.2040"]
    ac.an$Sum.Sq.Pooled[ac.an$ID == "08AEB.2042"]<- ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"] + ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] + ac.an$Sum.Sq[ac.an$ID == "11AEBGC.2040"]
    ac.an$Mean.Sq.Pooled[ac.an$ID == "08AEB.2042"] <- ac.an$Sum.Sq.Pooled[ac.an$ID == "08AEB.2042"] / ac.an$Df.Pooled[ac.an$ID == "08AEB.2042"]
    ac.an$F.Value[ac.an$ID == "08AEB.2042"]<-ac.an$Mean.Sq[ac.an$ID == "08AEB.2042"] / ac.an$Mean.Sq.Pooled[ac.an$ID == "08AEB.2042"]
    ac.an$P.Value.Upper.Tail[ac.an$ID == "08AEB.2042"]<- stats::pf(ac.an$F.Value[ac.an$ID == "08AEB.2042"], ac.an$Df[ac.an$ID == "08AEB.2042"], ac.an$Df.Pooled[ac.an$ID == "08AEB.2042"], lower.tail = F)

    # 	F = (MS_P(Aft) xIC / MS_P(Aft) xCC), or
    di08AEB.2030<-cbind(Source.of.Variation = "B:E",
                        sym.F[sym.F$Source.of.Variation == "A:E:B",2:3] -
                          sym.B[sym.B$Source.of.Variation == "B:E",2:3] -
                          sym.C[sym.C$Source.of.Variation == "A:E:B",2:3] +
                          sym.CB[sym.CB$Source.of.Variation == "B:E",2:3])
    ai.an<-dplyr::bind_rows(ai.an,di08AEB.2030)
    ai.an$ID[ai.an$Source.of.Variation == "B:E"] <- "08AEB.2030"
    ai.an$Df.Formula[ai.an$ID == "08AEB.2030"]<-"(p-1)"
    ai.an$Calculated.From[ai.an$ID == "08AEB.2030"]<-"sym.F-sym.B-sym.C+sym.BC"
    ai.an$Mean.Sq[ai.an$ID == "08AEB.2030"]<-ai.an$Sum.Sq[ai.an$ID == "08AEB.2030"] / ai.an$Df[ai.an$ID == "08AEB.2030"]
    ai.an$Source.of.Variation[ai.an$ID == "08AEB.2030"] <-"P(Aft) x I"
    ai.an$Denominator[ai.an$ID == "08AEB.2030"]<-"P(Aft)xC"
    ai.an$F.Value[ai.an$ID == "08AEB.2030"]<-ai.an$Mean.Sq[ai.an$ID == "08AEB.2030"] / ac.an$Mean.Sq[ac.an$ID == "08AEB.2040"]
    ai.an$P.Value.Upper.Tail[ai.an$ID == "08AEB.2030"]<- stats::pf(ai.an$F.Value[ai.an$ID == "08AEB.2030"], ai.an$Df[ai.an$ID == "08AEB.2030"], ac.an$Df[ac.an$ID == "08AEB.2040"], lower.tail = F)

    #F = (MS_B xS(CC) / MS_P(Aft) xS(CC)), or
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E:G",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E:G"] <- "07AEG.0040"
    c.an$Df.Formula[c.an$ID == "07AEG.0040"]<-"(b-1)c(s-1)"
    c.an$Calculated.From[c.an$ID == "07AEG.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "07AEG.0040"] <- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] / c.an$Df[c.an$ID == "07AEG.0040"]
    c.an$Source.of.Variation[c.an$ID == "07AEG.0040"] <-"B x S(C)"
    c.an$Denominator[c.an$ID == "07AEG.0040"]<-"P(Aft)xS(C)"
    c.an$F.Value[c.an$ID == "07AEG.0040"]<-c.an$Mean.Sq[c.an$ID == "07AEG.0040"] / ac.an$Mean.Sq[ac.an$ID == "09AEBG.2040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "07AEG.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "07AEG.0040"], c.an$Df[c.an$ID == "07AEG.0040"], ac.an$Df[ac.an$ID == "09AEBG.2040"], lower.tail = F)

    # 	F = (MS_B xS(IC) / MS_B xS(CC)), or
    di07AEG.0030<-cbind(Source.of.Variation = "A:E:G",
                        sym.F[sym.F$Source.of.Variation == "A:E:G",2:3] -
                          sym.C[sym.C$Source.of.Variation == "A:E:G",2:3])
    i.an<-dplyr::bind_rows(i.an,di07AEG.0030)
    i.an$ID[i.an$Source.of.Variation == "A:E:G"] <- "07AEG.0030"
    i.an$Df.Formula[i.an$ID == "07AEG.0030"]<-"(b-1)(s-1)"
    i.an$Calculated.From[i.an$ID == "07AEG.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "07AEG.0030"]<-i.an$Sum.Sq[i.an$ID == "07AEG.0030"] / i.an$Df[i.an$ID == "07AEG.0030"]
    i.an$Source.of.Variation[i.an$ID == "07AEG.0030"] <-"B x S(I)"
    i.an$Denominator[i.an$ID == "07AEG.0030"]<-"BxS(C)"
    i.an$F.Value[i.an$ID == "07AEG.0030"]<-i.an$Mean.Sq[i.an$ID == "07AEG.0030"] / c.an$Mean.Sq[c.an$ID == "07AEG.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "07AEG.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "07AEG.0030"], i.an$Df[i.an$ID == "07AEG.0030"], c.an$Df[c.an$ID == "07AEG.0040"], lower.tail = F)

    # 	If varP(Aft) xCC <= 0 and varB xS(CC) > 0 then F = (MS_B xCC / MS_B xS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0040"
    c.an$Df.Formula[c.an$ID == "03AE.0040"]<-"(b-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03AE.0040"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03AE.0040"] <- c.an$Sum.Sq[c.an$ID == "03AE.0040"] / c.an$Df[c.an$ID == "03AE.0040"]
    c.an$Source.of.Variation[c.an$ID == "03AE.0040"]<-"B x C"
    c.an$Denominator[c.an$ID == "03AE.0040"]<-"BxS(C)"
    c.an$F.Value[c.an$ID == "03AE.0040"]<-c.an$Mean.Sq[c.an$ID == "03AE.0040"] / c.an$Mean.Sq[c.an$ID == "07AEG.0040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0040"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0040"], c.an$Df[c.an$ID == "03AE.0040"], c.an$Df[c.an$ID == "07AEG.0040"], lower.tail = F)

    # 	If varB xS(CC) <= 0 and varP(Aft) xCC > 0, then F = (MS_B xCC / MS_P(Aft) xCC)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0041"
    c.an$Df.Formula[c.an$ID == "03AE.0041"]<-"(b-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03AE.0041"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03AE.0041"] <- c.an$Sum.Sq[c.an$ID == "03AE.0041"] / c.an$Df[c.an$ID == "03AE.0041"]
    c.an$Source.of.Variation[c.an$ID == "03AE.0041"]<-"B x C"
    c.an$Denominator[c.an$ID == "03AE.0041"]<-"P(Aft)xC"
    c.an$F.Value[c.an$ID == "03AE.0041"]<-c.an$Mean.Sq[c.an$ID == "03AE.0041"] / ac.an$Mean.Sq[ac.an$ID == "08AEB.2040"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0041"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0041"], c.an$Df[c.an$ID == "03AE.0041"], ac.an$Df[ac.an$ID == "08AEB.2040"], lower.tail = F)

    # 	If varB xS(CC) <= 0, varP(Aft) xCC <= 0, varT(P(Aft)) xCC <= 0 and varP(Aft) xS(CC) > 0,
    # then F = (MS_B xCC / MS_pooled), where MS_pooled = (SS_B xS(CC) + SS_P(Aft) xCC + SS_P(Aft) xS(CC)) / (df_B xS(CC) + df_P(Aft) xCC + df_P(Aft) xS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0042"
    c.an$Df.Formula[c.an$ID == "03AE.0042"]<-"(b-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03AE.0042"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03AE.0042"] <- c.an$Sum.Sq[c.an$ID == "03AE.0042"] / c.an$Df[c.an$ID == "03AE.0042"]
    c.an$Source.of.Variation[c.an$ID == "03AE.0042"]<-"B x C"
    c.an$Denominator[c.an$ID == "03AE.0042"]<-"BxS(C)+P(Aft)xC+P(Aft)xS(C)"
    c.an$Df.Pooled[c.an$ID == "03AE.0042"]<- c.an$Df[c.an$ID == "07AEG.0040"] + ac.an$Df[ac.an$ID == "08AEB.2040"] + ac.an$Df[ac.an$ID == "09AEBG.2040"]
    c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0042"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"]
    c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0042"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0042"] / c.an$Df.Pooled[c.an$ID == "03AE.0042"]
    c.an$F.Value[c.an$ID == "03AE.0042"]<-c.an$Mean.Sq[c.an$ID == "03AE.0042"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0042"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0042"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0042"], c.an$Df[c.an$ID == "03AE.0042"], c.an$Df.Pooled[c.an$ID == "03AE.0042"], lower.tail = F)

    # 	If varB xS(CC) <= 0, varP(Aft) xCC <= 0, varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xCC > 0,
    #then F = (MSB xCC / MSpooled), where MSpooled = (SS_B xS(CC) + SS_P(Aft) xCC + SS_T(P(Aft)) xCC) / (df_B xS(CC) + df_P(Aft) xCC + df_T(P(Aft)) xCC)
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0043"
    c.an$Df.Formula[c.an$ID == "03AE.0043"]<-"(b-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03AE.0043"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03AE.0043"] <- c.an$Sum.Sq[c.an$ID == "03AE.0043"] / c.an$Df[c.an$ID == "03AE.0043"]
    c.an$Source.of.Variation[c.an$ID == "03AE.0043"]<-"B x C"
    c.an$Denominator[c.an$ID == "03AE.0043"]<-"BxS(C)+P(Aft)xC+T(P(Aft))xC"
    c.an$Df.Pooled[c.an$ID == "03AE.0043"]<- c.an$Df[c.an$ID == "07AEG.0040"] + ac.an$Df[ac.an$ID == "08AEB.2040"] + ac.an$Df[ac.an$ID == "10AEBC.2040"]
    c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0043"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"]
    c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0043"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0043"] / c.an$Df.Pooled[c.an$ID == "03AE.0043"]
    c.an$F.Value[c.an$ID == "03AE.0043"]<-c.an$Mean.Sq[c.an$ID == "03AE.0043"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0043"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0043"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0043"], c.an$Df[c.an$ID == "03AE.0043"], c.an$Df.Pooled[c.an$ID == "03AE.0043"], lower.tail = F)

    # 	If varB xS(CC) <= 0, varP(Aft) xCC <= 0, varP(Aft) xS(CC) <= 0, varT(P(Aft)) xCC <= 0,
    #then F = (MSB xCC / MSpooled), where MSpooled = (SS_B xS(CC) + S_SP(Aft) xCC + SS_P(Aft) xS(CC) + SS_T(P(Aft)) xCC + SS_T(P(Aft)) xS(CC)) / (df_B xS(CC) + df_P(Aft) xCC + df_P(Aft) xS(CC) + df_T(P(Aft)) xCC + df_T(P(Aft)) xS(CC))
    c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
    c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0044"
    c.an$Df.Formula[c.an$ID == "03AE.0044"]<-"(b-1)(c-1)"
    c.an$Calculated.From[c.an$ID == "03AE.0044"]<-"sym.C"
    c.an$Mean.Sq[c.an$ID == "03AE.0044"] <- c.an$Sum.Sq[c.an$ID == "03AE.0044"] / c.an$Df[c.an$ID == "03AE.0044"]
    c.an$Source.of.Variation[c.an$ID == "03AE.0044"]<-"B x C"
    c.an$Denominator[c.an$ID == "03AE.0044"]<-"BxS(C)+P(Aft)xC+P(Aft)xS(C)+T(P(Aft))xC+T(P(Aft))xS(C)"
    c.an$Df.Pooled[c.an$ID == "03AE.0044"]<- c.an$Df[c.an$ID == "07AEG.0040"] + ac.an$Df[ac.an$ID == "08AEB.2040"] + ac.an$Df[ac.an$ID == "09AEBG.2040"] + ac.an$Df[ac.an$ID == "10AEBC.2040"] + ac.an$Df[ac.an$ID == "11AEBGC.2040"]
    c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0044"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"] + ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] + ac.an$Sum.Sq[ac.an$ID == "11AEBGC.2040"]
    c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0044"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0044"] / c.an$Df.Pooled[c.an$ID == "03AE.0044"]
    c.an$F.Value[c.an$ID == "03AE.0044"]<-c.an$Mean.Sq[c.an$ID == "03AE.0044"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0044"]
    c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0044"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0044"], c.an$Df[c.an$ID == "03AE.0044"], c.an$Df.Pooled[c.an$ID == "03AE.0044"], lower.tail = F)

    # 	F = (MS_B xIC / MS_B xCC), or
    di03AE.0030<-cbind(Source.of.Variation = "A:E",
                       sym.F[sym.F$Source.of.Variation == "A:E",2:3] -
                         sym.C[sym.C$Source.of.Variation == "A:E",2:3])
    i.an<-dplyr::bind_rows(i.an,di03AE.0030)
    i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0030"
    i.an$Df.Formula[i.an$ID == "03AE.0030"]<-"(b-1)"
    i.an$Calculated.From[i.an$ID == "03AE.0030"]<-"sym.F-sym.C"
    i.an$Mean.Sq[i.an$ID == "03AE.0030"]<-i.an$Sum.Sq[i.an$ID == "03AE.0030"] / i.an$Df[i.an$ID == "03AE.0030"]
    i.an$Source.of.Variation[i.an$ID == "03AE.0030"] <-"B x I"
    i.an$Denominator[i.an$ID == "03AE.0030"]<-"BxC"
    i.an$F.Value[i.an$ID == "03AE.0030"]<-i.an$Mean.Sq[i.an$ID == "03AE.0030"] / c.an$Mean.Sq[c.an$ID == "03AE.0040"]
    i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0030"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0030"], i.an$Df[i.an$ID == "03AE.0030"], c.an$Df[c.an$ID == "03AE.0040"], lower.tail = F)

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	If varT(P(Aft)) xS(CC) <= 0, then F = (MS_T(P(Aft)) xS(IC) / MS_pooled), where MS_pooled = (SS_T(P(Aft)) xS(CC) + SS_Res) / (df_T(P(Aft)) xS(CC) + df_Res)
      ai.an<-dplyr::bind_rows(ai.an,di11AEBGC.2030)
      ai.an$ID[ai.an$Source.of.Variation == "B:E:G:C"] <- "11AEBGC.2031"
      ai.an$Df.Formula[ai.an$ID == "11AEBGC.2031"]<-"p(t-1)(s-1)"
      ai.an$Calculated.From[ai.an$ID == "11AEBGC.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
      ai.an$Mean.Sq[ai.an$ID == "11AEBGC.2031"]<-ai.an$Sum.Sq[ai.an$ID == "11AEBGC.2031"] / ai.an$Df[ai.an$ID == "11AEBGC.2031"]
      ai.an$Source.of.Variation[ai.an$ID == "11AEBGC.2031"] <-"T(P(Aft)) x S(I)"
      ai.an$Denominator[ai.an$ID == "11AEBGC.2031"]<-"T(P(Aft))xS(C)+Residuals"
      ai.an$Df.Pooled[ai.an$ID == "11AEBGC.2031"]<- ac.an$Df[ac.an$ID == "11AEBGC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
      ai.an$Sum.Sq.Pooled[ai.an$ID == "11AEBGC.2031"]<- ac.an$Sum.Sq[ac.an$ID == "11AEBGC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
      ai.an$Mean.Sq.Pooled[ai.an$ID == "11AEBGC.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "11AEBGC.2031"] / ai.an$Df.Pooled[ai.an$ID == "11AEBGC.2031"]
      ai.an$F.Value[ai.an$ID == "11AEBGC.2031"]<-ai.an$Mean.Sq[ai.an$ID == "11AEBGC.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "11AEBGC.2031"]
      ai.an$P.Value.Upper.Tail[ai.an$ID == "11AEBGC.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "11AEBGC.2031"], ai.an$Df[ai.an$ID == "11AEBGC.2031"], ai.an$Df.Pooled[ai.an$ID == "11AEBGC.2031"], lower.tail = F)

      if(method == "eliminate") {
        # 	F = (MS_T(P(Aft)) xCC / MS_pooled), where MS_pooled = (SS_T(P(Aft)) xS(CC)) + SS_Res) / (df_T(P(Aft)) xS(CC)) + df_Res)
        ac.an<-dplyr::bind_rows(ac.an,di10AEBC.2040)
        ac.an$ID[ac.an$Source.of.Variation == "B:E:C"] <- "10AEBC.2041"
        ac.an$Df.Formula[ac.an$ID == "10AEBC.2041"]<-"p(t-1)(c-1)"
        ac.an$Calculated.From[ac.an$ID == "10AEBC.2041"]<-"sym.C-sym.CB"
        ac.an$Mean.Sq[ac.an$ID == "10AEBC.2041"] <- ac.an$Sum.Sq[ac.an$ID == "10AEBC.2041"] / ac.an$Df[ac.an$ID == "10AEBC.2041"]
        ac.an$Source.of.Variation[ac.an$ID == "10AEBC.2041"] <-"T(P(Aft)) x C"
        ac.an$Denominator[ac.an$ID == "10AEBC.2041"]<-"T(P(Aft))xS(C)+Residuals"
        ac.an$Df.Pooled[ac.an$ID == "10AEBC.2041"]<- ac.an$Df[ac.an$ID == "11AEBGC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ac.an$Sum.Sq.Pooled[ac.an$ID == "10AEBC.2041"]<- ac.an$Sum.Sq[ac.an$ID == "11AEBGC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ac.an$Mean.Sq.Pooled[ac.an$ID == "10AEBC.2041"] <- ac.an$Sum.Sq.Pooled[ac.an$ID == "10AEBC.2041"] / ac.an$Df.Pooled[ac.an$ID == "10AEBC.2041"]
        ac.an$F.Value[ac.an$ID == "10AEBC.2041"]<-ac.an$Mean.Sq[ac.an$ID == "10AEBC.2041"] / ac.an$Mean.Sq.Pooled[ac.an$ID == "10AEBC.2041"]
        ac.an$P.Value.Upper.Tail[ac.an$ID == "10AEBC.2041"]<- stats::pf(ac.an$F.Value[ac.an$ID == "10AEBC.2041"], ac.an$Df[ac.an$ID == "10AEBC.2041"], ac.an$Df.Pooled[ac.an$ID == "10AEBC.2041"], lower.tail = F)

        # 	If varT(P(Aft)) xCC <= 0, then F = (MS_T(P(Aft)) xIC / MS_pooled), where MS_pooled = (SS_T(P(Aft)) xCC + SS_Res) / (df_T(P(Aft)) xCC + df_Res)
        ai.an<-dplyr::bind_rows(ai.an,di10AEBC.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E:C"] <- "10AEBC.2031"
        ai.an$Df.Formula[ai.an$ID == "10AEBC.2031"]<-"p(t-1)"
        ai.an$Calculated.From[ai.an$ID == "10AEBC.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "10AEBC.2031"]<-ai.an$Sum.Sq[ai.an$ID == "10AEBC.2031"] / ai.an$Df[ai.an$ID == "10AEBC.2031"]
        ai.an$Source.of.Variation[ai.an$ID == "10AEBC.2031"] <-"T(P(Aft)) x I"
        ai.an$Denominator[ai.an$ID == "10AEBC.2031"]<-"T(P(Aft))xC+Residuals"
        ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"]<- ac.an$Df[ac.an$ID == "10AEBC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEBC.2031"]<- ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEBC.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEBC.2031"] / ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"]
        ai.an$F.Value[ai.an$ID == "10AEBC.2031"]<-ai.an$Mean.Sq[ai.an$ID == "10AEBC.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEBC.2031"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEBC.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEBC.2031"], ai.an$Df[ai.an$ID == "10AEBC.2031"], ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"], lower.tail = F)

        # 	F = (MS_P(Aft) xS(CC) / MS_pooled), where MS_pooled = (SS_T(P(Aft)) xS(CC) + SS_Res) / (df_T(P(Aft)) xS(CC) + df_Res)
        ac.an<-dplyr::bind_rows(ac.an,di09AEBG.2040)
        ac.an$ID[ac.an$Source.of.Variation == "B:E:G"] <- "09AEBG.2041"
        ac.an$Df.Formula[ac.an$ID == "09AEBG.2041"]<-"(p-1)c(s-1)"
        ac.an$Calculated.From[ac.an$ID == "09AEBG.2041"]<-"sym.C-sym.CB"
        ac.an$Mean.Sq[ac.an$ID == "09AEBG.2041"] <- ac.an$Sum.Sq[ac.an$ID == "09AEBG.2041"] / ac.an$Df[ac.an$ID == "09AEBG.2041"]
        ac.an$Source.of.Variation[ac.an$ID == "09AEBG.2041"] <-"P(Aft) x S(C)"
        ac.an$Denominator[ac.an$ID == "09AEBG.2041"]<-"T(P(Aft))xS(C)+Residuals"
        ac.an$Df.Pooled[ac.an$ID == "09AEBG.2041"]<- ac.an$Df[ac.an$ID == "11AEBGC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ac.an$Sum.Sq.Pooled[ac.an$ID == "09AEBG.2041"]<- ac.an$Sum.Sq[ac.an$ID == "11AEBGC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ac.an$Mean.Sq.Pooled[ac.an$ID == "09AEBG.2041"] <- ac.an$Sum.Sq.Pooled[ac.an$ID == "09AEBG.2041"] / ac.an$Df.Pooled[ac.an$ID == "09AEBG.2041"]
        ac.an$F.Value[ac.an$ID == "09AEBG.2041"]<-ac.an$Mean.Sq[ac.an$ID == "09AEBG.2041"] / ac.an$Mean.Sq.Pooled[ac.an$ID == "09AEBG.2041"]
        ac.an$P.Value.Upper.Tail[ac.an$ID == "09AEBG.2041"]<- stats::pf(ac.an$F.Value[ac.an$ID == "09AEBG.2041"], ac.an$Df[ac.an$ID == "09AEBG.2041"], ac.an$Df.Pooled[ac.an$ID == "09AEBG.2041"], lower.tail = F)

        # 	If varP(Aft) xS(CC) <= 0, then F = (MS_P(Aft) xS(IC) / MSpooled), where MSpooled = (SS_P(Aft) xS(CC) + SS_Res) / (df_P(Aft) xS(CC) + df_Res)
        ai.an<-dplyr::bind_rows(ai.an,di09AEBG.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E:G"] <- "09AEBG.2031"
        ai.an$Df.Formula[ai.an$ID == "09AEBG.2031"]<-"(p-1)(s-1)"
        ai.an$Calculated.From[ai.an$ID == "09AEBG.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "09AEBG.2031"]<-ai.an$Sum.Sq[ai.an$ID == "09AEBG.2031"] / ai.an$Df[ai.an$ID == "09AEBG.2031"]
        ai.an$Source.of.Variation[ai.an$ID == "09AEBG.2031"] <-"P(Aft) x S(I)"
        ai.an$Denominator[ai.an$ID == "09AEBG.2031"]<-"P(Aft)xS(C)+Residuals"
        ai.an$Df.Pooled[ai.an$ID == "09AEBG.2031"]<- ac.an$Df[ac.an$ID == "09AEBG.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "09AEBG.2031"]<- ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "09AEBG.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "09AEBG.2031"] / ai.an$Df.Pooled[ai.an$ID == "09AEBG.2031"]
        ai.an$F.Value[ai.an$ID == "09AEBG.2031"]<-ai.an$Mean.Sq[ai.an$ID == "09AEBG.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "09AEBG.2031"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "09AEBG.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "09AEBG.2031"], ai.an$Df[ai.an$ID == "09AEBG.2031"], ai.an$Df.Pooled[ai.an$ID == "09AEBG.2031"], lower.tail = F)

        # 	If varT(P(Aft)) xCC <= and varP(Aft) xS(CC) > 0, then F = (MS_P(Aft) xCC / MSpooled), where MSpooled = (SS_P(Aft) xS(CC) + SS_Res) / (df_P(Aft) xS(CC) + df_Res)
        ac.an<-dplyr::bind_rows(ac.an,di08AEB.2040)
        ac.an$ID[ac.an$Source.of.Variation == "B:E"] <- "08AEB.2043"
        ac.an$Df.Formula[ac.an$ID == "08AEB.2043"]<-"(p-1)(c-1)"
        ac.an$Calculated.From[ac.an$ID == "08AEB.2043"]<-"sym.C-sym.CB"
        ac.an$Mean.Sq[ac.an$ID == "08AEB.2043"] <- ac.an$Sum.Sq[ac.an$ID == "08AEB.2043"] / ac.an$Df[ac.an$ID == "08AEB.2043"]
        ac.an$Source.of.Variation[ac.an$ID == "08AEB.2043"]<-"P(Aft) x C"
        ac.an$Denominator[ac.an$ID == "08AEB.2043"]<-"P(Aft)xS(C)+Residuals"
        ac.an$Df.Pooled[ac.an$ID == "08AEB.2043"]<- ac.an$Df[ac.an$ID == "09AEBG.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ac.an$Sum.Sq.Pooled[ac.an$ID == "08AEB.2043"]<- ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ac.an$Mean.Sq.Pooled[ac.an$ID == "08AEB.2043"] <- ac.an$Sum.Sq.Pooled[ac.an$ID == "08AEB.2043"] / ac.an$Df.Pooled[ac.an$ID == "08AEB.2043"]
        ac.an$F.Value[ac.an$ID == "08AEB.2043"]<-ac.an$Mean.Sq[ac.an$ID == "08AEB.2043"] / ac.an$Mean.Sq.Pooled[ac.an$ID == "08AEB.2043"]
        ac.an$P.Value.Upper.Tail[ac.an$ID == "08AEB.2043"]<- stats::pf(ac.an$F.Value[ac.an$ID == "08AEB.2043"], ac.an$Df[ac.an$ID == "08AEB.2043"], ac.an$Df.Pooled[ac.an$ID == "08AEB.2043"], lower.tail = F)

        # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xCC > 0, then F = (MS_P(Aft) xCC / MS_pooled), where MSpooled = (SS_T(P(Aft)) xCC + SS_Res) / (df_T(P(Aft)) xCC + df_Res)
        ac.an<-dplyr::bind_rows(ac.an,di08AEB.2040)
        ac.an$ID[ac.an$Source.of.Variation == "B:E"] <- "08AEB.2044"
        ac.an$Df.Formula[ac.an$ID == "08AEB.2044"]<-"(p-1)(c-1)"
        ac.an$Calculated.From[ac.an$ID == "08AEB.2044"]<-"sym.C-sym.CB"
        ac.an$Mean.Sq[ac.an$ID == "08AEB.2044"] <- ac.an$Sum.Sq[ac.an$ID == "08AEB.2044"] / ac.an$Df[ac.an$ID == "08AEB.2044"]
        ac.an$Source.of.Variation[ac.an$ID == "08AEB.2044"]<-"P(Aft) x C"
        ac.an$Denominator[ac.an$ID == "08AEB.2044"]<-"T(P(Aft))xC+Residuals"
        ac.an$Df.Pooled[ac.an$ID == "08AEB.2044"]<- ac.an$Df[ac.an$ID == "10AEBC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ac.an$Sum.Sq.Pooled[ac.an$ID == "08AEB.2044"]<- ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ac.an$Mean.Sq.Pooled[ac.an$ID == "08AEB.2044"] <- ac.an$Sum.Sq.Pooled[ac.an$ID == "08AEB.2044"] / ac.an$Df.Pooled[ac.an$ID == "08AEB.2044"]
        ac.an$F.Value[ac.an$ID == "08AEB.2044"]<-ac.an$Mean.Sq[ac.an$ID == "08AEB.2044"] / ac.an$Mean.Sq.Pooled[ac.an$ID == "08AEB.2044"]
        ac.an$P.Value.Upper.Tail[ac.an$ID == "08AEB.2044"]<- stats::pf(ac.an$F.Value[ac.an$ID == "08AEB.2044"], ac.an$Df[ac.an$ID == "08AEB.2044"], ac.an$Df.Pooled[ac.an$ID == "08AEB.2044"], lower.tail = F)

        # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xCC <= 0, then F = (MSP(Aft) xCC / MSpooled), where MSpooled = (SS_P(Aft) xS(CC) + SS_T(P(Aft)) xCC + SS_T(P(Aft)) xS(CC) + SSRes) / (dfP(Aft) xS(CC) + dfT(P(Aft)) xCC + dfT(P(Aft)) xS(CC) + dfRes)
        ac.an<-dplyr::bind_rows(ac.an,di08AEB.2040)
        ac.an$ID[ac.an$Source.of.Variation == "B:E"] <- "08AEB.2045"
        ac.an$Df.Formula[ac.an$ID == "08AEB.2045"]<-"(p-1)(c-1)"
        ac.an$Calculated.From[ac.an$ID == "08AEB.2045"]<-"sym.C-sym.CB"
        ac.an$Mean.Sq[ac.an$ID == "08AEB.2045"] <- ac.an$Sum.Sq[ac.an$ID == "08AEB.2045"] / ac.an$Df[ac.an$ID == "08AEB.2045"]
        ac.an$Source.of.Variation[ac.an$ID == "08AEB.2045"]<-"P(Aft) x C"
        ac.an$Denominator[ac.an$ID == "08AEB.2045"]<-"P(Aft)xS(C)+T(P(Aft))xC+T(P(Aft))xS(C)+Residuals"
        ac.an$Df.Pooled[ac.an$ID == "08AEB.2045"]<- ac.an$Df[ac.an$ID == "09AEBG.2040"] + ac.an$Df[ac.an$ID == "10AEBC.2040"] + ac.an$Df[ac.an$ID == "11AEBGC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ac.an$Sum.Sq.Pooled[ac.an$ID == "08AEB.2045"]<- ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"] + ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] + ac.an$Sum.Sq[ac.an$ID == "11AEBGC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ac.an$Mean.Sq.Pooled[ac.an$ID == "08AEB.2045"] <- ac.an$Sum.Sq.Pooled[ac.an$ID == "08AEB.2045"] / ac.an$Df.Pooled[ac.an$ID == "08AEB.2045"]
        ac.an$F.Value[ac.an$ID == "08AEB.2045"]<-ac.an$Mean.Sq[ac.an$ID == "08AEB.2045"] / ac.an$Mean.Sq.Pooled[ac.an$ID == "08AEB.2045"]
        ac.an$P.Value.Upper.Tail[ac.an$ID == "08AEB.2045"]<- stats::pf(ac.an$F.Value[ac.an$ID == "08AEB.2045"], ac.an$Df[ac.an$ID == "08AEB.2045"], ac.an$Df.Pooled[ac.an$ID == "08AEB.2045"], lower.tail = F)

        # 	If varP(Aft) xCC <= 0, then F = (MS_P(Aft) xIC / MS_pooled), where MS_pooled = (SS_P(Aft) xCC + SS_Res) / (df_P(Aft) xCC + df_Res)
        ai.an<-dplyr::bind_rows(ai.an,di08AEB.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E"] <- "08AEB.2031"
        ai.an$Df.Formula[ai.an$ID == "08AEB.2031"]<-"(p-1)"
        ai.an$Calculated.From[ai.an$ID == "08AEB.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "08AEB.2031"]<-ai.an$Sum.Sq[ai.an$ID == "08AEB.2031"] / ai.an$Df[ai.an$ID == "08AEB.2031"]
        ai.an$Source.of.Variation[ai.an$ID == "08AEB.2031"] <-"P(Aft) x I"
        ai.an$Denominator[ai.an$ID == "08AEB.2031"]<-"P(Aft)xC+Residuals"
        ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"]<- ac.an$Df[ac.an$ID == "08AEB.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2031"]<- ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2031"] / ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"]
        ai.an$F.Value[ai.an$ID == "08AEB.2031"]<-ai.an$Mean.Sq[ai.an$ID == "08AEB.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2031"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "08AEB.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "08AEB.2031"], ai.an$Df[ai.an$ID == "08AEB.2031"], ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"], lower.tail = F)

        # 	F = (MS_B xS(CC) / MS_pooled), where MS_pooled = (SS_P(Aft) xS(CC) + SS_Res) / (df_P(Aft) xS(CC) + df_Res)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E:G",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E:G"] <- "07AEG.0041"
        c.an$Df.Formula[c.an$ID == "07AEG.0041"]<-"(b-1)c(s-1)"
        c.an$Calculated.From[c.an$ID == "07AEG.0041"]<-"sym.C-sym.CB"
        c.an$Mean.Sq[c.an$ID == "07AEG.0041"] <- c.an$Sum.Sq[c.an$ID == "07AEG.0041"] / c.an$Df[c.an$ID == "07AEG.0041"]
        c.an$Source.of.Variation[c.an$ID == "07AEG.0041"] <-"B x S(C)"
        c.an$Denominator[c.an$ID == "07AEG.0041"]<-"P(Aft)xS(C)+Residuals"
        c.an$Df.Pooled[c.an$ID == "07AEG.0041"]<- ac.an$Df[ac.an$ID == "09AEBG.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "07AEG.0041"]<- ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "07AEG.0041"] <- c.an$Sum.Sq.Pooled[c.an$ID == "07AEG.0041"] / c.an$Df.Pooled[c.an$ID == "07AEG.0041"]
        c.an$F.Value[c.an$ID == "07AEG.0041"]<-c.an$Mean.Sq[c.an$ID == "07AEG.0041"] / c.an$Mean.Sq.Pooled[c.an$ID == "07AEG.0041"]
        c.an$P.Value.Upper.Tail[c.an$ID == "07AEG.0041"]<- stats::pf(c.an$F.Value[c.an$ID == "07AEG.0041"], c.an$Df[c.an$ID == "07AEG.0041"], c.an$Df.Pooled[c.an$ID == "07AEG.0041"], lower.tail = F)

        #If varB xS(CC) <= 0, then F = (MS_B xS(IC) / MS_pooled), where MSpooled = (SS_B xS(CC) + SS_Res) / (df_B xS(CC) + df_Res)
        i.an<-dplyr::bind_rows(i.an,di07AEG.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E:G"] <- "07AEG.0031"
        i.an$Df.Formula[i.an$ID == "07AEG.0031"]<-"(b-1)(s-1)"
        i.an$Calculated.From[i.an$ID == "07AEG.0031"]<-"sym.F-sym.B-sym.C+sym.BC"
        i.an$Mean.Sq[i.an$ID == "07AEG.0031"]<-i.an$Sum.Sq[i.an$ID == "07AEG.0031"] / i.an$Df[i.an$ID == "07AEG.0031"]
        i.an$Source.of.Variation[i.an$ID == "07AEG.0031"] <-"B x S(I)"
        i.an$Denominator[i.an$ID == "07AEG.0031"]<-"BxS(C)+Residuals"
        i.an$Df.Pooled[i.an$ID == "07AEG.0031"]<- c.an$Df[c.an$ID == "07AEG.0040"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0031"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0031"] / i.an$Df.Pooled[i.an$ID == "07AEG.0031"]
        i.an$F.Value[i.an$ID == "07AEG.0031"]<-i.an$Mean.Sq[i.an$ID == "07AEG.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "07AEG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "07AEG.0031"], i.an$Df[i.an$ID == "07AEG.0031"], i.an$Df.Pooled[i.an$ID == "07AEG.0031"], lower.tail = F)

        # 	If varP(Aft) xCC <= 0 and varB xS(CC) > 0 then F = (MS_B xCC / MSpooled), where MSpooled = (SS_B xS(CC) + SS_Res) / (df_B xS(CC) + df_Res)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0045"
        c.an$Df.Formula[c.an$ID == "03AE.0045"]<-"(b-1)(c-1)"
        c.an$Calculated.From[c.an$ID == "03AE.0045"]<-"sym.C"
        c.an$Mean.Sq[c.an$ID == "03AE.0045"] <- c.an$Sum.Sq[c.an$ID == "03AE.0045"] / c.an$Df[c.an$ID == "03AE.0045"]
        c.an$Source.of.Variation[c.an$ID == "03AE.0045"]<-"B x C"
        c.an$Denominator[c.an$ID == "03AE.0045"]<-"BxS(C)+Residuals"
        c.an$Df.Pooled[c.an$ID == "03AE.0045"]<- c.an$Df[c.an$ID == "07AEG.0040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0045"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0045"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0045"] / c.an$Df.Pooled[c.an$ID == "03AE.0045"]
        c.an$F.Value[c.an$ID == "03AE.0045"]<-c.an$Mean.Sq[c.an$ID == "03AE.0045"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0045"]
        c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0045"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0045"], c.an$Df[c.an$ID == "03AE.0045"], c.an$Df.Pooled[c.an$ID == "03AE.0045"], lower.tail = F)

        # 	If varB xS(CC) <= 0 and varP(Aft) xCC > 0, then F = (MS_B xCC / MSpooled), where MSpooled = (SS_P(Aft) xCC + SS_Res) / (df_P(Aft) xCC + df_Res)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0046"
        c.an$Df.Formula[c.an$ID == "03AE.0046"]<-"(b-1)(c-1)"
        c.an$Calculated.From[c.an$ID == "03AE.0046"]<-"sym.C"
        c.an$Mean.Sq[c.an$ID == "03AE.0046"] <- c.an$Sum.Sq[c.an$ID == "03AE.0046"] / c.an$Df[c.an$ID == "03AE.0046"]
        c.an$Source.of.Variation[c.an$ID == "03AE.0046"]<-"B x C"
        c.an$Denominator[c.an$ID == "03AE.0046"]<-"P(Aft)xC+Residuals"
        c.an$Df.Pooled[c.an$ID == "03AE.0046"]<- ac.an$Df[ac.an$ID == "08AEB.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0046"]<- ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0046"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0046"] / c.an$Df.Pooled[c.an$ID == "03AE.0046"]
        c.an$F.Value[c.an$ID == "03AE.0046"]<-c.an$Mean.Sq[c.an$ID == "03AE.0046"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0046"]
        c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0046"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0046"], c.an$Df[c.an$ID == "03AE.0046"], c.an$Df.Pooled[c.an$ID == "03AE.0046"], lower.tail = F)

        # 	varB xS(CC) <= 0, varP(Aft) xCC <= 0, varT(P(Aft)) xCC <= 0 and varP(Aft) xS(CC) > 0, then F = (MS_B xCC / MSpooled),
        # where MSpooled = (SS_B xS(CC) + SS_P(Aft) xCC + SS_P(Aft) xS(CC) + SS_Res) / (df_B xS(CC) + df_P(Aft) xCC + df_P(Aft) xS(CC) + df_Res)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0047"
        c.an$Df.Formula[c.an$ID == "03AE.0047"]<-"(b-1)(c-1)"
        c.an$Calculated.From[c.an$ID == "03AE.0047"]<-"sym.C"
        c.an$Mean.Sq[c.an$ID == "03AE.0047"] <- c.an$Sum.Sq[c.an$ID == "03AE.0047"] / c.an$Df[c.an$ID == "03AE.0047"]
        c.an$Source.of.Variation[c.an$ID == "03AE.0047"]<-"B x C"
        c.an$Denominator[c.an$ID == "03AE.0047"]<-"BxS(C)+P(Aft)xC+P(Aft)xS(C)+Residuals"
        c.an$Df.Pooled[c.an$ID == "03AE.0047"]<- c.an$Df[c.an$ID == "07AEG.0040"] + ac.an$Df[ac.an$ID == "08AEB.2040"] + ac.an$Df[ac.an$ID == "09AEBG.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0047"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0047"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0047"] / c.an$Df.Pooled[c.an$ID == "03AE.0047"]
        c.an$F.Value[c.an$ID == "03AE.0047"]<-c.an$Mean.Sq[c.an$ID == "03AE.0047"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0047"]
        c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0047"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0047"], c.an$Df[c.an$ID == "03AE.0047"], c.an$Df.Pooled[c.an$ID == "03AE.0047"], lower.tail = F)

        # 	If varB xS(CC) <= 0, varP(Aft) xCC <= 0, varP(Aft) xS(CC) <= 0, varT(P(Aft)) xCC > 0, then F = (MSB xCC / MSpooled),
        # where MSpooled = (SS_B xS(CC) + SS_P(Aft) xCC + SS_T(P(Aft)) xCC + SS_Res) / (df_B xS(CC) + df_P(Aft) xCC + df_T(P(Aft)) xCC + df_Res)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0048"
        c.an$Df.Formula[c.an$ID == "03AE.0048"]<-"(b-1)(c-1)"
        c.an$Calculated.From[c.an$ID == "03AE.0048"]<-"sym.C"
        c.an$Mean.Sq[c.an$ID == "03AE.0048"] <- c.an$Sum.Sq[c.an$ID == "03AE.0048"] / c.an$Df[c.an$ID == "03AE.0048"]
        c.an$Source.of.Variation[c.an$ID == "03AE.0048"]<-"B x C"
        c.an$Denominator[c.an$ID == "03AE.0048"]<-"BxS(C)+P(Aft)xC+T(P(Aft))xC+Residuals"
        c.an$Df.Pooled[c.an$ID == "03AE.0048"]<- c.an$Df[c.an$ID == "07AEG.0040"] + ac.an$Df[ac.an$ID == "08AEB.2040"] + ac.an$Df[ac.an$ID == "10AEBC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0048"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0048"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0048"] / c.an$Df.Pooled[c.an$ID == "03AE.0048"]
        c.an$F.Value[c.an$ID == "03AE.0048"]<-c.an$Mean.Sq[c.an$ID == "03AE.0048"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0048"]
        c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0048"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0048"], c.an$Df[c.an$ID == "03AE.0048"], c.an$Df.Pooled[c.an$ID == "03AE.0048"], lower.tail = F)

        # 	If varB xS(CC) <= 0, varP(Aft) xCC <= 0, varP(Aft) xS(CC) <= 0, varT(P(Aft)) xCC <= 0, then F = (MSB xCC / MSpooled),
        #where MSpooled = (SS_B xS(CC) + SS_P(Aft) xCC + SS_P(Aft) xS(CC) + SS_T(P(Aft)) xCC + SS_T(P(Aft)) xS(CC) + SS_Res) / (df_B xS(CC) + df_P(Aft) xCC + df_P(Aft) xS(CC) + df_T(P(Aft)) xCC+ df_T(P(Aft)) xS(CC) + df_Res)
        c.an<-dplyr::bind_rows(c.an,sym.C[sym.C$Source.of.Variation == "A:E",1:3])
        c.an$ID[c.an$Source.of.Variation == "A:E"]<-"03AE.0049"
        c.an$Df.Formula[c.an$ID == "03AE.0049"]<-"(b-1)(c-1)"
        c.an$Calculated.From[c.an$ID == "03AE.0049"]<-"sym.C"
        c.an$Mean.Sq[c.an$ID == "03AE.0049"] <- c.an$Sum.Sq[c.an$ID == "03AE.0049"] / c.an$Df[c.an$ID == "03AE.0049"]
        c.an$Source.of.Variation[c.an$ID == "03AE.0049"]<-"B x C"
        c.an$Denominator[c.an$ID == "03AE.0049"]<-"BxS(C)+P(Aft)xC+P(Aft)xS(C)+T(P(Aft))xC+T(P(Aft))xS(C)+Residuals"
        c.an$Df.Pooled[c.an$ID == "03AE.0049"]<- c.an$Df[c.an$ID == "07AEG.0040"] + ac.an$Df[ac.an$ID == "08AEB.2040"] + ac.an$Df[ac.an$ID == "09AEBG.2040"] + ac.an$Df[ac.an$ID == "10AEBC.2040"] + ac.an$Df[ac.an$ID == "11AEBGC.2040"] + f.an$Df[f.an$ID == "12R.0000"]
        c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0049"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"] + ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] + ac.an$Sum.Sq[ac.an$ID == "11AEBGC.2040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0049"] <- c.an$Sum.Sq.Pooled[c.an$ID == "03AE.0049"] / c.an$Df.Pooled[c.an$ID == "03AE.0049"]
        c.an$F.Value[c.an$ID == "03AE.0049"]<-c.an$Mean.Sq[c.an$ID == "03AE.0049"] / c.an$Mean.Sq.Pooled[c.an$ID == "03AE.0049"]
        c.an$P.Value.Upper.Tail[c.an$ID == "03AE.0049"]<- stats::pf(c.an$F.Value[c.an$ID == "03AE.0049"], c.an$Df[c.an$ID == "03AE.0049"], c.an$Df.Pooled[c.an$ID == "03AE.0049"], lower.tail = F)

        # 	If varB xCC <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SS_B xCC + SS_Res) / (df_B xCC + df_Res)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0031"
        i.an$Df.Formula[i.an$ID == "03AE.0031"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0031"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0031"]<-i.an$Sum.Sq[i.an$ID == "03AE.0031"] / i.an$Df[i.an$ID == "03AE.0031"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0031"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0031"]<-"BxC+Residuals"
        i.an$Df.Pooled[i.an$ID == "03AE.0031"]<- c.an$Df[c.an$ID == "03AE.0040"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"] / i.an$Df.Pooled[i.an$ID == "03AE.0031"]
        i.an$F.Value[i.an$ID == "03AE.0031"]<-i.an$Mean.Sq[i.an$ID == "03AE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0031"], i.an$Df[i.an$ID == "03AE.0031"], i.an$Df.Pooled[i.an$ID == "03AE.0031"], lower.tail = F)
      }
      if(method == "pool") {
        # 	If varT(P(Aft)) xCC <= 0 and varT(P(Aft)) xS(IC) > 0, then F = (MS_T(P(Aft)) xIC / MS_pooled),
        # where MS_pooled = (SS_T(P(Aft)) xCC + SS_T(P(Aft)) xS(IC)) / (df_T(P(Aft)) xCC + df_T(P(Aft)) xS(IC)) = SS1/df1
        ai.an<-dplyr::bind_rows(ai.an,di10AEBC.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E:C"] <- "10AEBC.2031"
        ai.an$Df.Formula[ai.an$ID == "10AEBC.2031"]<-"p(t-1)"
        ai.an$Calculated.From[ai.an$ID == "10AEBC.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "10AEBC.2031"]<-ai.an$Sum.Sq[ai.an$ID == "10AEBC.2031"] / ai.an$Df[ai.an$ID == "10AEBC.2031"]
        ai.an$Source.of.Variation[ai.an$ID == "10AEBC.2031"] <-"T(P(Aft)) x I"
        ai.an$Denominator[ai.an$ID == "10AEBC.2031"]<-"T(P(Aft))xC+T(P(Aft))xS(I)"
        ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"]<- ac.an$Df[ac.an$ID == "10AEBC.2040"] + ai.an$Df[ai.an$ID == "11AEBGC.2030"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEBC.2031"]<- ac.an$Sum.Sq[ac.an$ID == "10AEBC.2040"] + ai.an$Sum.Sq[ai.an$ID == "11AEBGC.2030"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEBC.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEBC.2031"] / ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"]
        ai.an$F.Value[ai.an$ID == "10AEBC.2031"]<-ai.an$Mean.Sq[ai.an$ID == "10AEBC.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEBC.2031"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEBC.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEBC.2031"], ai.an$Df[ai.an$ID == "10AEBC.2031"], ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"], lower.tail = F)

        # 	If varT(P(Aft)) xCC <= 0 and varT(P(Aft)) xS(IC) <= 0, then F = (MST(P(Aft)) xIC / MSpooled), where MSpooled = SS1 + SSRes / df1 + dfRes)
        ai.an<-dplyr::bind_rows(ai.an,di10AEBC.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E:C"] <- "10AEBC.2032"
        ai.an$Df.Formula[ai.an$ID == "10AEBC.2032"]<-"p(t-1)"
        ai.an$Calculated.From[ai.an$ID == "10AEBC.2032"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "10AEBC.2032"]<-ai.an$Sum.Sq[ai.an$ID == "10AEBC.2032"] / ai.an$Df[ai.an$ID == "10AEBC.2032"]
        ai.an$Source.of.Variation[ai.an$ID == "10AEBC.2032"] <-"T(P(Aft)) x I"
        ai.an$Denominator[ai.an$ID == "10AEBC.2032"]<-"T(P(Aft))xC+T(P(Aft))xS(I)+Residuals"
        ai.an$Df.Pooled[ai.an$ID == "10AEBC.2032"]<- ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"] + f.an$Df[f.an$ID == "12R.0000"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEBC.2032"]<-  ai.an$Df.Pooled[ai.an$ID == "10AEBC.2031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEBC.2032"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "10AEBC.2032"] / ai.an$Df.Pooled[ai.an$ID == "10AEBC.2032"]
        ai.an$F.Value[ai.an$ID == "10AEBC.2032"]<-ai.an$Mean.Sq[ai.an$ID == "10AEBC.2032"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "10AEBC.2032"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "10AEBC.2032"]<- stats::pf(ai.an$F.Value[ai.an$ID == "10AEBC.2032"], ai.an$Df[ai.an$ID == "10AEBC.2032"], ai.an$Df.Pooled[ai.an$ID == "10AEBC.2032"], lower.tail = F)

        # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xS(IC) > 0, then F = (MS_P(Aft) xS(IC) / MS_pooled),
        #where MSpooled = (SS_P(Aft) xS(CC) + SS_T(P(Aft)) xS(IC)) / (df_P(Aft) xS(CC) + df_T(P(Aft)) xS(IC)) = SS1/df1
        ai.an<-dplyr::bind_rows(ai.an,di09AEBG.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E:G"] <- "09AEBG.2031"
        ai.an$Df.Formula[ai.an$ID == "09AEBG.2031"]<-"(p-1)(s-1)"
        ai.an$Calculated.From[ai.an$ID == "09AEBG.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "09AEBG.2031"]<-ai.an$Sum.Sq[ai.an$ID == "09AEBG.2031"] / ai.an$Df[ai.an$ID == "09AEBG.2031"]
        ai.an$Source.of.Variation[ai.an$ID == "09AEBG.2031"] <-"P(Aft) x S(I)"
        ai.an$Denominator[ai.an$ID == "09AEBG.2031"]<-"P(Aft)xS(C)+T(P(Aft))xS(I)"
        ai.an$Df.Pooled[ai.an$ID == "09AEBG.2031"]<- ac.an$Df[ac.an$ID == "09AEBG.2040"] + ai.an$Df[ai.an$ID == "11AEBGC.2030"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "09AEBG.2031"]<- ac.an$Sum.Sq[ac.an$ID == "09AEBG.2040"] + ai.an$Sum.Sq[ai.an$ID == "11AEBGC.2030"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "09AEBG.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "09AEBG.2031"] / ai.an$Df.Pooled[ai.an$ID == "09AEBG.2031"]
        ai.an$F.Value[ai.an$ID == "09AEBG.2031"]<-ai.an$Mean.Sq[ai.an$ID == "09AEBG.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "09AEBG.2031"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "09AEBG.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "09AEBG.2031"], ai.an$Df[ai.an$ID == "09AEBG.2031"], ai.an$Df.Pooled[ai.an$ID == "09AEBG.2031"], lower.tail = F)

        # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xS(IC) <= 0, then F = (MS_P(Aft) xS(IC) / MS_pooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        ai.an<-dplyr::bind_rows(ai.an,di09AEBG.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E:G"] <- "09AEBG.2032"
        ai.an$Df.Formula[ai.an$ID == "09AEBG.2032"]<-"(p-1)(s-1)"
        ai.an$Calculated.From[ai.an$ID == "09AEBG.2032"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "09AEBG.2032"]<-ai.an$Sum.Sq[ai.an$ID == "09AEBG.2032"] / ai.an$Df[ai.an$ID == "09AEBG.2032"]
        ai.an$Source.of.Variation[ai.an$ID == "09AEBG.2032"] <-"P(Aft) x S(I)"
        ai.an$Denominator[ai.an$ID == "09AEBG.2032"]<-"P(Aft)xS(C)+T(P(Aft))xS(I)+Residuals)"
        ai.an$Df.Pooled[ai.an$ID == "09AEBG.2032"]<-  ai.an$Df.Pooled[ai.an$ID == "09AEBG.2031"] + f.an$Df[f.an$ID == "12R.0000"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "09AEBG.2032"]<- ai.an$Sum.Sq.Pooled[ai.an$ID == "09AEBG.2031"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "09AEBG.2032"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "09AEBG.2032"] / ai.an$Df.Pooled[ai.an$ID == "09AEBG.2032"]
        ai.an$F.Value[ai.an$ID == "09AEBG.2032"]<-ai.an$Mean.Sq[ai.an$ID == "09AEBG.2032"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "09AEBG.2032"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "09AEBG.2032"]<- stats::pf(ai.an$F.Value[ai.an$ID == "09AEBG.2032"], ai.an$Df[ai.an$ID == "09AEBG.2032"], ai.an$Df.Pooled[ai.an$ID == "09AEBG.2032"], lower.tail = F)

        # 	If varP(Aft) xCC <= 0, varT(P(Aft)) xIC <= and varP(Aft) xS(IC) > 0, then F = (MS_P(Aft) xIC / MSpooled), where MSpooled = (SS_P(Aft) xCC + SS_P(Aft) xS(IC)) / (df_P(Aft) xCC + df_P(Aft) xS(IC))
        ai.an<-dplyr::bind_rows(ai.an,di08AEB.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E"] <- "08AEB.2031"
        ai.an$Df.Formula[ai.an$ID == "08AEB.2031"]<-"(p-1)"
        ai.an$Calculated.From[ai.an$ID == "08AEB.2031"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "08AEB.2031"]<-ai.an$Sum.Sq[ai.an$ID == "08AEB.2031"] / ai.an$Df[ai.an$ID == "08AEB.2031"]
        ai.an$Source.of.Variation[ai.an$ID == "08AEB.2031"] <-"P(Aft) x I"
        ai.an$Denominator[ai.an$ID == "08AEB.2031"]<-"P(Aft)xC+P(Aft)xS(I)"
        ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"]<- ac.an$Df[ac.an$ID == "08AEB.2040"] + ai.an$Df[ai.an$ID == "09AEBG.2030"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2031"]<- ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + ai.an$Sum.Sq[ai.an$ID == "09AEBG.2030"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2031"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2031"] / ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"]
        ai.an$F.Value[ai.an$ID == "08AEB.2031"]<-ai.an$Mean.Sq[ai.an$ID == "08AEB.2031"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2031"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "08AEB.2031"]<- stats::pf(ai.an$F.Value[ai.an$ID == "08AEB.2031"], ai.an$Df[ai.an$ID == "08AEB.2031"], ai.an$Df.Pooled[ai.an$ID == "08AEB.2031"], lower.tail = F)

        # 	If varP(Aft) xCC <= 0, varP(Aft) xS(IC) <= 0 and varT(P(Aft)) xIC > 0, then F = (MSP(Aft) xIC / MSpooled), where MSpooled = (SS_P(Aft) xCC + SS_T(P(Aft)) xIC) / (df_P(Aft) xCC + df_T(P(Aft)) xIC)
        ai.an<-dplyr::bind_rows(ai.an,di08AEB.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E"] <- "08AEB.2032"
        ai.an$Df.Formula[ai.an$ID == "08AEB.2032"]<-"(p-1)"
        ai.an$Calculated.From[ai.an$ID == "08AEB.2032"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "08AEB.2032"]<-ai.an$Sum.Sq[ai.an$ID == "08AEB.2032"] / ai.an$Df[ai.an$ID == "08AEB.2032"]
        ai.an$Source.of.Variation[ai.an$ID == "08AEB.2032"] <-"P(Aft) x I"
        ai.an$Denominator[ai.an$ID == "08AEB.2032"]<-"P(Aft)xC+T(P(Aft))xI"
        ai.an$Df.Pooled[ai.an$ID == "08AEB.2032"]<- ac.an$Df[ac.an$ID == "08AEB.2040"] + ai.an$Df[ai.an$ID == "10AEBC.2030"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2032"]<- ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + ai.an$Sum.Sq[ai.an$ID == "10AEBC.2030"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2032"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2032"] / ai.an$Df.Pooled[ai.an$ID == "08AEB.2032"]
        ai.an$F.Value[ai.an$ID == "08AEB.2032"]<-ai.an$Mean.Sq[ai.an$ID == "08AEB.2032"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2032"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "08AEB.2032"]<- stats::pf(ai.an$F.Value[ai.an$ID == "08AEB.2032"], ai.an$Df[ai.an$ID == "08AEB.2032"], ai.an$Df.Pooled[ai.an$ID == "08AEB.2032"], lower.tail = F)

        # 	If varP(Aft) xCC <= 0, varP(Aft) xS(IC) <= 0, varT(P(Aft)) xIC <= 0, and varT(P(Aft)) xS(IC) > 0, then F = (MSP(Aft) xIC / MSpooled), where MSpooled = (SS_P(Aft) xCC + SS_P(Aft) xS(IC) + SS_T(P(Aft)) xIC + SS_T(P(Aft)) xS(IC)) / (dfP(Aft) xCC + dfP(Aft) xS(IC) + dfT(P(Aft)) xIC + dfT(P(Aft)) xS(IC)) = SS1/df1
        ai.an<-dplyr::bind_rows(ai.an,di08AEB.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E"] <- "08AEB.2033"
        ai.an$Df.Formula[ai.an$ID == "08AEB.2033"]<-"(p-1)"
        ai.an$Calculated.From[ai.an$ID == "08AEB.2033"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "08AEB.2033"]<-ai.an$Sum.Sq[ai.an$ID == "08AEB.2033"] / ai.an$Df[ai.an$ID == "08AEB.2033"]
        ai.an$Source.of.Variation[ai.an$ID == "08AEB.2033"] <-"P(Aft) x I"
        ai.an$Denominator[ai.an$ID == "08AEB.2033"]<-"P(Aft)xC+P(Aft)xS(I)+T(P(Aft))xI+T(P(Aft))xS(I)"
        ai.an$Df.Pooled[ai.an$ID == "08AEB.2033"]<- ac.an$Df[ac.an$ID == "08AEB.2040"] + ai.an$Df[ai.an$ID == "09AEBG.2030"] + ai.an$Df[ai.an$ID == "10AEBC.2030"] + ai.an$Df[ai.an$ID == "11AEBGC.2030"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2033"]<- ac.an$Sum.Sq[ac.an$ID == "08AEB.2040"] + ai.an$Sum.Sq[ai.an$ID == "09AEBG.2030"] + ai.an$Sum.Sq[ai.an$ID == "10AEBC.2030"] + ai.an$Sum.Sq[ai.an$ID == "11AEBGC.2030"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2033"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2033"] / ai.an$Df.Pooled[ai.an$ID == "08AEB.2033"]
        ai.an$F.Value[ai.an$ID == "08AEB.2033"]<-ai.an$Mean.Sq[ai.an$ID == "08AEB.2033"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2033"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "08AEB.2033"]<- stats::pf(ai.an$F.Value[ai.an$ID == "08AEB.2033"], ai.an$Df[ai.an$ID == "08AEB.2033"], ai.an$Df.Pooled[ai.an$ID == "08AEB.2033"], lower.tail = F)

        # 	If varP(Aft) xCC <= 0, varT(P(Aft)) xIC <= 0, varP(Aft) xS(IC) <= 0) and varT(P(Aft)) xS(IC) <= 0, then F = (MSP(Aft) xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        ai.an<-dplyr::bind_rows(ai.an,di08AEB.2030)
        ai.an$ID[ai.an$Source.of.Variation == "B:E"] <- "08AEB.2034"
        ai.an$Df.Formula[ai.an$ID == "08AEB.2034"]<-"(p-1)"
        ai.an$Calculated.From[ai.an$ID == "08AEB.2034"]<-"sym.F-sym.B-sym.C+sym.BC"
        ai.an$Mean.Sq[ai.an$ID == "08AEB.2034"]<-ai.an$Sum.Sq[ai.an$ID == "08AEB.2034"] / ai.an$Df[ai.an$ID == "08AEB.2034"]
        ai.an$Source.of.Variation[ai.an$ID == "08AEB.2034"] <-"P(Aft) x I"
        ai.an$Denominator[ai.an$ID == "08AEB.2034"]<-"P(Aft)xC+P(Aft)xS(I)+T(P(Aft))xI+T(P(Aft))xS(I)+Residuals"
        ai.an$Df.Pooled[ai.an$ID == "08AEB.2034"]<- ai.an$Df.Pooled[ai.an$ID == "08AEB.2033"] + f.an$Df[f.an$ID == "12R.0000"]
        ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2034"]<- ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2033"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2034"] <- ai.an$Sum.Sq.Pooled[ai.an$ID == "08AEB.2034"] / ai.an$Df.Pooled[ai.an$ID == "08AEB.2034"]
        ai.an$F.Value[ai.an$ID == "08AEB.2034"]<-ai.an$Mean.Sq[ai.an$ID == "08AEB.2034"] / ai.an$Mean.Sq.Pooled[ai.an$ID == "08AEB.2034"]
        ai.an$P.Value.Upper.Tail[ai.an$ID == "08AEB.2034"]<- stats::pf(ai.an$F.Value[ai.an$ID == "08AEB.2034"], ai.an$Df[ai.an$ID == "08AEB.2034"], ai.an$Df.Pooled[ai.an$ID == "08AEB.2034"], lower.tail = F)

        # 	If varB xS(CC) <= 0 and varP(Aft) xS(IC) > 0, then F = (MSB xS(IC) / MSpooled), where MSpooled = (SSB xS(CC) + SS_P(Aft) xS(IC)) / (dfB xS(CC) + dfP(Aft) xS(IC)) = SS1/df1
        i.an<-dplyr::bind_rows(i.an,di07AEG.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E:G"] <- "07AEG.0031"
        i.an$Df.Formula[i.an$ID == "07AEG.0031"]<-"(b-1)(s-1)"
        i.an$Calculated.From[i.an$ID == "07AEG.0031"]<-"sym.F-sym.B-sym.C+sym.BC"
        i.an$Mean.Sq[i.an$ID == "07AEG.0031"]<-i.an$Sum.Sq[i.an$ID == "07AEG.0031"] / i.an$Df[i.an$ID == "07AEG.0031"]
        i.an$Source.of.Variation[i.an$ID == "07AEG.0031"] <-"B x S(I)"
        i.an$Denominator[i.an$ID == "07AEG.0031"]<-"BxS(C)+P(Aft)xS(I)"
        i.an$Df.Pooled[i.an$ID == "07AEG.0031"]<- c.an$Df[c.an$ID == "07AEG.0040"] + ai.an$Df[ai.an$ID == "09AEBG.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0031"]<- c.an$Sum.Sq[c.an$ID == "07AEG.0040"] + ai.an$Sum.Sq[ai.an$ID == "09AEBG.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0031"] / i.an$Df.Pooled[i.an$ID == "07AEG.0031"]
        i.an$F.Value[i.an$ID == "07AEG.0031"]<-i.an$Mean.Sq[i.an$ID == "07AEG.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "07AEG.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "07AEG.0031"], i.an$Df[i.an$ID == "07AEG.0031"], i.an$Df.Pooled[i.an$ID == "07AEG.0031"], lower.tail = F)

        # 	If varB xS(CC) <= 0, varP(Aft) xS(IC) <= 0 and varT(P(Aft)) xS(IC) > 0, then F = (MSB xS(IC) / MSpooled), where MSpooled = (SS1 + SS_T(P(Aft)) xS(IC)) / (df_1 + dfT(P(Aft)) xS(IC)) = SS2/df2
        i.an<-dplyr::bind_rows(i.an,di07AEG.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E:G"] <- "07AEG.0032"
        i.an$Df.Formula[i.an$ID == "07AEG.0032"]<-"(b-1)(s-1)"
        i.an$Calculated.From[i.an$ID == "07AEG.0032"]<-"sym.F-sym.B-sym.C+sym.BC"
        i.an$Mean.Sq[i.an$ID == "07AEG.0032"]<-i.an$Sum.Sq[i.an$ID == "07AEG.0032"] / i.an$Df[i.an$ID == "07AEG.0032"]
        i.an$Source.of.Variation[i.an$ID == "07AEG.0032"] <-"B x S(I)"
        i.an$Denominator[i.an$ID == "07AEG.0032"]<-"BxS(C)+P(Aft)xS(I)+T(P(Aft))xS(I)"
        i.an$Df.Pooled[i.an$ID == "07AEG.0032"]<- i.an$Df.Pooled[i.an$ID == "07AEG.0031"] + ai.an$Df[ai.an$ID == "11AEBGC.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0032"]<- i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0031"] + ai.an$Sum.Sq[ai.an$ID == "11AEBGC.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0032"] / i.an$Df.Pooled[i.an$ID == "07AEG.0032"]
        i.an$F.Value[i.an$ID == "07AEG.0032"]<-i.an$Mean.Sq[i.an$ID == "07AEG.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0032"]
        i.an$P.Value.Upper.Tail[i.an$ID == "07AEG.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "07AEG.0032"], i.an$Df[i.an$ID == "07AEG.0032"], i.an$Df.Pooled[i.an$ID == "07AEG.0032"], lower.tail = F)

        # 	If varB xS(CC) <= 0, varP(Aft) xS(IC) <= 0 and varT(P(Aft)) xS(IC) <= 0, then F = (MSB xS(IC) / MSpooled), where MSpooled = (SS2 + SSRes) / (df2 + dfRes)
        i.an<-dplyr::bind_rows(i.an,di07AEG.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E:G"] <- "07AEG.0033"
        i.an$Df.Formula[i.an$ID == "07AEG.0033"]<-"(b-1)(s-1)"
        i.an$Calculated.From[i.an$ID == "07AEG.0033"]<-"sym.F-sym.B-sym.C+sym.BC"
        i.an$Mean.Sq[i.an$ID == "07AEG.0033"]<-i.an$Sum.Sq[i.an$ID == "07AEG.0033"] / i.an$Df[i.an$ID == "07AEG.0033"]
        i.an$Source.of.Variation[i.an$ID == "07AEG.0033"] <-"B x S(I)"
        i.an$Denominator[i.an$ID == "07AEG.0033"]<-"BxS(C)+P(Aft)xS(I)+T(P(Aft))xS(I)+Residuals"
        i.an$Df.Pooled[i.an$ID == "07AEG.0033"]<- i.an$Df.Pooled[i.an$ID == "07AEG.0032"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0033"]<- i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0032"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0033"] <- i.an$Sum.Sq.Pooled[i.an$ID == "07AEG.0032"] / i.an$Df.Pooled[i.an$ID == "07AEG.0033"]
        i.an$F.Value[i.an$ID == "07AEG.0033"]<-i.an$Mean.Sq[i.an$ID == "07AEG.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "07AEG.0033"]
        i.an$P.Value.Upper.Tail[i.an$ID == "07AEG.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "07AEG.0033"], i.an$Df[i.an$ID == "07AEG.0033"], i.an$Df.Pooled[i.an$ID == "07AEG.0033"], lower.tail = F)

        # 	If varB xCC <= 0, varP(Aft) xIC <= 0 and varB xS(IC) > 0 then F = (MSB xIC / MSpooled),
        # where MSpooled = (SS_B xCC + SS_B xS(IC)) / (df_B xCC + df_B xS(IC))
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0031"
        i.an$Df.Formula[i.an$ID == "03AE.0031"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0031"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0031"]<-i.an$Sum.Sq[i.an$ID == "03AE.0031"] / i.an$Df[i.an$ID == "03AE.0031"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0031"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0031"]<-"BxC+BxS(I)"
        i.an$Df.Pooled[i.an$ID == "03AE.0031"]<- c.an$Df[c.an$ID == "03AE.0040"] + i.an$Df[i.an$ID == "07AEG.0030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + i.an$Sum.Sq[i.an$ID == "07AEG.0030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0031"] / i.an$Df.Pooled[i.an$ID == "03AE.0031"]
        i.an$F.Value[i.an$ID == "03AE.0031"]<-i.an$Mean.Sq[i.an$ID == "03AE.0031"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0031"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0031"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0031"], i.an$Df[i.an$ID == "03AE.0031"], i.an$Df.Pooled[i.an$ID == "03AE.0031"], lower.tail = F)

        # 	If varB xCC <= 0, varB xS(IC) <= 0 and varP(Aft) xIC > 0, then F = (MSB xIC / MSpooled),
        # where MSpooled = (SS_B xCC + SS_P(Aft) xIC) / (dfB xCC + dfP(Aft) xIC)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0032"
        i.an$Df.Formula[i.an$ID == "03AE.0032"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0032"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0032"]<-i.an$Sum.Sq[i.an$ID == "03AE.0032"] / i.an$Df[i.an$ID == "03AE.0032"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0032"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0032"]<-"BxC+P(Aft)xI"
        i.an$Df.Pooled[i.an$ID == "03AE.0032"]<- c.an$Df[c.an$ID == "03AE.0040"] + ai.an$Df[ai.an$ID == "08AEB.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0032"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + ai.an$Sum.Sq[ai.an$ID == "08AEB.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0032"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0032"] / i.an$Df.Pooled[i.an$ID == "03AE.0032"]
        i.an$F.Value[i.an$ID == "03AE.0032"]<-i.an$Mean.Sq[i.an$ID == "03AE.0032"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0032"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0032"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0032"], i.an$Df[i.an$ID == "03AE.0032"], i.an$Df.Pooled[i.an$ID == "03AE.0032"], lower.tail = F)

        # 	If varB xCC <= 0, varB xS(IC) <= 0, varP(Aft) xIC <= 0, varP(Aft) xS(IC) <= 0 and varT(P(Aft)) xIC > 0, then F = (MSB xIC / MSpooled),
        # where MSpooled = (SS_B xCC + SS_B xS(IC) + SS_P(Aft) xIC + SS_P(Aft) xS(IC)) / (df_B xCC + df_B xS(IC) + df_P(Aft) xIC + df_P(Aft) xS(IC))
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0033"
        i.an$Df.Formula[i.an$ID == "03AE.0033"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0033"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0033"]<-i.an$Sum.Sq[i.an$ID == "03AE.0033"] / i.an$Df[i.an$ID == "03AE.0033"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0033"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0033"]<-"BxC+BxS(I)+P(Aft)xI+P(Aft)xS(I)"
        i.an$Df.Pooled[i.an$ID == "03AE.0033"]<- c.an$Df[c.an$ID == "03AE.0040"] + i.an$Df[i.an$ID == "07AEG.0030"] + ai.an$Df[ai.an$ID == "08AEB.2030"] + ai.an$Df[ai.an$ID == "09AEBG.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0033"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + i.an$Sum.Sq[i.an$ID == "07AEG.0030"] + ai.an$Sum.Sq[ai.an$ID == "08AEB.2030"] + ai.an$Sum.Sq[ai.an$ID == "09AEBG.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0033"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0033"] / i.an$Df.Pooled[i.an$ID == "03AE.0033"]
        i.an$F.Value[i.an$ID == "03AE.0033"]<-i.an$Mean.Sq[i.an$ID == "03AE.0033"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0033"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0033"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0033"], i.an$Df[i.an$ID == "03AE.0033"], i.an$Df.Pooled[i.an$ID == "03AE.0033"], lower.tail = F)

        # 	If varB xCC <= 0, varB xS(IC) <= 0, varP(Aft) xIC <= 0, varT(P(Aft)) xIC <= 0 and varP(Aft) xS(IC) > 0, then F = (MSB xIC / MSpooled),
        # where MSpooled = (SS_B xCC + SS_B xS(IC) + SS_P(Aft) xIC + SS_T(P(Aft)) xIC) / (df_B xCC + df_B xS(IC) + df_P(Aft) xIC + df_T(P(Aft)) xIC)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0034"
        i.an$Df.Formula[i.an$ID == "03AE.0034"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0034"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0034"]<-i.an$Sum.Sq[i.an$ID == "03AE.0034"] / i.an$Df[i.an$ID == "03AE.0034"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0034"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0034"]<-"BxC+BxS(I)+P(Aft)xI+T(P(Aft))xI"
        i.an$Df.Pooled[i.an$ID == "03AE.0034"]<- c.an$Df[c.an$ID == "03AE.0040"] + i.an$Df[i.an$ID == "07AEG.0030"] + ai.an$Df[ai.an$ID == "08AEB.2030"] + ai.an$Df[ai.an$ID == "10AEBC.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0034"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + i.an$Sum.Sq[i.an$ID == "07AEG.0030"] + ai.an$Sum.Sq[ai.an$ID == "08AEB.2030"] + ai.an$Sum.Sq[ai.an$ID == "10AEBC.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0034"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0034"] / i.an$Df.Pooled[i.an$ID == "03AE.0034"]
        i.an$F.Value[i.an$ID == "03AE.0034"]<-i.an$Mean.Sq[i.an$ID == "03AE.0034"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0034"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0034"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0034"], i.an$Df[i.an$ID == "03AE.0034"], i.an$Df.Pooled[i.an$ID == "03AE.0034"], lower.tail = F)

        # 	If varB xCC <= 0, varB xS(IC) <= 0, varP(Aft) xIC <= 0, varP(Aft) xS(IC) <= 0, varT(P(Aft)) xIC <= 0 and varT(P(Aft)) xS(IC) > 0, then F = (MSB xIC / MSpooled),
        # where MSpooled = (SS_B xCC + SS_B xS(IC) + SS_P(Aft) xIC + SS_P(Aft) xS(IC) + SS_T(P(Aft)) xIC + SS_T(P(Aft)) xS(IC)) / (dfB xCC + dfB xS(IC) + dfP(Aft) xIC + dfP(Aft) xS(IC) + dfT(P(Aft)) xIC + dfT(P(Aft)) xS(IC)) = SS1/df1
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0035"
        i.an$Df.Formula[i.an$ID == "03AE.0035"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0035"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0035"]<-i.an$Sum.Sq[i.an$ID == "03AE.0035"] / i.an$Df[i.an$ID == "03AE.0035"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0035"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0035"]<-"BxC+BxS(I)+P(Aft)xI+P(Aft)xS(I)+T(P(Aft))xI+T(P(Aft))xS(I)"
        i.an$Df.Pooled[i.an$ID == "03AE.0035"]<- c.an$Df[c.an$ID == "03AE.0040"] + i.an$Df[i.an$ID == "07AEG.0030"] + ai.an$Df[ai.an$ID == "08AEB.2030"] + ai.an$Df[ai.an$ID == "09AEBG.2030"] + ai.an$Df[ai.an$ID == "10AEBC.2030"] + ai.an$Df[ai.an$ID == "11AEBGC.2030"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0035"]<- c.an$Sum.Sq[c.an$ID == "03AE.0040"] + i.an$Sum.Sq[i.an$ID == "07AEG.0030"] + ai.an$Sum.Sq[ai.an$ID == "08AEB.2030"] + ai.an$Sum.Sq[ai.an$ID == "09AEBG.2030"] + ai.an$Sum.Sq[ai.an$ID == "10AEBC.2030"] + ai.an$Sum.Sq[ai.an$ID == "11AEBGC.2030"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0035"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0035"] / i.an$Df.Pooled[i.an$ID == "03AE.0035"]
        i.an$F.Value[i.an$ID == "03AE.0035"]<-i.an$Mean.Sq[i.an$ID == "03AE.0035"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0035"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0035"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0035"], i.an$Df[i.an$ID == "03AE.0035"], i.an$Df.Pooled[i.an$ID == "03AE.0035"], lower.tail = F)

        # 	If varB xCC <= 0, varB xS(IC) <= 0, varP(Aft) xIC <= 0, varP(Aft) xS(IC) <= 0, varT(P(Aft)) xIC <= 0 and varT(P(Aft)) xS(IC) <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        i.an<-dplyr::bind_rows(i.an,di03AE.0030)
        i.an$ID[i.an$Source.of.Variation == "A:E"] <- "03AE.0036"
        i.an$Df.Formula[i.an$ID == "03AE.0036"]<-"(b-1)"
        i.an$Calculated.From[i.an$ID == "03AE.0036"]<-"sym.F-sym.C"
        i.an$Mean.Sq[i.an$ID == "03AE.0036"]<-i.an$Sum.Sq[i.an$ID == "03AE.0036"] / i.an$Df[i.an$ID == "03AE.0036"]
        i.an$Source.of.Variation[i.an$ID == "03AE.0036"] <-"B x I"
        i.an$Denominator[i.an$ID == "03AE.0036"]<-"BxC+BxS(I)+P(Aft)xI+P(Aft)xS(I)+T(P(Aft))xI+T(P(Aft))xS(I)+Residuals"
        i.an$Df.Pooled[i.an$ID == "03AE.0036"]<- i.an$Df.Pooled[i.an$ID == "03AE.0035"] + f.an$Df[f.an$ID == "12R.0000"]
        i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0036"]<- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0035"] + f.an$Sum.Sq[f.an$ID == "12R.0000"]
        i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0036"] <- i.an$Sum.Sq.Pooled[i.an$ID == "03AE.0036"] / i.an$Df.Pooled[i.an$ID == "03AE.0036"]
        i.an$F.Value[i.an$ID == "03AE.0036"]<-i.an$Mean.Sq[i.an$ID == "03AE.0036"] / i.an$Mean.Sq.Pooled[i.an$ID == "03AE.0036"]
        i.an$P.Value.Upper.Tail[i.an$ID == "03AE.0036"]<- stats::pf(i.an$F.Value[i.an$ID == "03AE.0036"], i.an$Df[i.an$ID == "03AE.0036"], i.an$Df.Pooled[i.an$ID == "03AE.0036"], lower.tail = F)
      }
    }
  }


  ###################################################################################
  ### NO IF CHAMA A COMBINACAO LM ORIGINAL, MAS NA ANOVA PARCIAL ESTA SEM O "A"


  ###--------------------------------- AMONG Before and Controls = From sym.CB
  ##--------------------- code xx.1040

  bc.an<-baci[F,]
  ###NO IF CHAMA A COMBINACAO LM ORIGINAL, MAS NA ANOVA PARCIAL ESTA SEM O "A"

  if("AEB" %in% aav.model.lm == T) {
    bc.an<-dplyr::bind_rows(bc.an,sym.CB[sym.CB$Source.of.Variation == "B:E",1:3])
    bc.an$ID[bc.an$Source.of.Variation == "B:E"] <- "08AEB.1040"
    bc.an$Calculated.From[bc.an$ID == "08AEB.1040"]<-"sym.CB"
    bc.an$Df.Formula[bc.an$ID == "08AEB.1040"]<-"(p-1)(c-1)"
    bc.an$Mean.Sq[bc.an$ID == "08AEB.1040"]<-bc.an$Sum.Sq[bc.an$ID == "08AEB.1040"] / bc.an$Df[bc.an$ID == "08AEB.1040"]
    bc.an$Source.of.Variation[bc.an$ID == "08AEB.1040"] <-"P(Bef) x C"
    # 2-tail
    # 	Two-tailed test F = (MS_P(Bef) xCC / MS_P(Aft) xCC)
    bc.an$Denominator[bc.an$ID == "08AEB.1040"]<-"P(Aft)xC"
    bc.an$F.Value[bc.an$ID == "08AEB.1040"]<-bc.an$Mean.Sq[bc.an$ID == "08AEB.1040"] / ac.an$Mean.Sq[ac.an$ID == "08AEB.2040"]
    bc.an$P.Value.Upper.Tail[bc.an$ID == "08AEB.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "08AEB.1040"], bc.an$Df[bc.an$ID == "08AEB.1040"], ac.an$Df[ac.an$ID == "08AEB.2040"], lower.tail = F)
    bc.an$P.Value.Lower.Tail[bc.an$ID == "08AEB.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "08AEB.1040"], bc.an$Df[bc.an$ID == "08AEB.1040"], ac.an$Df[ac.an$ID == "08AEB.2040"], lower.tail = T)
  }
  if("AEC" %in% aav.model.lm == T) {
    bc.an<-dplyr::bind_rows(bc.an,sym.CB[sym.CB$Source.of.Variation == "C:E",1:3])
    bc.an$ID[bc.an$Source.of.Variation == "C:E"] <- "10AEC.1040"
    bc.an$Calculated.From[bc.an$ID == "10AEC.1040"]<-"sym.CB"
    bc.an$Df.Formula[bc.an$ID == "10AEC.1040"]<-"(t-1)(c-1)"
    bc.an$Mean.Sq[bc.an$ID == "10AEC.1040"]<-bc.an$Sum.Sq[bc.an$ID == "10AEC.1040"] / bc.an$Df[bc.an$ID == "10AEC.1040"]
    bc.an$Source.of.Variation[bc.an$ID == "10AEC.1040"] <-"T(Bef) x C"
    # 2-tail
    # 	Two-tailed test F = (MS_T(Bef) xCC / MS_T(Aft) xCC)
    bc.an$Denominator[bc.an$ID == "10AEC.1040"]<-"T(Aft)xC"
    bc.an$F.Value[bc.an$ID == "10AEC.1040"]<-bc.an$Mean.Sq[bc.an$ID == "10AEC.1040"] / ac.an$Mean.Sq[ac.an$ID == "10AEC.2040"]
    bc.an$P.Value.Upper.Tail[bc.an$ID == "10AEC.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "10AEC.1040"], bc.an$Df[bc.an$ID == "10AEC.1040"], ac.an$Df[ac.an$ID == "10AEC.2040"], lower.tail = F)
    bc.an$P.Value.Lower.Tail[bc.an$ID == "10AEC.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "10AEC.1040"], bc.an$Df[bc.an$ID == "10AEC.1040"], ac.an$Df[ac.an$ID == "10AEC.2040"], lower.tail = T)
  }
  if("AEBG" %in% aav.model.lm == T) {
    bc.an<-dplyr::bind_rows(bc.an,sym.CB[sym.CB$Source.of.Variation == "B:E:G",1:3])
    bc.an$ID[bc.an$Source.of.Variation == "B:E:G"] <- "09AEBG.1040"
    bc.an$Calculated.From[bc.an$ID == "09AEBG.1040"]<-"sym.CB"
    bc.an$Df.Formula[bc.an$ID == "09AEBG.1040"]<-"(t-1)c(s-1)"
    bc.an$Mean.Sq[bc.an$ID == "09AEBG.1040"]<-bc.an$Sum.Sq[bc.an$ID == "09AEBG.1040"] / bc.an$Df[bc.an$ID == "09AEBG.1040"]
    bc.an$Source.of.Variation[bc.an$ID == "09AEBG.1040"] <-"P(Bef) x S(C)"
    # 2-tail
    # 	Two-tailed test F = (MS_P(Bef) xS(CC) / MS_P(Aft) xS(CC))
    bc.an$Denominator[bc.an$ID == "09AEBG.1040"]<-"P(Aft)xS(C)"
    bc.an$F.Value[bc.an$ID == "09AEBG.1040"]<-bc.an$Mean.Sq[bc.an$ID == "09AEBG.1040"] / ac.an$Mean.Sq[ac.an$ID == "09AEBG.2040"]
    bc.an$P.Value.Upper.Tail[bc.an$ID == "09AEBG.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "09AEBG.1040"], bc.an$Df[bc.an$ID == "09AEBG.1040"], ac.an$Df[ac.an$ID == "09AEBG.2040"], lower.tail = F)
    bc.an$P.Value.Lower.Tail[bc.an$ID == "09AEBG.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "09AEBG.1040"], bc.an$Df[bc.an$ID == "09AEBG.1040"], ac.an$Df[ac.an$ID == "09AEBG.2040"], lower.tail = T)
  }
  if("AEBC" %in% aav.model.lm == T) {
    bc.an<-dplyr::bind_rows(bc.an,sym.CB[sym.CB$Source.of.Variation == "B:E:C",1:3])
    bc.an$ID[bc.an$Source.of.Variation == "B:E:C"] <- "10AEBC.1040"
    bc.an$Calculated.From[bc.an$ID == "10AEBC.1040"]<-"sym.CB"
    bc.an$Df.Formula[bc.an$ID == "10AEBC.1040"]<-"p(t-1)(c-1)"
    bc.an$Mean.Sq[bc.an$ID == "10AEBC.1040"]<- bc.an$Sum.Sq[bc.an$ID == "10AEBC.1040"] / bc.an$Df[bc.an$ID == "10AEBC.1040"]
    bc.an$Source.of.Variation[bc.an$ID == "10AEBC.1040"] <-"T(P(Bef)) x C"
    # 2-tail
    # 	Two-tailed test F = (MS_T(P(Bef)) xCC / MS_T(P(Aft)) xCC)
    bc.an$Denominator[bc.an$ID == "10AEBC.1040"]<-"T(P(Aft))xC"
    bc.an$F.Value[bc.an$ID == "10AEBC.1040"]<-bc.an$Mean.Sq[bc.an$ID == "10AEBC.1040"] / ac.an$Mean.Sq[ac.an$ID == "10AEBC.2040"]
    bc.an$P.Value.Upper.Tail[bc.an$ID == "10AEBC.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "10AEBC.1040"], bc.an$Df[bc.an$ID == "10AEBC.1040"], ac.an$Df[ac.an$ID == "10AEBC.2040"], lower.tail = F)
    bc.an$P.Value.Lower.Tail[bc.an$ID == "10AEBC.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "10AEBC.1040"], bc.an$Df[bc.an$ID == "10AEBC.1040"], ac.an$Df[ac.an$ID == "10AEBC.2040"], lower.tail = T)
  }
  if("AEGC" %in% aav.model.lm == T) {
    bc.an<-dplyr::bind_rows(bc.an,sym.CB[sym.CB$Source.of.Variation == "C:E:G",1:3])
    bc.an$ID[bc.an$Source.of.Variation == "C:E:G"] <- "11AEGC.1040"
    bc.an$Calculated.From[bc.an$ID == "11AEGC.1040"]<-"sym.CB"
    bc.an$Df.Formula[bc.an$ID == "11AEGC.1040"]<-"(t-1)c(s-1)"
    bc.an$Mean.Sq[bc.an$ID == "11AEGC.1040"]<-bc.an$Sum.Sq[bc.an$ID == "11AEGC.1040"] / bc.an$Df[bc.an$ID == "11AEGC.1040"]
    bc.an$Source.of.Variation[bc.an$ID == "11AEGC.1040"] <-"T(Bef) x S(C)"
    # 2-tail
    # 	Two-tailed test F = (MST(Bef) xS(CC) / MST(Aft) xS(CC))
    bc.an$Denominator[bc.an$ID == "11AEGC.1040"]<-"T(Aft)xS(C)"
    bc.an$F.Value[bc.an$ID == "11AEGC.1040"]<-bc.an$Mean.Sq[bc.an$ID == "11AEGC.1040"] / ac.an$Mean.Sq[ac.an$ID == "11AEGC.2040"]
    bc.an$P.Value.Upper.Tail[bc.an$ID == "11AEGC.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "11AEGC.1040"], bc.an$Df[bc.an$ID == "11AEGC.1040"], ac.an$Df[ac.an$ID == "11AEGC.2040"], lower.tail = F)
    bc.an$P.Value.Lower.Tail[bc.an$ID == "11AEGC.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "11AEGC.1040"], bc.an$Df[bc.an$ID == "11AEGC.1040"], ac.an$Df[ac.an$ID == "11AEGC.2040"], lower.tail = T)
  }
  if("AEBGC" %in% aav.model.lm == T) {
    bc.an<-dplyr::bind_rows(bc.an,sym.CB[sym.CB$Source.of.Variation == "B:E:G:C",1:3])
    bc.an$ID[bc.an$Source.of.Variation == "B:E:G:C"] <- "11AEBGC.1040"
    bc.an$Calculated.From[bc.an$ID == "11AEBGC.1040"]<-"sym.CB"
    bc.an$Df.Formula[bc.an$ID == "11AEBGC.1040"]<-"p(t-1)c(s-1)"
    bc.an$Mean.Sq[bc.an$ID == "11AEBGC.1040"]<- bc.an$Sum.Sq[bc.an$ID == "11AEBGC.1040"] / bc.an$Df[bc.an$ID == "11AEBGC.1040"]
    bc.an$Source.of.Variation[bc.an$ID == "11AEBGC.1040"] <-"T(P(Bef)) x S(C)"
    # 2-tail
    # 	Two-tailed test F = (MS_T(P(Bef)) xS(CC) / MS_T(P(Aft)) xS(CC))
    bc.an$Denominator[bc.an$ID == "11AEBGC.1040"]<-"T(P(Aft))xS(C)"
    bc.an$F.Value[bc.an$ID == "11AEBGC.1040"]<-bc.an$Mean.Sq[bc.an$ID == "11AEBGC.1040"] / ac.an$Mean.Sq[ac.an$ID == "11AEBGC.2040"]
    bc.an$P.Value.Upper.Tail[bc.an$ID == "11AEBGC.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "11AEBGC.1040"], bc.an$Df[bc.an$ID == "11AEBGC.1040"], ac.an$Df[ac.an$ID == "11AEBGC.2040"], lower.tail = F)
    bc.an$P.Value.Lower.Tail[bc.an$ID == "11AEBGC.1040"]<- stats::pf(bc.an$F.Value[bc.an$ID == "11AEBGC.1040"], bc.an$Df[bc.an$ID == "11AEBGC.1040"], ac.an$Df[ac.an$ID == "11AEBGC.2040"], lower.tail = T)
  }


  ###---------------------------- IMPACTs (BEFORE) = From sym.B - sym.CB
  ##--------------------- codigo xx.1030

  bi.an<-baci[F,]
  ###NO IF CHAMA A COMBINACAO LM ORIGINAL, MAS NA ANOVA PARCIAL ESTA SEM O "A"

  if("AEB" %in% aav.model.lm == T) {
    di08AEB.1030<-cbind(Source.of.Variation = "B:E",
                        sym.B[sym.B$Source.of.Variation == "B:E",2:3] -
                          sym.CB[sym.CB$Source.of.Variation == "B:E",2:3])
    bi.an<-dplyr::bind_rows(bi.an,di08AEB.1030)
    bi.an$ID[bi.an$Source.of.Variation == "B:E"] <- "08AEB.1030"
    bi.an$Calculated.From[bi.an$ID == "08AEB.1030"]<-"sym.B-sym.CB"
    bi.an$Df.Formula[bi.an$ID == "08AEB.1030"]<-"(p-1)"
    bi.an$Mean.Sq[bi.an$ID == "08AEB.1030"]<-bi.an$Sum.Sq[bi.an$ID == "08AEB.1030"] / bi.an$Df[bi.an$ID == "08AEB.1030"]
    bi.an$Source.of.Variation[bi.an$ID == "08AEB.1030"] <-"P(Bef) x I"
    # 2-tail
    # 	Two-tailed test F = (MS_P(Bef) xIC / MS_P(Aft) xIC)
    bi.an$Denominator[bi.an$ID == "08AEB.1030"]<-"P(Aft)xI"
    bi.an$F.Value[bi.an$ID == "08AEB.1030"]<-bi.an$Mean.Sq[bi.an$ID == "08AEB.1030"] / ai.an$Mean.Sq[ai.an$ID == "08AEB.2030"]
    bi.an$P.Value.Upper.Tail[bi.an$ID == "08AEB.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "08AEB.1030"], bi.an$Df[bi.an$ID == "08AEB.1030"], ai.an$Df[ai.an$ID == "08AEB.2030"], lower.tail = F)
    bi.an$P.Value.Lower.Tail[bi.an$ID == "08AEB.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "08AEB.1030"], bi.an$Df[bi.an$ID == "08AEB.1030"], ai.an$Df[ai.an$ID == "08AEB.2030"], lower.tail = T)
  }
  if("AEBG" %in% aav.model.lm == T) {
    di09AEBG.1030<-cbind(Source.of.Variation = "B:E:G",
                         sym.B[sym.B$Source.of.Variation == "B:E:G",2:3] -
                           sym.CB[sym.CB$Source.of.Variation == "B:E:G",2:3])
    bi.an<-dplyr::bind_rows(bi.an,di09AEBG.1030)
    bi.an$ID[bi.an$Source.of.Variation == "B:E:G"] <- "09AEBG.1030"
    bi.an$Calculated.From[bi.an$ID == "09AEBG.1030"]<-"sym.B-sym.CB"
    bi.an$Df.Formula[bi.an$ID == "09AEBG.1030"]<-"(p-1)(s-1)"
    bi.an$Mean.Sq[bi.an$ID == "09AEBG.1030"]<- bi.an$Sum.Sq[bi.an$ID == "09AEBG.1030"] /  bi.an$Df[bi.an$ID == "09AEBG.1030"]
    bi.an$Source.of.Variation[bi.an$ID == "09AEBG.1030"] <-"P(Bef) x S(I)"
    # 2-tail
    # Two-tailed test F = (MS_P(Bef) xS(IC) / MS_P(Aft) xS(IC))
    bi.an$Denominator[bi.an$ID == "09AEBG.1030"]<-"P(Aft)xS(I)"
    bi.an$F.Value[bi.an$ID == "09AEBG.1030"]<-bi.an$Mean.Sq[bi.an$ID == "09AEBG.1030"] / ai.an$Mean.Sq[ai.an$ID == "09AEBG.2030"]
    bi.an$P.Value.Upper.Tail[bi.an$ID == "09AEBG.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "09AEBG.1030"], bi.an$Df[bi.an$ID == "09AEBG.1030"], ai.an$Df[ai.an$ID == "09AEBG.2030"], lower.tail = F)
    bi.an$P.Value.Lower.Tail[bi.an$ID == "09AEBG.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "09AEBG.1030"], bi.an$Df[bi.an$ID == "09AEBG.1030"], ai.an$Df[ai.an$ID == "09AEBG.2030"], lower.tail = T)
  }
  if("AEBC" %in% aav.model.lm == T) {
    di10AEBC.1030<-cbind(Source.of.Variation = "B:E:C",
                         sym.B[sym.B$Source.of.Variation == "B:E:C",2:3] -
                           sym.CB[sym.CB$Source.of.Variation == "B:E:C",2:3])
    bi.an<-dplyr::bind_rows(bi.an,di10AEBC.1030)
    bi.an$ID[bi.an$Source.of.Variation == "B:E:C"] <- "10AEBC.1030"
    bi.an$Calculated.From[bi.an$ID == "10AEBC.1030"]<-"sym.B-sym.CB"
    bi.an$Df.Formula[bi.an$ID == "10AEBC.1030"]<-"p(t-1)"
    bi.an$Mean.Sq[bi.an$ID == "10AEBC.1030"]<- bi.an$Sum.Sq[bi.an$ID == "10AEBC.1030"] / bi.an$Df[bi.an$ID == "10AEBC.1030"]
    bi.an$Source.of.Variation[bi.an$ID == "10AEBC.1030"] <-"T(P(Bef) x I"
    # 2-tail
    # 	Two-tailed test F = (MS_T(P(Bef)) xIC / MS_T(P(Aft)) xIC)
    bi.an$Denominator[bi.an$ID == "10AEBC.1030"]<-"T(P(Aft)xI"
    bi.an$F.Value[bi.an$ID == "10AEBC.1030"]<- bi.an$Mean.Sq[bi.an$ID == "10AEBC.1030"] / ai.an$Mean.Sq[ai.an$ID == "10AEBC.2030"]
    bi.an$P.Value.Upper.Tail[bi.an$ID == "10AEBC.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "10AEBC.1030"], bi.an$Df[bi.an$ID == "10AEBC.1030"], ai.an$Df[ai.an$ID == "10AEBC.2030"], lower.tail = F)
    bi.an$P.Value.Lower.Tail[bi.an$ID == "10AEBC.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "10AEBC.1030"], bi.an$Df[bi.an$ID == "10AEBC.1030"], ai.an$Df[ai.an$ID == "10AEBC.2030"], lower.tail = T)
  }
  if("AEC" %in% aav.model.lm == T) {
    di10AEC.1030<-cbind(Source.of.Variation = "C:E",
                        sym.B[sym.B$Source.of.Variation == "C:E",2:3] -
                          sym.CB[sym.CB$Source.of.Variation == "C:E",2:3])
    bi.an<-dplyr::bind_rows(bi.an,di10AEC.1030)
    bi.an$ID[bi.an$Source.of.Variation == "C:E"] <- "10AEC.1030"
    bi.an$Calculated.From[bi.an$ID == "10AEC.1030"]<-"sym.B-sym.CB"
    bi.an$Df.Formula[bi.an$ID == "10AEC.1030"]<-"(t-1)"
    bi.an$Mean.Sq[bi.an$ID == "10AEC.1030"]<- bi.an$Sum.Sq[bi.an$ID == "10AEC.1030"] /  bi.an$Df[bi.an$ID == "10AEC.1030"]
    bi.an$Source.of.Variation[bi.an$ID == "10AEC.1030"] <-"T(Bef) x I"
    # 2-tail
    # 	Two-tailed test F = (MS_T(Bef) xIC / MS_T(Aft) xIC)
    bi.an$Denominator[bi.an$ID == "10AEC.1030"]<-"T(Aft)xI"
    bi.an$F.Value[bi.an$ID == "10AEC.1030"]<-bi.an$Mean.Sq[bi.an$ID == "10AEC.1030"] / ai.an$Mean.Sq[ai.an$ID == "10AEC.2030"]
    bi.an$P.Value.Upper.Tail[bi.an$ID == "10AEC.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "10AEC.1030"], bi.an$Df[bi.an$ID == "10AEC.1030"], ai.an$Df[ai.an$ID == "10AEC.2030"], lower.tail = F)
    bi.an$P.Value.Lower.Tail[bi.an$ID == "10AEC.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "10AEC.1030"], bi.an$Df[bi.an$ID == "10AEC.1030"], ai.an$Df[ai.an$ID == "10AEC.2030"], lower.tail = T)
  }
  if("AEGC" %in% aav.model.lm == T) {
    di11AEGC.1030<-cbind(Source.of.Variation = "C:E:G",
                         sym.B[sym.B$Source.of.Variation == "C:E:G",2:3] -
                           sym.CB[sym.CB$Source.of.Variation == "C:E:G",2:3])
    bi.an<-dplyr::bind_rows(bi.an,di11AEGC.1030)
    bi.an$ID[bi.an$Source.of.Variation == "C:E:G"] <- "11AEGC.1030"
    bi.an$Calculated.From[bi.an$ID == "11AEGC.1030"]<-"sym.B-sym.CB"
    bi.an$Df.Formula[bi.an$ID == "11AEGC.1030"]<-"(t-1)(s-1)"
    bi.an$Mean.Sq[bi.an$ID == "11AEGC.1030"]<- bi.an$Sum.Sq[bi.an$ID == "11AEGC.1030"] /  bi.an$Df[bi.an$ID == "11AEGC.1030"]
    bi.an$Source.of.Variation[bi.an$ID == "11AEGC.1030"] <-"T(Bef) x S(I)"
    # 2-tail
    # 	Two-tailed test F = (MS_T(Bef) xS(IC) / MS_T(Aft) xS(IC))
    bi.an$Denominator[bi.an$ID == "11AEGC.1030"]<-"T(Aft)xS(I)"
    bi.an$F.Value[bi.an$ID == "11AEGC.1030"]<-bi.an$Mean.Sq[bi.an$ID == "11AEGC.1030"] / ai.an$Mean.Sq[ai.an$ID == "11AEGC.2030"]
    bi.an$P.Value.Upper.Tail[bi.an$ID == "11AEGC.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "11AEGC.1030"], bi.an$Df[bi.an$ID == "11AEGC.1030"], ai.an$Df[ai.an$ID == "11AEGC.2030"], lower.tail = F)
    bi.an$P.Value.Lower.Tail[bi.an$ID == "11AEGC.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "11AEGC.1030"], bi.an$Df[bi.an$ID == "11AEGC.1030"], ai.an$Df[ai.an$ID == "11AEGC.2030"], lower.tail = T)
  }
  if("AEBGC" %in% aav.model.lm == T) {
    di11AEBGC.1030<-cbind(Source.of.Variation = "B:E:G:C",
                          sym.B[sym.B$Source.of.Variation == "B:E:G:C",2:3] -
                            sym.CB[sym.CB$Source.of.Variation == "B:E:G:C",2:3])
    bi.an<-dplyr::bind_rows(bi.an,di11AEBGC.1030)
    bi.an$ID[bi.an$Source.of.Variation == "B:E:G:C"] <- "11AEBGC.1030"
    bi.an$Calculated.From[bi.an$ID == "11AEBGC.1030"]<-"sym.B-sym.CB"
    bi.an$Df.Formula[bi.an$ID == "11AEBGC.1030"]<-"p(t-1)(s-1)"
    bi.an$Mean.Sq[bi.an$ID == "11AEBGC.1030"]<- bi.an$Sum.Sq[bi.an$ID == "11AEBGC.1030"] /  bi.an$Df[bi.an$ID == "11AEBGC.1030"]
    bi.an$Source.of.Variation[bi.an$ID == "11AEBGC.1030"] <-"T(P(Bef)) x S(I)"
    # 2-tail
    # Two-tailed test F = (MS T(P(Bef)) xS(IC) / MS_T(P(Aft)) xS(IC))
    bi.an$Denominator[bi.an$ID == "11AEBGC.1030"]<-"T(P(Aft))xS(I)"
    bi.an$F.Value[bi.an$ID == "11AEBGC.1030"]<-bi.an$Mean.Sq[bi.an$ID == "11AEBGC.1030"] / ai.an$Mean.Sq[ai.an$ID == "11AEBGC.2030"]
    bi.an$P.Value.Upper.Tail[bi.an$ID == "11AEBGC.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "11AEBGC.1030"], bi.an$Df[bi.an$ID == "11AEBGC.1030"], ai.an$Df[ai.an$ID == "11AEBGC.2030"], lower.tail = F)
    bi.an$P.Value.Lower.Tail[bi.an$ID == "11AEBGC.1030"]<- stats::pf(bi.an$F.Value[bi.an$ID == "11AEBGC.1030"], bi.an$Df[bi.an$ID == "11AEBGC.1030"], ai.an$Df[ai.an$ID == "11AEBGC.2030"], lower.tail = T)
  }


  ###----------------------------------------------------------------------- End

  ###----------------------------------------- #################################
  ###----------------------------------------- ### Table: Asymmetrical ANOVA ###
  ###----------------------------------------- #################################

  baci<-rbind(f.an,c.an,i.an,b.an,a.an,bc.an,bi.an,ac.an,ai.an)
  baci<-dplyr::arrange(baci,ID)

  ###------------------------------------------------------------- Significant
  ## The F-ratio for some terms of the aav.models AEBGC, AEGC
  ## should be made by pooling procedures. Otherwise there is no test!

  # agora considerando P/2 para os casos de 2-tail
  # * se um deles for menor que P.avv/2; * se o unico for menor que P.avv;
  # ns se o unico for menor que p.pooling; eliminate se o unico for maior que p.pooling;
  # ns para todos os outros (se os dois forem > P.avv/2)
  if(pooling != T | (pooling == T & p.aav < p.pooling)) {
    baci$Significant<-ifelse(baci$P.Value.Upper.Tail <= (p.aav/2) & !is.na(baci$P.Value.Upper.Tail) |
           baci$P.Value.Lower.Tail <= (p.aav/2) & !is.na(baci$P.Value.Lower.Tail), "*",
         ifelse(baci$P.Value.Upper.Tail <= p.aav & !is.na(baci$P.Value.Upper.Tail) & is.na(baci$P.Value.Lower.Tail), "*",
                ifelse(baci$P.Value.Upper.Tail <= p.pooling & !is.na(baci$P.Value.Upper.Tail) & is.na(baci$P.Value.Lower.Tail), "ns",
                       ifelse(baci$P.Value.Upper.Tail > p.pooling & !is.na(baci$P.Value.Upper.Tail) & is.na(baci$P.Value.Lower.Tail),"ns (Can be Eliminated)",
                              "ns"))))

#    baci$Significant<-ifelse(baci$P.Value.Upper.Tail <= p.aav & !is.na(baci$P.Value.Upper.Tail) |
#                               baci$P.Value.Lower.Tail <= p.aav & !is.na(baci$P.Value.Lower.Tail), "*",
#                             ifelse(baci$P.Value.Upper.Tail < p.pooling & !is.na(baci$P.Value.Upper.Tail) |
#                                      baci$P.Value.Lower.Tail < p.pooling & !is.na(baci$P.Value.Lower.Tail), "ns",
#                                    ifelse(baci$P.Value.Upper.Tail >= p.pooling & !is.na(baci$P.Value.Upper.Tail) |
#                                             baci$P.Value.Lower.Tail >= p.pooling & !is.na(baci$P.Value.Lower.Tail),"ns (Can be Eliminated)", NA)))
  }
  if(pooling == T) {
    if(p.aav > p.pooling)   stop("p.aav must be equal or less than p.pooling")
    if(p.aav == p.pooling) {
      baci$Significant<-ifelse(baci$P.Value.Upper.Tail <= (p.aav/2) & !is.na(baci$P.Value.Upper.Tail) |
               baci$P.Value.Lower.Tail <= (p.aav/2) & !is.na(baci$P.Value.Lower.Tail), "*",
             ifelse(baci$P.Value.Upper.Tail <= p.aav & !is.na(baci$P.Value.Upper.Tail) & is.na(baci$P.Value.Lower.Tail), "*",
                    ifelse(baci$P.Value.Upper.Tail > p.pooling & !is.na(baci$P.Value.Upper.Tail) & is.na(baci$P.Value.Lower.Tail),"ns (Can be Eliminated)",
                           "ns")))

#      baci$Significant<-ifelse(baci$P.Value.Upper.Tail <= p.aav & !is.na(baci$P.Value.Upper.Tail) |
#                                 baci$P.Value.Lower.Tail <= p.aav & !is.na(baci$P.Value.Lower.Tail), "*",
#                               ifelse(baci$P.Value.Upper.Tail > p.pooling & !is.na(baci$P.Value.Upper.Tail) |
#                                        baci$P.Value.Lower.Tail > p.pooling & !is.na(baci$P.Value.Lower.Tail),"ns (Can be Eliminated)", NA))
    }
  }

  ###-------------------------------------------------------- Post.Hoc.Pooling

  if(aav.model == "AEC") {
    # 	F = (MS_T(Aft) xCC / MS_Res)
    baci$Post.Hoc.Pooling[baci$ID == "10AEC.2040"]<-"Cannot be pooled"

    # 	F = (MSB xCC / MST(Aft) xCC)
    baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]<-"Cannot be pooled"

    if(pooling != T) {
      # 	F = (MST(Aft) xIC / MST(Aft) xCC)
      baci$Post.Hoc.Pooling[baci$ID == "10AEC.2030"]<-"Cannot be pooled"

      # 	F = (MSB xIC / MSB xCC), or
      baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"
    }

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	F = (MST(Aft) xIC / MST(Aft) xCC)
      if(baci$P.Value.Upper.Tail[baci$ID == "10AEC.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEC.2030"]<-"Cannot be pooled"

      # 	If varT(Aft) xCC <= 0, then F = (MST(Aft) xIC / MSpooled), where MSpooled = (SST(Aft) xCC + SSRes) / (dfT(Aft) xCC + dfRes)
      if(baci$P.Value.Upper.Tail[baci$ID == "10AEC.2040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEC.2031"]<-"Can be pooled"
      AEC10<- c("10AEC.2030","10AEC.2031")
      AEC10n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEC.2030"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEC.2031"])) == F)

      if("eliminate" %in% method == T) {
        # 	F = (MSB xCC / MSpooled), where MSpooled = (SST(Aft) xCC + SSRes) / (dfT(Aft) xCC + dfRes)
        baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]<-NA
        baci$Post.Hoc.Pooling[baci$ID == "03AE.0041"]<-"Can be pooled"

        # 	If varB xCC <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SSRes) / (dfB xCC + dfRes)
        baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Can be pooled"
      }
      if("pool" %in% method == T) {
        # 	F = (MSB xIC / MSB xCC), or
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"

        # 	If varB xCC <= 0 and varT(Aft) xIC > 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SST(Aft) xIC) / (dfB xCC + dfT(Aft) xIC) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Can be pooled"

        # 	If varB xCC <= 0 and varT(Aft) xIC <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0032"]<-"Can be pooled"
      }
    }
  }



  ###---------------- INTERACTIONS with THREE temporal and ONE spatial factors

  if(aav.model == "AEBC") {
    # 	F = (MS_T(P(Aft)) xCC / MS_Res)
    baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2040"]<-"Cannot be pooled"

    # 	F = (MS_P(Aft) xCC / MS_T(P(Aft)) xCC), or
    baci$Post.Hoc.Pooling[baci$ID == "08AEB.2040"]<-"Cannot be pooled"

    # 	F = (MS_B xCC / MS_P(Aft) xCC), or
    baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]<-"Cannot be pooled"

    if(pooling != T) {
      # 	F = (MS_T(P(Aft)) xIC / MS_T(P(Aft)) xCC), or
      baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2030"]<-"Cannot be pooled"

      # 	F = (MS_P(Aft) xIC / MS_P(Aft) xCC), or
      baci$Post.Hoc.Pooling[baci$ID == "08AEB.2030"]<-"Cannot be pooled"

      #F = (MSB xIC / MSB xCC), or
      baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"
    }

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	F = (MS_T(P(Aft)) xIC / MS_T(P(Aft)) xCC), or
      if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2030"]<-"Cannot be pooled"

      # 	If varT(P(Aft)) xCC <= 0, then F = (MS_T(P(Aft)) xIC / MSpooled), where MSpooled = (SST(P(Aft)) xCC + SSRes) / (dfT(P(Aft)) xCC + dfRes)
      if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2031"]<-"Can be pooled"
      AEC10<- c("10AEBC.2030","10AEBC.2031")
      AEC10n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2030"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2031"])) == F)

      if("eliminate" %in% method == T) {
        # 	F = (MS_P(Aft) xCC / MS_pooled), where MS_pooled = (SS_T(P(Aft)) xCC + SS_Res) / (df_T(P(Aft)) xCC + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "08AEB.2040"]<-NA
        baci$Post.Hoc.Pooling[baci$ID == "08AEB.2041"]<-"Can be pooled"

        # 	If varP(Aft) xCC <= 0, then F = (MS_P(Aft) xIC / MS_pooled), where MS_pooled = (SS_P(Aft) xCC + SS_Res) / (df_P(Aft) xCC + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "08AEB.2030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == "08AEB.2041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2030"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == "08AEB.2041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2031"]<-"Can be pooled"

        # 	F = (MS_B xCC / MS_pooled), where MS_pooled = (SS_P(Aft) xCC + SS_Res) / (df_P(Aft) xCC + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]<-NA
        baci$Post.Hoc.Pooling[baci$ID == "03AE.0041"]<-"Can be pooled"

        # 	If varB xCC <= 0, then F = (MS_B xIC / MSpooled), where MSpooled = (SSB xCC + SSRes) / (dfB xCC + dfRes)
        baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Can be pooled"
      }

      if("pool" %in% method == T) {
        # 	F = (MS_P(Aft) xIC / MS_P(Aft) xCC), or
        if(baci$P.Value.Upper.Tail[baci$ID == "08AEB.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2030"]<-"Cannot be pooled"
        # 	If varP(Aft) xCC <= 0 and varT(P(Aft)) xIC > 0, then F = (MSP(Aft) xIC / MSpooled), where MSpooled = (SSP(Aft) xCC + SST(P(Aft)) xIC) / (dfP(Aft) xCC + dfT(P(Aft)) xIC) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == "08AEB.2040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2031"]<-"Can be pooled"
        # 	If varP(Aft) xCC <= 0 and varT(P(Aft)) xIC <= 0, then F = (MSP(Aft) xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "08AEB.2040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2032"]<-"Can be pooled"
        AEB08<- c("08AEB.2030","08AEB.2031","08AEB.2032")
        AEB08n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2030"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2031"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2032"])) == F)

        #F = (MSB xIC / MSB xCC), or
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"
        # 	If varB xCC <= 0 and varP(Aft) xIC > 0, then F = (MS_B xIC / MS_pooled),
        # where MSpooled = (SS_B xCC + SS_P(Aft) xIC) / (df_B xCC + df_P(Aft) xIC) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Can be pooled"
        # 	If varB xCC <= 0, varP(Aft) xIC <= 0 and varT(P(Aft)) xIC > 0, then F = (MS_B xIC / MS_pooled),
        # where MS_pooled = (SS1 + SS_T(P(Aft)) xIC) / (df1 + df_T(P(Aft)) xIC) = SS2/df2
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0032"]<-"Can be pooled"
        # 	If varB xCC <= 0, varP(Aft) xIC <= 0 and varT(P(Aft)) xIC <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SS2 + SSRes) / (df2 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "03AE.0040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0033"]<-"Can be pooled"
      }
    }
  }


  if(aav.model == "AEGC") {
    #	F = (MS_T(Aft) xS(CC) / MS_Res)
    baci$Post.Hoc.Pooling[baci$ID == "11AEGC.2040"]<-"Cannot be pooled"

    #  	F = (MST(Aft) xCC / MS T(Aft) xS(CC)), or
    baci$Post.Hoc.Pooling[baci$ID == "10AEC.2040"]<-"Cannot be pooled"

    #  	F = (MSB xS(CC) / MST(Aft) xS(CC)), or
    baci$Post.Hoc.Pooling[baci$ID == "07AEG.0040"]<-"Cannot be pooled"

    # If the two terms cannot be eliminated (* or ns)
    if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] < p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "10AEC.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]<-"Cannot be pooled"
    # 	If varB xS(CC) <= 0 and varT(Aft) xCC > 0, then F = (MSB xCC / MST(Aft) xCC)
    if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "10AEC.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]<-"Can be pooled"
    # 	If varT(Aft) xCC <= 0 and varB xS(CC) > 0, then F = (MSB xCC / MSB xS(CC))
    if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] < p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "10AEC.2040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0041"]<-"Can be pooled"
    # 	If varT(Aft) xCC <= 0, varB xS(CC) <= 0, then F = (MSB xCC / MSpooled),
    # where MSpooled = (SS_T(Aft) xCC + SS_B xS(CC) + SS_T(Aft) xS(CC)) / (df_T(Aft) xCC + df_B xS(CC) + df_T(Aft) xS(CC)), or
    if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "10AEC.2040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0042"]<-"Can be pooled"
    AE03<- c("03AE.0040","03AE.0041","03AE.0042")
    AE03n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]),
                   is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0041"]),
                   is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0042"])) == F)

    if(pooling != T) {
      #  	F = (MST(Aft) xS(IC) / MST(Aft) xS(CC)), or
      baci$Post.Hoc.Pooling[baci$ID == "11AEGC.2030"]<-"Cannot be pooled"

      #F = (MST(Aft) xIC / MST(Aft) xCC), or
      baci$Post.Hoc.Pooling[baci$ID == "10AEC.2030"]<-"Cannot be pooled"

      # 	F = (MSB xS(IC) / MSB xS(CC)), or
      baci$Post.Hoc.Pooling[baci$ID == "07AEG.0030"]<-"Cannot be pooled"

      # 	F = (MSB xIC / MSB xCC), or
      baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"
    }

    ### POOLING STRUCTURE
    if(pooling == T) {
      #  	F = (MST(Aft) xS(IC) / MST(Aft) xS(CC)), or
      if(baci$P.Value.Upper.Tail[baci$ID == "11AEGC.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "11AEGC.2030"]<-"Cannot be pooled"
      #  	If varT(Aft) xS(CC) <= 0, then F = (MST(Aft) xS(IC) / MSpooled), where MSpooled = (SST(Aft) xS(CC) + SSRes) / (dfT(Aft) xS(CC) + dfRes)
      if(baci$P.Value.Upper.Tail[baci$ID == "11AEGC.2040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "11AEGC.2031"]<-"Can be pooled"
      AEGC11<- c("11AEGC.2030","11AEGC.2031")
      AEGC11n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "11AEGC.2030"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "11AEGC.2031"])) == F)

      if("eliminate" %in% method == T) {
        # 	F = (MST(Aft) xCC / MSpooled), where MSpooled = (SST(Aft) xS(CC) + SSRes) / (dfT(Aft) xS(CC) + dfRes)
        baci$Post.Hoc.Pooling[baci$ID == "10AEC.2040"]<-NA
        baci$Post.Hoc.Pooling[baci$ID == "10AEC.2041"]<-"Can be pooled"

        # 	If varT(Aft) xCC <= 0, then F = (MST(Aft) xIC / MSpooled), where MSpooled = (SST(Aft) xCC + SSRes) / (dfT(Aft) xCC + dfRes)
        baci$Post.Hoc.Pooling[baci$ID == "10AEC.2030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEC.2041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEC.2031"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEC.2041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEC.2031"]<-"Can be pooled"

        # 	F = (MSB xS(CC) / MSpooled), where MSpooled = (SST(Aft) xS(CC) + SSRes) / (dfT(Aft) xS(CC) + dfRes)
        baci$Post.Hoc.Pooling[baci$ID == "07AEG.0040"]<-NA
        baci$Post.Hoc.Pooling[baci$ID == "07AEG.0041"]<-"Can be pooled"

        # 	If varB xS(CC) <= 0, then F = (MSB xS(IC) / MSpooled), where MSpooled = (SSB xS(CC) + SSRes) / (dfB xS(CC) + dfRes)
        baci$Post.Hoc.Pooling[baci$ID == "07AEG.0030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0031"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0031"]<-"Can be pooled"

        baci$Post.Hoc.Pooling[baci$ID == AE03[AE03n]]<-NA
        # If the two terms cannot be eliminated (* or ns)
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "10AEC.2041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == AE03[AE03n]]<-"Can be pooled"
        # 	If varB xS(CC) <= 0 and varT(Aft) xCC > 0, then F = (MSB xCC / MSpooled), where MSpooled = (SST(Aft) xCC + SSRes) / (dfT(Aft) xCC  + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "10AEC.2041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0043"]<-"Can be pooled"
        # 	If varT(Aft) xCC <= 0 and varB xS(CC) > 0, then F = (MSB xCC / MSpooled), where MSpooled = (SSB xS(CC) + SSRes) / (dfB xS(CC) + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "10AEC.2041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0044"]<-"Can be pooled"
        # 	If varT(Aft) xCC <= 0, varB xS(CC) <= 0, then F = (MSB xCC / MSpooled),
        # where MSpooled = (SST(Aft) xCC + SSB xS(CC) + SST(Aft) xS(CC) + SSRes) / (dfT(Aft) xCC + dfB xS(CC) + dfT(Aft) xS(CC) + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "10AEC.2041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0045"]<-"Can be pooled"
        AE03e<- c("03AE.0040","03AE.0041","03AE.0042","03AE.0043","03AE.0044","03AE.0045")
        AE03en<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0041"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0042"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0043"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0044"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0045"])) == F)

        # 	If varB xCC <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SSRes) / (dfB xCC + dfRes)
        baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == AE03e[AE03en]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == AE03e[AE03en]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Can be pooled"
      }

      if("pool" %in% method == T) {
        #F = (MST(Aft) xIC / MST(Aft) xCC), or
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEC.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEC.2030"]<-"Cannot be pooled"
        # 	If varT(Aft) xCC <= 0 and varT(Aft) xS(IC) > 0, then F = (MST(Aft) xIC / MSpooled), where MSpooled = (SST(Aft) xCC + SST(Aft) xS(IC) / (dfT(Aft) xCC + dfT(Aft) xS(IC) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEC.2040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEGC11[AEGC11n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEC.2031"]<-"Can be pooled"
        # 	If varT(Aft) xCC <= 0 and varT(Aft) xS(IC) <= 0, then F = (MST(Aft) xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEC.2040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEGC11[AEGC11n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEC.2032"]<-"Can be pooled"
        AEC10<- c("10AEC.2030","10AEC.2031","10AEC.2032")
        AEC10n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEC.2030"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEC.2031"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEC.2032"])) == F)

        # 	F = (MSB xS(IC) / MSB xS(CC)), or
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0030"]<-"Cannot be pooled"
        # 	If varB xS(CC) <= 0 and varT(Aft) xS(IC) > 0, then F = (MSB xS(IC) / MSpooled),
        # where MSpooled = (SSB xS(CC) + SST(Aft) xS(IC)) / (dfB xS(CC) + dfT(Aft) xS(IC)) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEGC11[AEGC11n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0031"]<-"Can be pooled"
        # 	If varB xS(CC) <= 0 and varT(Aft) xS(IC) <= 0, then F = (MSB xS(IC) / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEGC11[AEGC11n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0032"]<-"Can be pooled"
        AEG07<- c("07AEG.0030","07AEG.0031","07AEG.0032")
        AEG07n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0030"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0031"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0032"])) == F)

        # If the two terms cannot be eliminated (* or ns)
        # 	F = (MSB xIC / MSB xCC), or
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]]< p.pooling |
           (baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
            baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] < p.pooling &
            baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] < p.pooling)) baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"
        # 	If varB xCC <= 0, varB xS(IC) <= 0 and varT(Aft) xIC > 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SST(Aft) xIC) / (dfB xCC + dfT(Aft) xIC)
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Can be pooled"
        # 	If varB xCC <= 0, varT(Aft) xIC <= 0 and varB xS(IC) > 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SSB xS(IC)) / (dfB xCC + dfB xS(IC))
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0032"]<-"Can be pooled"
        # 	If varB xCC <= 0, varT(Aft) xIC <= 0, varB xS(IC) <= 0 and varT(Aft) xS(IC) > 0, then F = (MSB xIC / MSpooled), where MSpooled = (SSB xCC + SST(Aft) xIC + SSB xS(IC) + SST(Aft) xS(IC)) / (dfB xCC + dfT(Aft) xIC + dfB xS(IC) + dfT(Aft) xS(IC)) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEGC11[AEGC11n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0033"]<-"Can be pooled"
        # 	If varB xCC <= 0, varB xS(IC) <= 0, varT(B) xIC <= 0, varT(B) xS(IC) <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEC10[AEC10n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEGC11[AEGC11n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0034"]<-"Can be pooled"
      }
    }
  }


  if(aav.model == "AEBGC") {
    # 	F = (MS_T(P(Aft)) xS(CC) / MS_Res)
    baci$Post.Hoc.Pooling[baci$ID == "11AEBGC.2040"]<-"Cannot be pooled"

    # 	F = (MS_T(P(Aft)) xCC / MS_T(P(Aft)) xS(CC)), or
    baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2040"]<-"Cannot be pooled"

    # 	F = (MS_P(Aft) xS(CC) / MS_T(P(Aft)) xS(CC)), or
    baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2040"]<-"Cannot be pooled"

    # If the two terms cannot be eliminated (* or ns)
    if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] < p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2040"]<-"Cannot be pooled"
    # 	If varT(P(Aft)) xCC <= and varP(Aft) xS(CC) > 0, then F = (MS_P(Aft) xCC / MS_P(Aft) xS(CC))
    if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2040"]<-"Can be pooled"
    # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xCC > 0, then F = (MS_P(Aft) xCC / MS_T(P(Aft)) xCC)
    if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] < p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2041"]<-"Can be pooled"
    # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xCC <= 0, then F = (MS_P(Aft) xCC / MS_pooled), where MS_pooled = (SS_P(Aft) xS(CC) + SS_T(P(Aft)) xCC + SS_T(P(Aft)) xS(CC)) / (df_P(Aft) xS(CC) + df_T(P(Aft)) xCC + df_T(P(Aft)) xS(CC))
    if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2042"]<-"Can be pooled"
    AEB08<- c("08AEB.2040","08AEB.2041","08AEB.2042")
    AEB08n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2040"]),
                    is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2041"]),
                    is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2042"])) == F)

    #F = (MS_B xS(CC) / MS_P(Aft) xS(CC)), or
    baci$Post.Hoc.Pooling[baci$ID == "07AEG.0040"]<-"Cannot be pooled"

    # If the two terms cannot be eliminated (* or ns)
    if((baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] < p.pooling &
        baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] < p.pooling) |
       (baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] < p.pooling &
        baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] < p.pooling)) baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]<-"Cannot be pooled"
    # 	If varP(Aft) xCC <= 0 and varB xS(CC) > 0 then F = (MS_B xCC / MS_B xS(CC))
    if(baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]<-"Can be pooled"
    # 	If varB xS(CC) <= 0 and varP(Aft) xCC > 0, then F = (MS_B xCC / MS_P(Aft) xCC)
    if(baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] < p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0041"]<-"Can be pooled"
    # 	If varB xS(CC) <= 0, varP(Aft) xCC <= 0, varT(P(Aft)) xCC <= 0 and varP(Aft) xS(CC) > 0,
    # then F = (MS_B xCC / MS_pooled), where MS_pooled = (SS_B xS(CC) + SS_P(Aft) xCC + SS_P(Aft) xS(CC)) / (df_B xS(CC) + df_P(Aft) xCC + df_P(Aft) xS(CC))
    if(baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0042"]<-"Can be pooled"
    # 	If varB xS(CC) <= 0, varP(Aft) xCC <= 0, varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xCC > 0,
    #then F = (MSB xCC / MSpooled), where MSpooled = (SS_B xS(CC) + SS_P(Aft) xCC + SS_T(P(Aft)) xCC) / (df_B xS(CC) + df_P(Aft) xCC + df_T(P(Aft)) xCC)
    if(baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] < p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0043"]<-"Can be pooled"
    # 	If varB xS(CC) <= 0, varP(Aft) xCC <= 0, varP(Aft) xS(CC) <= 0, varT(P(Aft)) xCC <= 0,
    #then F = (MSB xCC / MSpooled), where MSpooled = (SS_B xS(CC) + S_SP(Aft) xCC + SS_P(Aft) xS(CC) + SS_T(P(Aft)) xCC + SS_T(P(Aft)) xS(CC)) / (df_B xS(CC) + df_P(Aft) xCC + df_P(Aft) xS(CC) + df_T(P(Aft)) xCC + df_T(P(Aft)) xS(CC))
    if(baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] >= p.pooling &
       baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0044"]<-"Can be pooled"
    AE03<- c("03AE.0040","03AE.0041","03AE.0042","03AE.0043","03AE.0044")
    AE03n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]),
                   is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0041"]),
                   is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0042"]),
                   is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0043"]),
                   is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0044"])) == F)

    if(pooling != T) {
      # 	F = (MS_T(P(Aft)) xS(IC) / MS_T(P(Aft)) xS(CC)), or
      baci$Post.Hoc.Pooling[baci$ID == "11AEBGC.2030"]<-"Cannot be pooled"

      # 	F = (MS_T(P(Aft)) xIC / MS_T(P(Aft)) xCC), or
      baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2030"]<-"Cannot be pooled"

      # 	F = (MS_P(Aft) xS(IC) / MS_P(Aft) xS(CC)), or
      baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2030"]<-"Cannot be pooled"

      # 	F = (MS_P(Aft) xIC / MS_P(Aft) xCC)
      baci$Post.Hoc.Pooling[baci$ID == "08AEB.2030"]<-"Cannot be pooled"

      # 	F = (MS_B xS(IC) / MS_B xS(CC)), or
      baci$Post.Hoc.Pooling[baci$ID == "07AEG.0030"]<-"Cannot be pooled"

      # 	F = (MS_B xIC / MS_B xCC), or
      baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"
    }

    ### POOLING STRUCTURE
    if(pooling == T) {
      # 	F = (MS_T(P(Aft)) xS(IC) / MS_T(P(Aft)) xS(CC)), or
      if(baci$P.Value.Upper.Tail[baci$ID == "11AEBGC.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "11AEBGC.2030"]<-"Cannot be pooled"
      # 	If varT(P(Aft)) xS(CC) <= 0, then F = (MS_T(P(Aft)) xS(IC) / MS_pooled), where MS_pooled = (SS_T(P(Aft)) xS(CC) + SS_Res) / (df_T(P(Aft)) xS(CC) + df_Res)
      if(baci$P.Value.Upper.Tail[baci$ID == "11AEBGC.2040"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "11AEBGC.2031"]<-"Can be pooled"
      AEBGC11<- c("11AEBGC.2030","11AEBGC.2031")
      AEBGC11n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "11AEBGC.2030"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "11AEBGC.2031"])) == F)

      if(method == "eliminate") {
        # 	F = (MS_T(P(Aft)) xCC / MS_pooled), where MS_pooled = (SS_T(P(Aft)) xS(CC)) + SS_Res) / (df_T(P(Aft)) xS(CC)) + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2040"]<-NA
        baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2041"]<-"Can be pooled"

        # 	If varT(P(Aft)) xCC <= 0, then F = (MS_T(P(Aft)) xIC / MS_pooled), where MS_pooled = (SS_T(P(Aft)) xCC + SS_Res) / (df_T(P(Aft)) xCC + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2031"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2031"]<-"Can be pooled"

        # 	F = (MS_P(Aft) xS(CC) / MS_pooled), where MS_pooled = (SS_T(P(Aft)) xS(CC) + SS_Res) / (df_T(P(Aft)) xS(CC) + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2040"]<-NA
        baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2041"]<-"Can be pooled"

        # 	If varP(Aft) xS(CC) <= 0, then F = (MS_P(Aft) xS(IC) / MSpooled), where MSpooled = (SS_P(Aft) xS(CC) + SS_Res) / (df_P(Aft) xS(CC) + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2031"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2031"]<-"Can be pooled"

        baci$Post.Hoc.Pooling[baci$ID == AEB08[AEB08n]]<-NA
        # If both terms cannot be eliminated (* or ns)
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2041"] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == AEB08[AEB08n]]<-"Cannot be pooled"
        # 	If varT(P(Aft)) xCC <= and varP(Aft) xS(CC) > 0, then F = (MS_P(Aft) xCC / MSpooled), where MSpooled = (SS_P(Aft) xS(CC) + SS_Res) / (df_P(Aft) xS(CC) + df_Res)
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2041"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2043"]<-"Can be pooled"
        # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xCC > 0, then F = (MS_P(Aft) xCC / MS_pooled),
        # where MSpooled = (SS_T(P(Aft)) xCC + SS_Res) / (df_T(P(Aft)) xCC + df_Res)
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2041"] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2044"]<-"Can be pooled"
        # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xCC <= 0, then F = (MSP(Aft) xCC / MSpooled),
        # where MSpooled = (SS_P(Aft) xS(CC) + SS_T(P(Aft)) xCC + SS_T(P(Aft)) xS(CC) + SSRes) / (dfP(Aft) xS(CC) + dfT(P(Aft)) xCC + dfT(P(Aft)) xS(CC) + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2041"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2045"]<-"Can be pooled"
        AEB08e<- c("08AEB.2040","08AEB.2041","08AEB.2042","08AEB.2043","08AEB.2044","08AEB.2045")
        AEB08en<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2040"]),
                         is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2041"]),
                         is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2042"]),
                         is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2043"]),
                         is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2044"]),
                         is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2045"])) == F)

        # 	If varP(Aft) xCC <= 0, then F = (MS_P(Aft) xIC / MS_pooled), where MS_pooled = (SS_P(Aft) xCC + SS_Res) / (df_P(Aft) xCC + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "08AEB.2030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08e[AEB08en]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2031"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08e[AEB08en]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2031"]<-"Can be pooled"

        # 	F = (MS_B xS(CC) / MS_pooled), where MS_pooled = (SS_P(Aft) xS(CC) + SS_Res) / (df_P(Aft) xS(CC) + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "07AEG.0040"]<-NA
        baci$Post.Hoc.Pooling[baci$ID == "07AEG.0041"]<-"Can be pooled"

        #If varB xS(CC) <= 0, then F = (MS_B xS(IC) / MS_pooled), where MSpooled = (SS_B xS(CC) + SS_Res) / (df_B xS(CC) + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "07AEG.0030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0031"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0031"]<-"Can be pooled"

        baci$Post.Hoc.Pooling[baci$ID == AE03[AE03n]]<-NA
        # If both terms cannot be eliminated (* or ns)
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08e[AEB08en]] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == AE03[AE03n]]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08e[AEB08en]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2041"] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == AE03[AE03n]]<-"Cannot be pooled"
        # 	If varP(Aft) xCC <= 0 and varB xS(CC) > 0 then F = (MS_B xCC / MSpooled), where MSpooled = (SS_B xS(CC) + SS_Res) / (df_B xS(CC) + df_Res)
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08e[AEB08en]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0045"]<-"Can be pooled"
        # 	If varB xS(CC) <= 0 and varP(Aft) xCC > 0, then F = (MS_B xCC / MSpooled), where MSpooled = (SS_P(Aft) xCC + SS_Res) / (df_P(Aft) xCC + df_Res)
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08e[AEB08en]] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0046"]<-"Can be pooled"
        # 	varB xS(CC) <= 0, varP(Aft) xCC <= 0, varT(P(Aft)) xCC <= 0 and varP(Aft) xS(CC) > 0, then F = (MS_B xCC / MSpooled),
        # where MSpooled = (SS_B xS(CC) + SS_P(Aft) xCC + SS_P(Aft) xS(CC) + SS_Res) / (df_B xS(CC) + df_P(Aft) xCC + df_P(Aft) xS(CC) + df_Res)
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08e[AEB08en]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2041"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2041"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0047"]<-"Can be pooled"
        # 	If varB xS(CC) <= 0, varP(Aft) xCC <= 0, varP(Aft) xS(CC) <= 0, varT(P(Aft)) xCC > 0, then F = (MSB xCC / MSpooled),
        # where MSpooled = (SS_B xS(CC) + SS_P(Aft) xCC + SS_T(P(Aft)) xCC + SS_Res) / (df_B xS(CC) + df_P(Aft) xCC + df_T(P(Aft)) xCC + df_Res)
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08e[AEB08en]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2041"] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0048"]<-"Can be pooled"
        # 	If varB xS(CC) <= 0, varP(Aft) xCC <= 0, varP(Aft) xS(CC) <= 0, varT(P(Aft)) xCC <= 0, then F = (MSB xCC / MSpooled),
        #where MSpooled = (SS_B xS(CC) + SS_P(Aft) xCC + SS_P(Aft) xS(CC) + SS_T(P(Aft)) xCC + SS_T(P(Aft)) xS(CC) + SS_Res) / (df_B xS(CC) + df_P(Aft) xCC + df_P(Aft) xS(CC) + df_T(P(Aft)) xCC+ df_T(P(Aft)) xS(CC) + df_Res)
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08e[AEB08en]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "07AEG.0041"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2041"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2041"] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0049"]<-"Can be pooled"

        # 	If varB xCC <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SS_B xCC + SS_Res) / (df_B xCC + df_Res)
        baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-NA
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Can be pooled"
      }
      if(method == "pool") {
        # 	F = (MS_T(P(Aft)) xIC / MS_T(P(Aft)) xCC), or
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2030"]<-"Cannot be pooled"
        # 	If varT(P(Aft)) xCC <= 0 and varT(P(Aft)) xS(IC) > 0, then F = (MS_T(P(Aft)) xIC / MS_pooled),
        # where MS_pooled = (SS_T(P(Aft)) xCC + SS_T(P(Aft)) xS(IC)) / (df_T(P(Aft)) xCC + df_T(P(Aft)) xS(IC)) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBGC11[AEBGC11n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2031"]<-"Can be pooled"
        # 	If varT(P(Aft)) xCC <= 0 and varT(P(Aft)) xS(IC) <= 0, then F = (MST(P(Aft)) xIC / MSpooled), where MSpooled = SS1 + SSRes / df1 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "10AEBC.2040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBGC11[AEBGC11n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2032"]<-"Can be pooled"
        AEBC10<- c("10AEBC.2030","10AEBC.2031","10AEBC.2032")
        AEBC10n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2030"]),
                         is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2031"]),
                         is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2032"])) == F)

        # 	F = (MS_P(Aft) xS(IC) / MS_P(Aft) xS(CC)), or
        if(baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2030"]<-"Cannot be pooled"
        # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xS(IC) > 0, then F = (MS_P(Aft) xS(IC) / MS_pooled),
        #where MSpooled = (SS_P(Aft) xS(CC) + SS_T(P(Aft)) xS(IC)) / (df_P(Aft) xS(CC) + df_T(P(Aft)) xS(IC)) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBGC11[AEBGC11n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2031"]<-"Can be pooled"
        # 	If varP(Aft) xS(CC) <= 0 and varT(P(Aft)) xS(IC) <= 0, then F = (MS_P(Aft) xS(IC) / MS_pooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "09AEBG.2040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBGC11[AEBGC11n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2032"]<-"Can be pooled"
        AEBG09<- c("09AEBG.2030","09AEBG.2031","09AEBG.2032")
        AEBG09n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2030"]),
                         is.na(baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2031"]),
                         is.na(baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2032"])) == F)

        # If the two terms cannot be eliminated (* or ns)
        # 	F = (MS_P(Aft) xIC / MS_P(Aft) xCC)
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] < p.pooling |
           (baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
            baci$P.Value.Upper.Tail[baci$ID == AEBC10[AEBC10n]] < p.pooling &
            baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] < p.pooling)) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2030"]<-"Cannot be pooled"
        # 	If varP(Aft) xCC <= 0, varT(P(Aft)) xIC <= and varP(Aft) xS(IC) > 0, then F = (MS_P(Aft) xIC / MSpooled), where MSpooled = (SS_P(Aft) xCC + SS_P(Aft) xS(IC)) / (df_P(Aft) xCC + df_P(Aft) xS(IC))
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBC10[AEBC10n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2031"]<-"Can be pooled"
        # 	If varP(Aft) xCC <= 0, varP(Aft) xS(IC) <= 0 and varT(P(Aft)) xIC > 0, then F = (MSP(Aft) xIC / MSpooled), where MSpooled = (SS_P(Aft) xCC + SS_T(P(Aft)) xIC) / (df_P(Aft) xCC + df_T(P(Aft)) xIC)
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBC10[AEBC10n]] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2032"]<-"Can be pooled"
        # 	If varP(Aft) xCC <= 0, varP(Aft) xS(IC) <= 0, varT(P(Aft)) xIC <= 0, and varT(P(Aft)) xS(IC) > 0, then F = (MSP(Aft) xIC / MSpooled), where MSpooled = (SS_P(Aft) xCC + SS_P(Aft) xS(IC) + SS_T(P(Aft)) xIC + SS_T(P(Aft)) xS(IC)) / (dfP(Aft) xCC + dfP(Aft) xS(IC) + dfT(P(Aft)) xIC + dfT(P(Aft)) xS(IC)) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBC10[AEBC10n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBGC11[AEBGC11n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2033"]<-"Can be pooled"
        # 	If varP(Aft) xCC <= 0, varT(P(Aft)) xIC <= 0, varP(Aft) xS(IC) <= 0) and varT(P(Aft)) xS(IC) <= 0, then F = (MSP(Aft) xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == AEB08[AEB08n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBC10[AEBC10n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBGC11[AEBGC11n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "08AEB.2034"]<-"Can be pooled"
        AEB08.30<- c("08AEB.2030","08AEB.2031","08AEB.2032","08AEB.2033","08AEB.2034")
        AEB08.30n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2030"]),
                           is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2031"]),
                           is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2032"]),
                           is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2033"]),
                           is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2034"])) == F)

        # 	F = (MS_B xS(IC) / MS_B xS(CC)), or
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0030"]<-"Cannot be pooled"
        # 	If varB xS(CC) <= 0 and varP(Aft) xS(IC) > 0, then F = (MSB xS(IC) / MSpooled),
        # where MSpooled = (SSB xS(CC) + SS_P(Aft) xS(IC)) / (dfB xS(CC) + dfP(Aft) xS(IC)) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0031"]<-"Can be pooled"
        # 	If varB xS(CC) <= 0, varP(Aft) xS(IC) <= 0 and varT(P(Aft)) xS(IC) > 0, then F = (MSB xS(IC) / MSpooled),
        # where MSpooled = (SS1 + SS_T(P(Aft)) xS(IC)) / (df_1 + dfT(P(Aft)) xS(IC)) = SS2/df2
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBGC11[AEBGC11n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0032"]<-"Can be pooled"
        # 	If varB xS(CC) <= 0, varP(Aft) xS(IC) <= 0 and varT(P(Aft)) xS(IC) <= 0, then F = (MSB xS(IC) / MSpooled),
        # where MSpooled = (SS2 + SSRes) / (df2 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == "07AEG.0040"] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBGC11[AEBGC11n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "07AEG.0033"]<-"Can be pooled"
        AEG07<- c("07AEG.0030","07AEG.0031","07AEG.0032","07AEG.0033")
        AEG07n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0030"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0031"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0032"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0033"])) == F)

        # If the two terms cannot be eliminated (* or ns)
        # 	F = (MS_B xIC / MS_B xCC), or
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08.30[AEB08.30n]] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08.30[AEB08.30n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBC10[AEBC10n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]<-"Cannot be pooled"
        # 	If varB xCC <= 0, varP(Aft) xIC <= 0 and varB xS(IC) > 0 then F = (MSB xIC / MSpooled),
        # where MSpooled = (SS_B xCC + SS_B xS(IC)) / (df_B xCC + df_B xS(IC))
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08.30[AEB08.30n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]<-"Can be pooled"
        # 	If varB xCC <= 0, varB xS(IC) <= 0 and varP(Aft) xIC > 0, then F = (MSB xIC / MSpooled),
        # where MSpooled = (SS_B xCC + SS_P(Aft) xIC) / (dfB xCC + dfP(Aft) xIC)
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08.30[AEB08.30n]] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0032"]<-"Can be pooled"
        # 	If varB xCC <= 0, varB xS(IC) <= 0, varP(Aft) xIC <= 0, varP(Aft) xS(IC) <= 0 and varT(P(Aft)) xIC > 0, then F = (MSB xIC / MSpooled),
        # where MSpooled = (SS_B xCC + SS_B xS(IC) + SS_P(Aft) xIC + SS_P(Aft) xS(IC)) / (df_B xCC + df_B xS(IC) + df_P(Aft) xIC + df_P(Aft) xS(IC))
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08.30[AEB08.30n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] >= p.pooling  &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBC10[AEBC10n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0033"]<-"Can be pooled"
        # 	If varB xCC <= 0, varB xS(IC) <= 0, varP(Aft) xIC <= 0, varT(P(Aft)) xIC <= 0 and varP(Aft) xS(IC) > 0, then F = (MSB xIC / MSpooled),
        # where MSpooled = (SS_B xCC + SS_B xS(IC) + SS_P(Aft) xIC + SS_T(P(Aft)) xIC) / (df_B xCC + df_B xS(IC) + df_P(Aft) xIC + df_T(P(Aft)) xIC)
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08.30[AEB08.30n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] >= p.pooling  &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] < p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBC10[AEBC10n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0034"]<-"Can be pooled"
        # 	If varB xCC <= 0, varB xS(IC) <= 0, varP(Aft) xIC <= 0, varP(Aft) xS(IC) <= 0, varT(P(Aft)) xIC <= 0 and varT(P(Aft)) xS(IC) > 0, then F = (MSB xIC / MSpooled),
        # where MSpooled = (SS_B xCC + SS_B xS(IC) + SS_P(Aft) xIC + SS_P(Aft) xS(IC) + SS_T(P(Aft)) xIC + SS_T(P(Aft)) xS(IC)) / (dfB xCC + dfB xS(IC) + dfP(Aft) xIC + dfP(Aft) xS(IC) + dfT(P(Aft)) xIC + dfT(P(Aft)) xS(IC)) = SS1/df1
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08.30[AEB08.30n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] >= p.pooling  &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBC10[AEBC10n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBGC11[AEBGC11n]] < p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0035"]<-"Can be pooled"
        # 	If varB xCC <= 0, varB xS(IC) <= 0, varP(Aft) xIC <= 0, varP(Aft) xS(IC) <= 0, varT(P(Aft)) xIC <= 0 and varT(P(Aft)) xS(IC) <= 0, then F = (MSB xIC / MSpooled), where MSpooled = (SS1 + SSRes) / (df1 + dfRes)
        if(baci$P.Value.Upper.Tail[baci$ID == AE03[AE03n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEB08.30[AEB08.30n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEG07[AEG07n]] >= p.pooling  &
           baci$P.Value.Upper.Tail[baci$ID == AEBG09[AEBG09n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBC10[AEBC10n]] >= p.pooling &
           baci$P.Value.Upper.Tail[baci$ID == AEBGC11[AEBGC11n]] >= p.pooling) baci$Post.Hoc.Pooling[baci$ID == "03AE.0036"]<-"Can be pooled"
      }
    }
  }

  ###------------------------------------------ Next.Step &  Interpretation

  if(aav.model == "AEBGC") {
    AEBGC11<- c("11AEBGC.2030","11AEBGC.2031")
    AEBGC11n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "11AEBGC.2030"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "11AEBGC.2031"])) == F)

    # Sites vs Times
    # T(P(B))xS(I) = ns
    if(baci$Significant[baci$ID == AEBGC11[AEBGC11n]] != "*") {
      # T(P(B))xS(C) = ns; T(P(B))xS(I) = ns
      if(baci$Significant[baci$ID == "11AEBGC.2040"] != "*") {
        baci$Next.Step[baci$ID == AEBGC11[AEBGC11n]]<-baci$Next.Step[baci$ID == "11AEBGC.2040"]<- "Go to P(Aft)xS(C) and T(P(Aft))xC"
        baci$Interpretation[baci$ID == AEBGC11[AEBGC11n]]<-baci$Interpretation[baci$ID == "11AEBGC.2040"]<- "No short-term and small-scale interactions"
      }
      # T(P(B))xS(C) = *; T(P(B))xS(I) = ns
      if(baci$Significant[baci$ID == "11AEBGC.2040"] == "*") {
        baci$Next.Step[baci$ID == AEBGC11[AEBGC11n]]<-baci$Next.Step[baci$ID == "11AEBGC.2040"]<- "END"
        baci$Interpretation[baci$ID == AEBGC11[AEBGC11n]]<-baci$Interpretation[baci$ID == "11AEBGC.2040"]<- "Short-term and small-scale interactions: temporal trends are not associated with the disturbed sites"
      }
    }

    # T(P(B))xS(I) = *
    # T(P(B))xS(C) = ns; T(P(B))xS(I) = *   &       # T(P(B))xS(C) = *; T(P(B))xS(I) = ns
    if(baci$Significant[baci$ID == AEBGC11[AEBGC11n]] == "*") {
        baci$Next.Step[baci$ID == AEBGC11[AEBGC11n]]<-baci$Next.Step[baci$ID == "11AEBGC.2040"]<- "Go to T(P(Bef))xS(I) (2-tail test)"
        baci$Interpretation[baci$ID == AEBGC11[AEBGC11n]]<-baci$Interpretation[baci$ID == "11AEBGC.2040"]<- "Short-term and small-scale interactions"
        # 2-tail test
        baci$Next.Step[baci$ID == "11AEBGC.1030"]<-baci$Next.Step[baci$ID == "11AEBGC.1040"]<- "END"
        ifelse(baci$Significant[baci$ID == "11AEBGC.1030"] == "*" & baci$Significant[baci$ID =="11AEBGC.1040"] != "*",
               baci$Interpretation[baci$ID == "11AEBGC.1030"]<-baci$Interpretation[baci$ID == "11AEBGC.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and small-scale: impact sites differs from before and after the disturbance occurs",
               baci$Interpretation[baci$ID == "11AEBGC.1030"]<-baci$Interpretation[baci$ID == "11AEBGC.1040"]<- "No short-term and small-scale impact detected: the changes are not associated with the disturbed sites")
        #"Short-term IMPACT DETECTED: impact sites differs from before and after the disturbance occurs"
        #"Short-term IMPACT DETECTED: control sites do not differs from before and after the disturbance"
        #### Therefore effect is specific to the impacted sites and coincident with the start of the disturbance. IMPACT DETECTED. ###
      }

#    # Sites vs Times
#    # T(P(B))xS(C) = ns; T(P(B))xS(I) = ns
#    if(baci$Significant[baci$ID == "11AEBGC.2040"] != "*") {
#      if(baci$Significant[baci$ID == AEBGC11[AEBGC11n]] != "*") {
#        baci$Next.Step[baci$ID == AEBGC11[AEBGC11n]]<-baci$Next.Step[baci$ID == "11AEBGC.2040"]<- "Go to P(Aft)xS(C) and T(P(Aft))xC"
#        baci$Interpretation[baci$ID == AEBGC11[AEBGC11n]]<-baci$Interpretation[baci$ID == "11AEBGC.2040"]<- "No short-term and small-scale interactions"
#      }
#      # T(P(B))xS(C) = ns; T(P(B))xS(I) = *
#      if(baci$Significant[baci$ID == AEBGC11[AEBGC11n]] == "*") {
#        baci$Next.Step[baci$ID == AEBGC11[AEBGC11n]]<-baci$Next.Step[baci$ID == "11AEBGC.2040"]<- "Go to T(P(Bef))xS(I) (2-tail test)"
#        baci$Interpretation[baci$ID == AEBGC11[AEBGC11n]]<-baci$Interpretation[baci$ID == "11AEBGC.2040"]<- "Short-term and small-scale interactions"
#        # 2-tail test
#        baci$Next.Step[baci$ID == "11AEBGC.1030"]<-baci$Next.Step[baci$ID == "11AEBGC.1040"]<- "END"
#        ifelse(baci$Significant[baci$ID == "11AEBGC.1030"] == "*" & baci$Significant[baci$ID =="11AEBGC.1040"] != "*",
#               baci$Interpretation[baci$ID == "11AEBGC.1030"]<-baci$Interpretation[baci$ID == "11AEBGC.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and small-scale: impact sites differs from before and after the disturbance occurs",
#               baci$Interpretation[baci$ID == "11AEBGC.1030"]<-baci$Interpretation[baci$ID == "11AEBGC.1040"]<- "No short-term and small-scale impact detected: the changes are not associated with the disturbed sites")
#        #"Short-term IMPACT DETECTED: impact sites differs from before and after the disturbance occurs"
#        #"Short-term IMPACT DETECTED: control sites do not differs from before and after the disturbance"
#        #### Therefore effect is specific to the impacted sites and coincident with the start of the disturbance. IMPACT DETECTED. ###
#      }
#    }
#    # T(P(B))xS(C) = *; T(P(B))xS(I) = ns
#    if(baci$Significant[baci$ID == "11AEBGC.2040"] == "*") {
#      if(baci$Significant[baci$ID == AEBGC11[AEBGC11n]] != "*") {
#        baci$Next.Step[baci$ID == AEBGC11[AEBGC11n]]<-baci$Next.Step[baci$ID == "11AEBGC.2040"]<- "Go to P(Aft)xS(C) and T(P(Aft))xC"
#        baci$Interpretation[baci$ID == AEBGC11[AEBGC11n]]<-baci$Interpretation[baci$ID == "11AEBGC.2040"]<- "Short-term and small-scale interactions: temporal trends are not associated with the disturbed sites"
#      }
#      # T(P(B))xS(C) = *; T(P(B))xS(I) = *
#      if(baci$Significant[baci$ID == AEBGC11[AEBGC11n]] == "*") {
#        baci$Next.Step[baci$ID == AEBGC11[AEBGC11n]]<-baci$Next.Step[baci$ID == "11AEBGC.2040"]<- "Go to P(Aft)xS(C) and T(P(Aft))xC"
#        baci$Interpretation[baci$ID == AEBGC11[AEBGC11n]]<-baci$Interpretation[baci$ID == "11AEBGC.2040"]<- "Short-term and small-scale interactions: widespread temporal trends are associated with both, control and disturbed sites"
#      }
#    }


    # Sites vs Periods P(B) x S(L)
    AEBG09.40<- c("09AEBG.2040","09AEBG.2041")
    AEBG09.40n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2040"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2041"])) == F)
    AEBG09.30<- c("09AEBG.2030","09AEBG.2031","09AEBG.2032")
    AEBG09.30n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2030"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2031"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "09AEBG.2032"])) == F)

    # P(B)xS(I) = ns
    # P(B)xS(C) = ns; P(B)xS(I) = ns
    if(baci$Significant[baci$ID == AEBG09.30[AEBG09.30n]] != "*") {
        if(baci$Significant[baci$ID == AEBG09.40[AEBG09.40n]] != "*") {
          baci$Next.Step[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Next.Step[baci$ID == AEBG09.40[AEBG09.40n]]<- "Go to BxS(C) and P(Aft)xC"
        baci$Interpretation[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Interpretation[baci$ID == AEBG09.40[AEBG09.40n]]<- "No medium-term and small-scale interactions"
        }
      # P(B)xS(C) = *; P(B)xS(I) = ns
        if(baci$Significant[baci$ID == AEBG09.40[AEBG09.40n]] == "*") {
          baci$Next.Step[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Next.Step[baci$ID == AEBG09.40[AEBG09.40n]]<- "END"
          baci$Interpretation[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Interpretation[baci$ID == AEBG09.40[AEBG09.40n]]<- "Medium-term and small-scale interactions: temporal trends are not associated with the disturbed sites"
        }
      }
      # P(B)xS(I) = *
      # P(B)xS(C) = ns; P(B)xS(I) = *  &   P(B)xS(C) = *; P(B)xS(I) = *
      if(baci$Significant[baci$ID == AEBG09.30[AEBG09.30n]] == "*") {
        baci$Next.Step[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Next.Step[baci$ID == AEBG09.40[AEBG09.40n]]<- "Go to P(Bef)xS(I) (2-tail test)"
        baci$Interpretation[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Interpretation[baci$ID == AEBG09.40[AEBG09.40n]]<- "Medium-term and small-scale interactions"
        # 2-tail test
        baci$Next.Step[baci$ID == "09AEBG.1030"]<-baci$Next.Step[baci$ID == "09AEBG.1040"]<- "END"
        ifelse(baci$Significant[baci$ID == "09AEBG.1030"] == "*" & baci$Significant[baci$ID =="09AEBG.1040"] != "*",
               baci$Interpretation[baci$ID == "09AEBG.1030"]<-baci$Interpretation[baci$ID == "09AEBG.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from medium-term and small-scale: impact sites differs from before and after the disturbance occurs",
               baci$Interpretation[baci$ID == "09AEBG.1030"]<-baci$Interpretation[baci$ID == "09AEBG.1040"]<- "No medium-term and small-scale impact detected: the changes are not associated with the disturbed location")
      }

#    # P(B)xS(C) = ns; P(B)xS(I) = ns
#    if(baci$Significant[baci$ID == AEBG09.40[AEBG09.40n]] != "*") {
#      if(baci$Significant[baci$ID == AEBG09.30[AEBG09.30n]] != "*") {
#        baci$Next.Step[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Next.Step[baci$ID == AEBG09.40[AEBG09.40n]]<- "Go to BxS(C) and P(Aft)xC"
#        baci$Interpretation[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Interpretation[baci$ID == AEBG09.40[AEBG09.40n]]<- "No medium-term and small-scale interactions"
#      }
#      # P(B)xS(C) = ns; P(B)xS(I) = *
#      if(baci$Significant[baci$ID == AEBG09.30[AEBG09.30n]] == "*") {
#        baci$Next.Step[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Next.Step[baci$ID == AEBG09.40[AEBG09.40n]]<- "Go to P(B)xS(I) (2-tail test), BxS(C) and P(Aft)xC"
#        baci$Interpretation[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Interpretation[baci$ID == AEBG09.40[AEBG09.40n]]<- "Medium-term and small-scale interactions"
#        # 2-tail test
#        baci$Next.Step[baci$ID == "09AEBG.1030"]<-baci$Next.Step[baci$ID == "09AEBG.1040"]<- "END"
#        ifelse(baci$Significant[baci$ID == "09AEBG.1030"] == "*" & baci$Significant[baci$ID =="09AEBG.1040"] != "*",
#               baci$Interpretation[baci$ID == "09AEBG.1030"]<-baci$Interpretation[baci$ID == "09AEBG.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from medium-term and small-scale: impact sites differs from before and after the disturbance occurs",
#               baci$Interpretation[baci$ID == "09AEBG.1030"]<-baci$Interpretation[baci$ID == "09AEBG.1040"]<- "No medium-term and small-scale impact detected: the changes are not associated with the disturbed location")
#      }
#    }
#    # P(B)xS(C) = *; P(B)xS(I) = ns
#    if(baci$Significant[baci$ID == AEBG09.40[AEBG09.40n]] == "*") {
#      if(baci$Significant[baci$ID == AEBG09.30[AEBG09.30n]] != "*") {
#        baci$Next.Step[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Next.Step[baci$ID == AEBG09.40[AEBG09.40n]]<- "Go to BxS(C) and P(Aft)xC"
#        baci$Interpretation[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Interpretation[baci$ID == AEBG09.40[AEBG09.40n]]<- "Medium-term and small-scale interactions: temporal trends are not associated with the disturbed sites"
#      }
#      # P(B)xS(C) = *; P(B)xS(I) = *
#      if(baci$Significant[baci$ID == AEBG09.30[AEBG09.30n]] == "*") {
#        baci$Next.Step[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Next.Step[baci$ID == AEBG09.40[AEBG09.40n]]<- "Go to BxS(C) and P(Aft)xC"
#        baci$Interpretation[baci$ID == AEBG09.30[AEBG09.30n]]<-baci$Interpretation[baci$ID == AEBG09.40[AEBG09.40n]]<- "Medium-term and small-scale interactions: widespread temporal trends are associated with both, control and disturbed sites"
#      }
#    }
#  }
  }


  ### Sites vs Before and after - B x S(L)
  if(aav.model == "AEBGC" | aav.model == "AEGC") {
    AEG07.40<- c("07AEG.0040","07AEG.0041")
    AEG07.40n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0040"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0041"])) == F)
    AEG07.30<- c("07AEG.0030","07AEG.0031","07AEG.0032","07AEG.0033")
    AEG07.30n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0030"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0031"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0032"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "07AEG.0033"])) == F)

    # BxS(I) = ns
    # BxS(C) = ns; BxS(I) = ns
    if(baci$Significant[baci$ID == AEG07.30[AEG07.30n]] != "*") {
      if(baci$Significant[baci$ID == AEG07.40[AEG07.40n]] != "*") {
        baci$Next.Step[baci$ID == AEG07.30[AEG07.30n]]<-baci$Next.Step[baci$ID == AEG07.40[AEG07.40n]]<- "Go to BxC"
        baci$Interpretation[baci$ID == AEG07.30[AEG07.30n]]<-baci$Interpretation[baci$ID == AEG07.40[AEG07.40n]]<- "No long-term and small-scale interactions"
      }
      # BxS(C) = *; BxS(I) = ns
      if(baci$Significant[baci$ID == AEG07.40[AEG07.40n]] == "*") {
        baci$Next.Step[baci$ID == AEG07.30[AEG07.30n]]<-baci$Next.Step[baci$ID == AEG07.40[AEG07.40n]]<- "END"
        baci$Interpretation[baci$ID == AEG07.30[AEG07.30n]]<-baci$Interpretation[baci$ID == AEG07.40[AEG07.40n]]<- "Long-term and small-scale interactions: temporal trends are not associated with the disturbed sites"
        }
      }
    # BxS(I) = *
    # BxS(C) = ns; BxS(I) = *
      if(baci$Significant[baci$ID == AEG07.30[AEG07.30n]] == "*") {
        baci$Next.Step[baci$ID == AEG07.30[AEG07.30n]]<-baci$Next.Step[baci$ID == AEG07.40[AEG07.40n]]<- "END"
        baci$Interpretation[baci$ID == AEG07.30[AEG07.30n]]<-baci$Interpretation[baci$ID == AEG07.40[AEG07.40n]]<- "PRESS DISTURBANCE - IMPACT DETECTED from long-term and small-scale: Disturbed sites differs from controls sites"
      }

#    # BxS(C) = ns; BxS(I) = ns
#    if(baci$Significant[baci$ID == AEG07.40[AEG07.40n]] != "*") {
#      if(baci$Significant[baci$ID == AEG07.30[AEG07.30n]] != "*") {
#        baci$Next.Step[baci$ID == AEG07.30[AEG07.30n]]<-baci$Next.Step[baci$ID == AEG07.40[AEG07.40n]]<- "Go to BxC"
#        baci$Interpretation[baci$ID == AEG07.30[AEG07.30n]]<-baci$Interpretation[baci$ID == AEG07.40[AEG07.40n]]<- "No long-term and small-scale interactions"
#      }
#      # BxS(C) = ns; BxS(I) = *
#      if(baci$Significant[baci$ID == AEG07.30[AEG07.30n]] == "*") {
#        baci$Next.Step[baci$ID == AEG07.30[AEG07.30n]]<-baci$Next.Step[baci$ID == AEG07.40[AEG07.40n]]<- "Go to BxC"
#        baci$Interpretation[baci$ID == AEG07.30[AEG07.30n]]<-baci$Interpretation[baci$ID == AEG07.40[AEG07.40n]]<- "PRESS DISTURBANCE - IMPACT DETECTED from long-term and small-scale: Disturbed sites differs from controls sites"
#      }
#    }
#    # BxS(C) = *; BxS(I) = ns
#    if(baci$Significant[baci$ID == AEG07.40[AEG07.40n]] == "*") {
#      if(baci$Significant[baci$ID == AEG07.30[AEG07.30n]] != "*") {
#        baci$Next.Step[baci$ID == AEG07.30[AEG07.30n]]<-baci$Next.Step[baci$ID == AEG07.40[AEG07.40n]]<- "Go to BxC"
#        baci$Interpretation[baci$ID == AEG07.30[AEG07.30n]]<-baci$Interpretation[baci$ID == AEG07.40[AEG07.40n]]<- "Long-term and small-scale interactions: temporal trends are not associated with the disturbed sites"
#      }
#      # BxS(C) = *; BxS(I) = *
#      if(baci$Significant[baci$ID == AEG07.30[AEG07.30n]] == "*") {
#        baci$Next.Step[baci$ID == AEG07.30[AEG07.30n]]<-baci$Next.Step[baci$ID == AEG07.40[AEG07.40n]]<- "Go to BxC"
#        baci$Interpretation[baci$ID == AEG07.30[AEG07.30n]]<-baci$Interpretation[baci$ID == AEG07.40[AEG07.40n]]<- "Long-term and small-scale interactions: widespread temporal trends are associated with both, control and disturbed sites"
#      }
#    }
#  }
  }


  if(aav.model == "AEBGC" | aav.model == "AEBC") {

    ### Location vs Period in Time - T(P(B)) x L
    AEBC10.40<- c("10AEBC.2040","10AEBC.2041")
    AEBC10.40n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2040"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2041"])) == F)
    AEBC10.30<- c("10AEBC.2030","10AEBC.2031","10AEBC.2032")
    AEBC10.30n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2030"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2031"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEBC.2032"])) == F)

    # T(P(B))xI = ns
    # T(P(B))xC = ns; T(P(B))xI = ns
      if(baci$Significant[baci$ID == AEBC10.30[AEBC10.30n]] != "*") {
        baci$Next.Step[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Next.Step[baci$ID == AEBC10.40[AEBC10.40n]]<- "Go to P(Aft)xC"
        baci$Interpretation[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Interpretation[baci$ID == AEBC10.40[AEBC10.40n]]<- "No short-term and local-scale interactions"
      }
    # T(P(B))xC = *; T(P(B))xI = ns
    if(baci$Significant[baci$ID == AEBC10.40[AEBC10.40n]] == "*") {
        baci$Next.Step[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Next.Step[baci$ID == AEBC10.40[AEBC10.40n]]<- "END"
        baci$Interpretation[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Interpretation[baci$ID == AEBC10.40[AEBC10.40n]]<- "Short-term and local-scale interactions: temporal trends are not associated with the disturbed location"
      }
    # T(P(B))xI = *
    # T(P(B))xC = ns; T(P(B))xI = *
      if(baci$Significant[baci$ID == AEBC10.30[AEBC10.30n]] == "*") {
        baci$Next.Step[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Next.Step[baci$ID == AEBC10.40[AEBC10.40n]]<- "Go to T(P(Bef))xI (2-tail test)"
        baci$Interpretation[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Interpretation[baci$ID == AEBC10.40[AEBC10.40n]]<- "Short-term and local-scale interactions"
        # 2-tail test
        baci$Next.Step[baci$ID == "10AEBC.1030"]<-baci$Next.Step[baci$ID == "10AEBC.1040"]<- "END"
        ifelse(baci$Significant[baci$ID == "10AEBC.1030"] == "*" & baci$Significant[baci$ID =="10AEBC.1040"] != "*",
               baci$Interpretation[baci$ID == "10AEBC.1030"]<-baci$Interpretation[baci$ID == "10AEBC.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and local-scale: impact location differs from before and after the disturbance occurs",
               baci$Interpretation[baci$ID == "10AEBC.1030"]<-baci$Interpretation[baci$ID == "10AEBC.1040"]<- "No short-term and local-scale impact detected: the changes are not associated with the disturbed location")
      }

#      # T(P(B))xC = ns; T(P(B))xI = ns
#      if(baci$Significant[baci$ID == AEBC10.40[AEBC10.40n]] != "*") {
#        if(baci$Significant[baci$ID == AEBC10.30[AEBC10.30n]] != "*") {
#          baci$Next.Step[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Next.Step[baci$ID == AEBC10.40[AEBC10.40n]]<- "Go to P(Aft)xC"
#          baci$Interpretation[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Interpretation[baci$ID == AEBC10.40[AEBC10.40n]]<- "No short-term and local-scale interactions"
#        }
#        # T(P(B))xC = ns; T(P(B))xI = *
#        if(baci$Significant[baci$ID == AEBC10.30[AEBC10.30n]] == "*") {
#          baci$Next.Step[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Next.Step[baci$ID == AEBC10.40[AEBC10.40n]]<- "Go to T(P(Bef))xI (2-tail test) and P(Aft)xC"
#          baci$Interpretation[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Interpretation[baci$ID == AEBC10.40[AEBC10.40n]]<- "Short-term and local-scale interactions"
#          # 2-tail test
#          baci$Next.Step[baci$ID == "10AEBC.1030"]<-baci$Next.Step[baci$ID == "10AEBC.1040"]<- "END"
#          ifelse(baci$Significant[baci$ID == "10AEBC.1030"] == "*" & baci$Significant[baci$ID =="10AEBC.1040"] != "*",
#                 baci$Interpretation[baci$ID == "10AEBC.1030"]<-baci$Interpretation[baci$ID == "10AEBC.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and local-scale: impact location differs from before and after the disturbance occurs",
#                 baci$Interpretation[baci$ID == "10AEBC.1030"]<-baci$Interpretation[baci$ID == "10AEBC.1040"]<- "No short-term and local-scale impact detected: the changes are not associated with the disturbed location")
#        }
#      }
#      # T(P(B))xC = *; T(P(B))xI = ns
#      if(baci$Significant[baci$ID == AEBC10.40[AEBC10.40n]] == "*") {
#        if(baci$Significant[baci$ID == AEBC10.30[AEBC10.30n]] != "*") {
#          baci$Next.Step[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Next.Step[baci$ID == AEBC10.40[AEBC10.40n]]<- "Go to P(Aft)xC"
#          baci$Interpretation[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Interpretation[baci$ID == AEBC10.40[AEBC10.40n]]<- "Short-term and local-scale interactions: temporal trends are not associated with the disturbed location"
#        }
#        # T(P(B))xC = *; T(P(B))xI = *
#        if(baci$Significant[baci$ID == AEBC10.30[AEBC10.30n]] == "*") {
#          baci$Next.Step[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Next.Step[baci$ID == AEBC10.40[AEBC10.40n]]<- "Go to P(Aft)xC"
#          baci$Interpretation[baci$ID == AEBC10.30[AEBC10.30n]]<-baci$Interpretation[baci$ID == AEBC10.40[AEBC10.40n]]<- "Short-term and local-scale interactions: widespread temporal trends are associated with both, control and disturbed locations"
#        }



    ### Location vs Period - P(B) x L
    AEB08.40<- c("08AEB.2040","08AEB.2041","08AEB.2042","08AEB.2043","08AEB.2044","08AEB.2045")
    AEB08.40n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2040"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2041"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2042"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2043"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2044"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2045"])) == F)
    AEB08.30<- c("08AEB.2030","08AEB.2031","08AEB.2032","08AEB.2033","08AEB.2034")
    AEB08.30n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2030"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2031"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2032"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2033"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "08AEB.2034"])) == F)

    # P(B)xI = ns
    # P(B)xC = ns; P(B)xI = ns
    if(baci$Significant[baci$ID == AEB08.30[AEB08.30n]] != "*") {
      if(baci$Significant[baci$ID == AEB08.40[AEB08.40n]] != "*") {
        baci$Next.Step[baci$ID == AEB08.30[AEB08.30n]]<-baci$Next.Step[baci$ID == AEB08.40[AEB08.40n]]<- "Go to BxC"
        baci$Interpretation[baci$ID == AEB08.30[AEB08.30n]]<-baci$Interpretation[baci$ID == AEB08.40[AEB08.40n]]<- "No medium-term and local-scale interactions"
      }
      # P(B)xC = *; P(B)xI = ns
      if(baci$Significant[baci$ID == AEB08.40[AEB08.40n]] == "*") {
          baci$Next.Step[baci$ID == AEB08.30[AEB08.30n]]<-baci$Next.Step[baci$ID == AEB08.40[AEB08.40n]]<- "END"
          baci$Interpretation[baci$ID == AEB08.30[AEB08.30n]]<-baci$Interpretation[baci$ID == AEB08.40[AEB08.40n]]<- "Medium-term and local-scale interactions: temporal trends are not associated with the disturbed location"
      }
    }
      # P(B)xC = ns; P(B)xI = *
      if(baci$Significant[baci$ID == AEB08.30[AEB08.30n]] == "*") {
        baci$Next.Step[baci$ID == AEB08.30[AEB08.30n]]<-baci$Next.Step[baci$ID == AEB08.40[AEB08.40n]]<- "Go to P(Bef)xI (2-tail test)"
        baci$Interpretation[baci$ID == AEB08.30[AEB08.30n]]<-baci$Interpretation[baci$ID == AEB08.40[AEB08.40n]]<- "Medium-term and local-scale interactions"
        # 2-tail test
        baci$Next.Step[baci$ID == "08AEB.1030"]<-baci$Next.Step[baci$ID == "08AEB.1040"]<- "END"
        ifelse(baci$Significant[baci$ID == "08AEB.1030"] == "*" & baci$Significant[baci$ID =="08AEB.1040"] != "*",
               baci$Interpretation[baci$ID == "08AEB.1030"]<-baci$Interpretation[baci$ID == "08AEB.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from Medium-term and local-scale: impact location differs from before and after the disturbance occurs",
               baci$Interpretation[baci$ID == "08AEB.1030"]<-baci$Interpretation[baci$ID == "08AEB.1040"]<- "No medium-term and local-scale impact detected: the changes are not associated with the disturbed location")
      }

#  # P(B)xC = ns; P(B)xI = ns
#    if(baci$Significant[baci$ID == AEB08.40[AEB08.40n]] != "*") {
#      if(baci$Significant[baci$ID == AEB08.30[AEB08.30n]] != "*") {
#        baci$Next.Step[baci$ID == AEB08.30[AEB08.30n]]<-baci$Next.Step[baci$ID == AEB08.40[AEB08.40n]]<- "Go to BxC"
#        baci$Interpretation[baci$ID == AEB08.30[AEB08.30n]]<-baci$Interpretation[baci$ID == AEB08.40[AEB08.40n]]<- "No medium-term and local-scale interactions"
#      }
#      # P(B)xC = ns; P(B)xI = *
#      if(baci$Significant[baci$ID == AEB08.30[AEB08.30n]] == "*") {
#        baci$Next.Step[baci$ID == AEB08.30[AEB08.30n]]<-baci$Next.Step[baci$ID == AEB08.40[AEB08.40n]]<- "Go to P(Bef)xI (2-tail test) and BxC"
#        baci$Interpretation[baci$ID == AEB08.30[AEB08.30n]]<-baci$Interpretation[baci$ID == AEB08.40[AEB08.40n]]<- "Medium-term and local-scale interactions"
#        # 2-tail test
#        baci$Next.Step[baci$ID == "08AEB.1030"]<-baci$Next.Step[baci$ID == "08AEB.1040"]<- "END"
#        ifelse(baci$Significant[baci$ID == "08AEB.1030"] == "*" & baci$Significant[baci$ID =="08AEB.1040"] != "*",
#               baci$Interpretation[baci$ID == "08AEB.1030"]<-baci$Interpretation[baci$ID == "08AEB.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from Medium-term and local-scale: impact location differs from before and after the disturbance occurs",
#               baci$Interpretation[baci$ID == "08AEB.1030"]<-baci$Interpretation[baci$ID == "08AEB.1040"]<- "No medium-term and local-scale impact detected: the changes are not associated with the disturbed location")
#      }
#    }
#    # P(B)xC = *; P(B)xI = ns
#    if(baci$Significant[baci$ID == AEB08.40[AEB08.40n]] == "*") {
#      if(baci$Significant[baci$ID == AEB08.30[AEB08.30n]] != "*") {
#        baci$Next.Step[baci$ID == AEB08.30[AEB08.30n]]<-baci$Next.Step[baci$ID == AEB08.40[AEB08.40n]]<- "Go to BxC"
#        baci$Interpretation[baci$ID == AEB08.30[AEB08.30n]]<-baci$Interpretation[baci$ID == AEB08.40[AEB08.40n]]<- "Medium-term and local-scale interactions: temporal trends are not associated with the disturbed location"
#      }
#      # P(B)xC = *; P(B)xI = *
#      if(baci$Significant[baci$ID == AEB08.30[AEB08.30n]] == "*") {
#        baci$Next.Step[baci$ID == AEB08.30[AEB08.30n]]<-baci$Next.Step[baci$ID == AEB08.40[AEB08.40n]]<- "Go to BxC"
#        baci$Interpretation[baci$ID == AEB08.30[AEB08.30n]]<-baci$Interpretation[baci$ID == AEB08.40[AEB08.40n]]<- "Medium-term and local-scale interactions: widespread temporal trends are associated with both, control and disturbed locations"
#      }
#    }

  }


  if(aav.model == "AEBGC" | aav.model == "AEBC" | aav.model == "AEGC" | aav.model == "AEC") {

    # Location vc Before and After B x L
    AE03.40<- c("03AE.0040","03AE.0041","03AE.0042","03AE.0043","03AE.0044","03AE.0045","03AE.0046","03AE.0047","03AE.0048","03AE.0049")
    AE03.40n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0040"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0041"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0042"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0043"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0044"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0045"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0046"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0047"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0048"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0049"])) == F)

    AE03.30<- c("03AE.0030","03AE.0031","03AE.0032","03AE.0033","03AE.0034","03AE.0035","03AE.0036")
    AE03.30n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0030"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0031"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0032"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0033"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0034"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0035"]),
                      is.na(baci$Post.Hoc.Pooling[baci$ID == "03AE.0036"])) == F)

    # BxI = ns
    # BxC = ns; BxI = ns
    if(baci$Significant[baci$ID == AE03.30[AE03.30n]] != "*") {
      if(baci$Significant[baci$ID == AE03.40[AE03.40n]] != "*") {
        baci$Next.Step[baci$ID == AE03.30[AE03.30n]]<-baci$Next.Step[baci$ID == AE03.40[AE03.40n]]<- "END"
        baci$Interpretation[baci$ID == AE03.30[AE03.30n]]<-baci$Interpretation[baci$ID == AE03.40[AE03.40n]]<- "No long-term and local-scale interactions"
      }
      # BxC = *; BxI = ns
      if(baci$Significant[baci$ID == AE03.40[AE03.40n]] == "*") {
          baci$Next.Step[baci$ID == AE03.30[AE03.30n]]<-baci$Next.Step[baci$ID == AE03.40[AE03.40n]]<- "END"
          baci$Interpretation[baci$ID == AE03.30[AE03.30n]]<-baci$Interpretation[baci$ID == AE03.40[AE03.40n]]<- "Long-term and local-scale interactions: temporal trends are not associated with the disturbed location"
        }
    }
    # BxI = *
    # BxC = ns; BxI = *  &  BxC = *; BxI = *
    if(baci$Significant[baci$ID == AE03.30[AE03.30n]] == "*") {
        baci$Next.Step[baci$ID == AE03.30[AE03.30n]]<-baci$Next.Step[baci$ID == AE03.40[AE03.40n]]<- "END"
        baci$Interpretation[baci$ID == AE03.30[AE03.30n]]<-baci$Interpretation[baci$ID == AE03.40[AE03.40n]]<- "PRESS DISTURBANCE - IMPACT DETECTED from long-term and local-scale: Disturbed location differs from controls locations"
    }

#    # BxC = ns; BxI = ns
#    if(baci$Significant[baci$ID == AE03.40[AE03.40n]] != "*") {
#      if(baci$Significant[baci$ID == AE03.30[AE03.30n]] != "*") {
#        baci$Next.Step[baci$ID == AE03.30[AE03.30n]]<-baci$Next.Step[baci$ID == AE03.40[AE03.40n]]<- "END"
#        baci$Interpretation[baci$ID == AE03.30[AE03.30n]]<-baci$Interpretation[baci$ID == AE03.40[AE03.40n]]<- "No long-term and local-scale interactions"
#      }
#      # BxC = ns; BxI = *
#      if(baci$Significant[baci$ID == AE03.30[AE03.30n]] == "*") {
#        baci$Next.Step[baci$ID == AE03.30[AE03.30n]]<-baci$Next.Step[baci$ID == AE03.40[AE03.40n]]<- "END"
#        baci$Interpretation[baci$ID == AE03.30[AE03.30n]]<-baci$Interpretation[baci$ID == AE03.40[AE03.40n]]<- "PRESS DISTURBANCE - IMPACT DETECTED from long-term and local-scale: Disturbed location differs from controls locations"
#      }
#    }
#    # BxC = *; BxI = ns
#    if(baci$Significant[baci$ID == AE03.40[AE03.40n]] == "*") {
#      if(baci$Significant[baci$ID == AE03.30[AE03.30n]] != "*") {
#        baci$Next.Step[baci$ID == AE03.30[AE03.30n]]<-baci$Next.Step[baci$ID == AE03.40[AE03.40n]]<- "END"
#        baci$Interpretation[baci$ID == AE03.30[AE03.30n]]<-baci$Interpretation[baci$ID == AE03.40[AE03.40n]]<- "Long-term and local-scale interactions: temporal trends are not associated with the disturbed location"
#      }
#      # BxC = *; BxI = *
#      if(baci$Significant[baci$ID == AE03.30[AE03.30n]] == "*") {
#        baci$Next.Step[baci$ID == AE03.30[AE03.30n]]<-baci$Next.Step[baci$ID == AE03.40[AE03.40n]]<- "END"
#        baci$Interpretation[baci$ID == AE03.30[AE03.30n]]<-baci$Interpretation[baci$ID == AE03.40[AE03.40n]]<- "Long-term and local-scale interactions: widespread temporal trends are associated with both, control and disturbed locations"
#      }
#    }
  }



  if(aav.model == "AEC") {

    # Location vc Time T(B) x L
    AEC10.30<- c("10AEC.2030","10AEC.2031")
    AEC10.30n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEC.2030"]),
                       is.na(baci$Post.Hoc.Pooling[baci$ID == "10AEC.2031"])) == F)

    # T(B)xI = ns
    # T(B)xC = ns; T(B)xI = ns
    if(baci$Significant[baci$ID == AEC10.30[AEC10.30n]] != "*") {
      if(baci$Significant[baci$ID == "10AEC.2040"] != "*") {
        baci$Next.Step[baci$ID == AEC10.30[AEC10.30n]]<-baci$Next.Step[baci$ID == "10AEC.2040"]<- "Go to BxC"
        baci$Interpretation[baci$ID == AEC10.30[AEC10.30n]]<-baci$Interpretation[baci$ID == "10AEC.2040"]<- "No short-term and local-scale interactions"
      }
      # T(B)xC = *; T(B)xI = ns
      if(baci$Significant[baci$ID == "10AEC.2040"] == "*") {
          baci$Next.Step[baci$ID == AEC10.30[AEC10.30n]]<-baci$Next.Step[baci$ID == "10AEC.2040"]<- "END"
          baci$Interpretation[baci$ID == AEC10.30[AEC10.30n]]<-baci$Interpretation[baci$ID == "10AEC.2040"]<- "Short-term and local-scale interactions: temporal trends are not associated with the disturbed location"
      }
    }
      # T(B)xI = *
      # T(B)xC = ns; T(B)xI = *  &  T(B)xC = *; T(B)xI = *
      if(baci$Significant[baci$ID == AEC10.30[AEC10.30n]] == "*") {
        baci$Next.Step[baci$ID == AEC10.30[AEC10.30n]]<-baci$Next.Step[baci$ID == "10AEC.2040"]<- "Go to T(Bef)xI (2-tail test)"
        baci$Interpretation[baci$ID == AEC10.30[AEC10.30n]]<-baci$Interpretation[baci$ID == "10AEC.2040"]<- "Short-term and local-scale interactions"
        # 2-tail test
        baci$Next.Step[baci$ID == "10AEC.1030"]<-baci$Next.Step[baci$ID == "10AEC.1040"]<- "END"
        ifelse(baci$Significant[baci$ID == "10AEC.1030"] == "*" & baci$Significant[baci$ID =="10AEC.1040"] != "*",
               baci$Interpretation[baci$ID == "10AEC.1030"]<-baci$Interpretation[baci$ID == "10AEC.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and local-scale: impact location differs from before and after the disturbance occurs",
               baci$Interpretation[baci$ID == "10AEC.1030"]<-baci$Interpretation[baci$ID == "10AEC.1040"]<- "No short-term and local-scale impact detected: the changes are not associated with the disturbed location")
        #"Short-term IMPACT DETECTED: impact sites differs from before and after the disturbance occurs"
        #"Short-term IMPACT DETECTED: control sites do not differs from before and after the disturbance"
        #### Therefore effect is specific to the impacted sites and coincident with the start of the disturbance. IMPACT DETECTED. ###
      }

#    # T(B)xC = ns; T(B)xI = ns
#    if(baci$Significant[baci$ID == "10AEC.2040"] != "*") {
#      if(baci$Significant[baci$ID == AEC10.30[AEC10.30n]] != "*") {
#        baci$Next.Step[baci$ID == AEC10.30[AEC10.30n]]<-baci$Next.Step[baci$ID == "10AEC.2040"]<- "Go to BxC"
#        baci$Interpretation[baci$ID == AEC10.30[AEC10.30n]]<-baci$Interpretation[baci$ID == "10AEC.2040"]<- "No short-term and local-scale interactions"
#      }
#      # T(B)xC = ns; T(B)xI = *
#      if(baci$Significant[baci$ID == AEC10.30[AEC10.30n]] == "*") {
#        baci$Next.Step[baci$ID == AEC10.30[AEC10.30n]]<-baci$Next.Step[baci$ID == "10AEC.2040"]<- "Go to T(Bef)xI (2-tail test) and BxC"
#        baci$Interpretation[baci$ID == AEC10.30[AEC10.30n]]<-baci$Interpretation[baci$ID == "10AEC.2040"]<- "Short-term and local-scale interactions"
#        # 2-tail test
#        baci$Next.Step[baci$ID == "10AEC.1030"]<-baci$Next.Step[baci$ID == "10AEC.1040"]<- "END"
#        ifelse(baci$Significant[baci$ID == "10AEC.1030"] == "*" & baci$Significant[baci$ID =="10AEC.1040"] != "*",
#               baci$Interpretation[baci$ID == "10AEC.1030"]<-baci$Interpretation[baci$ID == "10AEC.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and local-scale: impact location differs from before and after the disturbance occurs",
#               baci$Interpretation[baci$ID == "10AEC.1030"]<-baci$Interpretation[baci$ID == "10AEC.1040"]<- "No short-term and local-scale impact detected: the changes are not associated with the disturbed location")
#        #"Short-term IMPACT DETECTED: impact sites differs from before and after the disturbance occurs"
#        #"Short-term IMPACT DETECTED: control sites do not differs from before and after the disturbance"
#        #### Therefore effect is specific to the impacted sites and coincident with the start of the disturbance. IMPACT DETECTED. ###
#      }
#    }
#    # T(B)xC = *; T(B)xI = ns
#    if(baci$Significant[baci$ID == "10AEC.2040"] == "*") {
#      if(baci$Significant[baci$ID == AEC10.30[AEC10.30n]] != "*") {
#        baci$Next.Step[baci$ID == AEC10.30[AEC10.30n]]<-baci$Next.Step[baci$ID == "10AEC.2040"]<- "Go to BxC"
#        baci$Interpretation[baci$ID == AEC10.30[AEC10.30n]]<-baci$Interpretation[baci$ID == "10AEC.2040"]<- "Short-term and local-scale interactions: temporal trends are not associated with the disturbed location"
#      }
#      # T(B)xC = *; T(B)xI = *
#      if(baci$Significant[baci$ID == AEC10.30[AEC10.30n]] == "*") {
#        baci$Next.Step[baci$ID == AEC10.30[AEC10.30n]]<-baci$Next.Step[baci$ID == "10AEC.2040"]<- "Go to BxC"
#        baci$Interpretation[baci$ID == AEC10.30[AEC10.30n]]<-baci$Interpretation[baci$ID == "10AEC.2040"]<- "Short-term and local-scale interactions: widespread temporal trends are associated with both, control and disturbed locations"
#      }
#    }
  }


  if(aav.model == "AEGC") {

    # site vc Time T(B) x S(L)
    AEGC11.30<- c("11AEGC.2030","11AEGC.2031")
    AEGC11.30n<-which(c(is.na(baci$Post.Hoc.Pooling[baci$ID == "11AEGC.2030"]),
                        is.na(baci$Post.Hoc.Pooling[baci$ID == "11AEGC.2031"])) == F)

    # T(B)xS(I) = ns
    # T(B)xS(C) = ns; T(B)xS(I) = ns
    if(baci$Significant[baci$ID == AEGC11.30[AEGC11.30n]] != "*") {
      if(baci$Significant[baci$ID == "11AEGC.2040"] != "*") {
        baci$Next.Step[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Next.Step[baci$ID == "11AEGC.2040"]<- "Go to BxS(C) and T(Aft)xC"
        baci$Interpretation[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Interpretation[baci$ID == "11AEGC.2040"]<- "No short-term and small-scale interactions"
      }
      # T(B)xS(C) = *; T(B)xS(I() = ns
      if(baci$Significant[baci$ID == "11AEGC.2040"] == "*") {
          baci$Next.Step[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Next.Step[baci$ID == "11AEGC.2040"]<- "END"
          baci$Interpretation[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Interpretation[baci$ID == "11AEGC.2040"]<- "Short-term and small-scale interactions: temporal trends are not associated with the disturbed sites"
      }
    }
    # T(B)xS(I) = *
    # T(B)xS(C) = ns; T(B)xS(I) = *
    if(baci$Significant[baci$ID == AEGC11.30[AEGC11.30n]] == "*") {
      baci$Next.Step[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Next.Step[baci$ID == "11AEGC.2040"]<- "Go to T(B)xS(I) (2-tail test)"
      baci$Interpretation[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Interpretation[baci$ID == "11AEGC.2040"]<- "Short-term and small-scale interactions"
      # 2-tail test
      baci$Next.Step[baci$ID == "11AEGC.1030"]<-baci$Next.Step[baci$ID == "11AEGC.1040"]<- "END"
      ifelse(baci$Significant[baci$ID == "11AEGC.1030"] == "*" & baci$Significant[baci$ID =="11AEGC.1040"] != "*",
             baci$Interpretation[baci$ID == "11AEGC.1030"]<-baci$Interpretation[baci$ID == "11AEGC.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and small-scale: impact sites differs from before and after the disturbance occurs",
             baci$Interpretation[baci$ID == "11AEGC.1030"]<-baci$Interpretation[baci$ID == "11AEGC.1040"]<- "No short-term impact detected: the changes are not associated with the disturbed sites")
    }

#    # T(B)xS(C) = ns; T(B)xS(I) = ns
#    if(baci$Significant[baci$ID == "11AEGC.2040"] != "*") {
#      if(baci$Significant[baci$ID == AEGC11.30[AEGC11.30n]] != "*") {
#        baci$Next.Step[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Next.Step[baci$ID == "11AEGC.2040"]<- "Go to BxS(C) and T(Aft)xC"
#        baci$Interpretation[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Interpretation[baci$ID == "11AEGC.2040"]<- "No short-term and small-scale interactions"
#      }
#      # T(B)xS(C) = ns; T(B)xS(I) = *
#      if(baci$Significant[baci$ID == AEGC11.30[AEGC11.30n]] == "*") {
#        baci$Next.Step[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Next.Step[baci$ID == "11AEGC.2040"]<- "Go to T(B)xS(I) (2-tail test), BxS(C) and T(Aft)xC"
#        baci$Interpretation[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Interpretation[baci$ID == "11AEGC.2040"]<- "Short-term and small-scale interactions"
#        # 2-tail test
#        baci$Next.Step[baci$ID == "11AEGC.1030"]<-baci$Next.Step[baci$ID == "11AEGC.1040"]<- "END"
#        ifelse(baci$Significant[baci$ID == "11AEGC.1030"] == "*" & baci$Significant[baci$ID =="11AEGC.1040"] != "*",
#               baci$Interpretation[baci$ID == "11AEGC.1030"]<-baci$Interpretation[baci$ID == "11AEGC.1040"]<- "PULSE DISTURBANCE - IMPACT DETECTED from short-term and small-scale: impact sites differs from before and after the disturbance occurs",
#               baci$Interpretation[baci$ID == "11AEGC.1030"]<-baci$Interpretation[baci$ID == "11AEGC.1040"]<- "No short-term impact detected: the changes are not associated with the disturbed sites")
#      }
#    }
#    # T(B)xS(C) = *; T(B)xS(I() = ns
#    if(baci$Significant[baci$ID == "11AEGC.2040"] == "*") {
#      if(baci$Significant[baci$ID == AEGC11.30[AEGC11.30n]] != "*") {
#        baci$Next.Step[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Next.Step[baci$ID == "11AEGC.2040"]<- "Go to BxS(C) and T(Aft)xC"
#        baci$Interpretation[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Interpretation[baci$ID == "11AEGC.2040"]<- "Short-term and small-scale interactions: temporal trends are not associated with the disturbed sites"
#      }
#      # T(B)xS(C) = *; T(B)xS(I) = *
#      if(baci$Significant[baci$ID == AEGC11.30[AEGC11.30n]] == "*") {
#        baci$Next.Step[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Next.Step[baci$ID == "11AEGC.2040"]<- "Go to BxS(C) and T(Aft)xC"
#        baci$Interpretation[baci$ID == AEGC11.30[AEGC11.30n]]<-baci$Interpretation[baci$ID == "11AEGC.2040"]<- "Short-term and small-scale interactions: widespread temporal trends are associated with both, control and disturbed sites"
#      }
#    }
  }


  ###---------------------------------------------------------------  Final.Table

  baci$Final.Table<-ifelse(is.na(baci$P.Value.Upper.Tail) == T, "Yes",
                           ifelse(!is.na(baci$Post.Hoc.Pooling) == T, "Yes",
                                  ifelse(!is.na(baci$P.Value.Lower.Tail) &
                                           !is.na(baci$Interpretation) == T, "Yes","No")))

  ###-------------------------------------------------------------------------

  return(baci)
}

#----------------------------------------------------------------------- E N D
#----------------------------------------------------------------------- E N D
