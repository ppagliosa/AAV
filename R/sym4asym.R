
#' @title Symmetrical ANOVAs to Run an Asymmetrical ANOVA
#'
#' @description
#' Calculate the symmetrical ANOVAs to run an asymmetrical ANOVA
#' beyond BACI or post-impact design.
#'
#' Asymmetrical ANOVAs are essentially multi-factorial.
#' Beyond BACI design can be used for models with 2 or 3 temporal factors
#' and 1 or 2 spatial factors. Post-impact design with only one disturbed location
#' and multiple control locations can be used for models with 1 or 2 temporal factors
#' and 1 or 2 spatial factors. Post-impact design with more than one place with
#' one disturbed location and multiple control locations can be used for models
#' with 2 or 3 spatial factors.
#'
#' @param data a data frame with temporal factors, spatial factors and variable(s).
#' @param aav.design \code{baci} for before/after-control/impact design or \code{aci} for post-impact design.
#' @param n.ftemp number of temporal factors.
#' @param n.fspac number of spatial factors.
#' @param names.impact name of the asymmetrical level of the Control/Impact factor (e.g. "impact").
#' @param names.before only for \code{baci} design. Name of the before level of the Before/After factor (e.g. "before").
#' @param aav.file A file path to save the AAV output files. Default \code{NULL} uses current name and working directory.
#'
#' @details
#' The construction of asymmetrical analyses of variance consists of combining
#' the Sum of Squares values from separate (symmetrical) ANOVAs.
#' In beyond BACI designs, four fully orthogonal ANOVAs need to be done: on all
#' the data, on "control" data, on "before" data, and on "control" and "before" data.
#' In post-impact designs with only one disturbed location two fully orthogonal
#' ANOVAs need to be done: on all the data and on "control" data.
#' In post-impact designs with more than one place with one disturbed location
#' two fully orthogonal and one fully nested ANOVAs need to be done:
#' on all data (orthogonal) and on "control" data (orthogonal and nested).
#' \code{sym4asy} use \code{\link[GAD]{GAD}} (Sandrini-Neto et al. 2024) to perform symmetrical ANOVAs.
#'
#' The \code{sym4asym} analysis fits mixed ANOVA designs with a specific
#' combination of orthogonal/nested and fixed/random factors:
#' The  beyond BACI design runs for five-factors (3 temporal and 2 spatial),
#' four-factors (3 temporal and 1 spatial or 2 temporal and 2 spatial), and
#' three-factors (2 temporal and 1 spatial).
#' The post-impact design with only one disturbed location runs for
#' four-factors (2 temporal and 2 spatial), three-factors (2 temporal and 1 spatial or
#' 1 temporal and 2 spatial), and two-factors (1 temporal and 1 spatial).
#' The post-impact design with more than one place with one disturbed location runs for
#' three-factors (0 temporal and 3 spatial) and two-factors (0 temporal and 2 spatial).
#'
#' In all models the Location factor is the asymmetric one (e.g., one "impact"
#' group and two or more "control" groups), meaning there are different numbers
#' of levels of this nested factor within each different level of its upper-level
#' factor (e.g., Before vs After, Period, Time or Place factors).
#'
#' In all models the temporal and spatial factors are orthogonal.
#' Within each factor type (temporal or spatial) the factors are hierarchically nested.
#' In the beyond BACI and post-impact design with only one disturbed location models,
#' the first factor and the Location factor are fixed and all others are random.
#' In the post-impact design with more than one place with one disturbed location,
#' the first factor is fixed and all other factors are random.
#'
#' The data frame must follow a specific structure:
#' 1- first the temporal factors, then the spatial factors, and then the variable(s);
#' 2- each factor type must appear in a hierarchical order.
#'
#' For beyond BACI design the order for temporal factors is Before_after, Period and Time (3 factors)
#' or Before_after and Time (2 factors).
#' The order for spatial factors is Location and Site (2 factors)
#' or just Location (1 factor).
#' For post-impact design with only one disturbed location and multiple control locations,
#' the order for temporal factors must be Period and Time (2 factors) or just Time (1 factor), and
#' the order for spatial factors must be Location and Site (2 factors) or just Location (1 factor).
#' For post-impact design with more than one place with one disturbed location and multiple control locations,
#' the order for spatial factors must be Place, Location and Site (3 factors) or Place and Location (2 factors).
#'
#' If the names of the factors are different from the above names they will be changed.
#'
#' When the data frame have more than one variable the univariate ANOVAs are run
#' for each one separated.
#'
#' @returns
#' \code{Excel} files containing the symmetrical ANOVA tables used to
#' run the asymmetrical ANOVAs.
#' When analyzing more than one variable at once, each sheet have one ANOVA table.
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
#' Underwood, A.J. 1997. \emph{Experiments in Ecology: Their Logical Design and Interpretation Using
#' Analysis of Variance}. Cambridge University Press, Cambridge. \doi{10.1017/CBO9780511806407}.
#'
#' @seealso
#' \code{\link[AAV]{aav}},
#' \code{\link[AAV]{asym.aci}},
#' \code{\link[AAV]{asym.baci}},
#' \code{\link[GAD]{gad}}
#'
#' @examples
#' # Example for BACI design containing 3 temporal factors and 2 spatial factor:
#' data(baci.3t2s)
#' ex.baci.3t2s<-sym4asym(data = baci.3t2s, aav.design = "baci", n.ftemp = 3, n.fspac = 2,
#' names.impact = "impact", names.before = "before")
#' ex.baci.3t2s
#'
#' # Example for post-impact design with only one disturbed location
#' # containing 2 temporal factors and 2 spatial factor:
#' data(aci.2t2s)
#' ex.aci.2t2s<-sym4asym(data = aci.2t2s, aav.design = "aci", n.ftemp = 2, n.fspac = 2,
#'     names.impact = "impact")
#' ex.aci.2t2s
#'
#' # Example for post-impact design with more than one place with one
#' # disturbed location containing 0 temporal factors and 3 spatial factor:
#'  data(aci.0t3s)
#' ex.aci.0t3s<-sym4asym(data = aci.0t3s, aav.design = "aci", n.ftemp = 0, n.fspac = 3,
#'     names.impact = "impact")
#' ex.aci.0t3s
#'
#' @import GAD
#' @import openxlsx
#'
#' @importFrom stats lm
#'
#' @export

sym4asym<- function(data, aav.design, n.ftemp, n.fspac, names.impact,
                    names.before = NULL, aav.file = NULL) {

  ###--------------------------------------------------------------------  A C I
  ###--------------------------------------------------------------------  A C I
  ###--------------------------------------------------------------------  A C I

  if(aav.design %in% c("baci","aci") != T) stop("The 'aav.design' object does not exist. Check your command line.")
  if(names.impact %in% unique(data$location) != T) stop("The 'names.impact' object does not exist. Check your command line.")

  # lista para armazenar as ANOVAS simetricas
  symANOVA<-list()

  if(aav.design == "aci") {
    # Stops
    if(n.ftemp >= 3 | n.fspac == 0 | n.fspac >= 4) stop("Select different values for 'n.ftemp' or 'n.fspac'.")

    # Factors names
    if(n.ftemp == 2 & n.fspac == 2) {
      names(data)[1:4]<-c("period", "time", "location", "site")
      aav.model <-"BEGC"
      }
    if(n.ftemp == 2 & n.fspac == 1) {
      names(data)[1:3]<-c("period", "time", "location")
      aav.model <-"BEC"
      }
    if(n.ftemp == 1 & n.fspac == 2) {
      names(data)[1:3]<-c("time", "location", "site")
      aav.model <-"CEG"
      }
    if(n.ftemp == 1 & n.fspac == 1) {
      names(data)[1:2]<-c("time", "location")
      aav.model <-"CE"
      }
    if(n.ftemp == 0 & n.fspac == 2) {
      names(data)[1:2]<-c("place", "location")
      aav.model <-"DE"
      }
    if(n.ftemp == 0 & n.fspac == 3) {
      names(data)[1:3]<-c("place", "location", "site")
      aav.model <-"DEG"
    }

    ###------------------------------------------------- fULL = all data; full model
    B<-C<-D<-E<-G<-NULL
    ifelse("B" %in% substring(aav.model,1:4,1:4) == T, B <- as.fixed(data$period), B <-NA)
    if("C" %in% substring(aav.model,1:4,1:4) != T) C <-NA
    if("C" %in% substring(aav.model,1:4,1:4) == T & "B" %in% substring(aav.model,1:4,1:4) == T) C <- as.random(data$time)
    if("C" %in% substring(aav.model,1:4,1:4) == T & "B" %in% substring(aav.model,1:4,1:4) != T) C <- as.fixed(data$time)
    if("D" %in% substring(aav.model,1:4,1:4) != T) {
      D <-NA
      E <- as.fixed(data$location)
    }
    if("D" %in% substring(aav.model,1:4,1:4) == T) {
      D <- as.fixed(data$place)
      E <- as.random(data$location)
    }
    ifelse("G" %in% substring(aav.model,1:4,1:4) == T, G <- as.random(data$site), G <-NA)

    md <- data.frame(data[,(n.ftemp+n.fspac+1):length(data)])
    if(length(md) == 1) names(md)<-colnames(data[length(data)])

    # Additional Stops
    if(0 %in% colSums(md)) stop("At least one of the variable responses has only 0 values. Check your dataset.")
    md.ca<-data.frame(data[data$location != names.impact,(n.ftemp+n.fspac+1):length(data)])
    if(length(md.ca) == 1) names(md.ca)<-colnames(data[length(data)])
    if(0 %in% colSums(md.ca)) stop("At least one variable response has only 0 values in the 'Controls' location. Check your dataset.")

    # Running the symmetrical ANOVAs
    Ltemp<-list()
    Lfull.d<-list()
    for (i in 1:length(md)) {
      if(aav.model == "BEGC") {
        model.fd <- lm(md[,i] ~ B + E + B*E + G%in%E + C%in%B +
                         B*(G%in%E) + E*(C%in%B) +
                         (C%in%B)*(G%in%E), data = md)
      }
      if(aav.model == "BEC") {
        model.fd <- lm(md[,i] ~ B + E + B*E + C%in%B + E*(C%in%B), data = md)
      }
      if(aav.model == "CEG") {
        model.fd <- lm(md[,i] ~ C + E + C*E + G%in%E + C*(G%in%E), data = md)
      }
      if(aav.model == "CE") {
        model.fd <- lm(md[,i] ~ C + E + C*E, data = md)
      }
      if(aav.model == "DE") {
        model.fd <- lm(md[,i] ~ D + E + D*E, data = md)
      }
      if(aav.model == "DEG") {
        model.fd <- lm(md[,i] ~ D + E + D*E + G%in%E + D*(G%in%E), data = md)
      }
      #estimates(model.fd)
      Ltemp[[i]]<-gad(model.fd)
      Lfull.d[[i]]<-Ltemp[[i]][[1]]
    }
    names(Lfull.d)<-colnames(md)
    symANOVA[[1]]<-Lfull.d
    names(symANOVA)[[1]]<-"full_design"

    ###--------------------------- soC = excluding the disturbed ('Impact') location
    B<-C<-D<-E<-G<-NULL
    dataC<-data[data$location != names.impact,]

    ifelse("B" %in% substring(aav.model,1:4,1:4) == T, B <- as.fixed(dataC$period), B <-NA)
    if("C" %in% substring(aav.model,1:4,1:4) != T) C <-NA
    if("C" %in% substring(aav.model,1:4,1:4) == T & "B" %in% substring(aav.model,1:4,1:4) == T) C <- as.random(dataC$time)
    if("C" %in% substring(aav.model,1:4,1:4) == T & "B" %in% substring(aav.model,1:4,1:4) != T) C <- as.fixed(dataC$time)
    if("D" %in% substring(aav.model,1:4,1:4) != T) {
      D <-NA
      E <- as.fixed(dataC$location)
    }
    if("D" %in% substring(aav.model,1:4,1:4) == T) {
      D <- as.fixed(dataC$place)
      E <- as.random(dataC$location)
    }
    ifelse("G" %in% substring(aav.model,1:4,1:4) == T, G <- as.random(dataC$site), G <-NA)

    mdC <- data.frame(dataC[,(n.ftemp+n.fspac+1):length(dataC)])
    if(length(mdC) == 1) {names(mdC)<-colnames(dataC[length(dataC)])}

    Ltemp<-list()
    LC.d<-list()
    LCNtemp<-list()
    LCN.d<-list()
    for (i in 1:length(mdC)) {
      if(aav.model == "BEGC") {
        modelC <- lm(mdC[,i] ~ B + E + B*E + G%in%E + C%in%B +
                       B*(G%in%E) + E*(C%in%B) +
                       (C%in%B)*(G%in%E), data = mdC)
      }
      if(aav.model == "BEC") {
        modelC <- lm(mdC[,i] ~ B + E + B*E + C%in%B + E*(C%in%B), data = mdC)
      }
      if(aav.model == "CEG") {
        modelC <- lm(mdC[,i] ~ C + E + C*E + G%in%E + C*(G%in%E), data = mdC)
      }
      if(aav.model == "CE") {
        modelC <- lm(mdC[,i] ~ C + E + C*E, data = mdC)
      }
      if(aav.model == "DE") {
        modelC <- lm(mdC[,i] ~ D + E + D*E, data = mdC)
        modelCN <- lm(mdC[,i] ~ D + E%in%D, data = mdC)
      }
      if(aav.model == "DEG") {
        modelC <- lm(mdC[,i] ~ D + E + D*E + G%in%E + D*(G%in%E), data = mdC)
        modelCN <- lm(mdC[,i] ~ D + E%in%D + G%in%E%in%D, data = mdC)
      }
      #estimates(modelC)
      Ltemp[[i]]<-gad(modelC)
      LC.d[[i]]<-Ltemp[[i]][[1]]

      if (exists("modelCN") == T) {
        #estimates(modelCN)
        LCNtemp[[i]]<-gad(modelCN)
        LCN.d[[i]]<-LCNtemp[[i]][[1]]
      }
    }
    names(LC.d)<-colnames(mdC)
    symANOVA[[2]]<-LC.d
    names(symANOVA)[[2]]<-"control_design"

    if (exists("modelCN") == T) {
      names(LCN.d)<-colnames(mdC)
      symANOVA[[3]]<-LCN.d
      names(symANOVA)[[3]]<-"nested_control_design"
    }
  }


  ###--------------------------------------------------------------------  B A C I
  ###--------------------------------------------------------------------  B A C I
  ###--------------------------------------------------------------------  B A C I

  if(aav.design == "baci") {
    # Stops
    if(names.before %in% unique(data$before_after) != T) stop("The 'names.before' object does not exist. Check your command line")
    if(n.ftemp <= 1 | n.ftemp >= 4 | n.fspac == 0 | n.fspac >= 3) stop("Select different values for 'n.ftemp' or 'n.fspac'.")

    # Factors names
    if(n.ftemp == 3 & n.fspac == 2) {
      names(data)[1:5]<-c("before_after", "period", "time", "location", "site")
      aav.model <-"AEBGC"
    }
    if(n.ftemp == 3 & n.fspac == 1) {
      names(data)[1:4]<-c("before_after", "period", "time", "location")
      aav.model <-"AEBC"
    }
    if(n.ftemp == 2 & n.fspac == 2) {
      names(data)[1:4]<-c("before_after", "time", "location", "site")
      aav.model <-"AEGC"
    }
    if(n.ftemp == 2 & n.fspac == 1) {
      names(data)[1:3]<-c("before_after", "time", "location")
      aav.model <-"AEC"
    }

    ###------------------------------------------------- fULL = all data; full model
    A<-B<-C<-E<-G<-NULL
    ifelse("A" %in% substring(aav.model,1:5,1:5) == T, A <- as.fixed(data$before_after), A <-NA)
    ifelse("B" %in% substring(aav.model,1:5,1:5) == T, B <- as.random(data$period), B <-NA)
    ifelse("C" %in% substring(aav.model,1:5,1:5) == T, C <- as.random(data$time), C <-NA)
    ifelse("E" %in% substring(aav.model,1:5,1:5) == T, E <- as.fixed(data$location), E <-NA)
    ifelse("G" %in% substring(aav.model,1:5,1:5) == T, G <- as.random(data$site), G <-NA)

    md <- data.frame(data[,(n.ftemp+n.fspac+1):length(data)])
    if(length(md) == 1) names(md)<-colnames(data[length(data)])

    # Additional Stops
    if(0 %in% colSums(md)) stop("At least one of the variable responses has only 0 values. Check your dataset.")
    md.ca<-data.frame(data[data$location != names.impact & data$before_after != names.before,(n.ftemp+n.fspac+1):length(data)])
    if(length(md.ca) == 1) names(md.ca)<-colnames(data[length(data)])
    if(0 %in% colSums(md.ca)) stop("At least one variable response has only 0 values in the 'Controls' location after the disturbance has occurred. Check your dataset.")

    # Running the symmetrical ANOVAs
    Ltemp<-list()
    Lfull.d<-list()
    for (i in 1:length(md)) {
      if(aav.model == "AEBGC") {
        model.fd <- lm(md[,i] ~ A + E + A*E + B%in%A + G%in%E +
                         C%in%B%in%A + A*(G%in%E) + E*(B%in%A) +
                         (B%in%A)*(G%in%E) + E*(C%in%B%in%A) +
                         (C%in%B%in%A)*(G%in%E), data = md)
      }
      if(aav.model == "AEBC") {
        model.fd <- lm(md[,i] ~ A + E + A*E + B%in%A +
                         E*(B%in%A) + C%in%B%in%A +
                         E*(C%in%B%in%A), data = md)
      }
      if(aav.model == "AEGC") {
        model.fd <- lm(md[,i] ~ A + E + A*E + G%in%E + C%in%A +
                         A*(G%in%E) + E*(C%in%A) +
                         (C%in%A)*(G%in%E), data = md)
      }
      if(aav.model == "AEC") {
        model.fd <- lm(md[,i] ~ A + E + A*E + C%in%A + E*(C%in%A), data = md)
      }
      #estimates(model.fd)
      Ltemp[[i]]<-gad(model.fd)
      Lfull.d[[i]]<-Ltemp[[i]][[1]]
    }
    names(Lfull.d)<-colnames(md)
    symANOVA[[1]]<-Lfull.d
    names(symANOVA)[[1]]<-"full_design"


    ###--------------------------- soC = excluding the disturbed ('Impact') location
    A<-B<-C<-E<-G<-NULL
    dataC<-data[data$location != names.impact,]

    ifelse("A" %in% substring(aav.model,1:5,1:5) == T, A <- as.fixed(dataC$before_after), A <-NA)
    ifelse("B" %in% substring(aav.model,1:5,1:5) == T, B <- as.random(dataC$period), B <-NA)
    ifelse("C" %in% substring(aav.model,1:5,1:5) == T, C <- as.random(dataC$time), C <-NA)
    ifelse("E" %in% substring(aav.model,1:5,1:5) == T, E <- as.fixed(dataC$location), E <-NA) #random pq so tem controle? nao pq pra decompor tem q ser igual ao FULL
    ifelse("G" %in% substring(aav.model,1:5,1:5) == T, G <- as.random(dataC$site), G <-NA)

    mdC <- data.frame(dataC[,(n.ftemp+n.fspac+1):length(dataC)])
    if(length(mdC) == 1) {names(mdC)<-colnames(dataC[length(dataC)])}

    Ltemp<-list()
    LC.d<-list()
    for (i in 1:length(mdC)) {
      if(aav.model == "AEBGC") {
        modelC <- lm(mdC[,i] ~ A + E + A*E + B%in%A + G%in%E +
                       C%in%B%in%A + A*(G%in%E) + E*(B%in%A) +
                       (B%in%A)*(G%in%E) + E*(C%in%B%in%A) +
                       (C%in%B%in%A)*(G%in%E), data = mdC)
      }
      if(aav.model == "AEBC") {
        modelC <- lm(mdC[,i] ~ A + E + A*E + B%in%A +
                       E*(B%in%A) + C%in%B%in%A +
                       E*(C%in%B%in%A), data = mdC)
      }
      if(aav.model == "AEGC") {
        modelC <- lm(mdC[,i] ~ A + E + A*E + G%in%E + C%in%A +
                       A*(G%in%E) + E*(C%in%A) +
                       (C%in%A)*(G%in%E), data = mdC)
      }
      if(aav.model == "AEC") {
        modelC <- lm(mdC[,i] ~ A + E + A*E + C%in%A + E*(C%in%A), data = mdC)
      }
      #estimates(modelC)
      Ltemp[[i]]<-gad(modelC)
      LC.d[[i]]<-Ltemp[[i]][[1]]
    }
    names(LC.d)<-colnames(mdC)
    symANOVA[[2]]<-LC.d
    names(symANOVA)[[2]]<-"control_design"


    ###---------------------------------- soB = all locations before the disturbance
    A<-B<-C<-E<-G<-NULL
    dataB<-data[data$before_after == names.before,]

    #   ifelse("A" %in% substring(aav.model,1:5,1:5) == T, A <- as.fixed(dataB$before_after), A <-NA)
    ifelse("B" %in% substring(aav.model,1:5,1:5) == T, B <- as.random(dataB$period), B <-NA)   #fixed pq e o primeiro? nao pq pra decompor tem q ser igual ao FULL
    ifelse("C" %in% substring(aav.model,1:5,1:5) == T, C <- as.random(dataB$time), C <-NA)     #fixed qndo e o primeiro? nao pq pra decompor tem q ser igual ao FULL
    ifelse("E" %in% substring(aav.model,1:5,1:5) == T, E <- as.fixed(dataB$location), E <-NA)
    ifelse("G" %in% substring(aav.model,1:5,1:5) == T, G <- as.random(dataB$site), G <-NA)

    mdB <- data.frame(dataB[,(n.ftemp+n.fspac+1):length(dataB)])
    if(length(mdB) == 1) {names(mdB)<-colnames(dataB[length(dataB)])}

    Ltemp<-list()
    LB.d<-list()
    for (i in 1:length(mdB)) {
      if(aav.model == "AEBGC") {
        modelB <- lm(mdB[,i] ~ B + E + B*E + G%in%E + C%in%B +
                       B*(G%in%E) + E*(C%in%B) + (C%in%B)*(G%in%E), data = mdB)
      }
      if(aav.model == "AEBC") {
        modelB <- lm(mdB[,i] ~ B + E + B*E + C%in%B + E*(C%in%B), data = mdB)
      }
      if(aav.model == "AEGC") {
        modelB <- lm(mdB[,i] ~ C + E + C*E + G%in%E + C*(G%in%E), data = mdB)
      }
      if(aav.model == "AEC") {
        modelB <- lm(mdB[,i] ~ C + E + C*E, data = mdB)
      }
      #estimates(modelB)
      Ltemp[[i]]<-gad(modelB)
      LB.d[[i]]<-Ltemp[[i]][[1]]
    }
    names(LB.d)<-colnames(mdB)
    symANOVA[[3]]<-LB.d
    names(symANOVA)[[3]]<-"before_design"

    ###---------------- soBC = before the disturbance, excluding the impact location
    A<-B<-C<-E<-G<-NULL
    dataBC<-data[data$before_after == names.before,]
    dataBC<-dataBC[dataBC$location != names.impact,]

    #    ifelse("A" %in% substring(aav.model,1:5,1:5) == T, A <- as.fixed(dataBC$before_after), A <-NA)
    ifelse("B" %in% substring(aav.model,1:5,1:5) == T, B <- as.random(dataBC$period), B <-NA)  #fixed pq e o primeiro? nao pq pra decompor tem q ser igual ao FULL
    ifelse("C" %in% substring(aav.model,1:5,1:5) == T, C <- as.random(dataBC$time), C <-NA)    #fixed qndo e o primeiro? nao pq pra decompor tem q ser igual ao FULL
    ifelse("E" %in% substring(aav.model,1:5,1:5) == T, E <- as.fixed(dataBC$location), E <-NA) #random pq so tem Controles? nao pq pra decompor tem q ser igual ao FULL
    ifelse("G" %in% substring(aav.model,1:5,1:5) == T, G <- as.random(dataBC$site), G <-NA)

    mdBC <- data.frame(dataBC[,(n.ftemp+n.fspac+1):length(dataBC)])
    if(length(mdBC) == 1) {names(mdBC)<-colnames(dataBC[length(dataBC)])}

    Ltemp<-list()
    LBC.d<-list()
    for (i in 1:length(mdBC)) {
      if(aav.model == "AEBGC") {
        modelBC <- lm(mdBC[,i] ~ B + E + B*E + G%in%E + C%in%B +
                        B*(G%in%E) + E*(C%in%B) + (C%in%B)*(G%in%E), data = mdBC)
      }
      if(aav.model == "AEBC") {
        modelBC <- lm(mdBC[,i] ~ B + E + B*E + C%in%B + E*(C%in%B), data = mdBC)
      }
      if(aav.model == "AEGC") {
        modelBC <- lm(mdBC[,i] ~ C + E + C*E + G%in%E + C*(G%in%E), data = mdBC)
      }
      if(aav.model == "AEC") {
        modelBC <- lm(mdBC[,i] ~ C + E + C*E, data = mdBC)
      }
      #estimates(modelBC)
      Ltemp[[i]]<-gad(modelBC)
      LBC.d[[i]]<-Ltemp[[i]][[1]]
    }
    names(LBC.d)<-colnames(mdBC)
    symANOVA[[4]]<-LBC.d
    names(symANOVA)[[4]]<-"control_before_design"
  }
  return(symANOVA)
}
#----------------------------------------------------------------------- E N D
#----------------------------------------------------------------------- E N D



