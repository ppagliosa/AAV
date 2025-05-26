#' @title Dataset For Testing a Beyond BACI Design
#'
#' @docType data
#'
#' @description Dataset with 3 temporal (Before-after, Period, Time) and
#'     2 spatial (Location, Site) factors and 3 response variables simulating
#'     no impact, press disturbance and pulse disturbance.
#'
#' @format A dataframe with 288 rows and 8 collumns.
#'
#' @details
#' * Before-after: a fixed factor with two levels (before, after).
#' * Period: a random factor with two levels (p1, p2) nested in Before-after.
#' * Time: a random factor with two levels (t1, t2) nested in Period.
#' * Location: a random factor with three levels (control1, control2, impact).
#' * Site: a random factor with two levels (s1, s2) nested in Location.
#' * Var.random: response variable with random numbers; no impact simulation.
#' * Var.L.impact: response variable; press impact - a reduction in numbers at the impact location is simulated.
#' * Var.S.impact: response variable; pulse impact - a reduction in numbers at one site of the impact location is simulated.
#'
#' @source Created in-house to serve as an example.
#'
#' @examples
#'  data(baci.3t2s)
#'  baci.3t2s
"baci.3t2s"


#' @title Dataset For testing a Post-Impact Design
#'     With Only One Disturbed Location and Multiple Control Locations
#'
#' @docType data
#'
#' @description Dataset with 2 temporal (Period, Time) and
#'     2 spatial (Location, Site) factors and 3 response variables simulating
#'     no impact, press disturbance and pulse disturbance.
#'
#' @format A dataframe with 288 rows and 7 collumns.
#'
#' @details
#' * Period: a fixed factor with two levels (p1, p2).
#' * Time: a random factor with two levels (t1, t2) nested in Period.
#' * Location: a random factor with three levels (control1, control2, impact).
#' * Site: a random factor with two levels (s1, s2) nested in Location.
#' * Var.random: response variable with random numbers; no impact simulation.
#' * Var.L.impact: response variable; press impact - a reduction in numbers at the impact location is simulated.
#' * Var.S.impact: response variable; pulse impact - a reduction in numbers at one site of the impact location is simulated.
#'
#' @source Created in-house to serve as an example.
#'
#' @examples
#'  data(aci.2t2s)
#'  aci.2t2s
"aci.2t2s"


#' @title Dataset For Testing a Post-Impact Design
#'     With More Than One Place With One Disturbed Location and Multiple
#'     Control Locations
#'
#' @docType data
#'
#' @description Dataset with 0 temporal and 3 spatial (Place, Location, Site) factors
#'     and 3 response variables simulating no impact, press disturbance
#'     and pulse disturbance.
#'
#' @format A dataframe with 288 rows and 6 collumns.
#'
#' @details
#' * Place: a fixed factor with two levels (pl1, pl2).
#' * Location: a random factor with three levels (control1, control2, impact) nested in Places.
#' * Site: a random factor with two levels (s1, s2) nested in Location.
#' * Var.random: response variable with random numbers; no impact simulation.
#' * Var.L.impact: response variable; press impact - a reduction in numbers at the impact location is simulated.
#' * Var.S.impact: response variable; pulse impact - a reduction in numbers at one site of the impact location is simulated.
#'
#' @source Created in-house to serve as an example.
#'
#' @examples
#'  data(aci.0t3s)
#'  aci.0t3s
"aci.0t3s"


#' @title ANOVA Tables for a Beyond BACI Design
#'
#' @docType data
#'
#' @description There are four analysis of variance tables: one with the full dataset;
#'     one with the dataset from the control locations only;
#'     one with the data from before the disturbance occurred only;
#'     and one with the dataset from the control locations only from before the disturbance occurred.
#'     The ANOVAS are from 3 temporal factors (A, B and C) and 2 spatial factors (E and G).
#'     Source of variation codes in AAV are A (Before-after), B (Period), C (Time), D (Place), E (Location), G (Site).
#'
#' @format A list with 4 sub lists.
#'
#' @details
#' The ANOVAS are from 3 temporal factors (A, B and C) and 2 spatial factors (E and G).
#' Source of variation codes in AAV are A (Before-after), B (Period), C (Time), D (Place), E (Location), G (Site).
#' Each table is provided with:
#' * Source.of.variation: Source of variation.
#' * Df: Degree of freedom.
#' * Sum.Sq: Sum of squares.
#' * Mean.Sq: Mean square.
#' * F.value: F-test.
#' * Pr(>F): P-value.
#'
#' @source Created in-house to serve as an example.
#'
#' @examples
#'  data(symANOVA_baci.3t2s)
#'  symANOVA_baci.3t2s
"symANOVA_baci.3t2s"


#' @title ANOVA Tables for a Post-Impact Design With Only One Disturbed Location and Multiple Control Locations
#'
#' @docType data
#'
#' @description There are two analysis of variance tables: one with the full dataset,
#'     and one with only the control locations dataset.
#'     Analysis of variance table from 2 temporal factors (B and C) and 2 spatial factors (E and G).
#'     Source of variation codes in AAV are A (Before-after), B (Period), C (Time), D (Place), E (Location), G (Site).
#'
#' @format A list with 2 sub lists.
#'
#' @details
#' The ANOVAS are from 2 temporal factors (B and C) and 2 spatial factors (E and G).
#' Source of variation codes in AAV are A (Before-after), B (Period), C (Time), D (Place), E (Location), G (Site).
#' Each table is provided with:
#' * Source.of.variation: Source of variation.
#' * Df: Degree of freedom.
#' * Sum.Sq: Sum of squares.
#' * Mean.Sq: Mean square.
#' * F.value: F-test.
#' * Pr(>F): P-value.
#'
#' @source Created in-house to serve as an example.
#'
#' @examples
#'  data(symANOVA_aci.2t2s)
#'  symANOVA_aci.2t2s
"symANOVA_aci.2t2s"


#' @title ANOVA Tables for a Post-Impact Design With More Than One Place With One Disturbed Location and Multiple Control Locations
#'
#' @docType data
#'
#' @description There are three analysis of variance tables: one based on the full dataset (orthogonal),
#'     and two based on the dataset of only the control locations (one is orthogonal and the other is nested).
#'     Analysis of variance table from 0 temporal factors and 3 spatial factors (D, E and G).
#'     Source of variation codes in AAV are A (Before-after), B (Period), C (Time), D (Place), E (Location), G (Site).
#'
#' @format A list with 3 sub lists.
#'
#' @details
#' The ANOVAS are from 0 temporal factors and 3 spatial factors (D, E and G).
#' Source of variation codes in AAV are A (Before-after), B (Period), C (Time), D (Place), E (Location), G (Site).
#' Each table is provided with:
#' * Source.of.variation: Source of variation.
#' * Df: Degree of freedom.
#' * Sum.Sq: Sum of squares.
#' * Mean.Sq: Mean square.
#' * F.value: F-test.
#' * Pr(>F): P-value.
#'
#' @source Created in-house to serve as an example.
#'
#' @examples
#'  data(symANOVA_aci.0t3s)
#'  symANOVA_aci.0t3s
"symANOVA_aci.0t3s"
