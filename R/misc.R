# cohort = 5 <-> age = 3 <-> year = 7
# cohort = 5 <-> age = 2 <-> year = 6 
# cohort = 5 <-> age = 1 <-> year = 5
#
# ......-
# .....-.
# ....-..
#

#' Convert indices between age, cohort, and year.
#' 
#' If exactly one of the arguments is missing, the supplied indices will be used
#'     to calculate the missing index.
#' 
#' @param age The age of interest.
#' @param cohort The cohort of interest.
#' @param year The year of interest.
#' 
#' @return The value of the missing index.
#' 
#' @export
convert_ayc<- function(age, cohort, year) {
    if( missing(age) ) return(.co_y_to_a(cohort, year))
    if( missing(cohort) ) return(.a_y_to_co(age, year))
    if( missing(year) ) return(.a_co_to_y(age, cohort))
}
.a_co_to_y<- function(age, cohort) cohort + age - 1
.a_y_to_co<- function(age, year) year - age + 1
.co_y_to_a<- function(cohort, year) year - cohort + 1


#' Convert a stars object between using year or cohort indices.
#' 
#' @param object A univariate stars object with either a "cohort" or "year" 
#'     dimension and an "age" dimension.
#' 
#' @return A copy of the input with the year index switched to a cohort index,
#'     or the cohort index switched to a year index.
#' 
#' @export
exchange_cohort_year<- function(object) {
    if( "cohort" %in% dimnames(object) ) return(.exchange_c2y(object))
    if( "year" %in% dimnames(object) ) return(.exchange_y2c(object))
}

.exchange_y2c<- function(object) {
    olddim<- dim(object)
    n_age<- olddim["age"]
    n_year<- olddim["year"]
    n_cohort<- convert_ayc(age = n_age, year = n_year)
    newdim<- olddim
    names(newdim)[names(newdim) == "year"]<- "cohort"
    newdim["cohort"]<- n_cohort
    newobj<- array(NA, dim = newdim)

    oldidx<- olddim |>
        lapply(seq) |>
        do.call(expand.grid, args = _)
    oldidx$val<- object[[1]][as.matrix(oldidx)]
    oldidx<- oldidx[!is.na(oldidx$val), ]
    newidx<- oldidx
    colnames(newidx)[colnames(newidx) == "year"]<- "cohort"
    newidx$cohort<- convert_ayc(age = oldidx$age, year = oldidx$year)
    newobj[as.matrix(newidx[, colnames(newidx) != "val"])]<- newidx$val
    newdim<- stars::st_dimensions(object)
    names(newdim)[names(newdim) == "year"]<- "cohort"
    newdim$cohort<- stars::st_dimensions(cohort = seq(n_cohort))$cohort

    newobj<- stars::st_as_stars(
        list(x = newobj),
        dimensions = newdim
    )
    names(newobj)<- names(object)
    return(newobj)
}
.exchange_c2y<- function(object) {
    olddim<- dim(object)
    n_age<- olddim["age"]
    n_cohort<- olddim["cohort"]
    n_year<- convert_ayc(age = n_age, cohort = n_cohort)
    newdim<- olddim
    names(newdim)[names(newdim) == "cohort"]<- "year"
    newdim["year"]<- n_year
    newobj<- array(NA, dim = newdim)

    oldidx<- olddim |>
        lapply(seq) |>
        do.call(expand.grid, args = _)
    oldidx$val<- object[[1]][as.matrix(oldidx)]
    newidx<- oldidx
    colnames(newidx)[colnames(newidx) == "cohort"]<- "year"
    newidx$year<- convert_ayc(age = oldidx$age, cohort = oldidx$cohort)
    newobj[as.matrix(newidx[, colnames(newidx) != "val"])]<- newidx$val
    newdim<- stars::st_dimensions(object)
    names(newdim)[names(newdim) == "cohort"]<- "year"
    newdim$year<- stars::st_dimensions(year = seq(n_year))$year

    newobj<- stars::st_as_stars(
        list(x = newobj),
        dimensions = newdim
    )
    names(newobj)<- names(object)
    return(newobj)
}