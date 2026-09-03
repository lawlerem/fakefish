#' Plot growth parameters
#' 
#' @param growth
#'     The output of sim_growth, or a list with elements:
#'     * growth_matrix A n_class x n_class x n_age x n_year array.
#'     * growth_rate A n_class x n_age x n_year array.
#' @param age, year
#'     Either missing or an integer index. At least one must be missing.
#'     If both are missing, the growth rate adjustment is plotted.
#'     If age is missing and year is supplied, the growth transition matrix for
#'         the given year is plotted as a function of age.
#'     If age is supplied and year is missing, the growth transition matrix for
#'         the given age is plotted as a function of year.
#' 
#' @export
plot_growth<- function(
    growth,
    age,
    year
) { 
    n_class<- dim(growth$growth_matrix)[1]
    n_age<- dim(growth$growth_rate)[1]
    n_year<- dim(growth$growth_rate)[2]

    if( !missing(age) + !missing(year) == 2 ) {
        stop("At least one of the age or year arguments must be missing.")
    }

    if( missing(age) & missing(year) ) {
        opar<- par(
            mfrow = c(
                ceiling(sqrt(n_age)), 
                ceiling(n_age / ceiling(sqrt(n_age)))
            ),
            mar = c(1.1, 1.1, 1.1, 1.1)
        )
        ylim<- range(growth$growth_rate)
        for( a in seq_len(n_age) ) {
            plot(
                x = seq_len(n_year),
                y = growth$growth_rate[a, ],
                ylim = ylim,
                type = "l"
            )
        }
        par(opar)
    }

    if( !missing(age) ) {
        opar<- par(
            mfrow = c(n_class, n_class),
            mar = c(1.1, 1.1, 1.1, 1.1)
        )
        for( to in seq_len(n_class) ) {
        for( from in seq_len(n_class) ) {
            if( to < from ) {plot.new(); next}
            plot(
                x = seq_len(n_year),
                y = growth$growth_matrix[to, from, age, ],
                ylim = c(0, 1),
                type = "l"
            )
        }
        }
        par(opar)
    }

    if( !missing(year) ) {
        opar<- par(
            mfrow = c(n_class, n_class),
            mar = c(1.1, 1.1, 1.1, 1.1)
        )
        for( to in seq_len(n_class) ) {
        for( from in seq_len(n_class) ) {
            if( to < from ) {plot.new(); next}
            plot(
                x = seq_len(n_age),
                y = growth$growth_matrix[to, from, , year],
                ylim = c(0, 1),
                type = "l"
            )
        }
        }
        par(opar)
    }

    return(invisible())
}




#' Plot mortality parameters
#' 
#' @param mortality
#'     The output of sim_mortality, or a list with elements:
#'     * total_mortality A n_class x n_age x n_year array.
#'     * fishing_mortality A n_class x n_year array.
#'     * natural_mortality A n_class x n_age x n_year array.
#' @param type
#'     Either "total", "fishing" or "natural", determining which mortality type
#'     should be plotted.
#' 
#' @export
plot_mortality<- function(
    mortality,
    type = "total"
) {
    n_class<- dim(mortality$total_mortality)[1]
    n_age<- dim(mortality$total_mortality)[2]
    n_year<- dim(mortality$total_mortality)[3]

    mort<- switch(
        type,
        total = mortality$total_mortality,
        fishing = mortality$fishing_mortality,
        natural = mortality$natural_mortality
    )
    ylim<- range(mort)

    if( type == "fishing" ) {
        opar<- par(
            mfrow = c(n_class, 1),
            mar = c(1.1, 3.1, 1.1, 1.1)
        )
        for( c in rev(seq_len(n_class)) ) {
            plot(
                x = seq_len(n_year),
                y = mort[c, ],
                ylim = ylim,
                type = "l"
            )
        }
        par(opar)
        return(invisible())
    }

    opar<- par(
        mfrow = c(n_class, n_age),
        mar = c(1.1, 1.1, 1.1, 1.1)
    )
    for( c in rev(seq_len(n_class)) ) {
    for( a in seq_len(n_age) ) {
        plot(
            x = seq_len(n_year),
            y = mort[c, a, ],
            ylim = ylim,
            type = "l"
        )
    }
    }
    par(opar)

    return(invisible())
}



#' Plot abundance
#' 
#' @param abundance
#'     The output of sim_abundance, or a list with elements:
#'     * abundance A n_class x n_age x n_year array.
#' @param by
#'     Either "year" or "cohort" determining whether the stock should be tracked
#'     on a year-to-year basis or a cohort-to-cohort basis.
#' 
#' @export
plot_abundance<- function(
    abundance,
    by = "year"
) {
    n_class<- dim(abundance$abundance)[1]
    n_age<- dim(abundance$abundance)[2]
    n_year<- dim(abundance$abundance)[3]
    n_cohort<- convert_ayc(age = n_age, year = n_year)

    cuml_abundance<- abundance$abundance |>
        apply(2:3, cumsum)

    ylim<- range(cuml_abundance, na.rm = TRUE)

    if( by == "year" ) {
        opar<- par(
            mfrow = c(
                ceiling(sqrt(n_age)), 
                ceiling(n_age / ceiling(sqrt(n_age)))
            ),
            mar = c(1.1, 1.1, 1.1, 1.1)
        )
        for( a in seq_len(n_age) ) {
            plot(
                x = c(1, n_year),
                y = ylim,
                ylim = ylim,
                type = "n"
            )
            for( c in seq_len(n_class) ) {
                lines(
                    x = seq_len(n_year),
                    y = cuml_abundance[c, a, ]
                )
            }
        }
        par(opar)
    }
    if( by == "cohort" ) {
        opar<- par(
            mfrow = c(
                ceiling(sqrt(n_cohort)), 
                ceiling(n_cohort / ceiling(sqrt(n_cohort)))
            ),
            mar = c(1.1, 1.1, 1.1, 1.1)
        )
        for( co in seq_len(n_cohort) ) {
            plot(
                x = c(1, n_age),
                y = ylim,
                type = "n"
            )
            for( c in seq_len(n_class) ) {
                idx<- cbind(
                    c,
                    seq_len(n_age),
                    convert_ayc(age = seq_len(n_age), cohort = co)
                )
                lines(
                    x = seq_len(n_age),
                    y = cuml_abundance[idx],
                )
            }
        }
        par(opar)
    }

    return(invisible())
}