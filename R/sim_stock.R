#' Simulate recruitment for a fish stock
#' 
#' @param n_class
#'     The number of size classes in the stock.
#' @param n_cohort
#'     The number of cohorts tracked in the stock.
#' @param mean_recruitment_composition
#'     The typical composition of the recruit class. Entries will be rescaled to
#'      sum to  1.
#' @param mean_recruit_abundance
#'     The typical total abudance of the recruit class.
#' @param p_composition_cohort
#'     Parameters (height, stretch) governing yearly variability in composition.
#' @param p_abundance_cohort 
#'     Parameters (height, stretch) governing yearly variability in abundance.
#' 
#' @return
#'     A list with elements
#'     * abundance A n_class x n_cohort array giving the recruit abundance.
#'     * abundance_cohort A list with splines governing variability.
#'     
#' 
#' @export
sim_recruitment<- \(
    n_class, 
    n_cohort,
    mean_recruitment_composition = c(1, rep(0.01 / (n_class - 1), n_class - 1)),
    mean_recruit_abundance,
    p_composition_cohort = c(0.3, 6),
    p_abundance_cohort = c(0.3, 6)
) {
    delta0 = mean_recruitment_composition |>
        (\(x) x / sum(x))() |>
        elhelpers::icrl() |>
        c()
    composition_cohort_spline<- nnspline::create_nnspline(
            x = seq_len(n_cohort),
            parameters = p_composition_cohort
        ) |>
        nnspline::rspline(center = TRUE)
    abundance_cohort_spline<- nnspline::create_nnspline(
            x = seq_len(n_cohort),
            parameters = p_abundance_cohort
        ) |>
        nnspline::rspline(center = TRUE)

    delta<- outer(delta0, composition_cohort_spline$values, `+`)
    composition<- delta |> t() |> elhelpers::crl() |> t()

    recruits<- sweep(
        composition,
        2,
        mean_recruit_abundance * exp(abundance_cohort_spline$values),
        `*`
    )

    return(
        list(
            recruits = recruits,
            recruitment_splines = list(
                composition_cohort = composition_cohort_spline,
                abundance_cohort = abundance_cohort_spline
            )
        )
    )
}






#' Simulate mortality parameters for a fish stock
#' 
#' @param n_class
#'     The number of size classes in the stock.
#' @param n_age
#'     The number of ages tracked in the stock.
#' @param n_cohort
#'     The number of cohorts tracked in the stock.
#' @param mean_fishing_mortality
#'     The typical fishing mortality rate.
#' @param fishing_selectivity
#'     A n_class vector giving the gear selectivity of the fishery. Elements
#'         should be be between 0 and 1.
#' @param mean_natural_mortality
#'     The typical natural mortality rate.
#' @param p_fishing_year 
#'     Parameters (height, stretch) governing yearly variability in fishing 
#'     mortality.
#' @param p_natural_class 
#'     Parameters (height, stretch) governing class variability in natural 
#'     mortality.
#' @param p_natural_age 
#'     Parameters (height, stretch) governing age variability in natural 
#'     mortality.
#' @param p_natural_year 
#'     Parameters (height, stretch) governing yearly variability in natural 
#'     mortality.
#' 
#' @return
#'     A list with elements
#'     * total_mortality A n_class x n_age x n_year array giving the yearly 
#'         total mortality rate.
#'     * fishing_mortality A n_class x n_year array giving the yearly fishing
#'         mortality rate.
#'     * natural_mortality A n_class x n_age x n_year array giving the yearly 
#'         natural mortality rate.
#'     * growth_splines A list with splines governing variability.
#'     
#' 
#' @export
sim_mortality<- \(
    n_class,
    n_age,
    n_cohort,
    mean_fishing_mortality,
    fishing_selectivity,
    mean_natural_mortality,
    p_fishing_year = c(0.3, 6),
    p_natural_class = c(0.3, 2),
    p_natural_age = c(0.3, 2),
    p_natural_year = c(0.3, 2)
) {
    n_year<- convert_ayc(age = n_age, cohort = n_cohort)
    fishing_year_spline<- nnspline::create_nnspline(
            x = seq_len(n_year),
            parameters = p_fishing_year
        ) |>
        nnspline::rspline(center = TRUE)
    natural_class_spline<- nnspline::create_nnspline(
            x = seq_len(n_class),
            parameters = p_natural_class
        ) |>
        nnspline::rspline(center = TRUE)
    natural_age_spline<- nnspline::create_nnspline(
            x = seq_len(n_age),
            parameters = p_natural_age
        ) |>
        nnspline::rspline(center = TRUE)
    natural_year_spline<- nnspline::create_nnspline(
            x = seq_len(n_year),
            parameters = p_natural_year
        ) |>
        nnspline::rspline(center = TRUE)
   
    fishing_mortality<- outer(
        log(fishing_selectivity),
        mean_fishing_mortality * exp(fishing_year_spline$values),
        `+`
    )

    natural_mortality<- exp(natural_class_spline$values) |>
        outer(exp(natural_age_spline$values)) |>
        outer(exp(natural_year_spline$values)) |>
        (\(x) mean_natural_mortality * x)()
    
    total_mortality<- sweep(
        natural_mortality, 
        c(1, 3), 
        fishing_mortality,
        `+`
    )

    return(
        list(
            total_mortality = total_mortality,
            fishing_mortality = fishing_mortality,
            natural_mortality = natural_mortality,
            mortality_splines = list(
                fishing_year = fishing_year_spline,
                natural_class = natural_class_spline,
                natural_age = natural_age_spline,
                natural_year = natural_year_spline
            )
        )
    )
}






#' Simulate growth parameters for a fish stock
#' 
#' @param n_class
#'     The number of size classes in the stock.
#' @param n_age
#'     The number of ages tracked in the stock.
#' @param n_cohort
#'     The number of cohorts tracked in the stock.
#' @param mean_growth_rate
#'     The typical growth rate from one class to the next.
#' @param p_growth_class 
#'     Parameters (height, stretch) governing class variability in growth.
#' @param p_growth_age 
#'     Parameters (height, stretch) governing age variability in growth.
#' @param p_growth_year 
#'     Parameters (height, stretch) governing yearly variability in growth.
#' 
#' @return
#'     A list with elements
#'     * growth_matrix A n_class x n_class x n_age x n_year array giving the 
#'         yearly growth transition matrix.
#'     * class_growth_rate A n_class x n_class matrix giving the average growth
#'         rate from one class to the next.
#'     * growth_rate A n_age x n_year matrix giving adjustments to the average
#'         growth rate.
#'     * growth_splines A list with splines governing variability.
#'     
#' 
#' @export
sim_growth<- \(
    n_class,
    n_age,
    n_cohort,
    mean_growth_rate,
    p_growth_class = c(0.3, 2),
    p_growth_age = c(0.3, 2),
    p_growth_year = c(0.3, 2)
) {
    n_year<- convert_ayc(age = n_age, cohort = n_cohort)
    class_spline<- nnspline::create_nnspline(
            x = seq_len(n_class - 1),
            parameters = p_growth_class
        ) |>
        nnspline::rspline(center = TRUE)
    age_spline<- nnspline::create_nnspline(
            x = seq_len(n_age),
            parameters = p_growth_age
        ) |>
        nnspline::rspline(center = TRUE)
    year_spline<- nnspline::create_nnspline(
            x = seq_len(n_year),
            parameters = p_growth_year
        ) |>
        nnspline::rspline(center = TRUE)
    
    growth_rate<- outer(
        exp(age_spline$values),
        exp(year_spline$values)
    )

    class_growth_rate<- matrix(0, nrow = n_class, ncol = n_class)
    diag(class_growth_rate)<- c(
        -mean_growth_rate * exp(class_spline$values), 0
    )
    subdiagonal<- cbind(seq(n_class - 1) + 1, seq(n_class - 1))
    class_growth_rate[subdiagonal]<- -head(diag(class_growth_rate), -1)

    growth_matrix<- array(0, dim = c(n_class, n_class, n_age, n_year))
    for( a in seq_len(n_age) ) {
        for( y in seq_len(n_year) ) {
            growth_matrix[, , a, y]<- (growth_rate[a, y] * class_growth_rate) |>
                Matrix::expm() |> 
                as.matrix()
        }
    }

    return(
        list(
            growth_matrix = growth_matrix,
            class_growth_rate = class_growth_rate,
            growth_rate = growth_rate,
            growth_splines = list(
                class_spline = class_spline,
                age_spline = age_spline,
                year_spline = year_spline
            )
        )
    )
}




#' Simulate mortality parameters for a fish stock
#' 
#' @param recruitment
#'     The output of sim_recruitment, or a list with a n_class x n_cohort array 
#'     named "recruits" giving recruit abundance.
#' @param growth
#'     The output of sim_growth, or a list with a n_class x n_class x n_age x 
#'     n_age array named "growth_matrix" giving the yearly growth transition 
#'     matrix.
#' @param mortality
#'     The output of sim_mortality, or a list with a n_class x n_age x n_year 
#'     array named "total_mortality" giving the total mortality rate.
#' 
#' @return
#'     A list with elements
#'     * abundance A n_class x n_age x n_year array giving the stock abundance.
#' 
#' @export
sim_abundance<- \(
    recruitment,
    growth,
    mortality
) {
    n_class<- dim(mortality$total_mortality)[1]
    n_age<- dim(mortality$total_mortality)[2]
    n_year<- dim(mortality$total_mortality)[3]
    n_cohort<- convert_ayc(age = n_age, year = n_year)

    abundance<- array(NA, dim = c(n_class, n_age, n_year))
    abundance[, 1, seq_len(n_cohort)]<- recruitment$recruits
    for( co in seq_len(n_cohort) ) {
        for( a in seq_len(n_age - 1) ) {
            y<- convert_ayc(age = a, cohort = co)
            G<- growth$growth_matrix[, , a, y]
            Z<- mortality$total_mortality[, a, y]
            abundance[, a + 1, y + 1]<- G %*% diag(exp(-Z)) %*% abundance[, a, y]
        }
    }
    return(list(abundance = abundance))
}




#' Simulate a fish stock
#' 
#' @param n_class
#'     The number of size classes in the stock.
#' @param n_age
#'     The number of ages tracked in the stock.
#' @param n_cohort
#'     The number of cohorts tracked in the stock.
#' @param mean_recruitment_composition
#'     See ?sim_recruitment.
#' @param mean_recruit_abundance
#'     See ?sim_recruitment.
#' @param mean_growth_rate
#'     See ?sim_growth.
#' @param mean_fishing_mortality
#'     See ?sim_mortality.
#' @param fishing_selectivity
#'     See ?sim_mortality.
#' @param mean_natural_mortality
#'     See ?sim_mortality.
#' @param p_composition_cohort
#'     See ?sim_recruitment
#' @param p_abundance_cohort
#'     See ?sim_recruitment
#' @param p_growth_class
#'     See ?sim_growth
#' @param p_growth_age
#'     See ?sim_growth
#' @param p_growth_year
#'     See ?sim_growth
#' @param p_fishing_year
#'     See ?sim_mortality
#' @param p_natural_class
#'     See ?sim_mortality
#' @param p_natural_age
#'     See ?sim_mortality
#' @param p_natural_year
#'     See ?sim_mortality
#' 
#' @return
#'     A list concatenating the results of sim_recruitment, sim_growth,
#'     sim_mortality, and sim_abundance.
#' 
#' @export
sim_population<- function(
    n_class,
    n_age,
    n_cohort,
    mean_recruitment_composition,
    mean_recruit_abundance,
    mean_growth_rate,
    mean_fishing_mortality,
    fishing_selectivity,
    mean_natural_mortality,
    p_composition_cohort = c(0.3, 6),
    p_abundance_cohort = c(0.3, 6),
    p_growth_class = c(0.3, 2),
    p_growth_age = c(0.3, 2),
    p_growth_year = c(0.3, 2),
    p_fishing_year = c(0.3, 6),
    p_natural_class = c(0.3, 2),
    p_natural_age = c(0.3, 2),
    p_natural_year = c(0.3, 2)
) {
    recruitment<- sim_recruitment(
        n_class = n_class,
        n_cohort = n_cohort,
        mean_recruitment_composition = mean_recruitment_composition,
        mean_recruit_abundance = mean_recruit_abundance,
        p_composition_cohort = p_composition_cohort,
        p_abundance_cohort = p_abundance_cohort
    )
    growth<- sim_growth(
        n_class = n_class,
        n_age = n_age,
        n_cohort = n_cohort,
        mean_growth_rate = mean_growth_rate,
        p_growth_class = p_growth_class,
        p_growth_age = p_growth_age,
        p_growth_year = p_growth_year
    )
    mortality<- sim_mortality(
        n_class = n_class,
        n_age = n_age,
        n_cohort = n_cohort,
        mean_fishing_mortality = mean_fishing_mortality,
        fishing_selectivity = fishing_selectivity,
        mean_natural_mortality = mean_natural_mortality,
        p_fishing_year = p_fishing_year,
        p_natural_class = p_natural_class,
        p_natural_age = p_natural_age,
        p_natural_year = p_natural_year
    )
    abundance<- sim_abundance(
        recruitment = recruitment,
        growth = growth,
        mortality = mortality
    )
    return(
        c(
            recruitment,
            growth,
            mortality,
            abundance
        )
    )
}