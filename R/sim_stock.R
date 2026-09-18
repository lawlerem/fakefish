#' Simulate recruitment for a fish stock
#' 
#' @param n_class
#'     The number of size classes in the stock.
#' @param n_cohort
#'     The number of cohorts tracked in the stock.
#' @param geometry
#'     An sf object with polygon geometries describing the geometry of the stock
#'     area.
#' @param mean_recruitment_composition
#'     The typical composition of the recruit class. Entries will be rescaled to
#'      sum to  1.
#' @param mean_recruit_abundance
#'     The typical total abudance of the recruit class.
#' @param p_composition_cohort
#'     Parameters (height, stretch) governing yearly variability in composition.
#' @param p_composition_geometry
#'     Parameters (height, stretch) governing spatial variability in 
#'     composition.
#' @param p_abundance_cohort 
#'     Parameters (height, stretch) governing yearly variability in abundance.
#' @param p_abundance_geometry
#'     Parameters (height, stretch) governing spatial variability in abundance.
#' 
#' @return
#'     A list with elements
#'     * abundance A n_class x n_cohort x n_geometry stars array giving the 
#'           recruit abundance.
#'     * recruitment_splines A list with splines governing variability.
#'     
#' 
#' @export
sim_recruitment<- \(
    n_class, 
    n_cohort,
    geometry,
    mean_recruitment_composition = c(1, rep(0.01 / (n_class - 1), n_class - 1)),
    mean_recruit_abundance,
    p_composition_cohort = c(0.3, 6),
    p_composition_geometry = c(0.3, 6),
    p_abundance_cohort = c(0.3, 6),
    p_abundance_geometry = c(0.3, 6)
) {
    delta0 = mean_recruitment_composition |>
        (\(x) x / sum(x))() |>
        elhelpers::icrl() |>
        c()
    cohort_spline<- seq_len(n_cohort) |> nnspline::create_nnspline()
    geom_spline<- geometry |> sf::st_centroid() |> sf::st_coordinates() |>
        nnspline::create_nnspline()
    composition_cohort_spline<- cohort_spline |>
        nnspline::update_spline(parameters = p_composition_cohort) |>
        nnspline::rspline(center = TRUE)
    composition_geometry_spline<- geom_spline |>
        nnspline::update_spline(parameters = p_composition_geometry) |>
        nnspline::rspline(center = TRUE)
    abundance_cohort_spline<- cohort_spline |>
        nnspline::update_spline(parameters = p_abundance_cohort) |>
        nnspline::rspline(center = TRUE)
    abundance_geometry_spline<- geom_spline |>
        nnspline::update_spline(parameters = p_abundance_geometry) |>
        nnspline::rspline(center = TRUE)

    delta<- delta0 |>
        outer(composition_cohort_spline$values, `+`) |>
        outer(composition_geometry_spline$values, `+`)
    composition<- delta |> 
        aperm(c(2:3, 1)) |> 
        elhelpers::crl() |> 
        aperm(c(3, 1:2))

    abundance<- outer(
            abundance_cohort_spline$values, 
            abundance_geometry_spline$values,
            `+`
        ) |>
        exp() |>
        (\(x) mean_recruit_abundance * x)()
    recruits<- sweep(
        composition,
        2:3,
        abundance,
        `*`
    )

    return(
        list(
            recruits = stars::st_as_stars(
                list(recruits = recruits),
                dimensions = stars::st_dimensions(
                    class = seq_len(n_class),
                    cohort = seq_len(n_cohort),
                    geometry = geometry |> sf::st_geometry()
                )
            ),
            recruitment_splines = list(
                composition_cohort = composition_cohort_spline,
                composition_geometry = composition_geometry_spline,
                abundance_cohort = abundance_cohort_spline,
                abundance_geometry = abundance_geometry_spline
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
#' @param geometry
#'     An sf object with polygon geometries describing the geometry of the stock
#'     area.
#' @param mean_growth_rate
#'     The typical growth rate from one class to the next.
#' @param p_growth_class 
#'     Parameters (height, stretch) governing class variability in growth.
#' @param p_growth_age 
#'     Parameters (height, stretch) governing age variability in growth.
#' @param p_growth_year 
#'     Parameters (height, stretch) governing yearly variability in growth.
#' @param p_growth_geometry
#'     Parameters (height, stretch) governing spatial variability in growth.
#' 
#' @return
#'     A list with elements
#'     * growth_matrix A n_class x n_class x n_age x n_year x n_geom stars array
#'           giving the yearly growth transition matrix.
#'     * class_growth_rate A n_class x n_class matrix giving the average growth
#'         rate from one class to the next.
#'     * growth_rate A n_age x n_year x n_geom stars array giving adjustments to
#'           the average growth rate.
#'     * growth_splines A list with splines governing variability.
#'     
#' 
#' @export
sim_growth<- \(
    n_class,
    n_age,
    n_cohort,
    geometry,
    mean_growth_rate,
    p_growth_class = c(0.3, 2),
    p_growth_age = c(0.3, 2),
    p_growth_year = c(0.3, 2),
    p_growth_geometry = c(0.3, 2)
) {
    n_year<- convert_ayc(age = n_age, cohort = n_cohort)

    year_spline<- seq_len(n_year) |> nnspline::create_nnspline()
    class_spline<- seq_len(n_class - 1) |> nnspline::create_nnspline()
    age_spline<- seq_len(n_age) |> nnspline::create_nnspline()
    geom_spline<- geometry |> sf::st_centroid() |> sf::st_coordinates() |>
        nnspline::create_nnspline()

    class_spline<- class_spline |>
        nnspline::update_spline(parameters = p_growth_class) |>
        nnspline::rspline(center = TRUE)
    age_spline<- age_spline |>
        nnspline::update_spline(parameters = p_growth_age) |>
        nnspline::rspline(center = TRUE)
    year_spline<- year_spline |>
        nnspline::update_spline(parameters = p_growth_year) |>
        nnspline::rspline(center = TRUE)
    geom_spline<- geom_spline |>
        nnspline::update_spline(parameters = p_growth_geometry) |>
        nnspline::rspline(center = TRUE)
    
    growth_rate<- exp(age_spline$values) |>
        outer(exp(year_spline$values)) |>
        outer(exp(geom_spline$values))

    class_growth_rate<- matrix(0, nrow = n_class, ncol = n_class)
    diag(class_growth_rate)<- c(
        -mean_growth_rate * exp(class_spline$values), 0
    )
    subdiagonal<- cbind(seq(n_class - 1) + 1, seq(n_class - 1))
    class_growth_rate[subdiagonal]<- -head(diag(class_growth_rate), -1)

    growth_matrix<- array(
        0, 
        dim = c(n_class, n_class, n_age, n_year, nrow(geometry))
    )
    for( a in seq_len(n_age) ) {
        for( y in seq_len(n_year) ) {
            for( s in seq(nrow(geometry)) ) {
                growth_matrix[, , a, y, s]<- (growth_rate[a, y, s] * class_growth_rate) |>
                    Matrix::expm() |> 
                    as.matrix()
            }
        }
    }

    return(
        list(
            growth_matrix = stars::st_as_stars(
                list(growth_matrix = growth_matrix),
                dimensions = stars::st_dimensions(
                    to_class = seq_len(n_class),
                    from_class = seq_len(n_class),
                    age = seq_len(n_age),
                    year = seq_len(n_year),
                    geometry = geometry |> sf::st_geometry()
                )
            ),
            class_growth_rate = class_growth_rate,
            growth_rate = stars::st_as_stars(
                list(growth_rate = growth_rate),
                dimensions = stars::st_dimensions(
                    age = seq_len(n_age),
                    year = seq_len(n_year),
                    geometry = geometry |> sf::st_geometry()
                )
            ),
            growth_splines = list(
                class_spline = class_spline,
                age_spline = age_spline,
                year_spline = year_spline,
                geometry_spline = geom_spline
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
#' @param geometry
#'     An sf object with polygon geometries describing the geometry of the stock
#'     area.
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
#' @param p_fishing_geometry 
#'     Parameters (height, stretch) governing spatial variability in fishing 
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
#' @param p_natural_geometry 
#'     Parameters (height, stretch) governing spatial variability in natural 
#'     mortality.
#' 
#' @return
#'     A list with elements
#'     * total_mortality A n_class x n_age x n_year x n_geom stars array giving
#'           the yearly total mortality rate.
#'     * fishing_mortality A n_class x n_year x n_geom stars array giving the
#'           yearly fishing mortality rate.
#'     * natural_mortality A n_class x n_age x n_year x n_geom stars array
#'           giving the yearly natural mortality rate.
#'     * growth_splines A list with splines governing variability.
#'     
#' 
#' @export
sim_mortality<- \(
    n_class,
    n_age,
    n_cohort,
    geometry,
    mean_fishing_mortality,
    fishing_selectivity,
    mean_natural_mortality,
    p_fishing_year = c(0.3, 6),
    p_fishing_geometry = c(0.3, 6),
    p_natural_class = c(0.3, 2),
    p_natural_age = c(0.3, 2),
    p_natural_year = c(0.3, 2),
    p_natural_geometry = c(0.3, 2)
) {
    n_year<- convert_ayc(age = n_age, cohort = n_cohort)

    year_spline<- seq_len(n_year) |> nnspline::create_nnspline()
    class_spline<- seq_len(n_class) |> nnspline::create_nnspline()
    age_spline<- seq_len(n_age) |> nnspline::create_nnspline()
    geom_spline<- geometry |> sf::st_centroid() |> sf::st_coordinates() |>
        nnspline::create_nnspline()

    fishing_year_spline<- year_spline |>
        nnspline::update_spline(parameters = p_fishing_year) |>
        nnspline::rspline(center = TRUE)
    fishing_geometry_spline<- geom_spline |>
        nnspline::update_spline(parameters = p_fishing_geometry) |>
        nnspline::rspline(center = TRUE)
    natural_class_spline<- class_spline |>
        nnspline::update_spline(parameters = p_natural_class) |>
        nnspline::rspline(center = TRUE)
    natural_age_spline<- age_spline |>
        nnspline::update_spline(parameters = p_natural_age) |>
        nnspline::rspline(center = TRUE)
    natural_year_spline<- year_spline |>
        nnspline::update_spline(parameters = p_natural_year) |>
        nnspline::rspline(center = TRUE)
    natural_geometry_spline<- geom_spline |>
        nnspline::update_spline(parameters = p_natural_geometry) |>
        nnspline::rspline(center = TRUE)
   
    # Mortality<- log(selectivity) + exp(log(mu_mort) + year + geom)
    fishing_mortality<- fishing_year_spline$values |>
        outer(
            fishing_geometry_spline$values,
            `+`
        ) |>
        exp() |>
        (\(x) x * mean_fishing_mortality)() |>
        outer(
            X = log(fishing_selectivity),
            Y = _,
            `+`
        )

    natural_mortality<- exp(natural_class_spline$values) |>
        outer(exp(natural_age_spline$values)) |>
        outer(exp(natural_year_spline$values)) |>
        outer(exp(natural_geometry_spline$values)) |>
        (\(x) mean_natural_mortality * x)()
    
    total_mortality<- sweep(
        natural_mortality, 
        c(1, 3:4), 
        fishing_mortality,
        `+`
    )

    return(
        list(
            total_mortality = stars::st_as_stars(
                list(total_mortality = total_mortality),
                dimensions = stars::st_dimensions(
                    class = seq_len(n_class),
                    age = seq_len(n_age),
                    year = seq_len(n_year),
                    geometry = geometry |> sf::st_geometry()
                )
            ),
            fishing_mortality = stars::st_as_stars(
                list(fishing_mortality = fishing_mortality),
                dimensions = stars::st_dimensions(
                    class = seq_len(n_class),
                    year = seq_len(n_year),
                    geometry = geometry |> sf::st_geometry()
                )
            ),
            natural_mortality = stars::st_as_stars(
                list(natural_mortality = natural_mortality),
                dimensions = stars::st_dimensions(
                    class = seq_len(n_class),
                    age = seq_len(n_age),
                    year = seq_len(n_year),
                    geometry = geometry |> sf::st_geometry()
                )
            ),
            mortality_splines = list(
                fishing_year = fishing_year_spline,
                fishing_geometry = fishing_geometry_spline,
                natural_class = natural_class_spline,
                natural_age = natural_age_spline,
                natural_year = natural_year_spline,
                natural_geometry = natural_geometry_spline
            )
        )
    )
}




#' Simulate mortality parameters for a fish stock
#' 
#' @param recruitment
#'     The output of sim_recruitment, or a list with a n_class x n_cohort x
#'     n_geom array named "recruits" giving recruit abundance.
#' @param growth
#'     The output of sim_growth, or a list with a n_class x n_class x n_age x 
#'     n_age x n_geom stars array named "growth_matrix" giving the yearly growth
#'     transition matrix.
#' @param mortality
#'     The output of sim_mortality, or a list with a n_class x n_age x n_year x
#'     n_geom stars array named "total_mortality" giving the total mortality 
#'     rate.
#' @param fishing_mean_time
#'     The average number of steps spent in each polygon.
#' @param fishing_mean_duration
#'     The average duration of fishing trips.
#' @param fishing_mean_trips
#'     The average number of distinct fishing trips per year.
#' @param catchability_mean
#'     The average catchability.
#' @param catchability_var
#'     The variance of qlogis(catchability)
#' @param fishing_area_mean
#'     The average area fished per step.
#' @param fishing_area_var
#'     The variance of area fished per step.
#' 
#' @return
#'     A list with elements
#'     * abundance A n_class x n_age x n_year x n_geom stars array giving the
#'           stock abundance.
#' 
#' @export
sim_abundance<- \(
    recruitment,
    growth,
    mortality,
    fishing_mean_time,
    fishing_mean_duration,
    fishing_mean_trips,
    fishing_selectivity,
    catchability_mean,
    catchability_var,
    fishing_area_mean,
    fishing_area_var
) {
    n_class<- dim(mortality$total_mortality)[1]
    n_age<- dim(mortality$total_mortality)[2]
    n_year<- dim(mortality$total_mortality)[3]
    n_cohort<- convert_ayc(age = n_age, year = n_year)
    n_geom<- dim(mortality$total_mortality)[4]
    geom<- sf::st_geometry(mortality$total_mortality)

    fishing_track<-
        fishing_catch<-
        fishing_catch_by_geom<- 
        fishing_catchability<- 
        fishing_area<- vector("list", n_year)

    abundance<- array(NA, dim = c(n_class, n_age, n_year, 2, n_geom))
    abundance[, 1, seq_len(n_cohort), 1, ]<- recruitment$recruits$recruits
    for( y in seq_len(n_year) ) {
        # fish stuff out
        target<- abundance[ , , y, 1, ] |>
            sweep(1, fishing_selectivity, `*`) |>
            apply(3, sum, na.rm = TRUE)
        fishing_track[[y]]<- (1 + rpois(1, fishing_mean_trips)) |>
            seq_len() |>
            lapply(\(i) generate_fishing_track(
                geometry = geom,
                target = target,
                mean_time = fishing_mean_time * runif(1, 0.9, 1.1),
                duration = fishing_mean_duration * runif(1, 0.9, 1.1)
            )) |>
            do.call(c, args = _)
        catch<- fish_track(
            fishing_track = fishing_track[[y]],
            fish = abundance[, , y, 1, ] |> (\(x) {x[is.na(x)]<- 0; x})(),
            selectivity = fishing_selectivity,
            geometry = geometry,
            catchability_mean = catchability_mean,
            catchability_var = catchability_var,
            area_mean = fishing_area_mean,
            area_var = fishing_area_var
        )
        fishing_catch[[y]]<- catch$catch
        fishing_catch_by_geom[[y]]<- catch$catch_by_geometry
        fishing_catchability[[y]]<- catch$catchability
        fishing_area[[y]]<- catch$area
        abundance[, , y, 2, ]<- (abundance[, , y, 1, ] - catch$catch_by_geometry$catch) |>
            sapply(max, 0)

        for( a in seq_len(n_age - 1) ) {
            co<- convert_ayc(age = a, year = y)
            if( co > n_cohort ) next
            for( s in seq_len(n_geom) ) {
                G<- growth$growth_matrix$growth_matrix[, , a, y, s]
                M<- mortality$natural_mortality$natural_mortality[, a, y, s]
                abundance[, a + 1, y + 1, 1, s]<-
                    G %*% diag(exp(-M)) %*% abundance[, a, y, 2, s]
            }
        }
    }

    return(
        list(
            abundance = stars::st_as_stars(
                list(abundance = abundance),
                dimensions = stars::st_dimensions(
                    class = seq_len(n_class),
                    age = seq_len(n_age),
                    year = seq_len(n_year),
                    season = c("prefishing", "postfishing"),
                    geometry = mortality$total_mortality |> sf::st_geometry()
                )
            ),
            fishing_track = fishing_track,
            fishing_catch = fishing_catch,
            fishing_catch_by_geom = fishing_catch_by_geom,
            fishing_catchability = fishing_catchability,
            fishing_area = fishing_area
        )
    )
}



# #' @param mean_recruitment_composition
# #'     See ?sim_recruitment.
# #' @param mean_recruit_abundance
# #'     See ?sim_recruitment.
# #' @param mean_growth_rate
# #'     See ?sim_growth.
# #' @param mean_fishing_mortality
# #'     See ?sim_mortality.
# #' @param mean_natural_mortality
# #'     See ?sim_mortality.
# #' @param fishing_mean_time
# #'     See ?sim_abundance
# #' @param fishing_mean_duration
# #'     See ?sim_abundance
# #' @param fishing_mean_tracks
# #'     See ?sim_abundance
# #' @param fishing_selectivity
# #'     See ?sim_abundance
# #' @param catchability_mean
# #'     See ?sim_abundance
# #' @param catchability_var
# #'     See ?sim_abundance
# #' @param fishing_area_mean
# #'     See ?sim_abundance
# #' @param fishing_area_var
# #'     See ?sim_abundance
# #' @param p_composition_cohort
# #'     See ?sim_recruitment
# #' @param p_composition_geometry
# #'     See ?sim_recruitment
# #' @param p_abundance_cohort
# #'     See ?sim_recruitment
# #' @param p_abundance_geometry
# #'     See ?sim_recruitment
# #' @param p_growth_class
# #'     See ?sim_growth
# #' @param p_growth_age
# #'     See ?sim_growth
# #' @param p_growth_year
# #'     See ?sim_growth
# #' @param p_growth_geometry
# #'     See ?sim_growth
# #' @param p_fishing_year
# #'     See ?sim_mortality
# #' @param p_fishing_geometry
# #'     See ?sim_mortality
# #' @param p_natural_class
# #'     See ?sim_mortality
# #' @param p_natural_age
# #'     See ?sim_mortality
# #' @param p_natural_year
# #'     See ?sim_mortality
# #' @param p_natural_geometry
# #'     See ?sim_mortality



#' Simulate a fish stock
#' 
#' @param n_class
#'     The number of size classes in the stock.
#' @param n_age
#'     The number of ages tracked in the stock.
#' @param n_cohort
#'     The number of cohorts tracked in the stock.
#' @param geometry
#'     An sf object with polygon geometries describing the fishing area.
#' @inheritParams sim_recruitment
#' @inheritParams sim_growth
#' @inheritParams sim_mortality
#' @inheritParams sim_abundance
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
    geometry,
    mean_recruitment_composition,
    mean_recruit_abundance,
    mean_growth_rate,
    mean_fishing_mortality,
    mean_natural_mortality,
    fishing_mean_time,
    fishing_mean_duration,
    fishing_mean_trips,
    fishing_selectivity,
    catchability_mean,
    catchability_var,
    fishing_area_mean,
    fishing_area_var,
    p_composition_cohort = c(0.3, 6),
    p_composition_geometry = c(0.3, 6),
    p_abundance_cohort = c(0.3, 6),
    p_abundance_geometry = c(0.3, 6),
    p_growth_class = c(0.3, 2),
    p_growth_age = c(0.3, 2),
    p_growth_year = c(0.3, 2),
    p_growth_geometry = c(0.3, 2),
    p_fishing_year = c(0.3, 6),
    p_fishing_geometry = c(0.3, 6),
    p_natural_class = c(0.3, 2),
    p_natural_age = c(0.3, 2),
    p_natural_year = c(0.3, 2),
    p_natural_geometry = c(0.3, 2)
) {
    recruitment<- sim_recruitment(
        n_class = n_class,
        n_cohort = n_cohort,
        geometry = geometry,
        mean_recruitment_composition = mean_recruitment_composition,
        mean_recruit_abundance = mean_recruit_abundance,
        p_composition_cohort = p_composition_cohort,
        p_composition_geometry = p_composition_geometry,
        p_abundance_cohort = p_abundance_cohort,
        p_abundance_geometry = p_abundance_geometry
    )
    growth<- sim_growth(
        n_class = n_class,
        n_age = n_age,
        n_cohort = n_cohort,
        geometry = geometry,
        mean_growth_rate = mean_growth_rate,
        p_growth_class = p_growth_class,
        p_growth_age = p_growth_age,
        p_growth_year = p_growth_year,
        p_growth_geometry = p_growth_geometry
    )
    mortality<- sim_mortality(
        n_class = n_class,
        n_age = n_age,
        n_cohort = n_cohort,
        geometry = geometry,
        mean_fishing_mortality = mean_fishing_mortality,
        fishing_selectivity = fishing_selectivity,
        mean_natural_mortality = mean_natural_mortality,
        p_fishing_year = p_fishing_year,
        p_fishing_geometry = p_fishing_geometry,
        p_natural_class = p_natural_class,
        p_natural_age = p_natural_age,
        p_natural_year = p_natural_year,
        p_natural_geometry = p_natural_geometry
    )
    abundance<- sim_abundance(
        recruitment = recruitment,
        growth = growth,
        mortality = mortality,
        fishing_mean_time = fishing_mean_time,
        fishing_mean_duration = fishing_mean_duration,
        fishing_mean_trips = fishing_mean_trips,
        fishing_selectivity = fishing_selectivity,
        catchability_mean = catchability_mean,
        catchability_var = catchability_var,
        fishing_area_mean = fishing_area_mean,
        fishing_area_var = fishing_area_var 
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
