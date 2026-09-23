#' Simulate a scientific survey
#' 
#' @param target
#'     A n_class x n_age x n_year x n_geometry stars object giving the target
#'     surface to survey.
#' @param selectivity
#'     A numeric vector giving the survey selectivity by class.
#' @param stations_per_year
#'     The number of survey stations per year.
#' @param catchability_mean
#'     The mean catchability of the survey.
#' @param catchability_power
#'     The power parameter in a Tweedie distribution for survey catch. Must be
#'     between 1 and 2.
#' @param catchability_dispersion
#'     The dispersion parameter in a Tweedie distribution for survey catch. Must
#'     be positive.
#' @param composition_concentration
#'     The concentration parameter in a Dirichlet distribution for the survey
#'     catch composition. Must be positive.
#' 
#' @return 
#'     A list with elements
#'     * stations A data.frame giving the year and geometry for each survey 
#'           station.
#'     * survey A n_class x n_age x n_survey stars array giving the survey 
#'           count for each survey station.
#'     * survey_by_geometry A n_class x n_age x n_year x n_geometry stars array
#'           giving the combined survey count for each year.
#' 
#' @export
sim_survey<- \(
    target,
    selectivity,
    stations_per_year,
    catchability_mean,
    catchability_power,
    catchability_dispersion,
    composition_concentration
) {
    n_class<- dim(target)["class"]
    n_age<- dim(target)["age"]
    n_year<- dim(target)["year"]
    n_geom<- dim(target)["geometry"]
    
    selected_target<- target |> sweep(1, selectivity, `*`)
    combined_target<- selected_target |> stars::st_apply(3:4, sum, na.rm = TRUE)
    
    stations<- target |>
        sf::st_geometry() |>
        sf::st_area() |>
        (\(x) sample(
            length(x), 
            size = stations_per_year * n_year,
            replace = TRUE,
            prob = x
        ))() |>
        data.frame(
            year = rep(seq_len(n_year), each = stations_per_year),
            geometry = _
        )
    stations$n<- mgcv::rTweedie(
            catchability_mean * combined_target$sum[as.matrix(stations)],
            p = catchability_power,
            phi = catchability_dispersion
        ) |> 
        round()
    survey<- array(NA, dim = c(n_class, n_age, nrow(stations)))
    survey_by_geom<- array(0, dim = c(n_class, n_age, n_year, n_geom))
    for( i in seq_len(nrow(stations)) ) {
        stat<- stations[i, ]
        iseltarget<- selected_target[[1]][, , stat$year, stat$geometry]
        iseltarget[is.na(iseltarget)]<- 0
        comp<- elhelpers::rdirichlet(
            1,
            pi = iseltarget |> (\(x) x / sum(x))(),
            concentration = composition_concentration
        )
        survey[, , i]<- rmultinom(
            n = 1,
            size = stat$n,
            prob = comp
        )
        survey_by_geom[, , stat$year, stat$geometry]<-
            survey_by_geom[, , stat$year, stat$geometry] +
            survey[, , i]
    }


    survey<- stars::st_as_stars(
        list(count = survey),
        dimensions = stars::st_dimensions(
            class = seq_len(n_class),
            age = seq_len(n_age),
            station = seq_len(nrow(stations))
        )
    )
    survey_by_geom<- stars::st_as_stars(
        list(count = survey_by_geom),
        dimensions = stars::st_dimensions(
            class = seq_len(n_class),
            age = seq_len(n_age),
            year = seq_len(n_year),
            geometry = target |> sf::st_geometry()
        )
    )
    return(
        list(
            stations = stations[, 1:2],
            survey = survey,
            survey_by_geometry = survey_by_geom
        )
    )
}