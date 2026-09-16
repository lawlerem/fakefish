#' Simulate the path of a fishing vessel
#' 
#' A fishing path if simulated using a Metropolis-Hastings sampler using a 
#' target surface as the probability distribution from which to sample. 
#' Transitions between polygons of the fishing area are determined by the
#' Metropolis-Hastings sampler, while the lengths of stays within each polygon
#' are sampled from independent Poisson draws plus 1.
#' 
#' @param geometry
#'     An sf object with polygon geometries describing the fishing area.
#' @param target
#'     A strictly positive numeric vector giving the target surface to fish.
#' @param mean_time
#'     A non-negative number giving the average length of stay within each
#'     polygon.
#' @param duration
#'     The total length of the fishing track.
#' 
#' @return
#'     A numeric vector giving the geometry fished. The long-run average of time
#'     fished in each geometry is proportional to the target surface.
#' 
#' @export
generate_fishing_track<- \(
    geometry,
    target,
    mean_time = 1,
    duration = 10000
) {
    transition<- sf::st_intersects(geometry, remove_self = FALSE) |>
        lapply(\(x) {
            p<- x |> length() |> (\(n) rep(1/n, n))()
            names(p)<- x
            p
        })
    pmf<- target / sum(target)
    state_lengths<- rpois(n = duration, lambda = mean_time) + 1
    state_lengths<- state_lengths |> 
        head(
            which(cumsum(state_lengths) >= duration)[1]
        )


    state_values<- numeric(length(state_lengths))
    state_values[1]<- sample(length(pmf), size = 1, prob = pmf)
    for( t in seq_along(state_values) |> head(-1) ) {
        current<- state_values[t]
        prop<- transition[[current]] |> 
            (\(p) sample(names(p), 1, prob = p))() |>
            as.numeric()
        prob<- (transition[[prop]][paste(current)] * pmf[prop]) / 
            (transition[[current]][paste(prop)] * pmf[current])
        prob<- min(prob, 1)
        state_values[t + 1]<- sample(
                c(prop, current), 
                size = 1, 
                prob = c(prob, 1 - prob)
            )
    }
    return(inverse.rle(list(lengths = state_lengths, values = state_values)))
}

#' Simulate catch along a fishing track
#' 
#' @param fishing_track
#'     An integer vector giving the polygons fished.
#' @param fish
#'     A vector giving the amount of fish in each polygon.
#' @param geometry
#'     An sf object with polygon geometries describing the fishing area.
#' @param catchability_mean
#'     The average catchability
#' @param catchability_var
#'     The variance of qlogis(catchability).
#' @param area_mean
#'     The average area fished per step of the fishing track.
#' @param area_var
#'     The variance of area fished per step of the fishing track.
#' 
#' @return
#'     A list with vectors giving the catch, catchability, and area fished for 
#'     each step of the fishing track.
#' 
#' @export
fish_track<- \(
    fishing_track,
    fish,
    geometry,
    catchability_mean,
    catchability_var,
    area_mean,
    area_var
) {
    catchability<- rnorm(
            length(fishing_track),
            qlogis(catchability_mean),
            catchability_var
        ) |>
        plogis()
    catch<- numeric(length(fishing_track))
    area<- elhelpers::rgammamv(
        length(fishing_track),
        mean = area_mean,
        var = area_var
    )
    geo_area<- sf::st_area(geometry) |> units::drop_units()
    area<- (area / geo_area[fishing_track]) |> sapply(min, 1)
    for( i in seq_along(catch) ) {
        catch[i]<- area[i] * catchability[i] * fish[fishing_track[i]]
        fish[fishing_track[i]]<- fish[fishing_track[i]] - catch[i]
    }

    return(
        list(
            catch = catch,
            catchability = catchability,
            area = area
        )
    )
}