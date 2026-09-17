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
#'     A n_class x n_age x n_geometry stars array giving the abundance of fish
#'         at each polygon.
#' @param selectivity
#'     A vector of length n_class giving the gear selectivity by class.
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
#'     A list with elements
#'     * catch A n_class x n_age x n_track stars array giving the catch for each
#'           fishing step.
#'     * catch_by_geometry A n_class x n_age x n_geom stars array giving the
#'           catch for each polygon.
#'     * catchability A vector of actual catchabilities for each fishing step.
#'     * area_fished A vector of area towed for each fishing step.
#' 
#' @export
fish_track<- \(
    fishing_track,
    fish,
    selectivity,
    geometry,
    catchability_mean,
    catchability_var,
    area_mean,
    area_var
) {
    n_class<- dim(fish)[1]
    n_age<- dim(fish)[2]
    n_geom<- dim(fish)[3]

    catchability<- rnorm(
            length(fishing_track),
            qlogis(catchability_mean),
            catchability_var
        ) |>
        plogis()
    catch<- array(0, dim = c(n_class, n_age, length(fishing_track)))
    catch_by_geo<- array(0, dim = c(n_class, n_age, nrow(geometry)))
    area<- elhelpers::rgammamv(
        length(fishing_track),
        mean = area_mean,
        var = area_var
    )
    geo_area<- sf::st_area(geometry) |> units::drop_units()
    parea<- (area / geo_area[fishing_track]) |> sapply(min, 1)


    for( i in seq_along(fishing_track) ) {
        catch[, , i]<- (parea[i] * catchability[i] * fish[[1]][, , fishing_track[i]]) |>
            sweep(1, selectivity, `*`)
        catch_by_geo[, , fishing_track[i]]<- catch_by_geo[, , fishing_track[i]] +
            catch[, , i]
        fish[[1]][, , fishing_track[i]]<- fish[[1]][, , fishing_track[i]] - catch[, , i]
    }

    return(
        list(
            catch = stars::st_as_stars(
                list(catch = catch),
                dimensions = stars::st_dimensions(
                    class = seq_len(n_class),
                    age = seq_len(n_age),
                    tow = seq_along(fishing_track)
                )
            ),
            catch_by_geometry = stars::st_as_stars(
                list(catch = catch_by_geo),
                dimensions = stars::st_dimensions(
                    class = seq_len(n_class),
                    age = seq_len(n_age),
                    geometry = geometry |> sf::st_geometry()
                )
            ),
            catchability = catchability,
            area = area
        )
    )
}
