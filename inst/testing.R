# library(fakefish)
devtools::load_all()

n_class<- 3
foo<- sim_population(
    n_class = 3,
    n_age = 5,
    n_cohort = 10,
    geometry = mini_banquereau,
    mean_recruitment_composition = exp(-seq_len(n_class)),
    mean_recruit_abundance = 1000,
    mean_growth_rate = 0.3,
    mean_fishing_mortality = 0.08,
    fishing_selectivity = rep(1, n_class),
    mean_natural_mortality = 0.1
)


# plot_growth(foo)
# plot_growth(foo, age = 5)
# plot_growth(foo, year = 5)

# plot_mortality(foo, type = "total")
# plot_mortality(foo, type = "fishing")
# plot_mortality(foo, type = "natural")

# plot_abundance(foo, by = "year")
# plot_abundance(foo, by = "cohort")


fish<- foo$abundance[, , , 9, , drop = TRUE]

selectivity<- c(0.03, 0.3, 1)
target<- fish$abundance |> 
    sweep(1, selectivity, `/`) |>
    apply(3, sum)
track<- seq(10) |>
    lapply(\(i) generate_fishing_track(
        geometry = mini_banquereau,
        target = target,
        mean_time = 30,
        duration = 1000
    )) |>
    do.call(c, args = _)
fishing_effort<- track |> 
    factor(levels = seq(nrow(mini_banquereau))) |>
    table() |>
    as.numeric() |>
    sf::st_sf(
        effort = _,
        geometry = mini_banquereau |> sf::st_geometry()
    )
catch<- fish_track(
    fishing_track = track,
    fish = fish,
    selectivity = selectivity,
    geometry = mini_banquereau,
    catchability_mean = 0.1,
    catchability_var = 0.2,
    area_mean = 1,
    area_var = 1
)