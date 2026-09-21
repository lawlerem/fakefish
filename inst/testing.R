# library(fakefish)
library(tmap)
devtools::load_all()

n_class<- 3
n_age<- 5
n_cohort<- 10
stock<- sim_population(
    n_class = n_class,
    n_age = n_age,
    n_cohort = n_cohort,
    geometry = mini_banquereau,
    mean_recruitment_composition = exp(-seq_len(n_class)),
    mean_recruit_abundance = 1000,
    mean_growth_rate = 0.3,
    mean_natural_mortality = 0.1,
    fishing_mean_time = 40,
    fishing_mean_duration = 5000,
    fishing_mean_trips = 4,
    fishing_selectivity = seq(0.05, 0.95, length.out = n_class) |> round(2),
    catchability_mean = 0.4,
    catchability_var = 0.05,
    fishing_area_mean = 0.05,
    fishing_area_var = 0.01,
    p_composition_cohort = c(0.3, 6),
    p_composition_geometry = c(0.3, 6),
    p_abundance_cohort = c(0.3, 6),
    p_abundance_geometry = c(0.3, 6),
    p_growth_class = c(0.3, 2),
    p_growth_age = c(0.3, 2),
    p_growth_year = c(0.3, 2),
    p_growth_geometry = c(0.3, 2),
    p_natural_class = c(0.3, 2),
    p_natural_age = c(0.3, 2),
    p_natural_year = c(0.3, 2),
    p_natural_geometry = c(0.3, 2)
)



tm_shape(stock$abundance[, , , 10, 2, , drop = TRUE]) +
    tm_fill(fill = "abundance", fill.scale = tm_scale_continuous())

biomass<- stock$abundance |>
    sweep(1, c(0.2, 1, 2), `*`) |>
    (\(x) {names(x)<- "biomass"; x})()
cbiomass<- biomass |>
    stars::st_apply(3:5, sum) |>
    (\(x) {names(x)<- "biomass"; x})()

tm_shape(cbiomass[, , 2, , drop = TRUE]) +
    tm_fill(fill = "biomass", fill.scale = tm_scale_continuous())