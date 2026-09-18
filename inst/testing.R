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







n_class<- 3
n_age<- 5
n_cohort<- 10
geometry<- mini_banquereau

recruitment<- sim_recruitment(
    n_class = n_class,
    n_cohort = n_cohort,
    geometry = geometry,
    mean_recruitment_composition = exp(-seq_len(n_class)),
    mean_recruit_abundance = 1000,
    p_composition_cohort = c(0.2, 3),
    p_composition_geometry = c(0.2, 10),
    p_abundance_cohort = c(0.2, 3),
    p_abundance_geometry = c(0.2, 10)
)
growth<- sim_growth(
    n_class = n_class,
    n_age = n_age,
    n_cohort = n_cohort,
    geometry = geometry,
    mean_growth_rate = 0.3,
    p_growth_class = c(0.2, 1),
    p_growth_age = c(0.2, 3),
    p_growth_year = c(0.2, 3),
    p_growth_geometry = c(0.2, 10)
)
mortality<- sim_mortality(
    n_class = n_class,
    n_age = n_age,
    n_cohort = n_cohort,
    geometry = geometry,
    mean_fishing_mortality = 0.08,
    fishing_selectivity = seq(0.05, 0.95, length.out = n_class) |> round(2),
    mean_natural_mortality = 0.1,
    p_fishing_year = c(0.2, 3),
    p_fishing_geometry = c(0.2, 10),
    p_natural_class = c(0.2, 3),
    p_natural_age = c(0.2, 3),
    p_natural_year = c(0.2, 3),
    p_natural_geometry = c(0.2, 10)
)
abundance<- sim_abundance(
    recruitment = recruitment,
    growth = growth,
    mortality = mortality,
    fishing_mean_time = 30,
    fishing_mean_duration = 1000,
    fishing_mean_trips = 10,
    fishing_selectivity = seq(0.05, 0.95, length.out = n_class) |> round(2),
    catchability_mean = 0.1,
    catchability_var = 0.05,
    fishing_area_mean = 0.025,
    fishing_area_var = 0.01
)

