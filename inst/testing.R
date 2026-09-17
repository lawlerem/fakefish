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


plot_growth(foo)
plot_growth(foo, age = 5)
plot_growth(foo, year = 5)

plot_mortality(foo, type = "total")
plot_mortality(foo, type = "fishing")
plot_mortality(foo, type = "natural")

plot_abundance(foo, by = "year")
plot_abundance(foo, by = "cohort")