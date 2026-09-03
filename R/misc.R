# cohort = 5 <-> age = 3 <-> year = 7
# cohort = 5 <-> age = 2 <-> year = 6 
# cohort = 5 <-> age = 1 <-> year = 5
#
# ......-
# .....-.
# ....-..
#

convert_ayc<- function(age, cohort, year) {
    if( missing(age) ) return(.co_y_to_a(cohort, year))
    if( missing(cohort) ) return(.a_y_to_co(age, year))
    if( missing(year) ) return(.a_co_to_y(age, cohort))
}
.a_co_to_y<- function(age, cohort) cohort + age - 1
.a_y_to_co<- function(age, year) year - age + 1
.co_y_to_a<- function(cohort, year) year - cohort + 1