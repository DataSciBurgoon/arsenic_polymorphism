################################################################################
# polymorphism_gis_analysis.R
#
################################################################################
library(zipcode)
library(ggmap)
library(ggplot2)
library(tidyr)
library(rgeos)
library(sp)
library(parallel)

#Read in the Census data
setwd("../../census_data/ACS_14_5YR_DP05")
us_census_race_ethnicity_data <- read.csv("ACS_14_5YR_DP05.csv", header=TRUE, skip=1,
                                          check.names = FALSE)
setwd("../../gis_analysis/gis_analysis")

#Add in the lat/long values
data("zipcode")
us_census_race_ethnicity_data$Zip <- clean.zipcodes(us_census_race_ethnicity_data$Id2)
us_census_race_ethnicity_data <- merge(us_census_race_ethnicity_data, zipcode, by.x="Zip", by.y = "zip")


#Only want certain columns from the census data
estimate_columns <- c("Estimate; HISPANIC OR LATINO AND RACE - Total population - Hispanic or Latino (of any race) - Mexican",
                      "Estimate; HISPANIC OR LATINO AND RACE - Total population - Hispanic or Latino (of any race) - Puerto Rican",
                      "Estimate; HISPANIC OR LATINO AND RACE - Total population - Hispanic or Latino (of any race) - Cuban",
                      "Estimate; HISPANIC OR LATINO AND RACE - Total population - Hispanic or Latino (of any race) - Other Hispanic or Latino",
                      "Estimate; HISPANIC OR LATINO AND RACE - Total population - Not Hispanic or Latino - White alone",
                      "Estimate; HISPANIC OR LATINO AND RACE - Total population - Not Hispanic or Latino - Black or African American alone",
                      "Estimate; HISPANIC OR LATINO AND RACE - Total population - Not Hispanic or Latino - American Indian and Alaska Native alone",
                      "Estimate; HISPANIC OR LATINO AND RACE - Total population - Not Hispanic or Latino - Asian alone",
                      "Estimate; HISPANIC OR LATINO AND RACE - Total population - Not Hispanic or Latino - Native Hawaiian and Other Pacific Islander alone",
                      "Estimate; HISPANIC OR LATINO AND RACE - Total population - Not Hispanic or Latino - Some other race alone",
                      "Estimate; HISPANIC OR LATINO AND RACE - Total population - Not Hispanic or Latino - Two or more races")
all_necessary_columns <- c(estimate_columns, "latitude", "longitude")
us_census_race_ethnicity_data_trimmed <- us_census_race_ethnicity_data[, which(colnames(us_census_race_ethnicity_data) %in% all_necessary_columns)]

#Genotype frequencies in each population
genotype_frequencies <- read.table("rs11191439.txt", sep="\t", header=TRUE)
weighted_avg_genome_freqs <- by(genotype_frequencies, genotype_frequencies$Larger.Group, function(x) weighted.mean(x$C_Freq, x$Count), simplify=FALSE)
global_average <- weighted.mean(genotype_frequencies$C_Freq, genotype_frequencies$Count)

#Now I need to put these genotype frequencies into the right order, and use the global average when we have no other information
genotype_ordered <- c(weighted_avg_genome_freqs$Mexican, 
                      weighted_avg_genome_freqs$`Puerto Rican`,
                      global_average,
                      global_average,
                      weighted_avg_genome_freqs$White,
                      weighted_avg_genome_freqs$African,
                      global_average,
                      weighted_avg_genome_freqs$Asian,
                      global_average,
                      global_average,
                      global_average)

#census_genotype_freqs <- us_census_race_ethicity_data_trimmed[, 1:11] * t(as.matrix(genotype_ordered))

prod_fun <- function(x, y){
  x * y
}

t_census_genotype_freqs <- apply(as.matrix(us_census_race_ethnicity_data_trimmed[, 1:11]), 
                               1, 
                               prod_fun, 
                               y=t(as.matrix(genotype_ordered)))

#Want to keep this so that the rows are the zip codes
census_genotype_freqs <- t(t_census_genotype_freqs)

#Aggregate the number of genetically susceptible people by zipcode
agg_census_genotype_by_latlong <- rowSums(census_genotype_freqs)

#Aggregate the population for each zipcode
agg_census_total_population_by_latlong <- rowSums(as.matrix(us_census_race_ethnicity_data_trimmed[, 1:11]))

#Add back in the geocoordinates
agg_census_genotype_by_latlong <- cbind(susc_individuals = agg_census_genotype_by_latlong, 
                                        latitude = us_census_race_ethnicity_data_trimmed$latitude,
                                        longitude = us_census_race_ethnicity_data_trimmed$longitude)

agg_census_total_population_by_latlong <- cbind(population = agg_census_total_population_by_latlong, 
                                        latitude = us_census_race_ethnicity_data_trimmed$latitude,
                                        longitude = us_census_race_ethnicity_data_trimmed$longitude)

prop_at_risk_census_by_latlong <- data.frame(agg_census_genotype_by_latlong)$susc_individuals / data.frame(agg_census_total_population_by_latlong)$population

prop_at_risk_census_by_latlong <- cbind(proportion = prop_at_risk_census_by_latlong,
                                        latitude = us_census_race_ethnicity_data_trimmed$latitude,
                                        longitude = us_census_race_ethnicity_data_trimmed$longitude)

prop_at_risk_census_by_latlong <- as.data.frame(prop_at_risk_census_by_latlong)

#Let's map where these susceptible individuals live
us_map <- get_map("united states", zoom=4)
puerto_rico_map <- get_map("puerto rico", zoom=9)
nc_map <- get_map("north carolina", zoom=7)

#NC bounding box
#34.996, -84.33
#33.84, -78.54
#36.55, -75.85
#36.58, -81.68

agg_census_genotype_by_latlong <- as.data.frame(agg_census_genotype_by_latlong)

agg_census_genotype_by_latlong_threshold <- agg_census_genotype_by_latlong[which(agg_census_genotype_by_latlong$susc_individuals > 5000), ]

nc_data <- subset(agg_census_genotype_by_latlong, 
                  -84.33 <= longitude & longitude <= -75.85 &
                    33.84 <= latitude & latitude <= 36.58)

ggmap(us_map) + geom_point(
  aes(x=longitude, y=latitude, show_guide = TRUE, colour=susc_individuals), 
  data=as.data.frame(agg_census_genotype_by_latlong), alpha=.30, na.rm = T, size=.5)  + 
  scale_color_gradient(low="beige", high="dark red")

png("us_susceptible_individuals_map.png", height=700, width=700)
agg_census_genotype_by_latlong <- as.data.frame(agg_census_genotype_by_latlong)
agg_census_genotype_by_latlong_threshold <- agg_census_genotype_by_latlong[which(agg_census_genotype_by_latlong$susc_individuals > 5000), ]
ggmap(us_map) + geom_point(
  aes(x=longitude, y=latitude, show_guide = TRUE, colour=susc_individuals), 
  data=agg_census_genotype_by_latlong_threshold, alpha=.30, na.rm = T, size=3)  + 
  scale_color_gradient(low="red", high="dark red")
dev.off()


ggmap(puerto_rico_map) + geom_point(
  aes(x=longitude, y=latitude, show_guide = TRUE, colour=susc_individuals), 
  data=as.data.frame(agg_census_genotype_by_latlong), alpha=.50, na.rm = T, size=3)  + 
  scale_color_gradient(low="red", high="dark red")

ggmap(nc_map) + geom_point(
  aes(x=longitude, y=latitude, show_guide = TRUE, colour=susc_individuals), 
  data=nc_data, alpha=.50, na.rm = T, size=3)  + 
  scale_color_gradient(low="beige", high="dark red")

gg_us_map <- ggmap(us_map, extent='device')
gg_pr_map <- ggmap(puerto_rico_map, extent='device')
gg_nc_map <- ggmap(nc_map, extent='device')

ggmap(us_map, extent='device', maprange = FALSE) + 
  geom_density2d(data=agg_census_genotype_by_latlong_threshold,
                 aes(x=longitude, y=latitude),
                 size=0.3) + 
  stat_density2d(
    aes(x = longitude, y = latitude, fill = ..level.., alpha=..level..),
    size = 2, bins = 15, data = agg_census_genotype_by_latlong_threshold,
      geom = "polygon") +
  scale_fill_gradient(low="beige", high="blue") +
  scale_alpha(range = c(.4, .75), guide = FALSE) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))

ggmap(puerto_rico_map, extent='device', maprange = FALSE) + 
  geom_density2d(data=as.data.frame(agg_census_genotype_by_latlong),
                 aes(x=longitude, y=latitude),
                 size=0.3) + 
  stat_density2d(
    aes(x = longitude, y = latitude, fill = ..level.., alpha=..level..),
    size = 2, bins = 10, data = as.data.frame(agg_census_genotype_by_latlong),
    geom = "polygon") +
  scale_fill_gradient(low="beige", high="blue") +
  scale_alpha(range = c(.4, .75), guide = FALSE) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))

ggmap(nc_map, extent='device', maprange = FALSE) + 
  geom_density2d(data=nc_data,
                 aes(x=longitude, y=latitude),
                 size=0.3) + 
  stat_density2d(
    aes(x = longitude, y = latitude, fill = ..level.., alpha=..level..),
    size = 2, bins = 4, data = nc_data,
    geom = "polygon") +
  scale_fill_gradient(low="beige", high="blue") +
  scale_alpha(range = c(.4, .75), guide = FALSE) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))


# Based on Beebe-Dimmer, et al (http://ehjournal.biomedcentral.com/articles/10.1186/1476-069X-11-43)
# Odds ratio for bladder cancer increases 1.7x for rs11191439 per each 1ug/L 
# increase in arsenic in the water

# Bringing in the USGS data on arsenic in the groundwater through 2001 -- it's
# a bit dated, but it's the best data we have available to us at this time
usgs_arsenic_data <- read.table("arsenic_nov2001_usgs.txt", sep="\t", header=TRUE)
usgs_arsenic_data <- usgs_arsenic_data[, c(10:12)]
colnames(usgs_arsenic_data) <- c("concentration", "latitude", "longitude")
usgs_arsenic_data$longitude <- -1 * usgs_arsenic_data$longitude
usgs_geospatial_odds_ratio <- usgs_arsenic_data$concentration * 1.7
usgs_geospatial_odds_ratio <- cbind(odds_ratio = usgs_geospatial_odds_ratio, 
                                    latitude = usgs_arsenic_data$latitude,
                                    longitude = usgs_arsenic_data$longitude)

usgs_geospatial_odds_ratio <- as.data.frame(usgs_geospatial_odds_ratio)

usgs_latlong <- usgs_geospatial_odds_ratio[, 2:3]
census_latlong <- prop_at_risk_census_by_latlong[, 2:3]

#set1sp <- SpatialPoints(usgs_latlong)
#set2sp <- SpatialPoints(census_latlong)



#This next step takes a LONG time to run
#set1$nearest_in_set2 <- apply(gDistance(set1sp, set2sp, byid=TRUE), 1, which.min)


library(geosphere)

# create distance matrix
mat <- distm(usgs_geospatial_odds_ratio[,c("longitude", "latitude")], census_latlong[,c("longitude", "latitude")], fun=distCosine)


# assign the name to the point in list1 based on shortest distance in the matrix
#list1$locality <- list2$locality[apply(mat, 1, which.min)]
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
mat_min_row <- parRapply(cl, mat, which.min)
stopCluster(cl)

usgs_x_census_latitude <- census_latlong$latitude[mat_min_row]
usgs_x_census_longitude <- census_latlong$longitude[mat_min_row]
usgs_prop_at_risk <- prop_at_risk_census_by_latlong$proportion[mat_min_row]

usgs_concentration_x_census <- cbind(concentration = usgs_arsenic_data$concentration,
                                     latitude=usgs_x_census_latitude,
                                     longitude=usgs_x_census_longitude)

usgs_concentration_x_census <- as.data.frame(usgs_concentration_x_census)

png("usgs_arsenic_ground_water_concentrations.png", width=700, height=700)
ggmap(us_map) + geom_point(
  aes(x=longitude, y=latitude, show_guide = TRUE, colour=log10(concentration)), 
  data=usgs_concentration_x_census, alpha=0.8, na.rm = T)  + 
  scale_color_gradient(low="yellow", high="dark red")
dev.off()

#http://stats.stackexchange.com/questions/131416/converting-adjusted-odds-ratios-to-its-rr-counterpart
#Relative Risk=Odds Ratio/((1–p0)+(p0∗Odds Ratio))
#PAR: PAR = Pe*(RRe-1)/([1 + Pe*(RRe-1)])
p0 <- 0.437 #from Beebe-Dimmer, et al
usgs_geospatial_rr <- (usgs_geospatial_odds_ratio$odds_ratio) / ((1-p0)+(p0*usgs_geospatial_odds_ratio$odds_ratio))
usgs_geospatial_par <- (usgs_prop_at_risk * (usgs_geospatial_rr-1))/(1 + usgs_prop_at_risk * (usgs_geospatial_rr - 1))
usgs_geospatial_par_latlong <- cbind(par = usgs_geospatial_par,
                                     latitude = usgs_x_census_latitude,
                                     longitude = usgs_x_census_longitude)
                          
ggmap(us_map, extent='device', maprange = FALSE) + 
  geom_density2d(data=as.data.frame(usgs_geospatial_par_latlong),
                 aes(x=longitude, y=latitude),
                 size=0.3) + 
  stat_density2d(
    aes(x = longitude, y = latitude, fill = ..level.., alpha=..level..),
    size = 2, bins = 10, data = as.data.frame(usgs_geospatial_par_latlong),
    geom = "polygon") +
  scale_fill_gradient(low="beige", high="blue") +
  scale_alpha(range = c(.4, .75), guide = FALSE) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))

usgs_geospatial_par_latlong <- as.data.frame(usgs_geospatial_par_latlong)
usgs_geospatial_par_incidence_latlong <- cbind(par_incidence = usgs_geospatial_par_latlong$par * data.frame(agg_census_total_population_by_latlong)$population[mat_min_row],
                                               latitude = usgs_x_census_latitude,
                                               longitude = usgs_x_census_longitude)

#Population attributable risk incidence map
ggmap(us_map, extent='device', maprange = FALSE) + 
  geom_density2d(data=as.data.frame(usgs_geospatial_par_incidence_latlong),
                 aes(x=longitude, y=latitude),
                 size=0.3) + 
  stat_density2d(
    aes(x = longitude, y = latitude, fill = ..level.., alpha=..level..),
    size = 2, bins = 10, data = as.data.frame(usgs_geospatial_par_incidence_latlong),
    geom = "polygon") +
  scale_fill_gradient(low="beige", high="blue") +
  scale_alpha(range = c(.4, .75), guide = FALSE) +
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))


us_map2 <- get_map("united states", zoom=4, maptype="hybrid")
png("population_attributable_risk_incidence_cases_us-wide.png", width=3000, height=3000, res=300)
ggmap(us_map2) + geom_point(
  aes(x=longitude, y=latitude, show_guide = TRUE, colour=par_incidence, size=par_incidence), 
  data=as.data.frame(usgs_geospatial_par_incidence_latlong), alpha=0.8, na.rm = T)  + 
  scale_color_gradient(low="light blue", high="dark red")
dev.off()

png("population_attributable_risk_incidence_cases_us-wide_state_boundaries.png", width=3000, height=3000, res=300)
ggmap(us_map) + geom_point(
  aes(x=longitude, y=latitude, show_guide = TRUE, colour=par_incidence, size=par_incidence), 
  data=as.data.frame(usgs_geospatial_par_incidence_latlong), alpha=0.8, na.rm = T)  + 
  scale_color_gradient(low="light blue", high="dark red")
dev.off()


#Where is the PAR the highest?
hist(as.data.frame(usgs_geospatial_par_incidence_latlong)$par_incidence)
cutoff <- quantile(as.data.frame(usgs_geospatial_par_incidence_latlong)$par_incidence, probs=0.75, na.rm = TRUE)
usgs_geospatial_par_incidence_latlong <- as.data.frame(usgs_geospatial_par_incidence_latlong)
usgs_geospatial_par_incidence_latlong_threshold <- usgs_geospatial_par_incidence_latlong[which(usgs_geospatial_par_incidence_latlong$par_incidence >= cutoff), ]


ggmap(us_map) + geom_point(
  aes(x=longitude, y=latitude, show_guide = TRUE, colour=par_incidence, size=par_incidence), 
  data=usgs_geospatial_par_incidence_latlong_threshold, alpha=0.8, na.rm = T)  + 
  scale_color_gradient(low="light blue", high="dark red")

png("par_map_us_arsenic_groundwater.png", height=700, width=700)
ggmap(us_map) + geom_point(
  aes(x=longitude, y=latitude, show_guide = TRUE, colour=par_incidence), 
  data=usgs_geospatial_par_incidence_latlong_threshold, alpha=0.8, na.rm = T, size=3)  + 
  scale_color_gradient(low="orange", high="purple")
dev.off()

#Posterior probability of bladder cancer in adults with arsenic exposure > 3.72ppb:
#Note: the prior is the US bladder cancer incidence, which includes the entire US population
# thus it's likely to be an underestimate of the true prior for the genotype.
prior_prob_bladder_cancer <- 20.1/100000 #http://seer.cancer.gov/statfacts/html/urinb.html on August 29, 2016
p_arsenic_given_bladder_cancer <- 0.70  #http://ehjournal.biomedcentral.com/articles/10.1186/1476-069X-11-43
denominator <- (prior_prob_bladder_cancer * p_arsenic_given_bladder_cancer) + (0.30 * (1 - prior_prob_bladder_cancer))
posterior_bladder_cancer_given_arsenic <- (prior_prob_bladder_cancer * p_arsenic_given_bladder_cancer) / denominator
bayes_factor <- (posterior_bladder_cancer_given_arsenic / (1 - posterior_bladder_cancer_given_arsenic)) / (prior_prob_bladder_cancer/ (1 - prior_prob_bladder_cancer))

posterior_bladder_cancer_given_arsenic * 100000 #incidence per 100,000 people is 112.5

#Let's redo this posterior analysis, but this time we're going to add in some uncertainty
#And this will be for lifetime cancer risk:
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# THE MODEL.
bladder_cancer_modelString = "
data {
  int<lower=0> N;     //number of items
  int y;           // y number of successes
}
parameters {
  real <lower=0, upper=1> theta;
}
model {
  theta ~ beta(0.000201*100000, (1-0.000201)*100000);
  y ~ binomial(N, theta);
}
"

#So in this model, I decided that we actually DON'T know what the prior actually
#should be. I'm using a flat prior here.

# THE MODEL.
bladder_cancer_modelString = "
data {
  int<lower=0> N;     //number of items
  int y;           // y number of successes
}
parameters {
  real <lower=0, upper=1> theta;
  real<lower=0,upper=1> lambda; // prior mean chance of success
  real<lower=0.1> kappa; // prior count
}
transformed parameters {
  real<lower=0> alpha; // prior success count
  real<lower=0> beta; // prior failure count
  alpha <- lambda * kappa;
  beta <- (1 - lambda) * kappa;
}
model {
  lambda ~ uniform(0,1); // hyperprior
  kappa ~ pareto(0.1,1.5); // hyperprior
  theta ~ beta(alpha,beta);
  y ~ binomial(N, theta);
}
"

writeLines(bladder_cancer_modelString , con="TEMPmodel.txt" )
stanDso <- stan_model( model_code=bladder_cancer_modelString )

N <- 20
y <- 14
dataList <- list( y = y , N = N) 

disease_allele_stanFit <- sampling( object = stanDso , data = dataList , chains = 3 , iter = 5000 , 
                                    warmup = 200 , thin = 1,
                                    control=list(adapt_delta=0.99))
stan_hist(disease_allele_stanFit)

posterior_dist <- extract(disease_allele_stanFit)[[1]]
mean(posterior_dist)

quantile(posterior_dist, c(0.05)) #boundary on 95% HDI
max(posterior_dist) #upper boundary on 95% HDI

#95% HDI: [0.52, 0.96]; mean 0.69

#Let's do the same for the ancestral allele
N <- 102
y <- 41
dataList <- list( y = y , N = N) 

ancestor_allele_stanFit <- sampling( object = stanDso , data = dataList , chains = 3,
                                     iter = 5000 , warmup = 200 , thin = 1,
                                     control=list(adapt_delta=0.99))
stan_hist(ancestor_allele_stanFit)

posterior_dist <- extract(ancestor_allele_stanFit)[[1]]
mean(posterior_dist)

quantile(posterior_dist, c(0.05)) #boundary on 95% HDI
max(posterior_dist) #upper boundary on 95% HDI

#95% HDI: [0.32, 0.59]; mean 0.40

#Keep in mind that the posteriors are sensitive to differences in the N values. 
#If you have a larger N, then the prior is weighted less, and that has a huge
#influence. So I chose to keep the N values constant, and change the y values
#accordingly.

#Posterior odds ratio
#3.34
(.69/(1-.69))/(.40/(1-.40))


###################
#Bayes Analysis 2
# Going out on a limb here...based on a study from NCI
# it said 20% greater incidence in a New England sample
# when exposed to arsenic in their drinking water compared to US average
# http://jnci.oxfordjournals.org/content/108/9/djw099.abstract

# So in this model, I'm going to assume that the prior probability is like
# 22%. 

# THE MODEL.
ne_prior_bladder_cancer_modelString = "
data {
  int<lower=0> N;     //number of items
  int y;           // y number of successes
}
parameters {
  real <lower=0, upper=1> theta;
}
model {
  theta ~ beta(1.15, 4);
  y ~ binomial(N, theta);
}
"

writeLines(ne_prior_bladder_cancer_modelString , con="TEMPmodel.txt" )
ne_prior_stanDso <- stan_model(model_code=ne_prior_bladder_cancer_modelString )

N <- 20
y <- 14
dataList <- list( y = y , N = N) 

disease_allele_stanFit <- sampling( object = ne_prior_stanDso , data = dataList , chains = 3 , iter = 5000 , 
                                    warmup = 200 , thin = 1,
                                    control=list(adapt_delta=0.99))
stan_hist(disease_allele_stanFit)

posterior_dist <- extract(disease_allele_stanFit)[[1]]
mean(posterior_dist)

quantile(posterior_dist, c(0.05)) #boundary on 95% HDI
max(posterior_dist) #upper boundary on 95% HDI

#95% HDI: [0.44, 0.88]; mean 0.60

#Let's do the same for the ancestral allele
N <- 102
y <- 41
dataList <- list( y = y , N = N) 

ancestor_allele_stanFit <- sampling( object = ne_prior_stanDso , data = dataList , chains = 3,
                                     iter = 5000 , warmup = 200 , thin = 1,
                                     control=list(adapt_delta=0.99))
stan_hist(ancestor_allele_stanFit)

posterior_dist <- extract(ancestor_allele_stanFit)[[1]]
mean(posterior_dist)

quantile(posterior_dist, c(0.05)) #boundary on 95% HDI
max(posterior_dist) #upper boundary on 95% HDI

#95% HDI: [0.32, 0.56]; mean 0.39

#Keep in mind that the posteriors are sensitive to differences in the N values. 
#If you have a larger N, then the prior is weighted less, and that has a huge
#influence. So I chose to keep the N values constant, and change the y values
#accordingly.

#Posterior odds ratio
#2.35
(.60/(1-.60))/(.39/(1-.39))


#Number of wells in the US that have 3ppm or more arsenic based on USGS data 33%
length(which(usgs_arsenic_data$concentration >= 3)) / length(usgs_arsenic_data$concentration)

library(gRain)

yn <- c("yes", "no")
races <- c("mexican", "puerto_rican", "cuban", "other_latino", "white", "black",
           "native", "asian", "hawaiian_pacific", "other")
r1 <- cptable(~race, values=c(rep(.1, 10)), levels=races)
g1 <-  cptable(~genotype:race, values=c(.07, .93, 
                                        .18, .82,
                                        round(global_average,2), 1-round(global_average,2),
                                        round(global_average,2), 1-round(global_average,2),
                                        .10, .90,
                                        .13, .87,
                                        round(global_average,2), 1-round(global_average,2),
                                        round(weighted_avg_genome_freqs$Asian,2), 1-round(weighted_avg_genome_freqs$Asian,2),
                                        round(global_average,2), 1-round(global_average,2),
                                        round(global_average,2), 1-round(global_average,2)),
               levels=yn)
w1 <- cptable(~arsenic_water, values=c(.33, .67), levels=yn)
c1 <- cptable(~cancer|genotype:arsenic_water, 
              values=c(.70, .30, .40, .60, .41, .59, .36, .64),
              levels=yn)

plist <- compileCPT(list(r1, g1, w1, c1))
arsenic_cancer_bn <- grain(plist)

querygrain(setEvidence(arsenic_cancer_bn, evidence=list(race="asian", arsenic_water="yes")))

#Posterior probability of bladder cancer in adults with arsenic exposure > 3.72ppb:
#Note: the prior is the US bladder cancer incidence, which includes the entire US population
# thus it's likely to be an underestimate of the true prior for the genotype.
prior_prob_bladder_cancer <- 20.1/100000 #http://seer.cancer.gov/statfacts/html/urinb.html on August 29, 2016
p_arsenic_given_bladder_cancer <- 0.70  #http://ehjournal.biomedcentral.com/articles/10.1186/1476-069X-11-43
denominator <- (prior_prob_bladder_cancer * p_arsenic_given_bladder_cancer) + (0.30 * (1 - prior_prob_bladder_cancer))
posterior_bladder_cancer_given_arsenic <- (prior_prob_bladder_cancer * p_arsenic_given_bladder_cancer) / denominator
bayes_factor <- (posterior_bladder_cancer_given_arsenic / (1 - posterior_bladder_cancer_given_arsenic)) / (prior_prob_bladder_cancer/ (1 - prior_prob_bladder_cancer))

posterior_bladder_cancer_given_arsenic * 100000 #incidence per 100,000 people

