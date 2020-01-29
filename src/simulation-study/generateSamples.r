#' Generate the simulation dataset.
#' 
#' @param samples Number of samples to generate.
#' @param participants Number of subjects in each sample.
#' @param timePoints Vector of times at which measurements will be simulated.
#' @param beta Numeric vector of length 3 containing the true values of the fixed effect parameters \beta.
#' @param sigma Standard deviation of i.i.d. Gaussian response noise \sigma
#' @param alpha Numeric vector of length 3 containing the true values of missingness mechanism parameters \alpha.
#' @param gamma Parameter linking the missingness mechanism to the random effect \gamma.
#' 
#' @return Dataframe with 
generateSamples <- function(samples = 10,
                            participants = 50,
                            timePoints = seq(0, 3, len=5),
                            beta = c(-1.2, 0.5, -1.5),
                            sigma = sqrt(0.5),
                            alpha = c(1.6, -0.5, 2.5),
                            gamma = 1.5,
                            seed = 666) {
  #initialize random number generation
  set.seed(seed = seed)
  
  if(participants %% 2) flog.error('The simulation assumes an even number of participants. Results may be nonsensical.')
  
  #generate dataset
  #first, non-time-dependent part
  df <- data.frame(
    sample = rep(1:samples, each = participants),
    subject = rep(1:participants, times = samples),
    randIntercept = as.vector(t(rmvnorm(
      n = samples*participants/2,
      mean = c(1.35, -1.35),
      sigma = diag(c(0.2^2, 0.6^2))
    ))),
    treatment = c(0,1)
  ) %>%
    uncount(length(timePoints)) %>% #this duplicates each row to make space for the repeated observations
    mutate(time = rep(timePoints, times = samples*participants),
           y = beta[1] + beta[2]*time + beta[3]*treatment + randIntercept + rnorm(samples*participants*length(timePoints), sd=sigma),
           p = case_when(
             time == 0 ~ 1,
             TRUE      ~ plogis(alpha[1] + alpha[2]*time + alpha[3]*treatment + gamma*randIntercept)),
           rMNAR = rbinom(n = samples*participants*length(timePoints),
                          size = 1,
                          prob = p)
    ) %>%
    select(-randIntercept) %>% #drop the random intercept column
    
  
  df
}
