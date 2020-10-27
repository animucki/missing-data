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
#' @return Dataframe with result
generateSamples <- function(samples = 10,
                            participants = 100,
                            timePoints = seq(0, 3, len=5),
                            beta = c(-1.2, 0.5, -1.5),
                            sigma = sqrt(0.5),
                            sigma.b = 0.5,
                            rho.b = 0.9,
                            alpha = c(-1, -0.5, 2.5),
                            gamma = 1.5,
                            delta = 1) {
  
  flog.info(paste0('Generating ', samples, ' samples...'))

  argl <- as.list(match.call())[-1]
  flog.debug(paste(names(argl),unlist(argl),sep=": ", collapse = ", "))
  
  if(participants %% 2) flog.warn('The simulation assumes an even number of participants. Results may be nonsensical.')
  
  #generate dataset
  #first, non-time-dependent part
  # - of which, the random intercepts are Gaussian copula s.t. they are marginally
  #   bi1,bi2 ~iid~ U[ -sqrt(3)*sigma.b, sqrt(3)*sigma.b ]
  #   with forall i: corr(bi1, bi2) = rho.b

  # this creates a bivariate sample ~ U[0,1] with correlation rho.b using the Gaussian copula
  uMat <- pnorm(rmultinormal(n = samples*participants, sigma = c(1, rho.b, rho.b, 1)))

  # this is a scaling step
  bMat <- 2 * sqrt(3) * sigma.b * (uMat - 0.5)

  df <- data.frame(
    sample = as.factor(rep(1:samples, each = participants)),
    subject = as.factor(rep(1:participants, times = samples)),
    classIntercept = rmultinomial(n = samples*participants, size=1, prob=c(3,1,2)) %*% c(-1, 0, 1.5),
    randIntercept1 = bMat[,1],
    randIntercept2 = bMat[,2],
    treatment = c(0,1)
  ) %>%
    uncount(length(timePoints)) %>% #this duplicates each row to make space for the repeated observations
    mutate(time = rep(timePoints, times = samples*participants),
           y = beta[1] + beta[2]*time + beta[3]*treatment + randIntercept1 + classIntercept + rnorm(samples*participants*length(timePoints), sd=sigma),
           p = case_when(
             near(time, 0) ~ 1,
             time > 0 ~ pnorm(alpha[1] + alpha[2]*time + alpha[3]*treatment + gamma*randIntercept1 + delta*classIntercept)),
           rMNAR = rbinom(n = samples*participants*length(timePoints),
                          size = 1,
                          prob = p),
           p3 = case_when(
             near(time, 0) ~ 1,
             time > 0 ~ pnorm(alpha[1] + alpha[2]*time + alpha[3]*treatment + gamma*randIntercept2)),
           rMNAR3 = rbinom(n = samples*participants*length(timePoints),
                          size = 1,
                          prob = p3),
           p4 = case_when(
             near(time, 0) ~ 1,
             time > 0 ~ pnorm(alpha[1] + alpha[2]*time + alpha[3]*treatment + delta*classIntercept)),
           rMNAR4 = rbinom(n = samples*participants*length(timePoints),
                          size = 1,
                          prob = p4)
    )
  
  #for MAR indicator, match to the non-random proportion of missingness 
  propObserved <- df %>% filter(time>0) %>% pull(rMNAR) %>% mean
  df <- df %>% mutate(rMAR = case_when(
    near(time, 0) ~ 1L,
    time > 0 ~ rbinom(n = samples*participants*length(timePoints),
                      size = 1,
                      prob = propObserved)))

  #create new y variables
  df <- df %>% mutate(
    yMAR = case_when(
      rMAR == 1 ~ y,
      rMAR == 0 ~ NA_real_
    ),
    yMNAR = case_when(
      rMNAR == 1 ~ y,
      rMNAR == 0 ~ NA_real_
    ),
    yMNAR3 = case_when(
      rMNAR3 == 1 ~ y,
      rMNAR3 == 0 ~ NA_real_
    ),
    yMNAR4 = case_when(
      rMNAR4 == 1 ~ y,
      rMNAR4 == 0 ~ NA_real_
    )
  ) %>%
    select(!c(y, classIntercept, p, p3, p4, randIntercept1, randIntercept2)) #drop underlying values

  flog.info('Data generation complete.')
  df
}
