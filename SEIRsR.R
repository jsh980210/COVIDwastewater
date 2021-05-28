library(tidyverse)
library(pomp)
options(stringsAsFactors=FALSE)
stopifnot(packageVersion("pomp")>="3.0")
set.seed(1350254336)


library(tidyverse)
library(pomp)

#step function showing how members of the population move between S, I, R_S (recovered but shedding), and R (recovered and no longer shedding)
seirsr_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-mu_EI*dt));
  double dN_IRS = rbinom(I,1-exp(-mu_IRS*dt));
  double dN_RSR = rbinom(RS,1-exp(-mu_RSR*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IRS;
  RS += dN_IRS - dN_RSR;
  R += dN_RSR;
  H += dN_RSR;
")

seirsr_init <- Csnippet("
  S = nearbyint(eta*N);
  E = 0;
  I = 1;
  RS = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

seirsr_dmeas <- Csnippet("
  lik = dbinom(reports,H,rho,give_log);
  ")

seirsr_rmeas <- Csnippet("
  reports = rbinom(H,rho);
  ")

read_csv(file = "~/UVA research R file/new haven DAY cases.csv") %>%
  select(day,reports=cases) %>% #reports column is labelled "cases" in csv, hence why it is redefined here
  filter(day<=81) %>%
  pomp(
    times="day",t0=0, #we want our times to be in unit day, starting at 0
    rprocess=euler(seirsr_step,delta.t=1/24), #want delta.t to be smaller than times unit of day. 24 hours in a day
    rinit=seirsr_init,
    rmeasure=seirsr_rmeas,
    dmeasure=seirsr_dmeas,
    accumvars="H", 
    statenames=c("S","E","I","RS","R","H"),
    paramnames=c("Beta","mu_EI","mu_IRS","mu_RSR","eta","rho","N"),
    params=c(Beta=337.0430, mu_EI=0.2922718, mu_IRS=0.4079919,mu_RSR=0.03094503000,rho=0.5000,eta=0.0103,N=855000)
  ) -> measSEIRSR

#simulation code with estimated parameters
measSEIRSR %>%
  simulate(
    nsim=20,format="data.frame",include.data=TRUE
  ) -> sims

sims %>%
  ggplot(aes(x=day,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

