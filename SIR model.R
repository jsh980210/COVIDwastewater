library(tidyverse)
library(pomp)
options(stringsAsFactors=FALSE)
stopifnot(packageVersion("pomp")>="3.0")
set.seed(1350254336)


library(tidyverse)
library(pomp)

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")


sir_init <- Csnippet("
  S = nearbyint(eta*N);
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

sir_dmeas <- Csnippet("
  lik = dbinom(reports,H,rho,give_log);
  ")

sir_rmeas <- Csnippet("
  reports = rbinom(H,rho);
  ")


read.csv(file = "~/UVA research R file/new haven DAY cases.csv") %>%
  select(day,reports=cases) %>% #reports column is labelled "cases" in csv, hence why it is redefined here
  filter(day<=81) %>%
  pomp(
    times="day",t0=0, #we want our times to be in unit day, starting at 0
    rprocess=euler(sir_step,delta.t=1/24), #want delta.t to be smaller than times unit of day. 24 hours in a day
    rinit=sir_init,
    rmeasure=sir_rmeas,
    dmeasure=sir_dmeas,
    accumvars="H", 
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","mu_IR","eta","rho","N"),
    params=c(Beta=15,mu_IR=0.5,rho=0.5,eta=0.06,N=25000)
  ) -> measSIR

measSIR %>%
  simulate(
    nsim=20,format="data.frame",include.data=TRUE
  ) -> sims

sims %>%
  ggplot(aes(x=day,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

