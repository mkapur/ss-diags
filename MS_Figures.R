## Code to make figures 1:3 for Carvalho et al 
## kapurm@uw.edu Wi 2020

require(dplyr)
require(r4ss)
require(ggplot2)
require(patchwork) #devtools::install_github("thomasp85/patchwork")
require(reshape2)
require(tseries)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# baseDir <- "./Base_Case/"
baseDir <- "./_mod/atl_tuna_max/New_Ref"
base_case <- SS_output(baseDir)
## first I saved the Names as a csv
nm <- read.csv("./Names.csv", stringsAsFactors = FALSE)[,c(1,3)]


## for use with MAX
nm <- data.frame("Current.names" = base_case$FleetNames[1:4], 
                 New.names = c(paste0("Fishery",1:3),"Survey_1"))
## functions needed ----

runs.sig3 <- function(x,type="resid") {
  if(type=="resid"){mu = 0}else{mu = mean(x, na.rm = TRUE)} 
  # Average moving range
  mr  <- abs(diff(x - mu))
  amr <- mean(mr, na.rm = TRUE)
  # Upper limit for moving ranges
  ulmr <- 3.267 * amr
  # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
  mr  <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)
  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr / 1.128
  # Calculate control limits
  lcl <- mu - 3 * stdev
  ucl <- mu + 3 * stdev
  if(nlevels(factor(sign(x)))>1){ pvalue = round(tseries::runs.test(factor(sign(x)))$p.value,3)} else {
    pvalue = 0.001  
  }
  
  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}


get_ss3compsTA1.8 <- function(ss3rep,type=c('len','age','size','con'),fleet=NULL,
                              part = 0:2,pick.gender = 0:3,seas = NULL,method = NULL,
                              plotit = TRUE,maxpanel = 1000){
  
  # Check the type is correct and the pick.gender is correct
  is.in <- function (x, y)!is.na(match(x, y))
  if(!is.in(type[1],c('age','len','size','con'))){
    stop('Composition type incorrectly speficied')
  }else{
    if(sum(!is.in(pick.gender,c(0:3)))>0){
      stop('Unrecognised value for pick.gender')
    }
  }
  
  # Select the type of datbase
  dbase <- ss3rep[[paste0(type[1],'dbase')]]
  if(is.null(fleet)) fleet = unique(dbase$Fleet)
  sel <-  is.in(dbase$Fleet,fleet) & is.in(dbase$Part,part)
  if(type[1]!='con')sel <- sel & is.in(dbase$Sexes,pick.gender)
  if(type[1]=='size' & !is.null(method)) sel <- sel & is.in(dbase$method,method)
  if(sum(sel)==0) return()
  dbase <- dbase[sel,]
  if(is.null(seas)){
    seas <- 'comb'
    if(length(unique(dbase$Seas))>1)
      cat('Warning: combining data from multiple seasons\n')
  }
  # create label for partitions
  partitions <- sort(unique(dbase$Part)) # values are 0, 1, or 2
  partition.labels <- c("whole","discarded","retained")[partitions+1]
  partition.labels <- paste("(",paste(partition.labels,collapse="&")," catch)",sep="")
  gender.flag <- type[1]!='con' & max(tapply(dbase$Sexes,
                                             dbase$Fleet,function(x)length(unique(x))))>1
  indx <- paste(dbase$Fleet,dbase$Yr,if(type[1]=='con')dbase$'Lbin_lo' else
    '',if(seas=='sep')dbase$Seas else '')
  if(gender.flag)indx <- paste(indx,dbase$Sexes)
  method.flag <- if(type[1]=='size') length(unique(dbase$method))>1 else FALSE
  if(method.flag)
    indx <- paste(indx,dbase$method)
  uindx <- unique(indx)
  if(length(uindx)==1){
    # presumably the method is meaningless of there's only 1 point,
    # but it's good to be able to have the function play through
    cat('Warning: only one point to plot\n')
    return()
  }
  
  pldat <- matrix(0,length(uindx),10,
                  dimnames=list(uindx,
                                c('Obsmn','Obslo','Obshi','semn','Expmn','Std.res',
                                  'ObsloAdj','ObshiAdj','Fleet','Yr')))
  if(type[1]=='con')pldat <- cbind(pldat,Lbin=0)
  if(gender.flag)pldat <- cbind(pldat,pick.gender=0)
  if(method.flag)pldat <- cbind(pldat,method=0)
  
  # Find the weighting factor for this combination of factors
  for(i in 1:length(uindx)){  # each row of pldat is an individual comp
    subdbase <- dbase[indx==uindx[i],]
    xvar <- subdbase$Bin
    pldat[i,'Obsmn'] <- sum(subdbase$Obs*xvar)/sum(subdbase$Obs)
    pldat[i,'Expmn'] <- sum(subdbase$Exp*xvar)/sum(subdbase$Exp)
    pldat[i,'semn'] <- sqrt((sum(subdbase$Exp*xvar^2)/sum(subdbase$Exp)-
                               pldat[i,'Expmn']^2)/mean(subdbase$N))
    pldat[i,'Obslo'] <- pldat[i,'Obsmn']-2*pldat[i,'semn']
    pldat[i,'Obshi'] <- pldat[i,'Obsmn']+2*pldat[i,'semn']
    pldat[i,'Std.res'] <- (pldat[i,'Obsmn']-pldat[i,'Expmn'])/pldat[i,'semn']
    pldat[i,'Fleet'] <- mean(subdbase$Fleet)
    pldat[i,'Yr'] <- mean(if(seas=='comb')subdbase$Yr else subdbase$Yr.S)
    if(type[1]=='con')pldat[i,'Lbin'] <- mean(subdbase$'Lbin_lo')
    if(gender.flag)
      pldat[i,'pick.gender'] <- mean(subdbase$'Pick_gender')
    if(method.flag)
      pldat[i,'method'] <- mean(subdbase$method)
  }
  Nmult <- 1/var(pldat[,'Std.res'],na.rm=TRUE)
  
  # Find the adjusted confidence intervals
  for(i in 1:length(uindx)){
    pldat[i,'ObsloAdj'] <- pldat[i,'Obsmn']-2*pldat[i,'semn']/sqrt(Nmult)
    pldat[i,'ObshiAdj'] <- pldat[i,'Obsmn']+2*pldat[i,'semn']/sqrt(Nmult)
  }
  
  Nfleet <- length(unique(pldat[,'Fleet']))
  # make plot if requested
  if(plotit){
    plindx <- if(type[1]=='con'){
      paste(pldat[,'Fleet'],pldat[,'Yr'])
    }else{
      pldat[,'Fleet']
    }
    if(gender.flag)plindx <- paste(plindx,pldat[,'pick.gender'])
    if(method.flag)plindx <- paste(plindx,pldat[,'method'])
    uplindx <- unique(plindx)
    
    # Select number of panels
    Npanel <- length(uplindx)
    ## Ian T. 9/25/14: changing from having at least 4 panels to no minimum
    #NpanelSet <- max(4,min(length(uplindx),maxpanel))
    NpanelSet <- min(length(uplindx),maxpanel)
    Nr <- ceiling(sqrt(NpanelSet))
    Nc <- ceiling(NpanelSet/Nr)
    # save current graphical parameters
    par_current <- par()
    # set new parameters
    par(mfrow=c(Nr,Nc),mar=c(2,2,1,1)+0.1,mgp=c(0,0.5,0),oma=c(1.2,1.2,0,0),
        las=1)
    par(cex=1)
    for(i in 1:Npanel){
      subpldat <- pldat[plindx==uplindx[i],,drop=FALSE]
      x <- subpldat[,ifelse(type[1]=='con','Lbin','Yr')]
      plot(x,subpldat[,'Obsmn'],pch='-',
           xlim=if(length(x)>1)range(x) else c(x-0.5,x+0.5),
           ylim=range(subpldat[,c('Obslo','Obshi','ObsloAdj','ObshiAdj','Expmn')],
                      na.rm=TRUE),
           xlab='',ylab='')
      segments(x,subpldat[,'Obslo'],x,subpldat[,'Obshi'],lwd=3)
      arrows(x,subpldat[,'ObsloAdj'],x,subpldat[,'ObshiAdj'],lwd=1,
             length=0.04, angle=90, code=3)
      points(x,subpldat[,'Obsmn'],pch=21,bg='grey80')
      ord <- order(x)
      if(length(x)>1){
        lines(x[ord],subpldat[ord,'Expmn'],col=4)
      }else{
        lines(c(x-0.5,x+0.5),rep(subpldat[,'Expmn'],2),col=4)
      }
      # Lines
      fl <- ss3rep$FleetNames[subpldat[1,'Fleet']]
      yr <- paste(subpldat[1,'Yr'])
      lab <- if(type[1]=='con')ifelse(Nfleet>1,paste(yr,fl),yr) else fl
      if(gender.flag)lab <-
        paste(lab,ifelse(subpldat[1,'pick.gender']==0,'comb','sex'))
      if(method.flag)lab <- paste(lab,'meth',subpldat[1,'method'])
      lab <- paste(lab,partition.labels)
      mtext(lab,side=3,at=mean(x))
    }
    mtext(paste('Mean',ifelse(is.in(type[1],c('len','size')),'length','age')),
          side=2,las=0,outer=TRUE)
    mtext(ifelse(type[1]=='con','Length','Year'),side=1,outer=TRUE)
    # restore previous graphics parameters
    par(mfrow=par_current$mfrow, mar=par_current$mar, mgp=par_current$mgp,
        oma=par_current$oma, las=par_current$las)
  }
  tmp <- matrix(sample(pldat[,'Std.res'],1000*nrow(pldat),replace=TRUE),nrow(pldat))
  confint <- as.vector(quantile(apply(tmp,2,function(x)1/var(x,na.rm=TRUE)),
                                c(0.025,0.975),na.rm=TRUE))
  Output <- c(w=Nmult,lo=confint[1],hi=confint[2])
  Outs <- paste("Francis Weights - ", type[1], ": ", ss3rep$FleetNames[fleet],": ",
                round(Nmult,4), " (",round(confint[1],4),"-",round(confint[2],4),")",
                sep="")
  print(Outs)
  pldat=data.frame(pldat)
  yrs=pldat$Yr
  
  comps_out  = list(ss_out = pldat ,runs_dat = data.frame(Fleet=pldat$Fleet,Fleet_name=ss3rep$FleetNames[pldat$Fleet],Yr=yrs,
                                                          Obs=pldat$Obsmn,Exp=pldat$Expmn))
  
  
  # return(Output)
  return(comps_out)
}






## Re-organize data from raw version (as in GDrive 2020-02-19)

## Figure 1 Data plot ----
## We want this to have the correct fishery names and be in greyscale.



png("Fig1_dataplotA_MAX.png", res = 720, width = 6, height = 8,
    unit = 'in')
SSplotData(base_case, 
           fleetnames = c(nm$New.names),
           subplot = 1,
           fleetcol = grey(seq(0.1,0.5, length.out = length(base_case$FleetNames))))
dev.off()

png("Fig1_dataplotB.png", res = 720, width = 6, height = 8,
    unit = 'in')
SSplotData(base_case, fleetnames = c(nm$New.names), subplot = 2,
           fleetcol = grey(seq(0.1,0.5, length.out = length(nm$New.names))))
dev.off()

## Fig 2 Retrospective plots ----
## correct names and prettier than r4ss 

## make summary of retrospectives, package-free
retroSummary <- list.dirs("./Retrospective/", recursive = FALSE, full.names = TRUE)[grepl("retro*",list.dirs("./Retrospective/", recursive = FALSE, full.names = TRUE))] %>%
  lapply(.,SS_output) %>%
  SSsummarize() 

ssbvals       <- retroSummary$SpawnBio 
ssbvalsLower  <- retroSummary$SpawnBioLower
ssbvalsUpper  <- retroSummary$SpawnBioUpper


## manually remove the vals based on the retro year
for(i in 1:5){
  ssbvals[,paste0('model',i)][ssbvals$Label >= paste0("SSB_", base_case$endyr-i)] <- NA
  ssbvalsLower[,paste0('model',i)][ssbvalsLower$Label >= paste0("SSB_", base_case$endyr-i)] <- NA
  ssbvalsUpper[,paste0('model',i)][ssbvalsUpper$Label >= paste0("SSB_", base_case$endyr-i)] <- NA
}

## melt and bind
meltDat<- ssbvals %>% select(-Label) %>% melt(id = 'Yr') %>% 
  merge(.,  ssbvalsLower %>% select(-Label) %>% melt(id = 'Yr'), 
                                                         by = c('Yr','variable')) %>%
  merge(., ssbvalsUpper %>% select(-Label) %>% melt(id = 'Yr'), 
        by = c('Yr','variable'))

## rename the models
levels(meltDat$variable) <-  unlist(lapply(list.dirs("./Retrospective/", recursive = FALSE, full.names = TRUE)[grepl("retro*",list.dirs("./Retrospective/", recursive = FALSE, full.names = TRUE))], basename))
levels(meltDat$variable)[6] <- 'Base Case'
## custom order factor levels
meltDat$variable <- factor(meltDat$variable, levels(meltDat$variable)[c(6,1:5)])
## rename columns
names(meltDat)[3:5] <- c('SSB','LWR','UPR')


## raw recruits
ggplot(meltDat, aes(x = Yr, y = SSB, color = variable)) +
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  geom_line(lwd = 1.1) +
  geom_ribbon(aes(ymin = LWR, ymax = UPR, fill = variable), color = NA,  alpha = 0.15) +
  # scale_color_brewer(palette  = 'virid') +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  scale_x_continuous(limits = c(1985,2015), expand = c(0,0)) +
  # scale_y_continuous( limits = c(0,500), expand = c(0,0)) +
  labs(x = 'Year', y = 'SSB (mt)', fill = "", color = "")

ggsave(last_plot(), file = "./Fig2_RetroSSB_clean.png",
       dpi = 720, width = 8, height = 6, units = 'in'
       )



## Fig 3 Profile Panel ----

# for r0, comps and indices, with correct fishery names
## also had to save as CSV
r0vec <- seq(4.8,6.8,0.1) ## x axis

recProf <- read.table("./_/Profile.csv", skip = 1, sep = ",", header =TRUE)[,1:5]

compProf <- read.table("./Profile/Profile.csv", skip = 1, sep = ",", header =TRUE)[,7:11] %>% mutate(r0vec)
names(compProf) <- c(nm$New.names[nm$Current.names %in% names(compProf)], 'LNR0') ## swap names

idxProf <- read.table("./Profile/Profile.csv", skip = 1, sep = ",", header =TRUE)[,13:17] %>% mutate(r0vec)
names(idxProf) <- c(nm$New.names[nm$Current.names %in% names(idxProf)], 'LNR0') # ## swap names

p1 <- recProf %>%
  melt(id = 'SR_LN.R0.' ) %>%
ggplot(., aes(x = SR_LN.R0., y = value, col = variable)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = c(0.75,0.75),
        legend.text = element_text(size = 14),
        legend.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  geom_line(lwd = 1.1) +
  scale_color_viridis_d() +
  labs(x =expression("ln(R"[0]*")"), y = "Change in Log-Likelihood", col = "")

p2 <- compProf %>%
  melt(id = 'LNR0' ) %>%
  ggplot(., aes(x = LNR0, y = value, col = variable)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = c(0.75,0.75),
        legend.text = element_text(size = 14),
        legend.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  geom_line(lwd = 1.1) +
  scale_color_viridis_d() +
  labs(x =expression("ln(R"[0]*")"), y = "Change in Log-Likelihood", col = "")

p3 <- idxProf %>%
  melt(id = 'LNR0' ) %>%
  ggplot(., aes(x = LNR0, y = value, col = variable)) +
  theme_classic() + 
  theme(panel.grid = element_blank(), 
        legend.position = c(0.75,0.75),
        legend.text = element_text(size = 14),
        legend.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  geom_line(lwd = 1.1) +
  scale_color_viridis_d() +
  labs(x =expression("ln(R"[0]*")"), y = "Change in Log-Likelihood", col = "")

ggsave((p1 | p2  |p3), file = "./Fig3_LikelihoodPanel.png",
       width = 10, height = 8, dpi = 720, unit = 'in')

## FIG 3B Profile Panel MAX ----
r0vec <- seq(11.8,13.5,0.1) ## x axis

recProf <- read.table("./_mod/atl_tuna_max/Profile.csv", skip = 1, sep = ",", header =TRUE)[,c(1:3,5,6)]

agecompProf <- data.frame("OTB_ITA" = read.table("./_mod/atl_tuna_max/Profile.csv", skip = 1, sep = ",", header =TRUE)[,8]) %>% mutate(r0vec)
names(agecompProf) <- c(paste(nm$New.names[nm$Current.names %in% names(agecompProf)]), 'LNR0') ## swap names

lencompProf <-  data.frame(read.table("./_mod/atl_tuna_max/Profile.csv", skip = 1, sep = ",", header =TRUE)[,10:13]) %>% mutate(r0vec)
names(lencompProf)[4] <- "Medits15_16"
names(lencompProf) <- c(paste(nm$New.names[nm$Current.names %in% names(lencompProf)]), 'LNR0') ## swap names

idxProf <- data.frame("Medits15_16" = read.table("./_mod/atl_tuna_max/Profile.csv", skip = 1, sep = ",", header =TRUE)[,15])   %>% mutate(r0vec)
names(idxProf) <- c(paste(nm$New.names[nm$Current.names %in% names(idxProf)]), 'LNR0') # ## swap names


## I updated these so they are in GGPLOT and match the formatting of previous plots
# 1) In all plots update the name of the CPUE's for "CPUE_1" "CPUE_2" etc
# 2) On runs test panel add  Yaxis label "Log residuals" and X axis "Year"
# 3) Increase resolution of all residuals plots if possible.
p1 <- recProf %>%
  melt(id = 'SR_LN.R0.' ) %>%
  ggplot(., aes(x = SR_LN.R0., y = value, col = variable)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = c(0.25,0.75),
        legend.text = element_text(size = 14),
        legend.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  geom_line(lwd = 1.1) +
  scale_color_viridis_d() +
  labs(x =expression("ln(R"[0]*")"), y = "Change in Log-Likelihood", col = "")

p2 <- lencompProf %>%
  melt(id = 'LNR0' ) %>%
  ggplot(., aes(x = LNR0, y = value, col = variable)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = c(0.25,0.75),
        legend.text = element_text(size = 14),
        legend.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  geom_line(lwd = 1.1) +
  scale_color_viridis_d() +
  labs(x =expression("ln(R"[0]*")"), y = "Change in Log-Likelihood", col = "")

p3 <- idxProf %>%
  melt(id = 'LNR0' ) %>%
  ggplot(., aes(x = LNR0, y = value, col = variable)) +
  theme_classic() + 
  theme(panel.grid = element_blank(), 
        legend.position = c(0.25,0.75),
        legend.text = element_text(size = 14),
        legend.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  geom_line(lwd = 1.1) +
  scale_color_viridis_d() +
  labs(x =expression("ln(R"[0]*")"), y = "Change in Log-Likelihood", col = "")

ggsave((p1 | p2  |p3), 
       file = "./Fig3_LikelihoodPanel_Max.png",
       width = 10, height = 8, dpi = 720, unit = 'in')

## run code to generate objects
# source('./Residuals/ss3residualsdiag.R') ## calls ss3runsFun within it
lendat <- get_ss3compsTA1.8(base_case,type="len") # here size (default is length)

## Fig 4 Runs Test MeanL ----
## facet by fleet, point and line with residuals, green if passed overall,
## alpha rectangle indicating output of runs.sig3
mlRes <- lendat$runs_dat %>% 
  select(-Fleet) %>%
  merge(.,nm, by.x = "Fleet_name", by.y = "Current.names") %>%
  mutate(Fleet =  New.names) %>%
  select(Fleet, Yr, Obs, Exp, -Fleet_name,  -New.names) %>%
  mutate(logResidual =  log(Obs) - log(Exp))

mlLim <- mlRes %>%
  group_by(Fleet) %>%
  summarise(pruns = runs.sig3(logResidual)$p.runs, 
            lwr = runs.sig3(logResidual)$sig3lim[1],
            upr = abs(lwr),
            passFail = pruns > 0.05,
            medResidlt3 = median(abs(logResidual))<3,
            minYr = min(Yr), maxYr = max(Yr)) %>%
  filter(medResidlt3)
mlRunsDat <- merge(mlRes,mlLim,"Fleet") 
  
ggplot(mlRunsDat, aes(x = Yr, y = logResidual)) +
  theme_bw() +
  geom_rect(data =mlLim,  aes(x = minYr, y =upr, ## strangely need this to avoid overplotting
                                xmin = minYr, xmax = maxYr, 
                                ymin = -Inf, 
                                ymax = Inf),
            fill = ifelse(mlLim$passFail,
                          'seagreen3','red' ), alpha = 0.75)  +
  theme(panel.background = element_rect(fill = ifelse(mlRunsDat$medResidlt3,
                                                      'seagreen3','red' )),
        panel.grid = element_blank(), 
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title =  element_text(size = 12),
        legend.position = 'none') +
  
  geom_line(lwd = 0.6) + 
  geom_point(shape = 19, size = 2) +
  geom_rect(data =mlLim,  aes(x = minYr, y =upr, ## strangely need this to avoid overplotting
                              xmin = minYr, xmax = maxYr, ymin = lwr, ymax = upr),
            alpha = 0.5, fill = 'grey22')  +
  facet_wrap(~Fleet, scales = 'free_x', ncol = 5) +
  labs(x = 'Year', y = 'log residual')
ggsave(last_plot(), file = "./Fig4_meanLRunsTest.png",
       width = 10, height = 8, dpi = 720, unit = 'in')

## Figure 5 Runs Test CPUE ----
CPUERes <- base_case$cpue %>% 
  select(-Fleet) %>%
  merge(.,nm, by.x = "Fleet_name", by.y = "Current.names") %>%
  mutate(Fleet =  New.names) %>%
  select(Fleet, Yr, Obs, Exp, -Fleet_name,  -New.names) %>%
  mutate(logResidual =  log(Obs) - log(Exp))

CPUELim <- CPUERes %>%
  group_by(Fleet) %>%
  summarise(pruns = runs.sig3(logResidual)$p.runs, 
            lwr = runs.sig3(logResidual)$sig3lim[1],
            upr = abs(lwr),
            passFail = pruns > 0.05,
            medResidlt3 = median(abs(logResidual))<3,
            minYr = min(Yr), maxYr = max(Yr)) %>%
  filter(medResidlt3)

CPUERunsDat <- merge(CPUERes,CPUELim,"Fleet") 

ggplot(CPUERunsDat, aes(x = Yr, y = logResidual)) +
  theme_bw() + 
  geom_rect(data =CPUELim,  aes(x = minYr, y =upr, ## strangely need this to avoid overplotting
                                xmin = minYr, xmax = maxYr, 
                                ymin = -Inf, 
                                ymax = Inf),
                                fill = ifelse(CPUELim$passFail,
                                              'seagreen3','red' ), alpha = 0.75)  +
  
  theme( panel.grid = element_blank(),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 10),
        axis.title =  element_text(size = 12),
        legend.position = 'none') +
  
  geom_line(lwd = 0.6) + 
  geom_point(shape = 19, size = 2) +
  geom_rect(data =CPUELim,  aes(x = minYr, y =upr, ## strangely need this to avoid overplotting
                              xmin = minYr, xmax = maxYr, ymin = lwr, ymax = upr), 
            alpha = 0.5, fill = 'grey22')  +
  # scale_y_continuous(limits = c(-1.5,1.5)) +
  facet_wrap(~Fleet, scales = 'free_x', ncol = 3) +
  labs(x = 'Year', y = 'log residual')

ggsave(last_plot(), file = "./Fig5_CPUERunsTest.png",
       width = 10, height = 8, dpi = 720, unit = 'in')



## Figure X RecDev Test ----
## pausing likely drop due to KFJ comments

## Figure 6 Jabba Residuals CPUE ----
## no functions needed :)

CPUE_ResidRaw <- base_case$cpue %>%
  select(-Fleet) %>%
  merge(.,nm, by.x = "Fleet_name", by.y = "Current.names") %>%
  mutate(Fleet =  New.names) %>%
  select(-Fleet_name, -New.names) %>%
  mutate(logResidual = log(Obs) - log(Exp)) %>%
  filter(abs(logResidual) < 3) %>%
  select(Yr, Fleet, logResidual)

smooth.res <- CPUE_ResidRaw %>%
  mutate(smooth.res.all = predict(loess(logResidual~Yr))) %>%
  group_by(Yr) %>%
  dplyr::summarise(smoother = mean(smooth.res.all))



p1 <- ggplot(CPUE_ResidRaw, aes(x = Yr, y = logResidual)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
                     legend.text = element_text(size = 16),
                     axis.text = element_text(size = 12),
                     axis.title =  element_text(size = 12)) +
  geom_boxplot(aes(group = Yr)) +
  geom_point(aes(colour = Fleet), size = 2, shape = 19) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_line(data = smooth.res, aes(x = Yr, y = smoother, color = 'loess'),
            lwd = 0.8) +
  scale_color_manual(values = cbbPalette[c(2:7,1)],
                     guide = guide_legend(override.aes = list(
                       linetype = c(rep("blank", 6), "solid"),
                       shape = c(rep(19, 6), NA)))) +
  geom_text(aes(x = 2013, y = 0.75, label = paste0('RMSE = ',with(CPUE_ResidRaw,
                                                                  round(100*sqrt(sum(logResidual^2,na.rm =TRUE)/nrow(CPUE_ResidRaw)),1)),"%")), size = 5) +
  labs(x = 'Year', y = 'log residual', col = "")

ggsave(p1, file = "./Fig6_JABBAResid_CPUE.png",
       width = 10, height = 8, dpi = 720, unit = 'in')

## Fig 6B JABBA RESID max ----
CPUE_ResidRaw <- base_case$cpue %>%
  select(-Fleet) %>%
  merge(.,nm, by.x = "Fleet_name", by.y = "Current.names") %>%
  mutate(Fleet =  New.names) %>%
  select(-Fleet_name, -New.names) %>%
  mutate(logResidual = log(Obs) - log(Exp)) %>%
  filter(abs(logResidual) < 3) %>%
  select(Yr, Fleet, logResidual)

smooth.res <- CPUE_ResidRaw %>%
  mutate(smooth.res.all = predict(loess(logResidual~Yr))) %>%
  group_by(Yr) %>%
  dplyr::summarise(smoother = mean(smooth.res.all))



p1 <- ggplot(CPUE_ResidRaw, aes(x = Yr, y = logResidual)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title =  element_text(size = 12)) +
  geom_boxplot(aes(group = Yr)) +
  geom_point(aes(colour = Fleet), size = 2, shape = 19) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_line(data = smooth.res, aes(x = Yr, y = smoother, color = 'loess'),
            lwd = 0.8) +
  scale_color_manual(values = cbbPalette[c(1,2)])+
  geom_text(check_overlap = TRUE,aes(x = 2013, y = 0.75,
                label = paste0('RMSE = ',with(CPUE_ResidRaw,
                                                                  round(100*sqrt(sum(logResidual^2,na.rm =TRUE)/nrow(CPUE_ResidRaw)),1)),"%")), size = 5) +
  labs(x = 'Year', y = 'log residual', col = "")

ggsave(p1, file = "./Fig6_JABBAResid_CPUE_Max.png",
       width = 10, height = 8, dpi = 720, unit = 'in')

## use mlRes made for runs test plot
smooth.res <- mlRes %>%
  mutate(smooth.res.all = predict(loess(logResidual~Yr))) %>%
  group_by(Yr) %>%
  dplyr::summarise(smoother = mean(smooth.res.all))



p2 <- ggplot(mlRes, aes(x = Yr, y = logResidual)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title =  element_text(size = 12)) +
  geom_boxplot(aes(group = Yr)) +
  geom_point(aes(colour = Fleet), size = 2, shape = 19) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_line(data = smooth.res, aes(x = Yr, y = smoother, color = 'loess'),
            lwd = 0.8) +
  scale_color_manual(values = c(viridis::viridis(5),'black'),
                     guide = guide_legend(override.aes = list(
                       linetype = c(rep("blank", 5), "solid"),
                       shape = c(rep(19, 5), NA)))) +
  geom_text(aes(x = 2013, y = 0.25, 
                label = paste0('RMSE = ',with(mlRes,  round(100*sqrt(sum(logResidual^2,na.rm =TRUE)/nrow(mlRes)),1)),"%")), size = 5) +
  labs(x = 'Year', y = 'log residual', col = "")

ggsave(p2, file = "./Fig6_JABBAResid_meanL.png",
       width = 10, height = 8, dpi = 720, unit = 'in')


ggsave((p1  / p2), file = "./Fig6_JABBAResid_panel.png",
       width = 10, height = 8, dpi = 720, unit = 'in')


## Figure 7 Jitter panel plot ----
## two panels, one with likelihood and other with ssb trend from jittered runs
jit <- read.csv("./max/plots/Likelihoods.csv", sep = " ") %>% 
  filter(Label == "TOTAL") %>%
  reshape2::melt() %>%
  filter(!is.na(value)) %>%
  dplyr::mutate(idx0 = as.numeric(as.character(sub('replist',"",variable))))

# jitbase <- SS_output("./Jitter_Mako/02_SS_NEW_run/")
jitbase <- SS_output("./max/reference/")
p1 <-  jit %>%
   mutate(idx = 1:nrow(.)) %>%
   ggplot(., aes(x = idx, y = value)) +
     theme_bw() + 
   theme( panel.grid = element_blank(),
          strip.text = element_text(size = 16),
          axis.text = element_text(size = 10),
          axis.title =  element_text(size = 12),
          legend.position = 'none') +
   geom_hline(yintercept = jitbase$likelihoods_used[[1]][1], col = 'red', linetype = 'dashed') +
     geom_point() +
   # scale_y_continuous(limits = c(0,155)) +
   labs(x = "Jitter runs at a converged solution", y = "Total Likelihood")
 
  
jit = read.csv("./Jitter_Mako/plots/quants.csv", sep = " ")


## get values from converged runs
jitSummary <- SSgetoutput(keyvec=jit$idx0, getcomp=FALSE, 
                       dirvec="./Jitter_Mako/jitter/", getcovar=F) %>%
  SSsummarize(.)

ssbvals       <- jitSummary$SpawnBio 
ssbvalsLower  <- jitSummary$SpawnBioLower
ssbvalsUpper  <- jitSummary$SpawnBioUpper


BASE_ssb <- jitbase$derived_quants[grep("SSB",jitbase$derived_quants$Label),] %>%
  mutate(Yr = as.numeric(as.character(sub("SSB_","",Label)))) %>%
  filter(!is.na(Yr)) %>% 
  mutate(variable = "Base Case", 
         value.y = Value - 1.96*StdDev, 
         value.x = Value + 1.96*StdDev, value = Value) %>%
  select(Yr, variable, value, value.x,value.y )



## melt and bind
meltDat <- ssbvals %>% select(-Label) %>% melt(id = 'Yr') %>% 
  merge(.,  ssbvalsLower %>% select(-Label) %>% melt(id = 'Yr'), 
        by = c('Yr','variable')) %>%
  merge(., ssbvalsUpper %>% select(-Label) %>% melt(id = 'Yr'), 
        by = c('Yr','variable')) %>%
  bind_rows(., BASE_ssb) #%>%
  # mutate(showLeg = ifelse(variable == 'Base Case',T,F))

## rename columns
names(meltDat)[3:5] <- c('SSB','LWR','UPR')

meltDat$variable <-  factor(meltDat$variable)
# meltDat$variable <- factor(meltDat$variable, rev(levels(meltDat$variable)))
## raw SSB
p2 <-  ggplot(meltDat, aes(x = Yr, y = SSB, color = variable)) +
  theme_bw() + 
  theme( panel.grid = element_blank(),
         strip.text = element_text(size = 16),
         axis.text = element_text(size = 10),
         axis.title =  element_text(size = 12),
         legend.position = 'right') +  geom_line(lwd = 1.1) +
  
  geom_ribbon(aes(ymin = LWR, ymax = UPR, fill = variable), color = NA,  alpha = 0.1) +
  geom_line(data = subset(meltDat, variable == 'Base Case'), lwd = 1.1) +
  # scale_linetype_manual(breaks = 'Base Case', 
  #                       values = c('dashed', rep('solid',69)),
  #                       guide_legend(override.aes = list(colour = viri, 
  #                                                        fill = breaks = c("Base Case")))) +
  scale_color_viridis_d(breaks = c("Base Case")) +
  scale_fill_viridis_d(breaks = c("Base Case"))  +
  scale_x_continuous(limits = c(1950,2020), breaks = seq(1950,2010,10), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = 'Year', y = 'Spawning Biomass (mt)', fill = "", color = "")


  ggsave((p1 |p2), file = "./Fig7_JitterPanel.png",
         width = 10, height = 6, dpi = 720, unit = 'in')
  
  
  ## Fig 7B side by side jitter ----
  ## mako and max, just likelihood
## MAX JITTER
  maxDir <- "./_mod/atl_tuna_max/Reference/"
  max_base <- SS_output(maxDir)
    
  # jit <- read.csv("./_mod/atl_tuna_max/jit_like_Max.csv")
  # names(jit) <- c("idx", "value")
  #   
    ## METHOD for what is in plots/likelihoods
    jit <- read.csv("./_mod/atl_tuna_max/plots/Likelihoods.csv", sep = " ") %>%
    filter(Label == "TOTAL") %>%
    reshape2::melt() %>%
    filter(!is.na(value)) %>%
    dplyr::mutate(idx0 = as.numeric(as.character(sub('replist',"",variable))))

  p1 <-  
    jit %>%
    mutate(idx = 1:nrow(.)) %>%
    ggplot(., aes(x = idx, y = value)) +
    theme_bw() + 
    theme( panel.grid = element_blank(),
           strip.text = element_text(size = 16),
           axis.text = element_text(size = 10),
           axis.title =  element_text(size = 12),
           legend.position = 'none') +
    geom_hline(yintercept = max_base$likelihoods_used[[1]][1], 
               col = 'red', linetype = 'dashed') +
    geom_point() +
    scale_y_continuous(limits = c(0,3000)) +
    labs(x = "Jitter runs at a converged solution",
         y = "Total Likelihood")

  
  ## MAKO jitter
  
  makoDir <- "./_mod/Jitter_Mako/02_SS_NEW_run"
  mako_base <- SS_output(makoDir)

  jit <- read.csv("./_mod/Jitter_Mako/plots/Likelihoods.csv", sep = " ") %>% 
    filter(Label == "TOTAL") %>%
    reshape2::melt() %>%
    filter(!is.na(value)) %>%
    dplyr::mutate(idx0 = as.numeric(as.character(sub('replist',"",variable))))

  p2 <-  jit %>%
    mutate(idx = 1:nrow(.)) %>%
    ggplot(., aes(x = idx, y = value)) +
    theme_bw() + 
    theme( panel.grid = element_blank(),
           strip.text = element_text(size = 16),
           axis.text = element_text(size = 10),
           axis.title =  element_text(size = 12),
           legend.position = '') +
    geom_hline(yintercept = mako_base$likelihoods_used[[1]][1], 
               col = 'red', linetype = 'dashed') +
    geom_point() +
    scale_y_continuous(limits = c(0,155)) +
    labs(x = "Jitter runs at a converged solution",
         y = "Total Likelihood")
require(patchwork)
ggsave(p2 | p1,
       file = "./Fig7B_jitter_max_mako.png",
       width = 12, height = 10, dpi = 500, unit = 'in')



## Fig 10 ASPM panel ----
## load & fill plot list in a loop, first mako then max
plist = list(); idx = 1
for(i in c(1:2)){
  baseDir <- c("./_mod/Jitter_Mako/02_SS_NEW_run","./_mod/atl_tuna_max/ASPM")[i]
  aspmDir <- c("./_mod/ASPMII","./_mod/atl_tuna_max/New_Ref")[i] ## mako, max
  aspm <- SS_output(aspmDir, covar = FALSE)
  base_case <- SS_output(baseDir)
  
  if(i == 1){
    nm <- read.csv("./Names.csv", stringsAsFactors = FALSE)[,c(1,3)]
  } else{
    ## for use with MAX
    nm <- data.frame("Current.names" = base_case$FleetNames[1:4], 
                     New.names = c(paste0("Fishery",1:3),"Survey_1"))
  }
  



 
  

  ## CPUE 5 fits on single plot
  cp0 <- bind_rows(aspm$cpue %>% 
                     merge(.,nm, by.x = "Fleet_name", by.y = "Current.names") %>%
                     filter(Fleet_name == ifelse(i == 2, "Medits15_16", 'S5_EU_ESP_LL')) %>% ## S5 FOR PREVIOUS
                     mutate(model = 'ASPM'),
                   base_case$cpue %>% 
                     merge(.,nm, by.x = "Fleet_name", by.y = "Current.names") %>%
                     filter(Fleet_name == ifelse(i == 2, "Medits15_16", 'S5_EU_ESP_LL')) %>% ## S5 FOR PREVIOUS
                     mutate(model = 'Base Case')) %>%
    mutate(     lower_total = NA, #qnorm(.025, mean = Exp, sd = SE),
                upper_total = NA) %>%  #qnorm(.975, mean = Exp, sd = SE)) %>%
    select(Yr, Obs, Exp, SE, lower_total, upper_total, model)
  
  cp <- rbind(cp0,
              cp0 %>% filter(model == 'ASPM') %>% 
                mutate(model = 'Observed', 
                       lower_total = exp(log(Obs) + qt(0.025, df = 30) * 
                                           SE),
                       upper_total = exp(log(Obs) + qt(0.975, df = 30) * 
                                           SE)))
  
  
  plist[[idx]] <- ggplot(cp, aes(x = Yr, y = Exp, color = model)) +
    theme_bw() + 
    theme( panel.grid = element_blank(),
           strip.text = element_text(size = 16),
           axis.text = element_text(size = 10),
           axis.text.y = element_blank(),
           axis.title =  element_text(size = 12),
           legend.position = c('none')) +  
    scale_color_manual(values = c('dodgerblue2','seagreen4', 'black')) +
    guides(shape = FALSE, linetype = FALSE,
           colour = guide_legend(override.aes =  list(shape = c(NA,NA,1),
                                                      linetype = c("solid","solid","blank")))) +  
    geom_errorbar(data = subset(cp, model == 'Observed'),
                  aes(x = Yr, ymin = lower_total, 
                      ymax = upper_total,
                      group = model, color = model)) +
    geom_point(data = subset(cp, model == 'Observed'), shape = 1,
               aes(x = Yr, y = Obs, color = model)) +
    
    geom_line(data = subset(cp, model != 'Observed'),
              aes(x = Yr, y = Exp, color = model), lwd =1.1) +
    
    # scale_y_continuous(limits = c(0,80)) +
    scale_x_continuous(limits = c(1989,2018), 
                       breaks = seq(1990,2018,10)) +
    labs(x = 'Year', 
         y = paste0('Index ',
                    ifelse(i == 2, "Survey_1",  'CPUE_5')),  
         color = "")
  idx <- idx + 1
  
  ## Derived RelSSB with error and refpt line
  stdtable <- aspm$derived_quants[substring(aspm$derived_quants$Label,1,4)=="SSB_",]  %>%
    mutate(model = 'ASPM', SSBMSY =  aspm$derived_quants$Value[substring(aspm$derived_quants$Label,1,7)
                                                                    =="SSB_MSY"]) %>%
    bind_rows(.,
              base_case$derived_quants[substring(base_case$derived_quants$Label,1,4)=="SSB_",]  %>%
                mutate(model = 'Base Case', SSBMSY =  base_case$derived_quants$Value[substring(base_case$derived_quants$Label,1,7)
                                                                                     =="SSB_MSY"]))
  # stdtable <- stdtable[tolower(stdtable$Label)!="SSB_unfished",]
  # year as the part of the Label string starting with 6th character
  stdtable$Yr <- substring(stdtable$Label,5)
  # filling in Virgin and Initial years as 2 and 1 years prior to following years
  # stdtable$Yr[c(1,2,73,74)] <- as.numeric(stdtable$Yr[3])-(2:1) ## original
  
  stdtable$Yr <- as.numeric(stdtable$Yr)
  stdtable$val <- stdtable$Value
  # assume SSB have normal distribution
  plist[[idx]] <- 
    stdtable %>% 
    filter(!is.na(Yr) ) %>%
    group_by(Yr, model) %>%
    dplyr::summarise(Value = val, 
              SSBMSY = SSBMSY,
              lower = Value - 1.96*StdDev,
              upper = Value + 1.96*StdDev ) %>%
    ggplot(., aes(x = Yr, y = Value, fill = model, color = model)) +
    theme_bw() + 
    theme( panel.grid = element_blank(),
           strip.text = element_text(size = 16),
           axis.text = element_text(size = 10),
           axis.text.y = element_blank(),
           
           axis.title =  element_text(size = 12),
           legend.position = c('none')) +  
    scale_x_continuous(limits = c(1948,2016), breaks = seq(1950,2016,10)) +
    # scale_y_continuous(limits = c(0,2700)) +
    
    scale_color_manual(values = c('dodgerblue2','seagreen4')) +
    scale_fill_manual(values = c('dodgerblue2','seagreen4')) +
    geom_ribbon(aes(x = Yr, ymin = lower, ymax = upper, group = model, fill = model),
                alpha = 0.2 , color = NA) +
    geom_point() +
    
    geom_hline(aes(yintercept = SSBMSY, colour = model ), 
               linetype = 'dashed', lwd = 1.1, alpha = 0.2) +
    geom_text(aes(x = 1960, y = SSBMSY*1.1, label = 'SSB_MSY',color = model), 
              check_overlap = TRUE,
              alpha = 0.5) +
    guides( color = FALSE,
            fill = guide_legend(override.aes =  list(shape = c(19,19),
                                                     color = c('dodgerblue2','seagreen4')))) +
    labs(x = 'Year', y = 'Spawning Biomass (1000 mt)', fill = "")
  idx <- idx + 1
  
  
  ## Recruit TS on single plot
  stdtable <- aspm$derived_quants[substring(aspm$derived_quants$Label,1,5)=="Recr_",]  %>%
    mutate(model = 'ASPM') %>%
    bind_rows(.,
              base_case$derived_quants[substring(base_case$derived_quants$Label,1,5)=="Recr_",]  %>%
                mutate(model = 'Base Case'))
  stdtable <- stdtable[tolower(stdtable$Label)!="recr_unfished",]
  # year as the part of the Label string starting with 6th character
  stdtable$Yr <-substring(stdtable$Label,6)
  # filling in Virgin and Initial years as 2 and 1 years prior to following years
  stdtable$Yr[c(1,2,69,70)] <- as.numeric(stdtable$Yr[3])-(2:1)
  stdtable$Yr <- as.numeric(stdtable$Yr)
  
  # stdtable$Yr <- as.numeric(stdtable$Yr) + yrshift
  bioscale <- 1
  # scaling and calculation of confidence intervals
  v <- stdtable$Value * bioscale
  std <- stdtable$StdDev * bioscale
  # assume recruitments have log-normal distribution
  # from first principals (multiplicative survival probabilities)
  stdtable$logint <- sqrt(log(1+(std/v)^2))
  stdtable$lower <- exp(log(v) - 1.96*stdtable$logint)
  stdtable$upper <- exp(log(v) + 1.96*stdtable$logint)
  
  plist[[idx]] <-  ggplot(stdtable, aes(x = Yr, y = Value, color = model)) +
    theme_bw() + 
    theme( panel.grid = element_blank(),
           strip.text = element_text(size = 16),
           axis.text = element_text(size = 10),
           axis.text.y = element_blank(),
           
           axis.title =  element_text(size = 12),
           legend.position = 'bottom') +  
    scale_x_continuous(limits = c(base_case$startyr,2020),
                       breaks = seq(1950,2020,10)) +
    # scale_y_continuous(limits = c(0,600)) +
    scale_color_manual(values = c('dodgerblue2','seagreen4')) +
    # geom_ribbon(aes(x = Yr, ymin = lower, ymax = upper, group = model, fill = model),
    #             alpha = 0.2 ) +
    geom_errorbar(aes(x = Yr, ymin = lower, ymax = upper, group = model, color = model)) +
    geom_point() +
    labs(x = 'Year', y = 'Age-0 Recruits (1000s)', color = "")
  
  idx <- idx + 1
  
}

  
  ggsave( (plist[[1]]/plist[[2]]/plist[[3]]) | 
            (plist[[4]]/plist[[5]]/plist[[6]]),
         file = "./Fig10_ASPMPanel_mako_Max.png",
         width = 12, height = 12, dpi = 500, unit = 'in')
  
  
  
  