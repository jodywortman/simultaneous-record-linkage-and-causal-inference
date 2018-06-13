#Code for simulations- sections can be commented out or un-commented to run different simulation variations
library(RecordLinkage)
library(xtable)
library(dplyr)
set.seed(5)
data(RLdata10000)
data(RLdata500)
bigid=identity.RLdata10000

lg=cbind(RLdata10000, bigid)
lg=lg[order(bigid),]

oldtab=data.frame(table(lg$bigid))
names(oldtab)=c("bigid", "freq")

withdupes= (merge(oldtab, lg, by.x = "bigid", by.y = "bigid"))

#####################################Smaller large file
# largecleaned = subset(withdupes, freq>1)
# c(length(unique(largecleaned$bigid)), nrow(largecleaned))
# largeother= subset(withdupes, freq==1)
# largeothersamp=sample(largeother$bigid, 2000, replace = FALSE)
# largeother = largeother[largeother$bigid %in% largeothersamp,]
# 
# 
# c(length(unique(largeother$bigid)), nrow(largeother))
# 
# 
# largecleaned$tabid=1:nrow(largecleaned)
# small = largecleaned[duplicated(largecleaned$bigid), ]
# nrow(small)
# length(unique(small$bigid))
# 
# small2=largecleaned[!largecleaned$tabid %in% small$tabid, ]
# nrow(small2)
# length(unique(small2$bigid))
# 
# small=small[,1:9]
# large=rbind(small2[,1:9], largeother)
# 
# length(unique(large$bigid))
# nrow(large)
# nrow(small)
# length(unique(small$bigid))
# 
# largeind=sample(largeother$bigid, 1000, replace = FALSE)
# 
# smallextra = large[large$bigid %in% largeind,]
# small = rbind(small, smallextra)
# 
# large = large[!large$bigid %in% largeind,]
# c(nrow(small), length(unique(small$bigid)))
# c(nrow(large), length(unique(large$bigid)))
# 
# 
# # for(i in 1:nrow(small)){
# #   large$by[large$bigid==small$bigid[i]]=small$by[small$bigid==small$bigid[i]]
# # } ##too slow!
# 
# #faster way
# large$by[large$freq==2]=merge(large, small, by = 'bigid')$by.y
# 
# sum(merge(small, large, by="bigid")$by.x==merge(small, large, by="bigid")$by.y)
# 
# largebase = large
# smallbase = small
###################################################################


##################################original
largecleaned = subset(withdupes, freq>1)
c(length(unique(largecleaned$bigid)), nrow(largecleaned))
largeother= subset(withdupes, freq==1)
c(length(unique(largeother$bigid)), nrow(largeother))


largecleaned$tabid=1:nrow(largecleaned)
small = largecleaned[duplicated(largecleaned$bigid), ]
nrow(small)
length(unique(small$bigid))

small2=largecleaned[!largecleaned$tabid %in% small$tabid, ]
nrow(small2)
length(unique(small2$bigid))

small=small[,1:9]
large=rbind(small2[,1:9], largeother)

length(unique(large$bigid))
nrow(large)
nrow(small)
length(unique(small$bigid))

largeind=sample(largeother$bigid, 1000, replace = FALSE)

smallextra = large[large$bigid %in% largeind,]
small = rbind(small, smallextra)

large = large[!large$bigid %in% largeind,]
c(nrow(small), length(unique(small$bigid)))
c(nrow(large), length(unique(large$bigid)))


# for(i in 1:nrow(small)){
#   large$by[large$bigid==small$bigid[i]]=small$by[small$bigid==small$bigid[i]]
# } ##too slow!

#faster way
large$by[large$freq==2]=merge(large, small, by = 'bigid')$by.y

sum(merge(small, large, by="bigid")$by.x==merge(small, large, by="bigid")$by.y)

largebase = large
smallbase = small
########################################################################


#####################################FS LINKAGE#################
###see what the logdiffs will be
comparemat=merge(small, large, by="by")

comparemat$monthagree=as.numeric(comparemat$bm.x==comparemat$bm.y)
comparemat$dayagree=as.numeric(comparemat$bd.x==comparemat$bd.y)
comparemat$truematch=as.numeric(comparemat$bigid.x==comparemat$bigid.y)
comparemat$firstagree=jarowinkler(as.character(comparemat$fname_c1.x), as.character(comparemat$fname_c1.y))
comparemat$lastagree=jarowinkler(as.character(comparemat$lname_c1.x), as.character(comparemat$lname_c1.y))
#attach(comparemat)
comparemat$simplefirst=as.numeric(comparemat$firstagree>.95)
comparemat$simplelast=as.numeric(comparemat$lastagree>.95)
#attach(comparemat)

comp=comparemat

m=rep(.95, 4)
u=c(sum(comp$simplefirst)/length(comp$firstagree),
    sum(comp$simplelast/length(comp$lastagree)),
    sum(comp$monthagree)/length(comp$monthagree),
    sum(comp$dayagree)/length(comp$dayagree))


logdiff=log(m[1]/u[1], 2)*comp$simplefirst+
  log(m[2]/u[2], 2)*comp$simplelast+
  log((1-m[1])/(1-u[1]), 2)*(1-comp$simplefirst)+
  log((1-m[2])/(1-u[2]), 2)*(1-comp$simplelast)+
  log(m[1]/u[3], 2)*comp$monthagree+
  log(m[2]/u[4], 2)*comp$dayagree+
  log((1-m[1])/(1-u[3]), 2)*(1-comp$monthagree)+
  log((1-m[2])/(1-u[4]), 2)*(1-comp$dayagree)

#####now we can do the rest

###Now causal part
#create simulated data

threshold=sort(unique(logdiff))
threshold = threshold[threshold>0]
nt=length(threshold)
##############################################################



################################JW LINKAGE
# threshold=c(.99, 0.985, .975, .925, .9, .8)
# threshold = sort(threshold)
# nt=length(threshold)
##############################################

nsim=100
ttau=tvarhat=cminvartau=minvartau=dtau=dvarhat=ybar0=ybar1=rep(NA, nsim)
cminvarvarhat=minvarvarhat=tctau=tcvarhat=dctau=dcvarhat=rep(NA, nsim)
taus=varhats=ctaus=cvarhats=varss=matrix(NA, nrow=nsim, ncol=nt)



####################################################################################
for(k in 1:nsim){
  small = smallbase
  large = largebase
  
  bigmerge = merge(small, large, by="bigid", all = TRUE)
  
  # ##########Different level of tau
  # bigmerge$w=rbinom(nrow(bigmerge), 1, .5)
  # bigmerge$x1=rpois(nrow(bigmerge), (8-3*bigmerge$w))
  # bigmerge$x2 = rnorm(nrow(bigmerge), -bigmerge$w, 3)
  # truetau=10
  # bigmerge$y = 5*bigmerge$x1+3*bigmerge$x2+
  # truetau*bigmerge$w+rnorm(nrow(bigmerge), 5, 10)
  ###################################
  
  ##########Different level of sigma
  # bigmerge$w=rbinom(nrow(bigmerge), 1, .5)
  # bigmerge$x1=rpois(nrow(bigmerge), (8-3*bigmerge$w))
  # bigmerge$x2 = rnorm(nrow(bigmerge), -bigmerge$w, 3)
  # truetau=50
  # bigmerge$y = 5*bigmerge$x1+3*bigmerge$x2+
  #   truetau*bigmerge$w+rnorm(nrow(bigmerge), 5, 25)
  ###################################
  # 
  # ##########Non-linear response
  # bigmerge$w=rbinom(nrow(bigmerge), 1, .5)
  # bigmerge$x1=rpois(nrow(bigmerge), (8-3*bigmerge$w))
  # bigmerge$x2 = rnorm(nrow(bigmerge), -bigmerge$w, 3)
  # truetau=50
  # bigmerge$y = 5*(.2*bigmerge$x1)^2+exp(.7*bigmerge$x2)+
  #   truetau*bigmerge$w+rnorm(nrow(bigmerge), 5, 10)
  # ###################################
  # 
  ##########Diff levels of association
  # bigmerge$w=rbinom(nrow(bigmerge), 1, .5)
  # #bigmerge$x1=rpois(nrow(bigmerge), (8-4*bigmerge$w))
  # bigmerge$x1=rpois(nrow(bigmerge), (8-3*bigmerge$w))
  # bigmerge$x2 = rnorm(nrow(bigmerge), -bigmerge$w, 3)
  # truetau=50
  # bigmerge$y = 15*bigmerge$x1-7*bigmerge$x2+
  #   truetau*bigmerge$w+rnorm(nrow(bigmerge), 5, 10)
  ###################################
  #
  
  ##########Original
  bigmerge$w=rbinom(nrow(bigmerge), 1, .5)
  bigmerge$x1=rpois(nrow(bigmerge), (8-3*bigmerge$w))
  bigmerge$x2 = rnorm(nrow(bigmerge), -bigmerge$w, 3)
  truetau=50
  bigmerge$y = 5*bigmerge$x1+3*bigmerge$x2+
    truetau*bigmerge$w+rnorm(nrow(bigmerge), 5, 10)
  ###################################
  
  ######################Non constant
  # bigmerge$w=rbinom(nrow(bigmerge), 1, .5)
  # bigmerge$x1=rpois(nrow(bigmerge), (8-3*bigmerge$w))
  # bigmerge$x2 = rnorm(nrow(bigmerge), -bigmerge$w, 3)
  # truetau=50+.5*bigmerge$x1 +6.5*bigmerge$x2
  # bigmerge$y = 5*bigmerge$x1+3*bigmerge$x2+
  #   truetau*bigmerge$w+rnorm(nrow(bigmerge), 5, 10)
  #######################################################3
  
  ##########################Even smaller tau
  # bigmerge$w=rbinom(nrow(bigmerge), 1, .5)
  # bigmerge$x1=rpois(nrow(bigmerge), (8-3*bigmerge$w))
  # bigmerge$x2 = rnorm(nrow(bigmerge), -bigmerge$w, 3)
  # truetau=1
  # bigmerge$y = 5*bigmerge$x1+3*bigmerge$x2+
  #   truetau*bigmerge$w+rnorm(nrow(bigmerge), 5, 10)
  #########################################################
  
  ################Outcome dependent on birth year#######
  # b_month = coalesce(bigmerge$bm.x, bigmerge$bm.y) - 6
  # b_month = ifelse(b_month < -6 | b_month > 6, 0, b_month)
  # b_base = b_month + 6
  # w_probs = rbeta(length(b_base), b_base, 6)
  # bigmerge$w=rbinom(nrow(bigmerge), 1, w_probs)
  # bigmerge$x1=rpois(nrow(bigmerge), (8-3*bigmerge$w))
  # bigmerge$x2 = rnorm(nrow(bigmerge), -bigmerge$w, 3)
  # truetau= 50
  # bigmerge$y = 5*bigmerge$x1+3*bigmerge$x2+ 5*bmonth +
  #   truetau*bigmerge$w+rnorm(nrow(bigmerge), 5, 10)
  ######################################################
  
  
  largedat=subset(bigmerge, is.na(freq.y)==F)[,c(1, 21)]
  smalldat=subset(bigmerge, is.na(freq.x)==F)[,c(1, 18, 19, 20)]
  
  large = (merge(large, largedat, by = 'bigid'))
  small = merge(small, smalldat, by = 'bigid')
  
  truemerge=merge(small, large, by="bigid")
  #names(truemerge)
  #doublecheck
  #print(summary(lm(y~w+x1+x2, data=truemerge)))
  
  
  #Analyze with linked data
  #Propensity score model
  pmod=glm(w~x1+x2, data=small, family="binomial")
  small$pscore=predict(pmod, type="response")
  
  #split into subclasses based on quantiles
  determ=merge(small, large, by=c("fname_c1", "lname_c1", "by", "bm", "bd"))
  quantiles=quantile(determ$pscore, probs=c(.3, .4, .5, .7))
  #quantiles=quantile(determ$pscore, probs=c(.38, .46, .6, .72))
  small$subclass=rep(NA, nrow(small))
  small$subclass[small$pscore<quantiles[1]]=1
  small$subclass[small$pscore>quantiles[1]&small$pscore<=quantiles[2]]=2
  small$subclass[small$pscore>quantiles[2]&small$pscore<=quantiles[3]]=3
  small$subclass[small$pscore>quantiles[3]&small$pscore<=quantiles[4]]=4
  small$subclass[small$pscore>quantiles[4]]=5

  
    #check balance, readjust if necessary
  determ=merge(small, large, by=c("fname_c1", "lname_c1", "by", "bm", "bd"))
 #table(determ$subclass,determ$w)
  if(min(table(determ$subclass,determ$w)) < 3){
    print("Balance check 1 failed")
    quantiles=quantile(determ$pscore, probs=c(.4, .52, .65, .72))
    small$subclass=rep(NA, nrow(small))
    small$subclass[small$pscore<quantiles[1]]=1
    small$subclass[small$pscore>quantiles[1]&small$pscore<=quantiles[2]]=2
    small$subclass[small$pscore>quantiles[2]&small$pscore<=quantiles[3]]=3
    small$subclass[small$pscore>quantiles[3]&small$pscore<=quantiles[4]]=4
    small$subclass[small$pscore>quantiles[4]]=5
  }
  if(min(table(determ$subclass,determ$w)) >= 3){
    print("Balance check 1 passed")}
  
  #repeat balance check
  determ=merge(small, large, by=c("fname_c1", "lname_c1", "by", "bm", "bd"))
  if(min(table(determ$subclass,determ$w)) < 3){
    print("Balance check 2 failed")
    quantiles=quantile(determ$pscore, probs=c(.3, .45, .52, .72))
    small$subclass=rep(NA, nrow(small))
    small$subclass[small$pscore<quantiles[1]]=1
    small$subclass[small$pscore>quantiles[1]&small$pscore<=quantiles[2]]=2
    small$subclass[small$pscore>quantiles[2]&small$pscore<=quantiles[3]]=3
    small$subclass[small$pscore>quantiles[3]&small$pscore<=quantiles[4]]=4
    small$subclass[small$pscore>quantiles[4]]=5
  }
  
  determ=merge(small, large, by=c("fname_c1", "lname_c1", "by", "bm", "bd"))
  if(min(table(determ$subclass,determ$w)) < 3){
    print("Balance check 3 failed")
    quantiles=quantile(determ$pscore, probs=c(.3, .42, .55, .72))
    small$subclass=rep(NA, nrow(small))
    small$subclass[small$pscore<quantiles[1]]=1
    small$subclass[small$pscore>quantiles[1]&small$pscore<=quantiles[2]]=2
    small$subclass[small$pscore>quantiles[2]&small$pscore<=quantiles[3]]=3
    small$subclass[small$pscore>quantiles[3]&small$pscore<=quantiles[4]]=4
    small$subclass[small$pscore>quantiles[4]]=5
  }
  
  
  determ=merge(small, large, by=c("fname_c1", "lname_c1", "by", "bm", "bd"))
  if(min(table(determ$subclass,determ$w)) < 3){
    print("Balance check 4 failed")
    quantiles=quantile(determ$pscore, probs=c(.28, .38, .48, .68))
    small$subclass=rep(NA, nrow(small))
    small$subclass[small$pscore<quantiles[1]]=1
    small$subclass[small$pscore>quantiles[1]&small$pscore<=quantiles[2]]=2
    small$subclass[small$pscore>quantiles[2]&small$pscore<=quantiles[3]]=3
    small$subclass[small$pscore>quantiles[3]&small$pscore<=quantiles[4]]=4
    small$subclass[small$pscore>quantiles[4]]=5
  }
  
  if(min(table(determ$subclass,determ$w)) < 3){
    print("Balance check 5 failed")}
  
  subclass = small$subclass
  
  ########################################## FS Linkage
  comparemat=merge(small, large, by="by")

  comparemat$monthagree=as.numeric(comparemat$bm.x==comparemat$bm.y)
  comparemat$dayagree=as.numeric(comparemat$bd.x==comparemat$bd.y)
  comparemat$truematch=as.numeric(comparemat$bigid.x==comparemat$bigid.y)
  comparemat$firstagree=jarowinkler(as.character(comparemat$fname_c1.x), as.character(comparemat$fname_c1.y))
  comparemat$lastagree=jarowinkler(as.character(comparemat$lname_c1.x), as.character(comparemat$lname_c1.y))
  #attach(comparemat)
  comparemat$simplefirst=as.numeric(comparemat$firstagree>.95)
  comparemat$simplelast=as.numeric(comparemat$lastagree>.95)
  #attach(comparemat)

  comp=comparemat

  m=rep(.95, 4)
  u=c(sum(comp$simplefirst)/length(comp$firstagree),
      sum(comp$simplelast/length(comp$lastagree)),
      sum(comp$monthagree)/length(comp$monthagree),
      sum(comp$dayagree)/length(comp$dayagree))


  logdiff=log(m[1]/u[1], 2)*comp$simplefirst+
    log(m[2]/u[2], 2)*comp$simplelast+
    log((1-m[1])/(1-u[1]), 2)*(1-comp$simplefirst)+
    log((1-m[2])/(1-u[2]), 2)*(1-comp$simplelast)+
    log(m[1]/u[3], 2)*comp$monthagree+
    log(m[2]/u[4], 2)*comp$dayagree+
    log((1-m[1])/(1-u[3]), 2)*(1-comp$monthagree)+
    log((1-m[2])/(1-u[4]), 2)*(1-comp$dayagree)

  # hist(logdiff[truematch==1])
  # hist(logdiff[truematch==0])
  #
  # length(logdiff[logdiff>0])

  comparemat$weigh=logdiff
  #attach(comparemat)

  ####These two lines are the only difference b/w the two simulations#################
  maxvals=aggregate(weigh ~ bigid.x, data = comparemat, max)
  comparemat=merge(comparemat, maxvals, by='bigid.x')
  ################################################################################
  #############################################################################
  
  #######################################JW LINKAGE####################
  # comparemat=merge(small, large, by="by")
  # 
  # #attach(comparemat)
  # 
  # comparemat$monthagree=jarowinkler(as.character(comparemat$bm.x), as.character(comparemat$bm.y))
  # comparemat$dayagree=jarowinkler(as.character(comparemat$bd.x), as.character(comparemat$bd.y))
  # comparemat$truematch=as.numeric(comparemat$bigid.x==comparemat$bigid.y)
  # comparemat$firstagree=jarowinkler(as.character(comparemat$fname_c1.x), as.character(comparemat$fname_c1.y))
  # comparemat$lastagree=jarowinkler(as.character(comparemat$lname_c1.x), as.character(comparemat$lname_c1.y))
  # 
  # 
  # comp=comparemat
  # 
  # 
  # comparemat$weigh=.25*comparemat$monthagree+.25*comparemat$dayagree+
  #   .25*comparemat$firstagree + .25*comparemat$lastagree
  # #attach(comparemat)
  # maxvals=aggregate(weigh ~ bigid.x, data = comparemat, max)
  # comparemat=merge(comparemat, maxvals, by='bigid.x')

  ############################################################################
  
  
  tau_j=ctau_j=lj=n0=n1=s0=s1=cvar=rep(NA, length(unique(subclass)))
  tau=varhat=ctau=cvarhat=vars = rep(NA, nt)
  for (t in 1:nt){
    dat=subset(comparemat, weigh.x>=threshold[t]&weigh.x==weigh.y)
    #attach(dat)
    for(s in 1:length(unique(subclass)))
    {
      sdat=subset(dat, subclass==unique(subclass)[s])
      tau_j[s]=mean(sdat$y[sdat$w==1])-mean(sdat$y[sdat$w==0])
      lj[s]=nrow(sdat)/nrow(dat)
      n1[s]=sum(sdat$w)
      n0[s]=nrow(sdat)-n1[s]
      s0[s]=var(sdat$y[sdat$w==0])
      s1[s]=var(sdat$y[sdat$w==1])
      #regression correction
      lmod=lm(y~x1+x2+w, data=sdat)
      ctau_j[s]=coef(lmod)[4]
      cvar[s]=(summary(lmod)$coefficients[4,2])^2
    }
    tau[t]=sum(tau_j*lj)
    varhat[t]=sum((s0/n0+s1/n1)*lj^2)
    ctau[t]=sum(ctau_j*lj)
    cvarhat[t]=sum(cvar*lj^2)
    vars[t]=var(dat$y[dat$w==0])/(length(dat$w[dat$w==0])/length(dat$w))+
      var(dat$y[dat$w==1])/(length(dat$w[dat$w==1])/length(dat$w))
  }
  
  truemerge=merge(small, large, by="bigid")
  dat=truemerge
  for(s in 1:length(unique(subclass)))
  {
    sdat=subset(truemerge, subclass==unique(subclass)[s])
    tau_j[s]=mean(sdat$y[sdat$w==1])-mean(sdat$y[sdat$w==0])
    lj[s]=nrow(sdat)/nrow(dat)
    n1[s]=sum(sdat$w)
    n0[s]=nrow(sdat)-n1[s]
    s0[s]=var(sdat$y[sdat$w==0])
    s1[s]=var(sdat$y[sdat$w==1])
    lmod=lm(y~x1+x2+w, data=truemerge)
    ctau_j[s]=coef(lmod)[4]
    cvar[s]=(summary(lmod)$coefficients[4,2])^2
  }
  ttau[k]=sum(tau_j*lj)
  tvarhat[k]=sum((s0/n0+s1/n1)*lj^2)
  tctau[k]=sum(ctau_j*lj)
  tcvarhat[k]=sum(cvar*lj^2)
  
  determ=merge(small, large, by=c("fname_c1", "lname_c1", "by", "bm", "bd"))
  dat=determ
  for(s in 1:length(unique(subclass)))
  {
    sdat=subset(dat, subclass==unique(subclass)[s])
    tau_j[s]=mean(sdat$y[sdat$w==1])-mean(sdat$y[sdat$w==0])
    lj[s]=nrow(sdat)/nrow(dat)
    n1[s]=sum(sdat$w)
    n0[s]=nrow(sdat)-n1[s]
    s0[s]=var(sdat$y[sdat$w==0])
    s1[s]=var(sdat$y[sdat$w==1])
    lmod=lm(y~x1+x2+w, data=sdat)
    ctau_j[s]=coef(lmod)[4]
    cvar[s]=(summary(lmod)$coefficients[4,2])^2
  }
  
 
  dtau[k]=sum(tau_j*lj)
  dvarhat[k]=sum((s0/n0+s1/n1)*lj^2)
  dctau[k]=sum(ctau_j*lj)
  dcvarhat[k]=sum(cvar*lj^2)
  
  
  taus[k,]=tau
  varhats[k,]=varhat
  ctaus[k,]=ctau
  cvarhats[k,]=cvarhat
  varss[k,]=vars
  print(k)
  minvartau[k]=(tau[varhat==min(varhat)])
  cminvartau[k]=(ctau[cvarhat==min(cvarhat)])
  minvarvarhat[k]=(varhat[varhat==min(varhat)])
  cminvarvarhat[k]=(cvarhat[cvarhat==min(cvarhat)])
  ybar1[k]=mean(truemerge$y[truemerge$w==1])
  ybar0[k]=mean(truemerge$y[truemerge$w==0])
  
  
  
  
  print(c(round(k, 0), minvartau[k], threshold[cminvartau[k]==ctaus[k,]]))

 # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
#
####################################################################


chosenthreshold = rep(NA, nsim)
for(k in 1:nsim){
  chosenthreshold[k]=threshold[minvartau[k]==taus[k,]]
}
table(chosenthreshold)
mean(chosenthreshold[is.na(chosenthreshold)==F])
median(chosenthreshold)


chosencthreshold = rep(NA, nsim)
for(k in 1:nsim){
  chosencthreshold[k]=threshold[cminvartau[k]==ctaus[k,]]
}
table(chosencthreshold)
mean(chosenthreshold[is.na(chosenthreshold)==F])
median(chosenthreshold)
mean(chosencthreshold)

#dtau=dtau[is.na(dtau)==F]
#dvarhat=dvarhat[is.na(dvarhat)==F]


nt=length(threshold)
matchrate=units=duplicates=rep(NA, nt)
exampledat = NULL
for (t in 1:nt){
  dat=subset(comparemat, weigh.x>=threshold[t]&weigh.x==weigh.y)
  matchrate[t]=sum(dat$bigid.x==dat$bigid.y)/nrow(dat)
  units[t]=nrow(dat)
  duplicates[t]=nrow(dat) - length(unique(dat$bigid.x))
  newdat = subset(comparemat, weigh.x>=threshold[t]&weigh.x<threshold[t+1]&weigh.x==weigh.y)
  
  exampledat = rbind(exampledat, as.vector(newdat[sample(nrow(newdat), 1), ]))
  
}


localminthreshold = localmintauhat = localminvarhat = 
  cutminthreshold = cutmintauhat = cutminvarhat = 
  cutmincthreshold = cutminctauhat = cutmincvarhat = 
  etsrthreshold = etsrtauhat = etsrvarhat= 
  etsrthreshold2 = etsrtauhat2 = etsrvarhat2= 
  etsrthreshold1 = etsrtauhat1 = etsrvarhat1= 
  etsrthresholdpoint5 = etsrtauhatpoint5 = etsrvarhatpoint5= 
  etsrcthresholdpoint5 = etsrctauhatpoint5 = etsrcvarhatpoint5= 
  etsrthreshold3 = etsrtauhat3 = etsrvarhat3= 
  simpletau = simplethreshold = simplevarhat = 
  simplectau = simplecthreshold = simplecvarhat = 
  localmincthreshold = localminctauhat = localmincvarhat = rep(NA, nsim)
for (k in 1:nsim){
  simpletau[k]=(taus[k,][varss[k,]==min(varss[k,])])
  simplevarhat[k]=(varhats[k,][varss[k,]==min(varss[k,])])
  simplethreshold[k]=(threshold[varss[k,]==min(varss[k,])])
  
  simplectau[k]=(ctaus[k,][varss[k,]==min(varss[k,])])
  simplecvarhat[k]=(cvarhats[k,][varss[k,]==min(varss[k,])])
  simplecthreshold[k]=(threshold[varss[k,]==min(varss[k,])])
  
  ismin=as.numeric(varhats[k,2:nt]<varhats[k,1:nt-1])
  
  if(sum(ismin)==0){localminthreshold[k] = (threshold[1])}
  if(sum(ismin)>0){localminthreshold[k] =(max(threshold[2:nt][ismin==1]))}
  
  if(sum(ismin)==0){localmintauhat[k] = (taus[k,1])}
  if(sum(ismin)>0){localmintauhat[k] =((taus[k, 2:nt][threshold[2:nt]==localminthreshold[k]]))}
  
  if(sum(ismin)==0){localminvarhat[k] = (varhats[k,1])}
  if(sum(ismin)>0){localminvarhat[k] =((varhats[k, 2:nt][threshold[2:nt]==localminthreshold[k]]))}
  
  ismin=as.numeric(cvarhats[k,2:nt]<cvarhats[k,1:nt-1])
  if(sum(ismin)==0){localmincthreshold[k] = (threshold[1])}
  if(sum(ismin)>0){localmincthreshold[k] =(max(threshold[2:nt][ismin==1]))}
  
  if(sum(ismin)==0){localminctauhat[k] = (ctaus[k,1])}
  if(sum(ismin)>0){localminctauhat[k] =((ctaus[k, 2:nt][threshold[2:nt]==localmincthreshold[k]]))}
  
  if(sum(ismin)==0){localmincvarhat[k] = (cvarhats[k,1])}
  if(sum(ismin)>0){localmincvarhat[k] =((cvarhats[k, 2:nt][threshold[2:nt]==localmincthreshold[k]]))}
  
  cutmintauhat[k] = taus[k,][varhats[k,]==min(varhats[k,][taus[k,]>(taus[k,nt]-.1*sqrt(varhats[k,nt]))&
                                                            taus[k,]<(taus[k,nt]+.1*sqrt(varhats[k,nt]))])] 
  cutminthreshold[k] = threshold[varhats[k,]==min(varhats[k,][taus[k,]>(taus[k,nt]-.1*sqrt(varhats[k,nt]))&
                                                                taus[k,]<(taus[k,nt]+.1*sqrt(varhats[k,nt]))])] 
  cutminvarhat[k] = varhats[k,][varhats[k,]==min(varhats[k,][taus[k,]>(taus[k,nt]-.1*sqrt(varhats[k,nt]))&
                                                               taus[k,]<(taus[k,nt]+.1*sqrt(varhats[k,nt]))])] 
  #old version-compares each to h0
  #contributedtaus=c((taus[k,2:nt-1]*units[2:nt-1]-taus[k,nt]*units[nt])/(units[2:nt-1]-units[nt]), taus[k,nt])
  
  #new version- compares each to h-1
  contributedtaus=c((taus[k,1:(nt-1)]*units[1:(nt-1)]-taus[k,2:nt]*units[2:nt])/(units[1:(nt-1)]-units[2:nt]), taus[k,nt])
  
  # 
  # etsrtauhat[k] = taus[k,][varhats[k,]==min(varhats[k,][contributedtaus>taus[k,nt]-1.96*sqrt(varhats[k,nt])&
  #                                                         contributedtaus<taus[k,nt]+1.96*sqrt(varhats[k,nt])])]
  # etsrthreshold[k] = threshold[varhats[k,]==min(varhats[k,][contributedtaus>taus[k,nt]-1.96*sqrt(varhats[k,nt])&
  #                                                         contributedtaus<taus[k,nt]+1.96*sqrt(varhats[k,nt])])]
  # etsrvarhat[k] = varhats[k,][varhats[k,]==min(varhats[k,][contributedtaus>taus[k,nt]-1.96*sqrt(varhats[k,nt])&
  #                                                         contributedtaus<taus[k,nt]+1.96*sqrt(varhats[k,nt])])]

  etsrtauhat2[k] = taus[k,][varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-2*sqrt(varhats[k,nt])&
                                                           taus[k,]<taus[k,nt]+2*sqrt(varhats[k,nt])])]
  etsrthreshold2[k] = threshold[varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-2*sqrt(varhats[k,nt])&
                                                               taus[k,]<taus[k,nt]+2*sqrt(varhats[k,nt])])]
  etsrvarhat2[k] = varhats[k,][varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-2*sqrt(varhats[k,nt])&
                                                              taus[k,]<taus[k,nt]+2*sqrt(varhats[k,nt])])]
  
  etsrtauhat1[k] = taus[k,][varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-1*sqrt(varhats[k,nt])&
                                                           taus[k,]<taus[k,nt]+1*sqrt(varhats[k,nt])])]
  etsrthreshold1[k] = threshold[varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-1*sqrt(varhats[k,nt])&
                                                               taus[k,]<taus[k,nt]+1*sqrt(varhats[k,nt])])]
  etsrvarhat1[k] = varhats[k,][varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-1*sqrt(varhats[k,nt])&
                                                              taus[k,]<taus[k,nt]+1*sqrt(varhats[k,nt])])]
  
  etsrtauhatpoint5[k] = taus[k,][varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-.5*sqrt(varhats[k,nt])&
                                                                taus[k,]<taus[k,nt]+.5*sqrt(varhats[k,nt])])]
  etsrthresholdpoint5[k] = threshold[varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-.5*sqrt(varhats[k,nt])&
                                                                    taus[k,]<taus[k,nt]+.5*sqrt(varhats[k,nt])])]
  etsrvarhatpoint5[k] = varhats[k,][varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-.5*sqrt(varhats[k,nt])&
                                                                   taus[k,]<taus[k,nt]+.5*sqrt(varhats[k,nt])])]
  
  etsrtauhat3[k] = taus[k,][varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-3*sqrt(varhats[k,nt])&
                                                           taus[k,]<taus[k,nt]+3*sqrt(varhats[k,nt])])]
  etsrthreshold3[k] = threshold[varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-3*sqrt(varhats[k,nt])&
                                                               taus[k,]<taus[k,nt]+3*sqrt(varhats[k,nt])])]
  etsrvarhat3[k] = varhats[k,][varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-3*sqrt(varhats[k,nt])&
                                                              taus[k,]<taus[k,nt]+3*sqrt(varhats[k,nt])])]
  
  etsrtauhat[k] = taus[k,][varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-.1*sqrt(varhats[k,nt])&
                                                          taus[k,]<taus[k,nt]+.1*sqrt(varhats[k,nt])])]
  etsrthreshold[k] = threshold[varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-.1*sqrt(varhats[k,nt])&
                                                              taus[k,]<taus[k,nt]+.1*sqrt(varhats[k,nt])])]
  etsrvarhat[k] = varhats[k,][varhats[k,]==min(varhats[k,][taus[k,]>taus[k,nt]-.1*sqrt(varhats[k,nt])&
                                                             taus[k,]<taus[k,nt]+.1*sqrt(varhats[k,nt])])]
  

  
  cutminctauhat[k] = ctaus[k,][cvarhats[k,]==min(cvarhats[k,][ctaus[k,]>(ctaus[k,nt]-sqrt(cvarhats[k,nt]))&
                                                                ctaus[k,]<(ctaus[k,nt]+sqrt(cvarhats[k,nt]))])] 
  cutmincthreshold[k] = threshold[cvarhats[k,]==min(cvarhats[k,][ctaus[k,]>(ctaus[k,nt]-sqrt(cvarhats[k,nt]))&
                                                                   ctaus[k,]<(ctaus[k,nt]+sqrt(cvarhats[k,nt]))])] 
  cutmincvarhat[k] = cvarhats[k,][cvarhats[k,]==min(cvarhats[k,][ctaus[k,]>(ctaus[k,nt]-sqrt(cvarhats[k,nt]))&
                                                                   ctaus[k,]<(ctaus[k,nt]+sqrt(cvarhats[k,nt]))])]
  etsrctauhatpoint5[k] = ctaus[k,][cvarhats[k,]==min(cvarhats[k,][ctaus[k,]>ctaus[k,nt]-.5*sqrt(cvarhats[k,nt])&
                                                                ctaus[k,]<ctaus[k,nt]+.5*sqrt(cvarhats[k,nt])])]
  etsrcthresholdpoint5[k] = threshold[cvarhats[k,]==min(cvarhats[k,][ctaus[k,]>ctaus[k,nt]-.5*sqrt(cvarhats[k,nt])&
                                                                    ctaus[k,]<ctaus[k,nt]+.5*sqrt(cvarhats[k,nt])])]
  etsrcvarhatpoint5[k] = cvarhats[k,][cvarhats[k,]==min(cvarhats[k,][ctaus[k,]>ctaus[k,nt]-.5*sqrt(cvarhats[k,nt])&
                                                                   ctaus[k,]<ctaus[k,nt]+.5*sqrt(cvarhats[k,nt])])]
  
  
}

# ETSR for tau = 1
# ETSRpoint5 = c(mean(etsrtauhatpoint5), var(etsrtauhatpoint5), mean(etsrvarhatpoint5), mean((etsrtauhatpoint5-1)^2), mean(etsrthresholdpoint5))
# ETSR1 = c(mean(etsrtauhat1), var(etsrtauhat1), mean(etsrvarhat1), mean((etsrtauhat1-1)^2), mean(etsrthreshold1))
# ETSR = c(mean(etsrtauhat), var(etsrtauhat), mean(etsrvarhat), mean((etsrtauhat-1)^2), mean(etsrthreshold))
# ETSR2 = c(mean(etsrtauhat2), var(etsrtauhat2), mean(etsrvarhat2), mean((etsrtauhat2-1)^2), mean(etsrthreshold2))
# ETSR3 = c(mean(etsrtauhat3), var(etsrtauhat3), mean(etsrvarhat3), mean((etsrtauhat3-1)^2), mean(etsrthreshold3))
# 
# xtable(print(rbind(ETSR, ETSRpoint5, ETSR1,
#                    ETSR2, ETSR3)), digits=1)

#ETSR for tau = 10
# ETSRpoint5 = c(mean(etsrtauhatpoint5), var(etsrtauhatpoint5), mean(etsrvarhatpoint5), mean((etsrtauhatpoint5-10)^2), mean(etsrthresholdpoint5))
# ETSR1 = c(mean(etsrtauhat1), var(etsrtauhat1), mean(etsrvarhat1), mean((etsrtauhat1-10)^2), mean(etsrthreshold1))
# ETSR = c(mean(etsrtauhat), var(etsrtauhat), mean(etsrvarhat), mean((etsrtauhat-10)^2), mean(etsrthreshold))
# ETSR2 = c(mean(etsrtauhat2), var(etsrtauhat2), mean(etsrvarhat2), mean((etsrtauhat2-10)^2), mean(etsrthreshold2))
# ETSR3 = c(mean(etsrtauhat3), var(etsrtauhat3), mean(etsrvarhat3), mean((etsrtauhat3-10)^2), mean(etsrthreshold3))
# 
# xtable(print(rbind(ETSR, ETSRpoint5, ETSR1,
#                    ETSR2, ETSR3)), digits=1)

# ETSR for tau = 50
ETSRpoint5 = c(mean(etsrtauhatpoint5), var(etsrtauhatpoint5), mean(etsrvarhatpoint5), mean((etsrtauhatpoint5-50)^2), mean(etsrthresholdpoint5))
ETSR1 = c(mean(etsrtauhat1), var(etsrtauhat1), mean(etsrvarhat1), mean((etsrtauhat1-50)^2), mean(etsrthreshold1))
ETSR = c(mean(etsrtauhat), var(etsrtauhat), mean(etsrvarhat), mean((etsrtauhat-50)^2), mean(etsrthreshold))
ETSR2 = c(mean(etsrtauhat2), var(etsrtauhat2), mean(etsrvarhat2), mean((etsrtauhat2-50)^2), mean(etsrthreshold2))
ETSR3 = c(mean(etsrtauhat3), var(etsrtauhat3), mean(etsrvarhat3), mean((etsrtauhat3-50)^2), mean(etsrthreshold3))

xtable(print(rbind(ETSR, ETSRpoint5, ETSR1, 
                    ETSR2, ETSR3)), digits=1)


par(mfrow=c(1,1))
plot(tau~threshold, type='l')


par(mfrow=c(2, 1))
colnames(taus)=round(threshold, digits=3)
boxplot(taus[, rev(seq_len(ncol(taus)))],
        ylab="Estimated treatment effect",
        xlab="Threshold cutoff")
abline(h=mean(ttau), col="red")

colnames(varhats)=round(threshold, digits=3)
meanvarhats=apply(varhats[, rev(seq_len(ncol(varhats)))], 2, mean)
names(meanvarhats)=round(sort(threshold, decreasing = TRUE), digits = 3)

plot(meanvarhats, type='l',
     ylab="Estimated variance",
     xlab="Threshold cutoff",
     xaxt='n')
axis(1, at=1:(nt), labels=names(meanvarhats))


##numeric summaries

sum(cutminthreshold == localminthreshold)


#################tau = 1
# Known=c(mean(dtau), var(dtau), mean(dvarhat), mean((dtau-1)^2))
# MEV=c(mean(minvartau), var(minvartau), mean(minvarvarhat), mean((minvartau-1)^2), mean(chosenthreshold))
# Perfect=c(mean(ttau), var(ttau), mean(tvarhat), mean((ttau-1)^2))
# LMEV = c(mean(localmintauhat), var(localmintauhat), mean(localminvarhat), mean((localmintauhat-1)^2))
# MEDOV = c(mean(simpletau), var(simpletau), mean(simplevarhat), mean((simpletau-1)^2), mean(simplethreshold))
# #ETSR = c(mean(cutmintauhat), var(cutmintauhat), mean(cutminvarhat), mean((cutmintauhat-1)^2))
# ETSR = c(mean(etsrtauhatpoint5), var(etsrtauhatpoint5), mean(etsrvarhatpoint5), mean((etsrtauhatpoint5-1)^2), mean(etsrthresholdpoint5))

#################tau = 10
# Known=c(mean(dtau), var(dtau), mean(dvarhat), mean((dtau-10)^2))
# MEV=c(mean(minvartau), var(minvartau), mean(minvarvarhat), mean((minvartau-10)^2), mean(chosenthreshold))
# Perfect=c(mean(ttau), var(ttau), mean(tvarhat), mean((ttau-10)^2))
# LMEV = c(mean(localmintauhat), var(localmintauhat), mean(localminvarhat), mean((localmintauhat-10)^2))
# MEDOV = c(mean(simpletau), var(simpletau), mean(simplevarhat), mean((simpletau-10)^2), mean(simplethreshold))
# #ETSR = c(mean(cutmintauhat), var(cutmintauhat), mean(cutminvarhat), mean((cutmintauhat-1)^2))
# ETSR = c(mean(etsrtauhatpoint5), var(etsrtauhatpoint5), mean(etsrvarhatpoint5), mean((etsrtauhatpoint5-10)^2), mean(etsrthresholdpoint5))


######## tau = 50
Known=c(mean(dtau), var(dtau), mean(dvarhat), mean((dtau-50)^2))
MEV=c(mean(minvartau), var(minvartau), mean(minvarvarhat), mean((minvartau-50)^2), mean(chosenthreshold))
Perfect=c(mean(ttau), var(ttau), mean(tvarhat), mean((ttau-50)^2))
LMEV = c(mean(localmintauhat), var(localmintauhat), mean(localminvarhat), mean((localmintauhat-50)^2))
MEDOV = c(mean(simpletau), var(simpletau), mean(simplevarhat), mean((simpletau-50)^2), mean(simplethreshold))
#ETSR = c(mean(cutmintauhat), var(cutmintauhat), mean(cutminvarhat), mean((cutmintauhat-1)^2))
ETSR = c(mean(etsrtauhatpoint5), var(etsrtauhatpoint5), mean(etsrvarhatpoint5), mean((etsrtauhatpoint5-50)^2), mean(etsrthresholdpoint5))

xtable(print(rbind(Perfect, Known, MEV, 
                   ETSR, MEDOV)), digits=1)


###Combined with regression
par(mfrow=c(2, 1))
colnames(ctaus)=round(threshold, digits=1)
boxplot(ctaus[, rev(seq_len(ncol(ctaus)))],
        ylab="Estimated treatment effect",
        xlab="Threshold cutoff")
abline(h=mean(truetau), col="red")

colnames(cvarhats)=round(threshold, digits=1)
meancvarhats=apply(cvarhats[, rev(seq_len(ncol(cvarhats)))], 2, mean)
names(meancvarhats)=round(sort(threshold, decreasing = TRUE), digits = 1)

plot(meancvarhats, type='l',
     ylab="Estimated variance",
     xlab="Threshold cutoff",
     xaxt='n')
axis(1, at=1:(nt), labels=names(meancvarhats))


chosencthreshold = rep(NA, nsim)
for(k in 1:nsim){
  chosencthreshold[k]=threshold[cminvartau[k]==ctaus[k,]]
}
table(chosenthreshold)


Known=c(mean(dctau), var(dctau), mean(dcvarhat), mean((dctau-50)^2))
MEV=c(mean(cminvartau), var(cminvartau), mean(cminvarvarhat), mean((cminvartau-50)^2), mean(chosencthreshold))
Perfect=c(mean(tctau), var(tctau), mean(tcvarhat), mean((tctau-50)^2))
LMEV = c(mean(localminctauhat), var(localminctauhat), mean(localmincvarhat), mean((localminctauhat-50)^2))
MEDOV = c(mean(simplectau), var(simplectau), mean(simplecvarhat), mean((simplectau-50)^2), mean(simplecthreshold))
ETSR = c(mean(etsrctauhatpoint5), var(etsrctauhatpoint5), mean(etsrcvarhatpoint5), mean((etsrctauhatpoint5-50)^2), mean(etsrcthresholdpoint5))

xtable(print(rbind(Perfect, Known, MEV, 
                   ETSR, MEDOV)), digits=1)



#####tau = 10
# Known=c(mean(dctau), var(dctau), mean(dcvarhat), mean((dctau-10)^2))
# MEV=c(mean(cminvartau), var(cminvartau), mean(cminvarvarhat), mean((cminvartau-10)^2), mean(chosencthreshold))
# Perfect=c(mean(tctau), var(tctau), mean(tcvarhat), mean((tctau-10)^2))
# LMEV = c(mean(localminctauhat), var(localminctauhat), mean(localmincvarhat), mean((localminctauhat-10)^2))
# MEDOV = c(mean(simplectau), var(simplectau), mean(simplecvarhat), mean((simplectau-10)^2), mean(simplecthreshold))
# ETSR = c(mean(etsrctauhatpoint5), var(etsrctauhatpoint5), mean(etsrcvarhatpoint5), mean((etsrctauhatpoint5-10)^2), mean(etsrcthresholdpoint5))
# 
# xtable(print(rbind(Perfect, Known, MEV, 
#                    ETSR, MEDOV)), digits=1)
########tau = 1
# Known=c(mean(dctau), var(dctau), mean(dcvarhat), mean((dctau-1)^2))
# MEV=c(mean(cminvartau), var(cminvartau), mean(cminvarvarhat), mean((cminvartau-1)^2), mean(chosencthreshold))
# Perfect=c(mean(tctau), var(tctau), mean(tcvarhat), mean((tctau-1)^2))
# LMEV = c(mean(localminctauhat), var(localminctauhat), mean(localmincvarhat), mean((localminctauhat-1)^2))
# MEDOV = c(mean(simplectau), var(simplectau), mean(simplecvarhat), mean((simplectau-1)^2), mean(simplecthreshold))
# ETSR = c(mean(etsrctauhatpoint5), var(etsrctauhatpoint5), mean(etsrcvarhatpoint5), mean((etsrctauhatpoint5-1)^2), mean(etsrcthresholdpoint5))
# 
# xtable(print(rbind(Perfect, Known, MEV, 
#                    ETSR, MEDOV)), digits=1)      


threshold = threshold[1:nt]
nt=length(threshold)
matchrate=units=duplicates=rep(NA, nt)
exampledat = NULL
for (t in 1:nt){
  dat=subset(comparemat, weigh.x>=threshold[t]&weigh.x==weigh.y)
  matchrate[t]=sum(dat$bigid.x==dat$bigid.y)/nrow(dat)
  units[t]=nrow(dat)
  duplicates[t]=nrow(dat) - length(unique(dat$bigid.x))
  newdat = subset(comparemat, weigh.x>=threshold[t]&weigh.x<threshold[t+1]&weigh.x==weigh.y)

  exampledat = rbind(exampledat, as.vector(newdat[sample(nrow(newdat), 1), ]))
  
}
newdat = subset(comparemat, weigh.x>=threshold[nt]&weigh.x==weigh.y)
exampledat = rbind(exampledat, as.vector(newdat[sample(nrow(newdat), 1), ]))

datasetsummary = data.frame(cbind(Threshold=round(threshold, 1), 
                                  MatchRate=paste(round(matchrate*100, 1), '%'), 
                                  Units =units, Duplicates = duplicates))
exampledata = (cbind(Threshold=round(threshold, 2), exampledat[,c(4, 17, 6, 19, 8, 9, 2, 21, 22, 2, 1, 15)]))
exampledata = cbind(exampledata[,1:11], ifelse(exampledata[,12]==exampledata[,13], "True", "False"))
colnames(exampledata)=c("Threshold", "F. Name:A", "F. Name:B", "L. Name:A", "L. Name:B",
                        "Mo:A", "D:A", "Y:A", "Mo:B", "D:B", "Y:B", "Status")
print(exampledata)
print(xtable(exampledata), include.rownames = F)


print(xtable(datasetsummary), include.rownames = F)
