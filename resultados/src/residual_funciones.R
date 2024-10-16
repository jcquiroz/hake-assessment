

ressd <- function(obs, pre, log = TRUE){
  obs[obs == 0] <- NA
  if (log){
    dife = log(obs) - log(pre)
  } else {
    dife = obs - pre
  }
  return(dife)
}

ssruns_sig3 <- function(x, type=NULL, mixing="less") {
  if(is.null(type)) type="resid"
  if(type=="resid"){
    mu = 0} else {mu = mean(x, na.rm = TRUE)}
  alternative=c("two.sided","left.sided")[which(c("two.sided", "less")%in%mixing)]
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
  if(nlevels(factor(sign(x)))>1){
    # Make the runs test non-parametric
    runstest = randtests::runs.test(x,threshold = 0,alternative = alternative)
    if(is.na(runstest$p.value)) p.value =0.001
    pvalue = round(runstest$p.value,3)} else {
      pvalue = 0.001
    }
  
  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}

#ss3rep=ss3sma
mixing="less"
subplots=c("cpue","len","age")[1]
plot=TRUE
print=FALSE
png=print
pdf=FALSE
indexselect = NULL
miny = 2 #######
pch=21
lty=1
lwd=2
tickEndYr=FALSE
xlim="default"
ylimAdj=1.4
xaxs="i"
yaxs="i"
xylabs=TRUE
type="o"
legend=TRUE
legendloc="top"
legendcex=1
pwidth=6.5
pheight=5.0
punits="in"
res=300
ptsize=10
cex.main=1
plotdir=NULL
filenameprefix=""
par=list(mar=c(5,4,1,1)+.1)
verbose=TRUE
new=TRUE
add=FALSE




SSplotRunstest <- function(ss3rep=ss3sma, mixing="less", subplots=c("cpue","len","age")[1],
                             plot=TRUE, print=FALSE, png=print, pdf=FALSE,
                             indexselect = NULL,
                             miny = 1,
                             pch=21, lty=1, lwd=2,
                             tickEndYr=FALSE,
                             xlim="default", ylimAdj=1.4,
                             xaxs="i", yaxs="i",
                             xylabs=TRUE,
                             type="o", 
                             legend=TRUE, legendloc="top",
                             legendcex=1,
                             pwidth=6.5,pheight=5.0,punits="in",res=300,ptsize=10,cex.main=1,
                             plotdir=NULL,
                             filenameprefix="",
                             par=list(mar=c(5,4,1,1)+.1),
                             verbose=TRUE,
                             new=TRUE,
                             add=FALSE){
  
    #-------------------------------------------
    # r4ss plotting functions and coding style
    #-------------------------------------------
    # subfunction to write png files
    if(!add) graphics.off()
    if(add){
      print=F
      png=F
    }
  
    subplots = subplots[2]
    #datatypes= c("Index","Mean length","Mean age")
    datatypes= c("Indice","Longitud Media","Edad Media")
    ylabel = datatypes[which(c("cpue","len","age")%in%subplots)]
    if(verbose) cat('\n',"Running Runs Test Diagnosics for",datatypes[which(c("cpue","len","age")%in%subplots)],'\n')
    if(subplots=="cpue"){
    #cpue = ss3rep$cpue
    #cpue$residuals = ifelse(is.na(cpue$Obs) | is.na(cpue$Like),NA,log(cpue$Obs)-log(cpue$Exp))
    #if(is.null(cpue$Fleet_name)){ # Deal with Version control
    #cpue$Fleet_name = cpue$Name}
    #Res = cpue
      Res = sout$pf_res ############ 
    }
    
    ### FIX for length or age residual
    if(subplots=="len" | subplots=="age"){
      #comps = SScompsTA1.8(ss3rep,fleet=NULL,type=subplots,plotit = FALSE)$runs_dat
      #comps$residuals = ifelse(is.na(comps$Obs),NA,log(comps$Obs)-log(comps$Exp))
      #if(is.null(comps$Fleet_name)){ # Deal with Version control
      #comps$Fleet_name = comps$Name}
      #Res = comps  
      Res = prop_m$pf_res ##############
    }  
      
    pngfun <- function(file){
    # if extra text requested, add it before extention in file name
    file <- paste0(filenameprefix, file)
    # open png file
    png(filename=file.path(plotdir,file),
        width=pwidth,height=pheight,units=punits,res=res,pointsize=ptsize)
    # change graphics parameters to input value
    par(par)
  }
  
  # subset if indexselect is specified
  if(is.null(indexselect) == F & is.numeric(indexselect)){
    iname =  unique(Res$Fleet_name)[indexselect]
    if(TRUE %in% is.na(iname)) stop("One or more index numbers exceed number of available indices")
    Res = Res[Res$Fleet_name%in%iname,]
  }

    # Define indices
    #indices = unique(Res$Fleet_name)
    indices = 'Talla Media' ####### Ajustar para nombres
    n.indices = length(indices)
    series = 1:n.indices
    
    
    
    if(png) print <- TRUE
    if(png & is.null(plotdir))
      stop("to print PNG files, you must supply a directory as 'plotdir'")
    
    # check for internal consistency
    if(pdf & png){
      stop("To use 'pdf', set 'print' or 'png' to FALSE.")
    }
    if(pdf){
      if(is.null(plotdir)){
        stop("to write to a PDF, you must supply a directory as 'plotdir'")
      }
      pdffile <- file.path(plotdir,
                           paste0(filenameprefix, "SSplotComparisons_",
                                  format(Sys.time(), '%d-%b-%Y_%H.%M' ), ".pdf"))
      pdf(file=pdffile, width=pwidth, height=pheight)
      if(verbose) cat("PDF file with plots will be:",pdffile,'\n')
      par(par)
    }
    
    #---------------------------------------
    ######################## 
    plot_runs <- function(resid){  
      resid = prop_ym
      labels=c("Year",             #1
               "Residuals",        #2
               "Log Indice")        #3
      
      # open new window if requested
      if(plot & png==FALSE){
        if(!add) dev.new(width=pwidth,height=pheight,pointsize=ptsize,record=TRUE)
        
      } else {
        
        if(!add) par(par)
      }
      
      # get quantities for plot
      # yr <- resid$Yr
      # ti <- resid$Time  
      
      yr <- resid$year   ################
      ti <- NULL
      #ylab= paste(ylabel,"residuals")
      ylab= paste(ylabel,"de residuales")
      miny = 2
      
      ### make plot of index fits
      
      # Do runs test
      # runstest = ssruns_sig3(x=as.numeric(resid$residuals),type="resid",mixing=mixing)
      runstest = ssruns_sig3(x=as.numeric(resid$res_c),type="resid",mixing=mixing)
      ################
      
      # if no values included in subset, then set ylim based on all values
      ylim=c(min(-miny,runstest$sig3lim[1]*ylimAdj),max(miny,runstest$sig3lim[2]*ylimAdj))
      
      if(xlim[1]=="default") xlim = c(floor(min(ti,yr)-.1),ceiling(max(ti,yr)+0.1))
      
        plot(0, type = "n", xlim = xlim, yaxs = yaxs, 
             ylim = ylim, xlab = ifelse(xylabs,"Año",""), ylab = ifelse(xylabs,ylab,""), axes = FALSE)
        
        
        lims = runstest$sig3lim
        cols =  c(rgb(1,0,0,0.5),rgb(0,1,0,0.5))[ifelse(runstest$p.runs<0.05,1,2)]
        #rect(min(resid$Yr-1),lims[1],max(resid$Yr+1),lims[2],col=cols,border=cols) # only show runs if RMSE >= 0.1
        rect(min(resid$year-1),lims[1],max(resid$year+1),lims[2],col=cols,border=cols)
        ################
        abline(h=0,lty=2)
        for(j in 1:length(resid$year)){
          # lines(c(resid$Time[j],resid$Time[j]),c(0,resid$residuals[j]))
          lines(c(resid$year[j],resid$year[j]),c(0,resid$res_c[j]))
        }
        points(resid$year,resid$res_c,pch=pch,bg=ifelse(resid$res_c < lims[1] | resid$res_c > lims[2],2,"white"),cex=1)
        
        if(legend){
        # legend(legendloc,paste(resid$Fleet_name[1]),bty="n",y.intersp = -0.2,cex=legendcex+0.1)
        legend(legendloc,paste('CRUCERO'),bty="n",y.intersp = -0.2,cex=legendcex+0.1) } ################
        
        #legend("topright", bty='n',
        #       c("Passed","Failed"),pch=15,col=c(rgb(0,1,0,0.5),rgb(1,0,0,0.5)),pt.cex=2,cex=legendcex,y.intersp=0.9 )
        
        axis(1, at=resid$year)
        if(tickEndYr) axis(1, at=max(resid$year))
        
        axis(2)
        box()
        
      return(runstest)
    } # End of plot_runs function  
    #------------------------------------------------------------
    
    
    if(verbose) cat("Plotting Residual Runs Tests \n")
    if(plot){ 
      # LOOP through fleets
      nfleets=n.indices
      if(print){
        
        runs = NULL
        for(fi in 1:nfleets){
          resid = Res[Res$Fleet_name==indices[fi],]
          pngfun(paste0("residruns_",indices[fi],".png",sep=""))
          par(par)
          if(nrow(resid)>3 & (max(resid$Time)-min(resid$Time))>3){
          get_runs = plot_runs(resid)    
          dev.off()
          runs = rbind(runs,c(get_runs$p.runs,get_runs$sig3lim))
          } else {
            runs = rbind(runs,c(NA,NA,NA))}
          
        } # End of Fleet Loop
      }
      
      
      runs = NULL
      for(fi in 1:nfleets){
        resid = Res[Res$Fleet_name==indices[fi],]
        #if(nrow(resid)>3 & (max(resid$Time)-min(resid$Time))>3){ 
        if(nrow(resid)>3 & (max(resid$Time)-min(resid$Time))>3){ #############
        if(!add)(par)
        get_runs = plot_runs(resid)    
        runs = rbind(runs,c(get_runs$p.runs,get_runs$sig3lim))
        # End of Fleet Loop
        } else {
        runs = rbind(runs,c(NA,NA,NA))
        }
    }   
    }
    
    runstable = data.frame(Index=indices,runs.p=as.matrix(runs)[,1],Test=ifelse(is.na(as.matrix(runs)[,1]),"Excluded",ifelse(as.matrix(runs)[,1]<0.05,"Failed","Passed")),sigma3.lo=as.matrix(runs)[,2],sigma3.hi=as.matrix(runs)[,3],type=subplots) 
    colnames(runstable) = c("Index","runs.p","test","sigma3.lo","sigma3.hi","type")
    if(verbose) cat(paste0("\n","Runs Test stats by ",datatypes[which(c("cpue","len","age")%in%subplots)],":","\n"))
    return(runstable)
} # end of SSplotRuns()
#-----------------------------------------------------------------------------------------


