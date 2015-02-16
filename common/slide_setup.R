require(knitr)
# default settings, # settings for presentation version
echo.val<-F
fig.height<-5
dpi<-80
cache<-T
fig.path<-"pres_figures/"
cache.path<-"pres_cache/"

if(exists("printed")) { # settings for printed version (if the variable exists)
  echo.val<-T # show code
  fig.height<-3
  dpi<-100
  cache<-T
  fig.path<-"print_figures/"
  cache.path<-"print_cache/"
}

opts_chunk$set(dpi=dpi,cache=cache,
               cache.path=cache.path,results='hide',
               warning=F,fig.align='center',echo=echo.val,
               fig.height=fig.height,fig.path=fig.path,message=F) #dev="CairoPNG"