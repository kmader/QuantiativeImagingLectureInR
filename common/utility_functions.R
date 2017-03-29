# ---- libraries ----
require(ggplot2) #for plots
require(lattice) # nicer scatter plots
require(plyr) # for processing data.frames
library(dplyr)
require(grid) # contains the arrow function
require(jpeg)
#require(doMC) # for parallel code
require(png) # for reading png images
require(gridExtra)

require(reshape2) # for the melt function

## To install EBImage

if (!require("EBImage")) { # for more image processing
  source("http://bioconductor.org/biocLite.R")
  biocLite("EBImage")
}

used.libraries<-c("ggplot2","lattice","plyr","reshape2","grid","gridExtra","biOps","png","EBImage")

# ---- common-functions ----
# start parallel environment
# registerDoMC()
# functions for converting images back and forth
im.to.df<-function(in.img,out.col="val") {
  out.im<-expand.grid(x=1:nrow(in.img),y=1:ncol(in.img))
  out.im[,out.col]<-as.vector(in.img)
  out.im
}
df.to.im<-function(in.df,val.col="val",inv=F) {
  in.vals<-in.df[[val.col]]
  if(class(in.vals[1])=="logical") in.vals<-as.integer(in.vals*255)
  if(inv) in.vals<-255-in.vals
  out.mat<-matrix(in.vals,nrow=length(unique(in.df$x)),byrow=F)
  attr(out.mat,"type")<-"grey"
  out.mat
}
ddply.cutcols<-function(...,cols=1) {
  # run standard ddply command 
  cur.table<-ddply(...)
  cutlabel.fixer<-function(oVal) {
    sapply(oVal,function(x) {
      cnv<-as.character(x)
      mean(as.numeric(strsplit(substr(cnv,2,nchar(cnv)-1),",")[[1]]))
    })
  }
  cutname.fixer<-function(c.str) {
    s.str<-strsplit(c.str,"(",fixed=T)[[1]]
    t.str<-strsplit(paste(s.str[c(2:length(s.str))],collapse="("),",")[[1]]
    paste(t.str[c(1:length(t.str)-1)],collapse=",")
  }
  for(i in c(1:cols)) {
    cur.table[,i]<-cutlabel.fixer(cur.table[,i])
    names(cur.table)[i]<-cutname.fixer(names(cur.table)[i])
  }
  cur.table
}

# since biOps isn't available use a EBImage based version
imgSobel<-function(img) {
  kern<-makeBrush(3)*(-1)
  kern[2,2]<-8
  filter2(img,kern)
}


show.pngs.as.grid<-function(file.list,title.fun,zoom=1) {
  preparePng<-function(x) rasterGrob(readPNG(x,native=T,info=T),width=unit(zoom,"npc"),interp=F)
  labelPng<-function(x,title="junk") (qplot(1:300, 1:300, geom="blank",xlab=NULL,ylab=NULL,asp=1)+
                                        annotation_custom(preparePng(x))+
                                        labs(title=title)+theme_bw(24)+
                                        theme(axis.text.x = element_blank(),
                                              axis.text.y = element_blank()))
  imgList<-llply(file.list,function(x) labelPng(x,title.fun(x)) )
  do.call(grid.arrange,imgList)
}
## Standard image processing tools which I use for visualizing the examples in the script
commean.fun<-function(in.df) {
  ddply(in.df,.(val), function(c.cell) {
    weight.sum<-sum(c.cell$weight)
    data.frame(xv=mean(c.cell$x),
               yv=mean(c.cell$y),
               xm=with(c.cell,sum(x*weight)/weight.sum),
               ym=with(c.cell,sum(y*weight)/weight.sum)
    )
  })
}

colMeans.df<-function(x,...) as.data.frame(t(colMeans(x,...)))

pca.fun<-function(in.df) {
  ddply(in.df,.(val), function(c.cell) {
    c.cell.cov<-cov(c.cell[,c("x","y")])
    c.cell.eigen<-eigen(c.cell.cov)
    
    c.cell.mean<-colMeans.df(c.cell[,c("x","y")])
    out.df<-cbind(c.cell.mean,
                  data.frame(vx=c.cell.eigen$vectors[1,],
                             vy=c.cell.eigen$vectors[2,],
                             vw=sqrt(c.cell.eigen$values),
                             th.off=atan2(c.cell.eigen$vectors[2,],c.cell.eigen$vectors[1,]))
    )
  })
}
vec.to.ellipse<-function(pca.df) {
  ddply(pca.df,.(val),function(cur.pca) {
    # assume there are two vectors now
    create.ellipse.points(x.off=cur.pca[1,"x"],y.off=cur.pca[1,"y"],
                          b=sqrt(5)*cur.pca[1,"vw"],a=sqrt(5)*cur.pca[2,"vw"],
                          th.off=pi/2-atan2(cur.pca[1,"vy"],cur.pca[1,"vx"]),
                          x.cent=cur.pca[1,"x"],y.cent=cur.pca[1,"y"])
  })
}

# test function for ellipse generation
# ggplot(ldply(seq(-pi,pi,length.out=100),function(th) create.ellipse.points(a=1,b=2,th.off=th,th.val=th)),aes(x=x,y=y))+geom_path()+facet_wrap(~th.val)+coord_equal()
create.ellipse.points<-function(x.off=0,y.off=0,a=1,b=NULL,th.off=0,th.max=2*pi,pts=36,...) {
  if (is.null(b)) b<-a
  th<-seq(0,th.max,length.out=pts)
  data.frame(x=a*cos(th.off)*cos(th)+b*sin(th.off)*sin(th)+x.off,
             y=-1*a*sin(th.off)*cos(th)+b*cos(th.off)*sin(th)+y.off,
             id=as.factor(paste(x.off,y.off,a,b,th.off,pts,sep=":")),...)
}
deform.ellipse.draw<-function(c.box) {
  create.ellipse.points(x.off=c.box$x[1],
                        y.off=c.box$y[1],
                        a=c.box$a[1],
                        b=c.box$b[1],
                        th.off=c.box$th[1],
                        col=c.box$col[1])                    
}
bbox.fun<-function(in.df) {
  ddply(in.df,.(val), function(c.cell) {
    c.cell.mean<-colMeans.df(c.cell[,c("x","y")])
    xmn<-emin(c.cell$x)
    xmx<-emax(c.cell$x)
    ymn<-emin(c.cell$y)
    ymx<-emax(c.cell$y)
    out.df<-cbind(c.cell.mean,
                  data.frame(xi=c(xmn,xmn,xmx,xmx,xmn),
                             yi=c(ymn,ymx,ymx,ymn,ymn),
                             xw=xmx-xmn,
                             yw=ymx-ymn
                  ))
  })
}

# since the edge of the pixel is 0.5 away from the middle of the pixel
emin<-function(...) min(...)-0.5
emax<-function(...) max(...)+0.5
extents.fun<-function(in.df) {
  ddply(in.df,.(val), function(c.cell) {
    c.cell.mean<-colMeans.df(c.cell[,c("x","y")])
    out.df<-cbind(c.cell.mean,data.frame(xmin=c(c.cell.mean$x,emin(c.cell$x)),
                                         xmax=c(c.cell.mean$x,emax(c.cell$x)),
                                         ymin=c(emin(c.cell$y),c.cell.mean$y),
                                         ymax=c(emax(c.cell$y),c.cell.mean$y)))
  })
}

common.image.path<-"../common"
qbi.file<-function(file.name) file.path(common.image.path,"figures",file.name)
qbi.data<-function(file.name) file.path(common.image.path,"data",file.name)

th_fillmap.fn<-function(max.val) scale_fill_gradientn(colours=rainbow(10),limits=c(0,max.val))