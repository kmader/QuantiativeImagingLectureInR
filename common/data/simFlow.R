## Generate a nice random flowing system
library(plyr)
# The basic setup
# Edge creation criteria 

#' Generate frames 
#' @author Kevin Mader (kevin.mader@gmail.com)
#' Generates flow with given object count, frame count and randomness
#' the box and crop are introduced to allow for objects entering and 
#' leaving the field of view
#'
#' @param n.objects Object count
#' @param n.frames Number of frames
#' @param box.size Size of initial box in x,y,z where objects are randomly placed
#' @param crop.size Size of cropping to be performed at the last step
#' @param flow.field The vector function used to calculate the flow field at each point x,y,z
#'                    takes (x0,y0,z0) and returns (vx,vy,vz) which will then be scaled by the mean
#' @param base.rand Base randomness used for the default values of the flow.x,y,z.rand
#' @param rand.fun The function to use to generate random numbers
#' 
generate.frames<-function(n.objects=50,n.frames=10,box.size=0.5,crop.size=c(-.5,.5),
                          base.rand=0.05,flow.x.rand=NA,flow.y.rand=NA,flow.z.rand=NA,rand.fun=runif,
                          flow.field=NA,flow.rate=0.1,obj.creation=0,obj.destruction=0) {
  if(is.na(flow.x.rand)) flow.x.rand<-base.rand
  if(is.na(flow.y.rand)) flow.y.rand<-base.rand
  if(is.na(flow.z.rand)) flow.z.rand<-base.rand
  if(suppressWarnings(is.na(flow.field))) flow.field<-flow.linear(0,0,flow.rate)
  prf<-function(n=n.objects) rand.fun(n,-box.size,box.size)
  start.objs<-data.frame(LACUNA_NUMBER=c(1:n.objects),
                         REAL_LACUNA_NUMBER=c(1:n.objects),
                         POS_X=prf(),
                         POS_Y=prf(),
                         POS_Z=prf(),
                         sample=1)
  out.frames<-list(start.objs)
  
  for(c.frame in c(2:n.frames)) {
    last.frame<-out.frames[[c.frame-1]]
    next.frame<-data.frame(last.frame[,!(names(last.frame) %in% c("sample"))],sample=c.frame)
    rf<-function() rand.fun(n.objects,-1,1)
    flow.v<-flow.field(last.frame$POS_X,last.frame$POS_Y,last.frame$POS_Z)
    next.frame$POS_X<-last.frame$POS_X+flow.v$vx+flow.x.rand*rf()
    next.frame$POS_Y<-last.frame$POS_Y+flow.v$vy+flow.y.rand*rf()
    next.frame$POS_Z<-last.frame$POS_Z+flow.v$vz+flow.z.rand*rf()
    next.frame$LACUNA_NUMBER<-order(rf()) # scramble id
    next.frame$REAL_LACUNA_NUMBER<-start.objs$LACUNA_NUMBER # keep original ID
    out.frames[[c.frame]]<-next.frame
  }
  #crop objects
  out.frames<-llply(out.frames,function(c.objs) {
    subset(c.objs,(POS_X>crop.size[1] & POS_X<crop.size[2]) &
             (POS_Y>crop.size[1] & POS_Y<crop.size[2]) & 
             (POS_Z>crop.size[1] & POS_Z<crop.size[2]))
  })
  # random obj destruction
  if(obj.destruction>0) {
    out.frames<-llply(out.frames,function(c.objs) {
      c.objs[sample(1:nrow(c.objs),round((1-obj.destruction)*nrow(c.objs))),]
    })
  }
  if(obj.creation>0) {
    out.frames<-llply(out.frames,function(c.objs) {
      n.cnt<-round(obj.creation*n.objects)
      rbind(c.objs,
            data.frame(LACUNA_NUMBER=n.objects+c(1:n.cnt),
                       REAL_LACUNA_NUMBER=-1,
                       POS_X=prf(n.cnt),
                       POS_Y=prf(n.cnt),
                       POS_Z=prf(n.cnt),
                       sample=c.objs$sample[1]))
    })
  }
  out.frames
}
#' Generate edges 
#' @author Kevin Mader (kevin.mader@gmail.com)
#' Generates edges in a given list of frames by calculating all objects
#' that are sufficiently close to other objects
#'
#' @param edge.max.length the maximum distance in order to form an edge
#' 
generate.edges<-function(flow.frames,edge.max.length=0.3) {
  llply(flow.frames,function(c.objs) {
    keep.cols<-c("LACUNA_NUMBER","POS_X","POS_Y","POS_Z")
    cross.objs<-merge(c.objs[,keep.cols],c.objs[,keep.cols],by=c())
    
    cross.objs$dist<-with(cross.objs,sqrt((POS_X.x-POS_X.y)^2+(POS_Y.x-POS_Y.y)^2+(POS_Z.x-POS_Z.y)^2))
    # filter the edges so only unique edges under a certain length are kept
    cross.objs<-subset(cross.objs,(dist<edge.max.length) & (LACUNA_NUMBER.x<LACUNA_NUMBER.y)) 
    if(nrow(cross.objs)<1) c.samp<-c()
    else c.samp<-c.objs[1,"sample"]
    data.frame(Component.1=cross.objs$LACUNA_NUMBER.x,Component.2=cross.objs$LACUNA_NUMBER.y,Voxels=cross.objs$dist/edge.max.length,sample=c.samp)
  })
}
#' Generate flow 
#' @author Kevin Mader (kevin.mader@gmail.com)
#' Generates frames and edges
#' @param ... forwards most of the parameters directly to the generate.frame code
#' @param edge.max.length the maximum distance in order to form an edge
#' 
generate.flow<-function(...,edge.max.length=0.3) {
  cur.flow<-generate.frames(...)
  list(frames=cur.flow,edges=generate.edges(cur.flow,edge.max.length=edge.max.length))
}
#' Linear Flow Simulation
#' @author Kevin Mader (kevin.mader@gmail.com)
#' The field for a constant linear flow
#' @param flow.z.mean the mean velocity in z
flow.linear<-function(flow.x=0.1,flow.y=0.0,flow.z=0.0) {
  function(x,y,z) {data.frame(vx=flow.x,vy=flow.y,vz=flow.z)}
}

#' Shear Flow Simulation
#' @author Kevin Mader (kevin.mader@gmail.com)
#' The field for a shear flow
#' @param shear.slope the slope of the shear
flow.shear.z<-function(shear.slope=0.1) {
  function(x,y,z) {data.frame(vx=0,vy=0,vz=shear.slope*x)}
}
#' Parallel Reversed Flows
#' @author Kevin Mader (kevin.mader@gmail.com)
#' two parallel linear flows in opposite directions
#' @param flow.z.mean the mean velocity in z
flow.parallel<-function(flow.z=0.1) {
  function(x,y,z) {data.frame(vx=0,vy=0,vz=ifelse(x>0,flow.z,-flow.z))}
}

#' Explosion / Implosion
#' @author Kevin Mader (kevin.mader@gmail.com)
#' flow is away from the center
#' @param exp.rate is the rate of explosion
flow.explosion<-function(exp.rate=0.1) {
  function(x,y,z) {
    r<-sqrt(x^2+y^2+z^2) # calculate unit vector then scale it
    data.frame(vx=x/r*exp.rate,vy=y/r*exp.rate,vz=z/r*exp.rate)
  }
}

#' Twist
#' @author Kevin Mader (kevin.mader@gmail.com)
#' radian step size
#' @param rot.rate
flow.twist<-function(rot.rate=0.1) {
  function(x,y,z) {
    data.frame(vx=cos(rot.rate)*x-sin(rot.rate)*z - x,
               vy=0,
               vz=cos(rot.rate)*z+sin(rot.rate)*x - z)
  }
}

flow.fields<-list(Linear=function(flow.z) flow.linear(flow.z=flow.z),
                  Shear=function(shear.slope) flow.shear.z(shear.slope),
                  Parallel=function(flow.z) flow.parallel(flow.z),
                  Explosion=function(exp.rate) flow.explosion(exp.rate),
                  Twist=function(rot.rate) flow.twist(rot.rate)
)


#library(ggplot2)
#cur.flow<-generate.flow()
#ggplot(do.call(rbind,cur.flow$frames),aes(x=POS_X,y=POS_Z,color=sample))+geom_point()
#ggplot(do.call(rbind,cur.flow$frames),aes(x=POS_X,y=POS_Z,color=as.factor(REAL_LACUNA_NUMBER),group=as.factor(REAL_LACUNA_NUMBER)))+geom_point()+geom_path()
