truncate.str <- function(s,l)
  substr(s,1,l)
fixwidth <- function(s, width, justify=c('left','right','center'))
  sapply(s,function(s0)
    switch(justify,
           left=paste0(s0,strrep(' ',max(width-nchar(s0),0))),
           right= paste0(strrep(' ',max(width-nchar(s0),0)),s0),
           center=paste0(strrep(' ',max(floor(width-nchar(s0),0)/2)),s0,strrep(' ',ceil(max(width-nchar(s0),0)/2)))
    ))
print.matrix <- function(x,nsmall=4,space=2,width=12){
  collab <- dimnames(x)[[2]]
  rowlab <- dimnames(x)[[1]]
  v <- as.vector(x)

  v <- round(v,nsmall)
  width <- max(width,max(sapply(v,nchar)))
  collab <- truncate.str(collab,width)
  rowlab <- truncate.str(rowlab,width)
  rlwidth <- max(nchar(rowlab))

  if((ncol(x)+1)*width+ncol(x)+space>options('width')$width){
    print.matrix(x[,1:floor(options('width')$width/width)])
    print.matrix(x[,(floor(options('width')$width/width)+1):ncol(x)])
  }else{
    v <- matrix(v,nrow=nrow(x),ncol=ncol(x))
    hdrrow <- paste0(c(paste0(rep(' ',rlwidth),collapse=''),
                       fixwidth(collab,width,justify='right')),
                     collapse=strrep(' ',space))
    cat(hdrrow)
    cat('\n')
    for(i in 1:nrow(x)){
      row <- paste0(c(fixwidth(rowlab[i],rlwidth,justify='right'),
                      fixwidth(format(v[i,]),width,justify='right')),
                    collapse=strrep(' ',space))
      cat(row)
      cat('\n')
    }
  }
}
print.coef.matrix <- function(x,se,nsmall=3,space=1,width=12){
  collab <- dimnames(x)[[2]]
  rowlab <- dimnames(x)[[1]]
  v <- as.vector(x)
  sev <- as.vector(se)
  tval <- abs(v/sev)
  sig <- c('   ','*  ','** ','***')[((tval>qnorm(0.975))*1 + (tval>qnorm(0.995))*1 + (tval>qnorm(0.999))*1)+1]
  sig <- truncate.str(sig,max(nchar(sig)))
  v <- round(v,nsmall)
  v <- paste0(format(v),sig)

  width <- max(width,max(sapply(v,nchar)))
  collab <- truncate.str(collab,width)
  rowlab <- truncate.str(rowlab,width)
  rlwidth <- max(nchar(rowlab))

  if(ncol(x)*width+rlwidth+ncol(x)*space>options('width')$width){
    j <- floor((options('width')$width-rlwidth)/(width+space))
    print.coef.matrix(x[,1:j],
                      se[,1:j],
                      nsmall=nsmall,space=space,width=width)
    cat('\n')
    print.coef.matrix(x[,(j+1):ncol(x)],
                      se[,(j+1):ncol(x)],nsmall=nsmall,space=space,width=width)
  }else{
    v <- matrix(v,nrow=nrow(x),ncol=ncol(x))
    hdrrow <- paste0(c(paste0(rep(' ',rlwidth),collapse=''),
                       fixwidth(collab,width,justify='left')),
                     collapse=strrep(' ',space))
    cat(hdrrow)
    cat('\n')
    for(i in 1:nrow(x)){
      row <- paste0(c(fixwidth(rowlab[i],rlwidth,justify='right'),
                      fixwidth(format(v[i,]),width,justify='right')),
                    collapse=strrep(' ',space))
      cat(row)
      cat('\n')
    }
  }
}
