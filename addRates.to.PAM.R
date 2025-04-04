#Function to add rates to PAM

addRates.to.PAM<- function(PAM, rates){
  
  map.names <- colnames(PAM)
  
  x <- intersect(rates$SPECIES,map.names) ## get names in common between both tables
  
  rates.sub <- rates[rates$SPECIES%in%x,] ## get rates for common species
  
  x <- c('Longitude.x','Latitude.y',x)
  
  PAM.sub <- PAM[,colnames(PAM)%in%x] ## make PAM only for common species
  
  n <- rates.sub$SPECIES ### save species names
  
  rates.matrix <- as.matrix(t(rates.sub[,-(1)])) ## put species as columns to match siteXspp
  
  colnames(rates.matrix) <- n
  
  x <- intersect(rates$SPECIES,map.names) ## get names in common between both tables
  
  mat = matrix(ncol=length(x), nrow=nrow(PAM.sub)) ## create empty matrix to fill in loop
  
  #### Loop to replace the values (1's) on the map by speciation rates ####
  for (i in 1:length(x)){
    col <- PAM.sub[,x[i]]
    col[col==1] <- rates.matrix[4,x[i]] ###This number three refers to the third row in the rates.matrix object where the rates are
    mat[, i] = col
  }
  
  colnames(mat) <- x #### put column names in the matrix
  
  PAM.sub <- cbind(PAM[,1:2],mat) ### put back the xy coordinates
  
  PAM.sub  <- as.data.frame(apply(PAM.sub, 2, as.numeric)) 
  
  PAM.sub [PAM.sub  == 0] <- NA ### replace zeros by NA
  
  return(PAM.sub)
}



