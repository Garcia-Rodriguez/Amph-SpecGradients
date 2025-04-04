##Function to calculate mean, sd, cv and richness across lines in a large PAM
df <- data.frame(
  col1 = c("A", "B", "C", "D"),
  col2 = c("J", "K", "L", "M"),
  col1 = c(10, 20, 30, 40),
  col2 = c(15, 25, "NA", 35),
  col3 = c(20, "NA", "NA", "NA"),
  col4 = c(40, 42, 48, 46)
  )

PAM_mean.sd.cv<- function (df){
  ##Coordinates
  coords<-df[, 1:2]
  
  ##Speciation values
  values<-df[, 3:length(df)]
  values <- sapply(values, as.numeric)
  
  ##Mean speciation per grid-cell
  mean<-rowMeans(values, na.rm = TRUE)
  
  ##Standard deviation per grid-cell
  sd <- apply(values, 1, function(row) sd(row, na.rm = TRUE))
  
  ##Coefficient of variation per grid-cell
  cv<- (sd / mean) * 100
  
  ##Calculate species richness
  rich<- apply(values, 1, function(row) sum(!is.na(row)))
  
  ##Extract max and min per grid cell
  max<- apply(values, 1, max, na.rm = TRUE)
  min<- apply(values, 1, min, na.rm = TRUE)
  
  ##Calculate speciation range
  range<- max-min
  range[range == 0] <- NA
  
  ##Merge all stats
  return(data.frame(cbind(coords, mean, sd, cv, rich,max, min, range)))  
}

test<-PAM_mean.sd.cv(df)
test
