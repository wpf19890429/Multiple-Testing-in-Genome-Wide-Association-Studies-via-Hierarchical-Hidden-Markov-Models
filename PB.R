PB<-function(Methods,Rep)
{

################################ Show Progress Bar
	library(progress)
	pb <- progress_bar$new(
         format = paste(Methods,":completed [:bar] :percent, Execute time::elapsed",sep=""),
         total = Rep, clear = FALSE, width= 60)
	return(pb)
}