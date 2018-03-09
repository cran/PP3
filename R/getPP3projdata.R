getPP3projdata <-
function(PP3, number){

if (missing(number))
	stop("You have to specify which set of projected data you want by specifying a number argument")
return(PP3$pdata.list[[number]])
}
