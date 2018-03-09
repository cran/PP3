plot.PP3 <-
function(x, number, main, nbig=10, colvec=1, lab=NULL, ...){

if (missing(number))	{
	if (missing(main))
	    main <- "Projection Index Value Histogram"

	ix3sizeplot(x, main=main, nbig=nbig)
	}
else	{
	if (missing(main))
		main=""
	pdataplot(x, number=number, colvec=colvec, lab=lab, main=main, ...) 
}
}
