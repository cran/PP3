pdataplot <-
function (PP3manyobj, number, colvec=1, lab=NULL, ...) 
{
	oldpar <- par(mfrow=c(2,2), pty="s")

	if (is.null(lab))	{
		plot(PP3manyobj$pdata.list[[number]][1,], PP3manyobj$pdata.list[[number]][2,], col=colvec, xlab="Direction 1", ylab="Direction 2", ...)
		plot(PP3manyobj$pdata.list[[number]][1,], PP3manyobj$pdata.list[[number]][3,], col=colvec, xlab="Direction 1", ylab="Direction 3", ...)
		plot(PP3manyobj$pdata.list[[number]][2,], PP3manyobj$pdata.list[[number]][3,], col=colvec, xlab="Direction 2", ylab="Direction 3", ...)
		}
	else	{
		plot(PP3manyobj$pdata.list[[number]][1,], PP3manyobj$pdata.list[[number]][2,], xlab="Direction 1", ylab="Direction 2", type="n", ...)
		text(PP3manyobj$pdata.list[[number]][1,], PP3manyobj$pdata.list[[number]][2,], col=colvec, lab=lab)

		plot(PP3manyobj$pdata.list[[number]][1,], PP3manyobj$pdata.list[[number]][3,], xlab="Direction 1", ylab="Direction 3", type="n", ...)
		text(PP3manyobj$pdata.list[[number]][1,], PP3manyobj$pdata.list[[number]][3,], col=colvec, lab=lab)

		plot(PP3manyobj$pdata.list[[number]][2,], PP3manyobj$pdata.list[[number]][3,], xlab="Direction 2", ylab="Direction 3", type="n", ...)
		text(PP3manyobj$pdata.list[[number]][2,], PP3manyobj$pdata.list[[number]][3,], col=colvec, lab=lab)
	}

	par(oldpar)
}
