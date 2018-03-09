PP3slowDF3 <-
function (Pvec, the.init, maxrow, k, maxcol, n, text) 
{

avec <- Pvec[1:maxrow]
bvec <- Pvec[(maxrow+1):(2*maxrow)]
cvec <- Pvec[(2*maxrow+1):(3*maxrow)]

answer <- PP3ix3dvsFromTU(the.init=the.init, avec=avec, bvec=bvec, cvec=cvec,
  maxrow=maxrow, k=k, maxcol=maxcol, n=n, text=text, type = "deriv")

return(answer)

}
