# modified from image.legend by Martin Maechler
image.legend <- function(x,y, zlim, at.z = NULL, col = heat.colors(12), legnd=NULL,
             lwd = max(3,32/length(col)), bg = NA, bty = "", ...)
  ## * kein y.i -- Benutzer soll rein ueber lwd steuern; sollte reichen.
  ## * legnd koennte interessant sein, falls Text geschrieben werden soll
  ##   (weiss mal wieder nicht, wie man aus legnd legend als option
  ##     macht)
  ## * lwd wird per default in Abh. von col gewaehlt.
{
    ## Purpose:
    ## Authors: Martin Maechler,   9 Jul 2001
    ##          Martin Schlather, 24 Jul 2001

  if (!is.null(legnd) && is.null(at.z))
      stop("at.z must be given if legnd is") ## falls legnd darf at.z
    ##                                nicht automatisch gewaehlt werden

    if(!is.numeric(zlim) || zlim[1] > zlim[2])
        stop("`zlim' must be numeric; zlim[1] <= zlim[2]")
    if(is.null(at.z)) {
        ## hier ein Versuch in Abhaengigkeit von n
        ## die Anzahl der labels zu bestimmen:
        n <- min(5, max(1,length(col)/10))
        at.z <- pretty(zlim,n=n,min.n=max(n %/% 3,1))

        ## es sieht nicht schoen aus, wenn pretty die letzte oder
        ## erste zahl weit ausserhalb des zlim legt.
        ## heuristisch nur 25%  (oder so) ueberschreitung bzgl
        ## intervalllaenge zulassen:
        tol <- diff(at.z)[1] / 4
        at.z <- at.z[(at.z>=zlim[1]-tol) & (at.z<=zlim[2]+tol)]
      }
    if(!is.numeric(at.z) || is.unsorted(at.z))
        stop("`at.z' must be numeric non-decreasing")
    n.at <- length(at.z)
    nc   <- length(col)
    if(n.at >= nc)
        stop("length(at.z) must be (much) smaller than length(col)")
    dz <- diff(zlim)
    ## The colors must run equidistantly from zlim[1] to zlim[2];
    ## col[i] is for z-interval zlim[1] + [i-1, i) * dz/nc  ; i = 1:nc
    ## i.e., an at.z[] value z0 is color i0 = floor(nc * (z0 - zlim[1])/dz)
    at.i <- floor(nc * (at.z - zlim[1])/dz )
    ## Possibly extend colors by `background' to the left and right
    bgC <- if(is.null(bg)) NA else bg
    if((xtra.l <- 1 - at.i[1]) > 0) {
        at.i <- at.i + xtra.l
        col <- c(rep(bgC, xtra.l), col)
    }
    if((xtra.r <- at.i[n.at] - nc) > 0)
        col <- c(col, rep(bgC, xtra.r))
    lgd <- character(length(col))

    ## folgende if-Anweisung ist neu:
    if (is.null(legnd)) lgd[at.i] <-format(at.z, dig = 3)
    else {
      if (length(legnd)!=length(at.z))
        stop("at.z and legnd must have the same length")
      lgd[at.i] <- legnd
    }
    if((V <- R.version)$major <= 1 && V$minor <= 3.0 && V$status == "")
{
        ## stop-gap fix around the bug that "NA" is not a valid color:
        if(is.na(bgC)) {
            lgd <- lgd[!is.na(col)]
            col <- col[!is.na(col)]
        }
    }
    legend(x,y, legend = rev(lgd), col = rev(col), y.i = lwd/16, bty = bty, lwd = lwd, bg = bg, ...)
}
