# n_figures describes how many figures you have next to one another, EG:
# n_figures=2 describes this type layout:
#        ____________      ____________
#       |            |    |            |
#       |  fig. Xa   |    |  fig. Xb   |
#       |____________|    |____________|
new_png <- function(filename, n_figures=1, w=6, h=w*(9/16)) {
  if (n_figures==2) {
    w=w/2
    h=h*1.3
  } else if (n_figures==3) {
    warning("n_figures=3 don't do anything")
  } else if (n_figures != 1) {
    stop("Won't allow more than 3 or fewer than 1 figures next to one anohter")
  }

  png(filename, width=w, height=h, units="in", res=300)

  margs <- par("mar")
  margs[1] = margs[1] - 2   # shrinks the bottom margin
  margs[2] = margs[2] - 1   # shrinks the left margin
  margs[3] = 0.75            # shrinks the top margin
  margs[4] = 0.75            # shrinks the right margin
  margs = margs

  # mgp setting moves axis labels in by 1 line, tick labels in by 1/2 line
  par(cex.lab=.9, cex.axis=.9, mar=margs, mgp=c(1.75*0.9,0.9/2,0))
}
