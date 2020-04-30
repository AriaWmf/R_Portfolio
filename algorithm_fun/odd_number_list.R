#odd number function

oddNumbers <- function(l, r) {
  kk <- c()
  for (i in l:r) {
    if ( (i %% 2) != 0) {
      kk <- append(kk,i)
    }
  }
  return(kk)
}

oddNumbers(3,10)
