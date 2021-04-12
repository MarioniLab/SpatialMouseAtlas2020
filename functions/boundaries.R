# helper function
subsetV = Vectorize(
  function(M, i, j) {
    M[i,j]
  },
  vectorize.args = c("i", "j")
)
