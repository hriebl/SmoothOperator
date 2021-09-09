set.seed(1337)

n <- 100
k <- 10

for (bs in supported) {
  prefix <- paste0("Basis '", bs, "':")

  df <- data.frame(x = runif(n))
  smt <- Smooth$new("x", df, bs, k)
  x1 <- smt$construct()

  smt$add_knots()
  x2 <- smt$construct()

  test_that(paste(prefix, "Adding knots works"), {
    expect_true(is.data.frame(smt$knots))
    expect_equal(names(smt$knots), "x")

    expect_equal(x2$design_matrix, x1$design_matrix)
    expect_equal(x2$penalty_matrices, x1$penalty_matrices)
    expect_equal(x2$ranks, x1$ranks)
  })

  smt$add_centering_constraint()
  x2 <- smt$construct()

  test_that(paste(prefix, "Adding centering constraint works"), {
    expect_true(is.matrix(smt$constraints))
    expect_equal(nrow(smt$constraints), 1)

    expect_equal(x2$design_matrix, x1$design_matrix)
    expect_equal(x2$penalty_matrices, x1$penalty_matrices)
    expect_equal(x2$ranks, x1$ranks)
  })

  smt$add_point_constraints(data.frame(x = 0.5))
  x2 <- smt$construct()

  test_that(paste(prefix, "Adding point constraint works"), {
    expect_true(is.matrix(smt$constraints))
    expect_equal(nrow(smt$constraints), 2)

    expect_equal(ncol(x2$design_matrix), ncol(x1$design_matrix) - 1)
  })
}
