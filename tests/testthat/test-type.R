y=seq(-50,100,length.out = 30)
ncells=3

type=type(y=y,ncells=ncells)
Datatype=type

test_that("right length", {
  expect_equal(length(Datatype), length(y))
})

test_that("right number of cells", {
  expect_equal(ncells, max(Datatype))
  expect_equal(1, min(Datatype))
})
