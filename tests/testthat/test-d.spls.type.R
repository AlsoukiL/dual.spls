test_that("d.spls.type works", {

y=seq(-50,100,length.out = 30)
ncells=3

type=d.spls.type(y=y,ncells=ncells)
Datatype=type

# right length
  expect_equal(length(Datatype), length(y))


# right number of cells"
  expect_equal(ncells, max(Datatype))
  expect_equal(1, min(Datatype))
})
