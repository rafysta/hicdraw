test_that("hic200 index offsets are cumulative", {
  idx <- hic200_index(pombe_chrom_sizes, res = 200)
  expect_equal(idx$bin_offset[1], 0L)
  expect_equal(idx$bin_offset[2], as.integer(ceiling(5579133 / 200)))
  expect_equal(idx$nbin[3], as.integer(ceiling(2452883 / 200)))
})

test_that("bins_for_region converts chrII:800k-900k correctly", {
  idx <- hic200_index(pombe_chrom_sizes, res = 200)
  bb <- bins_for_region(idx, "II", 800000, 900000, res = 200)
  off2 <- idx$bin_offset[idx$chr == "II"]
  expect_equal(bb[1], off2 + 800000L / 200L)
  expect_equal(bb[2], off2 + 900000L / 200L)
})

test_that("slice_hic200 symmetrises and maps back to bp", {
  # tiny synthetic hic200 object: two chromosomes of 1000 bp at 200 bp
  cs <- data.frame(chr = c("A", "B"), length = c(1000, 1000))
  idx <- hic200_index(cs, res = 200)
  # chrA occupies bins 0..4, chrB bins 5..9
  dt <- data.table::data.table(
    bin1 = c(0L, 0L, 1L),
    bin2 = c(0L, 1L, 2L),
    score = c(5, 3, 2)
  )
  obj <- structure(list(dt = dt, index = idx, res = 200L, file = "x"),
                   class = "hic200")
  sl <- slice_hic200(obj, "A", 0, 800)
  # diagonal (0,0) once; off-diagonal (0,1) and (1,2) duplicated
  expect_equal(nrow(sl), 5)
  expect_true(all(sl$x %% 200 == 0))
  expect_equal(sum(sl$counts), 5 + 3 + 3 + 2 + 2)
})
