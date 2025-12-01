library(testthat)

test_that("edge_length and triangle_area compute expected values", {
    p1 <- c(0, 0)
    p2 <- c(3, 0)
    p3 <- c(0, 4)

    expect_equal(AlphaShapePackage:::edge_length(p1, p2), 3)
    expect_equal(AlphaShapePackage:::edge_length(p1, p3), 4)
    expect_equal(AlphaShapePackage:::triangle_area(p1, p2, p3), 6)
})

test_that("circumcircle returns correct center and radius for right triangle", {
    p1 <- c(0, 0)
    p2 <- c(1, 0)
    p3 <- c(0, 1)

    cc <- AlphaShapePackage:::circumcircle(p1, p2, p3)
    # center should be (0.5, 0.5)
    expect_equal(round(cc$center, 6), c(0.5, 0.5))
    # radius should be sqrt(0.5)
    expect_equal(round(cc$r, 6), round(sqrt(0.5), 6))
})

test_that("polygon_area computes shoelace formula correctly", {
    square <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE)
    expect_equal(AlphaShapePackage:::polygon_area(square), 1)

    tri <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
    expect_equal(AlphaShapePackage:::polygon_area(tri), 0.5)
})

test_that("point_in_circle excludes boundary points (strict <)", {
    center <- c(0, 0)
    r <- 1
    inside <- c(0.5, 0)
    on_boundary <- c(1, 0)

    expect_true(AlphaShapePackage:::point_in_circle(inside, center, r))
    expect_false(AlphaShapePackage:::point_in_circle(on_boundary, center, r))
})
