#' Alpha Shape Data
#'
#' @name data
#' @rdname data
#' @importFrom dbscan dbscan kNNdist
#' @importFrom stats dist
NULL

#' visualize_alpha_shape_2D
#' 
#' 2D Alpha Shape Visualization
#' @param alpha_shape_result The result from \code{alpha_shape_2D}.
#' @param points The original set of points used to compute the alpha shape.
#' @param title Title for the plot.
#' @param r_outer Outer radius for true shape outline (if applicable).
#' @param r_inner Inner radius for true shape outline (if applicable).
#' @param add Logical indicating whether to add to an existing plot.
#' @param shape The type of true shape ("ring", "concentric", "squares").
#' @param num_rings Number of rings (for "concentric" or "squares" shapes).
#' @param spacing Spacing between rings (for "concentric" or "squares" shapes).
#' @param ring_width Width of each ring (for "concentric" or "squares" shapes).
#' @param save_path Optional file path to save the plot as a PNG.
#' 
#' @return NULL. The function produces a plot.
visualize_alpha_shape_2D <- function(alpha_shape_result, points, title = "Alpha Shape (2D)",
                                     r_outer = NULL, r_inner = NULL, add = FALSE, shape = "ring",
                                     num_rings = NULL, spacing = NULL, ring_width = NULL,
                                     save_path = NULL) {
    if (!is.null(save_path)) {
        png(filename = save_path, width = 800, height = 800)
        on.exit(dev.off())
    }

    plots <- alpha_shape_result
    loops <- plots$polygons

    # If not adding to an existing plot, initialize the plot window with no background
    if (!add) {
        plot(points,
            col = "black", pch = 19, main = title, asp = 1,
            xlab = "x", ylab = "y",
            bg = NA
        ) # transparent background
    } else {
        points(points, col = "black", pch = 19)
    }

    if (length(loops) == 0) {
        return(invisible(NULL))
    }

    # ---- Alpha-shape outlines (blue) ----
    for (i in seq_along(loops)) {
        loop <- loops[[i]]
        polygon(loop, border = "blue", lwd = 2, col = NA)
    }

    # ---- True shape outlines (red) ----
    if (shape == "ring" && !is.null(r_outer) && !is.null(r_inner)) {
        # outline between the two circles to represent the ring
        theta <- seq(0, 2 * pi, length.out = 400)
        lines(r_outer * cos(theta), r_outer * sin(theta),
            col = "red", lwd = 2
        )

        lines(r_inner * cos(theta), r_inner * sin(theta),
            col = "red", lwd = 2
        )
    } else if (shape == "concentric" && !is.null(num_rings) && !is.null(ring_width) && !is.null(spacing)) {
        # outline for multiple concentric rings
        r_inner <- spacing
        for (j in 1:num_rings) {
            r_outer <- r_inner + ring_width
            theta <- seq(0, 2 * pi, length.out = 400)
            lines(r_outer * cos(theta), r_outer * sin(theta), col = "red", lwd = 2)
            lines(r_inner * cos(theta), r_inner * sin(theta), col = "red", lwd = 2)
            r_inner <- r_outer + spacing
        }
    } else if (shape == "squares" && !is.null(spacing) && !is.null(ring_width) && !is.null(num_rings)) {

        for (i in 1:num_rings){

            square_width <- ring_width

            r_inner <- spacing + 2 * (i - 1) * (square_width + spacing)
            r_outer <- spacing + 2 * (i - 1) * (square_width + spacing) + 2 * square_width
            # outline between the two squares to represent the square ring
            square_outer <- rbind(
                c(-r_outer, -r_outer),
                c(r_outer, -r_outer),
                c(r_outer, r_outer),
                c(-r_outer, r_outer),
                c(-r_outer, -r_outer)
            )
            lines(square_outer, col = "red", lwd = 2)

            square_inner <- rbind(
                c(-r_inner, -r_inner),
                c(r_inner, -r_inner),
                c(r_inner, r_inner),
                c(-r_inner, r_inner),
                c(-r_inner, -r_inner)
            )
            lines(square_inner, col = "red", lwd = 2)

            
        }
        
        # outline between the two squares to represent the square ring
        square_outer <- rbind(
            c(-r_outer, -r_outer),
            c(r_outer, -r_outer),
            c(r_outer, r_outer),
            c(-r_outer, r_outer),
            c(-r_outer, -r_outer)
        )
        lines(square_outer, col = "red", lwd = 2)

        square_inner <- rbind(
            c(-r_inner, -r_inner),
            c(r_inner, -r_inner),
            c(r_inner, r_inner),
            c(-r_inner, r_inner),
            c(-r_inner, -r_inner)
        )
        lines(square_inner, col = "red", lwd = 2)
    } else {
        invisible(NULL)
    }


    invisible(NULL)
}

#' calculate_area_volume
#' 
#' Calculate the area (2D) or volume (3D) of an alpha shape
#' @param alpha_shape The output from \code{alpha_shape_2D} or \code{alpha_shape_3D}.
#' @param dim The dimension of the alpha shape (2 or 3).
#' 
#' @return The area (2D) or volume (3D) of the alpha shape.
calculate_area_volume <- function(alpha_shape, dim = 2) {
    if (dim == 2) {
        if (is.null(alpha_shape$polygons) || length(alpha_shape$polygons) == 0) {
            return(0)
        }

        # Calculate signed areas for each polygon
        signed_areas <- sapply(alpha_shape$polygons, function(poly) {
            if (nrow(poly) < 3) {
                return(0)
            }
            x <- poly[, 1]
            y <- poly[, 2]
            0.5 * sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1]))
        })

        if (all(signed_areas == 0)) {
            return(0)
        }

        # Use signed areas - this handles holes correctly when orientations are proper
        total_signed_area <- sum(signed_areas)

        # If total is very close to zero, try alternative calculation
        if (abs(total_signed_area) < 1e-10) {
            # Use absolute areas - assume largest is outer boundary
            abs_areas <- abs(signed_areas)
            if (length(abs_areas) == 1) {
                result <- abs_areas[1]
            } else {
                # Largest polygon is outer boundary, others are holes
                result <- max(0, abs_areas[1] - sum(abs_areas[-1]))
            }
        } else {
            result <- abs(total_signed_area)
        }

        return(result)
    } else {
        stop("3D alpha shape volume not yet implemented.")
    }
}

#' calculate_surface_perimeter
#' 
#' Calculate surface area (3D) or perimeter (2D) of the alpha shape
#' @param alpha_shape The output from \code{alpha_shape_2D} or \code{alpha_shape_3D}.
#' @param dim The dimension of the alpha shape (2 or 3).
#' 
#' @return The surface area (3D) or perimeter (2D) of the alpha shape.
calculate_surface_perimeter <- function(alpha_shape, dim = 2) {
    # Function to calculate surface area (3D) or perimeter (2D) of the alpha shape
    # alpha_shape: the output from alpha_shape_2D or alpha_shape_3D

    if (dim == 2) {
        return(polygon_perimeter(alpha_shape))
    } else {
        return(tetrahedron_surface_area(alpha_shape))
    }
}

#' torus_test_generation_3D
#' 
#' Generate test data: Torus in 3D
#' @param num_points Number of points to generate.
#' @param R Major radius.
#' @param r Minor radius.
#' @param noise Standard deviation of Gaussian noise to add.
#' 
#' @return A numeric matrix of points (n x 3).
torus_test_generation_3D <- function(num_points = 10000, R = 3, r = 1, noise = 0.05) {
    # Function to generate points on a torus
    # num_points: number of points to generate
    # R: major radius
    # r: minor radius

    u <- runif(num_points, 0, 2 * pi)
    v <- runif(num_points, 0, 2 * pi)

    rho <- r * sqrt(runif(num_points, 0, 1))

    x <- (R + rho * cos(v)) * cos(u)
    y <- (R + rho * cos(v)) * sin(u)
    z <- rho * sin(v)

    if (noise > 0) {
        x <- x + rnorm(num_points, 0, sd = noise)
        y <- y + rnorm(num_points, 0, sd = noise)
        z <- z + rnorm(num_points, 0, sd = noise)
    }

    points <- cbind(x, y, z)
    return(points)
}

#' ring_test_generation_2D
#' 
#' Generate test data: Ring in 2D
#' @param num_points Number of points to generate.
#' @param r_outer Outer radius.
#' @param r_inner Inner radius.
#' @param noise Proportion of noise points to add.
#' 
#' @return A numeric matrix of points (n x 2).
ring_test_generation_2D <- function(num_points = 10000, r_outer = 3, r_inner = 1, noise = 0.05) {
    # Function to generate points on a ring shape
    # num_points: number of points to generate
    # r_outer: outer radius
    # r_inner: inner radius

    n_noise <- floor(noise * num_points)
    n_ring <- num_points - n_noise

    theta <- runif(n_ring, 0, 2 * pi)
    rho <- r_inner + (r_outer - r_inner) * runif(n_ring, 0, 1)

    x_ring <- rho * cos(theta)
    y_ring <- rho * sin(theta)

    if (noise > 0) {
        x_noise <- runif(n_noise, -r_outer, r_outer)
        y_noise <- runif(n_noise, -r_outer, r_outer)

        x <- c(x_ring, x_noise)
        y <- c(y_ring, y_noise)
    } else {
        x <- x_ring
        y <- y_ring
    }

    points <- cbind(x, y)
    return(points)
}

#' concentric_rings_2D
#' 
#' Generate test data: Concentric rings in 2D
#' @param num_rings Number of rings.
#' @param spacing Spacing between rings.
#' @param ring_width Width of each ring.
#' @param total_points Total number of points to generate.
#' @param noise Proportion of noise points to add.
#' 
#' @return A numeric matrix of points (n x 2).
concentric_rings_2D <- function(num_rings = 3, spacing = 1, ring_width = 1.0, total_points = 1000, noise = 0.01) {
    outer_edge <- num_rings * spacing + ring_width + spacing
    n_noise <- floor(noise * total_points)
    n_ring_points <- total_points - n_noise

    # equal ratio of points per ring

    points_per_ring_ratios <- seq(1, num_rings) / sum(seq(1, num_rings))

    points <- matrix(0, nrow = total_points, ncol = 2)
    current_idx <- 1

    r_inner <- spacing

    # Generate ring points first
    for (i in 1:num_rings) {
        r_outer <- r_inner + ring_width

        points_per_ring <- points_per_ring_ratios[i] * n_ring_points

        r <- runif(points_per_ring, r_inner, r_outer)
        theta <- runif(points_per_ring, 0, 2 * pi)
        x <- r * cos(theta)
        y <- r * sin(theta)

        end_idx <- current_idx + points_per_ring - 1
        points[current_idx:end_idx, ] <- cbind(x, y)
        current_idx <- end_idx + 1

        r_inner <- r_outer + spacing
    }

    # Generate noise points in remaining slots
    if (n_noise > 0) {
        x_noise <- runif(n_noise, -outer_edge, outer_edge)
        y_noise <- runif(n_noise, -outer_edge, outer_edge)
        points[current_idx:(current_idx + n_noise - 1), ] <- cbind(x_noise, y_noise)
    }

    return(points)
}

#' alternating_squares_pointcloud_2D
#' 
#' Generate test data: Alternating squares in 2D
#' @param total_points Total number of points to generate.
#' @param num_squares Number of squares.
#' @param square_width Width of each square.
#' @param spacing Spacing between squares.
#' @param noise Proportion of noise points to add.
#' 
#' @return A numeric matrix of points (n x 2).
alternating_squares_pointcloud_2D <- function(total_points = 1000, num_squares = 5, square_width = 0.5, spacing = 1, noise = 0.05) {
    points <- matrix(0, nrow = 0, ncol = 2)

    total_points <- ceiling(total_points / 2)

    if (noise > 0) {
        noisy_points <- floor(noise * total_points)
        total_points <- total_points - noisy_points
    }

    points_per_square_ratio <- seq(1, num_squares) / sum(seq(1, num_squares))

    center_x <- 0
    center_y <- 0

    for (i in 1:num_squares) {
        inner_k <- spacing + 2 * (i - 1) * (square_width + spacing)
        outer_k <- spacing + 2 * (i - 1) * (square_width + spacing) + 2 * square_width

        num_square_points <- floor(points_per_square_ratio[i] * total_points)

        square_points <- matrix(0, nrow = 0, ncol = 2)

        # generate two numbers that are uniform in the hollow square and then generate a othrogonal point
        # in the actual bounds

        for (j in 1:num_square_points) {
            rand_x_1 <- runif(1, -outer_k, outer_k)

            y_bounds_1 <- runif(1, inner_k, outer_k)
            y_sign_1 <- sample(c(-1, 1), 1)
            rand_y_1 <- y_sign_1 * y_bounds_1

            rand_y_2 <- runif(1, -outer_k, outer_k)

            x_bounds_2 <- runif(1, inner_k, outer_k)
            x_sign_2 <- sample(c(-1, 1), 1)
            rand_x_2 <- x_sign_2 * x_bounds_2

            square_points <- rbind(square_points, cbind(rand_x_1, rand_y_1))
            square_points <- rbind(square_points, cbind(rand_x_2, rand_y_2))
        }

        points <- rbind(points, square_points)
    }

    # Add noise points
    if (noise > 0) {
        x_noise <- runif(noisy_points, min(points[, 1]) - spacing, max(points[, 1]) + spacing)
        y_noise <- runif(noisy_points, min(points[, 2]) - spacing, max(points[, 2]) + spacing)
        noise_points <- cbind(x_noise, y_noise)
        points <- rbind(points, noise_points)
    }

    return(points)
}

#' alternating_squares_area_2D
#' 
#' Calculate area of alternating squares shape in 2D
#' @param num_squares Number of squares.
#' @param square_width Width of each square.
#' @param spacing Spacing between squares.
#' @param center_x X-coordinate of the center.
#' @param center_y Y-coordinate of the center.
#' 
#' @return The total area of the alternating squares.
square_area_2D <- function(num_squares = 5, square_width = 0.5, spacing = 1, center_x = 0, center_y = 0) {
    total_area <- 0
    for (i in 0:(num_squares - 1)) {
        x_inner <- center_x + spacing / 2 + (i - 1) * (square_width + spacing)
        x_outer <- center_x + spacing / 2 + (i - 1) * (square_width + spacing) + square_width
        y_inner <- center_y + spacing / 2 + (i - 1) * (square_width + spacing)
        y_outer <- center_y + spacing / 2 + (i - 1) * (square_width + spacing) + square_width

        area_square <- (2 * x_outer * 2 * y_outer) - (2 * x_inner * 2 * y_inner)
        total_area <- total_area + area_square
    }
    return(total_area)
}

#' alternating_squares_polygon
#' 
#' Generate polygon structure for alternating squares in 2D
#' @param num_squares Number of squares.
#' @param square_width Width of each square.
#' @param spacing Spacing between squares.
#' @param center_x X-coordinate of the center.
#' @param center_y Y-coordinate of the center.
#' 
#' @return A list of matrices representing the polygons.
alternating_squares_polygon <- function(num_squares = 5,
                                        square_width = 0.5,
                                        spacing = 1,
                                        center_x = 0,
                                        center_y = 0) {
    polygon_list <- list()

    for (i in 1:num_squares) {
        inner_k <- spacing / 2 + (i - 1) * (square_width + spacing)
        outer_k <- inner_k + square_width

        # Outer boundary (CCW)
        outer <- rbind(
            c(center_x - outer_k, center_y + outer_k),
            c(center_x + outer_k, center_y + outer_k),
            c(center_x + outer_k, center_y - outer_k),
            c(center_x - outer_k, center_y - outer_k),
            c(center_x - outer_k, center_y + outer_k) # close
        )

        # Inner boundary (CW for a "hole")
        inner <- rbind(
            c(center_x - inner_k, center_y + inner_k),
            c(center_x - inner_k, center_y - inner_k),
            c(center_x + inner_k, center_y - inner_k),
            c(center_x + inner_k, center_y + inner_k),
            c(center_x - inner_k, center_y + inner_k) # close
        )
        inner <- inner[nrow(inner):1, ] # reverse to CW

        # Combine into a single polygon structure
        poly <- rbind(outer, inner)

        polygon_list[[i]] <- poly
    }

    return(polygon_list)
}
#' pitchfork_bifurcation_datacloud
#' 
#' Generate test data: Pitchfork bifurcation in 2D
#' @param num_points Number of points to generate.
#' @param noise Proportion of noise points to add.
#' @param direction Direction of the bifurcation ("up" or "down").
#' 
#' @return A numeric matrix of points (n x 2).
pitchfork_bifurcation_datacloud <- function(num_points = 1000, noise = 0.05, direction = "up") {
    # base equation: dx = r*x - x^3

    points <- matrix(0, nrow = num_points, ncol = 2)

    positive_points <- matrix(0, nrow = 0, ncol = 2)
    negative_points <- matrix(0, nrow = 0, ncol = 2)

    num_noise <- floor(noise * num_points)

    for (i in 1:num_points) {
        x <- runif(1, -2, 2)
        r <- runif(1, -2, 2)
        y <- r * x - x^3
        # assign (r,x) to (x,y) depending on sign(y)
        if (y >= 0) {
            positive_points <- rbind(positive_points, c(x, y))
        } else {
            negative_points <- rbind(negative_points, c(x, y))
        }
    }

    # add noise points
    if (num_noise > 0) {
        x_noise <- runif(num_noise, -2, 2)
        y_noise <- runif(num_noise, -2, 2)
        noise_points <- cbind(x_noise, y_noise)

        positive_points <- rbind(positive_points, noise_points)
        negative_points <- rbind(negative_points, noise_points)
    } else {
        positive_points <- positive_points
        negative_points <- negative_points
    }

    if (direction == "up") {
        return(positive_points)
    } else if (direction == "down") {
        return(negative_points)
    } else {
        return(points)
    }
}

#' polygon_pitchfork
#' 
#' Generate polygon representing pitchfork bifurcation shape in 2D
#' @param x_min Minimum x-value.
#' @param x_max Maximum x-value.
#' @param r_min Minimum r-value.
#' @param r_max Maximum r-value.
#' 
#' @return A matrix of polygon coordinates.
polygon_pitchfork <- function(x_min = -2, x_max = 2, r_min = -2, r_max = 2) {
    # create polygon representing the pitchfork bifurcation shape
    num_points <- 200
    x_seq <- seq(x_min, x_max, length.out = num_points)
    upper_branch <- cbind(x_seq, (r_max * x_seq - x_seq^3))
    lower_branch <- cbind(rev(x_seq), (r_min * rev(x_seq) - rev(x_seq)^3))

    polygon_coords <- rbind(upper_branch, lower_branch)
    return(polygon_coords)
}

#' area_concentric_rings_2D
#' 
#' Calculate area of concentric rings shape in 2D
#' @param num_rings Number of rings.
#' @param spacing Spacing between rings.
#' @param ring_width Width of each ring.
#' 
#' @return The total area of the concentric rings.
area_concentric_rings_2D <- function(num_rings = 3, spacing = 1, ring_width = 0.5) {
    total_area <- 0
    r_inner <- spacing
    for (i in 1:num_rings) {
        r_outer <- r_inner + ring_width
        area_ring <- pi * (r_outer^2 - r_inner^2)
        total_area <- total_area + area_ring
        r_inner <- r_outer + spacing
    }
    return(total_area)
}

#' heart_shape_datacloud
#' 
#' Generate test data: Heart shape in 2D
#' @param num_points Number of points to generate.
#' @param num_rings Number of rings.
#' @param spacing Spacing between rings.
#' @param ring_width Width of each ring.
#' @param noise Proportion of noise points to add.
#' 
#' @return A numeric matrix of points (n x 2).
heart_shape_datacloud <- function(num_points = 1000, num_rings = 3, spacing = 1, ring_width = 0.5, noise = 0.05) {
    points <- matrix(0, nrow = num_points, ncol = 2)

    if (noise > 0) {
        pass
    }
}

#' error_in_2D
#' 
#' Calculate percentage error between true area and estimated area in 2D
#' @param true_area The true area.
#' @param estimated_area The estimated area.
#' 
#' @return The percentage error.
error_in_2D <- function(true_area, estimated_area) {
    return(abs(true_area - estimated_area) / true_area * 100)
}

#' to_sf_poly
#' 
#' Convert list of rings to sf polygon
#' @param rings A list of matrices, each representing a ring (Mx2).
#' 
#' @return An sf polygon object.
to_sf_poly <- function(rings) {
    # rings is a list of matrices, each ring Mx2
    
    fix_ring <- function(mat) {
        mat <- as.matrix(mat)

        # must have at least 3 distinct points
        if (nrow(mat) < 3) stop("Ring has too few points")

        # force closure exactly
        if (!all(mat[1,] == mat[nrow(mat),])) {
            mat[nrow(mat),] <- mat[1,]
        }

        # ensure 4 rows at minimum (p1,p2,p3,p1)
        if (nrow(mat) < 4) {
            mat <- rbind(mat, mat[1,])
        }

        mat
    }

    closed_rings <- lapply(rings, fix_ring)

    st_sfc(st_polygon(closed_rings))
}

#' calculate_KL_divergence
#' 
#' Calculate KL divergence between true shape and estimated alpha shape
#' @param true_shape The true shape representation (e.g., polygon for 2D).
#' @param estimated_shape The estimated alpha shape representation.
#' @param grid A data frame of grid points with columns 'x' and 'y'.
#' @param dim The dimension (2 or 3).
#' @param grid_res Resolution of the grid for Monte Carlo integration.
#' @param num_rings Number of rings (for ring shapes).
#' @param spacing Spacing between rings.
#' @param width Width of each ring.
#' @param shape The shape type ("concentric" or "squares").
#' 
#' @return The KL divergence value.
calculate_KL_divergence <- function(true_shape, estimated_shape, grid, dim = 2, grid_res = 500, num_rings = 3, spacing = 3.0, width = 1.5, shape = "concentric") {
    # Function to calculate KL divergence between true shape and estimated alpha shape
    # true_shape: the true shape representation (e.g., polygon for 2D)
    # estimated_shape: the estimated alpha shape representation
    # dim: dimension (2 or 3)

    if (dim == 2) {
        if (length(true_shape$polygons) == 0 || length(estimated_shape$polygons) == 0) {
            return(NA)
        }

        
        true_sf <- to_sf_poly(true_shape$polygons)
        est_sf <- to_sf_poly(estimated_shape$polygons)


        # bounding box
        all_points <- rbind(do.call(rbind, true_shape$polygons), do.call(rbind, estimated_shape$polygons))
        x_range <- range(all_points[, 1])
        y_range <- range(all_points[, 2])

        x_seq <- seq(x_range[1], x_range[2], length.out = grid_res)
        y_seq <- seq(y_range[1], y_range[2], length.out = grid_res)

        # create a Monte Carlo grid

        # grid <- expand.grid(x = x_seq, y = y_seq)



        # helper to check if points are in any polygon
        points_in_sf <- function(points_df, sf_polygons) {
            library(sf)
            pts <- st_as_sf(points_df, coords = c("x", "y"), crs = NA)
            result <- st_intersects(pts, sf_polygons, sparse = FALSE)
            apply(result, 1, any)
        }

        if (shape == "squares") {
            inside_est <- points_in_sf(grid, est_sf)
            inside_true <- c()
            # for grid points, check if in any of the square rings using the math

            inside_radii <- c()
            for (j in 1:num_rings) {
                    r_inner <- spacing + 2 * (j - 1) * (width + spacing)
                    r_outer <- r_inner + 2 * width
                    inside_radii <- c(inside_radii, c(r_inner, r_outer))
                }

            for (i in 1:nrow(grid)) {
                point <- grid[i, ]
                is_inside <- FALSE

                if (length(inside_radii) >= 2) {
                    for (j in seq(1, length(inside_radii), by = 2)) {
                        r_inner <- inside_radii[j]
                        r_outer <- inside_radii[j + 1]
                        max_value <- max(abs(point$x), abs(point$y))
                        if (max_value >= r_inner && max_value <= r_outer) {
                            is_inside <- TRUE
                            break
                        }
                    }
                }
                
                inside_true <- c(inside_true, is_inside)
            }

            eps <- 1e-16

            P_raw <- inside_true + eps
            Q_raw <- inside_est + eps

            P <- P_raw / sum(P_raw)
            Q <- Q_raw / sum(Q_raw)

            # add small epsilon to avoid log(0)

          

            KL <- sum(P * log(P / Q))
            return(KL)
        }

        if (shape == "concentric") {
            inside_est <- points_in_sf(grid, est_sf)
            inside_true <- c()
            # for grid points, check if in any of the concentric rings using the math
            for (i in 1:nrow(grid)) {
                point <- grid[i, ]
                is_inside <- FALSE

                r_inner <- spacing

                for (j in 1:num_rings) {
                    r_outer <- r_inner + width
                    r_point <- sqrt(point$x^2 + point$y^2)
                    if (r_point >= r_inner && r_point <= r_outer) {
                        is_inside <- TRUE
                        break
                    }
                    r_inner <- r_outer + spacing
                }
                inside_true <- c(inside_true, is_inside)
            }

            eps <- 1e-16

            P_raw <- inside_true + eps
            Q_raw <- inside_est + eps

            P <- P_raw / sum(P_raw)
            Q <- Q_raw / sum(Q_raw)

            KL <- sum(P * log(P / Q))
            return(KL)
        }

        # true vs estimated occupancy
        inside_true <- points_in_sf(grid, true_sf)
        inside_est <- points_in_sf(grid, est_sf)

        P <- inside_true / sum(inside_true)
        Q <- inside_est / sum(inside_est)

        # add small epsilon to avoid log(0)
        eps <- 1e-16
        P <- P + eps
        Q <- Q + eps

        KL <- sum(P * log(P / Q))
        return(KL)
    } else {
        # Placeholder for 3D KL divergence calculation
        return(NA)
    }
}

#' IOU_score
#' 
#' Calculate Intersection over Union (IoU) score between true shape and estimated alpha shape
#' @param true_shape The true shape representation (e.g., polygon for 2D).
#' @param estimated_shape The estimated alpha shape representation.
#' @param grid A data frame of grid points with columns 'x' and 'y'.
#' @param dim The dimension (2 or 3).
#' @param grid_res Resolution of the grid for Monte Carlo integration.
#' @param shape The shape type ("concentric" or "squares").
#' @param num_rings Number of rings (for ring shapes).
#' @param spacing Spacing between rings.
#' @param width Width of each ring.
#' 
#' @return The IoU score.
IOU_score <- function(true_shape, estimated_shape, grid, dim = 2, grid_res = 500, shape = "concentric", num_rings = 3, spacing = 3.0, width = 1.5) {
    # Function to calculate Intersection over Union (IoU) score between true shape and estimated alpha shape
    # true_shape: the true shape representation (e.g., polygon for 2D)
    # estimated_shape: the estimated alpha shape representation
    # dim: dimension (2 or 3)

    if (dim == 2) {
        if (length(true_shape$polygons) == 0 || length(estimated_shape$polygons) == 0) {
            return(0)
        }

        # bounding box
        all_points <- rbind(do.call(rbind, true_shape$polygons), do.call(rbind, estimated_shape$polygons))
        
        true_sf <- to_sf_poly(true_shape$polygons)
        est_sf <- to_sf_poly(estimated_shape$polygons)

        # helper to check if points are in any polygon
        points_in_sf <- function(points_df, sf_polygons) {
            library(sf)
            pts <- st_as_sf(points_df, coords = c("x", "y"), crs = NA)
            result <- st_intersects(pts, sf_polygons, sparse = FALSE)
            apply(result, 1, any)
        }


        if (shape == "squares") {
            inside_est <- points_in_sf(grid, est_sf)
            inside_true <- c()
            # for grid points, check if in any of the square rings using the math
            for (i in 1:nrow(grid)) {
                point <- grid[i, ]
                is_inside <- FALSE

                for (j in 1:num_rings) {
                    r_inner <- spacing + 2 * (j - 1) * (width + spacing)
                    r_outer <- r_inner + 2 * width
                    max_value <- max(abs(point$x), abs(point$y))
                    if (max_value >= r_inner && max_value <= r_outer) {
                        is_inside <- TRUE
                        break
                    }
                }
                inside_true <- c(inside_true, is_inside)
            }

            intersection <- sum(inside_true & inside_est)
            union <- sum(inside_true | inside_est)

            IoU <- intersection / union
            return(IoU)
        }

        if (shape=="concentric") {
            inside_est <- points_in_sf(grid, est_sf)
            inside_true <- c()
            # for grid points, check if in any of the concentric rings using the math
            for (i in 1:nrow(grid)) {
                point <- grid[i, ]
                is_inside <- FALSE

                r_inner <- spacing

                for (j in 1:num_rings) {
                    r_outer <- r_inner + width
                    r_point <- sqrt(point$x^2 + point$y^2)
                    if (r_point >= r_inner && r_point <= r_outer) {
                        is_inside <- TRUE
                        break
                    }
                    r_inner <- r_outer + spacing
                }
                inside_true <- c(inside_true, is_inside)
            }

            intersection <- sum(inside_true & inside_est)
            union <- sum(inside_true | inside_est)

            IoU <- intersection / union
            return(IoU)
        }

        # true vs estimated occupancy
        inside_true <- points_in_sf(grid, true_sf)
        inside_est <- points_in_sf(grid, est_sf)

        intersection <- sum(inside_true & inside_est)
        union <- sum(inside_true | inside_est)

        IoU <- intersection / union
        return(IoU)
    } else {
        # Placeholder for 3D IoU calculation
        return(0)
    }
}

#' IoU_and_KL
#' 
#' Calculate both IoU score and KL divergence between true shape and estimated alpha shape
#' @param true_shape The true shape representation (e.g., polygon for 2D).
#' @param estimated_shape The estimated alpha shape representation.
#' @param grid A data frame of grid points with columns 'x' and 'y'.
#' @param dim The dimension (2 or 3).
#' @param shape The shape type ("concentric" or "squares").
#' @param num_rings Number of rings (for ring shapes).
#' @param spacing Spacing between rings.
#' @param width Width of each ring.
#' 
#' @return A list containing IoU score and KL divergence.
IoU_and_KL <- function(true_shape, estimated_shape, grid,
                       dim = 2,
                       shape = "concentric",
                       num_rings = 3,
                       spacing = 3.0,
                       width = 1.5) {

    if (dim != 2) return(list(IOU_score = NA, KL_divergence = NA))

    if (length(true_shape$polygons) == 0 || length(estimated_shape$polygons) == 0) {
        return(list(IOU_score = 0, KL_divergence = NA))
    }

    library(sf)

    # Convert polygons to sf
    true_sf <- to_sf_poly(true_shape$polygons)
    est_sf <- to_sf_poly(estimated_shape$polygons)

    # Vectorized polygon membership
    points_in_sf <- function(points_df, sf_polygons) {
        pts <- st_as_sf(points_df, coords = c("x", "y"), crs = NA)
        result <- st_intersects(pts, sf_polygons, sparse = FALSE)
        apply(result, 1, any)
    }

    n <- nrow(grid)

    if (shape == "concentric") {

        r_point <- sqrt(grid$x^2 + grid$y^2)

        # inner radii start, outer radii end
        r_inner <- spacing + (0:(num_rings-1)) * (width + spacing)
        r_outer <- r_inner + width

        # Broadcasting: (n grid points) x (num_rings)
        mat_inner <- matrix(r_inner, nrow = n, ncol = num_rings, byrow = TRUE)
        mat_outer <- matrix(r_outer, nrow = n, ncol = num_rings, byrow = TRUE)
        rpt_mat <- matrix(r_point, nrow = n, ncol = num_rings)

        inside_true <- apply(rpt_mat >= mat_inner & rpt_mat <= mat_outer, 1, any)

    } else if (shape == "squares") {

        max_val <- pmax(abs(grid$x), abs(grid$y))

        r_inner <- spacing + 2 * (0:(num_rings-1)) * (width + spacing)
        r_outer <- r_inner + 2 * width

        mat_inner <- matrix(r_inner, nrow = n, ncol = num_rings, byrow = TRUE)
        mat_outer <- matrix(r_outer, nrow = n, ncol = num_rings, byrow = TRUE)
        mv_mat <- matrix(max_val, nrow = n, ncol = num_rings)

        inside_true <- apply(mv_mat >= mat_inner & mv_mat <= mat_outer, 1, any)

    } else {
        # general polygons
        inside_true <- points_in_sf(grid, true_sf)
    }

    inside_est <- points_in_sf(grid, est_sf)


    intersection <- sum(inside_true & inside_est)
    union <- sum(inside_true | inside_est)
    IoU <- intersection / union


    eps <- 1e-16

    P_raw <- inside_true + eps
    Q_raw <- inside_est  + eps

    P <- P_raw / sum(P_raw)
    Q <- Q_raw / sum(Q_raw)

    KL_left <- sum(P * log(P / Q))
    KL_right <- sum(Q * log(Q / P))

    KL <- 0.5*(KL_left + KL_right)

    # Wasserstein distance (Earth Mover's Distance) approximation

    

    return(list(IOU_score = IoU, KL_divergence = KL))
}

#' main_test
#' 
#' Main test function for alpha shape visualization and metrics calculation
#' @param dim Dimension (2 or 3).
#' @param alpha Alpha parameter for alpha shape.
#' @param CL Confidence level for alpha shape.
#' @param num_points Number of points to generate.
#' @param noise_in_data Proportion of noise points to add.
#' @param title Title for the visualization.
#' @param r_inner Inner radius (for ring shapes).
#' @param r_outer Outer radius (for ring shapes).
#' @param shape Shape type ("ring", "concentric", "pitchfork", "squares").
#' @param spacing Spacing between rings or squares.
#' @param ring_width Width of each ring or square.
#' @param num_rings Number of rings or squares.
#' @param sampling_fraction Fraction of points to sample for alpha shape computation.
#' @param save_part Optional string to append to the saved filename.
#' 
#' @return None. Generates visualization and prints metrics.
main_test <- function(dim = 2, alpha = 0.2, CL = 1.0, num_points = 1000, noise_in_data = 0.1, title = "Alpha Shape Visualization", r_inner = 1, r_outer = 3, shape = "ring", spacing = 1, ring_width = 0.5, num_rings = 3, sampling_fraction = 0.8, save_part = NULL) {
    # generate test data

    if (dim == 2) {
        print("Generating 2D test data...")
        if (shape == "ring") {
            points <- ring_test_generation_2D(num_points = num_points, noise = noise_in_data, r_outer = r_outer, r_inner = r_inner)
            true_area <- pi * (r_outer^2 - r_inner^2)
            true_shape <- list(polygons = list(
                matrix(c(
                    r_outer * cos(seq(0, 2 * pi, length.out = 100)),
                    r_outer * sin(seq(0, 2 * pi, length.out = 100))
                ), ncol = 2),
                matrix(c(
                    r_inner * cos(seq(0, 2 * pi, length.out = 100)),
                    r_inner * sin(seq(0, 2 * pi, length.out = 100))
                ), ncol = 2)
            ))
        } else if (shape == "concentric") {
            points <- concentric_rings_2D(
                num_rings = num_rings, spacing = spacing,
                ring_width = ring_width, total_points = num_points, noise = noise_in_data
            )

            # Build true polygons with proper holes
            true_shape <- list(polygons = list())
            r_inner <- spacing
            for (i in 1:num_rings) {
                r_outer <- r_inner + ring_width
                # Outer polygon (counterclockwise)
                true_shape$polygons[[length(true_shape$polygons) + 1]] <- cbind(
                    r_outer * cos(seq(0, 2 * pi, length.out = 200)),
                    r_outer * sin(seq(0, 2 * pi, length.out = 200))
                )
                # Inner polygon (clockwise)
                true_shape$polygons[[length(true_shape$polygons) + 1]] <- cbind(
                    r_inner * cos(seq(0, -2 * pi, length.out = 200)),
                    r_inner * sin(seq(0, -2 * pi, length.out = 200))
                )
                r_inner <- r_outer + spacing
            }

            # True area sum
            true_area <- area_concentric_rings_2D(
                num_rings = num_rings,
                spacing = spacing, ring_width = ring_width
            )
        } else if (shape == "pitchfork") {
            points <- pitchfork_bifurcation_datacloud(num_points = num_points, noise = noise_in_data, direction = "up")
            true_shape <- list(polygons = list(
                polygon_pitchfork(x_min = -2, x_max = 2, r_min = -2, r_max = 2)
            ))
            true_area <- 8
        } else if (shape == "squares") {
            points <- alternating_squares_pointcloud_2D(total_points = num_points, num_squares = num_rings, square_width = ring_width, spacing = spacing, noise = noise_in_data)
            true_area <- square_area_2D(num_squares = num_rings, square_width = ring_width, spacing = spacing)
            true_shape <- list(polygons = alternating_squares_polygon(num_squares = num_rings, square_width = ring_width, spacing = spacing))
        } else {
            stop("Unknown shape type.")
        }

        print("Computing alpha shape...")

        alpha_shape_result <- alpha_shape_2D(points, alpha, CL, sampling_fraction = sampling_fraction)

        print("Calculating metrics...")

        # build a filename safely using base R; "%+%" is not a base R operator and
        # will error if the package that defines it isn't loaded
        if (!is.null(save_part)) {
            save_path <- paste0("alpha_shape_2D_alpha_", alpha, "_CL_", CL, "_noise_", noise_in_data, "_shape_", shape, "_numpoints_", num_points, "_", save_part, ".png")
        } else {
            save_path <- paste0("alpha_shape_2D_alpha_", alpha, "_CL_", CL, "_noise_", noise_in_data, "_shape_", shape, "_numpoints_", num_points, ".png")
        }
        visualize_alpha_shape_2D(alpha_shape_result, points, title, r_outer = r_outer, r_inner = r_inner, shape = shape, spacing = spacing, ring_width = ring_width, num_rings = num_rings, save_path = save_path)

        # estimated_area <- calculate_area_volume(alpha_shape_result, dim = 2)
        # error <- error_in_2D(true_area, estimated_area)
        # print(paste("Error in 2D area estimation:", error))
        # print(paste("Estimated Area:", estimated_area, "True Area:", true_area))

        # create grid for KL and IOU calculations
        all_points <- rbind(do.call(rbind, true_shape$polygons), do.call(rbind, alpha_shape_result$polygons))
        x_range <- range(all_points[, 1])
        y_range <- range(all_points[, 2])

        grid <- data.frame(x = runif(100000, min = x_range[1], max = x_range[2]), y = runif(100000, min = y_range[1], max = y_range[2]))

        
        # KL_divergence <- calculate_KL_divergence(
        #     true_shape = true_shape,
        #     estimated_shape = alpha_shape_result,
        #     grid = grid,
        #     dim = 2,
        #     shape = shape,
        #     num_rings = num_rings,
        #     spacing = spacing,
        #     width = ring_width
        # )
        # IOU_score_value <- IOU_score(
        #     true_shape = true_shape,
        #     estimated_shape = alpha_shape_result,
        #     grid = grid,
        #     dim = 2,
        #     shape = shape,
        #     num_rings = num_rings,
        #     spacing = spacing,
        #     width = ring_width
        # )
        

        metrics <- IoU_and_KL(
            true_shape = true_shape,
            estimated_shape = alpha_shape_result,
            grid = grid,
            dim = 2,
            shape = shape,
            num_rings = num_rings,
            spacing = spacing,
            width = ring_width
        )

        print(paste("KL Divergence:", metrics$KL_divergence, "IOU Score:", metrics$IOU_score))
    } else {
        points <- torus_test_generation_3D(num_points = num_points, noise = noise_in_data)
        print("3D alpha shape not yet implemented.")
    }

    # Return both values as a list so callers can unpack by name
    return(list(KL_divergence = metrics$KL_divergence, IOU_score = metrics$IOU_score, alpha_shape = alpha_shape_result))
}

#' alpha_variations
#' 
#' Run alpha shape tests over multiple parameter combinations
#' @param alpha_value_list List of alpha values to test.
#' @param dim Dimension (2 or 3).
#' @param CL_list List of confidence levels to test.
#' @param num_points Number of points to generate.
#' @param noise_in_data List of noise proportions to test.
#' @param r_outer Outer radius (for ring shapes).
#' @param r_inner Inner radius (for ring shapes).
#' @param shape Shape type ("ring", "concentric", "pitchfork", "squares").
#' @param spacing Spacing between rings or squares.
#' @param ring_width Width of each ring or square.
#' @param num_rings Number of rings or squares.
#' @param specific_parameters Optional list of specific parameter combinations to run.
#' @param save_path Optional path to save the results CSV.
#' @param save_part Optional string to append to the saved filename.
#'  
#' @return None. Saves results to CSV.
alpha_variations <- function(alpha_value_list = c(0.1, 0.2),
                             dim = 2,
                             CL_list = c(1.0, 0.95, 0.90),
                             num_points = 10000,
                             noise_in_data = c(0.1, 0.01),
                             r_outer = 5,
                             r_inner = 3,
                             shape = "concentric",
                             spacing = 3.0,
                             ring_width = 1.5,
                             num_rings = 3,
                             specific_parameters = NULL,
                             save_path = NULL,
                             save_part = NULL) {
    total_trials <- length(alpha_value_list) * length(CL_list) * length(noise_in_data)
    trial_count <- 0
    KL_divergence_results <- data.frame()
    IOU_score_results <- data.frame()

    score_results <- data.frame(
        alpha = numeric(),
        CL = numeric(),
        noise = numeric(),
        KL_divergence = numeric(),
        IOU_score = numeric(),
        shape = character(),
        time_taken = numeric()
    )

    if (!is.null(specific_parameters)) {
        print("Running specific parameter combinations.")
        total_trials <- length(specific_parameters)
        for (params in specific_parameters) {
            # params may come in as a character vector if constructed with c(..., "shape")
            # coerce the numeric entries explicitly to avoid non-numeric errors
            alpha <- as.numeric(params[1])
            CL <- as.numeric(params[2])
            noise <- as.numeric(params[3])
            shape_specific <- as.character(params[4])
            time_start <- Sys.time()
            print(paste("Trial", trial_count + 1, "of", total_trials))
            print(paste("Running for Alpha:", alpha, "CL:", CL, "Noise:", noise, "Shape:", shape_specific))
            title <- paste("Alpha Shape Visualization - Alpha:", alpha, "CL:", CL, "Noise:", noise, "Shape:", shape_specific)
            res <- main_test(dim = dim, alpha = alpha, CL = CL, num_points = num_points, noise_in_data = noise, title = title, r_outer = r_outer, r_inner = r_inner, shape = shape_specific, spacing = spacing, ring_width = ring_width, num_rings = num_rings, save_part = save_part)
            KL_divergence <- res$KL_divergence
            IOU_score_value <- res$IOU_score
            KL_divergence_results <- rbind(KL_divergence_results, data.frame(alpha = alpha, CL = CL, noise = noise, KL_divergence = KL_divergence))
            IOU_score_results <- rbind(IOU_score_results, data.frame(alpha = alpha, CL = CL, noise = noise, IOU_score = IOU_score_value))
            trial_count <- trial_count + 1
            time_end <- Sys.time()
            time_taken <- as.numeric(time_end - time_start, units = "secs")
            score_results <- rbind(score_results, data.frame(alpha = alpha, CL = CL, noise = noise, KL_divergence = KL_divergence, IOU_score = IOU_score_value, shape = shape_specific, time_taken = time_taken))
        }
    }
    if (is.null(specific_parameters)) {
        for (alpha in alpha_value_list) {
            for (CL in CL_list) {
                for (noise in noise_in_data) {
                    print(paste("Trial", trial_count + 1, "of", total_trials))
                    print(paste("Running for Alpha:", alpha, "CL:", CL, "Noise:", noise))
                    title <- paste("Alpha Shape Visualization - Alpha:", alpha, "CL:", CL, "Noise:", noise)
                    res <- main_test(dim = dim, alpha = alpha, CL = CL, num_points = num_points, noise_in_data = noise, title = title, r_outer = r_outer, r_inner = r_inner, shape = shape, spacing = spacing, ring_width = ring_width, num_rings = num_rings, save_part = save_part)
                    KL_divergence <- res$KL_divergence
                    IOU_score_value <- res$IOU_score
                    KL_divergence_results <- rbind(KL_divergence_results, data.frame(alpha = alpha, CL = CL, noise = noise, KL_divergence = KL_divergence))
                    IOU_score_results <- rbind(IOU_score_results, data.frame(alpha = alpha, CL = CL, noise = noise, IOU_score = IOU_score_value))
                    score_results <- rbind(score_results, data.frame(alpha = alpha, CL = CL, noise = noise, KL_divergence = KL_divergence, IOU_score = IOU_score_value))
                    trial_count <- trial_count + 1
                }
            }
        }
    }

    if (!is.null(save_path)) {
        save_path <- save_path
    } else {
        if (shape == "concentric") {
            temp <- paste0(str(spacing), "_", str(ring_width), "_", str(num_rings))
            save_path <- paste0("alpha_shape_KL_IOU_concentric_rings_", temp, "_results.csv")
        } else if (shape == "ring") {
            temp <- paste0(str(r_inner), "_", str(r_outer))
            save_path <- paste0("alpha_shape_KL_IOU_ring_", temp, "_results.csv")
        } else if (shape == "pitchfork") {
            save_path <- paste0("alpha_shape_KL_IOU_pitchfork_results.csv")
        } else if (shape == "squares") {
            temp <- paste0(str(ring_width), "_", str(spacing), "_", str(num_rings))
            save_path <- paste0("alpha_shape_KL_IOU_squares_", temp, "_results.csv")
        } else {
            save_path <- paste0("alpha_shape_KL_IOU_", shape, "_results.csv")
        }
    }

    # save results to csv
    write.csv(score_results, file = save_path, row.names = FALSE)
    print(paste("Results saved to", save_path))
    print("Alpha variations completed.")
}

#' datapoints_variations
#' 
#' Run alpha shape tests varying number of data points
#' @param num_points_list List of number of points to test. 
#' 
#' @return None. Saves timing results to CSV.
datapoints_variations <- function(
    num_points_list = c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6)) {
    time_taken_list <- c()
    for (num_points in num_points_list) {
        time_start <- Sys.time()
        print(paste("Running for Number of Points:", num_points))
        title <- paste("Alpha Shape Visualization - Number of Points:", num_points)
        main_test(dim = 2, alpha = 1.0, CL = 0.95, shape = "ring", num_points = num_points, noise_in_data = 0.01, title = title)
        time_end <- Sys.time()
        time_taken <- time_end - time_start
        time_taken_list <- c(time_taken_list, as.numeric(time_taken, units = "secs"))
        print(paste("Time taken (seconds):", as.numeric(time_taken, units = "secs")))
    }
    result <- data.frame(
        num_points = num_points_list,
        time_taken_seconds = time_taken_list
    )
    # save to csv
    write.csv(result, file = "alpha_shape_datapoints_timing.csv", row.names = FALSE)
    print("Data points variations completed.")
}

#' plot_datapoints_timing
#' 
#' Plot timing results for varying number of data points
#' @param csv_file CSV file containing timing results.
#' @param scale Scale type for the plot (1: linear-linear, 2: linear-log, 3: log-log, 4: n vs sqrt t, 5: log n vs log sqrt t).
#' 
#' @return None. Saves plot to PNG file.
plot_datapoints_timing <- function(csv_file = "alpha_shape_datapoints_timing.csv", scale = 1) {
    data <- read.csv(csv_file)

    if (scale == 1) {
        save_path <- "alpha_shape_timing_linear_linear.png"
        png(save_path)
        plot(data$num_points, data$time_taken_seconds,
            type = "b",
            xlab = "Number of Points", ylab = "Time Taken (seconds)",
            main = "Alpha Shape Computation Time vs Number of Points"
        )
        dev.off()
    } else if (scale == 2) {
        save_path <- "alpha_shape_timing_linear_log.png"
        png(save_path)
        plot(data$num_points, data$time_taken_seconds,
            log = "y", type = "b",
            xlab = "Number of Points", ylab = "Time Taken (seconds, log scale)",
            main = "Alpha Shape Computation Time vs Number of Points"
        )
        dev.off()
    } else if (scale == 3) {
        save_path <- "alpha_shape_timing_log_log.png"
        png(save_path)
        plot(data$num_points, data$time_taken_seconds,
            log = "xy", type = "b",
            xlab = "Number of Points (log scale)", ylab = "Time Taken (seconds, log scale)",
            main = "Alpha Shape Computation Time vs Number of Points"
        )
        dev.off()
    } else if (scale == 4) {
        save_path <- "alpha_shape_timing_n_vs_sqrt_t.png"
        png(save_path)
        plot(data$num_points, sqrt(data$time_taken_seconds),
            type = "b",
            xlab = "Number of Points", ylab = "Square Root of Time Taken (seconds)",
            # main = "Alpha Shape Computation Time vs Number of Points (n vs sqrt t scale)"
        )

        dev.off()
    } else if (scale == 5) { # log n vs log sqrt t
        save_path <- "alpha_shape_timing_log_n_vs_log_sqrt_t.png"
        png(save_path)
        plot(log(data$num_points), log(sqrt(data$time_taken_seconds)),
            type = "b",
            xlab = "Log(Number of Points)", ylab = "Log(Square Root of Time Taken (seconds))",
            main = "Alpha Shape Computation Time vs Number of Points (log-log sqrt scale)"
        )
        dev.off()
    }
}

#' sample_fraction_variations
#' 
#' Run alpha shape tests varying sampling fraction
#' @param fraction_list List of sampling fractions to test.
#' @param dim Dimension (2 or 3).
#' @param alpha Alpha parameter for alpha shape.
#' @param CL Confidence level for alpha shape.
#' @param num_points Number of points to generate.
#' @param noise_in_data Proportion of noise points to add.
#' @param shape Shape type ("ring", "concentric", "pitchfork", "squares").
#' @param r_outer Outer radius (for ring shapes).
#' @param r_inner Inner radius (for ring shapes).
#' @param num_rings Number of rings or squares.
#' @param spacing Spacing between rings or squares.
#' @param ring_width Width of each ring or square.
#' 
#' @return None. Saves results to CSV.
sample_fraction_variations <- function(
    fraction_list = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99),
    dim = 2,
    alpha = 1.0,
    CL = 0.95,
    num_points = 10000,
    noise_in_data = 0.01,
    shape = "ring",
    r_outer = 5,
    r_inner = 3,
    num_rings = 2,
    spacing = 3.0,
    ring_width = 1.5) {
    data_results <- data.frame()
    for (fraction in fraction_list) {
        result <- main_test(dim = dim, alpha = alpha, CL = CL, num_points = num_points, noise_in_data = noise_in_data, shape = shape, spacing = spacing, ring_width = ring_width, num_rings = num_rings, r_outer = r_outer, r_inner = r_inner, sampling_fraction = fraction, title = paste("Alpha Shape Visualization - Sampling Fraction:", fraction))
        print(paste("Completed for Sampling Fraction:", fraction))
        # save to csv
        data_results <- rbind(data_results, data.frame(sampling_fraction = fraction, KL_divergence = result$KL_divergence, IOU_score = result$IOU_score))
    }

    write.csv(data_results, file = "alpha_shape_sampling_fraction_results.csv", row.names = FALSE)
}

#' corner_variations
#' 
#' Run alpha shape tests varying corner shape parameters
#' @param s Spacing between squares.
#' @param w Width of each square.
#' @param num_rings Number of rings or squares.
#' @param alpha_list List of alpha parameters for alpha shape.
#' @param CL_list List of confidence levels for alpha shape.
#' @param num_points Number of points to generate.
#' @param noise_in_data_list List of proportions of noise points to add.
#' 
#' @return None. Saves results to CSV.
corner_variations <- function(s = 3.0, w = 1.5, num_rings = 3, alpha_list = c(1.0), CL_list = c(0.95), num_points = 10000, noise_in_data_list = c(0.01)) {
    result_dataframe <- data.frame()

    for (alpha in alpha_list) {
        for (CL in CL_list) {
            for (noise_in_data in noise_in_data_list) {
                print(paste("Running for Alpha:", alpha, "CL:", CL, "Noise:", noise_in_data))
                title <- paste("Alpha Shape Visualization - Alpha:", alpha, "CL:", CL, "Noise:", noise_in_data)
                result <- main_test(dim = 2, alpha = alpha, CL = CL, num_points = num_points, noise_in_data = noise_in_data, title = title, shape = "squares", spacing = s, ring_width = w, num_rings = num_rings)
                KL_divergence <- result$KL_divergence
                IOU_score_value <- result$IOU_score
                result_dataframe <- rbind(result_dataframe, data.frame(alpha = alpha, CL = CL, noise = noise_in_data, KL_divergence = KL_divergence, IOU_score = IOU_score_value))
            }
        }
    }

    write.csv(result_dataframe, file = "corner_variations_results.csv", row.names = FALSE)
}