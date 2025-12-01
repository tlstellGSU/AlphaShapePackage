#' Alpha Shape Utilities
#'
#' @name alpha_shape
#' @rdname alpha_shape
#' @importFrom dbscan dbscan kNNdist
#' @importFrom stats dist
NULL

#' Circumcircle
#'
#' Function to calculate the circumcircle of a triangle defined by points p1, p2, p3
#' @param p1 A vector of length 2 (x, y) coordinates of the first point.
#' @param p2 A vector of length 2 (x, y) coordinates of the second point.
#' @param p3 A vector of length 2 (x, y) coordinates of the third point.
#'
#' @return A list with elements 'center' (vector of length 2) and 'r' (radius).
circumcircle <- function(p1, p2, p3) {
    A <- p2 - p1
    B <- p3 - p1
    d <- 2 * (A[1] * B[2] - A[2] * B[1])

    if (abs(d) < 1e-16) {
        return(list(center = c(NA, NA), r = Inf))
    }

    ux <- ((B[2] * (A[1]^2 + A[2]^2) - A[2] * (B[1]^2 + B[2]^2)) / d) + p1[1]
    uy <- ((A[1] * (B[1]^2 + B[2]^2) - B[1] * (A[1]^2 + A[2]^2)) / d) + p1[2]
    center <- c(ux, uy)
    r <- sqrt((center[1] - p1[1])^2 + (center[2] - p1[2])^2)
    return(list(center = center, r = r))
}

#' point_in_circle
#'
#' Function to check if a point is inside a circle
#' @param p A vector of length 2 (x, y) coordinates of the point.
#' @param center A vector of length 2 (x, y) coordinates of the circle center.
#' @param r Radius of the circle.
#'
#' @return TRUE if the point is inside the circle, FALSE otherwise.
point_in_circle <- function(p, center, r) {
    # Function to check if a point is inside a circle
    # p: a vector of length 2 (x, y) coordinates of the point
    # center: a vector of length 2 (x, y) coordinates of the circle center
    # r: radius of the circle
    return(sqrt((p[1] - center[1])^2 + (p[2] - center[2])^2) < r)
}

#' circumsphere
#'
#' Function to calculate the circumsphere of a tetrahedron defined by points p1, p2, p3, p4
#' @param p1 A vector of length 3 (x, y, z) coordinates of the first point.
#' @param p2 A vector of length 3 (x, y, z) coordinates of the second point.
#' @param p3 A vector of length 3 (x, y, z) coordinates of the third point.
#' @param p4 A vector of length 3 (x, y, z) coordinates of the fourth point.
#'
#' @return A list with elements 'center' (vector of length 3) and 'r' (radius).
circumsphere <- function(p1, p2, p3, p4) {
    # Function to calculate the circumsphere of a tetrahedron defined by points p1, p2, p3, p4
    # p1, p2, p3, p4: vectors of length 3 (x, y, z) coordinates of the points

    A <- rbind(p2 - p1, p3 - p1, p4 - p1)
    B <- c(sum((p2 - p1)^2), sum((p3 - p1)^2), sum((p4 - p1)^2)) / 2

    center <- solve(A, B) + p1
    r <- sqrt(sum((center - p1)^2))
    return(list(center = center, r = r))
}

#' point_in_sphere
#'
#' Function to check if a point is inside a sphere
#' @param p A vector of length 3 (x, y, z) coordinates of the point.
#' @param center A vector of length 3 (x, y, z) coordinates of the sphere center.
#' @param r Radius of the sphere.
#'
#' @return TRUE if the point is inside the sphere, FALSE otherwise.
point_in_sphere <- function(p, center, r) {
    # Function to check if a point is inside a sphere
    # p: a vector of length 3 (x, y, z) coordinates of the point
    # center: a vector of length 3 (x, y, z) coordinates of the sphere center
    # r: radius of the sphere
    return(sqrt(sum((p - center)^2)) < r)
}

#' Delaunay triangulation
#'
#' Performs Delaunay triangulation for 2D points and Delaunay tetrahedralization for 3D points using the Bowyer-Watson algorithm.
#' @param points A numeric matrix of points (n x 2 for 2D, n x 3 for 3D).
#' @param dim Dimension of the points (2 or 3).
#' @return A matrix of indices representing the triangles (for 2D) or tetrahedra (for
#'  3D) formed by the Delaunay triangulation/tetrahedralization. Each row corresponds to a triangle (3 indices) or tetrahedron (4 indices).
deluanay_triangulation <- function(points, dim = 2) {
    if (dim != 2 && dim != 3) {
        stop("Dimension must be 2 or 3.")
    }
    if (dim == 2) {
        points <- points + matrix(rnorm(length(points), sd = 1e-8), ncol = 2)
        points <- as.matrix(points)
        if (ncol(points) != 2) {
            stop("Input points must have exactly two columns for x and y coordinates.")
        }

        min_x <- min(points[, 1])
        max_x <- max(points[, 1])
        min_y <- min(points[, 2])
        max_y <- max(points[, 2])

        dx <- max_x - min_x
        dy <- max_y - min_y

        delta_max <- max(dx, dy)
        mid_x <- (min_x + max_x) / 2
        mid_y <- (min_y + max_y) / 2

        super <- rbind(
            c(mid_x - 20 * delta_max, mid_y - delta_max),
            c(mid_x, mid_y + 20 * delta_max),
            c(mid_x + 20 * delta_max, mid_y - delta_max)
        )

        triangles <- matrix(NA, nrow = 1, ncol = 3)
        points_all <- rbind(points, super)

        n <- nrow(points)
        super_idx <- (n + 1):(n + 3)
        triangles[1, ] <- super_idx

        for (i in 1:n) {
            point <- points_all[i, ]
            bad_triangles <- c()
            circumcenters <- list()
            circumradii <- c()

            for (j in 1:nrow(triangles)) {
                tri <- triangles[j, ]
                p1 <- points_all[tri[1], ]
                p2 <- points_all[tri[2], ]
                p3 <- points_all[tri[3], ]

                cc <- circumcircle(p1, p2, p3)
                circumcenters[[j]] <- cc$center
                circumradii[j] <- cc$r

                if (point_in_circle(point, cc$center, cc$r)) {
                    bad_triangles <- c(bad_triangles, j)
                }
            }

            edge_count <- list()
            for (bt in bad_triangles) {
                tri <- triangles[bt, ]
                edges <- list(
                    paste(sort(c(tri[1], tri[2])), collapse = "-"),
                    paste(sort(c(tri[2], tri[3])), collapse = "-"),
                    paste(sort(c(tri[3], tri[1])), collapse = "-")
                )
                for (edge in edges) {
                    if (is.null(edge_count[[edge]])) {
                        edge_count[[edge]] <- 1
                    } else {
                        edge_count[[edge]] <- edge_count[[edge]] + 1
                    }
                }
            }

            triangles <- triangles[-bad_triangles, , drop = FALSE]

            for (edge in names(edge_count)) {
                if (edge_count[[edge]] == 1) {
                    verts <- as.numeric(unlist(strsplit(edge, "-")))
                    new_tri <- c(verts, i)
                    triangles <- rbind(triangles, new_tri)
                }
            }
        }

        triangles <- triangles[!apply(triangles, 1, function(tri) any(tri %in% super_idx)), , drop = FALSE]
        return(triangles)
    } else {
        points <- as.matrix(points)
        if (ncol(points) != 3) {
            stop("Input points must have exactly three columns for x, y, and z coordinates.")
        }

        n <- nrow(points)
        min_x <- min(points[, 1])
        max_x <- max(points[, 1])
        min_y <- min(points[, 2])
        max_y <- max(points[, 2])
        min_z <- min(points[, 3])
        max_z <- max(points[, 3])

        dx <- max_x - min_x
        dy <- max_y - min_y
        dz <- max_z - min_z
        delta_max <- max(dx, dy, dz)
        mid_x <- (min_x + max_x) / 2
        mid_y <- (min_y + max_y) / 2
        mid_z <- (min_z + max_z) / 2

        super <- rbind(
            c(mid_x - 20 * delta_max, mid_y - 20 * delta_max, mid_z - 20 * delta_max),
            c(mid_x + 20 * delta_max, mid_y + 20 * delta_max, mid_z - 20 * delta_max),
            c(mid_x + 20 * delta_max, mid_y - 20 * delta_max, mid_z + 20 * delta_max),
            c(mid_x - 20 * delta_max, mid_y + 20 * delta_max, mid_z + 20 * delta_max)
        )

        points_all <- rbind(points, super)
        super_idx <- (n + 1):(n + 4)

        tetrahedra <- matrix(NA, nrow = 1, ncol = 4)
        tetrahedra[1, ] <- super_idx

        for (i in 1:n) {
            point <- points_all[i, ]
            bad_tetrahedra <- c()
            circumcenters <- list()
            circumradii <- c()

            for (j in 1:nrow(tetrahedra)) {
                tet <- tetrahedra[j, ]
                p1 <- points_all[tet[1], ]
                p2 <- points_all[tet[2], ]
                p3 <- points_all[tet[3], ]
                p4 <- points_all[tet[4], ]

                cs <- circumsphere(p1, p2, p3, p4)
                circumcenters[[j]] <- cs$center
                circumradii[j] <- cs$r

                if (point_in_sphere(point, cs$center, cs$r)) {
                    bad_tetrahedra <- c(bad_tetrahedra, j)
                }
            }

            face_count <- list()
            for (bt in bad_tetrahedra) {
                tet <- tetrahedra[bt, ]
                faces <- list(
                    paste(sort(c(tet[1], tet[2], tet[3])), collapse = "-"),
                    paste(sort(c(tet[1], tet[2], tet[4])), collapse = "-"),
                    paste(sort(c(tet[1], tet[3], tet[4])), collapse = "-"),
                    paste(sort(c(tet[2], tet[3], tet[4])), collapse = "-")
                )
                for (face in faces) {
                    if (is.null(face_count[[face]])) {
                        face_count[[face]] <- 1
                    } else {
                        face_count[[face]] <- face_count[[face]] + 1
                    }
                }
            }

            tetrahedra <- tetrahedra[-bad_tetrahedra, , drop = FALSE]

            for (face in names(face_count)) {
                if (face_count[[face]] == 1) {
                    verts <- as.numeric(unlist(strsplit(face, "-")))
                    new_tet <- c(verts, i)
                    tetrahedra <- rbind(tetrahedra, new_tet)
                }
            }
        }

        tetrahedra <- tetrahedra[!apply(tetrahedra, 1, function(tet) any(tet %in% super_idx)), , drop = FALSE]

        return(tetrahedra)
    }
}

#' polygon_area
#'
#' Calculate the area of a polygon given its vertices
#' @param polygon A matrix of vertices (n x 2) representing the polygon.
#' @return The area of the polygon.
polygon_area <- function(polygon) {
    if (!is.matrix(polygon) || ncol(polygon) != 2) {
        stop("Input polygon must be a matrix with two columns (x and y coordinates).")
    }

    # Calculate the area using the shoelace formula
    x <- polygon[, 1]
    y <- polygon[, 2]

    n <- length(x)
    if (n < 3) {
        return(0) # Not a polygon
    }

    area <- 0.5 * abs(sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1])))
    return(area)
}

#' edge_length
#'
#' Calculate the length of an edge given its two endpoints
#' @param p1 A numeric vector representing the first point (x, y).
#' @param p2 A numeric vector representing the second point (x, y).
#' @return The length of the edge.
edge_length <- function(p1, p2) {
    return(sqrt(sum((p1 - p2)^2)))
}

#' triangle_area
#'
#' Calculate the area of a triangle given its three vertices
#' @param p1 A numeric vector representing the first point (x, y).
#' @param p2 A numeric vector representing the second point (x, y).
#' @param p3 A numeric vector representing the third point (x, y).
#' @return The area of the triangle.
triangle_area <- function(p1, p2, p3) {
    a <- edge_length(p1, p2)
    b <- edge_length(p2, p3)
    c <- edge_length(p3, p1)
    s <- (a + b + c) / 2
    area <- sqrt(s * (s - a) * (s - b) * (s - c))
    return(area)
}

#' tetrahedron_volume
#'
#' Calculate the volume of a tetrahedron given its four vertices
#' @param p1 A numeric vector representing the first point (x, y, z).
#' @param p2 A numeric vector representing the second point (x, y, z).
#' @param p3 A numeric vector representing the third point (x, y, z).
#' @param p4 A numeric vector representing the fourth point (x, y, z).
#' @return The volume of the tetrahedron.
tetrahedron_volume <- function(p1, p2, p3, p4) {
    mat <- rbind(
        c(1, p1),
        c(1, p2),
        c(1, p3),
        c(1, p4)
    )
    vol <- abs(det(mat)) / 6
    return(vol)
}

#' tetrahedron_surface_area
#'
#' Calculate the surface area of a tetrahedron given its four vertices
#' @param p1 A numeric vector representing the first point (x, y, z).
#' @param p2 A numeric vector representing the second point (x, y, z).
#' @param p3 A numeric vector representing the third point (x, y, z).
#' @param p4 A numeric vector representing the fourth point (x, y, z).
#' @return The surface area of the tetrahedron.
tetrahedron_surface_area <- function(p1, p2, p3, p4) {
    area1 <- triangle_area(p1, p2, p3)
    area2 <- triangle_area(p1, p2, p4)
    area3 <- triangle_area(p1, p3, p4)
    area4 <- triangle_area(p2, p3, p4)
    total_area <- area1 + area2 + area3 + area4
    return(total_area)
}

#' edges_to_loops
#'
#' Convert boundary edges to ordered loops (polygons)
#' @param boundary_edges A list of numeric vectors length 2 (vertex indices).
#' @param points A matrix of point coordinates (n x 2).
#' @return A list of polygons, each polygon is a matrix (m x 2) of coordinates (closed: last == first).
edges_to_loops <- function(boundary_edges, points) {
    # boundary_edges: list of numeric vectors length 2 (vertex indices)
    # points: original points matrix (n x 2)
    # returns: list of polygons, each polygon is a matrix (m x 2) of coordinates (closed: last == first)
    if (length(boundary_edges) == 0) {
        return(list())
    }

    # map edge key to logical used flag
    edge_keys <- sapply(boundary_edges, function(e) {
        paste(sort(e), collapse = "-")
    }, USE.NAMES = FALSE)

    edge_used <- setNames(rep(FALSE, length(edge_keys)), edge_keys)

    # adjacency: vertex -> neighbor indices (list of integer vectors)
    adj <- list()
    for (e in boundary_edges) {
        a <- as.character(e[1])
        b <- as.character(e[2])
        adj[[a]] <- unique(c(adj[[a]], e[2]))
        adj[[b]] <- unique(c(adj[[b]], e[1]))
    }

    loops <- list()

    # function to mark edge used by its key
    mark_edge_used <- function(u, v) {
        k <- paste(sort(c(u, v)), collapse = "-")
        edge_used[k] <<- TRUE
    }
    is_edge_used <- function(u, v) {
        k <- paste(sort(c(u, v)), collapse = "-")
        edge_used[k]
    }

    # iterate until all edges used
    while (any(!edge_used)) {
        # find first unused edge key and its vertices
        unused_keys <- names(edge_used)[!edge_used]
        key0 <- unused_keys[1]
        verts <- as.numeric(unlist(strsplit(key0, "-")))
        start <- verts[1]
        next_v <- verts[2]

        # start walking the loop
        loop_idx <- c(start, next_v)
        mark_edge_used(start, next_v)
        prev_v <- start
        cur_v <- next_v

        # walk until we return to start
        safety <- 0
        while (cur_v != start && safety < 10000) {
            safety <- safety + 1
            neighs <- adj[[as.character(cur_v)]]
            # choose neighbor not equal to prev_v and where edge is not already used, if possible
            candidate <- NA
            for (nv in neighs) {
                if (nv == prev_v) next
                if (!is.na(nv) && !is_edge_used(cur_v, nv)) {
                    candidate <- nv
                    break
                }
            }
            # if all remaining neighbors lead to used edges (rare due to boundary structure), choose any neighbor that is not prev_v
            if (is.na(candidate)) {
                for (nv in neighs) {
                    if (nv != prev_v) {
                        candidate <- nv
                        break
                    }
                }
            }
            if (is.na(candidate)) break

            loop_idx <- c(loop_idx, candidate)
            mark_edge_used(cur_v, candidate)
            prev_v <- cur_v
            cur_v <- candidate
        }
        # close the loop by ensuring first == last
        if (loop_idx[1] != loop_idx[length(loop_idx)]) loop_idx <- c(loop_idx, loop_idx[1])
        # convert to coordinates matrix
        poly_coords <- points[loop_idx, , drop = FALSE]
        loops[[length(loops) + 1]] <- poly_coords
    }

    return(loops)
}

#' build_adj_my
#'
#' build adjacency list from edges
#' @param edges A list of numeric vectors length 2 (vertex indices).
#' @return A named list where each name is a vertex and each element is a vector of adjacent vertices.
build_adj_my <- function(edges) {
    edges <- do.call(rbind, edges)
    nodes <- sort(unique(as.vector(edges)))

    adj <- vector("list", length(nodes))
    names(adj) <- nodes

    for (i in seq_len(nrow(edges))) {
        u <- as.character(edges[i, 1])
        v <- as.character(edges[i, 2])
        adj[[u]] <- c(adj[[u]], v)
        adj[[v]] <- c(adj[[v]], u)
    }
    adj
}

#' find_cycles_my
#'
#' find cycles in undirected graph
#' @param adj A named list where each name is a vertex and each element is a
#' vector of adjacent vertices.
#' @return A list of cycles, each cycle is a vector of vertex names.
find_cycles_my <- function(adj) {
    nodes <- names(adj)
    visited <- setNames(rep(FALSE, length(nodes)), nodes)
    parent <- setNames(rep(NA_character_, length(nodes)), nodes)

    cycles <- list()

    dfs <- function(u) {
        visited[[u]] <<- TRUE

        for (v in adj[[u]]) {
            if (!visited[[v]]) {
                parent[[v]] <<- u
                dfs(v)
            } else if (!identical(parent[[u]], v)) {
                cu <- u
                cv <- v
                path <- cu

                # climb parent chain until reaching v
                while (!is.na(parent[[cu]]) && !identical(parent[[cu]], cv)) {
                    cu <- parent[[cu]]
                    path <- c(path, cu)
                }

                if (!is.na(parent[[cu]]) && identical(parent[[cu]], cv)) {
                    path <- c(path, cv)

                    # canonical orientation to avoid duplicates
                    rpath <- rev(path)
                    canon <- if (paste(path, collapse = ",") < paste(rpath, collapse = ",")) path else rpath

                    cycles[[length(cycles) + 1]] <<- canon
                }
            }
        }
    }

    for (n in nodes) {
        if (!visited[[n]]) dfs(n)
    }

    unique(lapply(cycles, function(x) unique(x)))
}

#' auto_epsilon
#'
#' Estimate epsilon for DBSCAN using k-NN distance distribution
#' @param points A numeric matrix of points (n x d).
#' @param min_points The number of nearest neighbors to consider (k).
#' @param quantile_value The quantile of the k-NN distances to use as epsilon.
#' @return A numeric value representing the estimated epsilon.
#' @export
auto_epsilon <- function(points, min_points = 5, quantile_value = 0.90) {
    points <- as.matrix(points)
    k_distribution <- dbscan::kNNdist(points, k = min_points)
    k_distribution_vec <- as.numeric(k_distribution)
    epsilon_estimated <- as.numeric(quantile(k_distribution_vec, probs = quantile_value))

    return(epsilon_estimated)
}

#' fix_orientation
#'
#' Fix polygon orientation
#' @param mat A matrix of polygon vertices (n x 2).
#' @param clockwise Logical indicating if the polygon should be oriented clockwise.
#' @return A matrix of polygon vertices with the desired orientation.
fix_orientation <- function(mat, clockwise = FALSE) {
    x <- mat[, 1]
    y <- mat[, 2]
    # signed area (shoelace, ignoring duplicated endpoint)
    area <- sum(x[-nrow(mat)] * y[-1] - x[-1] * y[-nrow(mat)])
    if ((area > 0 && clockwise) || (area < 0 && !clockwise)) {
        mat <- mat[nrow(mat):1, , drop = FALSE]
    }
    mat
}

#' DBSCAN_filter
#'
#' DBSCAN-based filtering of points
#' @param points A numeric matrix of points (n x d).
#' @param min_points The minimum number of points to form a dense region.
#' @param epsilon The radius of the neighborhood to consider.
#' @param quantile_epsilon The quantile of the k-NN distances to use as epsilon if epsilon is NULL.
#' @param B The number of bootstrap samples.
#' @param sample_fraction The fraction of points to sample in each bootstrap.
#' @param CL The confidence level for filtering points.
#' @return A matrix of filtered points.
dbscan_filter <- function(points, min_points = 5, epsilon = NULL,
                          quantile_epsilon = 0.95, B = 2000,
                          sample_fraction = 0.95, CL = 0.95) {
    points <- as.matrix(points)
    n <- nrow(points)
    d <- ncol(points)

    if (is.null(epsilon)) {
        epsilon <- auto_epsilon(points, min_points, quantile_epsilon)
    }

    # m-out-of-n sample size
    m <- max(min(n - 1, floor(n * sample_fraction)), min_points + 1)

    # m-out-of-n density scaling:
    #
    #   epsilon_m = epsilon * (n/m)^(1/d)
    #
    eps_m <- epsilon * (n / m)^(1 / d)

    inclusion_counts <- integer(n)
    sampled_counts <- integer(n)

    for (b in 1:B) {
        idx <- sample.int(n, size = m, replace = FALSE)
        sampled_counts[idx] <- sampled_counts[idx] + 1L

        db <- dbscan::dbscan(points[idx, , drop = FALSE],
            eps = eps_m,
            minPts = min_points
        )

        clustered_local <- which(db$cluster != 0)

        inclusion_counts[idx[clustered_local]] <-
            inclusion_counts[idx[clustered_local]] + 1L
    }

    # Estimate per-point inclusion probability
    inclusion_prob <- inclusion_counts / pmax(sampled_counts, 1L)

    # CL now means PB( point is “clustered” | sampled ) >= CL
    keep <- which(inclusion_prob >= 1 - CL & sampled_counts > 0)

    points[keep, , drop = FALSE]
}

#' alpha_shape_2D
#'
#' Compute the 2D alpha shape of a set of points
#' @param points A numeric matrix of points (n x 2).
#' @param alpha The alpha parameter controlling the shape detail.
#' @param CL The confidence level for DBSCAN filtering (default is 1.0, no filtering).
#' @param sampling_fraction The fraction of points to sample in each bootstrap for DBSCAN filtering.
#' @return A list containing:
#' \item{polygons}{A list of polygons representing the alpha shape boundary. Each polygon is a matrix (m x 2) of coordinates.}
#' \item{points}{The filtered points used to compute the alpha shape.}
#'
#' @export
alpha_shape_2D <- function(points, alpha, CL = 1.0, sampling_fraction = 0.95) {
    points <- as.matrix(points)
    if (ncol(points) != 2) stop("points must be n x 2 for 2D alpha shape")

    if (CL < 1) {
        num_points <- max(4, floor(nrow(points) * CL))
        points <- dbscan_filter(points, min_points = 5, CL = CL, sample_fraction = sampling_fraction)
    }

    if (nrow(points) < 4) {
        warning("Not enough points to form an alpha shape after filtering.")
        return(list(polygons = list(), points = points))
    }

    triangles <- deluanay_triangulation(points, dim = 2)
    if (is.null(triangles) || nrow(triangles) == 0) {
        return(list())
    }

    edge_count <- list()
    keep_tri_idx <- c()
    for (i in 1:nrow(triangles)) {
        tri <- triangles[i, ]
        p1 <- points[tri[1], ]
        p2 <- points[tri[2], ]
        p3 <- points[tri[3], ]

        cc <- circumcircle(p1, p2, p3)
        # skip degenerate (collinear) triangles: cc$r may be Inf
        if (is.finite(cc$r) && cc$r < alpha) {
            # keep the triangle
            keep_tri_idx <- c(keep_tri_idx, i)
            # count its edges
            edges <- list(
                paste(sort(c(tri[1], tri[2])), collapse = "-"),
                paste(sort(c(tri[2], tri[3])), collapse = "-"),
                paste(sort(c(tri[3], tri[1])), collapse = "-")
            )
            for (e in edges) {
                if (is.null(edge_count[[e]])) edge_count[[e]] <- 1 else edge_count[[e]] <- edge_count[[e]] + 1
            }
        }
    }

    # boundary edges are those that appear exactly once among kept triangles
    boundary_keys <- names(edge_count)[sapply(edge_count, function(x) x == 1)]
    if (length(boundary_keys) == 0) {
        # No boundary edges found - this means either no triangles kept or all triangles are internal
        if (length(keep_tri_idx) == 0) {
            warning(paste("No triangles kept with alpha =", alpha, ". Try increasing alpha value."))
        }
        return(list(polygons = list(), points = points))
    }

    # convert keys back to pairs
    boundary_edges <- lapply(boundary_keys, function(k) as.numeric(unlist(strsplit(k, "-"))))

    adj <- build_adj_my(boundary_edges)
    cycles <- find_cycles_my(adj)

    # build ordered loops (polygons)
    loops <- lapply(cycles, function(idx) {
        poly <- points[as.numeric(idx), , drop = FALSE]

        # ensure closure
        if (!all(poly[1, ] == poly[nrow(poly), ])) {
            poly <- rbind(poly, poly[1, ])
        }

        poly
    })


    if (length(loops) == 0) {
        warning("No loops formed from boundary edges.")
        return(list(polygons = list(), points = points))
    }

    # sort loops by absolute area descending (outer loop is largest by area)
    areas <- sapply(loops, polygon_area)
    order_idx <- order(areas, decreasing = TRUE)
    loops <- loops[order_idx]

    # orientation: outer CCW, inner CW
    if (length(loops) >= 1) {
        loops[[1]] <- fix_orientation(loops[[1]], clockwise = FALSE) # outer ring
    }
    if (length(loops) >= 2) {
        for (i in 2:length(loops)) {
            loops[[i]] <- fix_orientation(loops[[i]], clockwise = TRUE) # holes
        }
    }


    return(list(polygons = loops, points = points))
}
