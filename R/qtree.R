get_parent_index <- function (level)
{
  n_tiles_current <- 4^level

  parent_index <- matrix (
    ncol = sqrt (n_tiles_current),
    nrow = sqrt (n_tiles_current)
  )
  ct <- 1
  nreps <- 2
  for (ii in seq (1, sqrt (n_tiles_current), by = nreps))
  {
    for (jj in seq (1, sqrt (n_tiles_current), by = nreps))
    {
      parent_index [ii + 0:(nreps - 1), jj + 0:(nreps - 1)] <- ct
      ct <- ct + 1
    }
  }
  parent_index <- t (parent_index)
  return (parent_index)
}

#' Builds the index structure based on the input data.
#'
#' @param dat_in A \code{sf} point object.
#' @param max_level Deepest level of the quadtree structure.
#' @param force_quadratic Should the bounding box of the data be extended to be
#' quadratic?
#' @param cutoff_function Custom function that defines a cutoff criteria. Must
#' return a boolean value.
#' @import sf
#' @export
build_quadtree <- function (dat_in, max_level, force_quadratic = TRUE,
                            cutoff_function)
{
  cells_out <- vector (mode = "list", length = 4^max_level)

  bb <- sf::st_bbox (dat_in)
  if (force_quadratic)
  {
    diffx <- as.numeric (abs (diff (bb [c (1, 3)])))
    diffy <- as.numeric (abs (diff (bb [c (2, 4)])))

    overlap <- abs (diffy - diffx) / 2

    if (diffx > diffy)
    {
      bb [2] <- bb [2] - overlap
      bb [4] <- bb [4] + overlap
    } else if (diffx < diffy)
    {
      bb [1] <- bb [1] - overlap
      bb [3] <- bb [3] + overlap
    }
  }

  l_dat_in <- nrow (dat_in)
  indices <- matrix (nrow = l_dat_in, ncol = max_level)

  for (i in 1:max_level)
  {
    max_tiles <- 4^i
    n_tiles_parent <- 4^(i - 1)

    # define tile outlines
    diffx_i <- as.numeric (abs (diff (bb [c (1, 3)]))) / sqrt (max_tiles)
    diffy_i <- as.numeric (abs (diff (bb [c (2, 4)]))) / sqrt (max_tiles)

    tile_counter <- 1
    for (jx in 1:sqrt (max_tiles))
    {
      for (jy in 1:sqrt (max_tiles))
      {
        parent_ids <- rep (TRUE, l_dat_in)
        parent_index <- get_parent_index (i)
        if (i > 1)
        {
          parent_id <- parent_index [tile_counter]

          parent_ids <- indices [, i - 1] == parent_id
          parent_ids [is.na (parent_ids)] <- FALSE
        }
        parent_index_flip <- apply (parent_index, 2, rev)
        ix_row <- row (parent_index_flip) [tile_counter]
        ix_col <- col (parent_index_flip) [tile_counter]
        bb_j <- bb
        bb_j [1] <- bb [1] + (diffx_i * (ix_col - 1))
        bb_j [3] <- bb [1] + (diffx_i * ix_col)
        bb_j [2] <- bb [2] + (diffx_i * (ix_row - 1))
        bb_j [4] <- bb [2] + (diffx_i * ix_row)


        dat_parent <- dat_in [parent_ids, ]
        bbj_sf <- sf::st_as_sfc (bb_j)
        bb_j_contains <- sf::st_within (dat_parent, bbj_sf)
        bb_j_contains <- !is.na (as.numeric (bb_j_contains))

        updatable_rows <- rep (TRUE, l_dat_in)
        if (i > 1)
          updatable_rows <- indices [, i - 1] == parent_id

        terminate <- FALSE
        if (!missing (cutoff_function) & any (bb_j_contains))
        {
          terminate <- cutoff_function (dat_in [updatable_rows, ])
        }
        sign <- ifelse (terminate, -1, 1)

        if (terminate | (i == max_level))
        {
          if (length (which (bb_j_contains)) > 0)
          {
            cell_id <- which (sapply (cells_out, is.null)) [1]
            bbj_sf <- sf::st_sf (bbj_sf)
            bbj_sf$level <- i
            bbj_sf$terminate <- terminate
            cells_out [[cell_id]] <- bbj_sf
          }
        }

        indices [which (updatable_rows) [bb_j_contains], i] <- tile_counter * sign
        tile_counter <- tile_counter + 1
      }
    }
  }

  res <- list ()
  res$index <- indices
  cells_out <- do.call (rbind, cells_out)
  res$cells <- cells_out
  return (res)
}
