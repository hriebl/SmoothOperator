# TODO: what about null.space.penalty?

supported <- c("tp", "ts", "ds", "cr", "cs", "cc", "ps", "cp", "gp")

#' R6 class for a penalized spline smooth
#'
#' @description
#' The penalized spline smooths are constructed using the
#' [mgcv][mgcv::mgcv.package] package. See also:
#'
#' - [`mgcv::s()`][mgcv::s()]
#' - [`mgcv::smooth.terms`][mgcv::smooth.terms]
#' - [`mgcv::smoothCon()`][mgcv::smoothCon()]
#'
#' @importFrom mgcv Predict.matrix s smooth.construct smoothCon
#' @importFrom R6 R6Class
#' @export

Smooth <- R6Class(
  classname = "Smooth",
  public = list(
    #' @field terms The names of the covariates as a character vector.
    terms = NULL,

    #' @field data The data frame from which the covariates are extracted.
    data = NULL,

    #' @field bs A basis abbreviation as defined in mgcv.
    #'           Defaults to `"tp"`.
    bs = "tp",

    #' @field k The basis dimension before constraints and penalization.
    #'          Defaults to 10.
    k = 10,

    #' @field m The order of the penalty.
    #'          Defaults to `NA` for auto-initialization.
    m = NA,

    #' @field xt Extra information to set up the basis.
    #'           Defaults to `NULL`.
    xt = NULL,

    #' @field knots The knots for the basis construction as a data frame with
    #'              the terms as names. Defaults to `NULL` for auto-initialization
    #'              with mgcv's basis-specific method.
    knots = NULL,

    #' @field constraints A matrix of (identifiability) constraints. Defaults
    #'                    to `NULL`, in which case only a centering constraint
    #'                    will be applied.
    constraints = NULL,

    #' @field absorb_constraints If the constraints should be absorbed into the
    #'                           basis. Defaults to `TRUE`.
    absorb_constraints = TRUE,

    #' @field scale_penalty If the penalty matrices should be scaled to match
    #'                      the inner product of the design matrix. Defaults to
    #'                      `TRUE`.
    scale_penalty = TRUE,

    #' @field identity_penalty If the smooth should be reparameterized to turn
    #'                         the penalty matrices into identity matrices. This
    #'                         is more or less equivalent to
    #'                         [`mgcv::smooth2random()`][mgcv::smooth2random()].
    #'                         Defaults to `FALSE`.
    identity_penalty = FALSE,

    #' @description
    #' Create a new smooth object.
    #' @param terms The names of the covariates as a character vector.
    #' @param data The data frame from which the covariates are extracted.
    #' @param bs A basis abbreviation as defined in mgcv.
    #'           Defaults to `"tp"`.
    #' @param k The basis dimension before constraints and penalization.
    #'          Defaults to 10.
    initialize = function(terms, data, bs = "tp", k = 10) {
      self$terms <- terms
      self$data <- data

      self$bs <- match.arg(bs, supported)

      self$k <- k
    },

    #' @description
    #' Initialize knots using mgcv's basis-specific method.
    initialize_knots = function() {
      smooth <- private$s()

      smooth <- smooth.construct(smooth, self$data, knots = NULL)

      knots <- switch(
        self$bs,
        tp = smooth$Xu,
        ts = smooth$Xu,
        ds = smooth$knt,
        cr = smooth$xp,
        cs = smooth$xp,
        cc = smooth$xp,
        ps = smooth$knots,
        cp = smooth$knots,
        gp = smooth$knt
      )

      shift <- smooth$shift

      if (is.null(shift)) {
        shift <- 0
      }

      width <- length(self$terms)
      knots <- matrix(knots, ncol = width)
      knots <- sweep(knots, 2, shift, "+")
      knots <- as.data.frame(knots)
      names(knots) <- self$terms

      self$knots <- knots
      invisible(self)
    },

    #' @description
    #' Initialize constraints with a centering constraint.
    initialize_constraints = function() {
      self$remove_all_constraints()
      self$add_centering_constraint()
      invisible(self)
    },

    #' @description
    #' Add a centering constraint.
    add_centering_constraint = function() {
      if (is.null(self$constraints)) {
        self$remove_all_constraints()
      }

      smooth <- private$s()
      smooth <- smooth.construct(smooth, self$data, self$knots)
      constraint <- matrix(colMeans(smooth$X), nrow = 1)

      private$add_constraints(constraint)
      invisible(self)
    },

    #' @description
    #' Add a point constraint where the smooth should pass through zero.
    #' @param data A data frame with the terms as names defining where the
    #'             smooth should pass through zero.
    add_point_constraints = function(data) {
      if (is.null(self$constraints)) {
        self$initialize_constraints()
      }

      smooth <- private$s()
      smooth <- smooth.construct(smooth, self$data, self$knots)
      constraint <- Predict.matrix(smooth, data)

      private$add_constraints(constraint)
      invisible(self)
    },

    #' @description
    #' Remove all constraints (including the centering constraint).
    remove_all_constraints = function() {
      width <- self$k

      if (self$bs == "cc") {
        width <- width - 1
      }

      self$constraints <- matrix(nrow = 0, ncol = width)
      invisible(self)
    },

    #' @description
    #' Construct the design matrix and the penalty matrices.
    #' @param data An optional data frame if different data should be used for
    #'             the construction of the design matrix than the basis.
    #' @return A list with the elements `design_matrix`, `penalty_matrices` and
    #'         `ranks`.
    construct = function(data = NULL) {
      if (is.null(data)) {
        data <- self$data
      } else {
        if (is.null(self$knots)) {
          self$initialize_knots()
        }

        if (is.null(self$constraints)) {
          self$initialize_constraints()
        }
      }

      smooth <- private$s()

      if (!is.null(self$constraints)) {
        smooth$C <- self$constraints
      }

      smooth <- smoothCon(
        object           = smooth,
        data             = data,
        knots            = self$knots,
        absorb.cons      = self$absorb_constraints,
        scale.penalty    = self$scale_penalty,
        diagonal.penalty = self$identity_penalty
      )

      list(
        design_matrix = smooth[[1]]$X,
        penalty_matrices = smooth[[1]]$S,
        ranks = smooth[[1]]$rank
      )
    }
  ),
  private = list(
    add_constraints = function(new) {
      old_rank <- qr(self$constraints)$rank
      new_matrix <- rbind(self$constraints, new)
      new_rank <- qr(new_matrix)$rank

      if (new_rank < old_rank + nrow(new)) {
        stop("New constraints are linearly dependent on existing ones")
      }

      self$constraints <- new_matrix
      invisible(self)
    },
    s = function() {
      args <- lapply(self$terms, as.name)

      args$k <- self$k
      args$bs <- self$bs
      args$m <- self$m
      args$xt <- self$xt

      do.call(s, args)
    }
  )
)
