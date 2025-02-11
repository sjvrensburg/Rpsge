UnorderedUniqueSet <- R6::R6Class(
  "UnorderedUniqueSet",
  private = list(
    elements = list(), # Still using a list to store elements (order doesn't matter now)
    lookup_env = new.env(hash = TRUE, parent = emptyenv()) # Environment for O(1) lookups
  ),
  public = list(
    #' Initialize the UnorderedUniqueSet
    #'
    #' @return An UnorderedUniqueSet object.
    initialize = function() {
      private$elements <- list()
      private$lookup_env <- new.env(hash = TRUE, parent = emptyenv())
    },

    #' Add an element to the UnorderedUniqueSet.
    #'
    #' If the element is already present, it is not added again (uniqueness).
    #' Order is NOT preserved (it's an unordered set).
    #'
    #' @param element The element to add.
    #' @return TRUE if the element was added (i.e., it was new), FALSE otherwise.
    add = function(element) {
      key_str <- digest::digest(element) # Use digest to create a hash key

      if (!exists(key_str, envir = private$lookup_env, inherits = FALSE)) {
        private$elements <- c(private$elements, list(element)) # Append to list (order doesn't matter)
        assign(key_str, TRUE, envir = private$lookup_env) # Store key in environment for lookup
        return(TRUE)
      }
      return(FALSE) # Element already exists
    },

    #' Check if an element is in the UnorderedUniqueSet.
    #'
    #' O(1) lookup using the environment.
    #'
    #' @param element The element to check for.
    #' @return TRUE if the element is in the set, FALSE otherwise.
    has = function(element) {
      key_str <- digest::digest(element) # Use digest to create a hash key
      exists(key_str, envir = private$lookup_env, inherits = FALSE)
    },

    #' Remove an element from the UnorderedUniqueSet.
    #'
    #' Order is not maintained, so removal is simpler.
    #'
    #' @param element The element to remove.
    #' @return TRUE if the element was removed, FALSE if it was not found.
    delete = function(element) {
      key_str <- digest::digest(element) # Use digest to create a hash key

      if (exists(key_str, envir = private$lookup_env, inherits = FALSE)) {
        rm(list = key_str, envir = private$lookup_env) # Remove from lookup environment

        # Remove from the elements list (still need to find it)
        # Since order doesn't matter, we can be slightly more efficient (though still O(n) in list removal)
        index_to_remove <- -1
        for (i in seq_along(private$elements)) {
          if (identical(private$elements[[i]], element)) {
            index_to_remove <- i
            break
          }
        }
        if (index_to_remove != -1) {
          # Just remove the element at that index (order doesn't matter)
          if (length(private$elements) == 1) {
            private$elements <- list() # Empty list if it was the only element
          } else {
            private$elements <- private$elements[-index_to_remove] # Remove element at index
          }
        }
        return(TRUE)
      }
      return(FALSE) # Element not found
    },

    #' Get all elements in the UnorderedUniqueSet (order is not guaranteed).
    #'
    #' @return A list of elements (order may vary).
    values = function() {
      private$elements
    },

    #' Get the number of elements in the UnorderedUniqueSet.
    #'
    #' @return The size of the set.
    size = function() {
      length(private$elements)
    },

    #' Print the UnorderedUniqueSet (for debugging/inspection)
    print = function() {
      cat("<UnorderedUniqueSet>\n")
      if (self$size() > 0) {
        cat("  Elements (unordered):\n") # Indicate unordered nature
        for (elem in private$elements) {
          cat("    - ", deparse(elem), "\n")
        }
      } else {
        cat("  (empty)\n")
      }
      invisible(self)
    }
  )
)

OrderedUniqueSet <- R6::R6Class(
  "OrderedUniqueSet",
  private = list(
    elements = list(), # Use a list to maintain insertion order
    lookup_env = new.env(hash = TRUE, parent = emptyenv()) # Environment for O(1) lookups
  ),
  public = list(
    #' Initialize the OrderedUniqueSet
    #'
    #' @return An OrderedUniqueSet object.
    initialize = function() {
      private$elements <- list()
      private$lookup_env <- new.env(hash = TRUE, parent = emptyenv())
    },

    #' Add an element to the OrderedUniqueSet.
    #'
    #' If the element is already present, it is not added again (uniqueness).
    #' Insertion order is preserved.
    #'
    #' @param element The element to add.
    #' @return TRUE if the element was added (i.e., it was new), FALSE otherwise.
    add = function(element) {
      key_str <- digest::digest(element) # Use digest to create a hash key

      if (!exists(key_str, envir = private$lookup_env, inherits = FALSE)) {
        private$elements <- c(private$elements, list(element)) # Append to list to maintain order
        assign(key_str, TRUE, envir = private$lookup_env) # Store key in environment for lookup
        return(TRUE)
      }
      return(FALSE) # Element already exists
    },

    #' Check if an element is in the OrderedUniqueSet.
    #'
    #' O(1) lookup using the environment.
    #'
    #' @param element The element to check for.
    #' @return TRUE if the element is in the set, FALSE otherwise.
    has = function(element) {
      key_str <- digest::digest(element) # Use digest to create a hash key
      exists(key_str, envir = private$lookup_env, inherits = FALSE)
    },

    #' Remove an element from the OrderedUniqueSet.
    #'
    #' Maintains insertion order of remaining elements.
    #'
    #' @param element The element to remove.
    #' @return TRUE if the element was removed, FALSE if it was not found.
    delete = function(element) {
      key_str <- digest::digest(element) # Use digest to create a hash key

      if (exists(key_str, envir = private$lookup_env, inherits = FALSE)) {
        rm(list = key_str, envir = private$lookup_env) # Remove from lookup environment

        # Remove from the ordered list (inefficient O(n) operation in the worst case for lists)
        # Find the index of the element (value comparison is needed)
        index_to_remove <- -1
        for (i in seq_along(private$elements)) {
          if (identical(private$elements[[i]], element)) { # Use identical for object comparison
            index_to_remove <- i
            break
          }
        }
        if (index_to_remove != -1) {
          if (index_to_remove == 1 && length(private$elements) == 1) {
            private$elements <- list() # Empty the list if it was the only element
          } else if (index_to_remove == 1) {
            private$elements <- private$elements[-1] # Remove first element
          } else if (index_to_remove == length(private$elements)) {
            private$elements <- private$elements[-length(private$elements)] # Remove last element
          } else {
            private$elements <- c(private$elements[1:(index_to_remove - 1)], private$elements[(index_to_remove + 1):length(private$elements)]) # Remove from middle
          }
        }

        return(TRUE)
      }
      return(FALSE) # Element not found
    },

    #' Get all elements in the order they were inserted.
    #'
    #' @return A list of elements in insertion order.
    values = function() {
      # Return the elements in their insertion order
      # This maintains the order of elements while providing access to them
      return(private$elements)
    },

    #' Get the number of elements in the OrderedUniqueSet.
    #'
    #' @return The size of the set.
    size = function() {
      length(private$elements)
    },

    #' Print the OrderedUniqueSet (for debugging/inspection)
    print = function() {
      # R equivalent of Python's __str__
      cat("OrderedUniqueSet(", "\n")

      if (length(private$elements) > 0) {
        # Print each element with proper formatting
        for (i in seq_along(private$elements)) {
          # Handle different types of elements
          elem <- private$elements[[i]]
          if (is.character(elem)) {
            elem_str <- paste0("'", elem, "'")
          } else if (is.null(elem)) {
            elem_str <- "NULL"
          } else {
            elem_str <- as.character(elem)
          }

          # Add comma if not last element
          if (i < length(private$elements)) {
            cat("  ", elem_str, ",\n", sep="")
          } else {
            cat("  ", elem_str, "\n", sep="")
          }
        }
      }

      cat(")")

      # Return invisibly for method chaining
      invisible(self)
    }
  )
)
