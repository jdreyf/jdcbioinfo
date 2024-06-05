#-----------------------------------------------------------------------------
#' Change the formal arguments of a function
#'
#' If the function is specified within a package, then update the function
#' within the namespace of the package.
#'
#' @param function_name character name of function
#' @param package_name character name of package. Set this if you want the function
#'                     to be updated within the package namespace.  If this is unset
#'                     then the changed function is placed in the specified environment
#' @param envir environment to place the function if package_name not set
#'
#' @return TRUE otherwise should throw an error
#' @export
#-----------------------------------------------------------------------------
update_function_arguments <- function(function_name, package_name=NULL, envir=parent.frame(), ...) {

  #---------------------------------------------------------------------------
  # rlang::quos() drops an empty argument if it's the last one.
  #---------------------------------------------------------------------------
  dots  <- rlang::list2(...)

  #---------------------------------------------------------------------------
  # Get the named function and its formal arguments.
  # If `package_name` is defined, then get the function from within that
  # package, otherwise
  #---------------------------------------------------------------------------
  if (is.null(package_name)) {
    func <- get(function_name, envir = envir)
  } else {
    func  <- getFromNamespace(function_name, ns=package_name)
  }
  fargs <- formals(func)

  #---------------------------------------------------------------------------
  # For each named item in ..., set the formal argument
  #---------------------------------------------------------------------------
  for (i in seq(dots)) {
    argument_name <- names(dots)[i]
    fargs[argument_name] <- dots[[i]]
  }

  #---------------------------------------------------------------------------
  # Update the function with new formal arguments
  #---------------------------------------------------------------------------
  formals(func) <- fargs

  #---------------------------------------------------------------------------
  # If no `package_name` is defined, then just assign function in .GlobalEnv.
  # Otherwise if `package_name` is set:
  #  - Get the package as an environment
  #  - Unlock the environment if it is locked
  #  - Assign the new function into the package
  #  - Re-lock the environment if it was initially locked
  #
  # This seems to mostly work.
  #---------------------------------------------------------------------------
  if (is.null(package_name)) {
    assign(function_name, func, envir=envir)
  } else {
    package_env <- as.environment(paste0('package:', package_name))

    locked <- bindingIsLocked(function_name, package_env)
    if (locked) { unlockBinding(function_name, package_env) }

    package_env[[function_name]] <- func

    if (locked) { lockBinding(function_name, package_env) }
  }

  #---------------------------------------------------------------------------
  # Return quietly
  #---------------------------------------------------------------------------
  invisible(TRUE)
}
