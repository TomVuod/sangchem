.load_globals <- function(.env){
  globals_list <- readRDS(file.path(system.file(package="sangchem"), "globals.rds"))
  lapply(ls(globals_list), function(obj.name) assign(obj.name, get(obj.name, pos = globals_list), envir = .env))
}

#' @import rlang
# .onLoad <- function(libname, pkgname){
#   run_on_load()
# }
#
# on_load({
#   .load_globals(sys.frame(sys.nframe()))
#   })

.onLoad <- function(libname, pkgname){
  .load_globals(parent.env(sys.frame(sys.nframe())))
}


#' @export
load_globals <- function(){
  .env <- parent.env(sys.frame(sys.nframe()))
  globals_list <- readRDS(file.path(system.file(package="sangchem"), "globals.rds"))
  lapply(ls(globals_list), function(obj.name) assign(obj.name, get(obj.name, pos = globals_list), envir = .env))
}
