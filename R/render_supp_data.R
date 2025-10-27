# function code based on bookdown:::merge_chapters

merge_chapters = function(files, to, before = NULL, after = NULL, orig = files) {
  # in the preview mode, only use some placeholder text instead of the full Rmd
  preview = FALSE
  input = NULL
  content = unlist(mapply(files, SIMPLIFY = FALSE, FUN = function(f) {
    x = xfun::read_utf8(f)
    # if a chapter is short enough (<= 30 lines), just include the full chapter for preview
    x = bookdown:::insert_code_chunk(x, before, after)
    c(x, '', paste0('<!--chapter:end:', f, '-->'), '')
  }))
  unlink(to)
  xfun::write_utf8(content, to)
  Sys.chmod(to, '644')
}

#' @export
render_supp_materials <- function(save_globals=FALSE, output_path = NULL){
  if(is.null(output_path)){
    output_path <- getwd()
  }
  chapters <- c("index.Rmd",
                "Species_identity.Rmd",
                "Callow_markers.Rmd",
                "CHCs_over_time.Rmd",
                "Non-parametric_tests.Rmd",
                "Separation_experiment.Rmd",
                "Slave_impact.Rmd",
                "Dummy_ants.Rmd",
                "Body_surface.Rmd",
                "Session_info.Rmd")
  chapters <- unlist(lapply(chapters, function(x) file.path(system.file("supp_materials/",package="sangchem"), x)))
  merge_chapters(chapters, file.path(system.file("supp_materials/",package="sangchem"), "Supplementary_materials.Rmd"))

  rmarkdown::render(system.file("supp_materials/Supplementary_materials.Rmd",package="sangchem"),
                    output_file = file.path(output_path, "Supplementary_materials.pdf"),
                    params=list(save_globals=save_globals))
}

