merge_chapters = function(files, to, before = NULL, after = NULL, orig = files) {
  # in the preview mode, only use some placeholder text instead of the full Rmd
  preview = FALSE
  input = NULL
  content = unlist(mapply(files, orig, SIMPLIFY = FALSE, FUN = function(f, o) {
    x = xfun::read_utf8(f)
    # if a chapter is short enough (<= 30 lines), just include the full chapter for preview
    x = if (preview && !(o %in% input)) create_placeholder(x) else {
      bookdown:::insert_code_chunk(x, before, after)
    }
    c(x, '', paste0('<!--chapter:end:', o, '-->'), '')
  }))
  if (preview && !(files[1] %in% input))
    content = c(fetch_yaml(xfun::read_utf8(files[1])), content)
  unlink(to)
  xfun::write_utf8(content, to)
  Sys.chmod(to, '644')
}


merge_chapters(c("index.Rmd",
                 #"Species_identity.Rmd",
                 "Callow_markers.Rmd","CHCs_over_time.Rmd",
                 "Non-parametric_tests.Rmd", "Separation_experiment.Rmd","Slave_impact.Rmd",
                 "Dummy_ants.Rmd", "Body_surface.Rmd"), "Supplementary_materials.Rmd")

rmarkdown::render("Supplementary_materials.Rmd", output_file = "Supplementary_materials.pdf")

