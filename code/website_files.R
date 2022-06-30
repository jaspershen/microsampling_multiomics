# distill::create_website(dir = "website_files/", gh_pages = TRUE)
###convert r file to rmarkdown
# knitr::spin("code/ensure_shake/metabolomics/milk_component_in_body.R", knit = FALSE)
# ##copy rmarkdown to website_files
# file.copy("code/ensure_shake/metabolomics/milk_component_in_body.Rmd",
#           to = "website_files/", overwrite = TRUE)
##add this in the begin of the rmarkdown file
# ```{r setup, include=FALSE}
# knitr::opts_chunk$set(echo = TRUE)
# knitr::opts_knit$set(root.dir = '../')
# ```

##render websites
rmarkdown::render_site(input = here::here("website_files/"))
##copy website docs to root folder
file.copy(
  from = here::here("website_files/docs/"),
  to = here::here(),
  recursive = TRUE,
  overwrite = TRUE
)
