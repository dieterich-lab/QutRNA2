library(ggplot2)

option_list <- list(
  optparse::make_option(c("-f", "--format"),
                        type = "character",
                        default = "pdf_document",
                        help = "Output format"),
  optparse::make_option(c("-p", "--params"),
                        type = "character",
                        help = "YAML file with params")
)
opts <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list),
  positional_arguments = TRUE
)

stopifnot(!is.null(opts$args))

params <- yaml::yaml.load_file(opts$options$params)
params <- list()
rmarkdown::render(input = opts$args,
                  output_format = opts$options$format,
                  params = params)
