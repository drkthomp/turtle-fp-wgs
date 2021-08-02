library(rmarkdown)
library(tidyverse)

types <- c("within", "any")
file_list <- as.data.frame(str_split_fixed(list.files(paste0("partials/cnv/"), recursive = TRUE), "/", 4)) %>%
  mutate(
    minReadCounts = parse_number(str_split_fixed(V4, "-", 4)[, 2]),
    V4 = str_split_fixed(V4, "-", 4)[, 1]
  ) %>%
  distinct() %>%
  na.omit()
names(file_list) <- c("type", "bp", "normType", "sizeFactor", "minReadCount")



apply(file_list, MARGIN = 1, function(x) {
  x <- as.data.frame(t(x))
  dir.create(paste("03-StringDB", x$type, sep = "/"), showWarnings = FALSE)
  dir.create(paste("03-StringDB", x$type, x$bp, sep = "/"), showWarnings = FALSE)
  dir.create(paste("03-StringDB", x$type, x$bp, x$normType, sep = "/"), showWarnings = FALSE)
  # rmarkdown::render("30-StringDB.Rmd",
  #                  params = list(type=x$type, bp=x$bp,minReadCount=x$minReadCount,normType=x$normType,sizeFactor=x$sizeFactor),
  #                  output_file = paste("03-StringDB",x$type,x$bp,x$normType,paste0(x$sizeFactor,"-",x$minReadCount, "_Gene-Analysis"),sep="/"))
  tryCatch(rmarkdown::render("30-StringDB.Rmd",
    params = list(type = x$type, bp = x$bp, minReadCount = x$minReadCount, normType = x$normType, sizeFactor = x$sizeFactor),
    output_file = paste("03-StringDB", x$type, x$bp, x$normType, paste0(x$sizeFactor, "-", x$minReadCount, "_Gene-Analysis"), sep = "/")
  ), error = function(cond) {
    message(cond)
  })
})
