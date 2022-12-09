## FOR USE WITH GI MAPPING PIPELINE

## Rmd options to hold messages
## from: https://stackoverflow.com/questions/55831383/r-markdown-results-hold-but-for-messages
# global store for current chunk messages
current_messages <- NULL

# override default message hook
# if message option is hold, copy messages to the global store and print nothing
# otherwise let knitr do its thing
knitr::knit_hooks$set(message = function(x, options) {
  if (options$message == "hold") {
    current_messages <<- c(current_messages, x)
    return(NULL)
  } else {
    return(knitr::hooks_html()$message(x, options))
  }
})

# override chunk hook
# if message option is hold, clear global store, and append prepared messages
# otherwise let knitr do its thing
knitr::knit_hooks$set(chunk = function(x, options) {
  if (options$message == "hold") {
    collapsed_messages <- paste(current_messages, collapse = "")
    options$message <- TRUE
    messages <- knitr::hooks_html()$message(collapsed_messages, options)
    current_messages <<- NULL
    paste(x, messages, collapse = "")
  } else {
    knitr::hooks_html()$chunk(x, options)
  }
})


## applying results = hold and message = hold to all chunks:
## see: https://github.com/workflowr/workflowr/issues/137
knitr::opts_chunk$set(results = "hold", message = "hold")

