if(exists("echo.val")){ 
  if(echo.val) { 
    cat("---
      bibliography: /Users/mader/Documents/library.bib
    ---")
    library(knitcitations)
    bibliography()
  }
}