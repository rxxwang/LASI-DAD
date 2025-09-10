library(magick)

biomarker = snakemake@params[["biomarker"]]
add_title_above <- function(img,
                            title_text,
                            title_height = 100,
                            font_size = 50) {
  width <- image_info(img)$width
  title_canvas <- image_blank(width = width,
                              height = title_height,
                              color = "white")
  title_canvas <- image_annotate(
    title_canvas,
    title_text,
    gravity = "center",
    size = font_size,
    color = "black",
    font = "Arial"
  )
  image_append(c(title_canvas, img), stack = TRUE)
}

img = list()
title = c(paste0("A. Model 1: log(", biomarker, ")"),
          paste0("B. Model 2: log(", biomarker, ") + PCs"),
          paste0("C. Model 3: winsorized log(", biomarker, ")"),
          paste0("D. Model 4: winsorized log(", biomarker, ") + PCs"),
          paste0("E. Model 5: inverse normalized ", biomarker),
          paste0("F. Model 6: inverse normalized ", biomarker, " + PCs"))
for(i in 1:6){
  img[[i]] = image_read(snakemake@input[["qqplot"]][i])
  img[[i]] <- add_title_above(img[[i]], title[i])
}

row1 <- image_append(c(img[[1]], img[[3]], img[[5]]))
row2 <- image_append(c(img[[2]], img[[4]], img[[6]]))

combined <- image_append(c(row1, row2), stack = TRUE)
image_write(combined,
            path = snakemake@output[["combined_qqplot"]])
