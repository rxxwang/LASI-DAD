library(magick)

adbio = c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")
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
title = adbio
for(i in 1:5){
  img[[i]] = image_read(snakemake@input[["qqplot"]][i])
  img[[i]] <- add_title_above(img[[i]], title[i])
}

row1 <- image_append(c(img[[1]], img[[2]], img[[3]]))
row2 <- image_append(c(img[[4]], img[[5]]))

combined <- image_append(c(row1, row2), stack = TRUE)
image_write(combined,
            path = snakemake@output[["combined_qqplot"]])
