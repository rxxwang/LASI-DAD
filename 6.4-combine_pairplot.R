library(magick)

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3")
adbio <- c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")
biomarker = c("GFAP", "NfL", "pTau181", "Totaltau", "Abeta42/40")
file = paste0("pairplot/scatterplot_", biomarker,"_old.png")

img1 <- image_read(paste0("pairplot/scatterplot_", adbio[1], "_old.png"))
img2 <- image_read(paste0("pairplot/scatterplot_", adbio[2], "_old.png"))
img3 <- image_read(paste0("pairplot/scatterplot_", adbio[3], "_old.png"))
img4 <- image_read(paste0("pairplot/scatterplot_", adbio[4], "_old.png"))
img5 <- image_read(paste0("pairplot/scatterplot_", adbio[5], "_old.png"))


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

# Add subtitles
img1 <- add_title_above(img1, biomarker[1])
img2 <- add_title_above(img2, biomarker[2])
img3 <- add_title_above(img3, biomarker[3])
img4 <- add_title_above(img4, biomarker[4])
img5 <- add_title_above(img5, biomarker[5])

row1 <- image_append(c(img1, img3, img5))
row2 <- image_append(c(img2, img4))

combined <- image_append(c(row1, row2), stack = TRUE)
image_write(combined,
            path = paste0("pairplot/scatterplot_compare_old.png"))