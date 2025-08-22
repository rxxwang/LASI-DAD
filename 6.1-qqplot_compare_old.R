library(magick)

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250702_EWAS2/qqplot/")
adbio <- c("r1gfap", "r1nfl", "r1ptau181", "r1totaltau", "abeta_ratio")

img1 <- image_read(paste0("qqplot_PCA_ld_", adbio[1], "_log_", adbio[1], ".png"))
img2 <- image_read(paste0("qqplot_PCA_ld_", adbio[2], "_log_", adbio[2], ".png"))
img3 <- image_read(paste0("qqplot_PCA_ld_", adbio[3], "_log_", adbio[3], ".png"))
img4 <- image_read(paste0("qqplot_PCA_ld_", adbio[4], "_log_", adbio[4], ".png"))
img5 <- image_read(paste0("qqplot_PCA_ld_", adbio[5], "_log_", adbio[5], ".png"))

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3")
adbio <- c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")
biomarker = c("GFAP", "NfL", "pTau181", "Totaltau", "Abeta42/40")

img6 <- image_read(paste0("qqplot/qqplot_model2_", adbio[1], "_log_", adbio[1],".png"))
img7 <- image_read(paste0("qqplot/qqplot_model2_", adbio[2], "_log_", adbio[2],".png"))
img8 <- image_read(paste0("qqplot/qqplot_model2_", adbio[3], "_log_", adbio[3],".png"))
img9 <- image_read(paste0("qqplot/qqplot_model2_", adbio[4], "_log_", adbio[4],".png"))
img10 <- image_read(paste0("qqplot/qqplot_model2_", adbio[5], "_log_", adbio[5],".png"))

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
img1 <- add_title_above(img1, paste0("old data:", biomarker[1]))
img2 <- add_title_above(img2, paste0("old data:", biomarker[2]))
img3 <- add_title_above(img3, paste0("old data:", biomarker[3]))
img4 <- add_title_above(img4, paste0("old data:", biomarker[4]))
img5 <- add_title_above(img5, paste0("old data:", biomarker[5]))

img6 <- add_title_above(img6, paste0("New data:", biomarker[1]))
img7 <- add_title_above(img7, paste0("New data:", biomarker[2]))
img8 <- add_title_above(img8, paste0("New data:", biomarker[3]))
img9 <- add_title_above(img9, paste0("New data:", biomarker[4]))
img10 <- add_title_above(img10, paste0("New data:", biomarker[5]))

row1 <- image_append(c(img1, img2, img3, img4, img5))
row2 <- image_append(c(img6, img7, img8, img9, img10))

combined <- image_append(c(row1, row2), stack = TRUE)
image_write(combined,
            path = paste0("qqplot/qqplot_compare_old.png"))