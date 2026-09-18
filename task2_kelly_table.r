grid <- read.csv2(
  "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Kelly_grid_expanded.csv",
  skip = 2,
  row.names = 1,
  na.strings = "NA",
  check.names = FALSE
)
str(grid)
sapply(grid, class)

kl_values <- 2:12

wss_lloyd <- numeric(11)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  km.res <- kmeans(
    x = grid,
    centers = kl,
    iter.max = 100,
    algorithm = "Hartigan-Wong",
    nstart = 50
  )
  wss_lloyd[i] <- km.res$tot.withinss 
}

plot(kl_values, wss_lloyd, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Внутрішньо-кластерна сума", main = "Метод Хартігана-Вонга", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/HartiganWong_wss_characters_part2.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()





grid_t <- t(grid)
kl_values <- 2:10

wss_lloyd <- numeric(9)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  km.res <- kmeans(
    x = grid_t,
    centers = kl,
    iter.max = 100,
    algorithm = "Hartigan-Wong",
    nstart = 50
  )
  wss_lloyd[i] <- km.res$tot.withinss 
}

plot(kl_values, wss_lloyd, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Внутрішньо-кластерна сума", main = "Метод Хартігана-Вонга", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/HartiganWong_wss_conctructs_part2.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()



dist_euclidean <- dist(grid, method = "euclidean")
kl_values <- 2:12

sil_euclidean <- numeric(11)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  eucl.res <- pam(dist_euclidean, k = kl, diss = TRUE)
  sil_euclidean[i] <- eucl.res$silinfo$avg.width 
}

plot(kl_values, sil_euclidean, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Середнє значення ширини силуету", main = "PAM з евклідовою відстанню", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Euclidean_sil_characters_part2.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()



dist_euclidean <- dist(grid_t, method = "euclidean")
kl_values <- 2:10

sil_euclidean <- numeric(9)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  eucl.res <- pam(dist_euclidean, k = kl, diss = TRUE)
  sil_euclidean[i] <- eucl.res$silinfo$avg.width 
}

plot(kl_values, sil_euclidean, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Середнє значення ширини силуету", main = "PAM з евклідовою відстанню", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Euclidean_sil_constructs_part2.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()

# Кластеризація конструктів методом Хартігана-Вонга

op_kl_constructs <- 4

hartigan_constructs.res <- kmeans(
  x = grid_t,
  centers = op_kl_constructs,
  iter.max = 100,
  algorithm = "Hartigan-Wong",
  nstart = 50
)

# Кластеризація конструктів методом PAM з евклідовою відстанню

dist_euclidean_constructs <- dist(grid_t, method = "euclidean")

pam_constructs.res <- pam(
  dist_euclidean_constructs,
  k = op_kl_constructs,
  diss = TRUE
)

pal <- c(
  "magenta",
  "red",
  "blue",
  "yellow"
)

# Матричні діаграми розсіювання конструктів

construct_sets <- list(
  grid_t[, 1:8],
  grid_t[, 9:16],
  grid_t[, 17:24],
  grid_t[, 25:32],
  grid_t[, 33:ncol(grid_t)]
)

# Матричні діаграми розсіювання конструктів методом Хартігана-Вонга

for (i in seq_along(construct_sets)) {
  
  filename <- paste0(
    "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/",
    "Hartigan_constructs_set_", i, ".png"
  )
  
  png(
    filename,
    width = 1800,
    height = 1800,
    res = 180
  )
  
  pairs(
    construct_sets[[i]],
    pch = 16,
    cex = 0.6,
    col = pal[hartigan_constructs.res$cluster],
    main = paste("Метод Хартігана-Вонга: конструкти, набір", i)
  )
  
  dev.off()
}

# Матричні діаграми розсіювання конструктів методом PAM з евклідовою відстанню

for (i in seq_along(construct_sets)) {
  
  filename <- paste0(
    "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/",
    "PAM_constructs_set_", i, ".png"
  )
  
  png(
    filename,
    width = 1800,
    height = 1800,
    res = 180
  )
  
  pairs(
    construct_sets[[i]],
    pch = 16,
    cex = 0.6,
    col = pal[pam_constructs.res$clustering],
    main = paste("PAM з евклідовою відстанню: конструкти, набір", i)
  )
  
  dev.off()
}