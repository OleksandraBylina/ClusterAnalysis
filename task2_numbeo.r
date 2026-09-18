# Частина 2
# Зчитування даних Numbeo

library(cluster)

prices <- read.csv(
  "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/numbeo_prices_eur.csv",
  sep = ";",
  dec = ".",
  row.names = 1,
  na.strings = "NA",
  check.names = FALSE
)

str(prices)
sapply(prices, class)

# Стандартизація даних

prices_scaled <- scale(prices)

# Розбиття на кластери методом Хартігана-Вонга

dev.set(3)

kl_values <- 2:12

wss_hartigan <- numeric(11)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  km.res <- kmeans(
    x = prices_scaled,
    centers = kl,
    iter.max = 100,
    algorithm = "Hartigan-Wong",
    nstart = 50
  )
  wss_hartigan[i] <- km.res$tot.withinss
}

plot(kl_values, wss_hartigan, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Внутрішньо-кластерна сума", main = "Метод Хартігана-Вонга", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/HartiganWong_wss_numbeo.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()

# Розбиття на кластери методом PAM з евклідовою відстанню

dev.set(3)

dist_euclidean <- dist(
  prices_scaled,
  method = "euclidean"
)

kl_values <- 2:12

sil_euclidean <- numeric(11)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  pam.res <- pam(
    dist_euclidean,
    k = kl,
    diss = TRUE
  )
  sil_euclidean[i] <- pam.res$silinfo$avg.width
}

plot(kl_values, sil_euclidean, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Середнє значення ширини силуету", main = "PAM з евклідовою відстанню", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Euclidean_sil_numbeo.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()

# Розбиття на кластери методом Хартігана-Вонга з уже визначеною кількістю кластерів

op_kl <- 5
set.seed(12345)

hartigan.res <- kmeans(
  x = prices_scaled,
  centers = op_kl,
  iter.max = 100,
  algorithm = "Hartigan-Wong",
  nstart = 50
)

# Розбиття на кластери методом PAM з уже визначеною кількістю кластерів

euclidean.res <- pam(
  dist_euclidean,
  k = op_kl,
  diss = TRUE
)

# Визначення індекса Ренда для двох розбиттів

rand.index <- function(cl1, cl2){
  n <- length(cl1)
  a <- 0
  b <- 0
  c <- 0
  d <- 0
  for (i in 1:(n - 1)){
    for (j in (i + 1):n){
      same1 <- cl1[i] == cl1[j]
      same2 <- cl2[i] == cl2[j]
      if (same1 && same2) {
        a <- a + 1
      } else if (!same1 && !same2) {
        b <- b + 1
      } else if (same1 && !same2) {
        c <- c + 1
      } else {
        d <- d + 1
      }
    }
  }
  return((a + b) / (a + b + c + d))
}

# Визначення індекса Ренда для розбиттів

rand_hartigan_pam <- rand.index(
  hartigan.res$cluster,
  euclidean.res$clustering
)

rand_hartigan_pam

# Таблиця спряженості для розбиттів

tab_hartigan_pam <- table(
  hartigan.res$cluster,
  euclidean.res$clustering
)

tab_hartigan_pam

# Функція для визначення виправленого індекса Ренда для двох розбиттів

adjusted_rand_index <- function(cl1, cl2) {
  
  tab <- table(cl1, cl2)
  
  comb2 <- function(x) {
    x * (x - 1) / 2
  }
  
  sum_ij <- sum(comb2(tab))
  row_sums <- rowSums(tab)
  col_sums <- colSums(tab)
  sum_rows <- sum(comb2(row_sums))
  sum_cols <- sum(comb2(col_sums))
  n <- sum(tab)
  total_pairs <- comb2(n)
  expected <- (sum_rows * sum_cols) / total_pairs
  max_index <- 0.5 * (sum_rows + sum_cols)
  ari <- (sum_ij - expected) / (max_index - expected)
  
  return(ari)
}

# Визначення виправленого індекса Ренда для розбиттів

ari_hartigan_pam <- adjusted_rand_index(
  hartigan.res$cluster,
  euclidean.res$clustering
)

ari_hartigan_pam

# Матричні діаграми розсіювання

pal <- c(
  "magenta",
  "red",
  "blue",
  "yellow",
  "purple",
  "green",
  "orange",
  "brown",
  "cyan",
  "grey",
  "black"
)

# Матрична діаграма розсіювання методом Хартігана-Вонга

dev.set(3)

pairs(
  prices_scaled,
  pch = 16,
  cex = 0.6,
  col = pal[hartigan.res$cluster],
  main = "Метод Хартігана-Вонга"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Hartigan_pairs_numbeo.png",
  width = 1800,
  height = 1800,
  res = 180
)

dev.off()

# Матрична діаграма розсіювання методом PAM з евклідовою відстанню

dev.set(3)

pairs(
  prices_scaled,
  pch = 16,
  cex = 0.6,
  col = pal[euclidean.res$clustering],
  main = "PAM з евклідовою відстанню"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/PAM_pairs_numbeo.png",
  width = 1800,
  height = 1800,
  res = 180
)

dev.off()