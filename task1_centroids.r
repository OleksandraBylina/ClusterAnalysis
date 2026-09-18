# Частина 1
# Розбиття на кластери методом Ллойда

samp<-read.table("/Users/oleksandrabylina/Documents/University/ClusterAnalysis/mult1.txt", header=TRUE)

kl_values <- 2:20

wss_lloyd <- numeric(19)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  km.res <- kmeans(
    x = samp,
    centers = kl,
    iter.max = 100,
    algorithm = "Lloyd",
    nstart = 50
  )
  wss_lloyd[i] <- km.res$tot.withinss 
}

plot(kl_values, wss_lloyd, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Внутрішньо-кластерна сума", main = "Метод Ллойда", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Lloyd_wss.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()

# Розбиття на кластери методом Макквіна

wss_macqueen <- numeric(19)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  km.res <- kmeans(
    x = samp,
    centers = kl,
    iter.max = 100,
    algorithm = "MacQueen",
    nstart = 50
  )
  wss_macqueen[i] <- km.res$tot.withinss 
}

plot(kl_values, wss_macqueen, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Внутрішньо-кластерна сума", main = "Метод Макквіна", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/MacQueen_wss.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()

# Розбиття на кластери методом Хартігана-Вонга

wss_hartigan_wong <- numeric(19)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  km.res <- kmeans(
    x = samp,
    centers = kl,
    iter.max = 100,
    algorithm = "Hartigan-Wong",
    nstart = 50
  )
  wss_hartigan_wong[i] <- km.res$tot.withinss 
}

plot(kl_values, wss_hartigan_wong, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Внутрішньо-кластерна сума", main = "Метод Хартігана-Вонга", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/HartiganWong_wss.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()

# Розбиття на кластери всіма центроїдними методами з уже визначеною кількістю кластерів

op_kl <- 3
set.seed(12345)
lloyd.res <- kmeans(
  x = samp,
  centers = op_kl,
  iter.max = 100,
  algorithm = "Lloyd",
  nstart = 50
)
macqueen.res <- kmeans(
  x = samp,
  centers = op_kl,
  iter.max = 100,
  algorithm = "MacQueen",
  nstart = 50
)
hartigan.res <- kmeans(
  x = samp,
  centers = op_kl,
  iter.max = 100,
  algorithm = "Hartigan-Wong",
  nstart = 50
)

# Функція для визначення індекса Ренда для двох розбиттів

rand_index <- function(cl1, cl2){
  n <- length(cl1)
  a <- 0
  b <- 0
  c <- 0
  d <- 0
  for (i in 1:(n-1)){
    for (j in (i + 1):(n)){
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

# Визначення індекса Ренда для всіх пар розбиттів

rand_lloyd_macqueen <- rand_index(lloyd.res$cluster, macqueen.res$cluster)
rand_lloyd_hartigan <- rand_index(lloyd.res$cluster, hartigan.res$cluster)
rand_macqueen_hartigan <- rand_index(macqueen.res$cluster, hartigan.res$cluster)

rand_lloyd_macqueen
rand_lloyd_hartigan
rand_macqueen_hartigan

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

# Таблиці спряженості для всіх пар розбиттів

tab_lloyd_macqueen <- table(
  lloyd.res$cluster,
  macqueen.res$cluster
)

tab_lloyd_hartigan <- table(
  lloyd.res$cluster,
  hartigan.res$cluster
)

tab_macqueen_hartigan <- table(
  macqueen.res$cluster,
  hartigan.res$cluster
)

tab_lloyd_macqueen
tab_lloyd_hartigan
tab_macqueen_hartigan


# Визначення виправленого індекса Ренда для всіх пар розбиттів

ari_lloyd_macqueen <- adjusted_rand_index(
  lloyd.res$cluster,
  macqueen.res$cluster
)

ari_lloyd_hartigan <- adjusted_rand_index(
  lloyd.res$cluster,
  hartigan.res$cluster
)

ari_macqueen_hartigan <- adjusted_rand_index(
  macqueen.res$cluster,
  hartigan.res$cluster
)

ari_lloyd_macqueen
ari_lloyd_hartigan
ari_macqueen_hartigan

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

# Матрична діаграма розсіювання методом Ллойда

dev.set(3)

pairs(
  samp,
  pch = 16,
  cex = 0.5,
  col = pal[lloyd.res$cluster],
  main = "Метод Ллойда"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Lloyd_pairs.png",
  width = 1400,
  height = 1400,
  res = 150
)

dev.off()


# Матрична діаграма розсіювання методом Макквіна

dev.set(3)

pairs(
  samp,
  pch = 16,
  cex = 0.5,
  col = pal[macqueen.res$cluster],
  main = "Метод Мак-Квіна"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/MacQueen_pairs.png",
  width = 1400,
  height = 1400,
  res = 150
)

dev.off()


# Матрична діаграма розсіювання методом Хартігана-Вонга

dev.set(3)

pairs(
  samp,
  pch = 16,
  cex = 0.5,
  col = pal[hartigan.res$cluster],
  main = "Метод Хартігана-Вонга"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/HartiganWong_pairs.png",
  width = 1400,
  height = 1400,
  res = 150
)

dev.off()