# Частина 1
# Розбиття на кластери методом PAM з евклідовою відстанню

library(cluster)

samp<-read.table("/Users/oleksandrabylina/Documents/University/ClusterAnalysis/mult1.txt", header=TRUE)

dist_euclidean <- dist(samp, method = "euclidean")
kl_values <- 2:20

sil_euclidean <- numeric(19)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  eucl.res <- pam(
    dist_euclidean,
    k = kl,
    diss = TRUE
  )
  sil_euclidean[i] <- eucl.res$silinfo$avg.width
}

plot(kl_values, sil_euclidean, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Середнє значення ширини силуету", main = "PAM з евклідовою відстанню", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Euclidean_sil.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()

# Розбиття на кластери методом PAM з відстанню Манхеттена

dist_manhattan <- dist(samp, method = "manhattan")
kl_values <- 2:20

sil_manhattan <- numeric(19)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  man.res <- pam(
    dist_manhattan,
    k = kl,
    diss = TRUE
  )
  sil_manhattan[i] <- man.res$silinfo$avg.width
}

plot(kl_values, sil_manhattan, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Середнє значення ширини силуету", main = "PAM з відстанню Манхеттена", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Manhattan_sil.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()

# Розбиття на кластери методом PAM з максимальною відстанню

dist_max <- dist(samp, method = "maximum")
kl_values <- 2:20

sil_max <- numeric(19)
set.seed(12345)
for (i in seq_along(kl_values)){
  kl <- kl_values[i]
  max.res <- pam(
    dist_max,
    k = kl,
    diss = TRUE
  )
  sil_max[i] <- max.res$silinfo$avg.width
}

plot(kl_values, sil_max, type ="o", pch = 16, cex = 0.7, lwd = 1.5,col = "steelblue",xlab = "Кількість кластерів", ylab = "Середнє значення ширини силуету", main = "PAM з максимальною відстанню", xaxt = "n")
axis(1, at = kl_values, labels = kl_values, cex.axis = 0.7)
dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Max_sil.png",
  width = 1200,
  height = 800,
  res = 150
)
dev.off()

# Розбиття на кластери всіма медоїдними методами з уже визначеною кількістю кластерів

op_kl <- 11

euclidean.res <- pam(
  dist_euclidean,
  k = op_kl,
  diss = TRUE
)

manhattan.res <- pam(
  dist_manhattan,
  k = op_kl,
  diss = TRUE
)

maximum.res <- pam(
  dist_max,
  k = op_kl,
  diss = TRUE
)

# Визначення індекса Ренда для всіх пар розбиттів

rand_euclidean_manhattan <- rand.index(
  euclidean.res$clustering,
  manhattan.res$clustering
)

rand_euclidean_maximum <- rand.index(
  euclidean.res$clustering,
  maximum.res$clustering
)

rand_manhattan_maximum <- rand.index(
  manhattan.res$clustering,
  maximum.res$clustering
)

rand_euclidean_manhattan
rand_euclidean_maximum
rand_manhattan_maximum

# Таблиці спряженості для всіх пар розбиттів

tab_euclidean_manhattan <- table(
  euclidean.res$clustering,
  manhattan.res$clustering
)

tab_euclidean_maximum <- table(
  euclidean.res$clustering,
  maximum.res$clustering
)

tab_manhattan_maximum <- table(
  manhattan.res$clustering,
  maximum.res$clustering
)

tab_euclidean_manhattan
tab_euclidean_maximum
tab_manhattan_maximum

# Визначення виправленого індекса Ренда для всіх пар розбиттів

ari_euclidean_manhattan <- adjusted_rand_index(
  euclidean.res$clustering,
  manhattan.res$clustering
)

ari_euclidean_maximum <- adjusted_rand_index(
  euclidean.res$clustering,
  maximum.res$clustering
)

ari_manhattan_maximum <- adjusted_rand_index(
  manhattan.res$clustering,
  maximum.res$clustering
)

ari_euclidean_manhattan
ari_euclidean_maximum
ari_manhattan_maximum

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

# Матрична діаграма розсіювання PAM з евклідовою відстанню

dev.set(3)

pairs(
  samp,
  pch = 16,
  cex = 0.5,
  col = pal[euclidean.res$clustering],
  main = "PAM: евклідова відстань"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Euclidean_pairs.png",
  width = 1400,
  height = 1400,
  res = 150
)

dev.off()

# Матрична діаграма розсіювання PAM з відстанню Манхеттена

dev.set(3)

pairs(
  samp,
  pch = 16,
  cex = 0.5,
  col = pal[manhattan.res$clustering],
  main = "PAM: відстань Манхеттена"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Manhattan_pairs.png",
  width = 1400,
  height = 1400,
  res = 150
)

dev.off()

# Матрична діаграма розсіювання PAM з максимальною відстанню

dev.set(3)

pairs(
  samp,
  pch = 16,
  cex = 0.5,
  col = pal[maximum.res$clustering],
  main = "PAM: максимальна відстань"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Maximum_pairs.png",
  width = 1400,
  height = 1400,
  res = 150
)

dev.off()

# Порівняння найкращого центроїдного і медоїдного розбиттів

rand_centroid_pam <- rand.index(
  lloyd.res$cluster,
  euclidean.res$clustering
)

ari_centroid_pam <- adjusted_rand_index(
  lloyd.res$cluster,
  euclidean.res$clustering
)

rand_centroid_pam
ari_centroid_pam

# Таблиця спряженості для центроїдного і медоїдного розбиттів

tab_centroid_pam <- table(
  lloyd.res$cluster,
  euclidean.res$clustering
)

tab_centroid_pam


# Варіант з k = 3
# Розбиття на кластери всіма медоїдними методами при кількості кластерів 3

op_kl2 <- 3

euclidean2.res <- pam(
  dist_euclidean,
  k = op_kl2,
  diss = TRUE
)

manhattan2.res <- pam(
  dist_manhattan,
  k = op_kl2,
  diss = TRUE
)

maximum2.res <- pam(
  dist_max,
  k = op_kl2,
  diss = TRUE
)

# Визначення індекса Ренда для всіх пар розбиттів

rand_euclidean_manhattan2 <- rand.index(
  euclidean2.res$clustering,
  manhattan2.res$clustering
)

rand_euclidean_maximum2 <- rand.index(
  euclidean2.res$clustering,
  maximum2.res$clustering
)

rand_manhattan_maximum2 <- rand.index(
  manhattan2.res$clustering,
  maximum2.res$clustering
)

rand_euclidean_manhattan2
rand_euclidean_maximum2
rand_manhattan_maximum2

# Таблиці спряженості для всіх пар розбиттів

tab_euclidean_manhattan2 <- table(
  euclidean2.res$clustering,
  manhattan2.res$clustering
)

tab_euclidean_maximum2 <- table(
  euclidean2.res$clustering,
  maximum2.res$clustering
)

tab_manhattan_maximum2 <- table(
  manhattan2.res$clustering,
  maximum2.res$clustering
)

tab_euclidean_manhattan2
tab_euclidean_maximum2
tab_manhattan_maximum2

# Визначення виправленого індекса Ренда для всіх пар розбиттів

ari_euclidean_manhattan2 <- adjusted_rand_index(
  euclidean2.res$clustering,
  manhattan2.res$clustering
)

ari_euclidean_maximum2 <- adjusted_rand_index(
  euclidean2.res$clustering,
  maximum2.res$clustering
)

ari_manhattan_maximum2 <- adjusted_rand_index(
  manhattan2.res$clustering,
  maximum2.res$clustering
)

ari_euclidean_manhattan2
ari_euclidean_maximum2
ari_manhattan_maximum2

# Матричні діаграми розсіювання

# Матрична діаграма розсіювання PAM з евклідовою відстанню

dev.set(3)

pairs(
  samp,
  pch = 16,
  cex = 0.5,
  col = pal[euclidean2.res$clustering],
  main = "PAM: евклідова відстань, 3 кластери"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Euclidean_pairs2.png",
  width = 1400,
  height = 1400,
  res = 150
)

dev.off()

# Матрична діаграма розсіювання PAM з відстанню Манхеттена

dev.set(3)

pairs(
  samp,
  pch = 16,
  cex = 0.5,
  col = pal[manhattan2.res$clustering],
  main = "PAM: відстань Манхеттена, 3 кластери"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Manhattan_pairs2.png",
  width = 1400,
  height = 1400,
  res = 150
)

dev.off()

# Матрична діаграма розсіювання PAM з максимальною відстанню

dev.set(3)

pairs(
  samp,
  pch = 16,
  cex = 0.5,
  col = pal[maximum2.res$clustering],
  main = "PAM: максимальна відстань, 3 кластери"
)

dev.copy(
  png,
  filename = "/Users/oleksandrabylina/Documents/University/ClusterAnalysis/Maximum_pairs2.png",
  width = 1400,
  height = 1400,
  res = 150
)

dev.off()

# Порівняння центроїдного і медоїдного розбиттів при однаковій кількості кластерів

rand_centroid_pam2 <- rand.index(
  lloyd.res$cluster,
  euclidean2.res$clustering
)

ari_centroid_pam2 <- adjusted_rand_index(
  lloyd.res$cluster,
  euclidean2.res$clustering
)

rand_centroid_pam2
ari_centroid_pam2

# Таблиця спряженості для центроїдного і медоїдного розбиттів

tab_centroid_pam2 <- table(
  lloyd.res$cluster,
  euclidean2.res$clustering
)

tab_centroid_pam2


