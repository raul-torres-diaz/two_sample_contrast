######################################################################
blocks_class <- function(){
  setClass("blocks", representation(vertices = "matrix",
                                    elements = "matrix",
                                    stats = "numeric",
                                    leftChild = "list",
                                    rightChild = "list",
                                    prefChild = "list"))
}
######################### Split function ############################# 
split <- function(split_dat, cutting, permutation, pvertices, p = FALSE){
  num_data <- nrow(split_dat)
  fun_class <- list(new("blocks", vertices = pvertices, stats = 1 : (num_data + 1), elements = split_dat))
  final_class <- NULL
  for (n1 in 1 : num_data) {
    fun_perm <- permutation[n1]
    fun_cut <- cutting[permutation[n1]]
    for (n2 in 1 : length(fun_class)) {
      if (fun_perm %in% fun_class[[n2]]@stats & length(is.na(fun_class[[n2]]@leftChild)) == 0) {
        data <- fun_class[[n2]]@elements
        temp_dat <- matrix(data[order(data[, fun_cut]), ], ncol = 2)
        
        min_stat <- min(fun_class[[n2]]@stats)
        kth_order <- fun_perm - min_stat + 1
        value_stat <- temp_dat[kth_order, ] 
        cut_point <- value_stat[fun_cut]
        
        stats_left <- fun_class[[n2]]@stats[1 : kth_order]
        stats_right <- fun_class[[n2]]@stats[-(1 : kth_order)]
        
        elements_left <- NULL
        elements_right <- NULL
        
        if (fun_cut == 1) {
          vertices_left <- rbind(fun_class[[n2]]@vertices[1, ], c(cut_point, fun_class[[n2]]@vertices[2, 2]))
          vertices_right <- rbind(c(cut_point, fun_class[[n2]]@vertices[1, 2]), fun_class[[n2]]@vertices[2, ])
          for (l in 1 : nrow(fun_class[[n2]]@elements)) {
            if (temp_dat[l, fun_cut] < cut_point) {
              elements_left <- rbind(elements_left, temp_dat[l, ])
            } else if (temp_dat[l, fun_cut] > cut_point) {
              elements_right <- rbind(elements_right, temp_dat[l, ])
            } else {
              element_child <- temp_dat[l, ]
            }
          }
        } else if (fun_cut == 2) {
          vertices_left <- rbind(fun_class[[n2]]@vertices[1, ], c(fun_class[[n2]]@vertices[2, 1], cut_point))
          vertices_right <- rbind(c(fun_class[[n2]]@vertices[1, 1], cut_point), fun_class[[n2]]@vertices[2, ])
          for (l in 1 : nrow(fun_class[[n2]]@elements)) {
            if (temp_dat[l, fun_cut] < cut_point) {
              elements_left <- rbind(elements_left, temp_dat[l, ])
            } else if (temp_dat[l, fun_cut] > cut_point) {
              elements_right <- rbind(elements_right, temp_dat[l, ])
            } else {
              element_child <- temp_dat[l, ]
            }
          }
        }
        
        if (length(elements_left) == 0) {
          elements_left <- 0
        }
        if (length(elements_right) == 0) {
          elements_right <- 0
        }
        
        child_left <- new("blocks", vertices = vertices_left, 
                          stats = stats_left, elements = as.matrix(elements_left))
        child_right <- new("blocks", vertices = vertices_right, 
                           stats = stats_right, elements = as.matrix(elements_right))
        
        fun_class[[n2]]@leftChild <- list(child_left)
        fun_class[[n2]]@rightChild <- list(child_right)
        fun_class[[n2]]@prefChild <- list(element_child)
        
        fun_class <- append(fun_class, append(fun_class[[n2]]@leftChild, fun_class[[n2]]@rightChild))
        break
      }
    }
    if (p == TRUE){
      if (fun_cut == 1) {
        segments(x0 = as.numeric(cut_point), y0 = vertices_left[1, 2], 
                 x1 = as.numeric(cut_point), y1 = vertices_left[2, 2], 
                 col = 1, lwd = 1)
        
      } else if (fun_cut == 2) {
        segments(y0 = as.numeric(cut_point), x0 = vertices_left[1, 1], 
                 y1 = as.numeric(cut_point), x1 = vertices_left[2, 1], 
                 col = 1, lwd = 1)
      }
    }
  }
  
  for (final in 1 : length(fun_class)){
    if (length(is.na(fun_class[[final]]@leftChild)) == 0){
      final_class <- append(final_class, fun_class[[final]])
    }
  }
  final_class <- order_blocks(final_class)
  return(final_class)
}
######################### Order function ############################# 
order_blocks <- function(Blocks){
  temp_res <- NULL
  N_Blocks <- length(Blocks)
  Assigned_Blocks <- 0
  i <- 1
  while (Assigned_Blocks < (N_Blocks)) {
    if (Blocks[[(i %% N_Blocks) + 1]]@stats == (Assigned_Blocks + 1)) {
      temp_res <- append(temp_res, Blocks[[(i %% N_Blocks) + 1]])
      Assigned_Blocks <- Assigned_Blocks + 1
    }
    i <- i + 1
  }
  return(temp_res)
}
######################### Count function ############################# 
counting_blocks <- function(fun_Blocks, fun_dat){
  temp_res <-  numeric(length(fun_Blocks))  
  fun_dat <- as.data.frame(fun_dat)
  for (k in 1 : length(fun_Blocks)){
    i <- 1
    while (i <= nrow(fun_dat)) {
      f1 <- findInterval(fun_dat[i, 1], c(fun_Blocks[[k]]@vertices[1, 1], c(fun_Blocks[[k]]@vertices[2, 1])), rightmost.closed = TRUE, left.open = TRUE)
      f2 <- findInterval(fun_dat[i, 2], c(fun_Blocks[[k]]@vertices[1, 2], c(fun_Blocks[[k]]@vertices[2, 2])), rightmost.closed = TRUE, left.open = TRUE)
      if (f1 == 1 && f2 == 1) {
        temp_res[k] <- temp_res[k] + 1
        fun_dat <- fun_dat[-i,]
        i <- i - 1
      }
      i <- i + 1
    }
  }
  return(temp_res)
}
######################################################################
prob_null <- function(s, N, M){
  num <- choose(N + 1, s) * choose(M - 1, N - s)
  dem <- choose(M + N, N) 
  return(num/dem)
}
test_block <- function(fun_count, test_k = 0, lam = NA){
  count_m <- sum(fun_count)
  count_n <- length(fun_count) - 1
  count_s <- data.frame("frequency" = 0 : count_m,
                        "counts" = unlist(lapply(0 : count_m, function(tem) sum(fun_count == tem))))
  if (is.na(lam)){
    lam <- count_m/count_n
  }
  p <- unlist(lapply(0 : test_k, function(tem) (lam^tem)/((lam + 1)^(tem + 1))))
  mean_p <- count_n * p
  u <- sum(unlist(lapply(0 : test_k, function(tem) (count_s[(tem + 1), 2] - mean_p[(tem + 1)]) * (tem - lam - test_k - 1))))
  v <- sqrt(lam * (lam + 1)) * sum(unlist(lapply(0 : test_k, function(tem) count_s[(tem + 1), 2] - mean_p[(tem + 1)])))
  Q1 <- sum(unlist(lapply(0 : test_k, function(tem) ((count_s[(tem + 1), 2] - mean_p[(tem + 1)])^2)/mean_p[(tem + 1)])))
  Q2 <- (u^2 + v^2)/(count_n * (lam^2) * (1 + lam) * tail(x = p, n = 1))
  Q <- Q1 + Q2
  res_chi <- pchisq(q = Q, df = test_k + 1, lower.tail = F)
  return(res_chi)
}
test_block_empty <- function(fun_count, value, alpha = 0.05) {
  count_zero <- sum(fun_count == 0)
  pvalue <- value[(count_zero + 1)]
  return(list(as.numeric(pvalue <  alpha), pvalue))
}
test_good <- function(p1, p2){
  com <- 2/(1/p1 + 1/p2)
  return(com)
}
test_bonf<- function(p1, p2){
  com <- 2 * pmin(p1, p2)
  return(com)
}


######################################################################
######################################################################
######################################################################
library(MASS)
N <- 10000
nx = 32
ny = 32
set.seed(2020)
perm_dat_x <- sample(1 : nx, replace = FALSE)
cut_dat_x <- sample(c(1, 2), nx, replace = TRUE)
perm_dat_y <- sample(1 : ny, replace = FALSE)
cut_dat_y <- sample(c(1, 2), ny, replace = TRUE)
blocks_class()
m_X <- c(0, 0)
m_Y <- c(0, 0)
cov_X <- matrix(c(1, 0, 0, 1), 2, 2)
cov_Y <- matrix(c(1, 0, 0, 1), 2, 2)


######################################################################
prob_x <- numeric((nx + 1))
Prob_x <- numeric((nx + 1))
prob_y <- numeric((ny + 1))
Prob_y <- numeric((ny + 1))

for (i in 1 : (nx + 1)){
  prob_x[i] <- prob_null((i - 1), nx, ny)
}

for (i in 1 : (nx + 1)){
  Prob_x[i] <- sum(prob_x[i : (nx + 1)])
}

for (i in 1 : (ny + 1)){
  prob_y[i] <- prob_null((i - 1), ny, nx)
}

for (i in 1 : (ny + 1)){
  Prob_y[i] <- sum(prob_y[i : (ny + 1)])
}

######################################################################
Xmin <- min(m_X[1] - 8 * cov_X[1, 1], m_Y[1] - 8 * cov_Y[1, 1])
Xmax <- max(m_X[1] + 8 * cov_X[1, 1], m_Y[1] + 8 * cov_Y[1, 1])
Ymin <- min(m_X[2] - 8 * cov_X[2, 2], m_Y[2] - 8 * cov_Y[2, 2])
Ymax <- max(m_X[2] + 8 * cov_X[2, 2], m_Y[2] + 8 * cov_Y[2, 2])
plot_vertices <- matrix(c(Xmin, Xmax, Ymin, Ymax), 2, 2)

######################################################################
f_old_x <- f_old_y <- f_old_g <- f_old_b <- numeric(N)
f_0_x <- f_0_y <- f_0_g <- f_0_b <- numeric(N)
f_1_x <- f_1_y <- f_1_g <- f_1_b <- numeric(N)
f_2_x <- f_2_y <- f_2_g <- f_2_b <- numeric(N)
f_h <- numeric(N)
for (i in 1 : N){
  x <- mvrnorm(nx, m_X, cov_X)
  y <- mvrnorm(ny, m_Y, cov_Y)
  order_res_x <- split(split_dat = x, cutting = cut_dat_x, permutation = perm_dat_x, pvertices = plot_vertices,
                       p = F)
  count_res_x <- counting_blocks(order_res_x, y)
  order_res_y <- split(split_dat = y, cutting = cut_dat_y, permutation = perm_dat_y, pvertices = plot_vertices,
                       p = F)
  count_res_y <- counting_blocks(order_res_y, x)
  ############################# Empty old #################################  
  res1 <- test_block_empty(count_res_x, Prob_x)
  res2 <- test_block_empty(count_res_y, Prob_y)
  f_old_x[i] <- res1[[1]]
  f_old_y[i] <- res2[[1]]
  f_old_g[i] <- test_good(p1 = res1[[2]], p2 = res2[[2]])
  f_old_b[i] <- test_bonf(p1 = res1[[2]], p2 = res2[[2]])
  ############################# K = 0 #################################      
  res3 <- test_block(count_res_x, test_k = 0)
  res4 <- test_block(count_res_y, test_k = 0)
  f_0_x[i] <- res3
  f_0_y[i] <- res4
  f_0_g[i] <- test_good(p1 = res3, p2 = res4)
  f_0_b[i] <- test_bonf(p1 = res3, p2 = res4)
  ############################# K = 1 #################################      
  res5 <- test_block(count_res_x, test_k = 1)
  res6 <- test_block(count_res_y, test_k = 1)
  f_1_x[i] <- res5
  f_1_y[i] <- res6
  f_1_g[i] <- test_good(p1 = res5, p2 = res6)
  f_1_b[i] <- test_bonf(p1 = res5, p2 = res6)
  ############################# K = 2 #################################      
  res7 <- test_block(count_res_x, test_k = 2)
  res8 <- test_block(count_res_y, test_k = 2)
  f_2_x[i] <- res7
  f_2_y[i] <- res8
  f_2_g[i] <- test_good(p1 = res7, p2 = res8)
  f_2_b[i] <- test_bonf(p1 = res7, p2 = res8)
  ############################# Hotelling's #################################   
  res9 <- hotelling.test(x, y)
  f_h[i] <- res9$pval
  print(i)
}

list("Simulation" = c("X mean" = m_X, "Y mean" = m_Y, "X num" = nx, "YX num" = ny),
     "X cov" = cov_X,
     "Y cov" = cov_Y,
     "Old Empty" = c("X" = sum(f_old_x)/N,  "Y" = sum(f_old_y)/N, "Bonf" = sum(f_old_b < 0.05)/N,  "Good" = sum(f_old_g < 0.05)/N), 
     "K = 0" = c("X" = sum(f_0_x < 0.05)/N, "Y" = sum(f_0_y < 0.05)/N, "Bonf" = sum(f_0_b < 0.05)/N, "Good" = sum(f_0_g < 0.05)/N),
     "K = 1" = c("X" = sum(f_1_x < 0.05)/N, "Y" = sum(f_1_y < 0.05)/N, "Bonf" = sum(f_1_b < 0.05)/N, "Good" = sum(f_1_g < 0.05)/N),
     "K = 2" = c("X" = sum(f_2_x < 0.05)/N, "Y" = sum(f_2_y < 0.05)/N, "Bonf" = sum(f_2_b < 0.05)/N, "Good" = sum(f_2_g < 0.05)/N), 
     "Hotelling" = sum(f_h < 0.05)/N)

