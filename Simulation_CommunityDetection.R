rm(list=ls())
library(irlba)
library(aricode)  #NMI
library(mclust)  #ARI
library(RSpectra)
library(tidyverse)
library(clusterSim)
library(orca)
library(RColorBrewer)

source("HSC_function.R")


#Case 1
N<-400
K<-2
R<-2

p<-0.2
lambda<-0.5
a<-0.001
b<-1

m1<-1.2
m2<-1
sd<-1

Z<-matrix(0,N,K)

for (j in 1:K){
  Z[((j-1)*N/K+1):(j*N/K),j]<-1
}
B<-p*lambda*diag(K)+p*(1-lambda)*rep(1,K)%*%t(rep(1,K))
M<-(m1-m2)*diag(K)+m2*rep(1,K)%*%t(rep(1,K))
EX<-Z%*%M
EA<-Z%*%B%*%t(Z)

true_labels <- rep(1:K, each = N/K)
times<-2
result0<-replicate(times,simul(N,K,R,p,lambda,a,b,m1,m2,sd,Z,B,M,EX,EA,true_labels))
mean_result<-apply(result0,1,mean)
mean_result <- data.frame(t(mean_result))

MR_result<-data.frame(HSCwithCov=mean_result[[12]],
                      CASC=mean_result[[9]],
                      EdgeSC=mean_result[[3]],
                      MotifSC=mean_result[[4]],
                      CovC=mean_result[[1]])

NMI_result<-data.frame(HSCwithCov=mean_result[[27]],
                       CASC=mean_result[[24]],
                       EdgeSC=mean_result[[18]],
                       MotifSC=mean_result[[19]],
                       CovC=mean_result[[16]])

ARI_result<-data.frame(HSCwithCov=mean_result[[40]],
                       CASC=mean_result[[37]],
                       EdgeSC=mean_result[[31]],
                       MotifSC=mean_result[[32]],
                       CovC=mean_result[[29]])

#Case 2
N<-400
K<-2
R<-2

p<-0.18
lambda<-0.18
a<-0.001
b<-1

m1<-1.2
m2<-0.2
sd<-0.4

Z<-matrix(0,N,K)

for (j in 1:K){
  Z[((j-1)*N/K+1):(j*N/K),j]<-1
}
B<-p*lambda*diag(K)+p*(1-lambda)*rep(1,K)%*%t(rep(1,K))
M<-(m1-m2)*diag(K)+m2*rep(1,K)%*%t(rep(1,K))
EX<-Z%*%M
EA<-Z%*%B%*%t(Z)

true_labels <- rep(1:K, each = N/K)
times<-100

result0<-replicate(times,simul(N,K,R,p,lambda,a,b,m1,m2,sd,Z,B,M,EX,EA,true_labels))
mean_result<-apply(result0,1,mean)
mean_result


#Case 3
N<-600
K<-4
R<-4

p<-0.2
lambda<-0.4
a<-0.1
b<-1.1

m1<-1.5
m2<-0.5
sd<-0.6

R.seq <-seq(1,10,1)
mymat <- matrix(data=NA,nrow=length(R.seq),ncol=38)
i=1
for (R in R.seq){
  Z<-matrix(0,N,K)
  
  for (j in 1:K){
    Z[((j-1)*N/K+1):(j*N/K),j]<-1
  }
  B<-p*lambda*diag(K)+p*(1-lambda)*rep(1,K)%*%t(rep(1,K))
  M<-(m1-m2)*diag(K)+m2*rep(1,K)%*%t(rep(1,K))
  EX<-Z%*%M
  EA<-Z%*%B%*%t(Z)
  
  true_labels <- rep(1:K, each = N/K)
  times<-100
  result0<-replicate(times,simul(N,K,R,p,lambda,a,b,m1,m2,sd,Z,B,M,EX,EA,true_labels))
  mean_result<-apply(result0,1,mean)
  mymat[i,] <- c(R,mean_result)
  print(mymat[i,])
  i<-i+1
}

# Case 4: Simulation for complex case 
# (B, Cov, Initialization of K-means)

N<-600
K<-3
R<-3

p<-0.4
lambda<-0.4
#a<-0.1
#b<-1.1

m1<-1.2
m2<-0.6
sd<-0.4

Z<-matrix(0,N,K)

for (j in 1:K){
  Z[((j-1)*N/K+1):(j*N/K),j]<-1
}

# adjMat
adjMat <- generate_adjMat_complex(N,K,p,B1=0.6,B2=0.9,lambda,a,b,Z)

# Generate different types of covariate matrices
covMat_spherical <- generate_covMat(N=N, K=K, R=R, Z=Z, m1=m1, m2=m2, sd=sd,type = "spherical")
covMat_non_spherical <- generate_covMat(N=N, K=K, R=R, Z=Z, m1=m1, m2=m2, sd=sd,type = "non-spherical", var_ratio = 2)
covMat_moons <- generate_covMat(N=N, K=K, R=R, Z=Z, m1=m1, m2=m2, sd=sd,type = "non-convex", 
                               shape_params = list(shape_type = "moons", noise = 0.1))
covMat_circles <- generate_covMat(N=N, K=K, R=R, Z=Z, m1=m1, m2=m2, sd=sd,type = "non-convex",
                                 shape_params = list(shape_type = "circles", noise = 0.05))



# Clustering
L_spherical <- simul_complexCov_enhanced(adjMat, covMat_spherical$covMat, N, K, criterion = "wcss", clust_method = "kmeans",n_init = 20)
L_non_spherical <- simul_complexCov_enhanced(adjMat, covMat_non_spherical$covMat, N, K, criterion = "wcss", clust_method = "kmeans",n_init = 20)
L_moons <- simul_complexCov_enhanced(adjMat, covMat_moons$covMat, N, K, criterion = "wcss", clust_method = "kmeans",n_init = 20)
L_circles <- simul_complexCov_enhanced(adjMat, covMat_circles$covMat, N, K, criterion = "wcss", clust_method = "kmeans",n_init = 20)

result_spherical <- specClust_enhanced(L_spherical$Lmat, K, nIter = 20, verbose = TRUE, clust_method = "kmeans", n_init = 20)
result_non_spherical <- specClust_enhanced(L_non_spherical$Lmat, K, nIter = 20, verbose = TRUE, clust_method = "kmeans", n_init = 20)
result_moons <- specClust_enhanced(L_moons$Lmat, K, nIter = 20, verbose = TRUE, clust_method = "kmeans", n_init = 20)
result_circles <- specClust_enhanced(L_circles$Lmat, K, nIter = 20, verbose = TRUE, clust_method = "kmeans", n_init = 20)

# Visualize Covariates
plot1_cov <- visualize_embeddings(covMat_spherical$covMat, covMat_spherical$true_labels, "Spherical Gaussian")
plot2_cov <- visualize_embeddings(covMat_non_spherical$covMat, covMat_non_spherical$true_labels, "Non-Spherical Gaussian") %>% 
  plotly::layout(scene = list(camera = list(eye = list(x = 1.25, y = 1.25, z = 0.8))))
plot3_cov <- visualize_embeddings(covMat_moons$covMat, covMat_moons$true_labels, "Moon-shaped") %>% 
  plotly::layout(scene = list(camera = list(eye = list(x = 0.25, y = 0.25, z = 2))))
plot4_cov <- visualize_embeddings(covMat_circles$covMat, covMat_circles$true_labels, "Circle-shaped") %>%
  plotly::layout(scene = list(camera = list(eye = list(x = 0.25, y = 0.25, z = 2))))

# Visualize embeddings of the similarity matrix
plot1_L <- visualize_embeddings(result_spherical$eig_vector, covMat_spherical$true_labels, "Spherical Gaussian") %>%
  plotly::layout(scene = list(camera = list(eye = list(x = 2, y = 0.5, z = 0.5))))
plot2_L <- visualize_embeddings(result_non_spherical$eig_vector, covMat_non_spherical$true_labels, "Non-Spherical Gaussian") %>%
  plotly::layout(scene = list(camera = list(eye = list(x = 2, y = 0.5, z = 0.5))))
plot3_L <- visualize_embeddings(result_moons$eig_vector, covMat_moons$true_labels, "Moon-shaped") %>%
  plotly::layout(scene = list(camera = list(eye = list(x = 2, y = 0.5, z = 0.5))))
plot4_L <- visualize_embeddings(result_circles$eig_vector, covMat_circles$true_labels, "Circle-shaped") %>%
  plotly::layout(scene = list(camera = list(eye = list(x = 2, y = 0.5, z = 0.5))))

plot_files <- c(
  "result/spherical_gaussian_cov.png",
  "result/non_spherical_gaussian_cov.png", 
  "result/moon_shaped_cov.png",
  "result/circle_shaped_cov.png",
  "result/spherical_gaussian_embed.png",
  "result/non_spherical_gaussian_embed.png",
  "result/moon_shaped_embed.png",
  "result/circle_shaped_embed.png"
)
# Read all PNG files and convert to grobs
plot_grobs <- list()
for (file in plot_files) {
  if (file.exists(file)) {
    plot_grobs[[basename(file)]] <- grid::rasterGrob(png::readPNG(file))
  }
}

# Create the combined plot with proper row and column labels
combined_plot <- gridExtra::grid.arrange(
  # First row: Covariates label
  grid::textGrob(
    "Covariates",
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  ),
  
  # Second row: Column labels for covariate plots
  gridExtra::arrangeGrob(
    grid::textGrob("Spherical Gaussian", gp = grid::gpar(fontsize = 16)),
    grid::textGrob("Non-Spherical Gaussian", gp = grid::gpar(fontsize = 16)),
    grid::textGrob("Moon-shaped", gp = grid::gpar(fontsize = 16)),
    grid::textGrob("Circle-shaped", gp = grid::gpar(fontsize = 16)),
    nrow = 1,
    ncol = 4
  ),
  
  # Third row: Covariate plots
  gridExtra::arrangeGrob(
    plot_grobs[["spherical_gaussian_cov.png"]],
    plot_grobs[["non_spherical_gaussian_cov.png"]],
    plot_grobs[["moon_shaped_cov.png"]],
    plot_grobs[["circle_shaped_cov.png"]],
    nrow = 1,
    ncol = 4
  ),
  
  # Fourth row: Spectral Embeddings label
  grid::textGrob(
    "Spectral Embeddings",
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  ),
  
  # Fifth row: Column labels for embedding plots
  gridExtra::arrangeGrob(
    grid::textGrob("Spherical Gaussian", gp = grid::gpar(fontsize = 16)),
    grid::textGrob("Non-Spherical Gaussian", gp = grid::gpar(fontsize = 16)),
    grid::textGrob("Moon-shaped", gp = grid::gpar(fontsize = 16)),
    grid::textGrob("Circle-shaped", gp = grid::gpar(fontsize = 16)),
    nrow = 1,
    ncol = 4
  ),
  
  # Sixth row: Embedding plots
  gridExtra::arrangeGrob(
    plot_grobs[["spherical_gaussian_embed.png"]],
    plot_grobs[["non_spherical_gaussian_embed.png"]],
    plot_grobs[["moon_shaped_embed.png"]],
    plot_grobs[["circle_shaped_embed.png"]],
    nrow = 1,
    ncol = 4
  ),
  
  nrow = 6,
  heights = c(0.5, 0.5, 4, 0.5, 0.5, 4)  # Adjust heights for labels vs plots

  )


# Save combined plot
combined_filename <- "plot/covariate_embedding_comparison.png"
ggplot2::ggsave(
  combined_filename, 
  combined_plot, 
  width = 12, 
  height = 7,  # Slightly taller to accommodate extra labels
  dpi = 300,
  bg = "white"
)


##
run_comprehensive_simulation <- function(n_sim = 100) {
  # Fixed parameters
  N <- 600
  K <- 3
  R <- 3
  p <- 0.4
  lambda <- 0.4
  a <- 0.1
  b <- 1.1
  m1 <- 1.2
  m2 <- 0.6
  sd_val <- 0.4
  
  # Define all combinations to test
  covMat_types <- c("spherical", "non-spherical", "moons", "circles")
  criterion_methods <- list(
    c("wcss", "kmeans"),
    c("silhouette", "silhouette"), 
    c("dbindex", "dbindex")
  )
  
  # Initialize results dataframe
  results_list <- list()
  
  # Create membership matrix Z
  Z <- matrix(0, N, K)
  for (j in 1:K) {
    Z[((j-1)*N/K + 1):(j*N/K), j] <- 1
  }
  
  # Loop through all combinations
  combo_counter <- 1
  for (cov_type in covMat_types) {
    for (crit_method in criterion_methods) {
      criterion <- crit_method[1]
      clust_method <- crit_method[2]
      
      cat("Running combination:", combo_counter, "/12\n")
      cat("Covariate type:", cov_type, "| Criterion:", criterion, 
          "| Method:", clust_method, "\n")
      
      # Initialize storage for metrics across simulations
      all_metrics <- matrix(0, nrow = n_sim, ncol = 15)
      
      # Run n_sim simulations
      for (sim in 1:n_sim) {
        # Generate new adjacency matrix for each simulation
        adjMat <- generate_adjMat_complex(N, K, p, B1 = 0.6, B2 = 0.9, lambda, a, b, Z)
        
        # Generate covariate matrix based on current type
        covMat_info <- switch(cov_type,
                         "spherical" = generate_covMat(N = N, K = K, R = R, Z = Z, m1 = m1, m2 = m2, 
                                                       sd = sd_val, type = "spherical"),
                         "non-spherical" = generate_covMat(N = N, K = K, R = R, Z = Z, m1 = m1, m2 = m2, 
                                                           sd = sd_val, type = "non-spherical", var_ratio = 2),
                         "moons" = generate_covMat(N = N, K = K, R = R, Z = Z, m1 = m1, m2 = m2, 
                                                   sd = sd_val, type = "non-convex",
                                                   shape_params = list(shape_type = "moons", noise = 0.1)),
                         "circles" = generate_covMat(N = N, K = K, R = R, Z = Z, m1 = m1, m2 = m2, 
                                                     sd = sd_val, type = "non-convex",
                                                     shape_params = list(shape_type = "circles", noise = 0.05))
        )
        
        # Run simulation with current parameters
        result <- simul_complexCov_enhanced(
          adjMat, covMat_info$covMat, N, K,
          criterion = criterion,
          clust_method = clust_method,
          n_init = 20
        )
        
        # Store metrics
        all_metrics[sim, ] <- result$metric
      }
      
      # Calculate average metrics across simulations
      avg_metrics <- colMeans(all_metrics)
      
      # Create result row
      result_row <- data.frame(
        covMat_type = cov_type,
        criterion = criterion,
        clust_method = clust_method,
        X_only_error = avg_metrics[1],
        We_error = avg_metrics[2],
        Wm_error = avg_metrics[3],
        casc_We_X_error = avg_metrics[4],
        Mixed_W_X_error = avg_metrics[5],
        X_only_nmi = avg_metrics[6],
        We_nmi = avg_metrics[7],
        Wm_nmi = avg_metrics[8],
        casc_We_X_nmi = avg_metrics[9],
        Mixed_W_X_nmi = avg_metrics[10],
        X_only_ari = avg_metrics[11],
        We_ari = avg_metrics[12],
        Wm_ari = avg_metrics[13],
        casc_We_X_ari = avg_metrics[14],
        Mixed_W_X_ari = avg_metrics[15]
      )
      
      results_list[[combo_counter]] <- result_row
      combo_counter <- combo_counter + 1
      
      cat("Completed", n_sim, "simulations for current combination\n\n")
    }
  }
  
  # Combine all results
  final_results <- do.call(rbind, results_list)
  return(final_results)
}

# Run the comprehensive simulation with 100 replicates per combination
set.seed(123)  # For reproducibility
simulation_results <- run_comprehensive_simulation(n_sim = 100)

# Display results
print(simulation_results)

# Save results to RData file
save(simulation_results, file = "complex_simulation_results.RData")
write.csv(simulation_results, "comprehensive_simulation_results.csv", row.names = FALSE)

