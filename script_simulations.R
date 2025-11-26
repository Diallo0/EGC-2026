source("./main.R")


#----------------------------------- data simulées  -------------------------------------------------------------------
t = seq(0, 1, length.out = 20)

data_condI = sim_condition_I(n_per_clus = 100, t, f_list)
data_condII = sim_condition_II(n_per_clus = 100, t, f_list)
data_condIII = sim_condition_III(n_per_clus = 100, t, f_list)


#------------------------ Lissage MCP ------------------------------------------------

data_condI_liss = liss(data_condI[, -c(1,2)])

data_condII_liss = liss(data_condII[, -c(1,2)])

data_condIII_liss = liss(data_condIII[, -c(1,2)])

##--------------- visuels courbes reconstruites --------------------------------

fine_grid_t = seq(range(t)[1], range(t)[2] , length.out = 1000)

gg_discret_continue(data_condI, data_condI_liss, fine_grid_t)
gg_discret_continue(data_condII, data_condII_liss, fine_grid_t)
gg_discret_continue(data_condIII, data_condIII_liss, fine_grid_t)



##------------------------ visuels des fonctions et de leurs derivées (Figure 1, 2, 3) ----------------------------------------

gg_f_f_prime(data_condI, data_condI_liss, fine_grid_t)
gg_f_f_prime(data_condII, data_condII_liss, fine_grid_t)
gg_f_f_prime(data_condIII, data_condIII_liss, fine_grid_t)


# -------------------------- Simulation des scénarios --------------------------------------------------------------------

# simulation
Simulation = function(t){list(
  sim_condition_I(n_per_clus = 100, t, f_list),
  sim_condition_II(n_per_clus = 100, t, f_list),
  sim_condition_III(n_per_clus = 100, t, f_list)
)}




# ----------- k-means Riemann pour un seul scénario ----------------
kmeans_clus_Riemann_sim = function(list_dist, liss_list, fine_grid, class_true){
  # list_dist : liste des distance dans un scénarios pour les oméga : calcul de silhouette
  # liss_list : list des données lissées pour ce scénario
  
  # X et X_deriv sur grille fine
  X = t(sapply(1:length(liss_list), function(i) eval.fd(fine_grid, liss_list[[i]])))
  X_deriv = t(sapply(1:length(liss_list), function(i) eval.fd(fine_grid, liss_list[[i]], Lfdobj = 1)))
  
  omegas = seq(0, 1, length.out = length(list_dist))
  
  classes = lapply(seq_along(omegas), function(j){
    k_means_dw(
      X = X, 
      X_deriv = X_deriv, 
      K = length(unique(class_true)), 
      method = "riemann", 
      omega = omegas[j], 
      fine_grid = fine_grid
    )$cluster
  })
  
  # ARI
  ARI = sapply(seq_along(classes), function(j){
    mclust::adjustedRandIndex(classes[[j]], class_true)
  })
  
  # SIL
  SIL = sapply(seq_along(classes), function(j){
    cluster::silhouette(classes[[j]], dist = list_dist[[j]])[, 3] %>% mean()
  })
  
  return(list(ARI = ARI, SIL = SIL))
}

# ----------- k-means ACPF pour un seul scénario -------------------
kmeans_clus_ACPF_sim = function(list_dist, pca_res, class_true){
  # list_dist : pour silhouette 
  # pca_res : l'acpf
  omegas = seq(0, 1, length.out = length(list_dist))
  
  classes = lapply(seq_along(omegas), function(j){
    k_means_dw(
      X = pca_res$coords, 
      X_deriv = pca_res$coords_deriv, 
      K = length(unique(class_true)), 
      method = "ACPF", 
      omega = omegas[j]
    )$cluster
  })
  
  # ARI
  ARI = sapply(seq_along(classes), function(j){
    mclust::adjustedRandIndex(classes[[j]], class_true)
  })
  
  # SIL
  SIL = sapply(seq_along(classes), function(j){
    cluster::silhouette(classes[[j]], dist = list_dist[[j]])[, 3] %>% mean()
  })
  
  return(list(ARI = ARI, SIL = SIL))
}



# ------------------------- fonction d'execution pour une liste de 3 tableaux des scénarios ------------------------------------

Execution_distance_pond = function(t, sim){
  # Execution des différentes appproches pour les trois scénarios 
  # sim : est une liste des tableaux de données des 3 scénarios de simulations
  
  liss = lapply(1:3, function(i) liss(sim[[i]][, -c(1,2)]))
  
   
  #  ---------  CAH (Ward) : Approche de riemman -----------------------------------
  fine_grid_t = seq(range(t)[1], range(t)[2] , length.out = 1000)
  ## DO
  dist_D0 = lapply(1:3, function(i) calculate_D0_matrix_parallel(liss[[i]], fine_grid_t))
  dist_D1 = lapply(1:3, function(i) calculate_D1_matrix_parallel(liss[[i]], fine_grid_t))
  dist_Dp_condI = lapply(seq(0.2, 0.8, by = 0.2), function(i) calculate_Dp_matrix_byD0D1(dist_D0[[1]], 
                                              dist_D1[[1]], i))
  dist_Dp_condII = lapply(seq(0.2, 0.8, by = 0.2), function(i) calculate_Dp_matrix_byD0D1(dist_D0[[2]], 
                                                                                     dist_D1[[2]], i))
  dist_Dp_condIII = lapply(seq(0.2, 0.8, by = 0.2), function(i) calculate_Dp_matrix_byD0D1(dist_D0[[3]], 
                                                                                     dist_D1[[3]], i))
 
  # clustering via hclust avec ward 
  class_D0 = lapply(1:3, function(i) cutree(hclust(dist_D0[[i]]%>%as.dist(), method = "ward.D2"), k = 4))
  class_D1 = lapply(1:3, function(i) cutree(hclust(dist_D1[[i]]%>%as.dist(), method = "ward.D2"), k = 4))
  # Dp
  class_Dp_condI = lapply(1:length(dist_Dp_condI), function(i) cutree(hclust(dist_Dp_condI[[i]]%>%as.dist(),
                                                                             method = "ward.D2"), k = 4))
  class_Dp_condII = lapply(1:length(dist_Dp_condII), function(i) cutree(hclust(dist_Dp_condII[[i]]%>%as.dist(),
                                                                             method = "ward.D2"), k = 4))
  class_Dp_condIII = lapply(1:length(dist_Dp_condIII), function(i) cutree(hclust(dist_Dp_condIII[[i]]%>%as.dist(),
                                                                            method = "ward.D2"), k = 4))
  
  ## ------------------------ silhouette (riemann)-----------------------------------
  sil_dist_D0  = lapply(1:length(class_D0), function(i) silhouette(class_D0[[i]], dist = dist_D0[[i]])[, 3]%>%mean())
  sil_dist_D1  = lapply(1:length(class_D1), function(i) silhouette(class_D1[[i]], dist = dist_D1[[i]])[, 3]%>%mean())
  sil_Dp_condI = lapply(1:length(class_Dp_condI), function(i) silhouette(class_Dp_condI[[i]], dist = dist_Dp_condI[[i]])[, 3]%>%mean())
  sil_Dp_condII = lapply(1:length(class_Dp_condII), function(i) silhouette(class_Dp_condII[[i]], dist = dist_Dp_condII[[i]])[, 3]%>%mean())
  sil_Dp_condIII = lapply(1:length(class_Dp_condIII), function(i) silhouette(class_Dp_condIII[[i]], dist = dist_Dp_condIII[[i]])[, 3]%>%mean())
  
  sil_Riemann  = list(sil_dist_D0 = sil_dist_D0, sil_dist_D1 = sil_dist_D1, 
                      sil_Dp_condI = sil_Dp_condI, sil_Dp_condII = sil_Dp_condII, 
                      sil_Dp_condIII = sil_Dp_condIII)
  
  ## ------------------------------------ARI riemann----------------------------------
 
  ari_dist_D0 = lapply(1:length(class_D0), function (i) round(mclust::adjustedRandIndex(class_D0[[i]], sim[[1]]$class), 2)) # meme class
  ari_dist_D1 = lapply(1:length(class_D1), function (i) round(mclust::adjustedRandIndex(class_D1[[i]], sim[[1]]$class), 2))
  
  ari_dist_Dp_condI = lapply(1:length(class_Dp_condI), function (i) round(mclust::adjustedRandIndex(class_Dp_condI[[i]],
                                                                                                    sim[[1]]$class), 2))
  ari_dist_Dp_condII = lapply(1:length(class_Dp_condII), function (i) round(mclust::adjustedRandIndex(class_Dp_condII[[i]],
                                                                                                    sim[[1]]$class), 2))
  ari_dist_Dp_condIII = lapply(1:length(class_Dp_condIII), function (i) round(mclust::adjustedRandIndex(class_Dp_condIII[[i]],
                                                                                                    sim[[1]]$class), 2))
  
  Ari_Riemann = list(ari_dist_D0 = ari_dist_D0, ari_dist_D1 = ari_dist_D1, ari_dist_Dp_condI = ari_dist_Dp_condI, 
                     ari_dist_Dp_condII = ari_dist_Dp_condII, ari_dist_Dp_condIII = ari_dist_Dp_condIII)
  
  
  #  ------------------- CAH (Ward) :  ACPF --------------------------------------------------
  
  pca = lapply(1:3, function(i) pca_sans_reaffect(liss[[i]], t))
  
  dist_D0_sans_reaffect = lapply(1:3, function(i) dist(pca[[i]]$coords))
  dist_D1_sans_reaffect = lapply(1:3, function(i) dist(pca[[i]]$coords_deriv))
  
  # Dp sans reaffect
  dist_Dp_condI_sans_reaffect = lapply(seq(0.2, 0.8, by = 0.2), function(i) pca_Dp(pca[[1]], i))
  dist_Dp_condII_sans_reaffect = lapply(seq(0.2, 0.8, by = 0.2), function(i) pca_Dp(pca[[2]], i))
  dist_Dp_condIII_sans_reaffect = lapply(seq(0.2, 0.8, by = 0.2), function(i) pca_Dp(pca[[3]], i))
  
  
  class_D0_sans_reaffect = lapply(1:3, function(i) cutree(hclust(dist_D0_sans_reaffect[[i]], method = "ward.D2"), k = 4))
  class_D1_sans_reaffect = lapply(1:3, function(i) cutree(hclust(dist_D1_sans_reaffect[[i]], method = "ward.D2"), k = 4))
  
  
  class_Dp_sans_reaffect_condI = lapply(1:length(dist_Dp_condI_sans_reaffect), function(i) cutree(hclust(dist_Dp_condI_sans_reaffect[[i]],
                                                                             method = "ward.D2"), k = 4))
  class_Dp_sans_reaffect_condII = lapply(1:length(dist_Dp_condII_sans_reaffect), function(i) cutree(hclust(dist_Dp_condII_sans_reaffect[[i]],
                                                                               method = "ward.D2"), k = 4))
  class_Dp_sans_reaffect_condIII = lapply(1:length(dist_Dp_condIII_sans_reaffect), function(i) cutree(hclust(dist_Dp_condIII_sans_reaffect[[i]],
                                                                                               method = "ward.D2"), k = 4))
  
  ## ------------------------------------------ silhouette ACPF --------------------------------
  sil_ACPF_dist_D0  = lapply(1:length(class_D0_sans_reaffect), function(i) silhouette(class_D0_sans_reaffect[[i]], 
                                                                                      dist = dist_D0_sans_reaffect[[i]])[, 3]%>%mean())
  sil_ACPF_dist_D1  = lapply(1:length(class_D1_sans_reaffect), function(i) silhouette(class_D1_sans_reaffect[[i]], 
                                                                                      dist = dist_D1_sans_reaffect[[i]])[, 3]%>%mean())
  sil_ACPF_Dp_condI = lapply(1:length(class_Dp_sans_reaffect_condI), function(i) silhouette(class_Dp_sans_reaffect_condI[[i]], 
                                                                              dist = dist_Dp_condI_sans_reaffect[[i]])[, 3]%>%mean())
  sil_ACPF_Dp_condII = lapply(1:length(class_Dp_sans_reaffect_condII), function(i) silhouette(class_Dp_sans_reaffect_condII[[i]],
                                                                                dist = dist_Dp_condII_sans_reaffect[[i]])[, 3]%>%mean())
  sil_ACPF_Dp_condIII = lapply(1:length(class_Dp_sans_reaffect_condIII), function(i) silhouette(class_Dp_sans_reaffect_condIII[[i]],
                                                                                  dist = dist_Dp_condIII_sans_reaffect[[i]])[, 3]%>%mean())
  sil_ACPF = list(sil_ACPF_dist_D0 = sil_ACPF_dist_D0, sil_ACPF_dist_D1 = sil_ACPF_dist_D1, 
                  sil_ACPF_Dp_condI = sil_ACPF_Dp_condI, sil_ACPF_Dp_condII = sil_ACPF_Dp_condII, 
                  sil_ACPF_Dp_condIII = sil_ACPF_Dp_condIII)
  
 ## ------------------------------------  ARI ACPF ------------------------------------------------------
  
  ari_dist_D0_sans_reaffect = lapply(1:length(class_D0_sans_reaffect), function (i) round(mclust::adjustedRandIndex(class_D0_sans_reaffect[[i]], sim[[1]]$class), 2))# meme class
  ari_dist_D1_sans_reaffect = lapply(1:length(class_D1_sans_reaffect), function (i) round(mclust::adjustedRandIndex(class_D1_sans_reaffect[[i]], sim[[1]]$class), 2))
  
  ari_dist_Dp_sans_reaffect_condI = lapply(1:length(class_Dp_sans_reaffect_condI), function (i) round(mclust::adjustedRandIndex(class_Dp_sans_reaffect_condI[[i]],
                                                                                                                  sim[[1]]$class), 2))
  ari_dist_Dp_sans_reaffect_condII = lapply(1:length(class_Dp_sans_reaffect_condII), function (i) round(mclust::adjustedRandIndex(class_Dp_sans_reaffect_condII[[i]],
                                                                                                                   sim[[1]]$class), 2))
  ari_dist_Dp_sans_reaffect_condIII = lapply(1:length(class_Dp_sans_reaffect_condIII), function (i) round(mclust::adjustedRandIndex(class_Dp_sans_reaffect_condIII[[i]],
                                                                                                                   sim[[1]]$class), 2))
  
  
  Ari_ACPF = list(ari_dist_D0_sans_reaffect = ari_dist_D0_sans_reaffect, ari_dist_D1_sans_reaffect = ari_dist_D1_sans_reaffect, 
                  ari_dist_Dp_sans_reaffect_condI = ari_dist_Dp_sans_reaffect_condI, 
                  ari_dist_Dp_sans_reaffect_condII = ari_dist_Dp_sans_reaffect_condII, 
                  ari_dist_Dp_sans_reaffect_condIII = ari_dist_Dp_sans_reaffect_condIII)


 # --------------------------------------------- CAH : BASELINE --------------------------------------------------------------------------
 
  coefs_liss_condI = sapply(1:length(liss[[1]]), function(i) liss[[1]][[i]]$coefs)%>%t()
  coefs_liss_condII = sapply(1:length(liss[[2]]), function(i) liss[[2]][[i]]$coefs)%>%t()
  coefs_liss_condIII = sapply(1:length(liss[[3]]), function(i) liss[[3]][[i]]$coefs)%>%t()
  # clustering
  baseline_clus_condI = cutree(hclust(dist(coefs_liss_condI), method = "ward.D2"), k = 4)
  baseline_clus_condII = cutree(hclust(dist(coefs_liss_condII), method = "ward.D2"), k = 4)
  baseline_clus_condIII = cutree(hclust(dist(coefs_liss_condIII), method = "ward.D2"), k = 4)
  
  
  # silhouette baseline
  sil_baseline_condI = silhouette(baseline_clus_condI, dist(coefs_liss_condI))[, 3]%>%mean()
  sil_baseline_condII = silhouette(baseline_clus_condII, dist(coefs_liss_condII))[, 3]%>%mean()
  sil_baseline_condIII = silhouette(baseline_clus_condIII, dist(coefs_liss_condIII))[, 3]%>%mean()
  
  sil_baseline = list(sil_baseline_condI = sil_baseline_condI, sil_baseline_condII = sil_baseline_condII, 
                      sil_baseline_condIII = sil_baseline_condIII)
  
  # ARI baseline
  Ari_baseline_condI = mclust::adjustedRandIndex(baseline_clus_condI, sim[[1]]$class)
  Ari_baseline_condII = mclust::adjustedRandIndex(baseline_clus_condII, sim[[1]]$class)
  Ari_baseline_condIII = mclust::adjustedRandIndex(baseline_clus_condIII, sim[[1]]$class)
  
  Ari_baseline = list(Ari_baseline_condI = Ari_baseline_condI, Ari_baseline_condII = Ari_baseline_condII, 
                      Ari_baseline_condIII = Ari_baseline_condIII)
  
  
  # -------------------------- k-means fonctionnelle -------------------------------------
  
  ## ----------------- Préparer les listes de distances pour k-means Riemann --------------
  omega_seq = seq(0, 1, length.out = 6)
  # pour chaque condition, construire une liste de 6 distances dans le même ordre que omega_seq
  
  # Scénario I
  list_dist_km_Riemann_condI = list()
  list_dist_km_Riemann_condI[[1]] = dist_D0[[1]]         # omega = 0
  list_dist_km_Riemann_condI[[6]] = dist_D1[[1]]         # omega = 1
  for(j in 1:length(dist_Dp_condI)){
    list_dist_km_Riemann_condI[[1 + j]] = dist_Dp_condI[[j]]
  }
  
  # Scénario II
  list_dist_km_Riemann_condII = list()
  list_dist_km_Riemann_condII[[1]] = dist_D0[[2]]
  list_dist_km_Riemann_condII[[6]] = dist_D1[[2]]
  for(j in 1:length(dist_Dp_condII)){
    list_dist_km_Riemann_condII[[1 + j]] = dist_Dp_condII[[j]]
  }
  
  # Scénario III
  list_dist_km_Riemann_condIII = list()
  list_dist_km_Riemann_condIII[[1]] = dist_D0[[3]]
  list_dist_km_Riemann_condIII[[6]] = dist_D1[[3]]
  for(j in 1:length(dist_Dp_condIII)){
    list_dist_km_Riemann_condIII[[1 + j]] = dist_Dp_condIII[[j]]
  }
  
  
  ## ------------ Préparer les listes de distances pour k-means ACPF ----------------

  # scénario I
  list_dist_km_ACPF_condI = list()
  list_dist_km_ACPF_condI[[1]] = dist_D0_sans_reaffect[[1]]      # omega = 0
  list_dist_km_ACPF_condI[[6]] = dist_D1_sans_reaffect[[1]]      # omega = 1
  for(j in 1:length(dist_Dp_condI_sans_reaffect)){
    list_dist_km_ACPF_condI[[1 + j]] = dist_Dp_condI_sans_reaffect[[j]]
  }
  
  # scénario II
  list_dist_km_ACPF_condII = list()
  list_dist_km_ACPF_condII[[1]] = dist_D0_sans_reaffect[[2]]
  list_dist_km_ACPF_condII[[6]] = dist_D1_sans_reaffect[[2]]
  for(j in 1:length(dist_Dp_condII_sans_reaffect)){
    list_dist_km_ACPF_condII[[1 + j]] = dist_Dp_condII_sans_reaffect[[j]]
  }
  
  # scénario III
  list_dist_km_ACPF_condIII = list()
  list_dist_km_ACPF_condIII[[1]] = dist_D0_sans_reaffect[[3]]
  list_dist_km_ACPF_condIII[[6]] = dist_D1_sans_reaffect[[3]]
  for(j in 1:length(dist_Dp_condIII_sans_reaffect)){
    list_dist_km_ACPF_condIII[[1 + j]] = dist_Dp_condIII_sans_reaffect[[j]]
  }
  
  
  ## ---------------- k-means : Riemann  -------------------
  km_Riemann_condI = kmeans_clus_Riemann_sim(
    list_dist = list_dist_km_Riemann_condI,
    liss_list = liss[[1]],
    fine_grid = fine_grid_t,
    class_true = sim[[1]]$class
  )
  
  km_Riemann_condII = kmeans_clus_Riemann_sim(
    list_dist = list_dist_km_Riemann_condII,
    liss_list = liss[[2]],
    fine_grid = fine_grid_t,
    class_true = sim[[2]]$class
  )
  
  km_Riemann_condIII = kmeans_clus_Riemann_sim(
    list_dist = list_dist_km_Riemann_condIII,
    liss_list = liss[[3]],
    fine_grid = fine_grid_t,
    class_true = sim[[3]]$class
  )
  
  Ari_kmeans_Riemann = list(
    condI = km_Riemann_condI$ARI,
    condII = km_Riemann_condII$ARI,
    condIII = km_Riemann_condIII$ARI
  )
  
  sil_kmeans_Riemann = list(
    condI = km_Riemann_condI$SIL,
    condII = km_Riemann_condII$SIL,
    condIII = km_Riemann_condIII$SIL
  )
  
  ##  ---------------- k-means : ACPF  ----------------------
  km_ACPF_condI = kmeans_clus_ACPF_sim(
    list_dist = list_dist_km_ACPF_condI,
    pca_res = pca[[1]],
    class_true = sim[[1]]$class
  )
  
  km_ACPF_condII = kmeans_clus_ACPF_sim(
    list_dist = list_dist_km_ACPF_condII,
    pca_res = pca[[2]],
    class_true = sim[[2]]$class
  )
  
  km_ACPF_condIII = kmeans_clus_ACPF_sim(
    list_dist = list_dist_km_ACPF_condIII,
    pca_res = pca[[3]],
    class_true = sim[[3]]$class
  )
  
  Ari_kmeans_ACPF = list(
    condI = km_ACPF_condI$ARI,
    condII = km_ACPF_condII$ARI,
    condIII = km_ACPF_condIII$ARI
  )
  
  sil_kmeans_ACPF = list(
    condI = km_ACPF_condI$SIL,
    condII = km_ACPF_condII$SIL,
    condIII = km_ACPF_condIII$SIL
  )
  
  # ----------------------------- PAM --------------------------------------------------
  
  ## ------------------------ PAM : Riemann -----------------------------------
  
  # D0
  pam_Riemann_D0 = lapply(1:3, function(i)
    pam_clus(dist_D0[[i]], class_true = sim[[i]]$class))
  # D1
  pam_Riemann_D1 = lapply(1:3, function(i)
    pam_clus(dist_D1[[i]], class_true = sim[[i]]$class))
  
  # Dp 
  pam_Riemann_Dp_condI   = pam_clus(dist_Dp_condI,   class_true = sim[[1]]$class)
  pam_Riemann_Dp_condII  = pam_clus(dist_Dp_condII,  class_true = sim[[2]]$class)
  pam_Riemann_Dp_condIII = pam_clus(dist_Dp_condIII, class_true = sim[[3]]$class)
  

  Ari_PAM_Riemann = list(
    D0 = lapply(pam_Riemann_D0,  function(x) x$ARI),
    D1 = lapply(pam_Riemann_D1,  function(x) x$ARI),
    Dp = list(
      condI   = pam_Riemann_Dp_condI$ARI,
      condII  = pam_Riemann_Dp_condII$ARI,
      condIII = pam_Riemann_Dp_condIII$ARI
    )
  )
  
  sil_PAM_Riemann = list(
    D0 = lapply(pam_Riemann_D0,  function(x) x$SIL),
    D1 = lapply(pam_Riemann_D1,  function(x) x$SIL),
    Dp = list(
      condI   = pam_Riemann_Dp_condI$SIL,
      condII  = pam_Riemann_Dp_condII$SIL,
      condIII = pam_Riemann_Dp_condIII$SIL
    )
  )
  
  
  ## ------------------------ PAM ACPF -----------------------------------
  
  # D0 (coords)
  pam_ACPF_D0 = lapply(1:3, function(i)
    pam_clus(dist_D0_sans_reaffect[[i]], class_true = sim[[i]]$class))
  
  # D1 (coords_deriv)
  pam_ACPF_D1 = lapply(1:3, function(i)
    pam_clus(dist_D1_sans_reaffect[[i]], class_true = sim[[i]]$class))
  
  # Dp PCA (listes)
  pam_ACPF_Dp_condI   = pam_clus(dist_Dp_condI_sans_reaffect,   class_true = sim[[1]]$class)
  pam_ACPF_Dp_condII  = pam_clus(dist_Dp_condII_sans_reaffect,  class_true = sim[[2]]$class)
  pam_ACPF_Dp_condIII = pam_clus(dist_Dp_condIII_sans_reaffect, class_true = sim[[3]]$class)
  
  Ari_PAM_ACPF = list(
    D0 = lapply(pam_ACPF_D0, function(x) x$ARI),
    D1 = lapply(pam_ACPF_D1, function(x) x$ARI),
    Dp = list(
      condI   = pam_ACPF_Dp_condI$ARI,
      condII  = pam_ACPF_Dp_condII$ARI,
      condIII = pam_ACPF_Dp_condIII$ARI
    )
  )
  
  sil_PAM_ACPF = list(
    D0 = lapply(pam_ACPF_D0, function(x) x$SIL),
    D1 = lapply(pam_ACPF_D1, function(x) x$SIL),
    Dp = list(
      condI   = pam_ACPF_Dp_condI$SIL,
      condII  = pam_ACPF_Dp_condII$SIL,
      condIII = pam_ACPF_Dp_condIII$SIL
    )
  )
  
  
  
  return(list(
    Ari_Riemann   = Ari_Riemann,
    Ari_ACPF      = Ari_ACPF,
    sil_ACPF      = sil_ACPF,
    sil_Riemann   = sil_Riemann,
    sil_baseline  = sil_baseline,
    Ari_baseline  = Ari_baseline,
    # résultats k-means (si tu as suivi l’étape précédente)
    Ari_kmeans_Riemann = Ari_kmeans_Riemann,
    sil_kmeans_Riemann = sil_kmeans_Riemann,
    Ari_kmeans_ACPF    = Ari_kmeans_ACPF,
    sil_kmeans_ACPF    = sil_kmeans_ACPF,
    # nouveaux résultats PAM
    Ari_PAM_Riemann = Ari_PAM_Riemann,
    sil_PAM_Riemann = sil_PAM_Riemann,
    Ari_PAM_ACPF    = Ari_PAM_ACPF,
    sil_PAM_ACPF    = sil_PAM_ACPF
  ))
  
  
  
}




f_arranger_ARI = function(ex_dist_pond){
  # arrange les resultats sous forme de matrices 3 x 6
  # lignes : 3 conditions de simulation
  # colonnes : omega = (0, 0.2, 0.4, 0.6, 0.8, 1)
  
  ## ---------------- CAH : Riemann ---------------------------------
  Ari_Riemann = matrix(NA, ncol = 6 , nrow = 3)
  Ari_Riemann[, 1] = ex_dist_pond$Ari_Riemann$ari_dist_D0 %>% unlist()
  Ari_Riemann[, 6] = ex_dist_pond$Ari_Riemann$ari_dist_D1 %>% unlist()
  Ari_Riemann[1, 2:5] = ex_dist_pond$Ari_Riemann$ari_dist_Dp_condI   %>% unlist()
  Ari_Riemann[2, 2:5] = ex_dist_pond$Ari_Riemann$ari_dist_Dp_condII  %>% unlist()
  Ari_Riemann[3, 2:5] = ex_dist_pond$Ari_Riemann$ari_dist_Dp_condIII %>% unlist()
  
  ## ---------------- CAH : ACPF ------------------------------------
  Ari_ACPF = matrix(NA, ncol = 6 , nrow = 3)
  Ari_ACPF[, 1] = ex_dist_pond$Ari_ACPF$ari_dist_D0_sans_reaffect %>% unlist()
  Ari_ACPF[, 6] = ex_dist_pond$Ari_ACPF$ari_dist_D1_sans_reaffect %>% unlist()
  Ari_ACPF[1, 2:5] = ex_dist_pond$Ari_ACPF$ari_dist_Dp_sans_reaffect_condI   %>% unlist()
  Ari_ACPF[2, 2:5] = ex_dist_pond$Ari_ACPF$ari_dist_Dp_sans_reaffect_condII  %>% unlist()
  Ari_ACPF[3, 2:5] = ex_dist_pond$Ari_ACPF$ari_dist_Dp_sans_reaffect_condIII %>% unlist()
  
  ## ---------------- k-means : Riemann -----------------------------
  Ari_kmeans_Riemann = rbind(
    ex_dist_pond$Ari_kmeans_Riemann$condI,
    ex_dist_pond$Ari_kmeans_Riemann$condII,
    ex_dist_pond$Ari_kmeans_Riemann$condIII
  )
  
  ## ---------------- k-means : ACPF -------------------------------
  Ari_kmeans_ACPF = rbind(
    ex_dist_pond$Ari_kmeans_ACPF$condI,
    ex_dist_pond$Ari_kmeans_ACPF$condII,
    ex_dist_pond$Ari_kmeans_ACPF$condIII
  )
  
  ## ---------------- PAM : Riemann --------------------------------
  Ari_PAM_Riemann = matrix(NA, ncol = 6, nrow = 3)
  Ari_PAM_Riemann[, 1] = ex_dist_pond$Ari_PAM_Riemann$D0 %>% unlist()
  Ari_PAM_Riemann[, 6] = ex_dist_pond$Ari_PAM_Riemann$D1 %>% unlist()
  Ari_PAM_Riemann[1, 2:5] = ex_dist_pond$Ari_PAM_Riemann$Dp$condI   %>% unlist()
  Ari_PAM_Riemann[2, 2:5] = ex_dist_pond$Ari_PAM_Riemann$Dp$condII  %>% unlist()
  Ari_PAM_Riemann[3, 2:5] = ex_dist_pond$Ari_PAM_Riemann$Dp$condIII %>% unlist()
  
  ## ---------------- PAM : ACPF -----------------------------------
  Ari_PAM_ACPF = matrix(NA, ncol = 6, nrow = 3)
  Ari_PAM_ACPF[, 1] = ex_dist_pond$Ari_PAM_ACPF$D0 %>% unlist()
  Ari_PAM_ACPF[, 6] = ex_dist_pond$Ari_PAM_ACPF$D1 %>% unlist()
  Ari_PAM_ACPF[1, 2:5] = ex_dist_pond$Ari_PAM_ACPF$Dp$condI   %>% unlist()
  Ari_PAM_ACPF[2, 2:5] = ex_dist_pond$Ari_PAM_ACPF$Dp$condII  %>% unlist()
  Ari_PAM_ACPF[3, 2:5] = ex_dist_pond$Ari_PAM_ACPF$Dp$condIII %>% unlist()
  
  return(list(
    Ari_Riemann        = Ari_Riemann,
    Ari_ACPF           = Ari_ACPF,
    Ari_kmeans_Riemann = Ari_kmeans_Riemann,
    Ari_kmeans_ACPF    = Ari_kmeans_ACPF,
    Ari_PAM_Riemann    = Ari_PAM_Riemann,
    Ari_PAM_ACPF       = Ari_PAM_ACPF
  ))
}




f_arranger_silhouette = function(ex_dist_pond){
  # lignes  : 3 conditions
  # colonnes: omega = (0, 0.2, 0.4, 0.6, 0.8, 1)
  
  ## --------------- CAH : Riemann ---------------------------------
  sil_Riemann = matrix(NA, ncol = 6, nrow = 3)
  sil_Riemann[, 1] = ex_dist_pond$sil_Riemann$sil_dist_D0      %>% unlist()
  sil_Riemann[, 6] = ex_dist_pond$sil_Riemann$sil_dist_D1      %>% unlist()
  sil_Riemann[1, 2:5] = ex_dist_pond$sil_Riemann$sil_Dp_condI   %>% unlist()
  sil_Riemann[2, 2:5] = ex_dist_pond$sil_Riemann$sil_Dp_condII  %>% unlist()
  sil_Riemann[3, 2:5] = ex_dist_pond$sil_Riemann$sil_Dp_condIII %>% unlist()
  
  ## --------------- CAH : ACPF ------------------------------------
  sil_ACPF = matrix(NA, ncol = 6, nrow = 3)
  sil_ACPF[, 1] = ex_dist_pond$sil_ACPF$sil_ACPF_dist_D0      %>% unlist()
  sil_ACPF[, 6] = ex_dist_pond$sil_ACPF$sil_ACPF_dist_D1      %>% unlist()
  sil_ACPF[1, 2:5] = ex_dist_pond$sil_ACPF$sil_ACPF_Dp_condI   %>% unlist()
  sil_ACPF[2, 2:5] = ex_dist_pond$sil_ACPF$sil_ACPF_Dp_condII  %>% unlist()
  sil_ACPF[3, 2:5] = ex_dist_pond$sil_ACPF$sil_ACPF_Dp_condIII %>% unlist()
  
  ## --------------- k-means : Riemann -----------------------------
  # sil_kmeans_Riemann$condI/II/III sont déjà des vecteurs de longueur 6
  sil_kmeans_Riemann = rbind(
    ex_dist_pond$sil_kmeans_Riemann$condI,
    ex_dist_pond$sil_kmeans_Riemann$condII,
    ex_dist_pond$sil_kmeans_Riemann$condIII
  )
  
  ## --------------- k-means : ACPF -------------------------------
  sil_kmeans_ACPF = rbind(
    ex_dist_pond$sil_kmeans_ACPF$condI,
    ex_dist_pond$sil_kmeans_ACPF$condII,
    ex_dist_pond$sil_kmeans_ACPF$condIII
  )
  
  ## --------------- PAM : Riemann --------------------------------
  sil_PAM_Riemann = matrix(NA, ncol = 6, nrow = 3)
  sil_PAM_Riemann[, 1] = ex_dist_pond$sil_PAM_Riemann$D0 %>% unlist()
  sil_PAM_Riemann[, 6] = ex_dist_pond$sil_PAM_Riemann$D1 %>% unlist()
  sil_PAM_Riemann[1, 2:5] = ex_dist_pond$sil_PAM_Riemann$Dp$condI   %>% unlist()
  sil_PAM_Riemann[2, 2:5] = ex_dist_pond$sil_PAM_Riemann$Dp$condII  %>% unlist()
  sil_PAM_Riemann[3, 2:5] = ex_dist_pond$sil_PAM_Riemann$Dp$condIII %>% unlist()
  
  ## --------------- PAM : ACPF -----------------------------------
  sil_PAM_ACPF = matrix(NA, ncol = 6, nrow = 3)
  sil_PAM_ACPF[, 1] = ex_dist_pond$sil_PAM_ACPF$D0 %>% unlist()
  sil_PAM_ACPF[, 6] = ex_dist_pond$sil_PAM_ACPF$D1 %>% unlist()
  sil_PAM_ACPF[1, 2:5] = ex_dist_pond$sil_PAM_ACPF$Dp$condI   %>% unlist()
  sil_PAM_ACPF[2, 2:5] = ex_dist_pond$sil_PAM_ACPF$Dp$condII  %>% unlist()
  sil_PAM_ACPF[3, 2:5] = ex_dist_pond$sil_PAM_ACPF$Dp$condIII %>% unlist()
  
  return(list(
    sil_Riemann        = sil_Riemann,
    sil_ACPF           = sil_ACPF,
    sil_kmeans_Riemann = sil_kmeans_Riemann,
    sil_kmeans_ACPF    = sil_kmeans_ACPF,
    sil_PAM_Riemann    = sil_PAM_Riemann,
    sil_PAM_ACPF       = sil_PAM_ACPF
  ))
}




CLUS_RES = function(t, sim){
  
    #sim = Simulation(t)
    ex_dist_pond = Execution_distance_pond(t, sim = sim)
    
    ARI = f_arranger_ARI(ex_dist_pond) 
    SIL = f_arranger_silhouette(ex_dist_pond)
    
    # ----------------  Méthodes concurentes -----------------------
    ari_yu_condI   = sapply(seq(0, 1, by = 0.2), function(i) Ari_distance_yu(sim[[1]]%>%as.data.frame(), i))
    ari_yu_condII  = sapply(seq(0, 1, by = 0.2), function(i) Ari_distance_yu(sim[[2]]%>%as.data.frame(), i))
    ari_yu_condIII = sapply(seq(0, 1, by = 0.2), function(i) Ari_distance_yu(sim[[3]]%>%as.data.frame(), i))
    
    # clusterMLD
    ari_clusterMLD_condI = mclust::adjustedRandIndex(sim[[1]]$class, clusteringMLD(sim[[1]][,-c(1:2)], 
                                                                                   choix_class = "CH"))
    ari_clusterMLD_condII = mclust::adjustedRandIndex(sim[[2]]$class, clusteringMLD(sim[[2]][,-c(1:2)], 
                                                                                    choix_class = "CH"))
    ari_clusterMLD_condIII = mclust::adjustedRandIndex(sim[[3]]$class, clusteringMLD(sim[[3]][,-c(1:2)], 
                                                                                     choix_class = "CH"))
    
    Ari_yu = rbind(ari_yu_condI, ari_yu_condII, ari_yu_condIII)
    Ari_clusterMLD = rbind(ari_clusterMLD_condI, ari_clusterMLD_condII, ari_clusterMLD_condIII)
    
    # Baseline
    Ari_baseline = ex_dist_pond$Ari_baseline%>%unlist()
    sil_baseline = ex_dist_pond$sil_baseline%>%unlist()
    
    return(list(ARI = ARI, 
                SIL = SIL,
                Ari_yu = Ari_yu,
                Ari_clusterMLD = Ari_clusterMLD,
                Ari_baseline = Ari_baseline, 
                sil_baseline = sil_baseline))
    
}




# -------------------------- Application sur 5 simulations ----------------------------------------
# Voir le dossier Data (100 premiers) sont les données pour les resultats soumis à EGC2026
Tableau_sim = readRDS("./Data/List_data_simu.rds")
# 25 premiers listes
res = lapply(1:25, function(i) CLUS_RES(t, Tableau_sim[[i]]))
# ------------------- Synthétiser les résultats ----------------------------------------------------------

Extraction_des_res = function(res){
  # res : liste de résultats de CLUS_RES
  
  # ARI clusterMLD
  ARI_clusterMLD = do.call(cbind, lapply(res, function(x) x$Ari_clusterMLD))
  ARI_clusterMLD = data.frame(mean = rowMeans(ARI_clusterMLD), sd = apply(ARI_clusterMLD, 1, sd))
  
  # Baseline
  SIL_baseline = do.call(rbind, lapply(res, function(x) x$sil_baseline))
  SIL_baseline = data.frame(mean = colMeans(SIL_baseline), sd = apply(SIL_baseline, 2, sd))
  
  ARI_baseline = do.call(rbind, lapply(res, function(x) x$Ari_baseline))
  ARI_baseline = data.frame(mean = colMeans(ARI_baseline), sd = apply(ARI_baseline, 2, sd))
  
  # ARI yu
  array_3D_yu = simplify2array(lapply(res, function(x) x$Ari_yu))
  mean_yu =  apply(array_3D_yu, c(1, 2), mean)
  sd_yu   =  apply(array_3D_yu, c(1, 2), sd)
  ARI_yu  = list(mean_yu = mean_yu, sd_yu = sd_yu)
  
  # Riemann (CAH)
  array_3D_Riemann = simplify2array(lapply(res, function(x) x$ARI$Ari_Riemann))
  mean_Riemann = apply(array_3D_Riemann, c(1, 2), mean)
  sd_Riemann   = apply(array_3D_Riemann, c(1, 2), sd)
  Ari_Riemann  = list(mean_Riemann = mean_Riemann, sd_Riemann = sd_Riemann)
  
  # ACPF (CAH)
  array_3D_ACPF = simplify2array(lapply(res, function(x) x$ARI$Ari_ACPF))
  mean_ACPF = apply(array_3D_ACPF, c(1, 2), mean)
  sd_ACPF   = apply(array_3D_ACPF, c(1, 2), sd)
  Ari_ACPF  = list(mean_ACPF = mean_ACPF, sd_ACPF = sd_ACPF)
  
  # k-means Riemann
  array_3D_kmeans_Riemann = simplify2array(lapply(res, function(x) x$ARI$Ari_kmeans_Riemann))
  mean_kmeans_Riemann = apply(array_3D_kmeans_Riemann, c(1, 2), mean)
  sd_kmeans_Riemann   = apply(array_3D_kmeans_Riemann, c(1, 2), sd)
  Ari_kmeans_Riemann  = list(mean_kmeans_Riemann = mean_kmeans_Riemann, sd_kmeans_Riemann = sd_kmeans_Riemann)
  
  # k-means ACPF
  array_3D_kmeans_ACPF = simplify2array(lapply(res, function(x) x$ARI$Ari_kmeans_ACPF))
  mean_kmeans_ACPF = apply(array_3D_kmeans_ACPF, c(1, 2), mean)
  sd_kmeans_ACPF   = apply(array_3D_kmeans_ACPF, c(1, 2), sd)
  Ari_kmeans_ACPF  = list(mean_kmeans_ACPF = mean_kmeans_ACPF, sd_kmeans_ACPF = sd_kmeans_ACPF)
  
  # PAM Riemann
  array_3D_PAM_Riemann = simplify2array(lapply(res, function(x) x$ARI$Ari_PAM_Riemann))
  mean_PAM_Riemann = apply(array_3D_PAM_Riemann, c(1, 2), mean)
  sd_PAM_Riemann   = apply(array_3D_PAM_Riemann, c(1, 2), sd)
  Ari_PAM_Riemann  = list(mean_PAM_Riemann = mean_PAM_Riemann, sd_PAM_Riemann = sd_PAM_Riemann)
  
  # PAM ACPF
  array_3D_PAM_ACPF = simplify2array(lapply(res, function(x) x$ARI$Ari_PAM_ACPF))
  mean_PAM_ACPF = apply(array_3D_PAM_ACPF, c(1, 2), mean)
  sd_PAM_ACPF   = apply(array_3D_PAM_ACPF, c(1, 2), sd)
  Ari_PAM_ACPF  = list(mean_PAM_ACPF = mean_PAM_ACPF, sd_PAM_ACPF = sd_PAM_ACPF)
  
  # Silhouette Riemann
  array_3D_Riemann_sil = simplify2array(lapply(res, function(x) x$SIL$sil_Riemann))
  Sil_Riemann = list(mean_Riemann = apply(array_3D_Riemann_sil, c(1, 2), mean),
                     sd_Riemann   = apply(array_3D_Riemann_sil, c(1, 2), sd))
  
  # Silhouette ACPF
  array_3D_ACPF_sil = simplify2array(lapply(res, function(x) x$SIL$sil_ACPF))
  Sil_ACPF = list(mean_ACPF = apply(array_3D_ACPF_sil, c(1, 2), mean),
                  sd_ACPF   = apply(array_3D_ACPF_sil, c(1, 2), sd))
  
  # Silhouette k-means Riemann
  array_3D_kmeans_Riemann_sil = simplify2array(lapply(res, function(x) x$SIL$sil_kmeans_Riemann))
  Sil_kmeans_Riemann = list(mean_kmeans_Riemann = apply(array_3D_kmeans_Riemann_sil, c(1, 2), mean),
                            sd_kmeans_Riemann   = apply(array_3D_kmeans_Riemann_sil, c(1, 2), sd))
  
  # Silhouette k-means ACPF
  array_3D_kmeans_ACPF_sil = simplify2array(lapply(res, function(x) x$SIL$sil_kmeans_ACPF))
  Sil_kmeans_ACPF = list(mean_kmeans_ACPF = apply(array_3D_kmeans_ACPF_sil, c(1, 2), mean),
                         sd_kmeans_ACPF   = apply(array_3D_kmeans_ACPF_sil, c(1, 2), sd))
  
  # Silhouette PAM Riemann
  array_3D_PAM_Riemann_sil = simplify2array(lapply(res, function(x) x$SIL$sil_PAM_Riemann))
  Sil_PAM_Riemann = list(mean_PAM_Riemann = apply(array_3D_PAM_Riemann_sil, c(1, 2), mean),
                         sd_PAM_Riemann   = apply(array_3D_PAM_Riemann_sil, c(1, 2), sd))
  
  # Silhouette PAM ACPF
  array_3D_PAM_ACPF_sil = simplify2array(lapply(res, function(x) x$SIL$sil_PAM_ACPF))
  Sil_PAM_ACPF = list(mean_PAM_ACPF = apply(array_3D_PAM_ACPF_sil, c(1, 2), mean),
                      sd_PAM_ACPF   = apply(array_3D_PAM_ACPF_sil, c(1, 2), sd))
  
  return(list(
    Ari_Riemann        = Ari_Riemann,
    Ari_ACPF           = Ari_ACPF,
    Ari_kmeans_Riemann = Ari_kmeans_Riemann,
    Ari_kmeans_ACPF    = Ari_kmeans_ACPF,
    Ari_PAM_Riemann    = Ari_PAM_Riemann,
    Ari_PAM_ACPF       = Ari_PAM_ACPF,
    ARI_clusterMLD     = ARI_clusterMLD,
    ARI_yu             = ARI_yu,
    ARI_baseline       = ARI_baseline,
    Sil_Riemann        = Sil_Riemann,
    Sil_ACPF           = Sil_ACPF,
    Sil_kmeans_Riemann = Sil_kmeans_Riemann,
    Sil_kmeans_ACPF    = Sil_kmeans_ACPF,
    Sil_PAM_Riemann    = Sil_PAM_Riemann,
    Sil_PAM_ACPF       = Sil_PAM_ACPF,
    SIL_baseline       = SIL_baseline
  ))
}

# resultats 
Extraction_des_res(res)


















