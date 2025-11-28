#------------------------------ Données rélles labellisées -------------------------------

library(fds)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(fdapace) # important Yu et al. 2025 utilise fdapace::FPCAder
library(cluster)
library(parallel)
library(doParallel)
library(mclust) # ARI
library(fda)
library(fda.usc) # Exemples de données labellisées
library(clusterMLD)

# Signaux audio de phonèmes (5 classes)
data("phoneme")

# Spectres de viande (classes par teneur en gras : 2 classes)
data("tecator")

# Données de croissance (garçons vs filles)
data("growth")

# ------------ Transformation des données en un tableau d'individus ----------------------------------
# Les temps de mesure sont représentés par des colonnes
# Colonnes : identifiant, classe, instants de mesure

# Phonème
data_phoneme = cbind(id = 1:length(phoneme$classlearn), class = phoneme$classlearn, 
                     y = phoneme$learn$data)%>%as.data.frame()
colnames(data_phoneme) = c("id", "class", as.character(phoneme$learn$argvals))

# tecator avec 2 groupes selon ---> fat >= 20% 
data_tecator = cbind(id = 1:dim(tecator$y)[1], class = ifelse(tecator$y$Fat >= 20, 1, 0), 
                     y = tecator$absorp.fdata$data)%>%as.data.frame()
colnames(data_tecator) = c("id", "class", as.character(tecator$absorp.fdata$argvals))

# growth
data_growth = cbind(growth$hgtm, growth$hgtf)%>%t()
class = data_growth%>%rownames() 
class = ifelse(grepl("^boy", class), 1, 2)
data_growth = cbind(id = 1:length(class), class, data_growth)%>%as.data.frame()



#------------------- Simulations  -----------------------------------------------------------------


f_list = list(
  # 1) f1
  function(t) {
    0.9 * cos(2 * pi * t)
  },
  
  # 2) f2
  function(t) {
    0.85 * sin(4 * pi * t + 0.4)
  },
  
  # 3) f3
  function(t) {
    0.95 * exp(-((t - 0.3)^2) / (2 * 0.07^2)) - 0.15
  },
  
  # 4) f4
  function(t) {
    -0.9 + 1.8 * t
  }
)

##--------------------- Scénario I ----------------------------------------------------


sim_condition_I = function(n_per_clus, t, f_list, nb_cluster = 4){
  
  # Génération des données
  sim_data = do.call(rbind, lapply(1:nb_cluster, function(k) {
    do.call(rbind, lapply(1:n_per_clus, function(i) {
      data.frame(
        id = i,
        class = k,
        t = t,
        y = f_list[[k]](t) + rnorm(length(t), 0, 0.5) # sigma^2 = 0.5^2
      )
    }))
  }))
  
  sim_data_pivot = sim_data%>%
    pivot_wider(id_cols = c(id, class), names_from = t, values_from = y)
  
  sim_data_pivot$id = 1:(n_per_clus*nb_cluster)
  
  return(sim_data_pivot)
  
}


## ------------------ Scénario II (effets aléatoires identique pour tous les clusters) ------------------

# On considère un effet aléatoire identique pour tous les clusters
# Les individus d'un même cluster présentent une variabilité
# et cette variabilité intra-cluster est supposée identique pour tous les clusters




sim_condition_II = function(n_per_clus, t, f_list, nb_cluster = 4){
  # Matrice de covariance pour les effets aléatoires
  Sigma = matrix(c(0.09, 0.09, -0.045,
                    0.09, 0.25, -0.025,
                    -0.045, -0.025, 0.25), nrow = 3)
  
  # Génération des données
  sim_data = do.call(rbind, lapply(1:nb_cluster, function(k) {
    do.call(rbind, lapply(1:n_per_clus, function(i) {
      
      # Tirage aléatoire des coefficients (b0, b1, b2) ~ N(0, Sigma)
      b_i = mvrnorm(1, mu = c(0, 0, 0), Sigma = Sigma)
      
      # Calcul de r_i(t) pour chaque point de temps
      r_i = b_i[1] + b_i[2] * t + b_i[3] * t^2
      
      # Valeurs y = f_k(t) + r_i(t) + epsilon
      data.frame(
        id = i,
        class = k,
        t = t,
        y = f_list[[k]](t) + r_i + rnorm(length(t), 0, 0.5)
      )
    }))
  }))
  
  # Pivot pour mettre les trajectoires sur une ligne par individu
  sim_data_pivot = sim_data %>%
    pivot_wider(id_cols = c(id, class), names_from = t, values_from = y)
  
  # Réindexation des individus
  sim_data_pivot$id = 1:(n_per_clus * nb_cluster)
  
  return(sim_data_pivot)
  
}



##----------------------------------- Scénario III (effets aléatoires  clusters-spécifique) -----------------------

sim_condition_III = function(n_per_clus, t, f_list, nb_cluster = 4){
  
  # Définition des matrices de covariance spécifiques à chaque cluster
  Sigma_list = list(
    matrix(c(0.160, 0.144, -0.072,
             0.144, 0.360, -0.036,
             -0.072, -0.036, 0.360), nrow = 3, byrow = TRUE),
    
    matrix(c(0.090, 0.030, -0.120,
             0.030, 0.250, 0.025,
             -0.120, 0.025, 0.250), nrow = 3, byrow = TRUE),
    
    matrix(c(0.040, -0.048, 0.024,
             -0.048, 0.160, 0.016,
             0.024, 0.016, 0.160), nrow = 3, byrow = TRUE),
    
    matrix(c(0.010, -0.004, 0.016,
             -0.004, 0.040, -0.004,
             0.016, -0.004, 0.040), nrow = 3, byrow = TRUE)
  )
  
  
  # Génération des données
  sim_data = do.call(rbind, lapply(1:nb_cluster, function(k) {
    do.call(rbind, lapply(1:n_per_clus, function(i) {
      
      # Tirage des coefficients aléatoires 
      b_ik = mvrnorm(1, mu = c(0, 0, 0), Sigma = Sigma_list[[k]])
      
      # Calcul du terme aléatoire r_ik(t)
      r_ik = b_ik[1] + b_ik[2] * t + b_ik[3] * t^2
      
      # Valeurs simulées
      data.frame(
        id = i,
        class = k,
        t = t,
        y = f_list[[k]](t) + r_ik + rnorm(length(t), 0, 0.5)
      )
    }))
  }))
  
  # Pivot pour avoir une ligne par individu
  sim_data_pivot = sim_data %>%
    pivot_wider(id_cols = c(id, class), names_from = t, values_from = y)
  
  # Réindexation des individus
  sim_data_pivot$id = 1:(n_per_clus * nb_cluster)
  
  return(sim_data_pivot)
}


#---------------------  Lissage MCP -----------------------------------------------------

liss = function(df, l_grille = 10^seq(-5, 5, length.out = 1000), D = 50){
  # Lissage par B-splines avec des nœuds placés sur les quantiles des observations
  # df : matrice ou data.frame où chaque ligne représente une courbe
  # D  : nombre de fonctions de base
 
  t = as.numeric(colnames(df)) 
  q = quantile(t)
  d = (D-2)/4 # nbre de noeuds approximative à placer dans chaque intervalle de q
  noeuds = c(seq(q[1], q[2], length.out = d)[-d], 
             seq(q[2], q[3], length.out = d+1)[-(d+1)], 
             seq(q[3], q[4], length.out = d+2)[-(d+2)], 
             seq(q[4], q[5], length.out = d))
  
  base = create.bspline.basis(rangeval = range(t), nbasis = D, breaks = noeuds)
  
  # optimiser lambda
  gcv = function(lambda){ smooth.basis(t, df%>%t(), fdPar(base, Lfdobj = 2, lambda = lambda))$gcv}
  # lambda optimaux
  lambda_optimaux = l_grille[apply(sapply(l_grille, gcv), MARGIN = 1, FUN = function(x) which.min(x))]
  fd_list = list()
  for (i in 1:dim(df)[1]){fd_list[[i]] = smooth.basis(t, df[i, ]%>%unname()%>%unlist(), 
                                                      fdPar(base, 2, lambda = lambda_optimaux[i]))$fd}
  
  return(fd_list)
}



# ---------------- Approximation de Riemman (D0, D1, Dp) -----------------------------


##------------------ D0  -----------------------------------
calculate_D0_matrix_parallel = function(fd_list, fine_grid) {
  
  n = length(fd_list)
  D0_matrix = matrix(0, nrow = n, ncol = n)
  
  # Créer un cluster avec le nombre de cœurs disponibles
  no_cores = detectCores() - 1  # Laisser un cœur libre
  cl = makeCluster(no_cores)
  
  # Forcer assignation dans .GlobalEnv (Erreur dès fois sinon selon les pcs)
  assign("fd_list", fd_list, envir = .GlobalEnv)
  assign("fine_grid", fine_grid, envir = .GlobalEnv)
  
  # Exporter les variables nécessaires aux workers
  clusterExport(cl, varlist = c("fd_list", "fine_grid", "eval.fd"))
  
  # Fonction pour calculer une ligne de la matrice
  compute_row = function(i) {
    n = length(fd_list)
    row_values = numeric(n)
    
    # Recalculer delta à l'intérieur de la fonction (évite l'erreur)
    delta = diff(range(fine_grid)) / (length(fine_grid) - 1)
    
    for (j in (i+1):n) {
      fd_i_vals = eval.fd(fine_grid, fd_list[[i]])
      fd_j_vals = eval.fd(fine_grid, fd_list[[j]])
      
      diff_squared = (fd_i_vals - fd_j_vals)^2
      integral = sum(diff_squared) * delta
      
      row_values[j] = sqrt(integral)
    }
    return(row_values)
  }
  
  # Appliquer la fonction en parallèle
  result_list = parLapply(cl, 1:(n-1), compute_row)
  
  # Fermer le cluster
  stopCluster(cl)
  
  # Remplir la matrice de résultats
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      D0_matrix[i, j] = result_list[[i]][j]
      D0_matrix[j, i] = D0_matrix[i, j]  
    }
  }
  
  return(D0_matrix)
}



## -------------------------------- D1  -------------------------------------
calculate_D1_matrix_parallel = function(fd_list, fine_grid) {
  n = length(fd_list)
  D1_matrix = matrix(0, nrow = n, ncol = n)
  
  # Créer un cluster
  no_cores = detectCores() - 1
  cl = makeCluster(no_cores)
  
  # Forcer assignation dans .GlobalEnv (Erreur dès fois sinon selon les pcs)
  assign("fd_list", fd_list, envir = .GlobalEnv)
  assign("fine_grid", fine_grid, envir = .GlobalEnv)
  
  # Exporter les variables nécessaires
  clusterExport(cl, varlist = c("fd_list", "fine_grid", "eval.fd"))
  
  # Fonction pour calculer une ligne de la matrice
  compute_row = function(i) {
    fd_list_deriv_i = eval.fd(fine_grid, fd_list[[i]], Lfdobj = 1)
    row_values = numeric(n)
    
    # Recalculer delta à l'intérieur de la fonction (évite l'erreur)
    delta = diff(range(fine_grid)) / (length(fine_grid) - 1)
    
    for (j in 1:n) {
      if (i != j) {
        fd_list_deriv_j = eval.fd(fine_grid, fd_list[[j]], Lfdobj = 1)
        diff = fd_list_deriv_i - fd_list_deriv_j
        row_values[j] = sum(diff^2)* delta  # Approximation par somme de Riemann
      }
    }
    
    return(row_values)
  }
  
  # Appliquer en parallèle
  result_list = parLapply(cl, 1:n, compute_row)
  
  # Fermer le cluster
  stopCluster(cl)
  
  # Remplir la matrice
  for (i in 1:n) {
    D1_matrix[i, ] = result_list[[i]]
  }
  
  return(D1_matrix)
}



## -------------------------------- Dp -------------------------------------

# partir des de d1 et d0

calculate_Dp_matrix_byD0D1 = function(D0_matrix, D1_matrix, omega, STANDARDIZE = TRUE){
  if (!STANDARDIZE){
    return(sqrt((1-omega)*D0_matrix^2 + omega*D1_matrix^2))  
  } else {
    # normalisation important pour la pondération car echelle pouvant être différentes
    D0_min = min(D0_matrix)
    D0_max = max(D0_matrix)
    D0_standardized = (D0_matrix - D0_min) / (D0_max - D0_min)
    
    D1_min = min(D1_matrix)
    D1_max = max(D1_matrix)
    D1_standardized = (D1_matrix - D1_min) / (D1_max - D1_min)
    
    return(sqrt((1 - omega) * D0_standardized^2 + omega * D1_standardized^2))
  }
  
}




# -------------------- Distances pondérées (via ACPF) -------------------------------


## ------------------ sans re-affectation ------------------------------------------

pca_sans_reaffect = function(data_liss, t,  pve = 0.99){
  # data_liss : données lissées sous forme de liste de fd (données fonctionnelles)
  # t         : vecteur des temps de mesure
  # pve       : pourcentage de variabilité conservée par l'ACPF
  
  
  fdobj = fd(sapply( 1:length(data_liss), function (i) data_liss[[i]]$coefs), data_liss[[1]]$basis)
  fdobj_deriv = deriv.fd(fdobj, deriv = 1)
  
  pca_res = pca.fd(fdobj, nharm = 50, centerfns = TRUE)
  pca_res_deriv = pca.fd(fdobj_deriv, nharm = 50, centerfns = TRUE)
  
  # nombre de composantes gardées
  nb1 = which(cumsum(pca_res$varprop)>= pve)[1] 
  nb2 = which(cumsum(pca_res_deriv$varprop)>= pve)[1]
  
  return(list(coords = pca_res$scores[, 1:nb1], 
              coords_deriv = pca_res_deriv$scores[, 1:nb2]))
  
}


pca_Dp = function(pca_res, omega = 0.5){
  # Donne la matrice de dissimilarité 
  # pca_res : retour de la fonction pca_sans_reaffect
  # omega   : valeur de la pondération 
  coords = pca_res$coords
  coords_deriv = pca_res$coords_deriv
  
  # normalisation min-max 
  min_D0 = min(coords)
  max_D0 = max(coords)
  min_D1 = min(coords_deriv)
  max_D1 = max(coords_deriv)
  
  coords = (coords - min_D0) / (max_D0 - min_D0)
  coords_deriv = (coords_deriv - min_D1) / (max_D1 - min_D1)
  
  
  #Obtenir tous les couples (i, j) 
  indices = combn(nrow(pca_res$coords), 2, simplify = FALSE)
  
  nb_cores = detectCores() - 1
  cl = makeCluster(nb_cores)
  registerDoParallel(cl)
  
  results = foreach(pair = indices, 
                    .combine = rbind, 
                    .export = c("coords", "coords_deriv", "omega")) %dopar% {
                      i = pair[1]
                      j = pair[2]
                      
                      d = sqrt((1 - omega) * sqrt(sum((coords[i, ] - coords[j, ])^2)) +
                                 (omega * sqrt(sum((coords_deriv[i, ] - coords_deriv[j, ])^2))))
                      
                      c(i, j, d)  # retourne une ligne : i, j, distance
                    }
  
  stopCluster(cl)# arrêt
  tab = matrix(0, nrow = nrow(pca_res$coords), ncol = nrow(pca_res$coords))
  for (k in 1:nrow(results)) {
    i = results[k, 1]
    j = results[k, 2]
    d = results[k, 3]
    tab[i, j] = d
    tab[j, i] = d
  }
  
  
  
  return(as.dist(tab))
  
}



## ----------------  Re-affectation subspace --------------------------------------

# Voir la documentation sur les différents paramètres dans le dossier Yu-main

source("./Yu-main/SPFCder.R") 
source("./Yu-main/SPFConlyder.R")


Ari_distance_yu = function(data, wt = 0.5){
  # Calcule l'ARI de la méthode Yu et al. (2025)
  # data : tableau avec id et class dans les deux premières colonnes,
  #        suivies des colonnes de temps t
  # wt   : la pondération (=:omega)
  
  t = as.numeric(colnames(data)[-(1:2)])
  # calcul des derivées avec fda.usc::fdata
  fdata_obj = fdata(data[, -c(1:2)], argvals = t) 
  yder = fdata.deriv(fdata_obj, nderiv = 1)$data
  y = split(data[, -c(1:2)], seq(nrow(data)))
  y = lapply(y, function(x) as.numeric(x)) 
  yder = split(yder, seq(nrow(yder)))
  t  = lapply(1:dim(data)[1], function(i) t)
  # appliquons maintenant l'approche de Yu
  class = kCFCder(y, yder, t = t, wt = wt, k = data$class%>%unique()%>%length())$cluster
  # ARI caclul
  return(mclust::adjustedRandIndex(class, data$class))
}




# ---------------------- ClusterMLD (Zhou 2023)---------------------------------------------

clusteringMLD = function(data, choix_class = "CH", No.Class){
  # Renvoie la partition obtenue par la CAH fonctionnelle
  # choix_class : méthode de coupure du dendrogramme
  # No.Class    : nombre de clusters fixé a priori
  
  data_t = as.matrix(data%>%t())
  data_long = cbind(t = as.numeric(colnames(data)), data_t)%>%as.data.frame()
  colnames(data_long) = c("t", 1:dim(data)[1])
  data_long = pivot_longer(data_long, cols = colnames(data_long)[-1])
  
  id = as.numeric(data_long$name)
  x = data_long$t
  y = data_long$value
  
  cl = makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  clust_obj = LongDataCluster(x = x, Y = y , id = id, parallel = TRUE, No.Class = No.Class)
  stopCluster(cl)
  
  if(choix_class == "CH"){
    class_id  =  clust_obj$Dat.label %>%
      distinct(id, label.CH) %>%  
      arrange(id) 
    
    if(dim(class_id)[1]!= dim(data)[1]){
      missing_add = data.frame(id = setdiff(1:dim(data)[1], class_id$id), label.CH = 0)
      class_id_full = rbind(class_id[class_id$id %in% (1:dim(data)[1]), ],  missing_add)%>%arrange(id)
      class = class_id_full$label.CH
    } else{
      class = class_id$label.CH
    }                    
    
  } else {
    class_id  =  clust_obj$Dat.label %>%
      distinct(id, label.Gapb) %>%  
      arrange(id) 
    
    if(dim(class_id)[1]!= dim(data)[1]){
      missing_add = data.frame(id = setdiff(1:dim(data)[1], class_id$id), label.Gapb = 0)
      class_id_full = rbind(class_id[class_id$id %in% (1:dim(data)[1]), ],  missing_add)%>%arrange(id)
      class = class_id_full$label.Gapb
    } else{
      class = class_id$label.Gapb
    }   
  }
  
  return(class)
}



# ------------------------------ k-means integrant ces distances pondérées ---------------------------------------------------------

k_means_dw = function(X, X_deriv, K, method = "ACPF", omega = 0.5, max_iter = 100, nb_init = 20, fine_grid = NULL){

  # k-means intégrant une distance pondérée 
  # X        : scores de l'ACPF sur les fonctions lissées
  #            ou évaluations des fonctions lissées sur une grille fine
  # X_deriv  : scores de l'ACPF sur les dérivées
  #            ou évaluations des dérivées sur une grille fine
  # K        : nombre de classes
  # method   : approches pour la distance pondérée, "ACPF" ou "riemann"
  # omega    : valeur de la pondération
  # max_iter : nombre maximal d'itérations
  # nb_init  : nombre d'initialisation des centroïdes
  # fine_grid: grille fine d'évaluation
  
  
  # Normalisation min-max de chaque matrice
  X_min = min(X)
  X_max = max(X)
  X = (X - X_min) / (X_max - X_min)
  
  X_deriv_min = min(X_deriv)
  X_deriv_max = max(X_deriv)
  X_deriv = (X_deriv - X_deriv_min) / (X_deriv_max - X_deriv_min)
  
  n = nrow(X)
  p = ncol(X)
  q = ncol(X_deriv)
  
  
  delta = 1
  if(method == "riemann"){
    if(is.null(fine_grid)){
      stop("fine_grid can not be NULL")
    } else{
      delta = diff(range(fine_grid))/(p - 1) 
    }
  }
  
  
  
  # Faire une boucle pour stocker les clusters et les inerties
  list_inertie  = numeric(nb_init)
  list_clusters = vector("list", nb_init)
  
  for (init in 1:nb_init) {
    
    # Initialisation aléatoire des centres
    idx_init = sample(1:n, K)
    centres_X = X[idx_init, , drop = FALSE]
    centres_X_deriv = X_deriv[idx_init, , drop = FALSE]
    
    clusters = rep(0, n)
    
    for (iter in 1:max_iter) {
      
      # Matrice des distances pondérées
      dmat_sq = matrix(0, n, K)  
      for (k in 1:K) {
        dx       = sweep(X, 2, centres_X[k, ], "-")
        dx_deriv = sweep(X_deriv, 2, centres_X_deriv[k, ], "-")
        
        # distance^2 pondérée
        dmat_sq[, k] = (1 - omega) * rowSums(dx^2) * delta +
          omega      * rowSums(dx_deriv^2) * delta
      }
      
      # Pour l'affectation
      dmat = sqrt(dmat_sq)
      clusters_new = max.col(-dmat)
      
      # Mise à jour des centres
      new_centres_X = matrix(NA, K, p)
      new_centres_X_deriv = matrix(NA, K, q)
      
      for (k in 1:K) {
        members = which(clusters_new == k)
        if (length(members) == 0) {
          # réinitialisation si cluster vide
          rand_idx = sample(1:n, 1)
          new_centres_X[k, ]       = X[rand_idx, ]
          new_centres_X_deriv[k, ] = X_deriv[rand_idx, ]
        } else {
          new_centres_X[k, ]       = colMeans(X[members, , drop = FALSE])
          new_centres_X_deriv[k, ] = colMeans(X_deriv[members, , drop = FALSE])
        }
      }
      
      centres_X       = new_centres_X
      centres_X_deriv = new_centres_X_deriv
      
      if (all(clusters == clusters_new)) break
      clusters = clusters_new
    }
    
    # inertie de l'initialisation
    inertie_init = sum(dmat_sq[cbind(1:n, clusters)])
    list_inertie[init]      = inertie_init
    list_clusters[[init]]   = clusters
    
  }
  
  # choix inertie minimale
  best_idx = which.min(list_inertie) # idx de la meilleure partition
  
  
  return(list(cluster = list_clusters[[best_idx]], inertie = list_inertie[[best_idx]]))
  
}


# --------------------------------- PAM clustering -------------------------------------------------------------------------------------
pam_clus = function(list_dist, class_true){
  # fournit l'ARI et silhouette de pam
  if(is.list(list_dist)){
    classes = lapply(1:length(list_dist), function (i) cluster::pam(list_dist[[i]], k = unique(class_true)%>%length(), diss = T)$clustering)
    # ARI 
    ARI = sapply(1:length(classes), function(i) mclust::adjustedRandIndex(classes[[i]], class_true))
    # SIL
    SIL = sapply(1:length(classes), function(i) cluster::silhouette(classes[[i]], dist = list_dist[[i]])[, 3]%>%mean())
    
  } else{
    classe = cluster::pam(list_dist, k = unique(class_true)%>%length(), diss = T)$clustering
    # ARI
    ARI = mclust::adjustedRandIndex(classe, class_true)
    #SIL 
    SIL = cluster::silhouette(classe, dist = list_dist)[, 3]%>%mean()
  }
  
  return(list(ARI = ARI, SIL = SIL))
}


#--------------- Autres (visuels) ----------------------------------------------

##--------------- figures courbes brutes vs lissées -------------------------

gg_discret_continue = function(data, data_liss, fine_grid){
  
  data = cbind(type = rep("brutes", dim(data)[1]), data)
  mat_liss = sapply(1:dim(data)[1], function(i) eval.fd(fine_grid, data_liss[[i]]))%>%t()%>%data.frame()
  data_liss = cbind(type = rep("lissées", dim(data)[1]), id = 1:dim(data)[1], class = data$class, 
                    y = mat_liss) 
  colnames(data_liss) = c("type", "id", "class", as.character(fine_grid))
  
  data_pivot = pivot_longer(data, cols = -(1:3))
  data_liss_pivot = pivot_longer(data_liss, cols = -(1:3))
  
  gg = rbind(data_pivot, data_liss_pivot) %>%
    ggplot(aes(x = as.numeric(name), y = value, group = id, color = factor(class))) +
    
    # Lignes plus fines et alpha pour réduire la surcharge visuelle
    geom_line(linewidth = 0.8, alpha = 0.7) +
    
    # Facettes avec labels plus lisibles et un peu d'espace entre elles
    facet_wrap(~type, scales = "free_y", ncol = 2) +
    
    # Amélioration des titres et légendes
    labs(
      title = "Prétraitement par lissage MCP",
      subtitle = "Passage des données brutes à des données lissées",
      x = "Time",
      y = "Values",
      color = "Classe"
    ) +
    
    # Thème minimal et moderne
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 0, hjust = 1),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90")
    ) +
    
    # Palette de couleurs plus harmonieuse
    scale_color_brewer(palette = "Set1")
  
  
  return(gg)
  
  
}


##------------------ figures courbes reconstruites (fonctions) vs dérivées --------------------------------------------------

gg_f_f_prime = function(data, data_liss, fine_grid){
  
  mat_liss = sapply(1:dim(data)[1], function(i) eval.fd(fine_grid, data_liss[[i]]))%>%t()%>%data.frame()
  mat_liss_prime = sapply(1:dim(data)[1], function(i) eval.fd(fine_grid, data_liss[[i]], Lfdobj = 1))%>%t()%>%data.frame()
  data_liss = cbind(type = rep("fonctions", dim(data)[1]), id = 1:dim(data)[1], class = data$class, 
                    y = mat_liss) 
  colnames(data_liss) = c("type", "id", "class", as.character(fine_grid))
  data_liss_prime = cbind(type = rep("fonctions dérivées", dim(data)[1]), id = 1:dim(data)[1], class = data$class, 
                          y = mat_liss_prime) 
  colnames(data_liss_prime) = c("type", "id", "class", as.character(fine_grid))
  
  data_liss_pivot = pivot_longer(data_liss, cols = -(1:3))
  data_liss_prime_pivot = pivot_longer(data_liss_prime, cols = -(1:3))
  
  gg = rbind(data_liss_pivot, data_liss_prime_pivot) %>%
    ggplot(aes(x = as.numeric(name), y = value, group = id , color = factor(class)))+ 
    
    # Lignes plus fines et alpha pour réduire la surcharge visuelle
    geom_line(linewidth = 0.8, alpha = 0.7) +
    
    # Facettes avec labels plus lisibles et un peu d'espace entre elles
    facet_wrap(~factor(type, labels = c("(a)", "(b)")), scales = "free_y", ncol = 2) +
    
    # Amélioration des titres et légendes
    labs(
      #title = "Prétraitement par lissage MCP",
      #subtitle = "fonctions et fonctions dérivées évaluées sur une grille fine",
      x = "t",
      y = "Valeurs des courbes au temps t",
      color = "Classes"
    ) +
    
    # Thème minimal et moderne
    theme_bw(base_size = 14) +
    theme(
      #plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      #plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text      = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      legend.title   = element_text(size = 14, face = "bold"),      # titre légende
      legend.text    = element_text(size = 12, face = "bold"),      # texte légende
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90"),
      strip.background = element_rect(fill = "gray80", colour = "black"),  # fond facet
      strip.text       = element_text(size = 14, face = "bold")           # texte facet
    )+
    guides(color = guide_legend(override.aes = list(size = 3)))+
  
    # Palette de couleurs
    scale_color_brewer(palette = "Set1")
  
  return(gg)
  
}












































