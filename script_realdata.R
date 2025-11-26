source("./main.R")


# --------------------------- Lissage MCP -----------------------------------------
t_phoneme = as.numeric(colnames(data_phoneme[, -c(1:2)]))
t_tecator = as.numeric(colnames(data_tecator[, -c(1:2)]))
t_growth = as.numeric(colnames(data_growth[, -c(1:2)]))


phoneme_liss = liss(data_phoneme[, -c(1,2)])

tecator_liss = liss(data_tecator[, -c(1,2)])

growth_liss = liss(data_growth[, -c(1,2)])


##--------------- plot des données liss vs données normaux --------------------------------

fine_grid_phoneme = seq(min(t_phoneme), max(t_phoneme), length.out = 1000)
fine_grid_tecator = seq(min(t_tecator), max(t_tecator), length.out = 1000)
fine_grid_growth  = seq(min(t_growth), max(t_growth), length.out = 1000) 


gg_discret_continue(data_phoneme, phoneme_liss, fine_grid_phoneme)
gg_discret_continue(data_tecator, tecator_liss, fine_grid_tecator)
gg_discret_continue(data_growth, growth_liss, fine_grid_growth)



##------------------------ Plots des fonctions et de leurs derivées (Figure 5)-----------------------------------------

gg_f_f_prime(data_phoneme, phoneme_liss, fine_grid_phoneme)
gg_f_f_prime(data_tecator, tecator_liss, fine_grid_tecator)
gg_f_f_prime(data_growth, growth_liss, fine_grid_growth)


#-------------------- distances classiques (approximation de Riemman) ----------------------------------------

##-------------------------phoneme data ----------------------------------------------------------------------
dist_D0_phoneme = calculate_D0_matrix_parallel(phoneme_liss, fine_grid_phoneme)
dist_D1_phoneme = calculate_D1_matrix_parallel(phoneme_liss, fine_grid_phoneme)
dist_Dp_phoneme = lapply(seq(0, 1, by = 0.2), function(i) calculate_Dp_matrix_byD0D1(dist_D0_phoneme,
                                                         dist_D1_phoneme, i))

## ----------------------- tecator data ----------------------------------------------
dist_D0_tecator = calculate_D0_matrix_parallel(tecator_liss, fine_grid_tecator)
dist_D1_tecator = calculate_D1_matrix_parallel(tecator_liss, fine_grid_tecator)
dist_Dp_tecator = lapply(seq(0, 1, by = 0.2), function(i) calculate_Dp_matrix_byD0D1(dist_D0_tecator,
                                                                                     dist_D1_tecator, i))
## ------------------------ growth data -----------------------------------------------

dist_D0_growth = calculate_D0_matrix_parallel(growth_liss, fine_grid_growth)
dist_D1_growth = calculate_D1_matrix_parallel(growth_liss, fine_grid_growth)
dist_Dp_growth = lapply(seq(0, 1, by = 0.2), function(i) calculate_Dp_matrix_byD0D1(dist_D0_growth,
                                                                                    dist_D1_growth, i))

###------------------ clustering ----------------------------------------------------


clus = function(dissimilarite, nb_class = NULL, dist = TRUE){
  # Clus pratique une CAH avec linkage ward 
  # dissimilarité : matrice de distance, matrice des individus ou une liste de matrices de distances ou individus 
  # dist : indique si ce qu'on prends en entrer c'est une distance ou une liste de distance
  
  if(is.null(nb_class)){
    stop("Spécifier le nombre de classe")
  }
  
  if(dist){
    if(is.matrix(dissimilarite)){
      classe = cutree(hclust(dissimilarite%>%as.dist(), method = "ward.D2"), k = nb_class)
    } else{
      classe = list()
      for (i in 1:length(dissimilarite)){
        classe[[i]] = cutree(hclust(dissimilarite[[i]]%>%as.dist(), method = "ward.D2"), k = nb_class)
      }
    }
    
  } else{ # dist = FALSE 
    if(!is.list(dissimilarite)){
      classe = cutree(hclust(dist(dissimilarite), method = "ward.D2"), k = nb_class)
    }else{
      classe = list()
      for (i in 1:length(dissimilarite)){
          classe[[i]] = cutree(hclust(dissimilarite[[i]], method = "ward.D2"), k = nb_class)
        }
      } 
          
  }
   
  return(classe)
}
 
 


Ari_partition = function(classe, classe_true){
  # Calucl l'ARI de chaque partition
  if(!is.list(classe)){
    ari = round(mclust::adjustedRandIndex(classe, classe_true), 2)
  } else{
    ari = lapply(1:length(classe), function (i) round(mclust::adjustedRandIndex(classe[[i]], classe_true), 2))
  }
  
  return(ari)
  
}



# phoneme
Ari_partition(clus(dist_D0_phoneme, nb_class = 5), data_phoneme$class)
Ari_partition(clus(dist_D1_phoneme, nb_class = 5), data_phoneme$class)
Ari_partition(clus(dist_Dp_phoneme, nb_class = 5), data_phoneme$class)
# tecator
Ari_partition(clus(dist_D0_tecator, nb_class = 2), data_tecator$class)
Ari_partition(clus(dist_D1_tecator, nb_class = 2), data_tecator$class)
Ari_partition(clus(dist_Dp_tecator, nb_class = 2), data_tecator$class)
#growth
Ari_partition(clus(dist_D0_growth, nb_class = 2), data_growth$class)
Ari_partition(clus(dist_D1_growth, nb_class = 2), data_growth$class)
Ari_partition(clus(dist_Dp_growth, nb_class = 2), data_growth$class)




#-------------------- distances pondérées (ACPF) -----------------------------------------

## ------------------- sans réaffectation -------------------------------------------------------------
# Les ACPs fonctionnelles pour avoir les deux tableaux de scores
pca_res_phoneme = pca_sans_reaffect(phoneme_liss, fine_grid_phoneme)
pca_res_tecator = pca_sans_reaffect(tecator_liss, fine_grid_tecator)
pca_res_growth  = pca_sans_reaffect(growth_liss, fine_grid_growth)


# Les listes des distances D_p (sans reaffectation)
dist_Dp_sans_reaffect_phoneme = lapply(seq(0.2, 0.8, by = 0.2), function (i) pca_Dp(pca_res_phoneme, i))
dist_Dp_sans_reaffect_tecator = lapply(seq(0.2, 0.8, by = 0.2), function (i) pca_Dp(pca_res_tecator, i))
dist_Dp_sans_reaffect_growth =  lapply(seq(0.2, 0.8, by = 0.2), function (i) pca_Dp(pca_res_growth, i))


### --------------------------- phoneme -------------------------------------------------

Ari_partition(clus(pca_res_phoneme$coords, nb_class = 5, dist = FALSE),
              data_phoneme$class) # DO

Ari_partition(clus(pca_res_phoneme$coords_deriv, nb_class = 5, dist = FALSE),
              data_phoneme$class) # D1

Ari_partition(clus(dist_Dp_sans_reaffect_phoneme, nb_class = 5, dist = FALSE),
              data_phoneme$class) # Dp


### -------------------------- tecator -----------------------------------

Ari_partition(clus(pca_res_tecator$coords, nb_class = 2, dist = FALSE),
              data_tecator$class) # DO

Ari_partition(clus(pca_res_tecator$coords_deriv, nb_class = 2, dist = FALSE),
              data_tecator$class) # D1

Ari_partition(clus(dist_Dp_sans_reaffect_tecator, nb_class = 2, dist = FALSE),
              data_tecator$class) # Dp


### -------------------------- growth -------------------------------------

Ari_partition(clus(pca_res_growth$coords, nb_class = 2, dist = FALSE),
              data_growth$class) # DO

Ari_partition(clus(pca_res_growth$coords_deriv, nb_class = 2, dist = FALSE),
              data_growth$class) # D1

Ari_partition(clus(dist_Dp_sans_reaffect_growth, nb_class = 2, dist = FALSE),
              data_growth$class) # Dp


## ------------------- avec re-affectation (Yu et al. 2025) --------------------------------------------------------

ari_yu_phoneme = sapply(seq(0, 1, by = 0.2), function(i) Ari_distance_yu(data_phoneme, i))
ari_yu_tecator = sapply(seq(0, 1, by = 0.2), function(i) Ari_distance_yu(data_tecator, i))
ari_yu_growth  = sapply(seq(0, 1, by = 0.2), function(i) Ari_distance_yu(data_growth, i))

#----------------------- clusterMLD (Zhou 2023) -----------------------------------------------------------

ari_clusterMLD_phoneme = mclust::adjustedRandIndex(data_phoneme$class, clusteringMLD(data_phoneme[, -c(1:2)], 
                                                                               choix_class = "CH", No.Class = 5))
ari_clusterMLD_tecator = mclust::adjustedRandIndex(data_tecator$class, clusteringMLD(data_tecator[, -c(1:2)], 
                                                                                          choix_class = "CH", No.Class = 2))

ari_clusterMLD_growth = mclust::adjustedRandIndex(data_growth$class, clusteringMLD(data_growth[, -c(1:2)], 
                                                                                          choix_class = "CH", No.Class = 2))

# ---------------------------- baseline (clustering sur les coefficients estimés du lissage )---------------------------------------------------------------------
ari_baseline_phoneme = mclust::adjustedRandIndex(data_phoneme$class, cutree(hclust(dist(sapply(phoneme_liss, function(fdobj) fdobj$coefs)%>%t()),
                                                               method = "ward.D2"), k = 5)) 
                                                 
ari_baseline_tecator = mclust::adjustedRandIndex(data_tecator$class, cutree(hclust(dist(sapply(tecator_liss, function(fdobj) fdobj$coefs)%>%t()),
                                                                                   method = "ward.D2"), k = 2)) 

ari_baseline_growth = mclust::adjustedRandIndex(data_growth$class, cutree(hclust(dist(sapply(growth_liss, function(fdobj) fdobj$coefs)%>%t()),
                                                                                   method = "ward.D2"), k = 2)) 
  
                                                                                                                            


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


## ------------------------------- Riemann ------------------------------------------------------------------------------

# phoneme 
pam_clus(dist_Dp_phoneme, class_true = data_phoneme$class)
# tecator
pam_clus(dist_Dp_tecator, class_true = data_tecator$class)
# growth
pam_clus(dist_Dp_growth, class_true = data_growth$class)


## ------------------------------ ACPF --------------------------------------------------------------------------------

###----------------------------- phoneme ----------------------------------------------
# DO
pam_clus(dist(pca_res_phoneme$coords), class_true = data_phoneme$class)
#D1 
pam_clus(dist(pca_res_phoneme$coords_deriv), class_true = data_phoneme$class)
# Dp
pam_clus(dist_Dp_sans_reaffect_phoneme, class_true = data_phoneme$class)

# ---------------------------- tecator----------------------------------------------------
# DO
pam_clus(dist(pca_res_tecator$coords), class_true = data_tecator$class)
#D1 
pam_clus(dist(pca_res_tecator$coords_deriv), class_true = data_tecator$class)
# Dp
pam_clus(dist_Dp_sans_reaffect_tecator, class_true = data_tecator$class)
#------------------------------growth ----------------------------------------------------
# DO
pam_clus(dist(pca_res_growth$coords), class_true = data_growth$class)
#D1 
pam_clus(dist(pca_res_growth$coords_deriv), class_true = data_growth$class)
# Dp
pam_clus(dist_Dp_sans_reaffect_growth, class_true = data_growth$class)


# ----------------------------------------- k-means fonctionnelle --------------------------------------------


## ------------------------------------ Riemann ----------------------------------------------------

kmeans_clus_Riemann = function(list_dist, X_liss, fine_grid, class_true){
  
  X = t(sapply(1:length(X_liss), function(i) eval.fd(fine_grid, X_liss[[i]])))
  X_deriv = t(sapply(1:length(X_liss), function(i) eval.fd(fine_grid, X_liss[[i]], Lfdobj = 1)))
  
  classes = lapply(seq(0, 1, length.out = 6), function(i) k_means_dw(X, X_deriv, 
                                                                     K = class_true%>%unique()%>%length(), method = "riemann", 
                                                                     omega = i, fine_grid = fine_grid)$cluster)
  
  # ARI 
  ARI = sapply(1:length(classes), function(i) mclust::adjustedRandIndex(classes[[i]], class_true))
  # SIL
  SIL = sapply(1:length(classes), function(i) cluster::silhouette(classes[[i]], dist = list_dist[[i]])[, 3]%>%mean()) 
  
  return(list(ARI = ARI, SIL = SIL))
}


# phoneme
kmeans_clus_Riemann(dist_Dp_phoneme, phoneme_liss, fine_grid_phoneme, class_true = data_phoneme$class)
# tecator
kmeans_clus_Riemann(dist_Dp_tecator, tecator_liss, fine_grid_tecator, class_true = data_tecator$class)
# growth
kmeans_clus_Riemann(dist_Dp_growth, growth_liss, fine_grid_growth, class_true = data_growth$class)


## ------------------------------------ ACPF -------------------------------------------------------

kmeans_clus_ACPF = function(list_dist, pca_res_, class_true){
  classes = lapply(seq(0, 1, length.out = 6), function(i) k_means_dw(X = pca_res_$coords, X_deriv = pca_res_$coords_deriv, 
                                                                     K = class_true%>%unique()%>%length(), method = "ACPF", 
                                                                     omega = i)$cluster)
  
  # ARI 
  ARI = sapply(1:length(classes), function(i) mclust::adjustedRandIndex(classes[[i]], class_true))
  # SIL
  SIL = sapply(1:(length(classes)-2), function(i) cluster::silhouette(classes[[i]], dist = list_dist[[i]])[, 3]%>%mean()) 
  # sil omega = 0
  sil_0 = cluster::silhouette(classes[[1]], dist = dist(pca_res_$coords))[, 3]%>%mean()
  # sil omega = 1
  sil_1 = cluster::silhouette(classes[[6]], dist = dist(pca_res_$coords_deriv))[, 3]%>%mean()
  
  SIL = c(sil_0, SIL, sil_1)
  
  return(list(ARI = ARI, SIL = SIL))
}



# phoneme
kmeans_clus_ACPF(dist_Dp_phoneme, pca_res_phoneme,  class_true = data_phoneme$class)
# tecator
kmeans_clus_ACPF(dist_Dp_tecator, pca_res_tecator, class_true = data_tecator$class)
# growth
kmeans_clus_ACPF(dist_Dp_growth, pca_res_growth, class_true = data_growth$class)









