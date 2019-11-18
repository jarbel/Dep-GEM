#' #' #' @export
#' Y <- Y[1:5, 1:20]
#' output_gibbs <- gibbs(n.iter=1000, Y = Y, X = X_jitter[1:5])
#' names(output_gibbs)
#' output_gibbs$sigma_Z_store
#' output_gibbs$lambda_store
#' output_gibbs$M_store
#' output_gibbs$d_store
#' output_gibbs$D_store
#' output_gibbs$log_a_store
#' 
#' 
#' output_gibbs$p_store
#' output_gibbs$V_store
#' 
#' 
#' zz <- output_gibbs$Z_store # I think this is it
#' zz[,,100] %>% apply(1, sum) # nope
#' 
#' zz <- output_gibbs$V_store # I think this is it
#' zz[,,100] %>% apply(1, sum) #nope
#' 
#' zz <- output_gibbs$p_store # I think this is it
#' zz[,,100] %>% apply(1, sum)  # yay
#' zz
#' 
#' get_estimates <- function(my_matrix) {
#'   my_matrix %>% apply(1, shannon)
#' }
#' apply(zz, 3, get_estimates) %>% apply(1, mean)
#' 
#' names(output_gibbs)
#' 
#' 
