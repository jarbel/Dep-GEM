#' @export
plot_diversity_fun = function(X, Xs, D_store, D_star, D_data){
  par(mar=c(4,4,.5,.3), mfrow = c(1,1), mgp = c(3, 1, 0))
  plot(X,D_data, las = 1,
       pch=16, ylim = c(0,max(D_data,apply(D_star,1, q_975))), xlim = c(0, max(Xs)), 
       axes = FALSE, xlab="TPH", ylab="Shannon diversity")
  grid()
  axis(1)
  axis(2, las=1)
  
  #   points(X,apply(D_store,1, mean), col=1:length(X),pch=16)
  lines(Xs,apply(D_star,1, mean), col=1, lwd = 1.5)
  lines(Xs, apply(D_star,1, q_025), col=1, lty = 2)
  lines(Xs, apply(D_star,1, q_975), col=1, lty = 2)
}

