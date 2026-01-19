# Funciones auxiliares ----

rho2lambda <- function(p){
  L <-  p / sqrt(1-p**2)
  return(L)
}

lambda2rho <- function(L){
  p <- L / sqrt(1 + L**2)
  return(p)
}

get_obs <- function(censo, out=c(8, 15), X=Xtmp, y=ym, reg=region){
  master <- which(!is.na(y))
  group <- as.numeric(reg)
  missing  <- which(is.na(y))
  observed <- which(!is.na(y))
  ytest <- NA
  idx <- NA
  n <- length(master)
  if(!censo){
    X <- X[master, ]
    group <- reg[master]
    idx <- which(region[master]%in%out)
    idx <- sort(idx, decreasing = F)
    y <- y[master]
    ytest <- y[idx]
    y[idx] <- NA
    missing  <- which(is.na(y))
    observed <- which(!is.na(y))
  }
  yobs <- y[observed]  
  ymis <- y[missing]
  Xobs <- X[observed, ] 
  Xmis <- X[missing, ]
  output <- list(missing=missing, observed=observed, y_obs=yobs, y_mis=ymis, 
                 X_obs=Xobs, X_mis=Xmis, group=group, idx=idx, master=master,
                 y_test=ytest, K=length(unique(reg)))
  return(output)
}

get_coefs <- function(a, intercept=FALSE){
  tmp <- coef(a)
  tmp[which(is.na(tmp))] <- 0.001
  if(!intercept){
    tmp <- tmp[-1]
  }
  return(tmp)
}

metrics <- function(fitted, actual){
  mae <- mean(abs(actual - fitted))
  mse <- mean((actual - fitted)**2)
  corr<- cor(actual, fitted)
  mape<- 1/length(actual) * sum(abs(fitted - actual)/actual)
  output <- c(corr, mae, sqrt(mse), mape)
  names(output) <- c('Corr.', 'MAE', 'RMSE', 'MAPE')
  return(output)
}

delta <- function(tau2, c2){
  d <- sqrt(tau2) * sqrt(2*c2*log(sqrt(c2)) / (c2 - 1))
  return(d)
}


## SSVS ----

posterior_prob <- function(VB_fit){
  posterior_inclusion <- VB_fit$draws("m_ind")
  Mout <- apply(posterior_inclusion, 1, paste, collapse = "") # Etiquetas únicas para cada configuración
  # Inicializar contenedores para los conteos y las configuraciones de modelos
  Nm <- numeric() # Conteos para cada configuración de modelo
  Mindex <- character() # Configuraciones únicas
  # Contar la ocurrencia de cada configuración única
  while (length(Mout) > 0) {
    Nm <- c(Nm, sum(Mout %in% Mout[1]))
    Mindex <- c(Mindex, Mout[1])
    Mout <- Mout[Mout != Mout[1]] # Eliminar configuraciones procesadas
  }#
  #Calcular las probabilidades posteriores
  Pm <- Nm / sum(Nm)
  # Ordenar por probabilidad posterior en orden descendente
  ord <- order(Pm, decreasing = TRUE)
  Pm <- Pm[ord]
  Mindex <- Mindex[ord]
  # Combinar en un data frame para visualizar mejor
  Mprob <- data.frame(Probabilidad_Posterior = round(Pm, 4), Configuracion_Modelo = Mindex)
  return(Mprob)
}

marginal_pi <- function(fit){
  #marginal posterior inclusion
  posterior_inclusion <- fit$draws("m_ind")
  ## each variable:
  nrow(posterior_inclusion)
  marginal <- round(colSums(posterior_inclusion) / nrow(posterior_inclusion), 3) * 100
  return(marginal)
}


