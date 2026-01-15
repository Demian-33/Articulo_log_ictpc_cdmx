# =============== Analisis a nivel hogares del Ing. cor. total per capita =================

library(dplyr)
library(tidyr)
library(ggplot2)
library(cmdstanr)
library(stringr)

# (0) Datos para el modelo ====
rm(list=ls())
ruta <- r'(C:\Maestría_CP\Tesis\Articulo)'
setwd(ruta)
load('./Xyregion-CDMX.Rdata')
source('./funciones-final.R')

muni <- c(
  'Azcapotzalco', 'Coyoacán', 'Cuajimalpa de Morelos', 'Gustavo A. Madero',
  'Iztacalco', 'Iztapalapa', 'La Magdalena Contreras', 'Milpa Alta',
  'Álvaro Obregón', 'Tláhuac', 'Tlalpan', 'Xochimilco', 
  'Benito Juárez','Cuauhtémoc', 'Miguel Hidalgo', 'Venustiano Carranza') %>% factor

## (0.1) Definir areas pequeñas ====
group <- as.numeric(region)
K <- length(unique(group))

## (0.2) Covariables categoricas (no ordenadas) a matriz diseño ====
## no se consideran en el analisis
# for(i in 1:ncol(Xcat)){
#   Xcat[, i] <- as.factor(as.integer(Xcat[, i]))
# }
# tail(tibble(Xcat)) # ya todas son factor
# Xcat <- model.matrix( ~ 0 + . , data=Xcat) # generar a matriz diseño

## (0.3) Remover combinaciones lineales y correlaciones ====
Xtmp <- cbind(Xcon, Xbin)
a <- caret::findLinearCombos(Xtmp)
Xtmp <- Xtmp[, -c(a$remove)]
b <- caret::findCorrelation(cov(Xtmp), verbose=T, exact=T, cutoff = 0.8)
Xtmp <- Xtmp[, -c(b)]

## (0.4) Componentes principales ====
# Eigen <- eigen(cov(Xcon))
# cumsum(Eigen$values) / ncol(Xcon) # 26: 95%
# Xpc <- Xcon %*% Eigen$vectors[, 1:26]
# Xtmp <- cbind(Xpc, Xbin)

# (1) Modelo para predicción ====
# Un intercepto puede suplir la necesidad de una a priori jerárquica en varios interceptos
file <- file.path('C:/Maestría_CP/Tesis/Articulo', 'cmdstanr', 'LogSN-rho01-SSVS-oneInt.stan')
mod <- cmdstan_model(file, compile=T) # compile for hpp

## (1.1) Seleccionar datos ====
## censo=F, entonces se hace pronóstico en las alcaldias 8 y 15
## los argumentos ya están por defecto, vea get_obs
datos <- get_obs(censo=F)

## datos para Stan
dat <- list(
  K = K, 
  N_obs = length(datos$y_obs),
  N_mis = length(datos$y_mis),
  p = ncol(datos$X_obs),
  y_obs = datos$y_obs,
  X_obs = datos$X_obs,
  X_mis = datos$X_mis,
  group = datos$group[datos$observed],
  group_mis = datos$group[datos$missing],
  # SSVS. Stan parametriza la normal con la desv. estándar
  tau=1/300,      # El 'pico' es N(0, sd=1/300)
  c=3000,         # La 'losa' es N(0, sd=100)
  epsilon=0.01    # aceptamos |beta_{i}| > 0.01
)

## (1.2) valores iniciales ====
## para beta: se usan minimos cuadrados
tmp <- data.frame(yobs=datos$y_obs, datos$X_obs)
r0 <- lm(yobs ~ ., data=tmp)
y_mishat <- predict(r0, newdata=data.frame(datos$X_mis))
## para rho, mu y sigma se usan valores iniciales de la normal sesgada
r1 <- coef(sn::selm(yobs ~ 1, data=tmp), param.type = 'DP')
mu0  <- r1[1]
sigma0  <- r1[2]
rho0 <- lambda2rho(r1[3])

## datos para Stan
set.seed(2)
init_list <- list(
  b = get_coefs(r0, F),
  sigma=sigma0,
  pr=rep(0.5, ncol(datos$X_mis)),
  mu=mu0,
  rho=rep(rho0, times=K),
  y_mis=y_mishat
)

## (1.3) Correr modelo con VB ====
t0 <- proc.time()
fit_vb <- mod$variational(
  data = dat,
  seed = 123,
  init = list(init_list),
  grad_samples = 1, # mas muestras, mas lento
  draws = 1000,
  adapt_iter=1000,
  adapt_engaged=F, # buen ajuste con eta=0.1 (tamaño de paso)
  eta=0.1,
  iter=75000,
  tol_rel_obj=1e-6,
  algorith='meanfield', # o fullrank, mas lento
)
tproc_vb <- proc.time()-t0 # Tiempo de proceso

## (1.4) Correr modelo con HMC ====
## es mas lento este ajuste, no correr
t0 <- proc.time()
fit_hmc <- mod$sample(
  output_dir = r'(C:\Maestría_CP\Tesis\Articulo-GitHub)',
  data = dat,
  seed = 123,
  init = list(init_list),
  chains = 1,
  thin = 2,
  iter_warmup = 10,
  iter_sampling = 5,
  refresh = 1000
)
tproc_hmc <- proc.time()-t0 # Tiempo de proceso
getwd()
save(fit_hmc, file = r'(C:\Maestría_CP\Tesis\Articulo-GitHub\test.Rdata)')

## (1.4) Precisión y tiempo (predicción) ====
tproc_vb; tproc_hmc
fitted_vb0 <- apply(exp(fit_vb$draws('y_mis')), 2, mean)
fitted_vb <- apply(fit_vb$draws('y_mis'), 2, mean)
x11(); plot(fitted_vb0, fitted_vb1); grid()
abline(a=0, b=1)
metrics(apply(fit_vb$draws('y_mis'), 2, mean), datos$y_test)
metrics(fitted_vb0, exp(datos$y_test))
# metrics(fitted_vb1, exp(datos$y_test))
metrics(exp(y_mishat + 0.524/2), exp(datos$y_test))

summary(r0)

fitted_hmc <- fitted_vb
# fitted_vb <- fitted_hmc
fitted_hmc <- apply(drop(exp(fit_hmc$draws('y_mis'))), 2, mean)
m0 <- c(metrics(fitted_vb0, exp(datos$y_test)), tproc_vb[3])
m1 <- c(metrics(fitted_hmc, exp(datos$y_test)), tproc_hmc[3])
metrics((fitted_vb), (datos$y_test))
metrics(exp(y_mishat), exp(datos$y_test))
metrics((y_mishat), (datos$y_test))

round(m0, 4)
round(m1, 4)

m <- c('Corr.', 'MAE', 'RMSE', 'MAPE', 'Tiempo (s)')
metricas_out <- tibble(m, m0, m1) %>%
  kable(format='latex', digits = 4)

writeLines(metricas_out, con = './Figuras/metricas-vb-hmc.txt')
metricas_out <- readLines(con = './Figuras/metricas-vb-hmc.txt')
metricas_out <- metricas_out[-c(1, 2, 3, 4, 5, 16)]
writeLines(metricas_out, con = './Figuras/metricas-vb-hmc.txt')

print(fit_vb$summary(c('rho', 'mu', 'mu_0', 'sigma_mu')), n=35)
print(fit_hmc$summary(c('rho', 'mu', 'mu_0', 'sigma_mu')), n=35)

## (1.5) Grafico observados vs ajustados =====
fitted_hmc <- rep(NA, length=dat$N_mis) # por si no se corre
fitted_vb <- rep(NA, length=dat$N_mis) # por si no se corre
fitted_hmc <- apply(drop(fit_hmc$draws('y_mis')), 2, mean)
cor(fitted_hmc, datos$y_test)
dftmp <- data.frame(
  Alcaldía=factor(muni[region[datos$master[datos$idx]]]),
  Observados=datos$y_test,
  Ajustados_vb=fitted_vb,
  Ajustados_hmc=fitted_hmc)

g0 <- ggplot(dftmp, aes(x=Ajustados_vb, y=Observados, color=Alcaldía, shape=Alcaldía, group=Alcaldía)) + 
  geom_abline(intercept = 0, slope=1, linewidth=1.2, color='gray5', lty=3) +
  geom_point(size=2) +
  scale_color_manual(values=c('Milpa Alta'='gray25', 'Miguel Hidalgo'='gray75'))+
  geom_smooth(method = 'lm', formula=y~x, se=F, color='gray25')+
  theme_minimal() + 
  theme(legend.position = 'top', text=element_text(family='serif', size=11)) + 
  labs(x='log-ICTPC ajustado', y='log-ICTPC observado')
x11(); print(g0)

g1 <- ggplot(dftmp, aes(x=Ajustados_hmc, y=Observados, color=Alcaldía, shape=Alcaldía, group=Alcaldía)) + 
  geom_abline(intercept = 0, slope=1, linewidth=1.2, color='gray5', lty=3) +
  geom_point(size=2) +
  scale_color_manual(values=c('Milpa Alta'='gray25', 'Miguel Hidalgo'='gray75'))+
  geom_smooth(method = 'lm', formula=y~x, se=F, color='gray25')+
  theme_minimal() + 
  theme(legend.position = 'top', text=element_text(family='serif', size=11)) + 
  labs(x='log-ICTPC ajustado', y='log-ICTPC observado')
x11(); print(g1)

ggsave('./Figuras/LogSN-CDMX-vb.png', plot=g0, width=15, heigh=15, units='cm')
ggsave('./Figuras/LogSN-CDMX-hmc.png', plot=g1, width=15, heigh=15, units='cm')

# (2) Modelo para estimación ====
# ¿unico o varios interceptos?
# ¿que valores iniciales tomar?
# ¿pr = 0.1 implica seleccionar menos covariables?
file <- file.path('C:/Maestría_CP/Tesis/Articulo', 'cmdstanr', 'LogSN-rho01-SSVS-genq.stan')
file <- file.path('C:/Maestría_CP/Tesis/Articulo', 'cmdstanr', 'LogSN-rho01-SSVS-genq-oneInt.stan')
mod <- cmdstan_model(file, compile=T) 
mod$print()

## (2.1) Seleccionar datos ====
caret::findLinearCombos(Xtmp)$remove
Xtmp <- Xtmp[, -caret::findLinearCombos(Xtmp)$remove]
Xtmp <- Xtmp[, -100] # PC
Xtmp <- Xtmp[, -110] # sin PC
dim(Xtmp)
datos <- get_obs(none=TRUE, y=ym, region=region, X=Xtmp) # solo observados
dat <- list(
  K = K,
  N_obs = length(datos$y_obs),
  N_mis = length(datos$y_mis),
  p = ncol(datos$X_obs),
  y_obs = datos$y_obs,
  X_obs = datos$X_obs,
  X_mis = datos$X_mis,
  X_all = Xtmp,
  group = datos$group[datos$observed],
  groupmis = datos$group[datos$missing],
  group_all = region,
  # SSVS; Stan parametriza la normal con la desv. estándar
  tau=1/300,  
  c=sqrt(10*300**2),     
  epsilon=0.01
)

# delta(dat$tau**2, dat$c**2)

## (2.2) Valores iniciales ====
## usando lm
## expande algunos betas
tmp <- data.frame(y_obs=datos$y_obs, datos$X_obs)
caret::findLinearCombos(tmp[, -1])$remove
# tmp <- tmp[, -caret::findLinearCombos(tmp)$remove]
ncol(tmp[, -1]); qr(tmp[, -1])$rank

dim(tmp)
r0 <- lm(y_obs ~ ., data=tmp)
# r0 <- step(r0)
## usando sn::selm
## parece que no expande mucho...
x0 <- as.matrix(cbind(1, tmp[, c(2:45)]))
r0 <- sn::selm.fit(x0, y=tmp[, 1])
max(r0$param$dp[-1]); min(r0$param$dp)
x11(); plot(r0$param$dp, pch=19); grid()
set.seed(1)
# eta1 <- list(list(X=tmp[, -1], model='BRR'))
# r1 <- BGLR(y=tmp[, 1], ETA=eta1, nIter=50000, burnIn=25000, thin=5, verbose = F)
# b01 <- r1$ETA[[1]][["b"]]
eta2 <- list(list(X=tmp[, -1], model='FIXED'))
r2 <- BGLR(y=tmp[, 1], ETA=eta2, nIter=50000, burnIn=25000, thin=5, verbose = F)
b02 <- r2$ETA[[1]][["b"]]
min(b02); max(b02)
min(coef(r0), na.rm = T); max(coef(r0)[-1], na.rm = T)
x11()
plot(b02, pch=19, ylim=c(-5, 5), col='gray0')
points(get_coefs(r0, F), col='gray50', pch=19)
points(fit_vb$summary('b')$mean, col='gray25', pch=25)
legend('topleft', legend=c('Bayes flat', 'OLS'), col=c('gray0', 'gray50'), pch=19)
grid()
# points(b02, col=4, pch=19)
# points(round(beta_out$Media, 2), pch=19, col='red')

mu0 <- c()
rho0 <- c()
dftmp <- tibble(ym, group)[datos$master, ]
# for(i in 1:K){
#   tmp <- subset(dftmp, group==i)
#   r0 <- sn::selm(ym ~ 1 , data=tmp)
#   mu0 <- c(mu0, coef(r0, param.type = 'DP')[1])
#   # rho0 <- c(rho0, lambda2rho(coef(r0, param.type = 'DP')[3]))
# }
r3 <- coef(sn::selm(ym ~ 1 , data=dftmp), param.type='DP')
mu0 <- r3[1]
sigma0 <- r3[2]
rho0 <- lambda2rho(r3[3])

min(b02); max(b02)
min(coef(r0), na.rm = T); max(coef(r0)[-1], na.rm = T)

set.seed(2)
init_list <- list(
  b = get_coefs(r0, F), # nice without PC
  # b=b02, # nice with PC
  # b = rnorm(dat$p, mean=0, sd=0.01),
  pr=rep(0.5, ncol(datos$X_mis)), # 0.5
  sigma=sigma0,
  mu=mu0,
  # mu=rep(mu0, times=K),
  rho=rep(rho0, times=K)
)

## (2.3) Correr el modelo con VB ====
t0 <- proc.time()
fit_vb <- mod$variational(
  data = dat,
  grad_samples = 1,
  eval_elbo = 100,
  seed = 123,
  init = list(init_list),
  draws = 1000,
  adapt_iter=5000,
  adapt_engaged=F, # se encontró óptimo eta=0.1
  eta=0.1, # 0.05
  iter=75000,
  tol_rel_obj=1e-6,
  algorith='meanfield', # fullrank
)
tproc_vb <- proc.time()-t0 # Tiempo de proceso
print(fit_vb$summary(c('rho', 'sigma', 'mu')), n=38)
print(fit_vb$summary(c('b')), n=38)
sort(round(fit_vb$summary(c('b'))$mean, 2))
x11(); boxplot(round(fit_vb$summary(c('b'))$mean, 2), col='gray50'); grid()

## (2.4) Correr el modelo con HMC ====
# t0 <- proc.time()
# fit_hmc <- mod$sample(
#   output_dir = r'(C:\Users\Lesau\Documents\EntornosR-Articulo)',
#   data = dat,
#   seed = 123,
#   init = list(init_list),
#   chains=1, 
#   thin=5,
#   iter_warmup = 45000,
#   iter_sampling= 5000,
# )
# tproc_hmc <- proc.time()-t0 # Tiempo de proceso: 89.8 minutos
# writeLines(paste0(tproc_hmc, sep='-'), con = 'Entornos-Articulo/Tiempo-hmc-censo20.txt')
# save(fit_hmc, file='Entornos-Articulo/fit_hmc-censo20_burnin_45k_iter5k_thin5_.Rdata')
# mod # K interceptos, no pool
# valores iniciales de sn::selm(y_obs ~ 1), DP
load('Entornos-Articulo/fit_hmc-censo20_burnin_45k_iter5k_thin5_.Rdata')
print(fit_hmc$summary(c('rho', 'sigma', 'mu')), n=38)
sort(round(fit_hmc$summary(c('b'))$mean, 2))

## (2.5) Precisión y tiempo (ajustados) ====
cat("\n","Tiempo de ejecución:",tproc_vb[3],"(elapsed)\n")
out <- which(!is.na(ym))
fitted_vb <- apply(fit_vb$draws('y_all'), 2, mean)
# fitted_hmc <- apply(fit_vb$draws('y_all'), 2, mean)
metrics(fitted_vb[out], datos$y_test)
metrics(fitted(r1), datos$y_test)

## (2.6) Gráfico para las estimaciones de rho, sigma, mu ====

# g_rho <- mcmc_intervals(fit_vb$draws('rho')) + 
#   scale_y_discrete(labels = parse(text = paste0("rho[", 1:16, "]")))
# g_sigma <- mcmc_hist(fit_vb$draws('sigma'), alpha = 0.2) + 
#   labs(x = parse(text = 'sigma'))
# g_mu <- mcmc_intervals(fit_vb$draws('mu')) +
#   scale_y_discrete(labels = parse(text = paste0("mu[", 1:16, "]")))

# x11(); print(g_rho)
# x11(); print(g_sigma)
# x11(); print(g_mu)

# ggsave(filename='int_rho.png', plot=g_rho, width=15, heigh=15, units='cm')
# ggsave(filename='hist_sigma.png', plot=g_sigma, width=15, heigh=15, units='cm')
# ggsave(filename='int_mu.png', plot=g_mu, width=15, heigh=15, units='cm')

# (3) Estimación a posteriori ====

## (3.1) Estimación de beta ====
# podemos hacer algo mas completo
# mostrar el modelo mas popular
# y las probabilidades marginales de aparicion:

beta_ <- fit_vb$draws("b")
Media <- apply(beta_, 2, mean)
Error <- apply(beta_, 2, sd)
Lo025 <- apply(beta_, 2, quantile, 0.025)
Up975 <- apply(beta_, 2, quantile, 0.975)
Frecuencia <- marginal_pi(fit_vb)
# Frecuencia <- 100*apply(fit_vb$draws("m_ind"), 2, mean)
# Frecuencia <- 100*apply(fit_vb$draws("pr"), 2, mean)
Parametro <- paste0("$\\beta_{", 1:ncol(dat$X_obs), "}$")

beta_out <- tibble(Parametro, Media, `Error est.`=Error, `0.025 \\%`=Lo025, `0.975 \\%`=Up975, `Frecuencia \\%`=Frecuencia)
print(beta_out, n=Inf)
beta_out %>% filter(`Frecuencia \\%` > 75) %>%
  print(n=Inf)

table(beta_out$`Frecuencia \\%`)
sum(table(beta_out$`Frecuencia \\%`[beta_out$`Frecuencia \\%`>75]))

beta_out %>% filter(`Frecuencia \\%` > 75) %>%
  kable(
    format='latex', escape=F, digits = 4,
    label='tab:estimaciones-beta',
    caption = 'Elaboración propia basada en la muestra \\underline{a posteriori}.'
  )

print(beta_out, n=Inf)
sort(round(beta_out$Media, 2))

# usando los valores de lm
pos <- c(42, 93, 69, 43, 25, 52, 83)
neg <- c(102, 109, 104, 80, 18, 16, 89)

beta_out[pos, ]
beta_out[neg, ]

colnames(Xtmp)[pos]
colnames(Xtmp)[neg]

table(Xtmp$sanit)
load('./Figuras/TablaCovariables.Rdata')
subset(tabla_cov, Variable%in%colnames(Xtmp)[pos])
subset(tabla_cov, Variable%in%colnames(Xtmp)[neg])

## (3.2) Covariables seleccionadas ====
## hay que tomar las 81 binarias que fueron seleccionadas
## (1) tomar el modelo mas popular
## (2) poner sus probabilidades marginales
load('./Figuras/TablaCovariables.Rdata')
ccon <- tabla_cov %>% filter(Tipo=='Continua') %>% select(Descripción) %>% unlist() %>% as.character()
# ccon <- paste0('{\\footnotesize ', 1:52, '. ', ccon, '} \\\\')
ccon <- paste0(1:52, '. ', ccon)
writeLines(ccon, con='./Figuras/lista_con.txt')

library(stringr)
idx0 <- 1:ncol(Xtmp) # todos
idx1 <- str_split(posterior_prob(fit_vb)[1, 2], '', simplify = T) %>% unlist() %>% as.numeric() %>% as.logical()
# idx2 <- which(beta_out$`Frecuencia \\%` > 75 & !str_starts(colnames(Xtmp), "[0-9]"))

# seleccionar los indices en idx1 tales que no sean componente principal
# -26 para mover la PC y +52 para recuperar la variable original
idx3 <- idx0[!str_starts(colnames(Xtmp), "[0-9]")] - 26 + 52

Parametro <- paste0("$\\beta_{", 1:length(idx0), "}$")
# no se convierte bien en pandoc
# Parametro <- paste0('\\makecell[l]{', '$\\beta_{', 1:length(idx0), '}$',
#                     '\\\\ (', sprintf('%.2f', beta_out$`Frecuencia \\%`[idx0]), '\\%)}')
Parametro <- paste0('$\\beta_{', 1:length(idx0), '}$ (', sprintf('%.2f', beta_out$`Frecuencia \\%`[idx0]), '\\%)')
Media_Error <- c()
Lo_Up <- c()
beta_out
for(k in 1:length(idx0)){
  a0 <- sprintf('%.4f', beta_out$Media[k])
  a1 <- sprintf('%.4f', beta_out$`Error est.`[k])
  # Media_Error[k] <- paste('\\makecell[l]{', a0, '\\\\', '(', a1, ')', '}', sep='')
  Media_Error[k] <- paste(a0, ' (', a1, ')', sep='')
  b0 <- sprintf('%.4f', beta_out$`0.025 \\%`[k])
  b1 <- sprintf('%.4f', beta_out$`0.975 \\%`[k])
  Lo_Up[k] <- paste('(', b0, ', ', b1, ')', sep='')
}

# param_out2 <- tibble(
#   `\\makecell[l]{Parámetro\\\\(Frecuencia)}`=Parametro,
#   `\\makecell[l]{Media\\\\(Error est.)}`= Media_Error,
#   `\\makecell[l]{Intervalo c. \\\\ (2.5 \\%, 97.5 \\%)}`= Lo_Up,
#   Descripción = c(rep('', times=26), tabla_cov$Descripción[idx3])
# )

param_out2 <- tibble(
  `Parámetro (Frecuencia)`=Parametro,
  `Media (Error est.)`= Media_Error,
  `Intervalo cred. (2.5 \\%, 97.5 \\%)`= Lo_Up,
  Descripción = tabla_cov$Descripción[idx0]
  # Descripción = c(rep('', times=26), tabla_cov$Descripción[idx3]) # PC
)

param_out2 <- param_out2[idx1, ]

print(param_out2, n=60)

param_out2_ <- kable(
  param_out2,
  format='latex', escape = F,
  label='covariables-beta',
  caption = 'Estimaciones y lista de covariables incluidas, únicamente
  se muestran aquellas covariables con frecuencias de aparición mayores al 75\\%.
  Fuente: elaboración propia basada en la muestra \\underline{a posteriori}.'
)

# write.table(file='covariables-beta.txt', param_out2_, row.names = F, col.names = F)
writeLines(param_out2_, con='./Figuras/covariables-beta.txt')
param_out2_ <- readLines('./Figuras/covariables-beta.txt')
param_out2_ <- param_out2_[-c(1:10, length(param_out2_)-2, length(param_out2_)-1, length(param_out2_))]
param_out2_[seq(2, 74, by=2)] <- '\\cline{1-3}'
writeLines(param_out2_, con='./Figuras/covariables-beta.txt')

## (3.2.1) Grafico de mosaico ====
all_pr <- posterior_prob(fit_vb)[, 1]
# all_pr <- posterior_prob(fit_hmc)[, 1]
all_pr[1:10]
cumsum(all_pr)[1:10]
x11(); plot(all_pr[1:10], col=4, type='b'); grid()
pr <- str_split(posterior_prob(fit_vb)[1, 2], '', simplify = T) %>% unlist() %>% as.numeric()
# pr <- str_split(posterior_prob(fit_hmc)[1, 2], '', simplify = T) %>% unlist() %>% as.numeric()
# colnames(Xtmp)[which(pr == 1)][-c(1:26)] # PC
colnames(Xtmp)[which(pr == 1)]
tabla_cov[which(pr == 1), ]

beta_est_nom <- tibble(
  fit_vb$summary('b')[which(pr == 1), 1:2],
  colnames(Xtmp)[which(pr == 1)],
  coef(r0)[-1][which(pr == 1)])
print(beta_est_nom%>%arrange(mean), n=Inf)


vb <- cbind(posterior_prob(fit_vb)[1:5, ], Metodo='VB')
# hmc <- cbind(posterior_prob(muestra_HMC), Metodo='HMC')
# vb <- cbind(posterior_prob(fit_hmc), Metodo='HMC')
# tmp <- rbind(vb, hmc)
tmp <- rbind(vb, vb)
colnames(tmp) <- c('Prob', 'Config', 'Método')

tmp <- tmp %>%
  mutate(Prob = factor(Prob), Config = str_split(Config, "")) %>%
  unnest_wider(Config, names_sep = "") %>%
  rename_with(~ paste0("beta_", seq_along(.)), starts_with("Config")) %>%
  pivot_longer(starts_with("beta_"), names_to = "beta", values_to = "included") %>%
  mutate(included = as.integer(included)) %>%
  complete(Prob, beta, fill = list(included=NA)) %>%
  mutate(Método=factor(Método, levels=c('VB', 'HMC')))

g0 <- ggplot(tmp, aes(x = beta, y = Prob, fill = factor(included))) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_wrap(~Método, scales = 'free_y') +
  scale_x_discrete(labels=parse(text = paste0("beta[", 1:ncol(Xtmp), "]"))) +
  theme_minimal() +
  geom_vline(xintercept = 21.5, linewidth=1.2, col='gray25')+
  labs(fill='Inclusión', y='Prob. a posteriori', x='') +
  scale_fill_discrete(labels=c('0'='No', '1'='Sí'))+
  theme(legend.position = 'top', text=element_text(family='serif'),
        plot.margin = margin(0, 0, 0, 0),  axis.ticks.length = unit(0, "pt"))
x11(); print(g0)
# no hay mucha diferencia en la eleccion de modelos
colnames(Xcon)

## (3.3) Estimacion de rho, sigma y mu ====
rho_ <- fit_vb$draws("rho")
Media_rho <- apply(rho_, 2, mean)
Error_rho <- apply(rho_, 2, sd)
Lo025_rho <- apply(rho_, 2, quantile, 0.025)
Up975_rho <- apply(rho_, 2, quantile, 0.975)

sigma_ <- fit_vb$draws("sigma")
Media_sigma <- apply(sigma_, 2, mean)
Error_sigma <- apply(sigma_, 2, sd)
Lo025_sigma <- apply(sigma_, 2, quantile, 0.025)
Up975_sigma <- apply(sigma_, 2, quantile, 0.975)

mu_ <- fit_vb$draws("mu")
Media_mu <- apply(mu_, 2, mean)
Error_mu <- apply(mu_, 2, sd)
Lo025_mu <- apply(mu_, 2, quantile, 0.025)
Up975_mu <- apply(mu_, 2, quantile, 0.975)

# Cuando se usa pooling
# mu0_ <- fit_vb$draws("mu_0")
# Media_mu0 <- apply(mu0_, 2, mean)
# Error_mu0 <- apply(mu0_, 2, sd)
# Lo025_mu0 <- apply(mu0_, 2, quantile, 0.025)
# Up975_mu0 <- apply(mu0_, 2, quantile, 0.975)
# 
# sigma_mu_ <- fit_vb$draws("sigma_mu")
# Media_sigma_mu <- apply(sigma_mu_, 2, mean)
# Error_sigma_mu <- apply(sigma_mu_, 2, sd)
# Lo025_sigma_mu <- apply(sigma_mu_, 2, quantile, 0.025)
# Up975_sigma_mu <- apply(sigma_mu_, 2, quantile, 0.975)

Parametro <- c(
  paste0("$\\rho_{", 1:16, "}$"),
  paste0("$\\sigma$"),
  # paste0("$\\mu_{", 1:16, "}$")
  paste0("$\\mu$")
  #paste0("$\\mu_{0}$"),      # pooling
  #paste0("$\\sigma_{\\mu}$") # pooling
)

param_out <- tibble(
  Parametro=Parametro,
  Media=c(Media_rho, Media_sigma, Media_mu), # , Media_mu0, Media_sigma_mu
  `Error est.`=c(Error_rho, Error_sigma, Error_mu), # , Error_mu0, Error_sigma_mu
  `0.025 \\%`=c(Lo025_rho, Lo025_sigma, Lo025_mu), # , Lo025_mu0, Lo025_sigma_mu
  `0.975 \\%`=c(Up975_rho, Up975_sigma, Up975_mu), # , Up975_mu0, Up975_sigma_mu
)

print(param_out, n=Inf)
param_out %>% kable(
  format='latex', escape = F, digits = 4,
  label='tab:estimaciones',
  caption = 'Elaboración propia basada en la muestra \\underline{a posteriori}.'
)

K <- 16
# par1 <- Parametro[1:(K+1)]
# par2 <- c(NA, Parametro[(K+1):length(Parametro)])
Media_Error <- c()
Lo_Up <- c()
for(k in 1:(K+1+1)){
  a0 <- sprintf('%.4f', param_out$Media[k])
  a1 <- sprintf('%.4f', param_out$`Error est.`[k])
  # Media_Error[k] <- paste('\\makecell[l]{', a0, '\\\\', '(', a1, ')', '}', sep='')
  Media_Error[k] <- paste(a0, ' (', a1, ')', sep='')
  b0 <- sprintf('%.4f', param_out$`0.025 \\%`[k])
  b1 <- sprintf('%.4f', param_out$`0.975 \\%`[k])
  Lo_Up[k] <- paste('(', b0, ', ', b1, ')', sep='')
}

param_out2 <- tibble(
  Parámetro=Parametro,
  `\\makecell[l]{Media\\\\(Error est.)}`= Media_Error,
  `(2.5 \\%, 97.5 \\%)`= Lo_Up
)

# kable(
#   list(param_out2[1:(1+K),], param_out2[(2+K):(2*K), ]),
#   format='latex', escape = F,
#   label='tab:estimaciones',
#   caption = 'Elaboración propia basada en la muestra \\underline{a posteriori}.'
# )

param_out2[(K+1+1), ] <- 'foo'
param_out2 <- cbind(param_out2[1:(K+1), ], param_out2[(K+2):(2*K+2), ])


param_out2 <- tibble(c(as.character(muni), NA, NA), param_out2)


param_out2 %>% kable(
  format='latex',
  escape='F'
)

## (3.4) Gini-Lorenz a posteriori ====
# 0: ámbito urbano
# 1: ámbito rural
library(ineq)
library(dineq)
out <- which(is.na(ym))
fitted_vb <- apply(exp(fit_vb$draws('y_all')), 2, mean)
fitted_vb <- apply(drop(exp(fit_hmc$draws('y_all'))), 2, mean)

dftmp <- data.frame(
  Ajustados=fitted_vb[out],
  Contexto=Xide$rururb[out],
  Factor=Xide$factor[out]
)

lc <- Lc(dftmp$Ajustados, dftmp$Factor)
lc_urb <- Lc(dftmp$Ajustados[dftmp$Contexto==0], dftmp$Factor[dftmp$Contexto==0])
lc_rur <- Lc(dftmp$Ajustados[dftmp$Contexto==1], dftmp$Factor[dftmp$Contexto==1])

round(Gini(dftmp$Ajustados),6)
round(gini.wtd((dftmp$Ajustados), dftmp$Factor),6)
gini_decomp((dftmp$Ajustados), dftmp$Contexto, dftmp$Factor)$gini_group$gini_group

len0 <- c(nrow(dftmp), table(dftmp$Contexto))+1L
len1 <- c()
for(k in seq_along(len0)){
  len1 <- c(len1, rep(c('Ambos', 'Urbano', 'Rural')[k], times=len0[k]))
}
len1 <- factor(len1)

LG_all <- data.frame(
  p=c(lc$p, lc_urb$p, lc_rur$p), L=c(lc$L, lc_urb$L, lc_rur$L), Ámbito=len1)

LG_plot <- ggplot(LG_all, aes(x=p, y=L, group=Ámbito, linetype=Ámbito, color=Ámbito, label=Ámbito)) + 
  geom_line(linewidth=1.2) + 
  # geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = c('Urbano'='gray75', 'Rural'='gray50', 'Ambos'='gray0'))+
  coord_fixed(ratio = 1) +          # ensures a square panel
  geom_abline(intercept = 0, slope=1, linewidth=1.2, color='gray75', lty=1) + 
  theme_bw() + 
  theme(text=element_text(family='serif', size=11), legend.position = 'top')

x11(); plot(LG_plot)
ggsave(filename='./Figuras/Lorenz-CDMX.png', plot=LG_plot, width=15, heigh=15, units='cm')

## (3.6) Porcentaje de personas con las LPI y LPEI ====

# indices <- which(!is.na(ym))
# sum((Xide$factor * Xide$tam_hog)[indices])  # provenientes de la ENIGH
# sum((Xide$factor * Xide$tam_hog)[-indices]) # provenientes del censo

LPI_urb <- 4564.97
LPI_rur <- 3296.92
LPEI_urb <- 2354.65
LPEI_rur <- 1800.55

out <- which(is.na(ym))
fitted_vb <- apply(exp(fit_vb$draws('y_all')), 2, mean)
# fitted_vb <- exp(apply(drop(fit_hmc$draws('y_all')), 2, mean))
# nuevo enfoque?, aqui podemos usar integracion MC

dftmp2 <- tibble(Ajustados=fitted_vb[out],
                 Ámbito=Xide$rururb[out],
                 Factor=Xide$factor[out],
                 Tamaño=Xide$tam_hog[out],
                 Alcaldia=muni[region][out])

dftmp2 <- dftmp2 %>%
  mutate(LPI = case_when(
    `Ámbito` == 0 & Ajustados < LPI_urb ~ 1L,
    `Ámbito` == 1 & Ajustados < LPI_rur ~ 1L,
    TRUE                                ~ 0L)) %>%
  mutate(LPEI = case_when(
    `Ámbito` == 0 & Ajustados < LPEI_urb ~ 1L,
    `Ámbito` == 1 & Ajustados < LPEI_rur ~ 1L,
    TRUE                                ~ 0L))

dftmp2 <- dftmp2 %>%
  group_by(Alcaldia) %>%
  summarise(Total = sum(Factor*Tamaño),
            Total_LPI = sum(Factor*LPI*Tamaño), Total_LPEI = sum(Factor*LPEI*Tamaño),
            Porcentaje_LPI = 100*Total_LPI/Total, Porcentaje_LPEI = 100*Total_LPEI/Total)


sum(dftmp2$Total_LPI) / sum(dftmp2$Total) # 25.7 % 
sum(dftmp2$Total_LPEI) / sum(dftmp2$Total) # 4.5 % 

# ordenar de acuerdo a muni
dftmp2 <- dftmp2[match(muni, dftmp2$Alcaldia), ]

dftmp2 %>%
  select(Alcaldia, starts_with('Porcentaje')) %>%
  kable(format = 'latex', digits=3,
        label='tab:porcentajes',
        caption='Elaboración propia basada en la muestra \\underline{a posteriori}.')

# save(muni, fit_vb, dftmp2, file='./fit_vb_CPV2020.RData')

## (3.8) Cuantificar sesgo ====
## solo puede calcularse entre valores observados y ajustados
idx <- which(!is.na(ym))
# fitted_vb <- apply(drop(fit_hmc$draws('y_all')), 2, mean) # hmc
fitted_vb <- apply(exp(fit_vb$draws('y_all')), 2, mean) # vb
fitted_vb[-idx] <- -99999 # lo que sea
fitted_vb[idx] <- exp(fitted(r0)) # ajuste de r0
fitted_vb[idx] <- exp(fitted(r0) + 0.524/2)
fitted_vb[idx] <- exp(r2$yHat) # ajuste de BGLR
fitted_vb[idx] <- exp(r2$yHat + r2$varE/2) # ajuste de BGLR

dftmp <- tibble(
  ICTPC_obs = exp(ym[idx]),
  ICTPC_fit=fitted_vb[idx],
  Alcaldia=muni[region][idx])

dftmp <- dftmp %>%
  group_by(Alcaldia) %>%
  summarise(mean_ICTPC_obs = mean(ICTPC_obs, na.rm=T), # cambiar por median, esta curioso
            mean_ICTPC_fit = mean(ICTPC_fit, na.rm=T),
            diff_ICTPC = abs(mean_ICTPC_obs-mean_ICTPC_fit))

mean(dftmp$diff_ICTPC)

x11(); plot(fitted_vb[idx], exp(ym[idx]))
points(fitted_vb[idx], exp(ym[idx]), col=4)
grid()
metrics(fitted_vb[idx], exp(dat$y_obs))

# ordenar de acuerdo a muni
dftmp <- dftmp[match(muni, dftmp$Alcaldia), ]

dftmp %>%
  kable(format = 'latex', digits=3,
        label='sesgo',
        caption='Elaboración propia basada en la muestra \\underline{a posteriori}.')

# (4) Descriptivos ====

## (4.1) Tablas de estadísticos ====

datos <- list()
datos$master <- which(!is.na(ym))
dftmp <- data.frame(Observados=ym[datos$master],
                    Ámbito=Xide$rururb[datos$master],
                    Factor=Xide$factor[datos$master],
                    Alcaldia=muni[region[datos$master]])

Res_ambito <- dftmp %>%
  group_by(Ámbito) %>%
  summarise(Minimo=min(exp(Observados)),
            Mediana = median(exp(Observados)),
            Media = mean(exp(Observados)),
            Maximo=max(exp(Observados)),
            `Desv. est.`=sd(exp(Observados)),
            Conteo=n(), 
            .groups = "drop")

Res_ambito %>%
  kable(digits=3, format='latex')

Ambito_agrupado <- dftmp %>%
  mutate(Ámbito=factor(Ámbito, labels=c('Urbano', 'Rural'))) %>%
  group_by(Ámbito) 

Ambito_resumen %>% kable(format = "latex",
                         booktabs = F,
                         caption = "Resumen por Alcaldía y Ámbito",
                         escape = TRUE,
                         linesep = "") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  collapse_rows(columns = 1, valign = "middle")

Ambito_resumen %>% kable(format = "latex",
                         booktabs = TRUE,
                         caption = "Resumen por Alcaldía y Ámbito",
                         escape = TRUE,
                         linesep = "") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  collapse_rows(columns = 1, valign = "middle")


dftmp %>% group_by(Alcaldia) %>%
  summarise(Minimo=min(exp(Observados)),
            Mediana = median(exp(Observados)),
            Media = mean(exp(Observados)),
            Maximo=max(exp(Observados)),
            `Desv. est.`=sd(exp(Observados)),
            Conteo=n(), 
            .groups = "drop")

table_out <- dftmp %>%
  mutate(Ámbito=factor(Ámbito, labels=c('Urbano', 'Rural')),
         Alcaldia=factor(Alcaldia, levels=muni)) %>%
  group_by(Alcaldia, Ámbito) %>%
  summarise(Minimo=round(min(exp(Observados)), 3),
            Mediana = round(median(exp(Observados)), 3),
            Media = round(mean(exp(Observados)), 3),
            Maximo=round(max(exp(Observados)), 3),
            `Desv. est.`=round(sd(exp(Observados)), 3),
            Conteo=n(), 
            .groups = "drop") %>%
  complete(Alcaldia, Ámbito,
           fill = list(
             Minimo = NA_real_,
             Mediana = NA_real_,
             Media = NA_real_,
             Maximo = NA_real_,
             `Desv. est.` = NA_real_,
             Conteo = 0)) 
# arrange(Alcaldia, `Ámbito`) %>% #Order rows using column values

# install.packages('kableExtra')
library(kableExtra)

latex_table <- table_out %>%
  kable(format = "latex",
        booktabs = F,
        caption = "Resumen por Alcaldía y Ámbito",
        escape = TRUE,
        linesep = "")
latex_table

latex_table <- table_out %>%
  kable(format = "latex",
        booktabs = TRUE,
        caption = "Resumen por Alcaldía y Ámbito",
        escape = TRUE,
        linesep = "") %>%
  kable_styling(latex_options = c("hold_position")) %>%
  collapse_rows(columns = 1, valign = "middle")

latex_table

## (4.2) Log-ICTPC grafico de violin ====

Ambito <- ggplot(Ambito_agrupado, aes(x = Ámbito, y = Observados, fill=Ámbito)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.25, outlier.shape = 2, outlier.size = 4) +
  geom_jitter(width = 0.12, alpha = 0.2, size = 1) +
  labs(x = "Ámbito", y = 'log - ingreso corriente total pér cápita', title = "") +
  theme_minimal() + 
  theme(legend.position = 'none', text = element_text(family = "serif", size=11))
x11(); plot(Ambito)
ggsave(filename='./Box-Ambito.png', plot=Ambito,
       width = 15, height = 15, bg = "transparent", units='cm')


Ing <- ggplot(dftmp, aes(x = Ámbito, y = Observados, fill=`Ámbito`)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.25, outlier.shape = 2) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1) +
  labs(x = "Alcaldia", y = 'log-ICTPC', title = "") +
  theme_minimal() + 
  theme(legend.position = 'none', text = element_text(family = "serif"))
x11(); print(Ing)

## (4.3) Lista de covariables ====

library(readxl)
library(tidyr)
library(stringr)

load('./Xyregion-CDMX.Rdata')
descripcion <- read_excel('descripcion_variables.xlsx', range='A1:D308')
head(descripcion)

descripcion <- descripcion %>%
  # borrar espacios invisibles y convertir "" a NA
  mutate(Variable = na_if(str_squish(Variable), "")) %>%
  fill(Variable, .direction = "down") %>%
  mutate(Valores=if_else(Valores=='(1…8)', '(0…999999999)', Valores)) %>%
  mutate(Valores=factor(Valores))
descripcion

nomcon <- descripcion %>%
  filter(Variable %in% colnames(Xcon)) %>%
  distinct(Variable, .keep_all = T)
dim(Xcon)
dim(nomcon)

nombin <- descripcion %>%
  filter(Variable %in% colnames(Xbin)) %>%
  distinct(Variable, .keep_all = T)
dim(Xbin)
dim(nombin)

nomcat <- descripcion %>%
  filter(Variable %in% colnames(Xcat)) %>%
  distinct(Variable, .keep_all = T)
dim(Xcat)
dim(nomcat)

nomall <- rbind(
  nomcon %>% arrange(Variable),
  nombin %>% arrange(Variable),
  nomcat %>% arrange(Variable)
)

# Tipo <- c(rep('Con.', ncol(Xcon)), rep('Bin.', ncol(Xbin)), rep('Cat.', ncol(Xcat)))
Tipo <- c(rep('Continua', ncol(Xcon)), rep('Binaria', ncol(Xbin)), rep('Categórica', ncol(Xcat)))
nomall <- cbind(nomall, Tipo)

tabla_cov <- nomall %>%
  select(!c('Valores', 'Etiquetas')) %>%
  relocate(Tipo) %>%
  tibble %>%
  filter(Tipo %in% c('Continua', 'Binaria'))
# %>% filter(Variable %in% colnames(Xtmp))

52+81+7
tabla_cov
dim(tabla_cov)
save(tabla_cov, file='./Figuras/TablaCovariables.Rdata')


## Mas graficos ====

x11()
mcmc_intervals(fit_vb$draws("pr"))
mcmc_intervals(fit_vb$draws("rho"))
mcmc_hist(fit_vb$draws("rho"))
mcmc_hist(fit_vb$draws("z"))
mcmc_hist(fit_vb$draws("mu"))
mcmc_hist(fit_vb$draws("sigma"))

# ajustados, todo el censo
fitted_vb <- apply(fit_vb$draws("y_mis"), 2, mean)[!is.na(ym)]
metrics(fitted_vb, datos$ytest)
metrics(fitted(r1), datos$ytest)
x11(); plot(fitted_vb, datos$ytest); points(fitted(r1), datos$ytest, col='red', pch=19)
#

fitted_vb <- apply(fit_vb$draws("ymis"), 2, mean)
metrics(fitted_vb, datos$ytest)
metrics(ymishat, datos$ytest)
x11(); plot(fitted_vb, datos$ytest); points(ymishat, datos$ytest, col='red', pch=19)

# write.csv(tproc_vb, file='./tiempo_VB_iter125k_grad1_elbo100_eta01.csv')
# save(muestra_VB, file='./muestraVB_iter125k_grad1_elbo100_eta01.Rdata')


## (4.5) Gráfico de distribución del ingreso (histogramas) ----
umbral <- 210 # filtrar dos o tres ingresos muy pequeños
umbral <- 500
tmp <- which(!is.na(ym) & ym>log(umbral)) # o solo 10
dftmp <- data.frame(y = ym[tmp], mun = muni[region[tmp]],
                    amb = factor(Xide$rururb[tmp], labels=c('Urbano', 'Rural')))
y_hist <- ggplot(dftmp, aes(x=y)) + 
  geom_histogram(aes(y = after_stat(density)), color='gray25', fill='gray90') +
  geom_density(linewidth=1.2, color='gray50', linetype=2) + 
  labs(x='log - ICTPC', y = 'Densidad estimada') + 
  theme_minimal() + 
  theme(text = element_text(family = "serif"))
x11(); plot(y_hist)
ggsave("./Figuras/dist-log-ict.png", plot = y_hist,
       # device = grDevices::cairo_pdf,
       width = 15, height = 15, units = "cm")

y_hist2 <- ggplot(dftmp, aes(x=y)) + 
  geom_histogram(aes(y = after_stat(density)), color='gray80', fill=4) +
  geom_density(linewidth=1, color='gray40') + 
  labs(x='log - ICTPC', y = 'Densidad estimada') + 
  theme_minimal() + 
  theme(text = element_text(family = "serif")) + 
  facet_wrap(~ mun, ncol = 4, nrow = 4)
x11(); plot(y_hist2)
ggsave("dist-log-ict-alcaldias.png", plot = y_hist2,
       # device = grDevices::cairo_pdf,
       width = 15, height = 15, units = "cm")

y_hist3 <- ggplot(dftmp, aes(x=y)) + 
  geom_histogram(aes(y = after_stat(density)), color='gray25', fill='gray90', bins = 30) +
  geom_density(linewidth=1, color='gray50', linetype=2) + 
  labs(x='log - ICTPC', y = 'Densidad estimada') + 
  theme_minimal() + 
  theme(text = element_text(family = "serif")) + 
  facet_wrap(~ amb, ncol = 2, nrow = 1, scales = 'fixed')
x11(); plot(y_hist3)
ggsave("./Figuras/dist-log-ict-ambito.png", plot = y_hist3,
       # device = grDevices::cairo_pdf,
       width = 15, height = 15, units = "cm")



## Gini-Lorenz antes del analisis posterior----
datos <- list()
datos$master <- which(!is.na(ym))
dftmp <- data.frame(Observados=ym[datos$master],
                    Ámbito=Xide$rururb[datos$master],
                    Factor=Xide$factor[datos$master],
                    Alcaldia=muni[region[datos$master]])

lc <- Lc(exp(dftmp$Observados), dftmp$Factor)
lc_urb <- Lc(exp(dftmp$Observados[dftmp$Ámbito==0]), dftmp$Factor[dftmp$Ámbito==0])
lc_rur <- Lc(exp(dftmp$Observados[dftmp$Ámbito==1]), dftmp$Factor[dftmp$Ámbito==1])
round(Gini(exp(dftmp$Observados)),4)
round(gini.wtd(exp(dftmp$Observados), dftmp$Factor),4)
gini_decomp(exp(dftmp$Observados), dftmp$Ámbito, dftmp$Factor)
# $gini_group
# $gini_group$gini_group
# 0         1 
# 0.4804783 0.3116281 



x11();
plot(lc,main="Indice de Gini y curva de Lorenz para EdoMex",lwd=2)
lines(lc_urb$p, lc_urb$L,col="red",lty=2,lwd=2)
lines(lc_rur$p, lc_rur$L,col="blue",lty=2,lwd=2)
legend(0.05,0.9,legend=c("Rural (0.3856)","Urbano (0.3598)","General (0.3656)"),
       col=c("Blue","Red","Black"),lty=c(2,2,1),lwd=2,bty="n",cex=1.1)
grid()



library(ggplot2)
library(scales)
factores <- ggplot(Xide, aes(x = as.factor(rururb), y = factor, fill=as.factor(rururb))) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.25, outlier.shape = 2, outlier.size = 4) +
  geom_jitter(width = 0.08, alpha = 0.05, size = 1) +
  labs(x = "Ámbito", y = bquote(log[10]- ~factor~de~expansión), title = "") +
  theme_minimal() + 
  theme(legend.position = 'none', text = element_text(family = "serif", size=11)) + 
  scale_y_log10(labels = label_log(digits = 2)) +
  annotation_logticks(sides='l')
x11();plot(factores)
ggsave(filename='factores_por_ambito.pdf', plot=factores, width=20, height = 20, units = 'cm')

# Para los mapas ====

output <- tibble(
  Municipio=stringr::str_pad(2:17, width=3, side = 'left', pad = 0),
  Porcentaje_LPI  = dftmp$Porcentaje_LPI,
  Porcentaje_LPEI = dftmp$Porcentaje_LPEI,
  Estimacion_rho  = apply(fit_vb$draws('rho'), 2, mean)
)
write.csv(output, file='Mapas/LPE_LPEI_rho_vb.csv', row.names = F)





