# =============== Analisis a nivel hogares del Ing. cor. total per capita =================

# manipulación y visualización
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
# ajuste del modelo con Stan
library(cmdstanr)
# para calcular Gini
library(ineq)
library(dineq)

# (0) Datos para el modelo ====
rm(list=ls())
gc()
ruta <- r'(C:\Maestria\Articulo-GitHub)'
setwd(ruta)
load('./codigo-r-stan-datos/Xyregion-CDMX.Rdata')
source('./codigo-r-stan-datos/funciones-final.R')

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

dim(Xtmp) # 110 covariables

## (0.4) Componentes principales ====
# Eigen <- eigen(cov(Xcon))
# cumsum(Eigen$values) / ncol(Xcon) # 26: 95%
# Xpc <- Xcon %*% Eigen$vectors[, 1:26]
# Xtmp <- cbind(Xpc, Xbin)

# (1) Modelo para predicción ====
# Un intercepto puede suplir la necesidad de una a priori jerárquica en varios interceptos
file <- file.path('codigo-r-stan-datos', 'LogSN-rho01-SSVS-oneInt.stan')
mod <- cmdstan_model(file, compile=T) # compilar genera el .exe

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
  eval_elbo = 100,  # mas muestras, mas lento
  draws = 1000,
  adapt_iter=1000,
  adapt_engaged=F, # buen resultado con eta=0.1
  eta=0.1, # (escala al) tamaño de paso
  iter=75000,
  tol_rel_obj=1e-6,
  algorith='meanfield', # o fullrank, mas lento
)
tproc_vb <- proc.time()-t0 # Tiempo de proceso

## (1.4) Correr modelo con HMC ====
## es mas lento este ajuste, no correr
## en cambio, puede cargarse un ajuste previo con load(...)
# t0 <- proc.time()
# fit_hmc <- mod$sample(
#   output_dir = r'(C:\Maestria\Articulo-GitHub)',
#   data = dat,
#   seed = 123,
#   init = list(init_list),
#   chains = 1,
#   thin = 2,
#   iter_warmup = 10000,
#   iter_sampling = 2000,
#   refresh = 1000
# )
# tproc_hmc <- proc.time()-t0 # Tiempo de proceso
load(r'(codigo-r-stan-datos\fit_hmc-oneInt_burnin_13k_iter2k_thin2_out8-15_noCPV.Rdata)')
tproc_hmc <- readLines(r'(codigo-r-stan-datos\Tiempo-hmc-oneInt.txt)')
tproc_hmc <- as.numeric(str_split(tproc_hmc, '-',simplify = T)[1:3, 1])

## (1.4) Precisión y tiempo (predicción) ====
fitted_vb0 <- apply(exp(fit_vb$draws('y_mis')), 2, mean)
fitted_hmc0 <- apply(drop(exp(fit_hmc$draws('y_mis'))), 2, mean)
m0 <- c(metrics(fitted_vb0, exp(datos$y_test)), tproc_vb[3])
m1 <- c(metrics(fitted_hmc0, exp(datos$y_test)), tproc_hmc[3])
## Cuadro 4. Métricas del pronóstico
round(m0, 4)
round(m1, 4)

## (1.5) Grafico observados vs ajustados =====
## se grafican en escala log. para mejor visualizacion
fitted_vb  <- apply(fit_vb$draws('y_mis'), 2, mean)
fitted_hmc <- apply(drop(fit_hmc$draws('y_mis')), 2, mean)
cor(fitted_vb, datos$y_test)
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
# Figura 4, (a) y (b). log-observados vs log-ajustados
# ggsave('./Figuras/LogSN-CDMX-vb.png', plot=g0, width=15, heigh=15, units='cm')
# ggsave('./Figuras/LogSN-CDMX-hmc.png', plot=g1, width=15, heigh=15, units='cm')

# (2) Modelo para estimación ====
file <- file.path('codigo-r-stan-datos', 'LogSN-rho01-SSVS-genq-oneInt.stan')
mod <- cmdstan_model(file, compile=T) 

## (2.1) Seleccionar datos ====
datos <- get_obs(censo=T)
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
  c=3000, # sqrt(10*300**2)
  epsilon=0.01
)

## (2.2) Valores iniciales ====
## usando lm para beta
tmp <- data.frame(y_obs=datos$y_obs, datos$X_obs)
r0 <- lm(y_obs ~ ., data=tmp)
## usando selm de sn para rho, sigma y mu
r1 <- coef(sn::selm(y_obs ~ 1 , data=tmp), param.type='DP')
mu0 <- r1[1]
sigma0 <- r1[2]
rho0 <- lambda2rho(r1[3])

set.seed(2)
init_list <- list(
  b = get_coefs(r0, F),
  pr=rep(0.5, ncol(datos$X_mis)),
  sigma=sigma0,
  mu=mu0,
  rho=rep(rho0, times=K)
)

## (2.3) Correr modelo con VB ====
t0 <- proc.time()
fit_vb <- mod$variational(
  data = dat,
  grad_samples = 1, # mas muestras, mas lento
  eval_elbo = 100,  # mas muestras, mas lento
  seed = 123,
  init = list(init_list),
  draws = 1000,
  adapt_iter=5000,
  adapt_engaged=F, # buen resultado con eta=0.1
  eta=0.1, # (escala al) tamaño de paso
  iter=75000,
  tol_rel_obj=1e-6,
  algorith='meanfield', # o fullrank, mas lento
)
tproc_vb <- proc.time()-t0 # Tiempo de proceso

# (3) Estimación a posteriori ====

## (3.1) Estimación y selección de beta ====
## se mostrará el modelo más popular 
## y sus probabilidades marginales de aparición
beta_ <- fit_vb$draws("b")
Media <- apply(beta_, 2, mean)
Error <- apply(beta_, 2, sd)
Lo025 <- apply(beta_, 2, quantile, 0.025)
Up975 <- apply(beta_, 2, quantile, 0.975)
Frecuencia <- marginal_pi(fit_vb)
Parametro <- paste0("$\\beta_{", 1:ncol(dat$X_obs), "}$")
idx0 <- 1:ncol(dat$X_obs)
idx1 <- str_split(posterior_prob(fit_vb)[1, 2], '', simplify = T) %>%
  unlist() %>% as.numeric() %>% as.logical()

load('./codigo-r-stan-datos/TablaCovariables.Rdata') # tabla de covariables procesada
beta_out <- tibble(
  Parametro,
  Freq=Frecuencia,
  Media,
  ErrorE=Error,
  Lo025=Lo025,
  Up975=Up975,
  Descripcion=tabla_cov$Descripción[idx0]
  )

beta_out <-beta_out[idx1, ]

## Cuadro 8. Estimaciones de beta y lista de covariables
print(beta_out, n=Inf)

## (3.2) Grafico de mosaico ====
## no incluido

tmp <- cbind(posterior_prob(fit_vb)[1:5, ], Metodo='VB')
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
  scale_fill_manual(labels=c('0'='No', '1'='Sí'), values=c('0'='gray80', '1'='gray50'))+
  geom_tile(color = "white", linewidth = 0.4) +
  facet_wrap(~Método, scales = 'free_y') +
  scale_x_discrete(labels=parse(text = paste0("beta[", 1:ncol(Xtmp), "]"))) +
  theme_minimal() +
  labs(fill='Inclusión', y='Prob. a posteriori', x='') +
  theme(legend.position = 'top', text=element_text(family='serif'),
        plot.margin = margin(0, 0, 0, 0),  axis.ticks.length = unit(0, "pt"))
x11(); print(g0)
# no hay mucha diferencia en la eleccion de modelos

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

Parametro <- c(
  paste0("$\\rho_{", 1:16, "}$"),
  paste0("$\\sigma$"),
  paste0("$\\mu$")
)

param_out <- tibble(
  Parametro=Parametro,
  Media=c(Media_rho, Media_sigma, Media_mu),
  ErrorE=c(Error_rho, Error_sigma, Error_mu),
  Lo025=c(Lo025_rho, Lo025_sigma, Lo025_mu),
  Lo975=c(Up975_rho, Up975_sigma, Up975_mu),
)

## Cuadro 5. Estimaciones de rho, sigma y mu
print(param_out, n=Inf)

## (3.4) Gini-Lorenz a posteriori ====
# 0: ámbito urbano, 1: ámbito rural
out <- which(is.na(ym))
fitted_vb <- apply(exp(fit_vb$draws('y_all')), 2, mean)

dftmp <- data.frame(
  Ajustados=fitted_vb[out],
  Contexto=Xide$rururb[out],
  Factor=Xide$factor[out]
)

lc <- Lc(dftmp$Ajustados, dftmp$Factor)
lc_urb <- Lc(dftmp$Ajustados[dftmp$Contexto==0], dftmp$Factor[dftmp$Contexto==0])
lc_rur <- Lc(dftmp$Ajustados[dftmp$Contexto==1], dftmp$Factor[dftmp$Contexto==1])

# Gini ponderado
round(gini.wtd((dftmp$Ajustados), dftmp$Factor),6)
# Gini ponderado por ámbito
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
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = c('Urbano'='gray75', 'Rural'='gray50', 'Ambos'='gray0'))+
  coord_fixed(ratio = 1) +
  geom_abline(intercept = 0, slope=1, linewidth=1.2, color='gray75', lty=1) + 
  theme_bw() + 
  theme(text=element_text(family='serif', size=11), legend.position = 'top')

# Figura 5. Curva de Lorenz
x11(); plot(LG_plot)
# ggsave(filename='./Figuras/Lorenz-CDMX.png', plot=LG_plot, width=15, heigh=15, units='cm')

## (3.5) Porcentaje de personas con las LPI y LPEI ====
## no se duplica la info: ambas fuentes reproducen la ciudad (casi) completa
indices <- which(!is.na(ym))
sum((Xide$factor * Xide$tam_hog)[indices])  # provenientes de la ENIGH
sum((Xide$factor * Xide$tam_hog)[-indices]) # provenientes del Censo 2020

LPI_urb <- 4564.97
LPI_rur <- 3296.92
LPEI_urb <- 2354.65
LPEI_rur <- 1800.55

out <- which(is.na(ym))
fitted_vb <- apply(exp(fit_vb$draws('y_all')), 2, mean)
dftmp2 <- tibble(
  Ajustados=fitted_vb[out],
  Ámbito=Xide$rururb[out],
  Factor=Xide$factor[out],
  Tamaño=Xide$tam_hog[out],
  Alcaldia=muni[region][out]
)

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

## Cuadro 6. Porcentajes bajo LPI y LPEI
print(dftmp2[, -c(2:4)], n=Inf)

## (3.6) Cuantificar sesgo ====
## solo puede calcularse entre valores observados y ajustados
idx <- which(!is.na(ym))
fitted_vb <- apply(exp(fit_vb$draws('y_all')), 2, mean)
fitted_vb[-idx] <- NA

dftmp3 <- tibble(
  ICTPC_obs = exp(ym[idx]),
  ICTPC_fit=fitted_vb[idx],
  Alcaldia=muni[region][idx])

dftmp3 <- dftmp3 %>%
  group_by(Alcaldia) %>%
  summarise(mean_ICTPC_obs = mean(ICTPC_obs, na.rm=T), # cambiar por median, esta curioso
            mean_ICTPC_fit = mean(ICTPC_fit, na.rm=T),
            diff_ICTPC = abs(mean_ICTPC_obs-mean_ICTPC_fit))

# ordenar de acuerdo a muni
dftmp3 <- dftmp3[match(muni, dftmp3$Alcaldia), ]

## Cuadro 7. Sesgo agrupado por alcaldía
print(dftmp3, n=Inf)
mean(dftmp3$diff_ICTPC)

# (4) Descriptivos ====

## (4.1) Tablas de estadísticos ====

datos <- list()
datos$master <- which(!is.na(ym))
dftmp <- data.frame(Observados=ym[datos$master],
                    Ámbito=Xide$rururb[datos$master],
                    Factor=Xide$factor[datos$master],
                    Alcaldia=muni[region[datos$master]])

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

## Cuadro 2. Estadísticos por alcaldía y ámbito
print(table_out, n=Inf)

Res_ambito <- dftmp %>%
  group_by(Ámbito) %>%
  summarise(Minimo=min(exp(Observados)),
            Mediana = median(exp(Observados)),
            Media = mean(exp(Observados)),
            Maximo=max(exp(Observados)),
            `Desv. est.`=sd(exp(Observados)),
            Conteo=n(), 
            .groups = "drop")

## Cuadro 3. Estadisticos por ámbito
print(Res_ambito)


## (4.2) Gráfico de distribución del ingreso (histogramas) ----
umbral <- 400 # filtrar tres ingresos muy pequeños
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

y_hist3 <- ggplot(dftmp, aes(x=y)) + 
  geom_histogram(aes(y = after_stat(density)), color='gray25', fill='gray90', bins = 30) +
  geom_density(linewidth=1, color='gray50', linetype=2) + 
  labs(x='log - ICTPC', y = 'Densidad estimada') + 
  theme_minimal() + 
  theme(text = element_text(family = "serif")) + 
  facet_wrap(~ amb, ncol = 2, nrow = 1, scales = 'fixed')
x11(); plot(y_hist3)

# Figura 3 (a) y (b)
# ggsave("./Figuras/dist-log-ict.png", plot = y_hist, width = 15, height = 15, units = "cm")
# ggsave("./Figuras/dist-log-ict-ambito.png", plot = y_hist3, width = 15, height = 15, units = "cm")



# (5) Mapas ----
# Requiere instalar y cargar varias librerías, además
# asume que se ha ejectuado todo el código previo
## (5.1) Cargar librerías ----
library(sf) # st_read
library(ggspatial) # annotation_map_tile

# Genera la lista de alcaldias
draw_key_cust <- function(data, params, size) {
  data_text <- data
  data_text[c("family")] <- 'serif'
  data_text[c("label")] <- key_label2[names(pal)[match(data$colour, pal)]]
  data_text[c("fill")] <- NULL
  data_text$colour <- "black" # color dentro del cuadrito
  data_text$alpha <- 1
  data_text$size <- 11 / .pt
  grid::grobTree(
    draw_key_rect(data, list()),
    draw_key_text(data_text, list())
  )
}

# cargar el mapa base
carto_lightnl_url <- "https://a.basemaps.cartocdn.com/light_nolabels/${z}/${x}/${y}.png" # no labels
# carto_light_url <- "https://a.basemaps.cartocdn.com/light_all/${z}/${x}/${y}.png" # con labels
# cargar el shapefile
shp <- st_read(r'(codigo-r-stan-datos\poligonos_alcaldias_cdmx\poligonos_alcaldias_cdmx.shp)')
# empezar en 001 y terminar en 016
shp$CVE_MUN <- str_pad(as.numeric(shp$CVE_MUN)-1, width=3, pad='0')
# estableces colores
pal <- scales::grey_pal(start=0, end=0.11)(16)
names(pal) <- shp$NOMGEO
key_label <- shp$NOMGEO
names(key_label) <- shp$CVE_MUN
key_label2 <- shp$CVE_MUN
names(key_label2) <- shp$NOMGEO
# poner cada estimación a la alcaldia correcta
orden <- match(dftmp2$Alcaldia, shp$NOMGEO) # el elemento i (i=1, 2, ...) de #1 esta en la posición indicada de #2
# dftmp2 contiene los porcentajes bajo LPI y LPEI
shp$LPI_porcentaje[orden] <- dftmp2$Porcentaje_LPI
shp$LPEI_porcentaje[orden] <- dftmp2$Porcentaje_LPEI

# Mapa LPI
g1 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPI_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPI (%)', title = "", x='', y='') + 
  theme(
    text = element_text(family = "serif", size=11),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g1)

g2 <- ggplot() +
  annotation_map_tile(type = carto_lightnl_url, zoomin = 0, alpha = 0.8, progress = "none") +
  geom_sf(data = shp, aes(fill= LPEI_porcentaje)) +
  geom_sf_label(data=shp, aes(label = CVE_MUN, colour=NOMGEO, family='serif'), size = 4,  key_glyph = "cust") + 
  scale_colour_manual(name="Alcaldías", values = pal, limits=key_label) + 
  annotation_scale(location = "bl", width_hint = 0.4, text_family='serif') +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) + 
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "Greys")) + 
  scale_fill_gradient(low = "gray100", high = "gray50", limits=c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(fill ='Población bajo\nLPEI (%)', title = "", x='', y='') + 
  theme(
    text = element_text(family = "serif", size=11),
    legend.position = 'right',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
x11(); print(g2)

getwd()
# ggsave(plot=g1, filename='./Figuras/mapa-lpi.png', width=17.53, height = 17.53, units='cm')
# ggsave(plot=g2, filename='./Figuras/mapa-lpei.png', width=17.53, height = 17.53, units='cm')





