## Cargamos librerías
library(knitr)
library(ggplot2)
library(caret)
library(class)
library(gmodels)
library(kableExtra)
library(neuralnet)
library(kernlab)
library(rpart)
library(rpart.plot)
library(randomForest)
library(rattle)
library(RColorBrewer)


# Preparamos directorios
workingDir <- getwd()
dataDir <- file.path(workingDir,"../Data/")
dir(path = dataDir)

SGA10 <- readRDS(paste0(dataDir, "SGA10_BremKruglyak2005_v1.rds"))
SGA16 <- readRDS(paste0(dataDir, "SGA16_BremKruglyak2005_v1.rds"))
allpaired <- readRDS(paste0(dataDir, "allpairsfeta.rds"))

# Visualizamos un análisis descriptivo de los datasets
dim(SGA10)
head(SGA10)
summary(SGA10)
str(SGA10)

dim(SGA16)
head(SGA16)
str(SGA16)

dim(allpaired)
head(allpaired)
tail(allpaired)
str(allpaired)

# Eliminamos las variables identificativas de los genes individuales
allpaired$V1 <- NULL
allpaired$V2 <- NULL

# Convertimos a tipo numérico todas las variables del dataset *allpaireid*
allpaired[] <- lapply(allpaired, as.numeric)
str(allpaired)

# Generamos dos nuevos data frames filtrando sólo los pares existentes en los SGAs
allpaired_10 <- allpaired[row.names(allpaired) %in% row.names(SGA10), ]
dim(allpaired_10)
# Vemos que en lugar de 227 registros hay 190, debido a que algunos genes han cambiado de identificador

allpaired_16 <- allpaired[row.names(allpaired) %in% row.names(SGA16), ]
dim(allpaired_16)
# Vemos que en lugar de 757 registros, hay 635, debido a que algunos genes han cambiado de identificador


# Generamos un nuevo dataset, uniendo los campos de SGA y allpaired, para cada uno de los experimientos
SGA10_full <- merge(SGA10, allpaired_10, by=0)
rownames(SGA10_full) <- SGA10_full$Row.names
SGA10_full$Row.names <- NULL
dim(SGA10_full)
head(SGA10_full)

SGA16_full <- merge(SGA16, allpaired_16, by=0)
rownames(SGA16_full) <- SGA16_full$Row.names
SGA16_full$Row.names <- NULL
dim(SGA16_full)
head(SGA16_full)


# Eliminamos las variables "descriptivas" de los genes
SGA10_full$GENE1 <- NULL
SGA10_full$GENE2 <- NULL

SGA16_full$GENE1 <- NULL
SGA16_full$GENE2 <- NULL

# Reproducimos el modelo de regresión lineal del TFG Aina Rill de la bibliografía, para verificar que partimos de la base correcta
# Modelo 1 (Tabla 1)
mod1 <- lm(SGAsco ~. -NPC -fcMF , SGA10_full)
summary(mod1)

(RMSE.lm1 <- sqrt(mean(mod1$residuals^2)))

# Modelo 2 (Tabla 2)
mod2 <- lm(SGAsco ~ Ohnology + Complexes + SameProtein +
                SameFunction + fcBP + fcCC + NPC, SGA10_full)
summary(mod2)

(RMSE.lm2 <- sqrt(mean(mod2$residuals^2)))

# Los resultados, si bien no son idénticos, son muy parecidos.

mod1.2 <- lm(SGAsco ~. -PCC -fcMF, SGA10_full)
summary(mod1.2)

(RMSE.lm1.2 <- sqrt(mean(mod1.2$residuals^2)))

mod2.2 <- lm(SGAsco ~ . -PCC -Homology -PhysicalInt -CommonReg
             -Colocalization - Phenotype -fcMF,SGA10_full)
summary(mod2.2)
(RMSE.lm2.2 <- sqrt(mean(mod2.2$residuals^2)))

mod2_16 <- lm(SGAsco ~ Ohnology + Complexes + SameProtein +
                SameFunction + fcBP + fcCC + NPC, SGA16_full)
summary(mod2_16)

(RMSE.lm2_16 <- sqrt(mean(mod2_16$residuals^2)))

mod2.2_16 <- lm(SGAsco ~ . -PCC -Homology -PhysicalInt -CommonReg
             -Colocalization - Phenotype -fcMF,SGA16_full)
summary(mod2.2_16)
(RMSE.lm2.2_16 <- sqrt(mean(mod2.2_16$residuals^2)))

# Algoritmos Machine Learning
## Prepararemos los datos tanto para el dataset SGA10 como SGA16, pero inicialmente realizaremos los cálculos sólo sobre los datos de SGA10, para simplificar la comparación, ya que es sobre el que tenemos la referencia de regresión lineal.

### Selección de datos para *training* y *test*
#### Dividimos los dataset en dos, para entrenamiento y validación, así como el vector con los valores a predecir (*labels*)
bound_10 <- floor(nrow(SGA10_full)*0.67)
bound_16 <- floor(nrow(SGA16_full)*0.67)

seed <- c(12345)

set.seed(seed)
row_train_10 <- sample(seq_len(nrow(SGA10_full)), size = bound_10)

set.seed(seed)
row_train_16 <- sample(seq_len(nrow(SGA16_full)), size = bound_16)

SGA10_train <- SGA10_full[row_train_10, ]
SGA10_test <- SGA10_full[-row_train_10, ]

SGA16_train <- SGA16_full[row_train_16, ]
SGA16_test <- SGA16_full[-row_train_16, ]


# Algoritmo k-NN
## Empezaremos por el algoritmo k-NN (Nearest Neighbors), que aunque se usa habitualmente en modelos de clasificación, también puede ofrecer buenos resulados de regresión.
## Transformación de los datos
#### No es necesario, al disponer ya de los datos normalizados y en datasets de *training* y *test* de pasos anteriores.

## Entrenamiento del modelo

# Para el coeficiente PCC:
set.seed(seed)
knn_model_PCC_1 <- train(SGAsco~. -NPC - fcMF, data=SGA10_train, method="knn")

plot(knn_model_PCC_1)

k_PCC_1 <- rownames(knn_model_PCC_1$bestTune)
knn_model_PCC_1$results[k_PCC_1,]

# Realizar la predicción
pred_PCC_1 <- knn_model_PCC_1 %>% predict(SGA10_test)

## Evaluar el rendimiento del modelo
## Generamos una tabla para poder visualizar un resumen de las diferentes combinaciones:
tabla_knn <- data.frame(Option = c("Default", "Param"), PCC = NA, NPC = NA)

(tabla_knn[1,2] <- RMSE(pred_PCC_1, SGA10_test$SGAsco))

# Para la variable NPC:
set.seed(seed)
knn_model_NPC_1 <- train(SGAsco~. -PCC, data=SGA10_train, method="knn")

plot(knn_model_NPC_1)
k_NPC_1 <- rownames(knn_model_NPC_1$bestTune)
knn_model_NPC_1$results[k_NPC_1,]

# Realizar la predicción
pred_NPC_1 <- knn_model_NPC_1 %>% predict(SGA10_test)

## Evaluar el rendimiento del modelo
(tabla_knn[1,3] <- RMSE(pred_NPC_1, SGA10_test$SGAsco))

# Mejorar el rendimiento del modelo
# Podemos probar cambiando alguno de los parámetros por defecto:
set.seed(seed)
knn_model_PCC_2 <- train(SGAsco~. -NPC, data=SGA10_train, method="knn",
                       trControl = trainControl("cv", number = 10),
                       preProcess = c("center","scale"),
                       tuneLength = 15)
plot(knn_model_PCC_2)
k_PCC_2 <- rownames(knn_model_PCC_2$bestTune)
knn_model_PCC_2$results[k_PCC_2,]
pred_PCC_2 <- knn_model_PCC_2 %>% predict(SGA10_test)
(tabla_knn[2,2] <- RMSE(pred_PCC_2, SGA10_test$SGAsco))

set.seed(seed)
knn_model_NPC_2 <- train(SGAsco~. -PCC, data=SGA10_train, method="knn",
                       trControl = trainControl("cv", number = 10),
                       preProcess = c("center","scale"),
                       tuneLength = 15)
plot(knn_model_NPC_2)
k_NPC_2 <- rownames(knn_model_NPC_2$bestTune)
knn_model_NPC_2$results[k_NPC_2,]
pred_NPC_2 <- knn_model_NPC_2 %>% predict(SGA10_test)
(tabla_knn[2,3] <- RMSE(pred_NPC_2, SGA10_test$SGAsco))

## Resumen resultados algoritmo SVM:
tabla_knn
# kable(tabla_ann, col.names=c("Option", "PCC", "NPC"),
#      align= c("l","c"))


# Algoritmo Artificial Neural Network
## Entrenamiento del modelo
set.seed(seed)
ann_model_1_PCC <- neuralnet(SGAsco ~. -NPC, data=SGA10_train)
plot(ann_model_1_PCC, rep="best")

set.seed(seed)
ann_model_1_NPC <- neuralnet(SGAsco ~. -PCC, data=SGA10_train)
plot(ann_model_1_NPC, rep="best")

## Evaluar el rendimiento del modelo
tabla_ann <- data.frame(hidden = c(1,2,3), PCC = NA, NPC = NA)

ann_results_1_PCC <- compute(ann_model_1_PCC, SGA10_test)
ann_pred_1_PCC <- ann_results_1_PCC$net.result
(tabla_ann[1,2] <- RMSE(ann_pred_1_PCC, SGA10_test$SGAsco))

ann_results_1_NPC <- compute(ann_model_1_NPC, SGA10_test)
ann_pred_1_NPC <- ann_results_1_NPC$net.result
(tabla_ann[1,3] <- RMSE(ann_pred_1_NPC, SGA10_test$SGAsco))


## Mejorar el rendimiento del modelo

#### Podemos intentar mejorar el rendimiento incrementando el número de *hidden nodes*.
#2
set.seed(seed)
ann_model_2_PCC <- neuralnet(SGAsco ~. -NPC, data=SGA10_train, hidden = 2)
plot(ann_model_2_PCC, rep="best")
ann_results_2_PCC <- compute(ann_model_2_PCC, SGA10_test)
ann_pred_2_PCC <- ann_results_2_PCC$net.result
(tabla_ann[2,2] <- RMSE(ann_pred_2_PCC, SGA10_test$SGAsco))

set.seed(seed)
ann_model_2_NPC <- neuralnet(SGAsco ~. -PCC, data=SGA10_train, hidden = 2)
plot(ann_model_2_NPC, rep="best")
ann_results_2_NPC <- compute(ann_model_2_NPC, SGA10_test)
ann_pred_2_NPC <- ann_results_2_NPC$net.result
(tabla_ann[2,3] <- RMSE(ann_pred_2_NPC, SGA10_test$SGAsco))

# 3
set.seed(seed)
ann_model_3_PCC <- neuralnet(SGAsco ~. -NPC, data=SGA10_train, hidden = 3)
plot(ann_model_3_PCC, rep="best")
ann_results_3_PCC <- compute(ann_model_3_PCC, SGA10_test)
ann_pred_3_PCC <- ann_results_3_PCC$net.result
(tabla_ann[3,2] <- RMSE(ann_pred_3_PCC, SGA10_test$SGAsco))

set.seed(seed)
ann_model_3_NPC <- neuralnet(SGAsco ~. -PCC, data=SGA10_train, hidden = 3)
plot(ann_model_3_NPC, rep="best")
ann_results_3_NPC <- compute(ann_model_3_NPC, SGA10_test)
ann_pred_3_NPC <- ann_results_3_NPC$net.result
(tabla_ann[3,3] <- RMSE(ann_pred_3_NPC, SGA10_test$SGAsco))

## Resumen resultados algoritmo SVM:
tabla_ann
# kable(tabla_ann, col.names=c("Option", "PCC", "NPC"),
#      align= c("l","c"))


# Algoritmo Support Vector Machine

## Entrenamiento del modelo
#### Empezamos entrenando el modelo con el *kernel* lineal (`vanilladot`)
set.seed(seed)
svm_model_PCC_1 <- ksvm(SGAsco ~. -NPC, data=SGA10_train, kernel = "vanilladot")

set.seed(seed)
svm_model_NPC_1 <- ksvm(SGAsco ~. -PCC, data=SGA10_train, kernel = "vanilladot")

## Evaluar el rendimiento del modelo
#Generamos la tabla en blanco para almacenar los diferentes resultados, como hemos hecho en el resto de algoritmos.
tabla_svm <- data.frame(kernel = c("vanilladot", "rbfdot", "laplacedot" ), PCC = NA, NPC = NA)

svm_pred_PCC_1 <- predict(svm_model_PCC_1, SGA10_test)
(tabla_svm[1,2] <- RMSE(svm_pred_PCC_1, SGA10_test$SGAsco))

svm_pred_NPC_1 <- predict(svm_model_NPC_1, SGA10_test)
(tabla_svm[1,3] <- RMSE(svm_pred_NPC_1, SGA10_test$SGAsco))

## Mejorar el rendimiento del modelo

#### Podemos intentar mejorar el rendimiento cambiando el tipo de *kernel* que utiliza el modelo.
#### Con el *kernel* `rbfdot` (*Radial Basis kernel "Gaussian"*)

set.seed(seed)
svm_model_PCC_2 <- ksvm(SGAsco ~. -NPC, data=SGA10_train, kernel = "rbfdot")
svm_pred_PCC_2 <- predict(svm_model_PCC_2, SGA10_test)
(tabla_svm[2,2] <- RMSE(svm_pred_PCC_2, SGA10_test$SGAsco))

set.seed(seed)
svm_model_NPC_2 <- ksvm(SGAsco ~. -PCC, data=SGA10_train, kernel = "rbfdot")
svm_pred_NPC_2 <- predict(svm_model_NPC_2, SGA10_test)
(tabla_svm[2,3] <- RMSE(svm_pred_NPC_2, SGA10_test$SGAsco))


#### Con el *kernel* `laplacedot` (*Laplacian kernel*)

set.seed(seed)
svm_model_PCC_3 <- ksvm(SGAsco~. -NPC, data=SGA10_train, kernel = "laplacedot")
svm_pred_PCC_3 <- predict(svm_model_PCC_3, SGA10_test)
(tabla_svm[3,2] <- RMSE(svm_pred_PCC_3, SGA10_test$SGAsco))

set.seed(seed)
svm_model_NPC_3 <- ksvm(SGAsco ~. -PCC, data=SGA10_train, kernel = "laplacedot")
svm_pred_NPC_3 <- predict(svm_model_NPC_3, SGA10_test)
(tabla_svm[3,3] <- RMSE(svm_pred_NPC_3, SGA10_test$SGAsco))

## Resumen resultados algoritmo SVM:
tabla_svm
# kable(tabla_svm, col.names=c("Option", "PCC", "NPC"),
#      align= c("l","c"))


# Algoritmo Árbol de Decisión
## Entrenamiento del Modelo

set.seed(seed)
ad_model_PCC_1 <- rpart(SGAsco ~. -NPC, data=SGA10_train,)
fancyRpartPlot(ad_model_PCC_1, caption = NULL)

set.seed(seed)
ad_model_NPC_1 <- rpart(SGAsco ~. -PCC, data=SGA10_train)
fancyRpartPlot(ad_model_NPC_1, caption = NULL)


## Evaluar el rendimiento del modelo
tabla_ad <- data.frame(Option = c("Default", "Param"), PCC = NA, NPC = NA)

ad_pred_PCC_1 <- predict(ad_model_PCC_1, SGA10_test)
(tabla_ad[1,2] <- RMSE(ad_pred_PCC_1, SGA10_test$SGAsco))

ad_pred_NPC_1 <- predict(ad_model_NPC_1, SGA10_test)
(tabla_ad[1,3] <- RMSE(ad_pred_NPC_1, SGA10_test$SGAsco))

## Mejorar el rendimiento del modelo
## Podemos probar a cambiar algún parámetro por defecto para comprobar si mejora el rendimiento del modelo

set.seed(seed)
ad_model_PCC_2 <- rpart(SGAsco ~. -NPC, data=SGA10_train,
                      minsplit = 2, minbucket = 1)
fancyRpartPlot(ad_model_PCC_2, caption = NULL)
ad_pred_PCC_2 <- predict(ad_model_PCC_2, SGA10_test)
(tabla_ad[2,2] <- RMSE(ad_pred_PCC_2, SGA10_test$SGAsco))


set.seed(seed)
ad_model_NPC_2 <- rpart(SGAsco ~. -PCC, data=SGA10_train,
                        minsplit = 2, minbucket = 1)
fancyRpartPlot(ad_model_NPC_2, caption = NULL)
ad_pred_NPC_2 <- predict(ad_model_NPC_2, SGA10_test)
(tabla_ad[2,3] <- RMSE(ad_pred_NPC_2, SGA10_test$SGAsco))

## Resumen resultados algoritmo Arbol de decisión:
tabla_ad
# kable(tabla_ad, col.names=c("Option", "PCC", "NPC"),
#      align= c("l","c"))


# Algoritmo Random Forest
## Entrenamiento del Modelo
set.seed(seed)
rf_model_PCC_1<- randomForest(SGAsco ~. -NPC, data=SGA10_train)

set.seed(seed)
rf_model_NPC_1<- randomForest(SGAsco ~. -PCC, data=SGA10_train)

#### Mostramos un gráfico del modelo
plot(rf_model_PCC_1)
plot(rf_model_NPC_1)

#### Observamos que el comportamiento es bastante estable a partir de 100 árboles aproximádamente.

#### Visualizamos la importancia de las variables
varImpPlot(rf_model_PCC_1)
varImpPlot(rf_model_NPC_1)

#### Observamos por ejemplo que la variable *NPC* gana importancia respecto *PCC*.

## Evaluar el rendimiento del modelo
tabla_rf <- data.frame(Option = c("Default", "Param"), PCC = NA, NPC = NA)

rf_pred_PCC_1 <- predict(rf_model_PCC_1, SGA10_test)
(tabla_rf[1,2] <- RMSE(rf_pred_PCC_1, SGA10_test$SGAsco))

rf_pred_NPC_1 <- predict(rf_model_NPC_1, SGA10_test)
(tabla_rf[1,3] <- RMSE(rf_pred_NPC_1, SGA10_test$SGAsco))

## Mejorar el rendimiento del modelo

#### Podemos intentar mejorar el rendimiento incrementando el número de árboles.
set.seed(seed)
rf_model_PCC_1000 <- randomForest(SGAsco ~. -NPC, data=SGA10_train, ntree=1000)
rf_pred_PCC_1000 <- predict(rf_model_PCC_1000, SGA10_test)
(tabla_rf[2,2] <- RMSE(rf_pred_PCC_1000, SGA10_test$SGAsco))

set.seed(seed)
rf_model_NPC_1000 <- randomForest(SGAsco ~. -PCC, data=SGA10_train, ntree=1000)
rf_pred_NPC_1000 <- predict(rf_model_NPC_1000, SGA10_test)
(tabla_rf[2,3] <- RMSE(rf_pred_NPC_1000, SGA10_test$SGAsco))

### Vemos que el modelo no ha mejorado aumentando el número de árboles del bosque a 1.000, como se intuía en el gráfico anterior.

## Resumen resultados algoritmo Arbol de decisión:
tabla_rf
# kable(tabla_rf, col.names=c("Option", "PCC", "NPC"),
#      align= c("l","c"))

### Como resumen de todos los modelos
RMSE.lm1
tabla_knn
tabla_ann
tabla_svm
tabla_ad
tabla_rf


# Vemos que ninguno mejora el valor del modelo de regresión lineal.
# El que ha obtenido el mejor valor es el algoritmo Arbol de Decisión parametrizado, y con la variable NPC
# Si realizamos los cálculos que con el dataset SGA16:

set.seed(seed)
ad_model_NPC_16 <- rpart(SGAsco ~. -PCC, data=SGA16_train,
                        minsplit = 2, minbucket = 1)
fancyRpartPlot(ad_model_NPC_16, caption = NULL)
ad_pred_NPC_16 <- predict(ad_model_NPC_16, SGA16_test)
RMSE(ad_pred_NPC_16, SGA16_test$SGAsco)

# Vemos que aunque incrementemos el número de observaciones, la predicción no mejora

