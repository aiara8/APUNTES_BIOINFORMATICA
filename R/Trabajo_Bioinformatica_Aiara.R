##########################################################################################
#####################                                             ########################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R      ########################
#####################                                             ########################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente.
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

datos <- read.table("archivo-entrada.txt",
                    header=TRUE,
                    sep="\t",
                    check.names=FALSE)

head(datos)

tail(datos)

dim(datos)


### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original.
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control,
### la media de Condición y la DE de Condición. (1,5 pto)

control <- datos[,c("C1-TPM","C2-TPM","C3-TPM")]

condicion <- datos[,c("1W5-TPM","1W6-TPM","1W5-COV")]

datos$Media_Control <- apply(control,1,mean)

datos$SD_Control <- apply(control,1,sd)

datos$Media_Condicion <- apply(condicion,1,mean)

datos$SD_Condicion <- apply(condicion,1,sd)

head(datos)


### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condición.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

par(mfrow=c(1,2))

hist(datos$Media_Control,
     col="skyblue",
     main="Media Control",
     xlab="TPM")

hist(datos$Media_Condicion,
     col="salmon",
     main="Media Condicion",
     xlab="TPM")

par(mfrow=c(1,1))


### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control.
### Añade el título, colores y una leyenda en la esquina superior izquierda. (1,5 pto)

boxplot(log10(datos$Media_Control + 1),
        log10(datos$Media_Condicion + 1),
        names=c("Control","Condicion"),
        col=c("skyblue","salmon"),
        main="Comparacion Control vs Condicion",
        ylab="log10(TPM+1)")

legend("topleft",
       legend=c("Control","Condicion"),
       fill=c("skyblue","salmon"))


### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente:
### Realiza un test de comparativa de medias para cada gen entre condición y control.
### Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original. (2 pto)

shapiro_control <- shapiro.test(datos$Media_Control)

shapiro_condicion <- shapiro.test(datos$Media_Condicion)

shapiro_control

shapiro_condicion

pvalues <- c()

for(i in 1:nrow(datos))
{
  test <- t.test(
    as.numeric(control[i,]),
    as.numeric(condicion[i,])
  )

  pvalues[i] <- test$p.value
}

datos$pvalue <- pvalues

head(datos)


### 6.- De los resultados del apartado 5, saca (mediante código) un listado
### de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué.
### En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

DEGs <- subset(datos,pvalue < 0.05)

DEGs

nrow(DEGs)

sobreexpresados <- subset(
  DEGs,
  Media_Condicion > Media_Control
)

reprimidos <- subset(
  DEGs,
  Media_Condicion < Media_Control
)

cat("Genes diferencialmente expresados:",
    nrow(DEGs), "\n")

cat("Sobreexpresados:",
    nrow(sobreexpresados), "\n")

cat("Reprimidos:",
    nrow(reprimidos), "\n")


### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

library(pheatmap)

DEGs_ordenados <- DEGs[order(DEGs$pvalue), ]

top50 <- DEGs_ordenados[1:50, ]

heatmap_data <- as.matrix(
top50[,c("C1-TPM",
          "C2-TPM",
          "C3-TPM",
          "1W5-TPM",
          "1W6-TPM",
          "1W5-COV")]
)

rownames(heatmap_data) <- top50[,1]

pheatmap(
 heatmap_data,
 scale="row",
 color=colorRampPalette(
   c("blue","white","red")
 )(100),
 main="Top 50 genes diferencialmente expresados"
)