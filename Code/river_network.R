library(stpca)
library(SSN)
library(rgdal)
library(GWmodel)
library(geosphere)
library(ggplot2)
library(ggdendro)
library(gasper)
library(grid)
library(gridBase)
library(gridExtra)
library(dplyr)
library(maps)
library(fields)

shape_Geum <- readOGR("Data/ProcessedData/Geum.shp")
# shape_miho <- readOGR("Data/ProcessedData/GeumMiho.shp")
# shape_miho_testR42 <- readOGR("Data/ProcessedData/GeumMiho.shp")
source("Code/source.R")
source("Code/utils.R")

binaryIDs_obj <- get_binaryIDs_stream(mouth_node=614 , shape_Geum) # defined in source.R, shape file needs RCH_ID, binaryID
# river mouth : node 614, called Geumgang Gapmun (금강갑문)

binaryIDs_obj_miho <- binaryIDs_obj
binaryIDs_obj_miho$rid <- 1:113
binaryIDs_obj_miho$binaryID


# get adjacency 
adjacency <- get_adjacency_stream(binaryIDs_obj) # defined in source.R 

shreve_order <- compute_shreve(adjacency) # defined in source.R 


# ========================== #
#  data load and preprocess  #
# ========================== #
ssn_Geum <- importSSN_stream_past(shreve_obj = shreve_order, location="Full", multipleplaces=TRUE) 
colnames(ssn_Geum@obspoints@SSNPoints[[1]]@point.data)[which(colnames(ssn_Geum@obspoints@SSNPoints[[1]]@point.data) == "위도.Degree.")] = "lat"
colnames(ssn_Geum@obspoints@SSNPoints[[1]]@point.data)[which(colnames(ssn_Geum@obspoints@SSNPoints[[1]]@point.data) == "경도.Degree.")] = "lon"

spatial.wt_geum <- createWeightS_from_ssn(ssndata=ssn_Geum, afvcol="shreve", adjacency=adjacency)


## load Geum River data
# start from 2018
geum2018 <- read.csv("Data/ProcessedData/geum2018.csv", 
                     header=TRUE, fileEncoding = "euc-kr")

geum2018 <- geum2018[geum2018$분류번호 %in% miho_ids,]

geum2018 <- geum2018[,c(1,2,3,8:14)]
colnames(geum2018) <- c("id","loc","date","DO","BOD","COD","TN","TP","TOC","SS")

geum2018 <- na.omit(geum2018)
geum2018$TP <- as.numeric(geum2018$TP)
geum2018$year <- sapply(strsplit(geum2018$date, "/"), function(x){x[1]})
geum2018$month <- sapply(strsplit(geum2018$date, "/"), function(x){x[2]})


geum2018_month <- geum2018 %>% group_by(loc,id,year,month) %>% 
  summarise(DO = mean(DO),BOD=mean(BOD),COD=mean(COD),
            TN=mean(TN), TP=mean(TP),TOC=mean(TOC),SS=mean(SS)) %>%
  arrange(location,year,month)


geum2018_month %>%
  group_by(year, month) %>%
  filter(n_distinct(loc) == length(unique(geum2018$loc))) %>%
  ungroup()


geum2019 <- read.csv("Data/ProcessedData/geum2019.csv", 
                     header=TRUE, fileEncoding = "euc-kr")

geum2019 <- geum2019[geum2019$분류번호 %in% miho_ids,]

geum2019 <- geum2019[,c(1,2,3,8:14)]
colnames(geum2019) <- c("id","loc","date","DO","BOD","COD","TN","TP","TOC","SS")

geum2019$year <- sapply(strsplit(geum2019$date, "/"), function(x){x[1]})
geum2019$month <- sapply(strsplit(geum2019$date, "/"), function(x){x[2]})

geum2019 <- na.omit(geum2019)
geum2019$TP <- as.numeric(geum2019$TP)

geum2019_month <- geum2019 %>% group_by(loc,id,year,month) %>% 
  summarise(DO = mean(DO),BOD=mean(BOD),COD=mean(COD),
            TN=mean(TN), TP=mean(TP),TOC=mean(TOC),SS=mean(SS)) %>%
  arrange(location,year,month)


geum2019_month %>%
  group_by(year, month) %>%
  filter(n_distinct(loc) == length(unique(geum2019$loc))) %>%
  ungroup()



geum2020 <- read.csv("Data/ProcessedData/geum2020.csv", 
                     header=TRUE, fileEncoding = "euc-kr")

geum2020 <- geum2020[geum2020$분류번호 %in% miho_ids,]

geum2020 <- geum2020[,c(1,2,3,8:14)]
colnames(geum2020) <- c("id","loc","date","DO","BOD","COD","TN","TP","TOC","SS")

geum2020$year <- sapply(strsplit(geum2020$date, "/"), function(x){x[1]})
geum2020$month <- sapply(strsplit(geum2020$date, "/"), function(x){x[2]})

geum2020 <- na.omit(geum2020)

geum2020_month <- geum2020 %>% group_by(loc,id,year,month) %>% 
  summarise(DO = mean(DO),BOD=mean(BOD),COD=mean(COD),
            TN=mean(TN), TP=mean(TP),TOC=mean(TOC),SS=mean(SS)) %>%
  arrange(location,year,month)



geum2021 <- read.csv("Data/ProcessedData/geum2021.csv", 
                     header=TRUE, fileEncoding = "euc-kr")

geum2021 <- geum2021[geum2021$분류번호 %in% miho_ids,]

geum2021 <- geum2021[,c(1,2,3,8:14)]
colnames(geum2021) <- c("id","loc","date","DO","BOD","COD","TN","TP","TOC","SS")

geum2021$year <- sapply(strsplit(geum2021$date, "/"), function(x){x[1]})
geum2021$month <- sapply(strsplit(geum2021$date, "/"), function(x){x[2]})

geum2021 <- na.omit(geum2021)

geum2021$TP <- as.numeric(geum2021$TP)


geum2021_month <- geum2021 %>% group_by(loc,id,year,month) %>% 
  summarise(DO = mean(DO),BOD=mean(BOD),COD=mean(COD),
            TN=mean(TN), TP=mean(TP),TOC=mean(TOC),SS=mean(SS)) %>%
  arrange(location,year,month)


geum2022 <- read.csv("Data/ProcessedData/geum2022.csv", 
                     header=TRUE, fileEncoding = "euc-kr")

geum2022 <- geum2022[geum2022$분류번호 %in% miho_ids,]

geum2022 <- geum2022[,c(1,2,3,8:14)]
colnames(geum2022) <- c("id","loc","date","DO","BOD","COD","TN","TP","TOC","SS")

geum2022$TP <- as.numeric(geum2022$TP)
geum2022$year <- sapply(strsplit(geum2022$date, "/"), function(x){x[1]})
geum2022$month <- sapply(strsplit(geum2022$date, "/"), function(x){x[2]})

geum2022 <- na.omit(geum2022)

geum2022_month <- geum2022 %>% group_by(loc,id,year,month) %>% 
  summarise(DO = mean(DO),BOD=mean(BOD),COD=mean(COD),
            TN=mean(TN), TP=mean(TP),TOC=mean(TOC),SS=mean(SS)) %>%
  arrange(location,year,month)


geum2023 <- read.csv("Data/ProcessedData/geum2023.csv", 
                     header=TRUE, fileEncoding = "euc-kr")

geum2023 <- geum2023[geum2023$분류번호 %in% miho_ids,]

geum2023 <- geum2023[,c(1,2,3,8:14)]
colnames(geum2023) <- c("id","loc","date","DO","BOD","COD","TN","TP","TOC","SS")

geum2023$TP <- as.numeric(geum2023$TP)
geum2023$year <- sapply(strsplit(geum2023$date, "/"), function(x){x[1]})
geum2023$month <- sapply(strsplit(geum2023$date, "/"), function(x){x[2]})

geum2023 <- na.omit(geum2023)

geum2023_month <- geum2023 %>% group_by(loc,id,year,month) %>% 
  summarise(DO = mean(DO),BOD=mean(BOD),COD=mean(COD),
            TN=mean(TN), TP=mean(TP),TOC=mean(TOC),SS=mean(SS)) %>%
  arrange(location,year,month)


geum_month <- rbind(geum2018_month, geum2019_month,geum2020_month,geum2021_month,
                    geum2022_month,geum2023_month)


geum_month <- geum_month %>% mutate(id = factor(id, levels = unique(geum_month$id)[order(match(unique(geum_month$id), miho_ids))])) %>%
  arrange(year, month, id,loc) 

# unique(geum_month$loc)[order(match(unique(geum_month$loc), miho_names))]

# match(unique(geum_month$id), ssn_Geum@obspoints@SSNPoints[[1]]@point.data$분류번호)
spatial.wt_miho_month <- spatial.wt_geum[match(unique(geum_month$id), ssn_Geum@obspoints@SSNPoints[[1]]@point.data$분류번호),
                                         match(unique(geum_month$id), ssn_Geum@obspoints@SSNPoints[[1]]@point.data$분류번호)]

spatial.wt_miho.sym_month <- spatial.wt_miho_month + t(spatial.wt_miho_month)

spatial.wt_miho.sym[lower.tri(spatial.wt_miho.sym)] <- t(spatial.wt_miho.sym)[lower.tri(spatial.wt_miho.sym)]
diag(spatial.wt_miho_month) <- 0
diag(spatial.wt_miho.sym_month) <- 0


colnames(spatial.wt_miho_month) <- unique(geum_month$loc)
rownames(spatial.wt_miho_month) <- unique(geum_month$loc)
colnames(spatial.wt_miho.sym_month) <- unique(geum_month$loc)
rownames(spatial.wt_miho.sym_month) <- unique(geum_month$loc)


# undirected graph approach (not used for analysis)
mihoriver_month <- list()
# location
mihoriver_month$xy <- mihoriver$xy[c(1:23,25),]
rownames(mihoriver_month$xy) <- unique(geum_month$loc)
# weight matrix
mihoriver_month$A <- spatial.wt_miho.sym_month

# sparse weight matrix
sp.wmat <- c()
for(i in 1:(nrow(spatial.wt_miho.sym_month)-1)){
  for(j in (i+1):nrow(spatial.wt_miho.sym_month)){
    if(spatial.wt_miho.sym_month[i,j]!=0){
      sp.wmat <- rbind(sp.wmat, c(rownames(spatial.wt_miho.sym_month)[i],rownames(spatial.wt_miho.sym_month)[j], 
                                  spatial.wt_miho.sym_month[i,j]))
    }
  }
}
sp.wmat <- as.data.frame(sp.wmat)
sp.wmat$V3 <- as.numeric(sp.wmat$V3)
mihoriver_month$sA <- sp.wmat


plot_graph(mihoriver_month)

# directed graph approach
mihoriver_month_directed <- mihoriver_month

mihoriver_month_directed$A <- spatial.wt_miho_month

# sparse weight matrix
sp.wmat <- c()
for(i in 1:nrow(spatial.wt_miho_month)){
  for(j in 1:nrow(spatial.wt_miho_month)){
    if(spatial.wt_miho_month[i,j]!=0){
      sp.wmat <- rbind(sp.wmat, c(rownames(spatial.wt_miho_month)[i],rownames(spatial.wt_miho_month)[j], 
                                  spatial.wt_miho_month[i,j]))
    }
  }
}
sp.wmat <- as.data.frame(sp.wmat)
sp.wmat$V3 <- as.numeric(sp.wmat$V3)
mihoriver_month_directed$sA <- sp.wmat


plot_graph_custom_directed(mihoriver_month_directed, e.size=1, v.size=6,  value="Passengers", ratio=0.9,
                           min=-2.15, max=1.92, mg=c(4,4,4,4), main.title.size = 20, signal=FALSE, directed=TRUE)


# Laplacian and GSO
L.miho_month <- laplacian_mat(mihoriver_month$A) # laplacian matrix
val1 <- eigensort(L.miho_month)
evalues.miho_month <- val1$evalues
evectors.miho_month <- val1$evectors
# largest eigenvalue
lmax.miho_month <- max(evalues.miho_month)

N.miho_month <- nrow(L.miho_month)


# tmpsvd <- svd(t(spatial.wt_miho_month))
tmpsvd <- svd(spatial.wt_miho_month)

S.miho_month <- tmpsvd$u %*% diag(c(rep(1,N.miho_month-1), det(tmpsvd$u%*%t(tmpsvd$v)))) %*% t(tmpsvd$v)
val2 <- eigen(S.miho_month)
evalues.miho_month2 <- val2$values
evectors.miho_month2 <- val2$vectors

# order complex-valued eigenvalues
ordering_tmp <- order(Mod(1-evalues.miho_month2))
evalues.miho_month2 <- evalues.miho_month2[ordering_tmp]
evectors.miho_month2 <- evectors.miho_month2[,ordering_tmp]

# preprocessing
p.miho_month <- 6
R.miho_month <- 72

geum_month$cnt <- rep(1:R.miho_month,each=N.miho_month)
geum_month <- as.data.frame(geum_month[,-11])


# mean signal
meansig_geum <- as.data.frame(geum_month %>% group_by(id,loc) %>% 
                                summarise(meanDO = mean(DO),meanBOD = mean(BOD),meanCOD = mean(COD),
                                          meanTN = mean(TN),meanTP = mean(TP),meanTOC = mean(TOC)) %>%
                                mutate(id = factor(id, levels = unique(geum_month$id)[order(match(unique(geum_month$id), miho_ids))])) %>%
                                arrange(id,loc))

head(meansig_geum,20)


meansig_geum_log <- as.data.frame(geum_month %>% group_by(id,loc) %>% 
                                    summarise(meanlog_DO = mean(log(1+DO)),meanlog_BOD = mean(log(1+BOD)),meanlog_COD = mean(log(1+COD)),
                                              meanlog_TN = mean(log(1+TN)),meanlog_TP = mean(log(1+TP)),meanlog_TOC = mean(log(1+TOC))) %>%
                                    mutate(id = factor(id, levels = unique(geum_month$id)[order(match(unique(geum_month$id), miho_ids))])) %>%
                                    arrange(id,loc))

head(meansig_geum_log,24)


# generate X
X.miho_month <- array(0, c(p.miho_month,N.miho_month,R.miho_month))

for(i in 1:p.miho_month){
  for(r in 1:R.miho_month){
    data.miho <- geum_month[geum_month$cnt==r,(4+i)]
    # data.metro_getin_line2_timediv_wo_centering <- log(1+data.metro_getin_line2_timediv)
    data.miho <- log(1+data.miho)
    data.miho <- data.miho - meansig_geum_log[,2+i]
    X.miho_month[i,,r] <- data.miho
  }
}


# Fourier coefficients
X.tilde.miho_month <- array(0, c(p.miho_month,N.miho_month,R.miho_month)) # for undirected
X.tilde.miho_month2 <- array(0, c(p.miho_month,N.miho_month,R.miho_month)) # for directed

for(i in 1:p.miho_month){
  X.tilde.miho_month[i,,] <- Conj(t(evectors.miho_month)) %*% X.miho_month[i,,]
  X.tilde.miho_month2[i,,] <- Conj(t(evectors.miho_month2)) %*% X.miho_month[i,,]
}

# reordering
X.ut.miho_month <- array(0,c(N.miho_month,R.miho_month,p.miho_month)) # n x R x p
X.ut.miho_month2 <- array(0,c(N.miho_month,R.miho_month,p.miho_month)) # n x R x p

for(l in 1:N.miho_month){
  X.ut.miho_month[l,,] <- t(X.tilde.miho_month[,l,])
  X.ut.miho_month2[l,,] <- t(X.tilde.miho_month2[,l,])
}


k.optimal.miho_month <- 2 # inferred from the scree plot below


f.bar.miho_month <- array(0, c(N.miho_month,R.miho_month,k.optimal.miho_month)) # n x R x k
B.bar.miho_month <- array(0, c(N.miho_month,p.miho_month,k.optimal.miho_month)) # n x p x k
f.bar.miho_month2 <- array(0, c(N.miho_month,R.miho_month,k.optimal.miho_month)) # n x R x k
B.bar.miho_month2 <- array(0, c(N.miho_month,p.miho_month,k.optimal.miho_month)) # n x p x k

eigenprops <- c() # for scree plot to determine factor number
eigenprops2 <- c() # for scree plot to determine factor number
for(l in 1:N.miho_month){
  eigenres <- eigen(X.ut.miho_month[l,,] %*% Conj(t(X.ut.miho_month[l,,])))
  evectors <- eigenres$vectors
  f.bar.miho_month[l,,] <- sqrt(R.miho_month)*evectors[,1:k.optimal.miho_month]
  B.bar.miho_month[l,,] <- t(X.ut.miho_month[l,,]) %*% Conj(f.bar.miho_month[l,,]) / R.miho_month
  eigenprops <- rbind(eigenprops, cumsum(eigenres$values / sum(eigenres$values))[1:10])
  
  eigenres2 <- eigen(X.ut.miho_month2[l,,] %*% Conj(t(X.ut.miho_month2[l,,])))
  evectors2 <- eigenres2$vectors
  f.bar.miho_month2[l,,] <- sqrt(R.miho_month)*evectors2[,1:k.optimal.miho_month]
  B.bar.miho_month2[l,,] <- t(X.ut.miho_month2[l,,]) %*% Conj(f.bar.miho_month2[l,,]) / R.miho_month
  eigenprops2 <- rbind(eigenprops2, cumsum(eigenres2$values / sum(eigenres2$values))[1:10])
}


# invert ordering
f.tilde.miho_month <- array(0, c(k.optimal.miho_month,N.miho_month,R.miho_month)) # k x n x R
B.miho_month <- array(0, c(p.miho_month,k.optimal.miho_month,N.miho_month))
f.tilde.miho_month2 <- array(0, c(k.optimal.miho_month,N.miho_month,R.miho_month)) # k x n x R
B.miho_month2 <- array(0, c(p.miho_month,k.optimal.miho_month,N.miho_month))

for(l in 1:N.miho_month){
  f.tilde.miho_month[,l,] <- t(f.bar.miho_month[l,,])
  B.miho_month[,,l] <- B.bar.miho_month[l,,]
  f.tilde.miho_month2[,l,] <- t(f.bar.miho_month2[l,,])
  B.miho_month2[,,l] <- B.bar.miho_month2[l,,]
}

# iGFT

# factors on vertex domain
f.miho_month <- array(0, c(k.optimal.miho_month,N.miho_month,R.miho_month)) # k x n x R
f.miho_month2 <- array(0, c(k.optimal.miho_month,N.miho_month,R.miho_month)) # k x n x R

for(j in 1:k.optimal.miho_month){
  f.miho_month[j,,] <- evectors.miho_month %*% f.tilde.miho_month[j,,]
  f.miho_month2[j,,] <- evectors.miho_month2 %*% f.tilde.miho_month2[j,,]
}


# clustering
sum(which(Im(f.miho_month2[3,,])!=0)) # factors are real-valued!
kmeans_res_miho_month2 <- kmeans(t(Re(f.miho_month2[1,,])), centers = 2, nstart = 40) # k-means clustering with factor 1
kmeans_res_miho_month2$cluster 


layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))  # 4:1 width ratio for plot and legend
# Plot area
par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot
plot(kmeans_res_miho_month2$cluster, col = rep(c("red", "red", "blue", "blue", "blue",
                                                 "green", "green", "green", "magenta",
                                                 "magenta", "magenta", "red"), 6), 
     pch = 16, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "cluster", cex.lab=1.4)

axis(2, at = c(1, 2), labels = c("1", "2"))
# Add vertical lines for year separation
abline(v = 0.5, lty = 2)
abline(v = 12.5, lty = 2)
abline(v = 24.5, lty = 2)
abline(v = 36.5, lty = 2)
abline(v = 48.5, lty = 2)
abline(v = 60.5, lty = 2)
abline(v = 72.5, lty = 2)

text_positions <- c(6.5, 18.5, 30.5, 42.5, 54.5, 66.5)  # Middle positions between lines
mtext(text = 2018:2023, side = 1, at = text_positions, line = 1)
par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", legend = c("spring", "summer", "autumn", "winter"),
       col = c("blue", "green", "magenta", "red"), pch = 16, cex = 1.2, bty = "n", y.intersp=1.3)


## loading plot for factor 1
layout(matrix(c(1, 2, 3), nrow = 1, byrow = TRUE), widths = c(3, 3, 1))

par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot
matplot(abs(B.bar.miho_month2[,,1]), type="l", lty=1, col=c("black","red","blue",
                                                            "green","magenta","cyan","orange"),
        xlab="graph frequency index", ylab="loading magnitude", cex.lab=1.6, lwd=1.5)

matplot(Arg(B.bar.miho_month2[,,1]), type="l", lty=1, col=c("black","red","blue",
                                                            "green","magenta","cyan","orange"),
        xlab="graph frequency index", ylab="loading argument", cex.lab=1.6, lwd=1.5)
abline(h=0, lty=2)

par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", inset=c(-0.4,0), legend=colnames(geum_month[,5:10]), 
       col=c("black","red","blue","green","magenta","cyan"), lty=1, bty = "n", lwd=1.5, cex=1.2)



## loading plot for factor 2
layout(matrix(c(1, 2, 3), nrow = 1, byrow = TRUE), widths = c(3, 3, 1))

par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot
matplot(abs(B.bar.miho_month2[,,2]), type="l", lty=1, col=c("black","red","blue",
                                                            "green","magenta","cyan","orange"),
        xlab="graph frequency index", ylab="loading magnitude", cex.lab=1.6, lwd=1.5)

matplot(Arg(B.bar.miho_month2[,,2]), type="l", lty=1, col=c("black","red","blue",
                                                            "green","magenta","cyan","orange"),
        xlab="graph frequency index", ylab="loading argument", cex.lab=1.6, lwd=1.5)
abline(h=0, lty=2)

par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", inset=c(-0.4,0), legend=colnames(geum_month[,5:10]), 
       col=c("black","red","blue","green","magenta","cyan"), lty=1, bty = "n", lwd=1.5, cex=1.2)



# arrow plot
layout(1)
par(mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))   
plot(0, 0, type="n", xlim=c(0, max(abs(B.bar.miho_month2[1,,1]))), 
     ylim=c(0, max(abs(B.bar.miho_month2[1,,2]))), 
     xlab="factor 1", 
     ylab="factor 2")

# Draw arrows from the origin to each point
for(i in 1:6) {
  arrows(0, 0, abs(B.bar.miho_month2[1,i,1]), abs(B.bar.miho_month2[1,i,2]), 
         col="black", lwd=2, length=0.1)  # length controls arrowhead size
  
  # Add text labels next to each arrow
  text(abs(B.bar.miho_month2[1,i,1]), abs(B.bar.miho_month2[1,i,2]), 
       labels=colnames(geum_month[,5:10])[i], pos=4)
}


# par(mfrow=c(2,3))
# plot(f.miho_month2[1,6,], f.miho_month2[2,3,], col=rep(c("red", "red", "blue", "blue", "blue",
#                                                          "green", "green", "green", "magenta",
#                                                          "magenta", "magenta", "red"), 6), pch=19)
# plot(f.miho_month2[1,13,], f.miho_month2[2,1,], col=rep(c("red", "red", "blue", "blue", "blue",
#                                                           "green", "green", "green", "magenta",
#                                                           "magenta", "magenta", "red"), 6), pch=19)
# plot(f.miho_month2[1,15,], f.miho_month2[2,1,], col=rep(c("red", "red", "blue", "blue", "blue",
#                                                           "green", "green", "green", "magenta",
#                                                           "magenta", "magenta", "red"), 6), pch=19)
# plot(f.miho_month2[1,18,], f.miho_month2[2,1,], col=rep(c("red", "red", "blue", "blue", "blue",
#                                                           "green", "green", "green", "magenta",
#                                                           "magenta", "magenta", "red"), 6), pch=19)
# plot(f.miho_month2[1,21,], f.miho_month2[2,1,], col=rep(c("red", "red", "blue", "blue", "blue",
#                                                           "green", "green", "green", "magenta",
#                                                           "magenta", "magenta", "red"), 6), pch=19)
# plot(f.miho_month2[1,24,], f.miho_month2[2,1,], col=rep(c("red", "red", "blue", "blue", "blue",
#                                                           "green", "green", "green", "magenta",
#                                                           "magenta", "magenta", "red"), 6), pch=19)




# if we use original signal of each dimension ???
# Set up a 2x3 grid for the plots and one column for the legend
layout(matrix(c(1, 2, 3, 7, 
                4, 5, 6, 7), nrow = 2, byrow = TRUE), widths = c(2, 2, 2, 1))

# Define parameters for the main plots
par(mar = c(5, 4, 4, 2) + 0.1, xpd = FALSE, oma = c(0, 1, 0, 0))  # Adjust margins for each plot

# Plot 6 times in the 2x3 grid
for (i in 1:6) {
  kmeans_res_miho_month3 <- kmeans(t(X.miho_month[i,,]), centers = 2, nstart = 40)
  plot(kmeans_res_miho_month3$cluster, col = rep(c("red", "red", "blue", "blue", "blue",
                                                   "green", "green", "green", "magenta",
                                                   "magenta", "magenta", "red"), 6),
       pch = 16, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "cluster", cex.lab = 1.4, main=paste("p =",i))
  
  # Add y-axis ticks and labels
  axis(2, at = c(1, 2), labels = c("1", "2"))
  
  # Add vertical lines to separate years
  abline(v = c(0.5, 12.5, 24.5, 36.5, 48.5, 60.5, 72.5), lty = 2)
  
  # Add year labels below the plot
  text_positions <- c(6.5, 18.5, 30.5, 42.5, 54.5, 66.5)
  mtext(text = 2018:2023, side = 1, at = text_positions, line = 1)
}

# Set parameters for the legend plot area
par(mar = c(0, 0, 0, 0))

# Add legend once in the designated legend area
plot.new()
legend("center", legend = c("spring", "summer", "autumn", "winter"),
       col = c("blue", "green", "magenta", "red"), pch = 16, cex = 1.4, bty = "n")


## what if we use averaged signal?
tmppp2 <- t(X.miho_month[1,,])
for(ii in 2:p.miho_month){
  tmppp2 <- tmppp2 + t(X.miho_month[ii,,])
}
kmeans_res_miho_month4 <- kmeans(tmppp2/6, centers = 2, nstart = 40)
kmeans_res_miho_month4$cluster # 1 1 2 2 1 2 / 13 14 15 16 17 18 


layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))  # 4:1 width ratio for plot and legend
# Plot area
par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot
plot(kmeans_res_miho_month4$cluster, col = rep(c("red", "red", "blue", "blue", "blue",
                                                 "green", "green", "green", "magenta",
                                                 "magenta", "magenta", "red"), 6), 
     pch = 16, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "cluster", cex.lab=1.4)

axis(2, at = c(1, 2), labels = c("1", "2"))
# Add vertical lines for year separation
abline(v = 0.5, lty = 2)
abline(v = 12.5, lty = 2)
abline(v = 24.5, lty = 2)
abline(v = 36.5, lty = 2)
abline(v = 48.5, lty = 2)
abline(v = 60.5, lty = 2)
abline(v = 72.5, lty = 2)

text_positions <- c(6.5, 18.5, 30.5, 42.5, 54.5, 66.5)  # Middle positions between lines
mtext(text = 2018:2023, side = 1, at = text_positions, line = 1)
par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", legend = c("spring", "summer", "autumn", "winter"),
       col = c("blue", "green", "magenta", "red"), pch = 16, cex = 1.2, bty = "n")




# draw figure
miho <- readRDS("Miho.RDS")
names(miho)
names(miho$TweedData)
dim(miho$model.mat)

dev.off()
par(oma=c(0.5,0.5,0.5,0.5))
plot(x=miho$TweedPredPoints$Long, y=miho$TweedPredPoints$Lat, 
     xlab = "Longitude", ylab = "Latitude",
     pch=16, cex=.3, cex.lab=1.4)

# ssn_miho@obspoints@SSNPoints[[1]]@point.data$분류번호[c(9,21,20,22,11,13,2,1,3,4,5,8,7,18,27,14,6,28,24,17,12,10,26,15,16)]
# unique(geum_month$id)[as.numeric(na.omit(match(ssn_miho@obspoints@SSNPoints[[1]]@point.data$분류번호, unique(geum_month$id))))]
# unique(geum_month$id)
points(ssn_miho@obspoints@SSNPoints[[1]]@point.data$경도.Degree.[c(9,21,20,22,11,13,2,1,3,4,5,8,7,18,27,14,6,28,24,17,12,10,26,16)], 
       ssn_miho@obspoints@SSNPoints[[1]]@point.data$위도.Degree.[c(9,21,20,22,11,13,2,1,3,4,5,8,7,18,27,14,6,28,24,17,12,10,26,16)],
       pch=16, cex=1.6, col="red")




# January, 2018
g1 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[2,,1], value="Value", ratio=0.9,
                                 min=-0.5, max=0.85, mg=c(4,4,4,4), title="BOD", main.title.size = 20)
g2 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[3,,1], value="Value", ratio=0.9,
                                 min=-0.45, max=0.35, mg=c(4,4,4,4), title="COD", main.title.size = 20)
g3 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[6,,1], value="Value", ratio=0.9,
                                 min=-0.4, max=0.6, mg=c(4,4,4,4), title="TOC", main.title.size = 20)
g4 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[1,,1], value="Value", ratio=0.9,
                                 min=0, max=0.45, mg=c(4,4,4,4), title="DO", main.title.size = 20)
g5 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[4,,1], value="Value", ratio=0.9,
                                 min=0, max=0.85, mg=c(4,4,4,4), title="TN", main.title.size = 20)
g6 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[5,,1], value="Value", ratio=0.9,
                                 min=-0.1, max=0.25, mg=c(4,4,4,4), title="TP", main.title.size = 20)
g7 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = Re(f.miho_month2[1,,1]), value="Value", ratio=0.9,
                                 min=-3.3, max=1.65, mg=c(4,4,4,4), title=expression(F[1]^"Miho"), main.title.size = 20)
g8 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = Re(f.miho_month2[2,,1]), value="Value", ratio=0.9,
                                 min=-1.35, max=3.9, mg=c(4,4,4,4), title=expression(F[2]^"Miho"), main.title.size = 20)
# g9 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = Re(f.miho_month2[3,,1]), value="Value", ratio=0.9,
#                          min=-0.85, max=1.75, mg=c(4,4,4,4), title=expression(F[3]^"Miho"), main.title.size = 20)
grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8, nrow=3)

# July, 2018
g10 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[2,,7], value="Value", ratio=0.9,
                                  min=-0.4, max=0.7, mg=c(4,4,4,4), title="BOD", main.title.size = 20)
g11 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[3,,7], value="Value", ratio=0.9,
                                  min=-0.3, max=0.45, mg=c(4,4,4,4), title="COD", main.title.size = 20)
g12 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[6,,7], value="Value", ratio=0.9,
                                  min=-0.25, max=0.6, mg=c(4,4,4,4), title="TOC", main.title.size = 20)
g13 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[1,,7], value="Value", ratio=0.9,
                                  min=-0.4, max=0.1, mg=c(4,4,4,4), title="DO", main.title.size = 20)
g14 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[4,,7], value="Value", ratio=0.9,
                                  min=-0.75, max=0.15, mg=c(4,4,4,4), title="TN", main.title.size = 20)
g15 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = X.miho_month[5,,7], value="Value", ratio=0.9,
                                  min=-0.1, max=0.25, mg=c(4,4,4,4), title="TP", main.title.size = 20)

g16 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = Re(f.miho_month2[1,,7]), value="Value", ratio=0.9,
                                  min=-1.8, max=1.95, mg=c(4,4,4,4), title=expression(F[1]^"Miho"), main.title.size = 20)
g17 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=8, vertex_color = Re(f.miho_month2[2,,7]), value="Value", ratio=0.9,
                                  min=-2.25, max=3.05, mg=c(4,4,4,4), title=expression(F[2]^"Miho"), main.title.size = 20)


grid.arrange(g10,g11,g12,g13,g14,g15,g16,g17, nrow=3)



# plot with seoulmetro data together
load("seoulmetro_line2.RData")


h1 <- plot_graph_custom_directed(mihoriver_month_directed, e.size=1.3, v.size=6, vertex_color = "black", 
                                 value="value", ratio=0.9, signal=FALSE)
h2 <- plot_graph_custom3(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = "black", 
                         value="value", ratio=0.6, signal=FALSE)
grid.arrange(h1,h2, nrow=1)
