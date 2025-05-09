# may be redundant 

library(RColorBrewer)
library(colourvalues)
library(grDevices)
# library(SDMTools)
library(network)
library(ggraph)
library(tidygraph)
library(dplyr)
library(gasper)
library(readxl)
library(forecast)
library(ggfortify)
library(Metrics)
library(GNAR)
library(DTWBI)
# library(vars)
library(geosphere)
# library(xlsx)
library(scales)
library(igraph)
library(pracma)
library(R.matlab)
library(geosphere)
library(grid)
library(gridBase)
library(gridExtra)
library(expm)
library(Hmisc)
library(cccd)
library(zoo)

source("utils.R")


trade2023 <- read.csv("BACI_HS92_V202501/BACI_HS92_Y2023_V202501.csv", header=TRUE)
trade2023$q <- as.numeric(trade2023$q)
trade2023 <- na.omit(trade2023)

trade2023$value <- trade2023$v * trade2023$q


identical(sort(unique(trade2023$i)), sort(unique(trade2023$j)))
head(trade2023)


countrycode <- read.csv("BACI_HS92_V202501/country_codes_V202501b.csv", header=TRUE)
countrycode_from_ind <- read.csv("economic/Indicators/agriculture.csv", skip=4)

countrycode3digits <- sort(intersect(countrycode$country_iso3, countrycode_from_ind$Country.Code))
countrylist_from_trade <- countrycode[countrycode$country_iso3 %in% countrycode3digits,]$country_name
countrylist_from_ind <- countrycode_from_ind[countrycode_from_ind$Country.Code %in% countrycode3digits,]$Country.Name

country_name_code_tbl <-cbind(countrylist_from_ind, countrycode3digits)
colnames(country_name_code_tbl) <- c("name", "iso3")
country_name_code_tbl <- as.data.frame(country_name_code_tbl)
country_name_code_tbl$lat <- 0
country_name_code_tbl$lon <- 0

coords <- read.csv("country_coord.csv", header=TRUE)
for(i in 1:nrow(country_name_code_tbl)){
  if(country_name_code_tbl$iso3[i]=="CUW"){
    country_name_code_tbl[i,"lat"] <- 12.2135
    country_name_code_tbl[i,"lon"] <- -68.9496
  } else if(country_name_code_tbl$iso3[i]=="SXM"){
    country_name_code_tbl[i,"lat"] <- 18.0347
    country_name_code_tbl[i,"lon"] <- -63.0681
  } else{
    country_name_code_tbl[i,"lat"] <-  coords[which(coords$Alpha.3.code==country_name_code_tbl$iso3[i]),5]
    country_name_code_tbl[i,"lon"] <-  coords[which(coords$Alpha.3.code==country_name_code_tbl$iso3[i]),6]
  }
}

countrycode <- countrycode[-c(18,81,203),] # remove Belgium-Luxembourg (...1998) / Fed. Rep. of Germany (...1990) / Sudan (...2011)
countrycode_used <- countrycode[countrycode$country_iso3 %in% country_name_code_tbl[,2],]


trade2023_used <- trade2023[(trade2023$i %in% countrycode_used$country_code) & (trade2023$j %in% countrycode_used$country_code),]



trade2023_final <- trade2023_used %>%
  # Create new variables i_new and j_new to treat (i, j) and (j, i) as the same pair
  mutate(i_new = pmin(i, j), j_new = pmax(i, j)) %>%
  # Group by these new variables
  group_by(i_new, j_new) %>%
  # Summarize by summing the total
  summarise(
    total = sum(v) / ifelse(sum(i != i_new | j != j_new) > 0, 2, 1),  # Divide by 2 if both directions exist
    .groups = 'drop'
  ) %>%
  # Rename columns back to original names if desired
  rename(i = i_new, j = j_new)


trade2023_final <- trade2023_final %>% arrange(desc(total))
trade2023_final <- trade2023_final[1:round(0.03*nrow(trade2023_final)),]

trade2023_final$total <- trade2023_final$total/10^6
trade2023_final$total_log <- log(trade2023_final$total)

countrycode_final <- countrycode_used[countrycode_used$country_code %in% union(trade2023_final$i, trade2023_final$j),]
countrycode_final_sort <- countrycode_final %>% arrange(country_iso3)

ITW <- list()
# location
country_name_code_tbl_final <- country_name_code_tbl[(country_name_code_tbl$iso3 %in% countrycode_final$country_iso3),]
ITW$xy <- country_name_code_tbl_final[,c("lon","lat")]
rownames(ITW$xy) <- country_name_code_tbl_final$iso3

# weight matrix
e.weight <- matrix(0,nrow=nrow(country_name_code_tbl_final), ncol=nrow(country_name_code_tbl_final))
colnames(e.weight) <- country_name_code_tbl_final$iso3
rownames(e.weight) <- country_name_code_tbl_final$iso3

e.sp.weight <- NULL
for(k in 1:nrow(trade2023_final)){
  i <- which(countrycode_final_sort$country_code==as.numeric(trade2023_final[k,1]))
  j <- which(countrycode_final_sort$country_code==as.numeric(trade2023_final[k,2]))
  e.weight[i,j] <- as.numeric(trade2023_final[k,4])
  e.weight[j,i] <- as.numeric(trade2023_final[k,4])
  e.sp.weight <- rbind(e.sp.weight, c(i,j,as.numeric(trade2023_final[k,4])))
  e.sp.weight <- rbind(e.sp.weight, c(j,i,as.numeric(trade2023_final[k,4])))
}


ITW$A <- e.weight

# sparse weight matrix
ITW$sA <- e.sp.weight[nrow(e.sp.weight):1,]

L.ITW <- laplacian_mat(ITW$A) # laplacian matrix
val1 <- eigensort(L.ITW)
evalues.ITW <- val1$evalues
evectors.ITW <- val1$evectors
# largest eigenvalue
lmax.ITW <- max(evalues.ITW)

N.ITW <- nrow(L.ITW)


plot_graph(ITW)



# load economic data
files <- list.files(path = "./economic/Indicators", pattern = ".csv")
p.ITW <- length(files)
R.ITW <- 34 # 1990~2023
X.ITW <- array(0, c(p.ITW,N.ITW,R.ITW))

tmp2 <- c()
for(i in 1:p.ITW){
  print(paste("i=",i))
  tmp <- read.csv(paste("economic/Indicators",files[i], sep="/"), skip=4)
  print((as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
                                              paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  tmp2 <- c(tmp2, (as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
                                                        paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  # for(l in 1:N.ITW){
  #   
  # }
  # X.ITW[i,,] <- as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
  #     paste("X", 1990:2023, sep="")]) na.approx(na.rm = F) %>% na.locf(na.rm = F) %>% na.locf(fromLast = T, na.rm = F)
}

variable_remove <- c(4,5,7,17,19,26,27,28,33,34,35,36,37) 

tmp3 <- c()
for(i in 1:p.ITW){
  if(i %in% variable_remove){
    next
  }
  tmp <- read.csv(paste("economic/Indicators",files[i], sep="/"), skip=4)
  print((as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
                                              paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  tmp3 <- c(tmp3, (as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
                                                        paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  # for(l in 1:N.ITW){
  #   
  # }
  # X.ITW[i,,] <- as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
  #     paste("X", 1990:2023, sep="")]) na.approx(na.rm = F) %>% na.locf(na.rm = F) %>% na.locf(fromLast = T, na.rm = F)
}

country_name_code_tbl_final[unique(tmp3),]
table(tmp3)

country_remove <- unique(tmp3)

ITW2 <- list()
ITW2$xy <- ITW$xy[-country_remove,]

countrycode_final_sort2 <- countrycode_final_sort[-country_remove, ]
rownames(countrycode_final_sort2) <- 1:N.ITW2
trade2023_final2 <- trade2023_final[(trade2023_final$i %in% countrycode_final_sort2$country_code),]
trade2023_final2 <- trade2023_final2[(trade2023_final2$j %in% countrycode_final_sort2$country_code),]
ITW2$A <- ITW$A[-country_remove, -country_remove] 
e.sp.weight <- NULL
for(k in 1:nrow(trade2023_final2)){
  i <- which(countrycode_final_sort2$country_code==as.numeric(trade2023_final2[k,1]))
  j <- which(countrycode_final_sort2$country_code==as.numeric(trade2023_final2[k,2]))
  # e.weight[i,j] <- as.numeric(trade2023_final2[k,4])
  # e.weight[j,i] <- as.numeric(trade2023_final2[k,4])
  e.sp.weight <- rbind(e.sp.weight, c(i,j,as.numeric(trade2023_final2[k,4])))
  e.sp.weight <- rbind(e.sp.weight, c(j,i,as.numeric(trade2023_final2[k,4])))
}

ITW2$sA <- e.sp.weight[nrow(e.sp.weight):1,]



G20_idx <- c(2,8,9,12,16,23,24,31,32,36,37,41,48,64,65,73,75,77)

countrycode_final_sort_G20 <- countrycode_final_sort2[G20_idx,]


ITW_G20 <- list()
ITW_G20$xy <- ITW2$xy[G20_idx,]

rownames(countrycode_final_sort_G20) <- 1:length(G20_idx)


ITW_G20_2023 <- ITW_G20
ITW_G20_2023$A <- matrix(0,0, nrow=N.ITW_G20, ncol=N.ITW_G20)
colnames(ITW_G20_2023$A) <- colnames(ITW_G20$A)
rownames(ITW_G20_2023$A) <- rownames(ITW_G20$A)
e.sp.weight <- NULL
for(k in 1:nrow(trade2023_final)){
  i <- which(countrycode_final_sort_G20$country_code==as.numeric(trade2023_final[k,1]))
  j <- which(countrycode_final_sort_G20$country_code==as.numeric(trade2023_final[k,2]))
  # e.weight[i,j] <- as.numeric(trade2023_final2[k,4])
  # e.weight[j,i] <- as.numeric(trade2023_final2[k,4])
  ITW_G20_2023$A[i,j] <- as.numeric(trade2023_final[k,4])
  ITW_G20_2023$A[j,i] <- as.numeric(trade2023_final[k,4])
  e.sp.weight <- rbind(e.sp.weight, c(min(i,j),max(i,j),as.numeric(trade2023_final[k,4])))
}

ITW_G20_2023$sA <- e.sp.weight[nrow(e.sp.weight):1,]

L.ITW_G20_2023 <- laplacian_mat(ITW_G20_2023$A) # laplacian matrix
val1 <- eigensort(L.ITW_G20_2023)
evalues.ITW_G20_2023 <- val1$evalues
evectors.ITW_G20_2023 <- val1$evectors
# largest eigenvalue
lmax.ITW_G20_2023 <- max(evalues.ITW_G20_2023)

N.ITW_G20 <- nrow(L.ITW_G20_2023)


p.ITW2 <- length(files) - length(variable_remove)
R.ITW <- 34 # 1990~2023
X.ITW2 <- array(0, c(p.ITW2,N.ITW_G20,R.ITW))


for(i in 1:p.ITW){
  if(i %in% variable_remove){
    next
  }
  tmp <- read.csv(paste("economic/Indicators",files[i], sep="/"), skip=4)
  print((as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% rownames(ITW$xy[-country_remove,]), 
                                              paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  # tmp4 <- c(tmp4, (as.numeric(which(apply(as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3, 
  #                                                       paste("X", 1990:2023, sep="")]), 1, function(x) all(is.na(x)))==TRUE))))
  # for(l in 1:N.ITW){
  #   X.ITW[i,,] <- as.matrix(tmp[tmp$Country.Code %in% country_name_code_tbl_final$iso3,
  #                               paste("X", 1990:2023, sep="")]) na.approx(na.rm = F) %>% na.locf(na.rm = F) %>% na.locf(fromLast = T, na.rm = F)
  # }
}


# tmp4 <- c()
j <- 0
for(i in 1:p.ITW){
  if(i %in% variable_remove){
    next
  }
  j <- j+1
  tmp <- read.csv(paste("economic/Indicators",files[i], sep="/"), skip=4)
  
  for(l in 1:N.ITW2){
    X.ITW2[j,l,] <- as.matrix(tmp[tmp$Country.Code %in% rownames(ITW$xy[-country_remove,]),
                                  paste("X", 1990:2023, sep="")])[l,] %>%  na.approx(na.rm = F) %>% na.locf(na.rm = F) %>% na.locf(fromLast = T, na.rm = F)
  }
}



variables_selected <- c(6,7,8,11,2,4,15,5,14,17,22)

var_eng <- c("GDP growth", "GDP per capita growth", "GDP per capita", "GNI per capita",
             "Current account balance", "Exports of goods and services", "Imports of goods and services",
             "Foreign direct investment (net inflows)", "Gross capital formation",
             "Inflation", "Price level ratio")

var_kor <- c("국내총생산 성장률 (연간 %)", "1인당 GDP 성장률 (연간 %)", "1인당 GDP (현재 미국 달러 기준)", "1인당 총국민소득 (구매력 평가 기준, 국제 달러)",
             "경상수지 (국제수지 기준, 현재 미국 달러 기준)", "상품 및 서비스 수출액 (GDP 대비 %)", "상품 및 서비스 수입액 (GDP 대비 %)",
             "해외직접투자 순유입 (현재 미국 달러 기준)", "총자본형성 (GDP 대비 %)",
             "물가상승률 (GDP 디플레이터 기준)", "구매력 기준 환율과 시장 환율의 비율")

p.ITW_G20.std.selected <- length(variables_selected)

X.ITW_G20_2023 <- X.ITW2[variables_selected,G20_idx,34]
for (i in 1:nrow(X.ITW_G20_2023)) {
  v <- X.ITW_G20_2023[i,]
  X.ITW_G20_2023[i,] <- (v - mean(v)) / sd(v)
}

X.ITW_G20_2023.list <- list()
for(i in 1:p.ITW_G20.std.selected){
  X.ITW_G20_2023.list[[i]] <- X.ITW_G20_2023[i,]
}


### Factor Analysis
num_window <- 50
X.ITW_G20_2023_forFA <- array(0, c(p.ITW_G20.std.selected,N.ITW_G20,num_window))
WB <- windowbank.random(N=length(evalues.ITW_G20_2023), M=num_window, V=evectors.ITW_G20_2023, sigma=0.1, seed=1)
for(i in 1:p.ITW_G20.std.selected){
  X.ITW_G20_2023_forFA[i,,1] <- X.ITW_G20_2023[i,]
  X.ITW_G20_2023_forFA[i,,2:num_window] <- (t(WB)*as.vector(X.ITW_G20_2023_forFA[i,,1]))[,2:num_window]
}


# Fourier coefficients
X.tilde.ITW_G20_2023_forFA <- array(0, c(p.ITW_G20.std.selected,N.ITW_G20,num_window))

for(i in 1:p.ITW_G20.std.selected){
  X.tilde.ITW_G20_2023_forFA[i,,] <- Conj(t(evectors.ITW_G20_2023)) %*% X.ITW_G20_2023_forFA[i,,]
}

# reordering
X.ut.ITW_G20_2023_forFA <- array(0,c(N.ITW_G20,num_window,p.ITW_G20.std.selected)) # n x R x p

for(l in 1:N.ITW_G20){
  X.ut.ITW_G20_2023_forFA[l,,] <- t(X.tilde.ITW_G20_2023_forFA[,l,])
}



# k.optimal.ITW_G20_2023_forFA <- 2
k.optimal.ITW_G20_2023_forFA <- 2

f.bar.ITW_G20_2023_forFA <- array(0, c(N.ITW_G20,num_window,k.optimal.ITW_G20_2023_forFA)) # n x R x k
B.bar.ITW_G20_2023_forFA <- array(0, c(N.ITW_G20,p.ITW_G20.std.selected,k.optimal.ITW_G20_2023_forFA)) # n x p x k

eigenprops <- c()
for(l in 1:N.ITW_G20){
  eigenres <- eigen(X.ut.ITW_G20_2023_forFA[l,,] %*% Conj(t(X.ut.ITW_G20_2023_forFA[l,,])))
  evectors <- eigenres$vectors
  f.bar.ITW_G20_2023_forFA[l,,] <- sqrt(num_window)*evectors[,1:k.optimal.ITW_G20_2023_forFA]
  B.bar.ITW_G20_2023_forFA[l,,] <- t(X.ut.ITW_G20_2023_forFA[l,,]) %*% Conj(f.bar.ITW_G20_2023_forFA[l,,]) / num_window
  eigenprops <- rbind(eigenprops, cumsum(eigenres$values / sum(eigenres$values))[1:10])
}
eigenprops


# invert ordering
f.tilde.ITW_G20_2023_forFA <- array(0, c(k.optimal.ITW_G20_2023_forFA,N.ITW_G20,num_window)) # k x n x R
B.ITW_G20_2023_forFA <- array(0, c(p.ITW_G20.std.selected,k.optimal.ITW_G20_2023_forFA,N.ITW_G20))

for(l in 1:N.ITW_G20){
  f.tilde.ITW_G20_2023_forFA[,l,] <- t(f.bar.ITW_G20_2023_forFA[l,,])
  
  B.ITW_G20_2023_forFA[,,l] <- B.bar.ITW_G20_2023_forFA[l,,]
}

# iGFT
f.ITW_G20_2023_forFA <- array(0, c(k.optimal.ITW_G20_2023_forFA,N.ITW_G20,num_window)) # k x n x R
for(j in 1:k.optimal.ITW_G20_2023_forFA){
  f.ITW_G20_2023_forFA[j,,] <- evectors.ITW_G20_2023 %*% f.tilde.ITW_G20_2023_forFA[j,,]
}


fa_result.ITW_G20_2023 <- fa(t(X.ITW_G20_2023_forFA[1:p.ITW_G20.std.selected,,1]), nfactors = k.optimal.ITW_G20_2023_forFA, rotate = "varimax", fm = "ml", covar=TRUE)



num_window <- 50
X.ITW_G20_2000_forFA <- array(0, c(p.ITW_G20.std.selected,N.ITW_G20,num_window))
WB <- windowbank.random(N=length(evalues.ITW_G20_2000), M=num_window, V=evectors.ITW_G20_2000, sigma=0.1, seed=1)
for(i in 1:p.ITW_G20.std.selected){
  X.ITW_G20_2000_forFA[i,,1] <- X.ITW_G20_2000[i,]
  X.ITW_G20_2000_forFA[i,,2:num_window] <- (t(WB)*as.vector(X.ITW_G20_2000_forFA[i,,1]))[,2:num_window]
}

# Fourier coefficients
X.tilde.ITW_G20_2000_forFA <- array(0, c(p.ITW_G20.std.selected,N.ITW_G20,num_window))

for(i in 1:p.ITW_G20.std.selected){
  X.tilde.ITW_G20_2000_forFA[i,,] <- Conj(t(evectors.ITW_G20_2000)) %*% X.ITW_G20_2000_forFA[i,,]
}

# reordering
X.ut.ITW_G20_2000_forFA <- array(0,c(N.ITW_G20,num_window,p.ITW_G20.std.selected)) # n x R x p

for(l in 1:N.ITW_G20){
  X.ut.ITW_G20_2000_forFA[l,,] <- t(X.tilde.ITW_G20_2000_forFA[,l,])
}



# k.optimal.ITW_G20_2000_forFA <- 3
k.optimal.ITW_G20_2000_forFA <- 1

f.bar.ITW_G20_2000_forFA <- array(0, c(N.ITW_G20,num_window,k.optimal.ITW_G20_2000_forFA)) # n x R x k
B.bar.ITW_G20_2000_forFA <- array(0, c(N.ITW_G20,p.ITW_G20.std.selected,k.optimal.ITW_G20_2000_forFA)) # n x p x k

eigenprops <- c()
for(l in 1:N.ITW_G20){
  eigenres <- eigen(X.ut.ITW_G20_2000_forFA[l,,] %*% Conj(t(X.ut.ITW_G20_2000_forFA[l,,])))
  evectors <- eigenres$vectors
  f.bar.ITW_G20_2000_forFA[l,,] <- sqrt(num_window)*evectors[,1:k.optimal.ITW_G20_2000_forFA]
  B.bar.ITW_G20_2000_forFA[l,,] <- t(X.ut.ITW_G20_2000_forFA[l,,]) %*% Conj(f.bar.ITW_G20_2000_forFA[l,,]) / num_window
  eigenprops <- rbind(eigenprops, cumsum(eigenres$values / sum(eigenres$values))[1:10])
}
eigenprops


# invert ordering
f.tilde.ITW_G20_2000_forFA <- array(0, c(k.optimal.ITW_G20_2000_forFA,N.ITW_G20,num_window)) # k x n x R
B.ITW_G20_2000_forFA <- array(0, c(p.ITW_G20.std.selected,k.optimal.ITW_G20_2000_forFA,N.ITW_G20))

for(l in 1:N.ITW_G20){
  f.tilde.ITW_G20_2000_forFA[,l,] <- t(f.bar.ITW_G20_2000_forFA[l,,])
  
  B.ITW_G20_2000_forFA[,,l] <- B.bar.ITW_G20_2000_forFA[l,,]
}

# iGFT
f.ITW_G20_2000_forFA <- array(0, c(k.optimal.ITW_G20_2000_forFA,N.ITW_G20,num_window)) # k x n x R
for(j in 1:k.optimal.ITW_G20_2000_forFA){
  f.ITW_G20_2000_forFA[j,,] <- evectors.ITW_G20_2000 %*% f.tilde.ITW_G20_2000_forFA[j,,]
}


fa_result.ITW_G20_2000 <- fa(t(X.ITW_G20_2000_forFA[1:p.ITW_G20.std.selected,,1]), nfactors = k.optimal.ITW_G20_2000_forFA, rotate = "varimax", fm = "ml")


fres1 <- plot_graph_custom4(ITW_G20_2000, e.size=1.3, v.size=6, vertex_color = fa_result.ITW_G20_2000$scores[,1], value="Value", ratio=0.6,
                            min=-1.7, max=2.2, mg=c(4,4,4,4), title="Factor 1 (FA, 2000)", main.title.size = 20)
fres2 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = fa_result.ITW_G20_2023$scores[,1], value="Value", ratio=0.6,
                            min=-1.7, max=2.2, mg=c(4,4,4,4), title="Factor 1 (FA, 2023)", main.title.size = 20)
fres3 <- plot_graph_custom4(ITW_G20_2000, e.size=1.3, v.size=6, vertex_color = f.ITW_G20_2000_forFA[1,,1], value="Value", ratio=0.6,
                            min=-1.7, max=2.2, mg=c(4,4,4,4), title="Factor 1 (Proposed, 2000)", main.title.size = 20)
fres4 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = f.ITW_G20_2023_forFA[1,,1], value="Value", ratio=0.6,
                            min=-1.7, max=2.2, mg=c(4,4,4,4), title="Factor 1 (Proposed, 2023)", main.title.size = 20)
fres5 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = f.ITW_G20_2023_forFA[1,,1], value="Value", ratio=0.6,
                            min=min(f.ITW_G20_2023_forFA[1,,1]), max=max(f.ITW_G20_2023_forFA[1,,1]), mg=c(4,4,4,4), title=expression(hat(F)[1]^"TN"), main.title.size = 20)
fres6 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = f.ITW_G20_2023_forFA[2,,1], value="Value", ratio=0.6,
                            min=min(f.ITW_G20_2023_forFA[2,,1]), max=max(f.ITW_G20_2023_forFA[2,,1]), mg=c(4,4,4,4), title=expression(hat(F)[2]^"TN"), main.title.size = 20)


grid.arrange(fres1,fres2,fres3,fres4, nrow=2)

grid.arrange(fres5, fres6, nrow=1)

grid.arrange(p1,p2,fres5, fres6, nrow=2)


## loading plot for factor 1
layout(matrix(c(1, 2), nrow = 1), widths = c(5, 3))  # 4:1 width ratio for plot and legend
par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot

matplot(B.bar.ITW_G20_2000_forFA[,,1], type="l", lty=1, col=c("red", "firebrick", "darkorange", "darkgoldenrod",      # red variations
                                                              "blue", "deepskyblue", "steelblue",          # blue variations
                                                              "green", "forestgreen",                # green variations
                                                              "black", "darkgray"),
        xlab="graph frequency index", ylab="loading", cex.lab=1.4, lwd=1.4)
abline(h=0, lty=2)

par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", inset=c(-0.4,0), legend=var_eng,
       col=c(c("red", "firebrick", "darkorange", "darkgoldenrod",      # red variations
               "blue", "deepskyblue", "steelblue",          # blue variations
               "green", "forestgreen",                # green variations
               "black", "darkgray")), lty=1, lwd=1.4, bty = "n", cex=1.1)


## loading plot for factor 1
layout(matrix(c(1, 2), nrow = 1), widths = c(5, 3))  # 4:1 width ratio for plot and legend
par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot

matplot(B.bar.ITW_G20_2023_forFA[,,1], type="l", lty=1, col=c("red", "firebrick", "darkorange", "darkgoldenrod",      # red variations
                                                              "blue", "deepskyblue", "steelblue",          # blue variations
                                                              "green", "forestgreen",                # green variations
                                                              "black", "darkgray"),
        xlab="graph frequency index", ylab="loading", cex.lab=1.4, lwd=1.4)
abline(h=0, lty=2)

par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", inset=c(-0.4,0), legend=var_eng,
       col=c(c("red", "firebrick", "darkorange", "darkgoldenrod",      # red variations
               "blue", "deepskyblue", "steelblue",          # blue variations
               "green", "forestgreen",                # green variations
               "black", "darkgray")), lty=1, lwd=1.4, bty = "n", cex=1.1)



fp1 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = evectors.ITW_G20_2023[,4], value="Value", ratio=0.6,
                          min=min(evectors.ITW_G20_2023[,4]), max=max(evectors.ITW_G20_2023[,4]), mg=c(4,4,4,4), title=expression(v[4]^"TN"), main.title.size = 20)

fp2 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = evectors.ITW_G20_2023[,17], value="Value", ratio=0.6,
                          min=min(evectors.ITW_G20_2023[,17]), max=max(evectors.ITW_G20_2023[,17]), mg=c(4,4,4,4), title=expression(v[17]^"TN"), main.title.size = 20)

grid.arrange(fp1,fp2, nrow=1)

gp1 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = X.ITW_G20_2023[1,], value="Value", ratio=0.6,
                          min=min(X.ITW_G20_2023[1,]), max=max(X.ITW_G20_2023[1,]), mg=c(4,4,4,4), title="GDP growth", main.title.size = 20)
gp2 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = X.ITW_G20_2023[5,], value="Value", ratio=0.6,
                          min=min(X.ITW_G20_2023[5,]), max=max(X.ITW_G20_2023[5,]), mg=c(4,4,4,4), title="Current account balance", main.title.size = 20)
gp3 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = X.ITW_G20_2023[8,], value="Value", ratio=0.6,
                          min=min(X.ITW_G20_2023[8,]), max=max(X.ITW_G20_2023[8,]), mg=c(4,4,4,4), title="Foreign direct investment", main.title.size = 20)
gp4 <- plot_graph_custom4(ITW_G20_2023, e.size=1.3, v.size=6, vertex_color = X.ITW_G20_2023[10,], value="Value", ratio=0.6,
                          min=min(X.ITW_G20_2023[10,]), max=max(X.ITW_G20_2023[10,]), mg=c(4,4,4,4), title="Inflation", main.title.size = 20)


grid.arrange(gp1,gp2,gp3,gp4,fres5,fres6,fp1,fp2, nrow=2)


# reconstruction error
tmp5 <- c()
for(pp in 1:p.ITW_G20.std.selected){
  tmp5 <- c(tmp5, X.ITW_G20_2000_forFA[pp,,1] - evectors.ITW_G20_2023 %*% diag(B.bar.ITW_G20_2023_forFA[,pp,1]) %*% t(Conj(evectors.ITW_G20_2023)) %*% f.ITW_G20_2023_forFA[1,,1] -
              evectors.ITW_G20_2023 %*% diag(B.bar.ITW_G20_2023_forFA[,pp,2]) %*% t(Conj(evectors.ITW_G20_2023)) %*% f.ITW_G20_2023_forFA[2,,1] - 
              evectors.ITW_G20_2023 %*% diag(B.bar.ITW_G20_2023_forFA[,pp,3]) %*% t(Conj(evectors.ITW_G20_2023)) %*% f.ITW_G20_2023_forFA[3,,1])
}
rmse.ITW_G20_2023_forFA <- norm(tmp5, type="2")

rmse.ITW_G20_2023_forFA
norm(t(X.ITW_G20_2000_forFA[,,1]) - fa_result.ITW_G20_2000$scores %*% t(fa_result.ITW_G20_2000$loadings), "F")




## loading plot for factor 1 (2023)
layout(matrix(c(1, 2), nrow = 1), widths = c(5, 3))  # 4:1 width ratio for plot and legend
par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot

matplot(B.bar.ITW_G20_2023_forFA[,,1], type="l", lty=1, col=c("red", "firebrick", "darkorange", "darkgoldenrod",      # red variations
                                                              "blue", "deepskyblue", "steelblue",          # blue variations
                                                              "green", "forestgreen",                # green variations
                                                              "black", "darkgray"),
        xlab="graph frequency index", ylab="loading", cex.lab=1.4, lwd=1.4)
abline(h=0, lty=2)

par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", inset=c(-0.4,0), legend=var_eng,
       col=c(c("red", "firebrick", "darkorange", "darkgoldenrod",      # red variations
               "blue", "deepskyblue", "steelblue",          # blue variations
               "green", "forestgreen",                # green variations
               "black", "darkgray")), lty=1, lwd=1.4, bty = "n", cex=1.1)


## loading plot for factor 2 (2023)
layout(matrix(c(1, 2), nrow = 1), widths = c(5, 3))  # 4:1 width ratio for plot and legend
par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot

matplot(B.bar.ITW_G20_2023_forFA[,,2], type="l", lty=1, col=c("red", "firebrick", "darkorange", "darkgoldenrod",      # red variations
                                                              "blue", "deepskyblue", "steelblue",          # blue variations
                                                              "green", "forestgreen",                # green variations
                                                              "black", "darkgray"),
        xlab="graph frequency index", ylab="loading", cex.lab=1.4, lwd=1.4)
abline(h=0, lty=2)

par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", inset=c(-0.4,0), legend=var_eng,
       col=c(c("red", "firebrick", "darkorange", "darkgoldenrod",      # red variations
               "blue", "deepskyblue", "steelblue",          # blue variations
               "green", "forestgreen",                # green variations
               "black", "darkgray")), lty=1, lwd=1.4, bty = "n", cex=1.1)



layout(matrix(c(1, 3, 2, 4), nrow = 2), widths = c(3,2))  # 2 plots and 2 legends stacked vertically

# --- Factor 1 Plot ---
par(mar = c(5, 5, 4, 1) + 0.1, xpd = FALSE)
matplot(B.bar.ITW_G20_2023_forFA[,,1], type = "l", lty = 1,
        col = c("red", "firebrick", "darkorange", "darkgoldenrod",    # red
                "blue", "deepskyblue", "steelblue",                   # blue
                "green", "forestgreen",                               # green
                "black", "darkgray"),                                 # neutral
        xlab = "graph frequency index", ylab = "loading",
        cex.lab = 1.4, lwd = 1.4)
abline(h = 0, lty = 2)

# --- Factor 1 Legend ---
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", inset = c(0, 0), legend = var_eng,
       col = c("red", "firebrick", "darkorange", "darkgoldenrod",
               "blue", "deepskyblue", "steelblue",
               "green", "forestgreen",
               "black", "darkgray"),
       lty = 1, lwd = 1.4, bty = "n", cex = 1.1)

# --- Factor 2 Plot ---
par(mar = c(5, 5, 4, 1) + 0.1, xpd = FALSE)
matplot(B.bar.ITW_G20_2023_forFA[,,2], type = "l", lty = 1,
        col = c("red", "firebrick", "darkorange", "darkgoldenrod",
                "blue", "deepskyblue", "steelblue",
                "green", "forestgreen",
                "black", "darkgray"),
        xlab = "graph frequency index", ylab = "loading",
        cex.lab = 1.4, lwd = 1.4)
abline(h = 0, lty = 2)

# --- Factor 2 Legend ---
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", inset = c(0, 0), legend = var_eng,
       col = c("red", "firebrick", "darkorange", "darkgoldenrod",
               "blue", "deepskyblue", "steelblue",
               "green", "forestgreen",
               "black", "darkgray"),
       lty = 1, lwd = 1.4, bty = "n", cex = 1.1)




## Arrow plot
# Set up a 2x3 grid for the plots and one column for the legend
layout(matrix(c(1, 2, 3, 4), nrow = 1, byrow = TRUE), widths = c(3, 3, 3, 3))


par(mar = c(5, 6, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot

for(l in c(5,10,15)){
  plot(0, 0, type="n", xlim=c(0, max(abs(B.bar.ITW_G20_2023_forFA[l,,1]))), 
       ylim=c(0, max(abs(B.bar.ITW_G20_2023_forFA[l,,2]))), 
       xlab="loading 1", 
       ylab="loading 2", cex.lab=2, lwd=1.2, main=bquote(lambda[.(l)]^"TN"), cex.main=2)
  
  # Draw arrows from the origin to each point
  for(i in 1:11) {
    arrows(0, 0, abs(B.bar.ITW_G20_2023_forFA[l,i,1]), abs(B.bar.ITW_G20_2023_forFA[l,i,2]), 
           col=c("red", "firebrick", "darkorange", "darkgoldenrod",
                 "blue", "deepskyblue", "steelblue",
                 "green", "forestgreen",
                 "black", "darkgray")[i], lwd=2, length=0.1)  # length controls arrowhead size
    
  }
}

par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", inset=c(-0.4,0), legend=var_eng,
       col=c("red", "firebrick", "darkorange", "darkgoldenrod",
             "blue", "deepskyblue", "steelblue",
             "green", "forestgreen",
             "black", "darkgray"), lty=1, lwd=2, bty = "n", cex=1.5)
