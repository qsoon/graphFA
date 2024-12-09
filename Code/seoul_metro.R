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
library(ggplot2)


source("Code/source.R")
source("Code/utils.R")
################
## graph load ##
################

# station Korean, English name
station.name <- read.csv("Data/seoulmetro/station_name.csv", 
                         header=TRUE, fileEncoding = "euc-kr")
colnames(station.name) <- c("code", "name_kor", "name", "line", "external_code")


# station latitude, longitude info
station.loc <- read.csv("Data/seoulmetro/station_location.csv", 
                        header=TRUE, fileEncoding = "euc-kr")
colnames(station.loc) <- c("ID", "name_kor", "line", "lon", "lat")

station.loc2 <- read.csv("Data/seoulmetro/station_location2.csv", 
                         header=TRUE, fileEncoding = "euc-kr")
station.loc2 <- station.loc2[,c(2:6)]
colnames(station.loc2) <- c("line", "station_num", "name_kor", "lon", "lat")


# distance between stations
station.distance <- read.csv("Data/seoulmetro/station_distance.csv", 
                             header=TRUE, fileEncoding = "euc-kr")

station.distance <- station.distance[, c(2:5)]
colnames(station.distance) <- c("line", "name_kor", "btwdist", "cumdist")

# 2021 hourly getting on/off info for each station
hourly.pplnum.2021 <- read.csv("Data/seoulmetro/hourly_pplnum_2021.csv", 
                               header=TRUE, fileEncoding = "euc-kr")
hourly.pplnum.2021 <- hourly.pplnum.2021[, -1]
colnames(hourly.pplnum.2021) <- c("date", "line", "station_num", "name_kor", "type",
                                  paste(rep("t",19), c(1:19), sep=""))


# We have getting on/off info only for line 1~8, target line = line 1~8
target_line <- c("01호선", "02호선", "03호선", "04호선",
                 "05호선", "06호선", "07호선", "08호선",
                 "1호선", "2호선", "3호선", "4호선",
                 "5호선", "6호선", "7호선", "8호선")

############################
#### Data preprocessing ####
############################

station.name$name_kor <- as.character(station.name$name_kor)
station.name <- station.name[station.name$name_kor != "이수",]
station.name$name_kor[which(station.name$name_kor == "4?19민주묘지")] <- "4.19민주묘지"
station.name <- station.name[station.name$line %in% target_line,]
station.name$line <- as.character(station.name$line)
station.name$line <- as.character(sapply(station.name$line,
                                         function(x) {strsplit(x, split="")[[1]][2]}))

station.loc$name_kor <- as.character(sapply(as.character(station.loc$name_kor), 
                                            function(x) {strsplit(x, split="\\(")[[1]][1]}))
station.loc <- station.loc[station.loc$line %in% target_line,]
station.loc <- as.data.frame(dplyr::select(station.loc, name_kor, lon, lat))

station.loc$name_kor[which(station.loc$name_kor=="이수")] <- "총신대입구"
# averaging location if there are several line passing the same station
station.loc <- as.data.frame(station.loc %>% group_by(name_kor) %>% 
                               summarise(lat = mean(lat), lon = mean(lon)))

# included in station.loc. lack of info 
station.loc2$name_kor <- as.character(station.loc2$name_kor)
station.loc2$name_kor[which(station.loc2$name_kor=="서울")] <- "서울역"

station.distance$name_kor <- as.character(station.distance$name_kor)
station.distance$name_kor[which(station.distance$name_kor=="이수")] <- "총신대입구"

# remove "(", ")"  in the station name
hourly.pplnum.2021$name_kor <- as.character(sapply(as.character(hourly.pplnum.2021$name_kor), 
                                                   function(x) {strsplit(x, split="\\(")[[1]][1]}))
hourly.pplnum.2021$name_kor[which(hourly.pplnum.2021$name_kor=="이수")] <- "총신대입구"
hourly.pplnum.2021 <- dplyr::select(hourly.pplnum.2021,-c("station_num"))
hourly.pplnum.2021 <- as.data.frame(hourly.pplnum.2021)

hourly.pplnum.2021_getin <- hourly.pplnum.2021[hourly.pplnum.2021$type=="승차",]
hourly.pplnum.2021_getout <- hourly.pplnum.2021[hourly.pplnum.2021$type=="하차",]

hourly.pplnum.2021_getin_line2 <- hourly.pplnum.2021_getin[hourly.pplnum.2021_getin$line=="2호선",]
hourly.pplnum.2021_getout_line2 <- hourly.pplnum.2021_getout[hourly.pplnum.2021_getout$line=="2호선",]

hourly.pplnum.2021_getin <- dplyr::select(hourly.pplnum.2021_getin,-c("line","type"))
hourly.pplnum.2021_getout <- dplyr::select(hourly.pplnum.2021_getout,-c("line","type"))
hourly.pplnum.2021_getin_line2 <- dplyr::select(hourly.pplnum.2021_getin_line2,-c("line","type"))
hourly.pplnum.2021_getout_line2 <- dplyr::select(hourly.pplnum.2021_getout_line2,-c("line","type"))


hourly.pplnum.2021_getin[,3:ncol(hourly.pplnum.2021_getin)] <- 
  sapply(hourly.pplnum.2021_getin[,3:ncol(hourly.pplnum.2021_getin)], as.numeric)

hourly.pplnum.2021_getout[,3:ncol(hourly.pplnum.2021_getout)] <- 
  sapply(hourly.pplnum.2021_getout[,3:ncol(hourly.pplnum.2021_getout)], as.numeric)

hourly.pplnum.2021_getin_line2[,3:ncol(hourly.pplnum.2021_getin_line2)] <- 
  sapply(hourly.pplnum.2021_getin_line2[,3:ncol(hourly.pplnum.2021_getin_line2)], as.numeric)

hourly.pplnum.2021_getout_line2[,3:ncol(hourly.pplnum.2021_getout_line2)] <- 
  sapply(hourly.pplnum.2021_getout_line2[,3:ncol(hourly.pplnum.2021_getout_line2)], as.numeric)

# summing the numbers of people across type, line
hourly.pplnum.2021_getin <- as.data.frame(hourly.pplnum.2021_getin %>% group_by(date, name_kor) %>% 
                                            summarise(across(everything(),sum)))

hourly.pplnum.2021_getout <- as.data.frame(hourly.pplnum.2021_getout %>% group_by(date, name_kor) %>% 
                                             summarise(across(everything(),sum)))

hourly.pplnum.2021_getin_line2 <- as.data.frame(hourly.pplnum.2021_getin_line2 %>% group_by(date, name_kor) %>% 
                                                  summarise(across(everything(),sum)))

hourly.pplnum.2021_getout_line2 <- as.data.frame(hourly.pplnum.2021_getout_line2 %>% group_by(date, name_kor) %>% 
                                                   summarise(across(everything(),sum)))


hourly.pplnum.2021_getin %>% group_by(date) %>% summarise(count = n())
hourly.pplnum.2021_getout %>% group_by(date) %>% summarise(count = n())
hourly.pplnum.2021_getin_line2 %>% group_by(date) %>% summarise(count = n())
hourly.pplnum.2021_getout_line2 %>% group_by(date) %>% summarise(count = n())

# stations with no info
station_removed <- c("신내", "강일", "하남검단산", "하남시청")

hourly.pplnum.2021_getin <- hourly.pplnum.2021_getin[!(hourly.pplnum.2021_getin$name_kor %in% station_removed), ]
hourly.pplnum.2021_getin %>% group_by(date) %>% summarise(count = n())

hourly.pplnum.2021_getout <- hourly.pplnum.2021_getout[!(hourly.pplnum.2021_getout$name_kor %in% station_removed), ]
hourly.pplnum.2021_getout %>% group_by(date) %>% summarise(count = n())

hourly.pplnum.2021_getin_line2 <- hourly.pplnum.2021_getin_line2[!(hourly.pplnum.2021_getin_line2$name_kor %in% station_removed), ]
hourly.pplnum.2021_getin_line2 %>% group_by(date) %>% summarise(count = n())

hourly.pplnum.2021_getout_line2 <- hourly.pplnum.2021_getout_line2[!(hourly.pplnum.2021_getout_line2$name_kor %in% station_removed), ]
hourly.pplnum.2021_getout_line2 %>% group_by(date) %>% summarise(count = n())


########################
#### aggregate info ####
########################

### Line 2 get in ###

# add location to getting on/off info
station.info.getin_line2 <- inner_join(hourly.pplnum.2021_getin_line2, station.loc, by='name_kor')

# add station's English name
station.info.getin_line2 <- inner_join(station.info.getin_line2, station.name[, c("name_kor", "name")][
  !duplicated(station.name[, c("name_kor", "name")]),], by='name_kor')

station.info.getin_line2$name <- as.character(station.info.getin_line2$name)

# column reordering
station.info.getin_line2 <- dplyr::select(station.info.getin_line2, 1,2, ncol(station.info.getin_line2):(ncol(station.info.getin_line2)-2),
                                          3:(ncol(station.info.getin_line2)-3))

# our target stations
target_station_line2 <- unique(station.info.getin_line2$name) # total 50 stations
target_station_kor_line2 <- unique(station.info.getin_line2$name_kor)

# fill NA values
tmp <- target_station_line2
tmp.kor <- target_station_kor_line2
tmp.loc <- sapply(station.info.getin_line2[!duplicated(station.info.getin_line2[,c("lon","lat")]),
                                           c("lon","lat")], as.numeric)

datenum <- length(unique(station.info.getin_line2$date))

tmp.mat <- as.data.frame(cbind(rep(tmp.kor, datenum), rep(tmp, datenum),
                               rep(tmp.loc[,1], datenum), rep(tmp.loc[,2], datenum)))

colnames(tmp.mat) <- c("name_kor", "name", "lon", "lat")
tmp.mat$name_kor <- as.character(tmp.mat$name_kor)
tmp.mat$name <- as.character(tmp.mat$name)
tmp.mat$lon <- as.numeric(as.character(tmp.mat$lon))
tmp.mat$lat <- as.numeric(as.character(tmp.mat$lat))
tmp.mat$date <- as.factor(as.character(rep(seq(as.Date("2021-01-01"), by = "day", length.out = datenum),
                                           each = length(target_station_line2))))

tmp2 <- left_join(tmp.mat, as.data.frame(station.info.getin_line2), 
                  by=c("date"="date","name_kor"="name_kor", 
                       "name"="name"))


station.info.getin_line2 <- dplyr::select(tmp2, 5,1,2,3,4,8:ncol(tmp2))
colnames(station.info.getin_line2)[4:5] <- c("lon", "lat")

station.info.getin_line2.timediv <- station.info.getin_line2
# station.info.getin_line2.timediv$T1 <- station.info.getin_line2.timediv$t1
# station.info.getin_line2.timediv$T2 <- rowSums(station.info.getin_line2.timediv[,7:10])
# station.info.getin_line2.timediv$T3 <- rowSums(station.info.getin_line2.timediv[,11:16])
# station.info.getin_line2.timediv$T4 <- rowSums(station.info.getin_line2.timediv[,17:20])
# station.info.getin_line2.timediv$T5 <- rowSums(station.info.getin_line2.timediv[,21:24])
# 
# station.info.getin_line2.timediv <- station.info.getin_line2.timediv[, c(1,2,3,4,5,25:29)]



### Line 2 get out ###

# add location to getting on/off info
station.info.getout_line2 <- inner_join(hourly.pplnum.2021_getout_line2, station.loc, by='name_kor')

# add station's English name
station.info.getout_line2 <- inner_join(station.info.getout_line2, station.name[, c("name_kor", "name")][
  !duplicated(station.name[, c("name_kor", "name")]),], by='name_kor')

station.info.getout_line2$name <- as.character(station.info.getout_line2$name)

# column reordering
station.info.getout_line2 <- dplyr::select(station.info.getout_line2, 1,2, ncol(station.info.getout_line2):(ncol(station.info.getout_line2)-2),
                                           3:(ncol(station.info.getout_line2)-3))

# our target stations
target_station_line2 <- unique(station.info.getout_line2$name) # total 50 stations
target_station_kor_line2 <- unique(station.info.getout_line2$name_kor)

# fill NA values
tmp <- target_station_line2
tmp.kor <- target_station_kor_line2
tmp.loc <- sapply(station.info.getout_line2[!duplicated(station.info.getout_line2[,c("lon","lat")]),
                                            c("lon","lat")], as.numeric)

datenum <- length(unique(station.info.getout_line2$date))

tmp.mat <- as.data.frame(cbind(rep(tmp.kor, datenum), rep(tmp, datenum),
                               rep(tmp.loc[,1], datenum), rep(tmp.loc[,2], datenum)))

colnames(tmp.mat) <- c("name_kor", "name", "lon", "lat")
tmp.mat$name_kor <- as.character(tmp.mat$name_kor)
tmp.mat$name <- as.character(tmp.mat$name)
tmp.mat$lon <- as.numeric(as.character(tmp.mat$lon))
tmp.mat$lat <- as.numeric(as.character(tmp.mat$lat))
tmp.mat$date <- as.factor(as.character(rep(seq(as.Date("2021-01-01"), by = "day", length.out = datenum),
                                           each = length(target_station_line2))))

tmp2 <- left_join(tmp.mat, as.data.frame(station.info.getout_line2), 
                  by=c("date"="date","name_kor"="name_kor", 
                       "name"="name"))


station.info.getout_line2 <- dplyr::select(tmp2, 5,1,2,3,4,8:ncol(tmp2))
colnames(station.info.getout_line2)[4:5] <- c("lon", "lat")

station.info.getout_line2.timediv <- station.info.getout_line2
# station.info.getout_line2.timediv$T1 <- station.info.getout_line2.timediv$t1
# station.info.getout_line2.timediv$T2 <- rowSums(station.info.getout_line2.timediv[,7:10])
# station.info.getout_line2.timediv$T3 <- rowSums(station.info.getout_line2.timediv[,11:16])
# station.info.getout_line2.timediv$T4 <- rowSums(station.info.getout_line2.timediv[,17:20])
# station.info.getout_line2.timediv$T5 <- rowSums(station.info.getout_line2.timediv[,21:24])

# station.info.getout_line2.timediv <- station.info.getout_line2.timediv[, c(1,2,3,4,5,25:29)]


#########################################
#### edge weight matrix construction ####
#########################################
station.distance <- as.data.frame(station.distance) 
# remove "(", ")" in the name
station.distance$name_kor <- as.character(sapply(as.character(station.distance$name_kor), 
                                                 function(x) {strsplit(x, split="\\(")[[1]][1]}))
# remove blank at the end of name
station.distance$name_kor <- as.character(sapply(as.character(station.distance$name_kor), 
                                                 function(x) {strsplit(x, split=" ")[[1]][1]}))
station.distance$name_kor[which(station.distance$name_kor=="신내역")] <- "신내"
# add station's English name
station.distance_line2 <- station.distance[station.distance$name_kor %in% target_station_kor_line2,]

station.distance_line2 <- inner_join(station.distance_line2, 
                                     station.info.getin_line2[!duplicated(station.info.getin_line2[,c("name_kor","name")]),
                                                              c("name_kor", "name")], by='name_kor')

station.distance_line2[station.distance_line2$name_kor=="산성",]$btwdist <-
  station.distance_line2[station.distance_line2$name_kor=="산성",]$btwdist + 1.5 # 남위례역이 hourly data 정보에 없어서 빠지기 때문에 간격 더해줌


station.distance_line2[station.distance_line2$name_kor=="미사",]$btwdist <-
  station.distance_line2[station.distance_line2$name_kor=="미사",]$btwdist + 0.8 # 강일이 hourly data 정보에 없어서 빠지기 때문에 간격 더해줌

station.distance_line2$name <- as.character(station.distance_line2$name)

e.weight <- matrix(0,nrow=length(target_station_line2), ncol=length(target_station_line2))
colnames(e.weight) <- target_station_line2
rownames(e.weight) <- target_station_line2

for(i in 2:2){
  tmp <-station.distance_line2[station.distance_line2$line==i,]
  if(i==2){
    n <- 44 # circular line. 54th line : City Hall again
    e.weight["Seongsu", "Yongdap"] <- tmp$btwdist[45]
    e.weight["Yongdap", "Seongsu"] <- tmp$btwdist[45]
    for(j in (n+1):47){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
    e.weight["Sindorim", "Dorimcheon"] <- tmp$btwdist[49]
    e.weight["Dorimcheon", "Sindorim"] <- tmp$btwdist[49]
    for(j in 49:50){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
  } else if(i==5){
    n <- 46 # Hanam Pungsan is the end station of one line in Line5
    e.weight["Gangdong", "Dunchon-dong"] <- tmp$btwdist[47]
    e.weight["Dunchon-dong", "Gangdong"] <- tmp$btwdist[47]
    for(j in (n+1):52){
      e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
      e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
    }
  } else{
    n <- nrow(tmp)
  }
  for(j in 1:(n-1)){
    e.weight[tmp$name[j], tmp$name[j+1]] <- tmp$btwdist[j+1]
    e.weight[tmp$name[j+1], tmp$name[j]] <- tmp$btwdist[j+1]
  }
}

e.weight.old <- e.weight
e.weight.old[e.weight.old!=0] <- exp(-e.weight.old[e.weight.old!=0])
e.weight[e.weight!=0] <- exp(-(e.weight[e.weight!=0])^2 / mean(e.weight[e.weight!=0])^2)

e.weight_circle <- e.weight[(!colnames(e.weight) %in% c("Yongdap","Yongdu","Sinjeongnegeori","Sinseol-dong","Dorimcheon","Yangcheon-gu Office","Sindap")),
                            (!colnames(e.weight) %in% c("Yongdap","Yongdu","Sinjeongnegeori","Sinseol-dong","Dorimcheon","Yangcheon-gu Office","Sindap"))]

e.sp.weight <- NULL
e.color <- c() # for line color
color.cand <- c("blue", "yellowgreen", "orangered", "cyan",
                "darkorchid", "chocolate3", "darkolivegreen", "hotpink")
for(i in 2:2){
  tmp <-station.distance_line2[station.distance_line2$line==i,]
  if(i==2){
    n <- 44 # circular line. 54th line : City Hall again
    e.sp.weight <- rbind(e.sp.weight, c("Seongsu", "Yongdap", tmp$btwdist[45]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, c("Yongdap", "Seongsu", tmp$btwdist[45]))
    e.color <- c(e.color, color.cand[i])
    for(j in (n+1):47){
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
    e.sp.weight <- rbind(e.sp.weight, c("Sindorim", "Dorimcheon", tmp$btwdist[49]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, c("Dorimcheon", "Sindorim", tmp$btwdist[49]))
    e.color <- c(e.color, color.cand[i])
    
    for(j in 49:50){
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, 
                           c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
  } else if(i==5){
    n <- 46 # Hanam Pungsan is the end station of one line in Line5
    e.sp.weight<- rbind(e.sp.weight, c("Gangdong", "Dunchon-dong", tmp$btwdist[47]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight<- rbind(e.sp.weight, c("Dunchon-dong", "Gangdong", tmp$btwdist[47]))
    e.color <- c(e.color, color.cand[i])
    for(j in (n+1):52){
      e.sp.weight <- rbind(e.sp.weight, c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
      e.sp.weight <- rbind(e.sp.weight, c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
      e.color <- c(e.color, color.cand[i])
    }
  } else{
    n <- nrow(tmp)
  }
  for(j in 1:(n-1)){
    e.sp.weight <- rbind(e.sp.weight, 
                         c(tmp$name[j], tmp$name[j+1], tmp$btwdist[j+1]))
    e.color <- c(e.color, color.cand[i])
    e.sp.weight <- rbind(e.sp.weight, 
                         c(tmp$name[j+1], tmp$name[j], tmp$btwdist[j+1]))
    e.color <- c(e.color, color.cand[i])
  }
}

tmp <- as.data.frame(e.sp.weight)
tmp$V1 <- as.character(tmp$V1)
tmp$V2 <- as.character(tmp$V2)
tmp2 <- exp(-as.numeric(as.character(tmp$V3)))
tmp$V3 <- exp(-(as.numeric(as.character(tmp$V3)))^2 / (mean(as.numeric(as.character(tmp$V3))))^2)

e.sp.weight <- tmp
e.sp.weight.old <- e.sp.weight
e.sp.weight.old[,3] <- tmp2



seoulmetro_line2 <- list()
# location
seoulmetro_line2$xy <- sapply(station.info.getin_line2[!duplicated(station.info.getin_line2[,c("lon","lat")]),
                                                       c("lon","lat")], as.numeric)
rownames(seoulmetro_line2$xy) <- target_station_line2

# weight matrix
seoulmetro_line2$A <- e.weight

# sparse weight matrix
seoulmetro_line2$sA <- e.sp.weight

# dist matrix
distmat <- e.weight.old
distmat[distmat!=0] <- -log(e.weight.old[e.weight.old!=0])

tmp <- e.sp.weight.old
tmp[,3] <- -log(tmp[,3])

seoulmetro_line2$dist <- distmat
seoulmetro_line2$sdist <- tmp

plot_graph(seoulmetro_line2)


L.metro_line2 <- laplacian_mat(seoulmetro_line2$A) # laplacian matrix
val1 <- eigensort(L.metro_line2)
evalues.metro_line2 <- val1$evalues
evectors.metro_line2 <- val1$evectors
# largest eigenvalue
lmax.metro_line2 <- max(evalues.metro_line2)

N.metro_line2 <- nrow(L.metro_line2)

#################
## signal load ##
#################
subway2019 <- read.csv("Data/seoulmetro/metro2019.csv", 
                       header=TRUE, fileEncoding = "euc-kr")
subway2019 <- subway2019[,c(1,2,4:25)]

colnames(subway2019) <- c("date","line","name_kor","type",paste("t",1:20,sep=""))
subway2019_line2 <- subway2019[subway2019$line=="2호선",]
subway2019_line2[subway2019_line2$name_kor=="동대문역사문화공원(DDP)",]$name_kor <- "동대문역사문화공원"
subway2019_line2[subway2019_line2$name_kor=="왕십리(성동구청)",]$name_kor <- "왕십리"
subway2019_line2[subway2019_line2$name_kor=="구의(광진구청)",]$name_kor <- "구의"
subway2019_line2[subway2019_line2$name_kor=="강변(동서울터미널)",]$name_kor <- "강변"
subway2019_line2[subway2019_line2$name_kor=="잠실(송파구청)",]$name_kor <- "잠실"
subway2019_line2[subway2019_line2$name_kor=="삼성(무역센터)",]$name_kor <- "삼성"
subway2019_line2[subway2019_line2$name_kor=="교대(법원.검찰청)",]$name_kor <- "교대"
subway2019_line2[subway2019_line2$name_kor=="서울대입구(관악구청)",]$name_kor <- "서울대입구"
subway2019_line2[subway2019_line2$name_kor=="대림(구로구청)",]$name_kor <- "대림"
subway2019_line2[subway2019_line2$name_kor=="충정로(경기대입구)",]$name_kor <- "충정로"
subway2019_line2[subway2019_line2$name_kor=="용두(동대문구청)",]$name_kor <- "용두"
subway2019_line2[subway2019_line2$name_kor=="동대문역사문화공원(DDP)(DDP)",]$name_kor <- "동대문역사문화공원"
subway2019_line2[subway2019_line2$name_kor=="낙성대(강감찬)",]$name_kor <- "낙성대"

sort(unique(subway2019_line2$name_kor))

subway2019_line2_getin <- subway2019_line2[subway2019_line2$type=="승차",]
subway2019_line2_getout <- subway2019_line2[subway2019_line2$type=="하차",]

subway2019_line2_getin <- subway2019_line2_getin[,c(1,3,5:24)]
subway2019_line2_getout <- subway2019_line2_getout[,c(1,3,5:24)]

for(i in 1:(nrow(subway2019_line2_getin)/50)){
  subway2019_line2_getin[(50*(i-1)+1):(50*i),] <-
    subway2019_line2_getin[50*(i-1)+order(unique(subway2019_line2$name_kor)),]
}

for(i in 1:(nrow(subway2019_line2_getout)/50)){
  subway2019_line2_getout[(50*(i-1)+1):(50*i),] <-
    subway2019_line2_getout[50*(i-1)+order(unique(subway2019_line2$name_kor)),]
}

subway2019_line2_getin$name <- station.info.getin_line2.timediv$name[1:50]
subway2019_line2_getin <- subway2019_line2_getin[,c(1,2,23,3:22)]

subway2019_line2_getout$name <- station.info.getout_line2.timediv$name[1:50]
subway2019_line2_getout <- subway2019_line2_getout[,c(1,2,23,3:22)]



subway2020 <- read.csv("Data/seoulmetro/metro2020.csv", 
                       header=TRUE, fileEncoding = "euc-kr")
subway2020 <- subway2020[,c(1,2,4:25)]

colnames(subway2020) <- c("date","line","name_kor","type",paste("t",1:20,sep=""))
subway2020_line2 <- subway2020[subway2020$line=="2호선",]
subway2020_line2[subway2020_line2$name_kor=="동대문역사문화공원(DDP)",]$name_kor <- "동대문역사문화공원"
subway2020_line2[subway2020_line2$name_kor=="왕십리(성동구청)",]$name_kor <- "왕십리"
subway2020_line2[subway2020_line2$name_kor=="구의(광진구청)",]$name_kor <- "구의"
subway2020_line2[subway2020_line2$name_kor=="강변(동서울터미널)",]$name_kor <- "강변"
subway2020_line2[subway2020_line2$name_kor=="잠실(송파구청)",]$name_kor <- "잠실"
subway2020_line2[subway2020_line2$name_kor=="삼성(무역센터)",]$name_kor <- "삼성"
subway2020_line2[subway2020_line2$name_kor=="교대(법원.검찰청)",]$name_kor <- "교대"
subway2020_line2[subway2020_line2$name_kor=="서울대입구(관악구청)",]$name_kor <- "서울대입구"
subway2020_line2[subway2020_line2$name_kor=="대림(구로구청)",]$name_kor <- "대림"
subway2020_line2[subway2020_line2$name_kor=="충정로(경기대입구)",]$name_kor <- "충정로"
subway2020_line2[subway2020_line2$name_kor=="용두(동대문구청)",]$name_kor <- "용두"
subway2020_line2[subway2020_line2$name_kor=="동대문역사문화공원(DDP)(DDP)",]$name_kor <- "동대문역사문화공원"
subway2020_line2[subway2020_line2$name_kor=="낙성대(강감찬)",]$name_kor <- "낙성대"

sort(unique(subway2020_line2$name_kor))

subway2020_line2_getin <- subway2020_line2[subway2020_line2$type=="승차",]
subway2020_line2_getout <- subway2020_line2[subway2020_line2$type=="하차",]

subway2020_line2_getin <- subway2020_line2_getin[,c(1,3,5:24)]
subway2020_line2_getout <- subway2020_line2_getout[,c(1,3,5:24)]

for(i in 1:(nrow(subway2020_line2_getin)/50)){
  subway2020_line2_getin[(50*(i-1)+1):(50*i),] <-
    subway2020_line2_getin[50*(i-1)+order(unique(subway2020_line2$name_kor)),]
}

for(i in 1:(nrow(subway2020_line2_getout)/50)){
  subway2020_line2_getout[(50*(i-1)+1):(50*i),] <-
    subway2020_line2_getout[50*(i-1)+order(unique(subway2020_line2$name_kor)),]
}

subway2020_line2_getin$name <- station.info.getin_line2.timediv$name[1:50]
subway2020_line2_getin <- subway2020_line2_getin[,c(1,2,23,3:22)]

subway2020_line2_getout$name <- station.info.getout_line2.timediv$name[1:50]
subway2020_line2_getout <- subway2020_line2_getout[,c(1,2,23,3:22)]


subway2021 <- read.csv("Data/seoulmetro/metro2021.csv", 
                       header=TRUE, fileEncoding = "euc-kr")
subway2021 <- subway2021[,c(2,3,5:26)]

colnames(subway2021) <- c("date","line","name_kor","type",paste("t",1:20,sep=""))
subway2021_line2 <- subway2021[subway2021$line=="2",]
subway2021_line2[subway2021_line2$name_kor=="동대문역사문화공원(DDP)",]$name_kor <- "동대문역사문화공원"
subway2021_line2[subway2021_line2$name_kor=="왕십리(성동구청)",]$name_kor <- "왕십리"
subway2021_line2[subway2021_line2$name_kor=="구의(광진구청)",]$name_kor <- "구의"
subway2021_line2[subway2021_line2$name_kor=="강변(동서울터미널)",]$name_kor <- "강변"
subway2021_line2[subway2021_line2$name_kor=="잠실(송파구청)",]$name_kor <- "잠실"
subway2021_line2[subway2021_line2$name_kor=="삼성(무역센터)",]$name_kor <- "삼성"
subway2021_line2[subway2021_line2$name_kor=="교대(법원.검찰청)",]$name_kor <- "교대"
subway2021_line2[subway2021_line2$name_kor=="서울대입구(관악구청)",]$name_kor <- "서울대입구"
subway2021_line2[subway2021_line2$name_kor=="대림(구로구청)",]$name_kor <- "대림"
subway2021_line2[subway2021_line2$name_kor=="충정로(경기대입구)",]$name_kor <- "충정로"
subway2021_line2[subway2021_line2$name_kor=="용두(동대문구청)",]$name_kor <- "용두"
# subway2021_line2[subway2021_line2$name_kor=="동대문역사문화공원(DDP)(DDP)",]$name_kor <- "동대문역사문화공원"
subway2021_line2[subway2021_line2$name_kor=="낙성대(강감찬)",]$name_kor <- "낙성대"

sort(unique(subway2021_line2$name_kor))

subway2021_line2_getin <- subway2021_line2[subway2021_line2$type=="승차",]
subway2021_line2_getout <- subway2021_line2[subway2021_line2$type=="하차",]

subway2021_line2_getin <- subway2021_line2_getin[,c(1,3,5:24)]
subway2021_line2_getout <- subway2021_line2_getout[,c(1,3,5:24)]

for(i in 1:(nrow(subway2021_line2_getin)/50)){
  subway2021_line2_getin[(50*(i-1)+1):(50*i),] <-
    subway2021_line2_getin[50*(i-1)+order(unique(subway2021_line2$name_kor)),]
}

for(i in 1:(nrow(subway2021_line2_getout)/50)){
  subway2021_line2_getout[(50*(i-1)+1):(50*i),] <-
    subway2021_line2_getout[50*(i-1)+order(unique(subway2021_line2$name_kor)),]
}

subway2021_line2_getin$name <- station.info.getin_line2.timediv$name[1:50]
subway2021_line2_getin <- subway2021_line2_getin[,c(1,2,23,3:22)]

subway2021_line2_getout$name <- station.info.getout_line2.timediv$name[1:50]
subway2021_line2_getout <- subway2021_line2_getout[,c(1,2,23,3:22)]


subway2022 <- read.csv("Data/seoulmetro/metro2022.csv", 
                       header=TRUE, fileEncoding = "euc-kr")
subway2022 <- subway2022[,c(2,3,5:26)]

colnames(subway2022) <- c("date","line","name_kor","type",paste("t",1:20,sep=""))
subway2022_line2 <- subway2022[subway2022$line=="2",]
subway2022_line2[subway2022_line2$name_kor=="동대문역사문화공원(DDP)",]$name_kor <- "동대문역사문화공원"
subway2022_line2[subway2022_line2$name_kor=="왕십리(성동구청)",]$name_kor <- "왕십리"
subway2022_line2[subway2022_line2$name_kor=="구의(광진구청)",]$name_kor <- "구의"
subway2022_line2[subway2022_line2$name_kor=="강변(동서울터미널)",]$name_kor <- "강변"
subway2022_line2[subway2022_line2$name_kor=="잠실(송파구청)",]$name_kor <- "잠실"
subway2022_line2[subway2022_line2$name_kor=="삼성(무역센터)",]$name_kor <- "삼성"
subway2022_line2[subway2022_line2$name_kor=="교대(법원.검찰청)",]$name_kor <- "교대"
subway2022_line2[subway2022_line2$name_kor=="서울대입구(관악구청)",]$name_kor <- "서울대입구"
subway2022_line2[subway2022_line2$name_kor=="대림(구로구청)",]$name_kor <- "대림"
subway2022_line2[subway2022_line2$name_kor=="충정로(경기대입구)",]$name_kor <- "충정로"
subway2022_line2[subway2022_line2$name_kor=="용두(동대문구청)",]$name_kor <- "용두"
# subway2022_line2[subway2022_line2$name_kor=="동대문역사문화공원(DDP)(DDP)",]$name_kor <- "동대문역사문화공원"
subway2022_line2[subway2022_line2$name_kor=="낙성대(강감찬)",]$name_kor <- "낙성대"

sort(unique(subway2022_line2$name_kor))

subway2022_line2_getin <- subway2022_line2[subway2022_line2$type=="승차",]
subway2022_line2_getout <- subway2022_line2[subway2022_line2$type=="하차",]

subway2022_line2_getin <- subway2022_line2_getin[,c(1,3,5:24)]
subway2022_line2_getout <- subway2022_line2_getout[,c(1,3,5:24)]

for(i in 1:(nrow(subway2022_line2_getin)/50)){
  subway2022_line2_getin[(50*(i-1)+1):(50*i),] <-
    subway2022_line2_getin[50*(i-1)+order(unique(subway2022_line2$name_kor)),]
}

for(i in 1:(nrow(subway2022_line2_getout)/50)){
  subway2022_line2_getout[(50*(i-1)+1):(50*i),] <-
    subway2022_line2_getout[50*(i-1)+order(unique(subway2022_line2$name_kor)),]
}

subway2022_line2_getin$name <- station.info.getin_line2.timediv$name[1:50]
subway2022_line2_getin <- subway2022_line2_getin[,c(1,2,23,3:22)]

subway2022_line2_getout$name <- station.info.getout_line2.timediv$name[1:50]
subway2022_line2_getout <- subway2022_line2_getout[,c(1,2,23,3:22)]



subway2023 <- read.csv("Data/seoulmetro/metro2023.csv", 
                       header=TRUE, fileEncoding = "euc-kr")
subway2023 <- subway2023[,c(2,3,5:26)]

colnames(subway2023) <- c("date","line","name_kor","type",paste("t",1:20,sep=""))
subway2023_line2 <- subway2023[subway2023$line=="2호선",]
subway2023_line2[subway2023_line2$name_kor=="동대문역사문화공원(DDP)",]$name_kor <- "동대문역사문화공원"
subway2023_line2[subway2023_line2$name_kor=="왕십리(성동구청)",]$name_kor <- "왕십리"
subway2023_line2[subway2023_line2$name_kor=="구의(광진구청)",]$name_kor <- "구의"
subway2023_line2[subway2023_line2$name_kor=="강변(동서울터미널)",]$name_kor <- "강변"
subway2023_line2[subway2023_line2$name_kor=="잠실(송파구청)",]$name_kor <- "잠실"
subway2023_line2[subway2023_line2$name_kor=="삼성(무역센터)",]$name_kor <- "삼성"
subway2023_line2[subway2023_line2$name_kor=="교대(법원.검찰청)",]$name_kor <- "교대"
subway2023_line2[subway2023_line2$name_kor=="서울대입구(관악구청)",]$name_kor <- "서울대입구"
subway2023_line2[subway2023_line2$name_kor=="대림(구로구청)",]$name_kor <- "대림"
subway2023_line2[subway2023_line2$name_kor=="충정로(경기대입구)",]$name_kor <- "충정로"
subway2023_line2[subway2023_line2$name_kor=="용두(동대문구청)",]$name_kor <- "용두"
# subway2023_line2[subway2023_line2$name_kor=="동대문역사문화공원(DDP)(DDP)",]$name_kor <- "동대문역사문화공원"
subway2023_line2[subway2023_line2$name_kor=="낙성대(강감찬)",]$name_kor <- "낙성대"

sort(unique(subway2023_line2$name_kor))

subway2023_line2_getin <- subway2023_line2[subway2023_line2$type=="승차",]
subway2023_line2_getout <- subway2023_line2[subway2023_line2$type=="하차",]

subway2023_line2_getin <- subway2023_line2_getin[,c(1,3,5:24)]
subway2023_line2_getout <- subway2023_line2_getout[,c(1,3,5:24)]

for(i in 1:(nrow(subway2023_line2_getin)/50)){
  subway2023_line2_getin[(50*(i-1)+1):(50*i),] <-
    subway2023_line2_getin[50*(i-1)+order(unique(subway2023_line2$name_kor)),]
}

for(i in 1:(nrow(subway2023_line2_getout)/50)){
  subway2023_line2_getout[(50*(i-1)+1):(50*i),] <-
    subway2023_line2_getout[50*(i-1)+order(unique(subway2023_line2$name_kor)),]
}

subway2023_line2_getin$name <- station.info.getin_line2.timediv$name[1:50]
subway2023_line2_getin <- subway2023_line2_getin[,c(1,2,23,3:22)]

subway2023_line2_getout$name <- station.info.getout_line2.timediv$name[1:50]
subway2023_line2_getout <- subway2023_line2_getout[,c(1,2,23,3:22)]


subway_line2_getin <- rbind(subway2019_line2_getin[,1:22],subway2020_line2_getin[,1:22],subway2021_line2_getin[,1:22],
                            subway2022_line2_getin[,1:22],subway2023_line2_getin[,1:22])
subway_line2_getout <- rbind(subway2019_line2_getin[,1:22],subway2020_line2_getout[,1:22],subway2021_line2_getout[,1:22],
                             subway2022_line2_getout[,1:22],subway2023_line2_getout[,1:22])

subway_line2_total <- subway_line2_getin
subway_line2_total[,4:22] <- subway_line2_getin[,4:22] + subway_line2_getout[,4:22]
subway_line2_total <- subway_line2_total[(50*6+1):(nrow(subway_line2_total)-50*303),]


# covid vs post-covid
# weekly data
p.subway <- 19
R.subway <- length(unique(subway_line2_total$date))/7 # 208 weeks

subway_line2_total_wk <- NULL
for(i in 1:R.subway){
  subway_line2_total_wk_tmp <- NULL
  for(j in 1:7){
    if(j==1){
      subway_line2_total_wk_tmp <- subway_line2_total[(50*(7*(i-1)+j-1)+1:50),4:22]
    } else{
      subway_line2_total_wk_tmp <- subway_line2_total_wk_tmp +
        subway_line2_total[(50*(7*(i-1)+j-1)+1:50),4:22]
    }
  }
  subway_line2_total_wk_tmp <- subway_line2_total_wk_tmp / 7
  if(i==1){
    subway_line2_total_wk <- subway_line2_total_wk_tmp
  } else{
    subway_line2_total_wk <- rbind(subway_line2_total_wk, subway_line2_total_wk_tmp)
  }
}


subway_line2_total_wk$weeknum <- rep(1:R.subway, each=50)
subway_line2_total_wk$name_kor <- subway_line2_total$name_kor[1:nrow(subway_line2_total_wk)]
subway_line2_total_wk$name <- subway_line2_total$name[1:nrow(subway_line2_total_wk)]

subway_line2_total_wk <- subway_line2_total_wk[,c(20:22,1:19)]
head(subway_line2_total_wk)

meansig_subway_log <- as.data.frame(subway_line2_total_wk %>% group_by(name_kor,name) %>% 
                                      summarise(meanlog_t1 = mean(log(1+t1)),meanlog_t2 = mean(log(1+t2)),meanlog_t3 = mean(log(1+t3)),
                                                meanlog_t4 = mean(log(1+t4)),meanlog_t5 = mean(log(1+t5)),meanlog_t6 = mean(log(1+t6)),
                                                meanlog_t7 = mean(log(1+t7)),meanlog_t8 = mean(log(1+t8)),meanlog_t9 = mean(log(1+t9)),
                                                meanlog_t10 = mean(log(1+t10)),meanlog_t11 = mean(log(1+t11)),meanlog_t12 = mean(log(1+t12)),
                                                meanlog_t13 = mean(log(1+t13)),meanlog_t14= mean(log(1+t14)),meanlog_t15 = mean(log(1+t15)),
                                                meanlog_t16 = mean(log(1+t16)),meanlog_t17 = mean(log(1+t17)),meanlog_t18 = mean(log(1+t18)),
                                                meanlog_t19 = mean(log(1+t19))) %>%
                                      arrange(name_kor))

head(meansig_subway_log,10)


X.subway_line2_total_wk <- array(0, c(p.subway,N.metro_line2,R.subway))

for(i in 1:p.subway){
  for(r in 1:R.subway){
    data.subway_total_wk <- subway_line2_total_wk[subway_line2_total_wk$weeknum == r, (i+3)]
    data.subway_total_wk <- log(1+data.subway_total_wk)
    data.subway_total_wk <- data.subway_total_wk - meansig_subway_log[,2+i]
    
    X.subway_line2_total_wk[i,,r] <- data.subway_total_wk
  }
}


# Fourier coefficients
X.tilde.subway.total.line2_wk <- array(0, c(p.subway,N.metro_line2,R.subway))

for(i in 1:p.subway){
  X.tilde.subway.total.line2_wk[i,,] <- Conj(t(evectors.metro_line2)) %*% X.subway_line2_total_wk[i,,]
}

# reordering
X.ut.subway.total.line2_wk <- array(0,c(N.metro_line2,R.subway,p.subway)) # n x R x p

for(l in 1:N.metro_line2){
  X.ut.subway.total.line2_wk[l,,] <- t(X.tilde.subway.total.line2_wk[,l,])
}



# resultant estimates
k.optimal.subway.total_wk <- 2 # inferred from the scree plot below


f.bar.subway.total_wk <- array(0, c(N.metro_line2,R.subway,k.optimal.subway.total_wk)) # n x R x k
B.bar.subway.total_wk <- array(0, c(N.metro_line2,p.subway,k.optimal.subway.total_wk)) # n x p x k

eigenprops <- c()
for(l in 1:N.metro_line2){
  eigenres <- eigen(X.ut.subway.total.line2_wk[l,,] %*% Conj(t(X.ut.subway.total.line2_wk[l,,])))
  evectors <- eigenres$vectors
  f.bar.subway.total_wk[l,,] <- sqrt(R.subway)*evectors[,1:k.optimal.subway.total_wk]
  B.bar.subway.total_wk[l,,] <- t(X.ut.subway.total.line2_wk[l,,]) %*% Conj(f.bar.subway.total_wk[l,,]) / R.subway
  
  eigenprops <- rbind(eigenprops, cumsum(eigenres$values / sum(eigenres$values))[1:10])
}


# invert ordering
f.tilde.subway.total_wk <- array(0, c(k.optimal.subway.total_wk,N.metro_line2,R.subway)) # k x n x R
B.subway.total_wk <- array(0, c(p.subway,k.optimal.subway.total_wk,N.metro_line2))

for(l in 1:N.metro_line2){
  f.tilde.subway.total_wk[,l,] <- t(f.bar.subway.total_wk[l,,])
  
  B.subway.total_wk[,,l] <- B.bar.subway.total_wk[l,,]
}

# iGFT
f.subway.total_wk <- array(0, c(k.optimal.subway.total_wk,N.metro_line2,R.subway)) # k x n x R
for(j in 1:k.optimal.subway.total_wk){
  f.subway.total_wk[j,,] <- evectors.metro_line2 %*% f.tilde.subway.total_wk[j,,]
}



## clustering
kmeans_res_total_wk <- kmeans(t(f.subway.total_wk[1,,]), centers = 2, nstart = 25)
kmeans_res_total_wk$cluster
table(kmeans_res_total_wk$cluster)

layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))  # 4:1 width ratio for plot and legend
# Plot area
par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot
plot(kmeans_res_total_wk$cluster, col = c(rep("blue",51), rep("green",53), rep("magenta",54), rep("red",53)), 
     pch = 16, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "cluster", cex.lab=1.4)

axis(2, at = c(1, 2), labels = c("1", "2"))
# Add vertical lines for year separation
abline(v = 0.5, lty = 2)
abline(v = 52.5, lty = 2)
abline(v = 104.5, lty = 2)
abline(v = 156.5, lty = 2)
abline(v = 208.5, lty = 2)

text_positions <- c(26.5, 78.5, 130.5, 182.5)  # Middle positions between lines
mtext(text = 2019:2022, side = 1, at = text_positions, line = 1)
par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", legend = c("2019", "2020", "2021", "2022"),
       col = c("blue", "green", "magenta", "red"), pch = 16, cex = 1.2, bty = "n")


# use original signals
for(ii in 1:p.subway){
  kmeans_res_total <- kmeans(t(X.subway_line2_total_wk[ii,,]), centers = 2, nstart = 40)

  print(ii)
  print(kmeans_res_total$cluster)
}

# use averaged signal
tmppp <- t(X.subway_line2_total_wk[1,,])
for(ii in 2:p.subway){
  tmppp <- tmppp + t(X.subway_line2_total_wk[ii,,])
}

kmeans_res_total <- kmeans(tmppp/19, centers = 2, nstart = 40)
kmeans_res_total$cluster
table(kmeans_res_total$cluster)


## loading plot for factor 1
layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))  # 4:1 width ratio for plot and legend
par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot

matplot(B.bar.subway.total_wk[,,1], type="l", lty=1, col=colorRampPalette(c("red","grey", "blue"))(19),
        xlab="graph frequency index", ylab="loading", cex.lab=1.4, lwd=1.4)
abline(h=0, lty=2)

par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", inset=c(-0.4,0), legend=c("before 6 a.m.","6 a.m. - 7 a.m.","7 a.m. - 8 a.m.","8 a.m. - 9 a.m.","9 a.m. - 10 a.m.",
                                           "10 a.m. - 11 a.m.","11 a.m. - 12 p.m.","12 p.m. - 1 p.m.","1 p.m. - 2 p.m.","2 p.m. - 3 p.m.",
                                           "3 p.m. - 4 p.m.","4 p.m. - 5 p.m.","5 p.m. - 6 p.m.","6 p.m. - 7 p.m.","7 p.m. - 8 p.m.",
                                           "8 p.m. - 9 p.m.","9 p.m. - 10 p.m.","10 p.m. - 11 p.m.","after 11 p.m."),
       col=colorRampPalette(c("red", "grey","blue"))(19), lty=1, lwd=1.4, bty = "n", cex=1.1)


## loading plot for factor 2
layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))  # 4:1 width ratio for plot and legend
par(mar = c(5, 4, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot

matplot(B.bar.subway.total_wk[,,2], type="l", lty=1, col=colorRampPalette(c("red", "yellow", "green", "blue", "purple"))(19),
        xlab="graph frequency index", ylab="loading", cex.lab=1.4, lwd=1.4)
abline(h=0, lty=2)

par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", inset=c(-0.4,0), legend=c("before 6 a.m.","6 a.m. - 7 a.m.","7 a.m. - 8 a.m.","8 a.m. - 9 a.m.","9 a.m. - 10 a.m.",
                                           "10 a.m. - 11 a.m.","11 a.m. - 12 p.m.","12 p.m. - 1 p.m.","1 p.m. - 2 p.m.","2 p.m. - 3 p.m.",
                                           "3 p.m. - 4 p.m.","4 p.m. - 5 p.m.","5 p.m. - 6 p.m.","6 p.m. - 7 p.m.","7 p.m. - 8 p.m.",
                                           "8 p.m. - 9 p.m.","9 p.m. - 10 p.m.","10 p.m. - 11 p.m.","after 11 p.m."),
       col=colorRampPalette(c("red","yellow", "green", "blue", "purple"))(19), lty=1, lwd=1.4, bty = "n", cex=1.1)



## Arrow plot
# Set up a 2x3 grid for the plots and one column for the legend
layout(matrix(c(1, 2, 3, 7, 
                4, 5, 6, 7), nrow = 2, byrow = TRUE), widths = c(2, 2, 2, 1))


par(mar = c(5, 6, 4, 1) + 0.1, xpd=FALSE, oma=c(0,1,0,0))  # Reset margins for the main plot

for(l in c(1,10,20,30,40,50)){
  plot(0, 0, type="n", xlim=c(0, max(abs(B.bar.subway.total_wk[l,,1]))), 
       ylim=c(0, max(abs(B.bar.subway.total_wk[l,,2]))), 
       xlab="loading 1", 
       ylab="loading 2", cex.lab=1.6, lwd=1.2, main=bquote(lambda[.(l)]^" metro"), cex.main=1.7)
  
  # Draw arrows from the origin to each point
  for(i in 1:19) {
    arrows(0, 0, abs(B.bar.subway.total_wk[l,i,1]), abs(B.bar.subway.total_wk[l,i,2]), 
           col=colorRampPalette(c("red","yellow", "green", "blue", "purple"))(19)[i], lwd=2, length=0.1)  # length controls arrowhead size
    
  }
}

par(mar = c(0, 0, 0, 0))  # Minimize margins in the legend area
plot.new()
legend("center", inset=c(-0.4,0), legend=c("before 6 a.m.","6 a.m. - 7 a.m.","7 a.m. - 8 a.m.","8 a.m. - 9 a.m.","9 a.m. - 10 a.m.",
                                           "10 a.m. - 11 a.m.","11 a.m. - 12 p.m.","12 p.m. - 1 p.m.","1 p.m. - 2 p.m.","2 p.m. - 3 p.m.",
                                           "3 p.m. - 4 p.m.","4 p.m. - 5 p.m.","5 p.m. - 6 p.m.","6 p.m. - 7 p.m.","7 p.m. - 8 p.m.",
                                           "8 p.m. - 9 p.m.","9 p.m. - 10 p.m.","10 p.m. - 11 p.m.","after 11 p.m."),
       col=colorRampPalette(c("red","yellow", "green", "blue", "purple"))(19), lty=1, lwd=2, bty = "n", cex=1.3)



## draw plot
g1 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = X.subway_line2_total_wk[1,,1], value="Value", ratio=0.6,
                         min=-0.6, max=1.05, mg=c(4,4,4,4), title="before 6 a.m.", main.title.size = 20)
g2 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = X.subway_line2_total_wk[4,,1], value="Value", ratio=0.6,
                         min=-1.8, max=0.6, mg=c(4,4,4,4), title="8 a.m. - 9 a.m.", main.title.size = 20)
g3 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = X.subway_line2_total_wk[11,,1], value="Value", ratio=0.6,
                         min=-0.05, max=0.55, mg=c(4,4,4,4), title="3 p.m. - 4 p.m.", main.title.size = 20)
g4 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = X.subway_line2_total_wk[17,,1], value="Value", ratio=0.6,
                         min=-0.9, max=0.95, mg=c(4,4,4,4), title="9 p.m. - 10 p.m.", main.title.size = 20)
g5 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = X.subway_line2_total_wk[19,,1], value="Value", ratio=0.6,
                         min=-1.1, max=1.25, mg=c(4,4,4,4), title="after 11 p.m.", main.title.size = 20)
g6 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = f.subway.total_wk[1,,1], value="Value", ratio=0.6,
                         min=-4.35, max=3.7, mg=c(4,4,4,4), title=expression(F[1]^"metro"), main.title.size = 20)
g7 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = f.subway.total_wk[2,,1], value="Value", ratio=0.6,
                         min=-4.05, max=1.8, mg=c(4,4,4,4), title=expression(F[2]^"metro"), main.title.size = 20)
g8 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = evectors.metro_line2[,1], value="Value", ratio=0.6,
                         min=0.1, max=0.15, mg=c(4,4,4,4), title=expression(v[1]^"metro"), main.title.size = 20)
g9 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = evectors.metro_line2[,8], value="Value", ratio=0.6,
                         min=-0.45, max=0.25, mg=c(4,4,4,4), title=expression(v[8]^"metro"), main.title.size = 20)

# g6 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = evectors.metro_line2[,1], value="Value", ratio=0.6,
#                          min=0.1, max=0.15, mg=c(4,4,4,4), title="v1", main.title.size = 20)
grid.arrange(g1,g2,g4,g5,g6,g7,g8,g9, nrow=3)


g10 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = X.subway_line2_total_wk[1,,104], value="Value", ratio=0.6,
                          min=-1.1, max=-0.05, mg=c(4,4,4,4), title="before 6 a.m.", main.title.size = 20)
g11 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = X.subway_line2_total_wk[4,,104], value="Value", ratio=0.6,
                          min=-0.7, max=0.2, mg=c(4,4,4,4), title="8 a.m. - 9 a.m.", main.title.size = 20)
g12 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = X.subway_line2_total_wk[11,,104], value="Value", ratio=0.6,
                          min=-0.9, max=-0.15, mg=c(4,4,4,4), title="3 p.m. - 4 p.m.", main.title.size = 20)
g13 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = X.subway_line2_total_wk[17,,104], value="Value", ratio=0.6,
                          min=-1.1, max=0, mg=c(4,4,4,4), title="9 p.m. - 10 p.m.", main.title.size = 20)
g14 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = X.subway_line2_total_wk[19,,104], value="Value", ratio=0.6,
                          min=-1.95, max=-0.6, mg=c(4,4,4,4), title="after 11 p.m.", main.title.size = 20)
g15 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = f.subway.total_wk[1,,104], value="Value", ratio=0.6,
                          min=-1.9, max=2.25, mg=c(4,4,4,4), title=expression(F[1]^"metro"), main.title.size = 20)
g16 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=8, vertex_color = f.subway.total_wk[2,,104], value="Value", ratio=0.6,
                          min=-3.15, max=2.05, mg=c(4,4,4,4), title=expression(F[2]^"metro"), main.title.size = 20)

grid.arrange(g10,g11,g13,g14,g15,g16,g8,g9, nrow=3)


g8 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = evectors.metro_line2[,1], value="Value", ratio=0.6,
                         min=0.1, max=0.15, mg=c(4,4,4,4), title=expression(v[1]^"metro"), main.title.size = 20)
g9 <- plot_graph_custom4(seoulmetro_line2, e.size=1.3, v.size=6, vertex_color = evectors.metro_line2[,8], value="Value", ratio=0.6,
                         min=-0.45, max=0.25, mg=c(4,4,4,4), title=expression(v[1]^"metro"), main.title.size = 20)


grid.arrange(g1,g2,g3,g4,g8,g9, nrow=2)
grid.arrange(g5,g6,g7,g4,g8,g9, nrow=2)

save(seoulmetro_line2, file="seoulmetro_line2.RData")


