# graphFA

This repo contains the source code for graph frequency domain factor analysis of multivariate graph signals.

## Description

- Code
  - `utils.R` is a code for functions used for graph frequency domain factor analysis.
  - `source.R` is a code for functions used for river network data analysis.
  - `simulation_karate.R`, `simulation_ush.R` are codes for simulation study.
  - `seoul_metro.R`, `river_network.R` are codes for real data analysis.


- Data
  - `seoulmetro` contains data of daily number of people getting on and off the Seoul Metropolitan Subway in South Korea.
  - `stationary` contains data of hourly temperature measurements recorded in Fahrenheit across the United States on August 1, 2010.
  - `KRF_3.0_Geumgang` contains shape files related to the Geum River such as catchment area shape, line shape, nodes. They can be downloaded from [Korean Reach File](http://water.nier.go.kr/web/gisKrf?pMENU_NO=89).
  - `ProcessedData` contains processed TOC data for Geum River. The original TOC data can be obtained from [Water Environment Information System](http://water.nier.go.kr/web/waterMeasure?pMENU_NO=2). 

- stpca_Rpackage
  - R package from "Flow-directed PCA (Gallacher et al., 2017)". 
  
## Code overview
We present a graph frequency domain factor model.

