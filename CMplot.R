library(readr)
library(tidyverse)
library(readxl)
Eastern_Bedik <- read_excel("C:/Users/M/Downloads/zasilka-WO476P3GXF33Y4H6/Eastern_Bedik_XPEHH.xlsx", 
                            col_names = TRUE)
Eastern_Chad <- read_excel("C:/Users/M/Downloads/zasilka-WO476P3GXF33Y4H6/Eastern_Chad_XPEHH.xlsx", 
                           col_names = TRUE)
Eastern_Yemen <- read_excel("C:/Users/M/Downloads/zasilka-WO476P3GXF33Y4H6/Eastern_Yemen_XPEHH.xlsx", 
                            col_names = TRUE)

data <- merge(Eastern_Bedik, Eastern_Chad, Eastern_Yemen,  by=c("CHR", "POSITION"), all.x=T, all.y=T, sort = F)
colnames(data) <- c("CHR", "POS","Bedik", "Chad", "Yemen")


library(CMplot)
p<-CMplot(abs(data),type="p",plot.type="c",LOG10=FALSE, r=1,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:5),sep=""),
          cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),
          ylab=expression(abs(p)), cir.band = 0.3, H = 0.5, 
          threshold.col=c("red","blue"),signal.line=1,
          signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
          bin.size=1e6,outward=TRUE,file="tiff",memo="",dpi=600,file.output=F,verbose=F,width=10,height=10, main = "Selection scans")
