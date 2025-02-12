#!/usr/bin/Rscript
library(ggplot2)


infilename=commandArgs()[5]
Max_len=commandArgs()[6]
Min_len=commandArgs()[7]
data=read.table(infilename)


color=c("darkorchid4", "seagreen", "orange3")
data[,4]=apply(data,1,function(x){b=as.numeric(x[1]);if(b%%3==0){return(0)}else if(b%%3==1){return(1)}else{return(2)}})
data_subset=data[data[,2]>= as.numeric(Min_len) & data[,2]<= as.numeric(Max_len),]
pdf_name=sub(".txt",paste0("_", Max_len, "_", Min_len, ".pdf"),infilename)
pdf(pdf_name,width=10,height=10)
ggplot(data_subset)+geom_bar(aes(x=V1+0.4,y=V3,fill=factor(V4)),stat = "identity")+facet_wrap(~V2,scales="free",ncol=2)+coord_cartesian(xlim=c(-96,6))+scale_fill_manual(name="Class",values=color)+theme_bw()+ylab("Read counts")+xlab("Distance to start codon")+theme(legend.position="non")+scale_x_continuous(breaks=seq(-96,6,by=6))
#ggplot(data_subset)+geom_line(aes(x=V1,y=V3,colour=factor(V4),group = 1),stat = "identity")+facet_wrap(~V2,scales="free",ncol=2)+coord_cartesian(xlim=c(-21,150))+theme_bw()+ylab("Read counts")+xlab("Distance to start codon")+scale_colour_manual(name="Class",values=color)+theme(legend.position="non")+scale_x_continuous(breaks=seq(-10,150,by=10))
data_subset_aggregate=aggregate(x=.~ V2+V4,data=data_subset,FUN=sum)

ggplot(data_subset_aggregate)+geom_bar(aes(x=V2,y=V3,fill=factor(V4,levels=c(0,1,2))),stat="identity")+theme_bw()+ylab("Read counts")+xlab("Read length")+theme(panel.grid.minor=element_blank())+scale_fill_manual(name="frame 5'end",values=color,labels=c("frame 0","frame 1","frame 2"))+scale_x_continuous(breaks = round(seq(min(data[,2]), max(data[,2]), by = 2)))
garbage <- dev.off()


