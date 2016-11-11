#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

# define the first function

Rest_Distri <-function(file_data, output) {

    data_in <- read.table(file_data, sep ="\t", header = TRUE)
    colnames(data_in) <- c("Fragment_Length", "In Silico", "Experiment")
    data_t <- melt(data_in, id = "Fragment_Length")
    data_t$variable <- as.factor(data_t$variable)
    p <- ggplot(data_t, aes(x=Fragment_Length, y=value,group = variable, fill=variable)) + 
    geom_bar(aes(color=variable),stat="identity", position="dodge") + 
    xlab("Fragment Length(bp)")+ylab("Number of Tags") + ggtitle("Test_Graph") + 
    labs(color = "Group", fill = "Group")
    outfile <- paste(output,"png",sep=".")
    ggsave(outfile, p, width=10,height=8, dpi=120) 
 
    }

#define the second function
Fr_Number_Distribution <- function(file_a, file_b, Out){
    
    data1 <- read.table(file_a, sep="\t", header=FALSE)
    colnames(data1) <- c("Length_range", "Num")
    data2 <- read.table(file_b, sep="\t", skip =1, header=FALSE, stringsAsFactors = FALSE)
    colnames(data2) <- c("Length_range", "Exp_Num")
    data_Tmp <- merge(data1, data2, by="Length_range", all=T)
    raw_data <- melt(data_Tmp, id = "Length_range")
    raw_data[,1] <- as.character(raw_data[,1])
    new_data <- raw_data
    new_data[,1] <- do.call(rbind, strsplit(new_data[,1], "-", fixed = TRUE))
    a <- data.frame(new_data[,1])
    new_data[,1] <- as.numeric(as.character(a[,1]))
    lab <- as.character(data1[,1])
    p <- ggplot(new_data, aes(x=Length_range, y = value, group = variable, fill = variable)) + geom_bar(aes(color = variable), stat= "identity", position = "dodge") + scale_x_continuous(breaks=seq(0,3800,200),labels= lab[seq(1, length(lab), 10)]) + theme(axis.text.x = element_text(face="bold", color="grey", size=5)) +
    labs(color = "Group", fill = "Group")
    print (Out)
    
    filename <- paste(Out,".png", sep="")
    ggsave(filename, p, width =10, height = 8, dpi=300)

    }

#define the third function
Skip_Firs_Line_Plot <- function(file_data, out) {
    src_data <- read.table(file_data, skip=1, header=FALSE, sep="\t", stringsAsFactors = FALSE)
    data_t <- melt(src_data, id="V1")

    New_Data <- data_t
    New_Data[,1] <- do.call(rbind, strsplit(New_Data[,1], '-', fixed = TRUE))
    a <- data.frame(New_Data[,1])
    lab <- as.character(data_t$V1)
    New_Data[,1] <- as.numeric(as.character(a[,1]))
    p <- ggplot(New_Data, aes(x=V1, y=value, group =variable,fill=variable)) + 
    geom_bar(aes(color=variable),stat="identity", position="dodge")+
    xlab("Fragment Length(bp)")+ylab("Number of Tags") + ggtitle("Test_Graph") + 
    labs(color ="Group", fill = "Group") +
    scale_x_continuous(breaks=seq(0,280,20), labels = lab) +
    scale_y_continuous(breaks=seq(0,1400000,100000), labels = seq(0,1400000,100000)/1000) +
    theme_bw() + theme(axis.text.x = element_text(face="bold", color="grey", size=6), plot.margin = unit(c(2,3,3,4),"cm"),panel.margin = unit(c(1,1,1,1), "cm") ,panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    outfile <- paste(out,"png",sep=".")
    ggsave(outfile, p, width=10,height=8)     
    }


args   = commandArgs(T)
f1      = args[1]

if(length(args) == 3){

    f2 = args[2]
    OutPic = args[3]
    
    Fr_Number_Distribution(f1,f2, OutPic)
} else if (is.na(args[2])) {
    OutPic <- "Test_pic"
    Skip_Firs_Line_Plot(f1, OutPic)
    #Rest_Distri(f1, OutPic)
} else {
    OutPic <- args[2]
    Skip_Firs_Line_Plot(f1, OutPic)
    #Rest_Distri(f1, OutPic)
}


