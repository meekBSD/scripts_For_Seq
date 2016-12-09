#!/usr/bin/env Rscript

if (!require("RODBC")) {
    install.packages("RODBC", dependencies = T, repos= "http://cran.us.r-project.org")
   }

library(RODBC)
args = commandArgs(T)

f_name = args[1]

a = strsplit(f_name, split="\\.")
len = length(a[[1]])

suffix = a[[1]][len]

if (suffix == "xlsx"){
    con.tab <- odbcConnectExcel2007(f_name)
    sheet.names <- sqlTables(con.tab)
    print (sheet.names$TABLE_NAME)
    odbcClose(con.tab)
} else if (suffix == "xls"){
    con.tab <- odbcConnect(f_name)   # check this function "odbcConnect" whether it exists in RODBC
    sheet.names <- sqlTables(con.tab)

    n<-length(sheet.names$TABLE_NAME)
    s.name=gsub('\\$','', sheet.names$TABLE_NAME)
    sname=gsub('\'','', s.name)
    sheets1 <- list()
    for(x in 1:n) {
        sheets1[[x]] <- sqlFetch(con, sname[x],colnames=FALSE,rownames=FALSE)
        sheets1[[x]] <- sheets1[[x]][,c(1,2,4)]
        write.table(sheets1[[x]], file = paste(x,"_result_data.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE)
    } 
    odbcClose(con.tab)
}
    
# use this script like this:
# Rscript test.R goldData.xls 


