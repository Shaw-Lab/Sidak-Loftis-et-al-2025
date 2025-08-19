
### Kaylee Vosbigian
### 08.09.25

# :::::::::::::::::::::::::::::::::::::: #
# TRAF2 binding motif search #
# ::::::::::::::::::::::::::::::::::::: #

#read in packages
library(seqinr)
library(stringr)

#select proteome file
mSeq<-read.fasta(file = "FILE NAME", seqtype = "AA")


#create empty list for results
output<-list()
count<-1

#detect protein motif in amino acid sequences
for (i in 1:length(mSeq)){
  s1<-paste0(mSeq[[i]],collapse = "")
  
  if (str_detect(s1,"[PSAT]...[QE]E")){
    start<-1
    myStr<-""
    res<-as.data.frame(str_locate_all(s1,"[PSAT]...[QE]E"))
    for (j in 1:nrow(res)) {
      end<-res[j,"start"]
      myStr<-paste0(myStr,substr(s1,start,end))
      myStr<-paste0(myStr,"*")
      start<-end
    }
    myStr<-paste0(myStr,substr(s1,start,nchar(s1)))
    attributes(myStr)<-attributes(mSeq[[i]])
    output[[count]]<-myStr
    count<-count+1
  }
}

#format data to include annotations
#create empty data frame with labeled columns
myOut<-data.frame(seq = character(nrows), name = character(nrows), Annot = character(nrows), gene = character(nrows),protein = character(nrows))

#add output data into data frame with annotations
for (i in 1:length(output)){
  myOut[i,"seq"]<-output[[i]]
  myOut[i,"name"]<-attributes(output[[i]])$name
  myOut[i,"Annot"]<-attributes(output[[i]])$Annot
  myOut[i,"gene"]<-str_extract( attributes(output[[i]])$Annot[1], "(?<=gene=)[^=]*(?=])")
  myOut[i,"protein"]<-str_extract( attributes(output[[i]])$Annot[1], "(?<=protein=)[^=]*(?=])")
}

# save the file as csv
write.csv(myOut,file ="TRAF2_CS_set_major_09_2022.csv") #replace file name in quotes with your own file name 
