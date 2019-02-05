setwd("/home/willy/RedoxChallenges/centerSelect/AA_utility")

prot_file = read.table("023HeadOff.pqr")


AAs = prot_file[,3:4]
AAs

levels(AAs[,2])
which(AAs$V4 == "ALA")
ALA = AAs[419:428,1]

which(AAs$V4 == "ARG")
ARG = AAs[100:123,1]

which(AAs$V4 == "ASN")
ASN = AAs[1215:1228,1]

which(AAs$V4 == "ASP")
ASP = AAs[1:14,1]

which(AAs$V4 == "CYS")
CYS = AAs[148:158,1]

which(AAs$V4 == "GLN")
GLN = AAs[124:140,1]

which(AAs$V4 == "GLU")
GLU = AAs[200:214,1]

which(AAs$V4 == "GLY")
GLY = AAs[141:147,1]

which(AAs$V4 == "ILE")
ILE = AAs[570:588,1]

which(AAs$V4 == "LEU")
LEU = AAs[15:33,1]

which(AAs$V4 == "MET")
MET = AAs[509:525,1]

which(AAs$V4 == "PHE")
PHE = AAs[69:88,1]

which(AAs$V4 == "PRO")
PRO = AAs[612:625,1]

which(AAs$V4 == "SER")
SER = AAs[89:99,1]

which(AAs$V4 == "THR")
THR = AAs[626:639,1]

which(AAs$V4 == "TRP")
TRP = AAs[1229:1252,1]

which(AAs$V4 == "TYR")
TYR = AAs[768:788,1]

which(AAs$V4 == "VAL")
VAL = AAs[34:49,1]


prot_file = read.table("009HeadOff.pqr")
AAs = prot_file[,3:4]
AAs

which(AAs$V4 == "LYS")
LYS = AAs[191:212,1]

which(AAs$V4 == "HIS")
HIS = AAs[969:985,1]

# prot_file = read.table("008HeadOff.pqr")
# AAs = prot_file[,3:4]
# length(levels(AAs$V4))
# which(AAs$V4 == "NA")

na.pad <- function(x,len){
  x[1:len]
}

makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}

df = makePaddedDataFrame(list(ALA=ALA,ARG=ARG,ASN=ASN,ASP=ASP,CYS=CYS,GLN=GLN,GLU=GLU,GLY=GLY,ILE=ILE,LEU=LEU,MET=MET,PHE=PHE,PRO=PRO,SER=SER,THR=THR,TRP=TRP,TYR=TYR,VAL=VAL,LYS=LYS, HIS=HIS))

# View(df)


write.csv2(df,"AAs.csv", row.names = FALSE)


df2 = read.csv2("AAs.csv")


get_length_of_AA <- function(AA_column){
  return(sum(!(is.na(AA_column)), na.rm = TRUE))
}

get_length_of_AA(df$ALA)


makeAALengthDataFrame <- function(l,...){
  print(sapply(l,length))
    #maxlen <- max(sapply(l,length))
  #data.frame(lapply(l,na.pad,len=maxlen),...)
  
  data.frame(t(sapply(l,length)))
}

AA_lengths = makeAALengthDataFrame(list(ALA=ALA,ARG=ARG,ASN=ASN,ASP=ASP,CYS=CYS,GLN=GLN,GLU=GLU,GLY=GLY,ILE=ILE,LEU=LEU,MET=MET,PHE=PHE,PRO=PRO,SER=SER,THR=THR,TRP=TRP,TYR=TYR,VAL=VAL, LYS=LYS, HIS=HIS))

write.csv2(AA_lengths,"AA_lengths.csv", row.names = FALSE)

df3 = read.csv2("AA_lengths.csv")

# View(df3)

