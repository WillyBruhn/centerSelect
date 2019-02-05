# wd = "/home/sysgen/Documents/LWB/EvaluateCluster/"
wd = "/home/willy/RedoxChallenges/centerSelect"

labels_file = "labels.csv"

setwd(wd)

labels = read.csv(labels_file, header = FALSE, na.strings=c("","NA"))

tab = labels[1:103,5:7]
names(tab) = c("protA", "label", "protB")

tab[2,]

which(is.na(tab$label))
is.na(tab$label)

levels(tab$protB)

tab$protB = replace(tab$protB, which("trx" == tab$protB), 'Trx')
tab$label = !is.na(tab$label)


tab


path_to_files = "/home/willy/RedoxChallenges/centerSelectTest/"
distances_file_neg = paste(path_to_files, "ListOfEMD_neg_100.csv", sep = "")
distances_file_pos = paste(path_to_files, "ListOfEMD_pos_100.csv", sep = "")

distances_neg = read.csv(distances_file_neg)
distances_pos = read.csv(distances_file_pos)










