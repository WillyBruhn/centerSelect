#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=7) {
  stop("need 7 aruments", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  #args[2] = "out.txt"
}

print("calling R")


path_to_centerSelect = args[1]
print(path_to_centerSelect)

folder=args[2]
print(folder)

prot_name=args[3]
print(prot_name)

# defualt: 30
boxSize=as.numeric(args[4])
print(boxSize)

# default: -5
depth=as.numeric(args[5])
print(depth)

# default: 0.2
eps=as.numeric(args[6])
print(eps)


# should the pts-file be changed such that only the points in the box remain in the file?
# SELECTBOX = TRUE

# default: 0.2
# SELECTBOX=as.numeric(args[7])
SELECTBOX=as.logical(args[7])
print(SELECTBOX)


library(plot3D)
# ??border3D

#----------------------------------------------------

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180.0)}

rotate_around_x <- function(theta, points){
  # X = x;
  # 
  # Y = y*cos(theta) - z*sin(theta);
  # 
  # Z = y*sin(theta) + z*cos(theta);
  
  if(ncol(points) != 3){
    print("Wrong format!")
  }
  
  ct = cos(deg2rad(theta))
  st = sin(deg2rad(theta))
  for(i in 1:nrow(points)){
    points[i,1] = points[i,1]
    points[i,2] = points[i,2]*ct - points[i,3]*st
    points[i,3] = points[i,2]*st + points[i,3]*ct
  }
  
  return(points)
}

rotate_around_y <- function(theta, points){
  #   X = x*cos(theta) + z*sin(theta);
  # 
  # Y = y;
  # 
  # Z = z*cos(theta) - x*sin(theta);
  
  if(ncol(points) != 3){
    print("Wrong format!")
  }
  
  ct = cos(deg2rad(theta))
  st = sin(deg2rad(theta))
  for(i in 1:nrow(points)){
    points[i,2] = points[i,2]
    points[i,1] = points[i,1]*ct + points[i,3]*st
    points[i,3] = points[i,3]*ct - points[i,1]*st
  }
  
  return(points)
}

rotate_around_z <- function(theta, points){
  #   X = x*cos(theta) - y*sin(theta);
  # 
  # Y = x*sin(theta) + y*cos(theta);
  # 
  # Z = z;
  
  if(ncol(points) != 3){
    print("Wrong format!")
  }
  
  ct = cos(deg2rad(theta))
  st = sin(deg2rad(theta))
  for(i in 1:nrow(points)){
    points[i,3] = points[i,3]
    points[i,1] = points[i,1]*ct - points[i,2]*st
    points[i,2] = points[i,1]*st + points[i,2]*ct
  }
  
  return(points)
}

plot_atoms_and_iso <- function(theta, phi, prot_file, iso, iso2, rx, ry, rz){
  gridCount = 129
  origin = c(-2.079555e+01, -5.045662e+01, -4.726355e+01)
  dx = c(8.244773e-01, 0.000000e+00, 0.000000e+00)
  dy = c(0.000000e+00, 6.956269e-01, 0.000000e+00)
  dz = c(0.000000e+00, 0.000000e+00, 7.680477e-01)
  
  lx = gridCount*dx[1]
  ly = gridCount*dy[2]
  lz = gridCount*dz[3]
  
  # box3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, border = FALSE)
  # box3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, border = FALSE)
  
  points3D(prot_file$V6, prot_file$V7, prot_file$V8, add = TRUE, col = "green", theta = theta, phi = phi, cex = 0.1)
  
  # prot_trans = rotate_araound_z(rz,prot_file[,6:8])
  
  # points3D(prot_trans[1], prot_trans[2], prot_trans[3], add = TRUE, col = "green", theta = theta)
  
  # iso_translated = translate_points(iso,c(lx/2,ly/2,lz/2))
  # iso_translated = translate_points(iso,c(lx,ly,lz))
  
  # iso_translated = rotate_araound_x(rx,iso_translated)
  # iso_translated = rotate_araound_y(ry,iso_translated)
  
  iso_translated = iso
  # iso_translated = translate_points(iso_translated,c(lx/2,ly/2,lz/2))
  # iso_translated = translate_points(iso_translated,c(lx/2.0,ly/2.0,lz/2.0))
  iso_translated = translate_points(iso_translated,c(lx,ly,lz))
  
  iso_translated = rotate_around_x(rx,iso_translated)
  iso_translated = rotate_around_y(ry,iso_translated)
  iso_translated = rotate_around_z(rz,iso_translated)
  
  plot_with_less_res(theta, phi, 1,iso_translated, TRUE, "blue")
  
  # iso2_translated = iso2
  # iso2_translated = translate_points(iso2,c(lx/2,ly/2,lz/2))
  # iso2_translated = rotate_around_z(rz,iso2_translated)
  # plot_with_less_res(theta,100,iso2_translated,TRUE, "red")
}


drawBox <- function(x_pos, y_pos, size){
  lines(x = c(x_pos,x_pos+size,x_pos+size, x_pos, x_pos), y = c(y_pos,y_pos,y_pos-size,y_pos-size, y_pos),col = "red")
}


plotVol <- function(box_size,x,y,z){
  
  x_range = range(x)
  y_range = range(y)
  z_range = range(z)
  
  
  par(mfrow=c(3,1))
  plot(x,y)
  drawBox(30,-20,box_size)
  
  plot(x,z)
  drawBox(30,10,box_size)
  
  plot(y,z)
  drawBox(-20,10,box_size)
  
}

cys_cand_valid <- function(cys_candidates){
  
  chain_length = 0
  cys_start = cys_candidates[1]
  cys_curr = cys_start
  for(i in 2:length(cys_candidates)){
    if(cys_candidates[i] == cys_curr+1){
      chain_length = chain_length + 1
      cys_curr = cys_candidates[i]
    }
    else {
      return(FALSE)
    }
    
    
    if(chain_length >= 6){
      return(TRUE)
    }
    
  }
  
  return(FALSE)
}


find_motive <- function(prot_file, motive, AA_lengths){
  
  # hits = c()
  
  hits_start = c()
  hits_end = c()
  
  # col.names = c("start", "end")
  
  # find start-anchors of motive
  start_ind = which(prot_file$V4 == motive[1])
  
  # print(paste("start-indices:", start_ind, sep = " "))
  
  last_start_AA = "NONE"
  last_start_ind = -101
  
  for(curr_start in start_ind){
    cand_valid = TRUE
    
    # print(paste("curr_start: ", curr_start, sep = ""))
    
    # check for rest of motive
    # for(i in 1:length(motive)){
    
    i = 1
    curr_ind = curr_start - 1 + i
    
    # print(paste(last_start_AA ,"=?=", as.character(prot_file$V4[curr_start]) , sep = " "))
    
    if(last_start_AA == as.character(prot_file$V4[curr_start]) && curr_start - last_start_ind < 20){
      cand_valid = FALSE
      
    }
    else {
      last_start_AA = as.character(prot_file$V4[curr_start])
      last_start_ind = curr_start
    }
    
    
    
    while(i <= length(motive) && cand_valid == TRUE){
      
      
      
      
      # print(paste(motive[i], "=?=",  prot_file$V4[curr_ind], curr_ind, "<=?", length(prot_file$V4), sep = " "))
      if(!(motive[i] == prot_file$V4[curr_ind] && curr_ind <= length(prot_file$V4))){
        # print("candidate not valid!")
        cand_valid = FALSE
        # print(paste("AA-chain length: ", i))
        # break
      }
      
      # check if the consecutive atoms also belong to the current AA
      while(prot_file$V4[curr_ind] == motive[i]){
        curr_ind = curr_ind + 1
        
        # print(paste(curr_ind, prot_file$V4[curr_ind], "==", motive[i], sep = " "))
      }
      
      # print("AA ended -------------------------------------------------")
      i = i+1
    
    }
      
    # }
    
    if(cand_valid == TRUE){
      print(paste("Valid candidate found at", curr_start, sep = " "))
      # hits = c(hits,curr_start)
      hits_start = c(hits_start,curr_start)
      hits_end = c(hits_end, curr_ind)

    }
  }
  
  hits = data.frame("start" = hits_start, "end" = hits_end)

  return(hits)
}

find_motive_with_reduced_chain <- function(reduced_chain, motive){
  
  print("find_motive_with_reduced_chain() ...")
  
  hits_start = c()
  hits_end = c()
  
  # find start-anchors of motive
  start_ind = which(reduced_chain$chain == motive[1])
  
  print(start_ind)
  
  for(curr_start in start_ind){
    valid = TRUE
    
    print(curr_start)
    for(j in 1:length(motive)){
      
      if(motive[j] != reduced_chain$chain[curr_start+j-1]){
        print(paste(motive[j], "!=", reduced_chain$chain[curr_start+j-1]))
        valid = FALSE
      }
    }
    if(valid == TRUE){
      print(paste("Valid candidate found at", curr_start, sep = " "))
      # hits = c(hits,curr_start)
      hits_start = c(hits_start,reduced_chain$start_indices[curr_start])
      hits_end = c(hits_end, reduced_chain$end_indices[curr_start + (length(motive)-1)])
      
    }
  }
  
  hits = data.frame("start" = hits_start, "end" = hits_end)
  return(hits)
}

reduce_AA_chain <- function(prot_file, AA_lengths){
  
  # remove water.
  # Should actually not contain water at this point...
  water2beRemoved = which(as.character(prot_file$V4) == "HOH")
  
  if(length(water2beRemoved) != 0){
    print(paste("removing water ", length(water2beRemoved)))
    
    prot_file = prot_file[-water2beRemoved,]
  }
  print(paste("length of prot_file is ", nrow(prot_file)))
  
  # return()

  valid = TRUE
  curr_AA = "NONE"
  last_AA = "NONE"
  
  chain = c()
  start_indices = c()
  end_indices = c()
  
  i = 1
  while(i < length(prot_file$V4)){
    curr_AA = as.character(prot_file$V4[i])
    print(curr_AA)
      
    for(j in 1:((AA_lengths[curr_AA][1,] - 1))){
        print(paste(i, j, AA_lengths[curr_AA][1,], curr_AA, prot_file$V4[i+j-1]))
        if(i+j-1 < length(prot_file$V4) && prot_file$V4[i+j-1] != curr_AA){
          
          # chain was shorter than the reference
          if(j < 7){
            print("chain was shorter than the reference")
            break
          }
          
          # chain was fine, and the next AA starts
          else if(abs(j - AA_lengths[curr_AA][1,]) < 4){
            print("chain was fine, and the next AA starts")

            # chain = c(chain, curr_AA)

            # i = i + j
            break
          }
          
          # chain has errors
          else {
            
            valid = FALSE
            print("pqr-file is not conform with the specified lengths of the AAs!")
            print(paste("expected lenght of", curr_AA, AA_lengths[curr_AA][1,], "; actual length: ", j, "at position ", i))
            
            return(FALSE)
          }

        }

        
    }
    
    # if(j == AA_lengths[curr_AA][1,]){
    #   chain = c(chain, curr_AA)
    # }
    
    # chain was fine, and the next AA starts
    if(abs(j - AA_lengths[curr_AA][1,]) < 4){
      print("chain was fine, and the next AA starts")
      
      chain = c(chain, curr_AA)
      start_indices = c(start_indices, i)
      end_indices = c(end_indices, i+j-1)
    } 
    
    
    # if(i > 300) return(FALSE)
    
    i = i + j
    
    # i = i + AA_lengths[curr_AA][1,]
  }
  
  l = list("chain" = chain, "start_indices" = start_indices, "end_indices" = end_indices)
  
  return(l)
  
}

display_found_motive <- function(motive, prot_file, hit){
  print(motive)
  prot_file$V4[hit:(hit+length(motive)-1)]
}

plot_vol <- function(box_size, center){
  library(plot3D)
  
  x = prot_file$V6
  y = prot_file$V7
  z = prot_file$V8
  
  par(mfrow=c(1,1))
  points3D(x,y,z)
  
  # border3D(x0 = c(center[1]), 
  #          y0 = c(center[2]),
  #          z0 = c(center[3]),
  #          x1 = c(center[1] + box_size),
  #          y1 = c(center[2] + box_size),
  #          z1 = c(center[3] + box_size), add = TRUE, col = "black")
  
  b = box_size/2.0
  
  border3D(x0 = c(center[1] - b), 
           y0 = c(center[2] - b),
           z0 = c(center[3] - b),
           x1 = c(center[1] + b),
           y1 = c(center[2] + b),
           z1 = c(center[3] + b), add = TRUE, col = "black")
  
  
  
  # points3D(x,y,z, add = TRUE)
  
  # ?border3D
  
}



get_pos_in_space <- function(line, prot_file){
  
  return(c(prot_file$V6[line],prot_file$V7[line],prot_file$V8[line]))
}

select_points_radial <- function(radius, center, volume){
  
  
  x = c()
  y = c()
  z = c()
 
  # select all points whose distance is small enough
  for(line in 1:nrow(volume)){
    
    point = volume[line,]
    
    # matrix
    
    if(dist(rbind(center, point), method = "euclidean") <= radius){
      # active_centre = rbind(active_centre,point)
      
      x = c(x,point$x)
      y = c(y,point$y)
      z = c(z,point$z)
    }
  }
  
  active_centre = data.frame("X" = x, "Y" = y, "Z" = z)
  return(active_centre)
}

select_points_parallel_to_z <- function(radius, resolution, lz, center, volume){
  
  
  x = c()
  y = c()
  z = c()
  
  # resolution = 10.0
  step_size = lz/resolution
  
  z_projections = matrix(0, nrow = resolution, ncol = 3)
  for(i in 1:resolution){
    z_projections[i,1] = center[1]
    z_projections[i,2] = center[2]
    z_projections[i,3] = center[3] + step_size
  }
  
  
  # select all points whose distance is small enough
  for(line in 1:nrow(volume)){
    
    point = volume[line,]
    
    # matrix
    
    inside_flag = FALSE
    
    for(l in nrow(z_projections)){
      if(dist(rbind(z_projections[l,], point), method = "euclidean") <= radius){
        # active_centre = rbind(active_centre,point)
        
        inside_flag = TRUE
        # break
      }
    }
    
    if(inside_flag == TRUE){
      x = c(x,point$x)
      y = c(y,point$y)
      z = c(z,point$z)
    }
  }
  
  active_centre = data.frame("X" = x, "Y" = y, "Z" = z)
  return(active_centre)
}

select_points_cubic <- function(edgeSize, center, prot_file){
  
  
  x = c()
  y = c()
  z = c()
  
  # select all points whose distance is small enough
  for(line in 1:length(prot_file$V7)){
    
    point = get_pos_in_space(line, prot_file)
    
    # matrix
    
    ?dist
    
    if(dist(rbind(center, point), method = "euclidean") <= radius){
      # active_centre = rbind(active_centre,point)
      
      x = c(x,point[1])
      y = c(y,point[2])
      z = c(z,point[3])
    }
  }
  
  active_centre = data.frame("X" = x, "Y" = y, "Z" = z)
  return(active_centre)
}

AA_to_list <- function(in_m){
  out_m = c()
  for(i in 1:ncol(in_m)){
    out_m = c(out_m, as.character(in_m[1,i]))
  }
  
  return(out_m)
}


translate_points <- function(points, delta){
  for(i in 1:nrow(points)){
    points[i,1] = points[i,1] - delta[1]
    points[i,2] = points[i,2] - delta[2]
    points[i,3] = points[i,3] - delta[3]
  }
  
  return(points)
}

# lower resolution for plotting
less_resolution <- function(n, points){
  less_res <- points[seq(1, nrow(points), n),]
  return(less_res)
}

plot_with_less_res <- function(theta, phi, n,points, add_flag = FALSE, col = "blue"){
  less_res = less_resolution(n,points)
  points3D(less_res[,1], less_res[,2], less_res[,3], add = add_flag, col = col, theta = theta, phi = phi, cex = 0.01)
}


select_box_parallel_to_z <- function(boxSize, depth, center, volume){
  
  print(center)
  dimX = c(min(volume[1]), max(volume[1]))
  dimY = c(min(volume[2]), max(volume[2]))
  dimZ = c(min(volume[3]), max(volume[3]))
  print(dimX)
  print(dimY)
  print(dimZ)
  
  
  # centerNotNormX = center[1]*(dimX[2]-dimX[1])+dimX[1]
  # centerNotNormY = center[2]*(dimY[2]-dimY[1])+dimY[1]
  # centerNotNormZ = center[3]*(dimZ[2]-dimZ[1])+dimZ[1]
  # print(centerNotNormX)
  # 
  # center = c(centerNotNormX,centerNotNormY,centerNotNormZ)
  # 
  # print(center)
  
  
  x = c()
  y = c()
  z = c()
  
  bs = boxSize/2.0
  
  leftBorder = center[1] - bs
  rightBorder = center[1] + bs
  
  bottomBorder = center[2] - bs
  topBorder = center[2] + bs
  
  # select all points whose distance is small enough
  for(line in 1:nrow(volume)){
    
    point = volume[line,]
    
    if(point[1] < rightBorder && point[1] > leftBorder && point[2] < topBorder && point[2] > bottomBorder && point[3] > center[3] - depth){
      x = c(x,point$x)
      y = c(y,point$y)
      z = c(z,point$z)
    }
    
  }
  
  active_centre = data.frame("X" = x, "Y" = y, "Z" = z)
  return(active_centre)
}

draw_sel_box <- function(boxSize, depth, center){
  bs = boxSize/2.0
  
  # center[3] + 40 should be going till the area of the box
  border3D(center[1]-bs,center[2]-bs,center[3] - depth, center[1] +bs, center[2]+bs, center[3]+40, add = TRUE, col = "black")
  
}

my_createisosurf <- function(dx_data, level, eps){
  
  v1 = as.numeric(dx_data$V1)
  v2 = as.numeric(dx_data$V2)
  v3 = as.numeric(dx_data$V3)
  
  merged = rep(0, 3*length(v1))
  
  for(i in 1:length(v1)){
    merged[1 + ((i-1)*3)] = v1[i]
    merged[1+1+((i-1)*3)] = v2[i]
    merged[1+2+((i-1)*3)] = v3[i]
  }
  
  # gridCount = 129
  # origin = c(-2.079555e+01, -5.045662e+01, -4.726355e+01)
  # dx = c(8.244773e-01, 0.000000e+00, 0.000000e+00)
  # dy = c(0.000000e+00, 6.956269e-01, 0.000000e+00)
  # dz = c(0.000000e+00, 0.000000e+00, 7.680477e-01)
  # 
  # lx = gridCount*dx[1]
  # ly = gridCount*dy[2]
  # lz = gridCount*dz[3]
  
  xVals = c()
  yVals = c()
  zVals = c()
  
  for(i in 1:gridCount){
    for(j in 1:gridCount){
      for(k in 1:gridCount){
        
        ind = ((i-1)*gridCount+(j-1))*gridCount+k
        # print(paste("testing ind ", ind ))
        if(abs(merged[ind] - level) < eps){
          # print(paste("found voxel with level = ", level, sep = ""))
          xVals = c(xVals, i*dx[1] + origin[1])
          yVals = c(yVals, j*dy[2] + origin[2])
          zVals = c(zVals, k*dz[3] + origin[3])
        }
      }
    }
  }
  
  df = data.frame(xVals, yVals, zVals)
  names(df) = c("x","y","z")
  return(df)
  
}


center_of_AA_chain <- function(start, end, prot_file, normalized = TRUE){
  x_m = mean(prot_file$V6[start:end])
  y_m = mean(prot_file$V7[start:end])
  z_m = mean(prot_file$V8[start:end])
  
  if(normalized){
    x_m = x_m/(max(prot_file$V6)-min(prot_file$V6))
    y_m = y_m/(max(prot_file$V7)-min(prot_file$V7))
    z_m = z_m/(max(prot_file$V8)-min(prot_file$V8))
  }
  
  return(c(x_m,y_m,z_m))
}

first_cystein_of_AA_chain <- function(start, end, prot_file){
  
  
  # print(prot_file$V4[start:end])
  
  # check if the first base is a cystein
  if(prot_file$V4[start] !="CYS"){
    return(0)
  } 
  
  # find end of the first cystein
  endInd = start
  while(prot_file$V4[endInd] =="CYS"){
    endInd = endInd + 1
  }
  endInd = endInd-1
  
  # print(endInd)
  # print(prot_file$V4[start:endInd])

  return(center_of_AA_chain(start, endInd, prot_file))
}
#####################################################################





#--------------------------------------------------------------------
# manual tests
#--------------------------------------------------------------------
# # 
# # 
# # 
# # 
# # fileName = paste(prot_name, "/", prot_name, ".pqr", sep = "")
# # 
# # path_to_centerSelect = "/home/sysgen/Documents/LWB/centerSelect/"
# # # folder = "/home/sysgen/Documents/LWB/centerSelectTest/Redox/Output/"
# folder = "/home/willy/PredictingProteinInteractions/data/additionalPDBS_1/Output/"
# prot_name="hTrx1_C73H"
# path_to_centerSelect = "/home/sysgen/Documents/LWB/centerSelect/"
# # folder = "/home/sysgen/Documents/LWB/centerSelectTest/Redox/Output/"
# folder = "/home/willy/RedoxChallenges/Redox_old/Output/"
# prot_name="016"

# # prot_name="084"
# boxSize=-1
# depth=10
# eps=0.3
# # #
# path_to_centerSelect = "/home/willy/PredictingProteinInteractions/PreProcessingProteins/centerSelect/"
# #
# path_to_centerSelect = "/home/willy/RedoxChallenges/centerSelect"
#--------------------------------------------------------------------
# folder = "/home/sysgen/Documents/LWB/BoxTest/Output/"
# prot_name="1aba"
# path_to_centerSelect = "/home/sysgen/Documents/LWB/centerSelect/"


setwd(folder)

dx_file = paste(prot_name, "/", prot_name, "_pot.dx", sep = "")
outPath = paste(prot_name, "/", sep ="")

fileName = paste(prot_name, "/", prot_name, "HeadOff.pqr", sep = "")
prot_file = read.table(fileName)
# prot_file <- head(prot_file,-1)

AA_utility_path = paste(path_to_centerSelect, "/" ,"AA_utility", sep = "")
fileNameMotives = "motifs3letterCode.txt"
motives = read.table(paste(AA_utility_path, "/", fileNameMotives, sep = ""))

AA_lengths = read.csv2(paste(AA_utility_path, "/", "AA_lengths.csv", sep = ""))

# AA_lengths$"ALA"

# remove water???
# which(prot_file$V3 == "H3")

# myData <- prot_file[-c(9), ]
# which(myData$V3 == "H3")

# prot_file <- prot_file[ !(prot_file$V3 %in% c("H2")), ]
# prot_file <- prot_file[ !(prot_file$V3 %in% c("H3")), ]


reduced_chain = reduce_AA_chain(prot_file, AA_lengths)

reduced_chain[1]

out = c()
for(i in 1:nrow(motives)){
# for(i in 63){
  
  # print(i)
  m = AA_to_list(motives[i,])
  print(m[1])
  
  # print(which(reduced_chain$chain == m[1]))

  # o = find_motive(prot_file, m, AA_lengths)
  o = find_motive_with_reduced_chain(reduced_chain, m)
  # print(o)
  
  out = c(out,o)
}

out

prot_file[(out$start:out$end),]

# motives[64,]
# 
# which("CYS" == reduced_chain$chain)
# which("GLY" == reduced_chain$chain)
# which("TYR" == reduced_chain$chain)
# which("THR" == reduced_chain$chain)

# 
# m = AA_to_list(motives[49,])
# 
# find_motive_with_reduced_chain(reduced_chain, m)



print(paste("writing to file ", outPath, "/", prot_name, "_active_center.csv ...", sep = ""))
write.csv2(out, file = paste(outPath,"/",prot_name,"_active_center.csv",sep=""), row.names = F)


center = c()

if(SELECTBOX){
  center = center_of_AA_chain(out$start,out$end, prot_file, normalized = FALSE)
}else {
  center = center_of_AA_chain(out$start,out$end, prot_file, normalized = TRUE)
}

# changed this on 5.2.19
# center1 = center_of_AA_chain(out$start,out$end, prot_file)

# commented this lien 5.11.19
# center = first_cystein_of_AA_chain(out$start,out$end, prot_file)

df_center = data.frame(t(center))
names(df_center) = c("x","y","z")

print(paste("writing to file ", outPath, "/", prot_name, "_active_center.pts ...", sep = ""))
#write.csv2(df_center, file = paste(outPath,"/",prot_name,"_active_center.pts",sep=""), dec = ".", row.names = F)
write.table(df_center, file = paste(outPath,"/",prot_name,"_active_center.pts",sep=""), dec = ".", row.names = F, sep = ";")


# boxSize = 20
# depth = 10
# eps = 0.3


# SELECTBOX = TRUE
# before was quit
if(SELECTBOX){

print(paste("reading dx-file ",folder, "/", dx_file, " ...", sep = ""))
dxData <- read.csv(file=paste(folder,"/",dx_file ,sep=""), sep=' ', skip= 11, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)

dxData <- head(dxData,-5)
# sometimes 5 sometimes 10!!!

dx_meta = read.csv(file=paste(folder,"/",dx_file ,sep=""), sep=' ', skip= 0, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)

gridCount = as.numeric(dx_meta$V6[5])
origin = as.numeric(c(dx_meta$V2[6], dx_meta$V3[6], dx_meta$V4[6]))
dx = as.numeric(c(dx_meta$V2[7], dx_meta$V3[7], dx_meta$V4[7]))
dy = as.numeric(c(dx_meta$V2[8], dx_meta$V3[8], dx_meta$V4[8]))
dz = as.numeric(c(dx_meta$V2[9], dx_meta$V3[9], dx_meta$V4[9]))

lx = gridCount*dx[1]
ly = gridCount*dy[2]
lz = gridCount*dz[3]


iso = my_createisosurf(dx_data = dxData, 1.0, eps = eps)
iso2 = my_createisosurf(dx_data = dxData, -1.0, eps = eps)





pos_surf_near_act_cent = 0
neg_surf_near_act_cent = 0

pos_abrev = "_pot_positive.pts"
neg_abrev = "_pot_negative.pts"
if(boxSize > 0){
  pos_surf_near_act_cent = select_box_parallel_to_z(boxSize, depth, center, iso)
  neg_surf_near_act_cent = select_box_parallel_to_z(boxSize, depth, center, iso2)
  
  # print(paste("writing to file ", outPath, "/", prot_name, "_active_center_", pos_abrev, " ...", sep = ""))
  # write.csv2(pos_surf_near_act_cent, file = paste(outPath,"/",prot_name,"_active_center_",pos_abrev,sep=""), row.names = F)
  # 
  # print(paste("writing to file ", outPath, "/", prot_name, "_active_center_", neg_abrev, " ...", sep = ""))
  # write.csv2(neg_surf_near_act_cent, file = paste(outPath,"/",prot_name, "_active_center_",neg_abrev,sep=""), row.names = F)
} else {
  pos_surf_near_act_cent = iso
  neg_surf_near_act_cent = iso2
  
  # print(paste("writing to file ", outPath, "/", prot_name, pos_abrev, " ...", sep = ""))
  # write.csv2(pos_surf_near_act_cent, file = paste(outPath,"/",prot_name,pos_abrev,sep=""), row.names = F)
  # 
  # print(paste("writing to file ", outPath, "/", prot_name, neg_abrev, " ...", sep = ""))
  # write.csv2(neg_surf_near_act_cent, file = paste(outPath,"/",prot_name,neg_abrev,sep=""), row.names = F)
}

# overwrite existing pts-file
print(paste("writing to file ", outPath, "/", prot_name, pos_abrev, " ...", sep = ""))
write.csv2(pos_surf_near_act_cent, file = paste(outPath,"/",prot_name,pos_abrev,sep=""), row.names = F)

print(paste("writing to file ", outPath, "/", prot_name,  neg_abrev, " ...", sep = ""))
write.csv2(neg_surf_near_act_cent, file = paste(outPath,"/",prot_name,neg_abrev,sep=""), row.names = F)



make_drawings <- function(){
  # svg(filename = paste(prot_name, "/active_center.svg", sep = ""),width = 8, height = 10.64)
  png(paste(prot_name, "/active_center.png", sep = ""))
  par(mfrow = c(3,2))
  
  theta = 0
  phi = 90
  
  border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
  plot_with_less_res(theta,phi,1,iso,TRUE)
  plot_with_less_res(theta,phi,1,iso2,TRUE, "red")
  
  # points3D(active_centre$X, active_centre$Y, active_centre$Z, add = TRUE, col = "black")
  
  
  border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
  plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
  plot_with_less_res(theta,phi,1,pos_surf_near_act_cent,TRUE, "blue")
  # plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
  
  draw_sel_box(boxSize, depth, center)
  
  
  theta = 0
  phi = 135
  
  border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
  plot_with_less_res(theta,phi,1,iso,TRUE)
  plot_with_less_res(theta,phi,1,iso2,TRUE, "red")
  
  # points3D(active_centre$X, active_centre$Y, active_centre$Z, add = TRUE, col = "black")
  
  
  border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
  plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
  plot_with_less_res(theta,phi,1,pos_surf_near_act_cent,TRUE, "blue")
  # plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
  
  draw_sel_box(boxSize, depth, center)
  
  theta = 0
  phi = 180
  
  border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
  plot_with_less_res(theta,phi,1,iso,TRUE)
  plot_with_less_res(theta,phi,1,iso2,TRUE, "red")
  
  # points3D(active_centre$X, active_centre$Y, active_centre$Z, add = TRUE, col = "black")
  
  
  border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
  plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
  plot_with_less_res(theta,phi,1,pos_surf_near_act_cent,TRUE, "blue")
  # plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
  
  draw_sel_box(boxSize, depth, center)
  
  dev.off()
}


df_center = data.frame(t(center))
names(df_center) = c("x","y","z")

print(paste("writing to file ", outPath, "/", prot_name, "_active_center.pts ...", sep = ""))
#write.csv2(df_center, file = paste(outPath,"/",prot_name,"_active_center.pts",sep=""), dec = ".", row.names = F)
write.table(df_center, file = paste(outPath,"/",prot_name,"_active_center.pts",sep=""), dec = ".", row.names = F, sep = ";")
make_drawings()

}


quit()

#---------------------------------------------------------------------------------------
# pre 22.4.2019
# we have a new way of generating the iso-surfaces
#---------------------------------------------------------------------------------------
if(FALSE){
  print(paste("reading dx-file ",folder, "/", dx_file, " ...", sep = ""))
  dxData <- read.csv(file=paste(folder,"/",dx_file ,sep=""), sep=' ', skip= 11, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)
  
  dxData <- head(dxData,-5)
  # sometimes 5 sometimes 10!!!
  
  dx_meta = read.csv(file=paste(folder,"/",dx_file ,sep=""), sep=' ', skip= 0, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)
  
  gridCount = as.numeric(dx_meta$V6[5])
  origin = as.numeric(c(dx_meta$V2[6], dx_meta$V3[6], dx_meta$V4[6]))
  dx = as.numeric(c(dx_meta$V2[7], dx_meta$V3[7], dx_meta$V4[7]))
  dy = as.numeric(c(dx_meta$V2[8], dx_meta$V3[8], dx_meta$V4[8]))
  dz = as.numeric(c(dx_meta$V2[9], dx_meta$V3[9], dx_meta$V4[9]))
  
  lx = gridCount*dx[1]
  ly = gridCount*dy[2]
  lz = gridCount*dz[3]
  
  
  iso = my_createisosurf(dx_data = dxData, 1.0, eps = eps)
  iso2 = my_createisosurf(dx_data = dxData, -1.0, eps = eps)
  
  
  
  # boxSize = 30
  # depth = -5
  pos_surf_near_act_cent = 0
  neg_surf_near_act_cent = 0
  
  pos_abrev = "_pot_positive.pts"
  neg_abrev = "_pot_negative.pts"
  if(boxSize > 0){
    pos_surf_near_act_cent = select_box_parallel_to_z(boxSize, depth, center, iso)
    neg_surf_near_act_cent = select_box_parallel_to_z(boxSize, depth, center, iso2)
  
    print(paste("writing to file ", outPath, "/", prot_name, "_active_center_", pos_abrev, " ...", sep = ""))
    write.csv2(pos_surf_near_act_cent, file = paste(outPath,"/",prot_name,"_active_center_",pos_abrev,sep=""), row.names = F)
  
    print(paste("writing to file ", outPath, "/", prot_name, "_active_center_", neg_abrev, " ...", sep = ""))
    write.csv2(neg_surf_near_act_cent, file = paste(outPath,"/",prot_name, "_active_center_",neg_abrev,sep=""), row.names = F)
  } else {
    pos_surf_near_act_cent = iso
    neg_surf_near_act_cent = iso2
  
    print(paste("writing to file ", outPath, "/", prot_name, pos_abrev, " ...", sep = ""))
    write.csv2(pos_surf_near_act_cent, file = paste(outPath,"/",prot_name,pos_abrev,sep=""), row.names = F)
  
    print(paste("writing to file ", outPath, "/", prot_name, neg_abrev, " ...", sep = ""))
    write.csv2(neg_surf_near_act_cent, file = paste(outPath,"/",prot_name,neg_abrev,sep=""), row.names = F)
  }
  
  
  make_drawings <- function(){
    # svg(filename = paste(prot_name, "/active_center.svg", sep = ""),width = 8, height = 10.64)
    png(paste(prot_name, "/active_center.png", sep = ""))
    par(mfrow = c(3,2))
  
    theta = 0
    phi = 90
  
    border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
    plot_with_less_res(theta,phi,1,iso,TRUE)
    plot_with_less_res(theta,phi,1,iso2,TRUE, "red")
  
    # points3D(active_centre$X, active_centre$Y, active_centre$Z, add = TRUE, col = "black")
  
  
    border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
    plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
    plot_with_less_res(theta,phi,1,pos_surf_near_act_cent,TRUE, "blue")
    # plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
  
    draw_sel_box(boxSize, depth, center)
  
  
    theta = 0
    phi = 135
  
    border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
    plot_with_less_res(theta,phi,1,iso,TRUE)
    plot_with_less_res(theta,phi,1,iso2,TRUE, "red")
  
    # points3D(active_centre$X, active_centre$Y, active_centre$Z, add = TRUE, col = "black")
  
  
    border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
    plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
    plot_with_less_res(theta,phi,1,pos_surf_near_act_cent,TRUE, "blue")
    # plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
  
    draw_sel_box(boxSize, depth, center)
  
    theta = 0
    phi = 180
  
    border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
    plot_with_less_res(theta,phi,1,iso,TRUE)
    plot_with_less_res(theta,phi,1,iso2,TRUE, "red")
  
    # points3D(active_centre$X, active_centre$Y, active_centre$Z, add = TRUE, col = "black")
  
  
    border3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, col = "black", theta = theta, phi = phi)
  
    plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
    plot_with_less_res(theta,phi,1,pos_surf_near_act_cent,TRUE, "blue")
    # plot_with_less_res(theta,phi,1,neg_surf_near_act_cent,TRUE, "red")
  
    draw_sel_box(boxSize, depth, center)
  
    dev.off()
  }
  
  
  
  make_drawings()
}







