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


find_motive <- function(prot_file, motive){
  
  # hits = c()
  
  hits_start = c()
  hits_end = c()
  
  # col.names = c("start", "end")
  
  # find start-anchors of motive
  start_ind = which(prot_file$V4 == motive[1])
  
  print(paste("start-indices:", start_ind, sep = " "))
  
  last_start_AA = "NONE"
  
  for(curr_start in start_ind){
    cand_valid = TRUE
    
    # check for rest of motive
    # for(i in 1:length(motive)){
    
    i = 1
    curr_ind = curr_start - 1 + i
    
    # print(paste(last_start_AA ,"=?=", as.character(prot_file$V4[curr_start]) , sep = " "))
    
    if(last_start_AA == as.character(prot_file$V4[curr_start])){
      cand_valid = FALSE
      
    }
    else {
      last_start_AA = as.character(prot_file$V4[curr_start])
    }
    
    while(i <= length(motive) && cand_valid == TRUE){
      
      # print(paste("AA-chain length: ", i))
      
      
      # print(paste(motive[i], "=?=",  prot_file$V4[curr_ind], curr_ind, "<=?", length(prot_file$V4), sep = " "))
      if(!(motive[i] == prot_file$V4[curr_ind] && curr_ind <= length(prot_file$V4))){
        print("candidate not valid!")
        cand_valid = FALSE
        break
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

select_points_radial <- function(radius, center, prot_file){
  
  
  x = c()
  y = c()
  z = c()
 
  # select all points whose distance is small enough
  for(line in 1:length(prot_file$V7)){
    
    point = get_pos_in_space(line, prot_file)
    
    # matrix
    
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
#####################################################################
path_to_centerSelect = "/home/willy/RedoxChallenges/centerSelect/"
setwd(path_to_centerSelect)
# setwd("/home/sysgen/Documents/LWB/centerSelect/")

# fileName = "000_ArsHeadOff.pqr"

# fileName = "000_Grx1HeadOff.pqr"

prot_name = "014"
dx_file = paste(prot_name, "/", prot_name, "_pot.dx", sep = "")

fileName = paste(prot_name, "/", prot_name, "HeadOff.pqr", sep = "")

outPath = paste(prot_name, "/", sep ="")

prot_file = read.table(fileName)


AA_utility_path = "AA_utility"
fileNameMotives = "motifs3letterCode.txt"
motives = read.table(paste(AA_utility_path, "/", fileNameMotives, sep = ""))

# library(plotly)

# install.packages("plotly")

# packageVersion('plotly')


# View(prot_file)


box_size = 10


plot_vol(box_size, c(10,10,10))

#####################################################################
# find CYS-motive
#####################################################################

# cys_candidates = which(prot_file$V4 == "CYS")
# 
# cys_cand_valid(cys_candidates)
# 
# 
# motive = rep("CYS",11)
# hits = find_motive(prot_file,motive)
# 
# 
# display_found_motive(motive, prot_file, hits[1])
# 
# 
# get_pos_in_space(hits[1],prot_file)
# 
# plot_vol(box_size, get_pos_in_space(hits[1],prot_file))
# 
# 
# plot_vol(box_size, get_pos_in_space(hits[2],prot_file))
# 
# 
# active_centre = select_points(10, get_pos_in_space(hits[1],prot_file), prot_file)
#   
# nrow(prot_file)
# 
# nrow(active_centre)
# 
# plot_vol(box_size, get_pos_in_space(hits[2],active_centre))

#####################################################################

out = c()
for(i in 1:nrow(motives)){
  
  print(i)
  m = AA_to_list(motives[i,])
  print(m)

  o = find_motive(prot_file, m)
  print(o)
  
  out = c(out,o)
}

out

center_of_AA_chain <- function(start, end, prot_file){
  x_m = mean(prot_file$V6)
  y_m = mean(prot_file$V7)
  z_m = mean(prot_file$V8)
  
  return(c(x_m,y_m,z_m))
}

center = center_of_AA_chain(out$start,out$end, prot_file)

radius = 5
plot_vol(radius, center)


active_centre = select_points_radial(radius = radius, center = center,prot_file = prot_file)

active_centre

points3D(active_centre$X, active_centre$Y, active_centre$Z, add = TRUE, col = "black")

#--------------------------------------------------------------------------------------

print(paste("reading dx-file ",path_to_centerSelect, "/", dx_file, " ...", sep = ""))
dxData <- read.csv(file=paste(path_to_centerSelect,"/",dx_file ,sep=""), sep=' ', skip= 11, header=F ,stringsAsFactors=FALSE,  check.names = FALSE)
dxData <- head(dxData,-5)
# sometimes 5 sometimes 10!!!

# v1 = as.numeric(dxData$V1)
# v2 = as.numeric(dxData$V2)
# v3 = as.numeric(dxData$V3)

v1 = as.numeric(dxData$V1)
v2 = as.numeric(dxData$V2)
v3 = as.numeric(dxData$V3)


size = 129
x <- c(1:size)
y <- c(1:size)
z <- c(1:size)

v_empty = rep(0,length(v3))
v_one = rep(1,length(v3))


merged <- as.vector(rbind(v1,v2,v3))
# merged <- as.vector(rbind(v_one,v_empty,v_empty))

V <- array(merged, c(size,size,size))

print("creating isosurface ...")
iso <- createisosurf(x, y, z, V, level = 1.0)
iso2 <- createisosurf(x, y, z, V, level = -1.0)


# lower resolution for plotting
less_resolution <- function(n, points){
  less_res <- points[seq(1, nrow(points), n),]
  return(less_res)
}

plot_with_less_res <- function(theta, phi, n,points, add_flag = FALSE, col = "blue"){
  less_res = less_resolution(n,points)
  points3D(less_res[,1], less_res[,2], less_res[,3], add = add_flag, col = col, theta = theta, phi = phi, cex = 0.01)
}


iso <- createisosurf(x, z, y, V, level = 1.0)

max(iso[,1])
min(iso[,1])

my_createisosurf <- function(dx_data, level){
  
  v1 = as.numeric(dxData$V1)
  v2 = as.numeric(dxData$V2)
  v3 = as.numeric(dxData$V3)
  
  merged = rep(0, 3*length(v1))
  
  for(i in 1:length(v1)){
    merged[1 + ((i-1)*3)] = v1[i]
    merged[1+1+((i-1)*3)] = v2[i]
    merged[1+2+((i-1)*3)] = v3[i]
  }
  
  gridCount = 129
  origin = c(-2.079555e+01, -5.045662e+01, -4.726355e+01)
  dx = c(8.244773e-01, 0.000000e+00, 0.000000e+00)
  dy = c(0.000000e+00, 6.956269e-01, 0.000000e+00)
  dz = c(0.000000e+00, 0.000000e+00, 7.680477e-01)
  
  lx = gridCount*dx[1]
  ly = gridCount*dy[2]
  lz = gridCount*dz[3]
  
  xVals = c()
  yVals = c()
  zVals = c()
  
  for(i in 1:gridCount){
    for(j in 1:gridCount){
      for(k in 1:gridCount){
        
        ind = ((i-1)*gridCount+(j-1))*gridCount+k
        # print(paste("testing ind ", ind ))
        if(abs(merged[ind] - level) < 0.1){
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

iso = my_createisosurf(dx_data = dxData, 1.0)
iso2 = my_createisosurf(dx_data = dxData, -1.0)


plot_with_less_res(50,30,1,iso,FALSE)
plot_with_less_res(50,30,1,iso2,TRUE, "red")

points3D(prot_file$V6, prot_file$V7, prot_file$V8, add = TRUE, col = "green", cex = 0.1)



# origin -2.079555e+01 -5.045662e+01 -4.726355e+01
# delta 8.244773e-01 0.000000e+00 0.000000e+00
# delta 0.000000e+00 6.956269e-01 0.000000e+00
# delta 0.000000e+00 0.000000e+00 7.680477e-01


th = 80
ph = 40
box3D(origin[1],origin[2],origin[3], origin[1] +lx, origin[2]+ly, origin[3]+lz, add = FALSE, border = FALSE, theta = th, phi = ph)
plot_atoms_and_iso(th, ph, prot_file, iso, iso2, 0, 0, 0)


r = rotate_around_z(90,iso)


View(r)

prot_trans = prot_file[1:10,6:8]

prot_trans = rotate_around_z(0,prot_trans)

t = matrix(c(0,0,0, 1,1,1, 0,1,0), ncol = 3, byrow = TRUE)

points3D(t[,1], t[,2], t[,3], add = FALSE, col = "green")
t = rotate_around_z(90,t)

points3D(t[,1], t[,2], t[,3], add = FALSE, col = "green")



print(paste("writing to file ", outPath, "/", prot_name, "_pos.pts ...", sep = ""))
write.csv2(iso, file = paste(outPath,"/",prot_name,"_pos.pts",sep=""), row.names = F)

print(paste("writing to file ", outPath, "/", fileName, "_neg.pts ...", sep = ""))
write.csv2(iso2, file = paste(outPath,"/",fileName,"_neg.pts",sep=""), row.names = F)





