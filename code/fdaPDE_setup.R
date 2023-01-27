library(fdaPDE)
library(rgl)
library(stringr)
library(sp)
library(rgdal)
library(minixml)

read_bnd_seg_nodes <- function(all_lines, idx_start){
   n_seg_nodes <- as.numeric(unlist(strsplit(all_lines[idx_start],' '))[1])
   idx_start <- idx_start + 1
   idx_end <- idx_start + n_seg_nodes -1 
   bnd_node_list <- as.numeric(all_lines[idx_start:idx_end])
   #list1 <- bnd_node_list[1:n_seg_nodes-1]
   #list2 <- bnd_node_list[2:n_seg_nodes]
   #return (data.frame(list1,list2))
  return(bnd_node_list)
}

build_bnd_nodes_pair <- function(bnd_nodes_ind) {
   n_seg_nodes <- length(bnd_nodes_ind)
   list1 <- bnd_nodes_ind[1:n_seg_nodes-1]
   list2 <- bnd_nodes_ind[2:n_seg_nodes]
   return (data.frame(list1,list2))
}

read_bnd_nodes <- function(all_lines, idx_start) {
    n_bnd_segments <- as.numeric(unlist(strsplit(all_lines[idx_start],' '))[1])
    idx_start <- idx_start + 2
    bnd_nodes_ind <- read_bnd_seg_nodes(all_lines, idx_start)
    bnd_nodes_pair <- build_bnd_nodes_pair(bnd_nodes_ind)
    n_seg_nodes <- length(bnd_nodes_ind)
    idx_start <- idx_start + n_seg_nodes + 1
    if (n_bnd_segments > 1 ) {
        for (i in 2:n_bnd_segments ) {
           temp_ind <- read_bnd_seg_nodes(all_lines, idx_start)
           temp_pair <- build_bnd_nodes_pair(temp_ind)
           n_seg_nodes <- length(temp_ind)
           idx_start <- idx_start + n_seg_nodes + 1
           bnd_nodes_pair <- rbind(bnd_nodes_pair, temp_pair)
           bnd_nodes_ind <- c(bnd_nodes_ind, temp_ind)
        }
    }
    return (list(ind=bnd_nodes_ind,pair=bnd_nodes_pair))
}

write_gr3_output_mod <- function(gr3_mesh, output_fname, desc_text, elv_smooth, i=1, target=1, scale=1){
  #target=1: smooth z, target=2: smooth dz
    nnode <- gr3_mesh@nnode
    idx_elm_start <- nnode + 2 + 1
    all_lines <- gr3_mesh@all_lines
    nodes <- gr3_mesh@nodes
    nodes[,1] <- nodes[,1]/scale
    nodes[,2] <- nodes[,2]/scale
    #fmt_nodes <- gr3_mesh@fmt_nodes
    nnode_last_line_idx <- idx_elm_start -1
    # node_lines <- matrix(all_lines[3:nnode_last_line_idx],nrow=nnode)
    # node_lines_sel <- substr(node_lines,1,fmt_nodes[4])
    zz <- file(output_fname,"w")
    title1 <- paste(all_lines[1],desc_text) 
    writeLines(title1, con = zz, sep = "\n")
    writeLines(all_lines[2], con = zz, sep = "\n")
    if (target == 2) {
        nodes[,3] <- (elv_smooth$fit.FEM$coeff[,i]+nodes$z)*-1.0
    } else {
        nodes[,3] <- (elv_smooth$fit.FEM$coeff[,i])*-1.0
    }
    write.table(format(nodes, width=18, nsmall=8), file=zz, sep="\t", eol="\n", row.names = TRUE, col.names = FALSE, quote = FALSE)
    write.table(all_lines[idx_elm_start:length(all_lines)], file=zz, sep="\t", eol="\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
    close(zz)
}

write_lxml_output <- function(elv_smooth3, output_filename){
   doc <- XMLDocument$new('LandXML'); 
   doc$update(xmlns="http://www.landxml.org/schema/LandXML-1.2",'xmlns:xsi' = "http://www.w3.org/2001/XMLSchema-instance", 
              'xsi:schemaLocation' ="http://www.landxml.org/schema/LandXML-1.2 http://www.landxml.org/schema/LandXML-1.2/LandXML-1.2.xsd", 
               date="2020-02-13", time="19:33:54", version="1.2", language="English", readOnly="false")
    # don't include Coordinate Reference System (CRS) information for now
    # need to modify code to support different CRS in addition to UTM10 in the future
    # doc$add('Units')$ add('Metric')$ update(areaUnit="squareMeter", linearUnit="meter", 
    #                 volumeUnit="cubicMeter", temperatureUnit="celsius", pressureUnit="mmHG", 
    #                 diameterUnit="meter" ,angularUnit="decimal degrees", directionUnit="decimal degrees")
    # doc$add('CoordinateSystem')$ update(desc="NAD83 NAD_1983_UTM_Zone_10N, Meter", epsgCode="26910",
    #                 ogcWktCode="PROJCS[&quot;NAD_1983_UTM_Zone_10N&quot;,GEOGCS[&quot;GCS_North_American_1983&quot;,DATUM[&quot;D_North_American_1983&quot;,SPHEROID[&quot;GRS1980&quot;,6378137.000,298.25722210]],PRIMEM[&quot;Greenwich&quot;,0.0],UNIT[&quot;Degree&quot;,0.017453292519943295]],PROJECTION[&quot;Transverse_Mercator&quot;],PARAMETER[&quot;false_easting&quot;,500000.0],PARAMETER[&quot;false_northing&quot;,0.0],PARAMETER[&quot;central_meridian&quot;,-123.0],PARAMETER[&quot;latitude_of_origin&quot;,0.0],PARAMETER[&quot;scale_factor&quot;,0.9996],UNIT[&quot;Meter&quot;,1.0]]", 
    #                 horizontalDatum="D_North_American_1983", horizontalCoordinateSystemName="NAD_1983_UTM_Zone_10N")
    elm <- doc$add('Surfaces') $add('Surface') $update(name=str_remove(output_filename,'.xml'), desc ='interpolated surface by fdaPDE') $add('Definition') $update(surfType="TIN")
    elm_pnts <- elm$add('Pnts')
    nnode <- dim(elv_smooth3[[1]][[2]]$mesh[[1]])[1]
    for (i in 1:nnode) {
      elm_pnts$add('P') $update(id=i) $update(paste(elv_smooth3[[1]][[2]]$mesh[[1]][i,2], 
                                      elv_smooth3[[1]][[2]]$mesh[[1]][i,1], elv_smooth3[[1]][[1]][i]))
    }
    elm_faces <- elm$add('Faces')
    nelm <- dim(elv_smooth3[[1]][[2]]$mesh[[4]])[1]
    for (i in 1:nelm) {
      elm_faces$add('F') $update(paste(elv_smooth3[[1]][[2]]$mesh[[4]][i,]))
    }
  doc$save(output_filename)
}

write_lxml_output_from_gr3 <- function(gr3_obj, output_filename){
  doc <- XMLDocument$new('LandXML'); 
  doc$update(xmlns="http://www.landxml.org/schema/LandXML-1.2",'xmlns:xsi' = "http://www.w3.org/2001/XMLSchema-instance", 
             'xsi:schemaLocation' ="http://www.landxml.org/schema/LandXML-1.2 http://www.landxml.org/schema/LandXML-1.2/LandXML-1.2.xsd", 
             date="2020-02-13", time="19:33:54", version="1.2", language="English", readOnly="false")
  doc$add('Units')$ add('Metric')$ update(areaUnit="squareMeter", linearUnit="meter", 
                                          volumeUnit="cubicMeter", temperatureUnit="celsius", pressureUnit="mmHG", 
                                          diameterUnit="meter" ,angularUnit="decimal degrees", directionUnit="decimal degrees")
  doc$add('CoordinateSystem')$ update(desc="NAD83 NAD_1983_UTM_Zone_10N, Meter", epsgCode="26910",
                                      ogcWktCode="PROJCS[&quot;NAD_1983_UTM_Zone_10N&quot;,GEOGCS[&quot;GCS_North_American_1983&quot;,DATUM[&quot;D_North_American_1983&quot;,SPHEROID[&quot;GRS1980&quot;,6378137.000,298.25722210]],PRIMEM[&quot;Greenwich&quot;,0.0],UNIT[&quot;Degree&quot;,0.017453292519943295]],PROJECTION[&quot;Transverse_Mercator&quot;],PARAMETER[&quot;false_easting&quot;,500000.0],PARAMETER[&quot;false_northing&quot;,0.0],PARAMETER[&quot;central_meridian&quot;,-123.0],PARAMETER[&quot;latitude_of_origin&quot;,0.0],PARAMETER[&quot;scale_factor&quot;,0.9996],UNIT[&quot;Meter&quot;,1.0]]", 
                                      horizontalDatum="D_North_American_1983", horizontalCoordinateSystemName="NAD_1983_UTM_Zone_10N")
  elm <- doc$add('Surfaces') $add('Surface') $update(name=gr3_obj@all_lines[1], desc ='converted from gr3 file') $add('Definition') $update(surfType="TIN")
  elm_pnts <- elm$add('Pnts')
  nnode <- gr3_obj@nnode
  for (i in 1:nnode) {
    elm_pnts$add('P') $update(id=i) $update(paste(gr3_obj@nodes[i,2], 
                                                  gr3_obj@nodes[i,1], gr3_obj@nodes[i,3]))
  }
  elm_faces <- elm$add('Faces')
  nelm <- gr3_obj@nelm
  for (i in 1:nelm) {
    elm_faces$add('F') $update(paste(gr3_obj@elms[i,3],gr3_obj@elms[i,4],gr3_obj@elms[i,5]))
  }
  doc$save(output_filename)
}


find_diff_coeff_matrix3 <- function (theta,Kx,Ky){
   mat<-matrix(c(1.0,0,0,1.0),nrow=2)
   if (is.na(theta) || Kx == Ky) {A <- mat*Kx}
   else {
      scale_factor <- 1.0 #1.0/sqrt(0.5*(Kx^2+Ky^2))
      diff_matrix_original <- matrix(c(Kx*scale_factor,0,0,Ky*scale_factor), nrow=2)
      R_matrix <- matrix(c(cos(theta),sin(theta),-1.0*sin(theta),cos(theta)),nrow=2)
      R_matrix_i <- matrix(c(cos(theta),sin(theta),-1.0*sin(theta),cos(theta)),nrow=2, byrow=T)
      A <- R_matrix %*% diff_matrix_original %*% R_matrix_i
    }
   return (A)
} 

find_diff_coeff_matrix4 <- function (theta,ratio,strength){
  mat<-matrix(c(1.0,0,0,1.0),nrow=2)
  if (is.na(theta) || ratio == 1) { A <- sqrt(0.5*strength^2)* mat } 
  else {
    scale_factor <- 1.0/sqrt(ratio^2+1^2)
    diff_matrix_original <- matrix(c(ratio*scale_factor,0,0,1.0*scale_factor), nrow=2)
    R_matrix <- matrix(c(cos(theta),sin(theta),-1.0*sin(theta),cos(theta)),nrow=2)
    R_matrix_i <- matrix(c(cos(theta),sin(theta),-1.0*sin(theta),cos(theta)),nrow=2, byrow=T)
    A <- strength*R_matrix %*% diff_matrix_original %*% R_matrix_i
  }
  return (A)
}

gr3_mesh <- setClass("gr3_mesh", 
    slots = c(nelm = "numeric", nnode = "numeric", nodes = "data.frame", elms = "matrix", 
              open_BC_indices = "numeric", land_BC_indices = "numeric", bnd_segments = "data.frame", 
              fmt_nodes = "numeric", all_lines = "character"))  
             

read_gr3_mod <- function (fname, scale=1) {
  # cannot hendle the case without land boundaries
  con <- file(fname,'r')
  all_lines <- readLines(con,n = -1)
  close(con)
  nelm <- as.numeric(unlist(strsplit(all_lines[2]," "))[1])
  nnode <- as.numeric(unlist(strsplit(all_lines[2]," "))[2])
  nnode_last_line_idx <- 3 + nnode -1
  node_lines <- matrix(all_lines[3:nnode_last_line_idx],nrow=nnode)
  #node_lines_sel <- substr(node_lines,1,fmt_nodes[4])
  x <- apply(node_lines,1,function(x) as.numeric(unlist(strsplit(x," +"))[2]))
  x <- x*scale
  y <- apply(node_lines,1,function(x) as.numeric(unlist(strsplit(x," +"))[3]))
  y <- y*scale
  z <- apply(node_lines,1,function(x) as.numeric(unlist(strsplit(x," +"))[4]))
  idx_elm_start <- nnode+2+1
  idx_elm_end <- idx_elm_start +nelm - 1
  elm_lines <- matrix(all_lines[idx_elm_start:idx_elm_end])
  elms <- t(matrix(as.numeric(apply(elm_lines,1,function(x) unlist(strsplit(x,'\\s+')))),nrow=5))
  idx_start_open_bnd_nodes <- nnode + nelm + 3
  n_open_bnd_segments <- as.numeric(unlist(strsplit(all_lines[idx_start_open_bnd_nodes],' '))[1])
  n_open_bnd_nodes_total <- as.numeric(unlist(strsplit(all_lines[idx_start_open_bnd_nodes+1],' '))[1])
  if (n_open_bnd_nodes_total> 0) {
    temp <- read_bnd_nodes(all_lines, idx_start_open_bnd_nodes)
    open_bnd_nodes <- temp[2]
    open_BC_indices <- unname(unlist(temp[1]))
  } else {
    open_BC_indices <- c(-1)
  }
  idx_start_land_bnd_nodes <- idx_start_open_bnd_nodes + n_open_bnd_segments + n_open_bnd_nodes_total +2
  n_land_bnd_nodes_total <- as.numeric(unlist(strsplit(all_lines[idx_start_land_bnd_nodes+1],' '))[1])
  if (n_land_bnd_nodes_total> 0) {
    temp <- read_bnd_nodes(all_lines, idx_start_land_bnd_nodes)
    land_bnd_nodes <- temp[2]
    land_BC_indices <- unname(unlist(temp[1]))
  } else {
    land_BC_indices <- c(-1)
  }
  nodes <- data.frame(x,y,z)
  nodes$z <- nodes$z*-1.0
  if ((open_BC_indices[1] > 0) & (land_BC_indices[1] > 0)) {
    bnd_segments <- rbind(data.frame(open_bnd_nodes),data.frame(land_bnd_nodes))
  } else if ((open_BC_indices[1] <= 0) & (land_BC_indices[1] > 0)) {
    bnd_segments <- data.frame(land_bnd_nodes)
  } else if ((open_BC_indices[1] > 0) & (land_BC_indices[1] <= 0)) {
    bnd_segments <- data.frame(open_bnd_nodes)
  } else {
    stop ("gr3 file doesn't have boundary information")
  }
    
  return (new("gr3_mesh", nelm = nelm, nnode = nnode, nodes = nodes,
              elms = elms, open_BC_indices = open_BC_indices, land_BC_indices = land_BC_indices, 
              bnd_segments = bnd_segments, all_lines = all_lines))
}

K_func5 <- function(points)
  #calculate K using three spatial arrays: angles, ratio, and magnitudes
{
  ndata <- nrow(points)
  points <- as.data.frame(points)
  colnames(points) <- c('x','y')
  points_org <- data.frame(points[,1]/scale_factor+x0, points[,2]/scale_factor+y0)
  #write.table(points_org, file=file_chk, append=TRUE, sep=',')
  directions <- extract(raster1,points)
  aniso_ratios <- extract(raster2,points)
  diff_coeffs <- extract(raster3,points)
  output <- array(0., c(2, 2, ndata))
  mat<-matrix(c(1.0,0,0,1.0),nrow=2)
  for (i in 1:ndata) {
    theta <- directions[i]
    ratio <- aniso_ratios[i]
    strength <- diff_coeffs[i]
    #print (paste(i, points$x[i], points$y[i], points$i_row[i],points$i_col[i], theta, ratio, strength))
    output[,,i] <- find_diff_coeff_matrix4(theta/180.0*pi, ratio, strength) 
    #print (as.numeric(output[,,i]))
    flush.console()
  }
  output
}

K_func6 <- function(points)
  #calculate K using three spatial arrays: angles, ratio, and magnitudes
  #allow aniso ratio and diff coeff to be constant 5/8/2020
{
  ndata <- nrow(points)
  points <- as.data.frame(points)
  colnames(points) <- c('x','y')
  points_org <- data.frame(points[,1]/scale_factor+x0, points[,2]/scale_factor+y0)
  #write.table(points_org, file=file_chk, append=TRUE, sep=',')
  directions <- extract(raster1,points)
  if (constant_aniso_ratio) {aniso_ratios<-rep(aniso_ratio,ndata)} 
      else {aniso_ratios <- extract(raster2,points)}
  if (constant_diff_coeff) {diff_coeffs<-rep(diff_coeff,ndata)} 
      else {diff_coeffs <- extract(raster3,points)}
  output <- array(0., c(2, 2, ndata))
  mat<-matrix(c(1.0,0,0,1.0),nrow=2)
  for (i in 1:ndata) {
    theta <- directions[i]
    ratio <- aniso_ratios[i]
    strength <- diff_coeffs[i]
    #print (paste(i, points$x[i], points$y[i], points$i_row[i],points$i_col[i], theta, ratio, strength))
    output[,,i] <- find_diff_coeff_matrix4(theta/180.0*pi, ratio, strength)
    #print (as.numeric(output[,,i]))
    flush.console()
  }
  output
}

b_func<-function(points)
{
  output <- array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 0
  output
}
c_func<-function(points)
{
  rep(c(0), nrow(points))
}
u_func<-function(points)
{
  rep(c(0), nrow(points))
}

create_output_file_names<-function(header, suffix)
{
  output_rds_file <- paste0(header,'_elv_smooth3.rds') #filename to save rds objct
  output_file_nodes <- paste0(header,'_node_output', suffix) #filename to export node elevations as cvs file
  output_file_target_points <- paste0(header, suffix) #filename to export elevations at target points as cvs file
  output_file_target_raster <- paste0(header, '.tif') #filename to export interpolated raster as tif file
  output_file_gr3 <- paste0(header, '.gr3') #filename to export fdaPDE results for FEM in gr3 format
  output_file_lxml <- paste0(header, '.xml') #filename to export fdaPDE results for FEM in LandXML format
  return (list(output_rds_file, output_file_nodes, output_file_target_points, 
               output_file_target_raster, output_file_gr3,output_file_lxml))
}

create_base_target_raster_layer <- function(target_locations, target_res, crs){
  target_locations_original <- matrix(c(target_locations[,1]/scale_factor+x0, 
                                        target_locations[,2]/scale_factor+y0),ncol=2)
  x_min <- min(target_locations_original[,1]) - 0.5*target_res
  y_min <- min(target_locations_original[,2]) - 0.5*target_res
  y_max <- max(target_locations_original[,2]) + 0.5*target_res
  x_max <- max(target_locations_original[,1]) + 0.5*target_res
  xsize <- abs((x_max-x_min)/target_res)
  ysize <- abs((y_max-y_min)/target_res)
  print(paste(x_min,y_min,x_max,y_max,xsize,ysize))
  flush.console()
  target_base_ras <- raster(ncol=xsize, nrow=ysize, xmn=x_min, xmx=x_max, ymn=y_min, ymx=y_max,crs=crs)
  return(target_base_ras)
}

write_results_to_files <- function(elv_smooth3, output_fnames, desc_text=NULL) {
  # save resulting R object on the disk for future use in case it is needed
  saveRDS(elv_smooth3, file=output_fnames[[1]])
  # write node elevations for checking purpose 
  xy <- elv_smooth3[[1]][[2]][[1]][[1]]
  z <- elv_smooth3[[1]][[1]]
  write.csv(cbind(xy,z), file=output_fnames[[2]], row.names=FALSE) 
  # derive elevations at target points and write (x,y,z) to a .csv file
  target_values <- eval.FEM(elv_smooth3$fit.FEM,target_locations)
  target_locations_original <- matrix(c(target_locations[,1]/scale_factor+x0, 
                                        target_locations[,2]/scale_factor+y0),ncol=2)
  fit_target <- cbind(target_locations_original,target_values)
  fit_target_mod <- fit_target[(!is.na(fit_target[,3])),]
  write.csv(fit_target_mod, file=output_fnames[[3]], row.names=FALSE)
  rasterize(target_locations_original,target_ras_base,target_values,fun=mean,update=TRUE,
            updateValue='all',filename=output_fnames[[4]], format="GTiff", overwrite=TRUE, crs=crs)
  if (is.null(desc_text)) { desc_text <- str_remove(output_fnames[[5]],'.gr3') }
  write_gr3_output_mod(gr3_obj1, output_fnames[[5]], desc_text, elv_smooth3)
  write_lxml_output(elv_smooth3, output_fnames[[6]])
}

execute_sp_var_smoothing<-function(case_identifier){
  print (paste0('Start processing spatial variable case: ',case_identifier,", ",date()))
  flush.console()
  PDE_parameters <- list(K = K_func6, b = b_func, c = c_func,u = u_func)
  # elv_smooth3 <- smooth.FEM(locations = ref_locations, observations = ref_elevations,
  #                                        FEMbasis = FEMbasis, lambda = 1, PDE_parameters = PDE_parameters,
  #                                        BC = BC,GCV = FALSE)
  
  elv_smooth3 <- smooth.FEM(locations = ref_locations, observations = ref_elevations,
                            FEMbasis = FEMbasis, lambda = 1, PDE_parameters = PDE_parameters,
                            BC = BC,GCV.inflation.factor = 1)
  print ('Exit fdaPDE smoothing function')
  flush.console()
  write_results_to_files(elv_smooth3, output_fnames, case_identifier)
  # calculate root-mean-square-error (rmse) for reference points
  fit_pt_values <- eval.FEM(elv_smooth3$fit.FEM,ref_locations)
  print(paste(i,j,n_aniso_ratio*(i-1)+j+i))
  case_rmse <- rmse(fit_pt_values[,1],ref_elevations, na_rm = TRUE)
  case_result <- elv_smooth3$fit.FEM$coeff
  print (paste0('Finish processing spatial variable case: ',case_identifier,", ",date()))
  return(list(case_rmse, case_result))
  flush.console()
}

