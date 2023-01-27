###############################################
#R script: bathymetry interpolation with fdaPDE
###############################################

library(raster)
library(yaml)
library(hydroGOF)

#read input parameters
yaml_file <- 'simple_fdaPDE.yaml'
input_parms <- read_yaml(yaml_file)

interpolation_script_directory <- input_parms$interpolation_script_directory


# Check errors in the input yaml file
pass_input_check <- TRUE
if(is.null(input_parms$interpolation_script_directory)) {
    print ('Error: missing input paramter <interpolation_script_directory> in input yaml file') 
    pass_input_check <- FALSE } 
if(is.null(input_parms$working_directory)) {
    print ('Error: missing input paramter <working_directory> in input yaml file') 
    pass_input_check <- FALSE}
if(is.null(input_parms$input_gr3_file)) {
    print ('Error: missing input paramter <input_gr3_file> in input yaml file') 
    pass_input_check <- FALSE}
if(is.null(input_parms$mesh_prefix)) {
    print ('Error: missing input paramter <mesh_prefix> in input yaml file') 
    pass_input_check <- FALSE}
if(is.null(input_parms$target_pts_fname)) {
    print ('Error: missing input paramter <target_pts_fname> in input yaml file') 
    pass_input_check <- FALSE}
if(is.null(input_parms$target_raster_resolution)) {
    print ('Error: missing input paramter <target_raster_resolution> in input yaml file') 
    pass_input_check <- FALSE}
if(is.null(input_parms$output_file_prefix)) {
    print ('Error: missing input paramter <output_file_prefix> in input yaml file') 
    pass_input_check <- FALSE}
if(is.null(input_parms$observation_points)) {
    print ('Error: missing input paramter <observation_points> in input yaml file') 
    pass_input_check <- FALSE}
if(is.null(input_parms$smooth_type)) {
    print ('Error: missing input paramter <smooth_type> in input yaml file')  
    pass_input_check <- FALSE
} else if (input_parms$smooth_type == "isotropic_stationary") {
    if (is.null(input_parms$lambda)){
        print ('Error: missing input paramter <lambda> in input yaml file') 
        pass_input_check <- FALSE}
} else if (input_parms$smooth_type == "anisotropic_stationary") {
    if (is.null(input_parms$lambda)){
        print ('Error: missing input paramter <lambda> in input yaml file') 
        pass_input_check <- FALSE}
    if (is.null(input_parms$aniso_ratio)){
        print ('Error: missing input paramter <aniso_ratio> in input yaml file') 
        pass_input_check <- FALSE}
    if (is.null(input_parms$aniso_direction)){
        print ('Error: missing input paramter <aniso_direction> in input yaml file') 
        pass_input_check <- FALSE}
} else if (input_parms$smooth_type == "anisotropic_nonstationary") {
    if (is.null(input_parms$case_id)){
        print ('Error: missing input paramter <case_id> in input yaml file') 
        pass_input_check <- FALSE }
    if (is.null(input_parms$direction_raster)){
        print ('Error: missing input paramter <direction_raster> in input yaml file') 
        pass_input_check <- FALSE }
#to be compatible with previous version without constant_aniso_ratio and constant_diff_coeff
#set them to FALSE if they don't exist in the input files.
    if(is.null(input_parms$constant_aniso_ratio)) { 
        constant_aniso_ratio <- FALSE } 
    else { constant_aniso_ratio <- input_parms$constant_aniso_ratio }
    if (constant_aniso_ratio & is.null(input_parms$aniso_ratio) ){
        print ('Error: missing input paramter <aniso_ratio> in input yaml file') 
        pass_input_check <- FALSE}
    if (!constant_aniso_ratio & is.null(input_parms$aniso_ratio_raster)){
        print ('Error: missing input paramter <aniso_ratio_raster> in input yaml file') 
        pass_input_check <- FALSE}
    if(is.null(input_parms$constant_diff_coeff)) { 
        constant_diff_coeff <- FALSE } 
    else { constant_diff_coeff <- input_parms$constant_diff_coeff }
    if (constant_diff_coeff & is.null(input_parms$diff_coeff) ){
        print ('Error: missing input paramter <diff_coeff> in input yaml file') 
        pass_input_check <- FALSE}
    if (!constant_diff_coeff & is.null(input_parms$diff_coeff_raster) ){
        print ('Error: missing input paramter <diff_coeff_raster> in input yaml file') 
        pass_input_check <- FALSE } 
} else {
    print ('Incorrect <smooth_type> in input yaml file') 
    pass_input_check <- FALSE 
}
if (!pass_input_check){stop ('Incomplete/incorrect input yaml file', call. = FALSE)}

source(paste0(interpolation_script_directory,'\\fdaPDE_setup.R'))
setwd(input_parms$working_directory)

#set Coordinate Reference System (crs), if not specified, default to utm10
if(is.null(input_parms$coordinate_reference_system)){
    crs <- '+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
} else {
    crs_str <- input_parms$coordinate_reference_system
    if (startsWith(crs_str, 'epsg:') || startsWith(crs_str, 'EPSG:')){
        crs <- crs(paste0("+init=",crs_str)) 
    } else {
        crs <- crs_str
    }
}

## Scaling of xy data seems to have some effects on fdaPDE results.
## This script set scale to 1 and origin to (0.0,0.0)
x0 <- 0.0 
y0 <- 0.0
scale_factor <- 1. 

input_gr3_file <- input_parms$input_gr3_file
target_pts_fname <- input_parms$target_pts_fname
output_file_suffix <- '.csv'
output_file_prefix <- input_parms$output_file_prefix


#for debug only
#file_chk <- paste0(mesh_prefix,output_file_prefix,'_K_func_points_chk.csv')

#set up FEM mesh
gr3_obj1 <- read_gr3_mod(input_gr3_file, scale_factor)
print (paste('gr3_obj1: ',object.size(gr3_obj1)/1000000.0,' MB'))
#adjust origin of the coordinates to (x0, y0)
gr3_obj1@nodes[,1] <- gr3_obj1@nodes[,1]-x0
gr3_obj1@nodes[,2] <- gr3_obj1@nodes[,2]-y0
mesh <- create.mesh.2D(nodes = gr3_obj1@nodes[,1:2], segments = gr3_obj1@bnd_segments, 
                       triangles = data.frame(gr3_obj1@elms[,3:5]), order = 1)
FEMbasis <- create.FEM.basis(mesh)
print (paste('mesh: ',object.size(mesh)/1000000.0,' MB'))
print (paste('FEMbasis: ',object.size(FEMbasis)/1000000.0,' MB'))

## set up observation points as the bases for interpolation
for (i in seq(1:length(input_parms$observation_points))) {
    fname <- input_parms$observation_points[[i]][[1]]
    print(paste('reading file:',fname))
    obs_pts <- read.csv(fname, header=T)
    print(dim(obs_pts))
    obs_pts <- cbind(obs_pts$POINT_X, obs_pts$POINT_Y, obs_pts$POINT_Z)
    colnames(obs_pts)<- c('x','y','z')
     if (i == 1) {
         ref_pts <- obs_pts }
     else {
         ref_pts <- rbind(ref_pts, obs_pts) }
}
print (paste('ref_pts: ',object.size(ref_pts)/1000000.0,' MB'))
ref_locations <- matrix(c(ref_pts[,1]*scale_factor-x0,ref_pts[,2]*scale_factor-y0),ncol=2)
ref_elevations <- ref_pts[,3]

## set up boundary condition
if (length(gr3_obj1@land_BC_indices) == 1 && gr3_obj1@land_BC_indices == -1) {
    BC <- NULL
} else {
    land_BC_indices <- gr3_obj1@land_BC_indices
    land_BC_indices <- unique(land_BC_indices)
    land_BC_values <- gr3_obj1@nodes[land_BC_indices,3]
    BC <- list(BC_indices=land_BC_indices, BC_values=land_BC_values)
}

## set up target points
target_pts <- read.csv(input_parms$target_pts_fname, header=T)
print(paste('reading target points:',input_parms$target_pts_fname))
x <- target_pts$POINT_X*scale_factor - x0
y <- target_pts$POINT_Y*scale_factor - y0
target_locations <- matrix(c(x,y),ncol=2)
print (paste('target_pts: ',object.size(target_pts)/1000000.0,' MB'))
remove(target_pts)

## set up base for target raster
target_raster_res <- input_parms$target_raster_resolution
target_ras_base <- create_base_target_raster_layer(target_locations, target_raster_res, crs)


n_lambda <- switch(input_parms$smooth_type, "isotropic_stationary" = length(input_parms$lambda), 
    "anisotropic_stationary" = length(input_parms$lambda), 
    "anisotropic_nonstationary" = ifelse(constant_diff_coeff,length(input_parms$diff_coeff),1))
n_aniso_ratio <- switch(input_parms$smooth_type, "isotropic_stationary" = 1, 
    "anisotropic_stationary" = length(input_parms$aniso_ratio), 
    "anisotropic_nonstationary" = ifelse(constant_aniso_ratio,length(input_parms$aniso_ratio),1))

# initialize a matrix to store node elevations of individual runs
# the first column store values from the input
result_summary <- matrix(rep(0), nrow=gr3_obj1@nnode, ncol=n_lambda*n_aniso_ratio+1)
result_summary[,1] <- gr3_obj1@nodes$z
# initialize a matrix to store root-mean-square-error (rmse) of individual runs
rmse_summary <- matrix(rep(0),nrow=n_aniso_ratio,ncol=n_lambda)
# initialize a list to store resulting objects from fdaPDE functions
results <- list()
i <- 1
j <- 1
## 
if (input_parms$smooth_type == 'isotropic_stationary') {
    lambda_list <- input_parms$lambda
    for (lambda in lambda_list) {
        print (paste('Start processing isotropic case:', 'lambda=',lambda, date()))
        flush.console()
        output_file_header <- paste0(output_file_prefix,'_iso_',format(lambda,scientific=F))
        output_fnames <- create_output_file_names(output_file_header, output_file_suffix)
        elv_smooth3 <- smooth.FEM(locations = ref_locations, observations = ref_elevations,
                                        FEMbasis, lambda = lambda,
                                        BC = BC, GCV = FALSE)
        print ('Exit fdaPDE smoothing function')
        flush.console() 
        write_results_to_files(elv_smooth3, output_fnames)
        # calculate root-mean-square-error (rmse) for reference points
        fit_pt_values <- eval.FEM(elv_smooth3$fit.FEM,ref_locations)
        rmse_summary[j,i] <- rmse(fit_pt_values[,1],ref_elevations, na_rm = TRUE)
        result_summary[,i+1] <- elv_smooth3$fit.FEM$coeff
        print (paste('Finish processing isotropic case:', 'lambda=',lambda, date()))
        flush.console()
        i <- i+1
    }
} else if (input_parms$smooth_type == 'anisotropic_stationary') {
     lambda_list <- input_parms$lambda
     aniso_ratio_list <- input_parms$aniso_ratio
     aniso_direction <- input_parms$aniso_direction
     for (aniso_ratio in aniso_ratio_list) {
         i <- 1
         for (lambda in lambda_list) {
             print (paste('Start processing case:','anisotropic ratio=',aniso_ratio,'lambda=',lambda,date()))
             flush.console()
             output_file_header <- paste0(output_file_prefix,'_aniso_',aniso_ratio,'_',format(lambda_list[i],scientific=F))
             output_fnames <- create_output_file_names(output_file_header, output_file_suffix)
             diffusion_matrix <- find_diff_coeff_matrix3(aniso_direction/180.*pi, aniso_ratio, 1)
             print(diffusion_matrix)
             flush.console()
             PDE_parameters <- list(K=diffusion_matrix, b=c(0,0), c=0)
             # elv_smooth3 <- smooth.FEM(locations = ref_locations, observations = ref_elevations,
             #                                     FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters,
             #                                     BC = BC, GCV = FALSE)
             elv_smooth3 <- smooth.FEM(locations = ref_locations, observations = ref_elevations,
                                       FEMbasis = FEMbasis, lambda = lambda, PDE_parameters = PDE_parameters,
                                       BC = BC, GCV.inflation.factor = 1.5)
             print ('Exit fdaPDE smoothing function')
             flush.console() 
             write_results_to_files(elv_smooth3, output_fnames)
             #calculate root-mean-square-error (rmse) for reference points
             fit_pt_values <- eval.FEM(elv_smooth3$fit.FEM,ref_locations)
             rmse_summary[j,i] <- rmse(fit_pt_values[,1],ref_elevations, na_rm = TRUE)
             result_summary[,i+1+(j-1)*2] <- elv_smooth3$fit.FEM$coeff
             print (paste('Finish processing case:','anisotropic ratio=',aniso_ratio,'lambda=',lambda,date()))
             flush.console()
             i <- i+1
        }
    j <- j+1 }
} else {
    output_file_header <- paste0(output_file_prefix,'_sp_var_',input_parms$case_id)
    print (paste('memory size before =', memory.size()))
    raster1 <- raster(input_parms$direction_raster)
    print (paste('direction raster: ',object.size(raster1)/1000000.0,'MB'))
    if (!constant_aniso_ratio & !constant_diff_coeff) {
        raster2 <- raster(input_parms$aniso_ratio_raster)
        print (paste('anisotropic ratio raster: ',object.size(raster2)/1000000.0,'MB'))
        raster3 <- raster(input_parms$diff_coeff_raster)
        print (paste('diffusion coefficient raster: ',object.size(raster2)/1000000.0,'MB'))
        print (paste('memory size after loading raster =', memory.size()/1000.0,'GB'))
        gc()
        print (paste('memory size after cleanup =', memory.size()/1000.0,'GB'))
        flush.console()
        output_fnames <- create_output_file_names(output_file_header, output_file_suffix) 
        print(output_fnames)
        print(paste(i,' ',j))
        result_list <- execute_sp_var_smoothing(input_parms$case_id)
        rmse_summary[i,j] <- unlist(result_list[1])
        result_summary[,n_lambda*(i-1)+j+1] <- unlist(result_list[2])
    } else if (constant_aniso_ratio & !constant_diff_coeff) {
        raster3 <- raster(input_parms$diff_coeff_raster)
        print (paste('diffusion coefficient raster: ',object.size(raster2)/1000000.0,'MB'))
        print (paste('memory size after loading raster =', memory.size()/1000.0,'GB'))
        gc()
        print (paste('memory size after cleanup =', memory.size()/1000.0,'GB'))
        flush.console()
        for (aniso_ratio in input_parms$aniso_ratio) {
            fname_ext <- paste0('_ar',aniso_ratio)
            output_fnames <- create_output_file_names(paste0(output_file_header,fname_ext), output_file_suffix)
            print(output_fnames)
            print(paste(i,' ',j))
            result_list <- execute_sp_var_smoothing(paste0(input_parms$case_id,fname_ext))
            rmse_summary[i,j] <- unlist(result_list[1])
            result_summary[,n_lambda*(i-1)+j+1] <- unlist(result_list[2])
            i<-i+1
        }
    } else if (!constant_aniso_ratio & constant_diff_coeff) {
        raster2 <- raster(input_parms$aniso_ratio_raster)
        print (paste('anisotropic ratio raster: ',object.size(raster2)/1000000.0,'MB'))
        print (paste('memory size after loading raster =', memory.size()/1000.0,'GB'))
        gc()
        print (paste('memory size after cleanup =', memory.size()/1000.0,'GB'))
        flush.console()
        for (diff_coeff in input_parms$diff_coeff) {
            fname_ext <- paste0('_dc',diff_coeff)
            output_fnames <- create_output_file_names(paste0(output_file_header,fname_ext), output_file_suffix)
            print(output_fnames)
            print(paste(i,' ',j))
            result_list <- execute_sp_var_smoothing(paste0(input_parms$case_id,fname_ext))
            rmse_summary[i,j] <- unlist(result_list[1])
            result_summary[,n_lambda*(i-1)+j+1] <- unlist(result_list[2])
            j<-j+1
        }
    } else {
        print (paste('memory size after loading raster =', memory.size()/1000.0,'GB'))
        gc()
        print (paste('memory size after cleanup =', memory.size()/1000.0,'GB'))
        flush.console()
        for (aniso_ratio in input_parms$aniso_ratio) {
            for (diff_coeff in input_parms$diff_coeff) {
                fname_ext <- paste0('_ar',aniso_ratio,'_dc',diff_coeff)
                output_fnames <- create_output_file_names(paste0(output_file_header,fname_ext), output_file_suffix)
                print(output_fnames)
                print(paste(i,' ',j))
                result_list <- execute_sp_var_smoothing(paste0(input_parms$case_id,fname_ext))
                rmse_summary[i,j] <- unlist(result_list[1])
                result_summary[,n_lambda*(i-1)+j+1] <- unlist(result_list[2])
                j<-j+1
            }
            i<-i+1
            j<-1
        }
    }
}
print ('summary of node elevations')
print (summary(result_summary))
print ('summary of rmse for reference points')
print (rmse_summary)














