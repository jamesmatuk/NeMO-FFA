library(pracma)

##### Supplemental functions #####
point_in_triangle <- function(triangle,query_point){
  v_point <- query_point - triangle[1,]
  v0 <- triangle[1,] - triangle[1,]
  v1 <- triangle[2,] - triangle[1,]
  v2 <- triangle[3,] - triangle[1,]
  
  a <- (det(cbind(v_point,v2)) - det(cbind(v0,v2)))/det(cbind(v1,v2))
  b <- -(det(cbind(v_point,v1)) - det(cbind(v0,v1)))/det(cbind(v1,v2))
  in_triangle <- FALSE
  if(a>=0 & b>=0 & (a + b)<=1){
    in_triangle <- TRUE
  }
  return(in_triangle)
}


# specify time grid to evaluate kernels on
t_grid <- 1:12/6
s_grid <- 0:12/6
#t_grid <- seq(1/6,12/6,length.out = 50)
M_t <- length(t_grid)
M_s <- length(s_grid)
st_grid <- meshgrid(s_grid,t_grid)

gridpoints <- cbind(c(st_grid$X),c(st_grid$Y))

# break observation domain into a square then triangle
centers_t <- seq(t_grid[1],t_grid[M_t],length.out = num_centers)
centers_s <- seq(s_grid[1],s_grid[M_s],length.out = num_centers)
spacing_t <- diff(unique(centers_t))[1]
spacing_s <- diff(unique(centers_s))[1]
centers_grid <- meshgrid(centers_s,centers_t)

triangle_diff <- array(0,dim= c(3,2,6))
triangle_diff[,,1] <- rbind(c(0,0),
                            c(spacing_s,0),
                            c(spacing_s,spacing_t))
triangle_diff[,,2] <- rbind(c(0,0),
                            c(0,spacing_t),
                            c(spacing_s,spacing_t))
triangle_diff[,,3] <- rbind(c(0,0),
                            c(0,spacing_t),
                            c(-spacing_s,0))
triangle_diff[,,4] <- rbind(c(0,0),
                            c(-spacing_s,0),
                            c(-spacing_s,-spacing_t))
triangle_diff[,,5] <- rbind(c(0,0),
                            c(0,-spacing_t),
                            c(-spacing_s,-spacing_t))
triangle_diff[,,6] <- rbind(c(0,0),
                            c(spacing_s,0),
                            c(0,-spacing_t))


# triangular grid where basis elements are centered
centers_triangle<- cbind(centers_grid$X[(centers_grid$X - min(centers_grid$X)) <= centers_grid$Y],
                         centers_grid$Y[(centers_grid$X- min(centers_grid$X)) <= centers_grid$Y])
num_vertices <- dim(centers_triangle)[1]

# define tent basis through the triangles
basis_funs <- array(0,dim = c(M_t,M_s,num_vertices))
for(v in 1:num_vertices){
  # what triangles could potentially make up the kernels
  potential_triangles <- array(0,dim = c(3,2,6))
  for(tri_number in 1:6){
    potential_triangles[,,tri_number] <- rep(1,3)%*%centers_triangle[v,,drop = F] + triangle_diff[,,tri_number]
  }
  # decide which triangular functions to keep
  min_centers <- min(centers_triangle)
  max_centers <- max(centers_triangle)
  keep_triangle <- rep(T,6)
  for(tri_number in 1:6){
    for(v2 in 2:3){
      if(sum(apply(round(centers_triangle,4) == rep(1,dim(centers_triangle)[1])%*%t(round(potential_triangles[v2,,tri_number],4)),1,sum) == 2)==0){
        keep_triangle[tri_number] = F
      }
    }
  }
  actual_triangle <- potential_triangles[,,keep_triangle,drop = F]
  
  triangle_area <- 1/2*spacing_t*spacing_s
  basis_fun <- array(0,dim = c(M_t,M_s))
  for(m1 in 1:M_t){
    for(m2 in 1:M_s){
      grid_point <- c(st_grid$X[m1,m2],st_grid$Y[m1,m2])
      for(tri_number in 1:sum(keep_triangle)){
        temp_triangle <- actual_triangle[,,tri_number]
        if(point_in_triangle(actual_triangle[,,tri_number],grid_point)){
          basis_fun[m1,m2] <- t(c(temp_triangle[2,1]*temp_triangle[3,2] - temp_triangle[3,1]*temp_triangle[2,2],
                                  temp_triangle[2,2] - temp_triangle[3,2],temp_triangle[3,1] - temp_triangle[2,1]))%*%
            c(1,grid_point[1],grid_point[2])/(2*triangle_area)
        }
      }
    }
  }
  basis_funs[,,v] <- abs(basis_fun)
}

