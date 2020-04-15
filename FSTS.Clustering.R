FSTS.clustering <- function(data, num.clust, part.coef, stop.crit, max.iter,
                            id.var = "Patient", time.vector, value.vars){
 
  time.vector <- c(-2, -1, time.vector)
  
  patient.id <- data[, id.var]
  data <- data[, value.vars]
  
  data.normalized <- normalization(data)
  
  # Matrix that indicates if the value is missing or not
  I <- is.na(data.normalized)
  
  # Initializes the missing values randomly
  data.normalized <- initialize.missing.values(data.normalized)
  
  # Add 2 extra columns corresponding to times -1 and -2 (required by FSTS)
  data.normalized <- cbind(matrix(0, nrow = length(patient.id), ncol = 2), data.normalized)
  
  U <- generate.partition.matrix(num.clust = num.clust, num.patients = length(patient.id))
  
  # Init cluster matrix
  v <- matrix(0, nrow = num.clust, ncol = length(time.vector))
  # Init distance matrix
  D <- matrix(0, nrow = num.clust, ncol = length(patient.id))
  
  obj.fcn <- rep(NA, max.iter)
  
  for(iter in 1:max.iter){
    
    # Store the last iteration clusters and partition matrix 
    v.ant <- v
    U.ant <- U
    
    # Calculate the distance matrix
    distances <- calculate.distance.and.centroids(data.normalized, I, v, U, m = part.coef, time.vector)
    
    D <- distances$D
    v <- distances$v
    
    U <- update.partition.matrix(D, part.coef)
    
    obj.fcn[iter] <- sum((D^2)*U^part.coef)
    
    data.normalized <- calculate.missing.values(data.normalized, I, U, part.coef, v)
    
    # If we have updated once the clusters
    if(iter > 1){
      
      # If the gain in the objective function is lower than the difference
      #  we end the optimization
      if(abs(obj.fcn[iter-1] - obj.fcn[iter]) < stop.crit){
        break
      }
      
    }
    
  }
  
  # Calculate the mean and standard deviation of each patient
  means <- apply(data, 1, FUN = mean, na.rm = TRUE)
  sd <- apply(data, 1, FUN = sd, na.rm = TRUE)
  
  data.not.normalized <- sweep(data.normalized, 1, sd, "*") # Multiply by standard deviation per row
  data.not.normalized <- sweep(data.not.normalized, 1, means, "+") # Add mean per patient
  
  data.not.normalized <- data.not.normalized[, -c(1:2)]
  colnames(data.not.normalized) <- value.vars
  data.not.normalized[, id.var] <- patient.id
  
  out <- list(clusters = v,
              cost = obj.fcn,
              data.imputed = data.not.normalized,
              U = U)
  
  out 
  
}

normalization <- function(data){
  
  # Calculate the mean and standard deviation of each patient
  means <- apply(data, 1, FUN = mean, na.rm = TRUE)
  sd <- apply(data, 1, FUN = sd, na.rm = TRUE)
  
  # Normalize the data
  
  data.normalized <- sweep(data, 1, means, "-") # Remove mean by row
  data.normalized <- sweep(data.normalized, 1, sd, "/") # Devide by standard deviation per row
  
  data.normalized
  
}

generate.partition.matrix <- function(num.clust, num.patients){
  
  set.seed(94)
  
  U <- matrix(runif(num.clust*num.patients, min = 0, max = 1), 
              nrow = num.clust, ncol = num.patients)
  
  # Normalization of the U matrix
  
  for(col in 1:ncol(U)){
    U[, col] <- U[, col]/sum(U[, col])
  }
  
  return(U)
}

initialize.missing.values <- function(data.normalized){
  
  set.seed(94)
  
  for(row in 1:nrow(data.normalized)){
    
    patient.max <- max(unlist(data.normalized[row, ]), na.rm = TRUE)
    patient.min <- min(unlist(data.normalized[row, ]), na.rm = TRUE)
    
    for(col in 1:ncol(data.normalized)){
      if(is.na(data.normalized[row, col])){
        data.normalized[row,col] <- runif(1) * (patient.max - patient.min) + patient.max
      }
    }
  }
  
  data.normalized
  
}

calculate.distance.and.centroids <- function(data.normalized, I, v, U, m, time.vector){
  
  data_1 = nrow(data.normalized)
  data_2 = ncol(data.normalized)

  ak = rep(1, data_2-1)
  ck = rep(1, data_2-1)
  
  for(k in 2:(data_2-1)){
    ak[k] = -(time.vector[k+1] - time.vector[k])^2
    ck[k] = -(time.vector[k] - time.vector[k-1])^2
  }
  
  bk = -(ak + ck)
  dk = -ak
  fk = -ck
  ek = -(dk + fk)
  
  mik = matrix(0, nrow = nrow(v), ncol = (data_2-1))
  aux_mik = 0
  aux_u = 0
  for( i in 1:nrow(v)) {
    for( k in 2:(data_2-1)) {
      for( j in 1:data_1){
        aux_mik = aux_mik + ( U[i,j]^m * ( ( dk[k]*data.normalized[j,k-1] ) + 
                                             ( ek[k]*data.normalized[j,k] ) +
                                             ( fk[k]*data.normalized[j,k+1] ) ) )
        aux_u = aux_u + U[i,j]^m
      }
      mik[i,k] = -( aux_mik/aux_u );
      aux_mik = 0;
      aux_u = 0;
    }
  }
  
  v = matrix(0, nrow = nrow(v), ncol = ncol(v));
  
  for (k in 2:(data_2-1)){
    for (i in 1:nrow(v)){
      v[i,k+1] = ( mik[i,k] - (ak[k]*v[i,k-1]) - (bk[k]*v[i,k]) )/( ck[k] );
    }
  }
  
  dist_STS = matrix(0, nrow = nrow(v), ncol = data_1)
  aux_dist = 0
  for(i in 1:nrow(v)){
    for(j in 1:data_1){
      for(k in 1:(data_2-1)){
        aux_dist = aux_dist + ( ( ( v[i,k+1]-v[i,k] )/( time.vector[k+1]-time.vector[k] ) ) -
                                  ( ( data.normalized[j,k+1]-data.normalized[j,k] )/
                                      ( time.vector[k+1]-time.vector[k] ) ) )^2
      }
      dist_STS[i,j] = aux_dist
      aux_dist = 0
    }
  }
  
  out = list(D = dist_STS, v = v)
  
  out
    
}

update.partition.matrix <- function(D, m){
  
  num.clusters <- nrow(D)
  
  D.m <- D^(-1/(m-1));
  U <- D.m/(matrix(1, nrow = num.clusters, ncol = 1) %*% colSums(D.m));
  
  U
}

calculate.missing.values <- function(data, I, U, m, v){
  
  n.patients <- nrow(data)
  n.time.points <- ncol(data)
  n.clusters <- nrow(v)
  
  new.data <- data
  
  for (k in 1:n.patients){
    for (j in 1:(n.time.points-2)){
      if (I[k,j] == TRUE){
        aux_OCS_n = 0
        for (i in 1:n.clusters){
          aux_OCS_n = aux_OCS_n + ( U[i,k]^m*v[i,j+2] );
        }
        new.data[k,j+2] = aux_OCS_n/(sum(U[, k]^m));
      }
    }
  }
  
  new.data
}

