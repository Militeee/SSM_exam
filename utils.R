library(markovchain)


fixation_times_CTX_238GS <- 48

from_table_to_Q <-  function(mt, model = c("SSWM", "PCTMC")){

  if(model == "SSWM"){

    Q_list <-  Q_SSWM(mt)

  } else if (model == "PCTMC"){

    Q_list <-  Q_PCTMC(mt)

  } else {

    print(paste0("Evolutionary model", model, "not implemented"))
    return(NULL)
  }




  return(Q_list)
}


calculate_fitness_from_growth_SSWM <-  function(g1,g2,d1,d2, fun = identity){

    if(g1 >= g2 | sum(d1 != d2) > 1 ){

      return(0)

    } else {

      rate <- (fixation_times_CTX_238GS * fun(g2/g1))
      rate <-  1/rate
    }

    return(rate)
}


Q_SSWM <- function(mt){

  drugs <- rownames(mt)

  states <-  colnames(mt)

  states_mat <- sapply(strsplit(states, split = ""), function(x) as.logical(as.numeric(x)))

  MC_list <-  apply(mt, 1, function(x)  generate_MC_SSWM(x,states, states_mat))

  names(MC_list) <-  drugs

  return(MC_list)

}

generate_q <- function(vec,states_mat){

  states <-  length(vec)
  Q <- matrix(nrow = ncol(states_mat), ncol = ncol(states_mat))
  for(i in seq_len(states)){
    for(j in seq_len(states)){
          Q[i,j] <- calculate_fitness_from_growth_SSWM(vec[i], vec[j], states_mat[,i], states_mat[,j])

    }
  }

  for(i in seq_len(states)){

    Q[i,i] <- sum(Q[i,]) * -1
  }

  return(Q)
}

generate_MC_SSWM <-  function(probs, states, states_mat){

  Q <- generate_q(probs,states_mat)
  CTMC <- new("ctmc", states = states, generator = Q)
  return(CTMC)

}


plot_ctmc <-  function(ctmc){

  igraph_obj <- graph_from_adjacency_matrix(ctmc@generator, weighted = T)
  igraph_obj <- simplify(igraph_obj)
  yrPal <- colorRampPalette(c('yellow','orange', 'red'))
  grad <- yrPal(150)
  E(igraph_obj)$color <- grad[round((E(igraph_obj)$weight - min(E(igraph_obj)$weight)) / max(E(igraph_obj)$weight / 100),3) + 10]
  V(igraph_obj)$color <- rep(x = "white", length(V(igraph_obj)))
  l <- cbind(c(7,4,6,8,10,2,4,6,8,10,12,4,6,8,10,7), c(0,2,2,2,2,4,4,4,4,4,4,6,6,6,6,8))
  plot(igraph_obj, layout = l, size = 20, vertex.label.color = "black")

  return(igraph_obj)

}


multi_sim_CTC <- function(ctmc, init ,N = 10^4, NC = 1000){

  res <-  lapply(1:N, function(x) simulation_CTMC(ctmc@generator,init,  NC) )
  names(res) <- paste0("run", 1:N)
  return(res)


}


simulation_CTMC <-  function(Q,p0,n_iter){
  ind = sample.int(length(p0), prob = p0,replace = T, size = 1)
  P = Q
  traj = vector(length = n_iter)
  traj[1] = ind
  time = vector(length = n_iter)
  time[1] = 0
  i = 1
  while(i < n_iter){

    exit_rate = - Q[ind,ind]
    if(exit_rate == 0) {
      break
    }
    pp <- P[ind,]
    pp[ind] <- 0
    ind = sample.int(length(p0), prob = pp, replace = T, size = 1)
    tau = -log(runif(1))/exit_rate
    i = i+ 1
    traj[i] = ind
    time[i] = time[i-1] + tau
  }

  return(cbind(traj[1:i], time[1:i]))

}

mean_abs_time_fromSM <- function(x,y, transl_states) {

 x2 <-  x %>% filter(state %in% y)
 tot <-  nrow(x2)
 x2 %>% group_by(state) %>%
    summarise(mean = mean(time), variance = var(time), prob = dplyr::n() / tot) %>% inner_join(.,transl_states, by = 'state') %>%
      select(mutations,mean,variance, prob)

}

calculate_commutativity <- function(final_probs, init_state){

  dims <- length(final_probs)
  res <- matrix(ncol = dims, nrow = dims)
  for(i in seq_len(dims)){
    for(j in seq_len(dims)){
      res[i,j] <- all(round(init_state %*% final_probs[[i]] %*% final_probs[[j]], 1) == round(init_state %*% final_probs[[j]] %*% final_probs[[i]], 1))
    }
  }
  colnames(res) <- names(final_probs)
  rownames(res) <-  colnames(res)
  return(res)
}


calculate_coll_sens_res <- function(final_probs, init_state, fitness_mat){

  dims <- length(final_probs)
  res <- matrix(ncol = dims, nrow = dims)
  for(i in seq_len(dims)){
    for(j in seq_len(dims)){
      res[i,j] <- (round(init_state %*% final_probs[[i]] %*% final_probs[[j]], 2) %*% fitness_mat[j,]) / (round(init_state %*% final_probs[[j]],2) %*% fitness_mat[j,])
    }
  }
  colnames(res) <- names(final_probs)
  rownames(res) <-  colnames(res)
  return(res)
}


calculate_best_drugs_1 <- function(final_probs, init_state, fitness_mat){

  dims <- length(final_probs)
  res <- as.data.frame(matrix(ncol = 2, nrow = dims))
  for(i in seq_len(dims)){
    res[i, 2] <- round(init_state %*% final_probs[[i]], 2) %*% fitness_mat[i,]
    res[i, 1] <-  rownames(data)[i]
  }
  colnames(res) <- c("drugs", "GR")
  res <- res[order(res[,2]),]
  return(res)
}


calculate_best_drugs_2 <- function(final_probs, init_state, fitness_mat){

  dims <- length(final_probs)
  res <- as.data.frame(matrix(ncol = 2, nrow = dims * dims - dims))
  for(i in seq_len(dims)){
    for(j in seq_len(dims)){
      if(i == j) next
      res[(i - 1) * dims + j, 2] <- round(init_state %*% final_probs[[i]] %*% final_probs[[j]], 2) %*% fitness_mat[j,]
      res[(i -1) * dims  + j, 1] <-  paste0(rownames(data)[i] , "->" , rownames(data)[j])
    }
  }
  colnames(res) <- c("drugs", "GR")
  res <- res[order(res[,2]),]
  return(res)
}

calculate_best_drugs_3 <-  function(final_probs, init_state, fitness_mat) {

  dims <- length(final_probs)
  res <- as.data.frame(matrix(ncol = 2, nrow = dims * dims* dims - dims * dims))
  for(i in seq_len(dims)){
    for(j in seq_len(dims)){
      if(i == j) next
      for(k in seq_len(dims)){
        if(j == k) next

        res[(i - 1) * dims * dims + (j - 1) * dims + k, 2] <- round(init_state %*% final_probs[[i]] %*% final_probs[[j]] %*% final_probs[[k]], 2) %*% fitness_mat[k,]
        res[(i -1) * dims * dims + (j - 1) * dims + k, 1] <-  paste0(rownames(data)[i] , "->" , rownames(data)[j], "->", rownames(data)[k])
      }
    }
  }
  colnames(res) <- c("drugs", "GR")
  res <- res[order(res[,2]),]
  return(res)
}


generate_mut_matrix <- function(genos, mut_rate = 10^-8){

  dims <-  length(genos)
  genos <- sapply(strsplit(genos, split = ""), function(x) as.logical(as.numeric(x)))

  res <- matrix(ncol = dims, nrow = dims )

  for(i in seq_len(dims)){
    for(j in seq_len(dims)){
      if(sum(genos[,i] != genos[,j]) == 1  ){
        res[i,j] <- mut_rate
      } else {
        res[i,j] <- 0
      }

    }
  }
  for(i in seq_len(dims)){
    res[i,i] <- sum(res[i,])
  }

  return(res)

}

generate_WF_probs <-  function(N_s,genotypes,mutation_matrix, fitness) {
  probs <-  vector(length = length(genotypes))
  for(i in seq_along(genotypes)){
    tmp1 <- ( N_s[i] * (1 + fitness[i]) * (1 - mutation_matrix[i,i]) ) + sum( N_s[-i]  * mutation_matrix[i,-i])
    probs[i] <- tmp1 / (N_s[i] * (1 + fitness[i]) + sum(N_s[-i]))

  }
  return(probs)
}

selection_coeff <- function(fitt) {

  yy <-  max(fitt)

  return((fitt / yy) - 0.6)
}

WF <- function(genotypes, mutation_matrix, fitness, N = 1000, P = 10^8, init_state = 1) {

  N_s <- matrix(0L, nrow = length(genotypes), ncol = N+1)
  P_s <- matrix(0L, nrow = length(genotypes), ncol = N+1)
  P_s[1,1] <- P
  N_s[1,1] <- 1
  for(i in seq_len(N) + 1){
   probs <-  generate_WF_probs(P_s[,i-1], genotypes, mutation_matrix, fitness)
   samples <-  rbinom(length(probs), size = P, prob = probs / sum(probs))
   #P <-  sum(samples)
   N_s[,i] <-  samples / P

   P_s[,i] <-  N_s[,i] * P
  }
  colnames(N_s) <- paste0(0:N)
  rownames(N_s) <- genotypes
  return(N_s)
}

