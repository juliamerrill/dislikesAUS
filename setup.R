library(tidyverse)
messagef <- function(...) message(sprintf(...))

setup_workspace2 <- function(fname = "data/dislikes_data.csv"){
  master <- readr::read_csv2(fname)
  bad_ids <- master %>% count(p_id) %>% filter(n != 49) %>% pull(p_id)
  bad_ids2 <- master %>% group_by(p_id) %>% summarise(m = sum(is.na(SMP.liking))) %>% filter(m == 49) %>% pull(p_id)
  master <- master %>% filter(!(p_id %in% union(bad_ids, bad_ids2)))
  master <- master %>%
    mutate(p_id = sprintf("%04d", as.integer(factor(p_id))),
           is_fan = case_when(is.na(SMP.liking) ~ F, SMP.liking > 4 ~ T, T ~ F))
  assign("master2", master, globalenv())

}

comp_style_pref <- function(fam1, lik1, fam2, lik2){
  w_lik1 <- fam1 * lik1
  w_lik2 <- fam2 * lik2
  cos_lik <- sum(w_lik1 * w_lik2, na.rm = T)
  w_lik1_norm <- sqrt(sum(w_lik1*w_lik1, na.rm = T))
  w_lik2_norm <- sqrt(sum(w_lik2*w_lik2, na.rm = T))
  cos_lik/(w_lik1_norm*w_lik2_norm)
}

comp_style_by_id <- function(data = master, id1, id2){
  if(id1 == id2){
    return(1)
  }
  fam1 <- master %>% filter(p_id == id1) %>% pull(SMP.familiarity)
  fam2 <- master %>% filter(p_id == id2) %>% pull(SMP.familiarity)
  lik1 <- master %>% filter(p_id == id1) %>% pull(SMP.liking)
  lik2 <- master %>% filter(p_id == id2) %>% pull(SMP.liking)
  comp_style_pref(fam1, lik1, fam2, lik2)

}

get_sim_matrix <- function(data = master, style = "Heavy Metal", n_sample = 30, as_matrix = F, as_dist = F){
  set.seed(555)
  metal_fans <- data %>% filter(style ==  !!style, is_fan) %>% pull(p_id)
  data <- data %>% filter(p_id %in% metal_fans)
  ids <- unique(data$p_id) %>% sort()

  if(!is.null(n_sample)){
    ids <- sample(ids, n_sample) %>% sort()
  }
  ret <- map_dfr(1:(length(ids)), function(i){
    fam1 <- data %>% filter(p_id == ids[i]) %>% pull(SMP.familiarity)
    lik1 <- data %>% filter(p_id == ids[i]) %>% pull(SMP.liking)
    map_dfr((i):length(ids), function(j){
      if(j %% 100 == 0) messagef("i =%d, j = %d", i, j)
      lik2 <- data %>% filter(p_id == ids[j]) %>% pull(SMP.liking)
      fam2 <- data %>% filter(p_id == ids[j]) %>% pull(SMP.familiarity)
      ret <- comp_style_pref(fam1, lik1, fam2, lik2)
      tibble(id1 = ids[i], id2 = ids[j], sim  = ret)
    })
  })
  ret <- ret %>%
    bind_rows(ret %>% select(id2 = id1, id1 = id2, sim)) %>%
    distinct() %>%
    arrange(id1, id2)
  if(as_dist){
    ret <- ret %>% mutate(sim = 1 - sim)
  }
  if(as_matrix){
    ret <- ret$sim %>% matrix(ncol = length(ids), nrow = length(ids))
    row.names(ret) <- ids
    colnames(ret) <- ids

    if(as_dist){
      ret <- as.dist(ret)
    }
  }
  ret
}
