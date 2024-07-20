# Functions for use with the Supplementary Materials for Moran & Lester (in prep)

library(tidyverse)
library(lmerTest)
library(effects)

# Randomly sample one inventory per language

rand.inv = function(df){
  set.seed(1)
  
  nas =  df %>%
    distinct(InventoryID, ISO6393) %>%
    filter(is.na(ISO6393)) %>%
    pull(InventoryID)
  
  inventory_ids_sampled_one_per_isocode = df %>%
    distinct(InventoryID, ISO6393) %>%
    drop_na(ISO6393) %>%
    group_by(ISO6393) %>%
    sample_n(1) %>%
    pull(InventoryID)
  
  
  inventory_ids_sampled_one_per_isocode <- c(inventory_ids_sampled_one_per_isocode, nas)
  
  output = df %>%
    filter(InventoryID %in% inventory_ids_sampled_one_per_isocode) 
  
  return(output)
  
}

# Function to find the mode
Mode <- function(x){
  mode <- names(sort(table(inv_sizes), decreasing = T))[1]
  return(as.numeric(mode))
}

# Function to build the inventories for a single run

# df: (dataframe-like object) a table containing PHOIBLE (e.g., phoible.csv)
# k: (integer) the number of random inventories to generate 
# inv.inf: (boolean) should the database be built using information about inventory sizes from df (PHOIBLE)?
# freq.inf: (boolean) should the database be built using information about relative frequency of phonemes?
# repl: (boolean, default = FALSE) should the random sample be constructed using replacement?

inventory.builder <- function(k, inv.size.types, pho.typs, inv.inf, freq.inf, repl=F, inv.size.probs=NULL, pho.probs=NULL){
  
  # Open a list to store the inventories
  inv.list <- vector(mode="list", length=k)
  
  # Open a list to store the iteration labels
  lab.list <- vector(mode="list", length=k)
  
  # Create k randomized databases
  for(i in 1:k){
    
    # Define inventory size(s):
    
    # Sample an inventory size based on its probability
    if(inv.inf){
      inv.size <- as.numeric(sample(inv.size.types, 1, prob=inv.size.probs))
    }
    
    # ... or simply take the mode (if inventory size info is ablated)
    else{
      inv.size = Mode(inv.size.typs)
    }
    
    # Sample sounds based on their probabilities
    if(freq.inf){
      inv.list[[i]]  <- sample(pho.typs,  inv.size, replace=repl, prob=pho.probs)
    }
    
    # ... or simply draw a uniform random set from the full set of phonemes
    else{
      inv.list[[i]]  = sample(pho.typs, inv.size, replace=repl)
    }
    
    # Update the iteration label vectors
    lab.list[[i]] <- rep(i, inv.size)
  }
  
  # Generate the randomzed database
  out.df <- data.frame("InventoryID" = unlist(lab.list), 
                       "Phoneme" = unlist(inv.list))
  
  return(out.df)
}


# Function to construct, process, and collate random inventories

# df: (dataframe-like object) a table containing PHOIBLE (e.g., phoible.csv)
# n: (integer): the number of random databases to generate
# k: (integer) the number of random inventories to generate
# inv.inf: (boolean) should the database(s) be built using information about inventory sizes from df (PHOIBLE)?
# freq.inf: (boolean) should the database(s) be built using information about relative frequency of phonemes?
# repl: (boolean, default = FALSE) should the random sample be constructed using replacement?

processing.iterator = function(df, n, k, inv.inf, freq.inf, repl=F){
  
  # Get the inventory sizes per languages
  inv.sizes = df %>%
    group_by(InventoryID) %>%
    summarize(inv_size = n()) %>%
    pull(inv_size)
  
  # Find probability distributions and sets for
  # inventory sizes and phonemes in df
  
  # .. if inventory information is preserved
  if(inv.inf){
    inv.size.probs = prop.table(table(inv.sizes))
    inv.size.typs = names(inv.size.probs)
  }
  
  # or just take the unique set of inventory sizes
  else{
    inv.probs = NULL
    inv.size.typs = unique(inv.sizes)
  }
  
  # ... if frequency information is preserved
  if(freq.inf){
    pho.probs = prop.table(table(df$Phoneme))
    pho.typs = names(pho.probs)
  }
  
  # or just take the unique set of phonemes
  else{
    pho.probs = NULL
    pho.typs = unique(df$Phoneme)
  }
  
  # Open a list to store the databases
  df.list <- vector(mode="list", length = n)
  
  # For as many n desired databases (copies of "full PHOIBLE" or beyond)
  #for(i in 1:n){
  df.list = lapply(1:n, function(x)
    inventory.builder(k, 
                      inv.size.typs, 
                      pho.typs, 
                      inv.inf, 
                      freq.inf, 
                      repl, 
                      inv.size.probs, 
                      pho.probs) %>%
      mutate(Iteration = x) %>%
      get_cumulative_segment_counts_df(.))
  
  # Combine all of the individual random databases
  all.dfs <- do.call(rbind, df.list)
  
  return(all.dfs)
}

# Arguments 
# - df: dataframe of segments (i.e., PHOIBLE or UPSID)
get_cumulative_segment_counts_df <- function(df) {
  output.df = df %>% 
    group_by(Iteration) %>%
    mutate(results = lengths(Reduce(function(x, y) unique(c(x, y)), Phoneme, acc = T))) %>%
    group_by(Iteration, InventoryID) %>%
    summarize(results = max(results)) %>%
    arrange(Iteration, as.numeric(results)) %>%
    group_by(Iteration) %>%
    mutate(InventoryID = seq(1, n()),
           Iteration = as.factor(Iteration))
  return(output.df)
}


## Function to plot density of inventory sizes

# phoible: (dataframe-like object) a dataframe containing PHOIBLE-like data (by segment)
# random.pho: (dataframe-like object) a dataframe containing a randomized version of PHOIBLE

inv.density.plotting = function(phoible, random.pho){
  pho.den = phoible %>%
    group_by(InventoryID) %>%
    summarize(sizes = n()) %>%
    mutate(Source = "PHOIBLE",
           InventoryID = as.factor(InventoryID))
  
  rand.den = random.pho %>%
    group_by(InventoryID) %>%
    summarize(sizes = n()) %>%
    mutate(Source = "Random",
           InventoryID = as.factor(InventoryID))
  
  den.dat = bind_rows(pho.den, rand.den)
  
  p = ggplot(den.dat, aes(x=sizes, fill = Source)) +
    geom_density(alpha=0.4) +
    theme_bw() + 
    ggtitle("Inventory size density")
  
  return(p)
}



## Compute the hypothetical Heaps' estimates based on the constructed databases


# df: (dataframe-like object) randomized form of PHOIBLE
heaps.preds = function(random.dataframe, mod){
  pho.inv.num = length(unique(random.dataframe$InventoryID))
  
  # Extract the predictions of the model with 95% confidence intervals
  results = as.data.frame(allEffects(mod, xlevels=pho.inv.num)[[1]]) %>%
    mutate(fit = exp(fit),
           lower = exp(lower),
           upper = exp(upper))
  
  
  # Create a label for the source of these estimates
  Source = as.factor(rep("Hyp. Heaps'", nrow(results)))
  
  # Create a dataframe that can be combined with the
  # "empirical" databases
  pred.df = data.frame("results" = results$fit, "InventoryID" = results$InventoryID, "Source" = Source, "Iteration" = "Fitted")
  
  return(pred.df)
  
}



## Function to combine "empirical" (true or random) and modeled results

# pho: (dataframe-like object) PHOIBLE or randomized version thereof
# modeled: (dataframe-like object) Herdan-Heaps' predictions based on data
plot.emp.vs.pred = function(pho.rand, modeled, pho.true) {
  
  # Combine the dataframes
  plot.dat <- bind_rows(pho.rand, modeled)
  
  # Remove repeated values from true data
  pho.true <- pho.true[!duplicated(pho.true$results),]
  
  # Plot empirical against modeled cumulative type counts
  p <- ggplot(plot.dat %>% filter(Iteration != "Fitted"), aes(x = InventoryID, y = results, color=Iteration)) +
    geom_line(lwd = 1.2, alpha=0.5) +
    geom_line(data=plot.dat %>% filter(Iteration == "Fitted"), aes(x=InventoryID, y=results), color = "black", lty=2) + 
    geom_line(data=pho.true, aes(x = InventoryID, y = results), color = "dodgerblue", alpha = 0.5, lwd = 1.5) +
    theme_bw()
  
  return(p)
}



exploratory.density.plots = function(df) {
  # How many sounds added per each additional language?
  p1 = ggplot(df, aes(x = InventoryID, y = diffs, color=Iteration)) +
    geom_line(size=0.1, alpha=0.5) +
    geom_smooth(method="gam", formula=y ~ s(x), color="black", size=.4, lty=2) +
    theme_bw() +
    ggtitle("Number of new phonemes per increment") 
  
  # Density of new sounds added per language
  p2 = ggplot(test, aes(x = diffs, color=Iteration, fill=Iteration)) +
    geom_density(alpha=0.2) +
    theme_bw() +
    xlab("Difference in new phonemes (lag 1)") + 
    ggtitle("Change in new phonemes per language") 
  
  return(c(p1, p2))
}





# Function to get frequencies and frequency ranks
freqs = function(df){
  out.df = df %>%
    group_by(Phoneme) %>%
    summarize(freq=n()) %>%
    arrange(-freq) %>%
    mutate(rank = seq(1, nrow(.)))
  return(out.df) 
}




## Plotting log-log fo random against true phoible

loglogplot = function(df){
  
  # Compute frequencies and ranks for PHOIBLE
  pho.freqs <- freqs(df) %>% mutate(Source="PHOIBLE")
  
  # Compute phoneme frequencies and ranks for the random database
  # and log them.
  rand.freqs <- freqs(rand_test) %>%
    mutate(Source = "Random",
           logFreq = log(freq),
           logRank = log(rank))
  
  # Combine with the PHOIBLE data
  freq.rank.dat <- bind_rows(rand.freqs,
                             pho.freqs %>%
                               mutate(logFreq = log(freq),
                                      logRank = log(rank)))
  
  # Plot the log-log (freq x rank) for true and random databases
  p <- ggplot(freq.rank.dat, aes(y = logFreq, x=logRank, color=Source)) +
    geom_line() +
    theme_bw()
  
  return(p)
}


## Function to compare and plot cumulative vocabulary growth rates

# pho.proc: (dataframe-like object) PHOIBLE with cumulative counts over k iterations of randomization
# rand.proc: (dataframe-like object) random samples of PHOIBLE with cumulative counts over k iterations of randomization
# rand.mod: (dataframe-like object) Herdan-Heaps' model predictions based on random sample (rand.proc)
rand.v.pho.plot = function(pho.proc, rand.proc, rand.mod){
  # Combine the PHOIBLE, random, and fitted-random data together
  rand.pho.dat <- bind_rows(pho.proc %>%
                              mutate(Source = "True"),
                            rand.proc,
                            rand.mod)
  
  # Just for the plot (reordering factor levels)
  rand.pho.dat$Source = as.factor(rand.pho.dat$Source)
  rand.pho.dat$Source = factor(rand.pho.dat$Source, levels = levels(rand.pho.dat$Source)[c(3, 2, 1)])
  
  p = ggplot(rand.pho.dat %>% 
               filter(Source!="Hyp. Heaps (Random)'") %>% 
               droplevels(.), 
             aes(y=results, x=InventoryID, color=Source)) +
    geom_line() +
    geom_smooth(data=rand.pho.dat %>%
                  filter(Source=="Random") %>%
                  droplevels(.),
                aes(y=results, x=InventoryID), color="darkgreen") +
    geom_line(data=rand.mod, aes(y=results, x=InventoryID), color="black", size=.8, lty=3) + 
    theme_bw() +
    ggtitle('Random (+sizes, +Zipf) vs. True PHOIBLE')
  
  return(p)
}

# Function to randomly segment PHOIBLE into training and testing samples, make predictions
cross.pred = function(df, n){
  
  df.list = lapply(1:n, function(x){
    
    set.seed(x)
    
    # Randomize PHOIBLE
    rand.invs = unique(df$InventoryID) %>% sample()
    
    # Training inventories
    train.invs = rand.invs[1:floor(.75*length(rand.invs))]
    
    # Test inventories
    test.invs = rand.invs[!rand.invs %in% train.invs]
    
    # Take first 75% of inventories as training data
    train.df = df %>% 
      filter(InventoryID %in% train.invs) %>% 
      arrange(match(InventoryID, train.invs)) %>%
      mutate(Iteration = n)
    
    # Take remaining 25% of inventories as test data
    test.df = df %>% 
      filter(InventoryID %in% test.invs) %>% 
      arrange(match(InventoryID, test.invs)) %>%
      mutate(Iteration = n)
    
    # Compute the cumulative frequencies for both
    train.df = get_cumulative_segment_counts_df(train.df) %>% mutate(InventoryID = seq(1, n()))
    
    test.df = get_cumulative_segment_counts_df(test.df) %>% mutate(InventoryID = seq(1, n()))
    
    # Train a model on the training data
    train.mod = lm(log(results) ~ log(InventoryID), train.df)
    
    # Get model predictions for test data
    preds = exp(predict(train.mod, newdata = test.df))
    
    out.df = data.frame(Iteration = x,
                        Empirical = test.df$results,
                        Predicted = preds,
                        InventoryID = seq(1:length(preds)))
    
    return(out.df)})
  
  return(do.call(rbind, df.list) %>% mutate(Iteration=as.factor(Iteration)))
}

# Function to plot output of cross.pred()
# df: (dataframe-like object) should be the output of cross.pred() function (see above)
# title: (character string) title for the plot
# plot.typ: (character string) one of two options:
# - plot.typ = "true.vs.pred" plots the correlation between true and predicted values
# - plot.typ = "heaps.pred" plots the predicted values
# - plot.typ = "heaps.emp" plots the fit of Herdan-Heaps' within the training set
cross.pred.plot = function(cross.pred.dat, mod.obj, title, plot.typ = "true.vs.pred"){
  
  if(plot.typ=="true.vs.pred"){
    
    # Model the relationship between predicted and empirical values
    
    # Extract the coefficient
    coeff = paste0("β = ", round(summary(mod.obj)$coefficients[2,1], 2))
    
    # Extract the p-value for Empirical
    pval = paste0("p = ", anova(mod.obj) %>% data.frame() %>% select(6) %>% pull() %>% round(2))
    
    # Combine the coefficient and p-value text, separated onto two lines, for plotting
    stat.sum = paste(coeff, pval, sep="\n")
    
    # Extract the model predictions
    eff.df = effect("Empirical", mod.obj, xlevels=100) %>% 
      data.frame() %>%
      rename(Predicted = fit) %>%
      mutate(Iteration = "Predictions")
    
    # Plotting
    p = ggplot(data=cross.pred.dat, aes(x=Empirical, y=Predicted, group = Iteration)) +
      geom_abline(slope=1, intercept=0, color="red", lwd=1.1) + 
      geom_line(alpha=0.15, lwd=1) +
      geom_line(data=eff.df, aes(x=Empirical, y=Predicted),  color="black", lwd=1.1) + 
      geom_ribbon(data=eff.df, aes(ymin=lower, ymax=upper), color = "darkblue", fill="dodgerblue", alpha=0.2, lty=2) +
      geom_text(data=eff.df, aes(x = max(Empirical)/4, y = (3*max(Predicted))/4, label=stat.sum), size=5, family="serif") +
      ggtitle(title) +
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5))
  }
  
  else if(plot.typ=="heaps.pred"){
    p = ggplot(data=cross.pred.dat, aes(x=InventoryID, y=Predicted, color=Iteration)) +
      geom_line(alpha=0.5) +
      geom_line(data=cross.pred.dat,
                aes(x=InventoryID, y=Empirical, color=Iteration), lty=2) + 
      ggtitle(title) +
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5)) + 
      guides(color = "none") + 
      xlab("Inventory") +
      ylab("Type count")
  }
  
  else if(plot.typ=="heaps.emp"){

    emp.eff.df = effect("log(InventoryID)", mod.obj, xlevels=100)
    
    emp.eff.df = emp.eff.df %>%
      data.frame() %>% 
      rename(Predicted = fit) %>% 
      mutate(Predicted = exp(Predicted),
             lower = exp(lower),
             upper = exp(upper))
    
    # Extract the coefficient
    coeff = paste0("β = ", round(exp(summary(mod.obj)$coefficients[2,1]), 2))
    
    # Extract the p-value for Empirical
    pval = paste0("p = ", anova(mod.obj) %>% data.frame() %>% select(6) %>% pull() %>% round(2))
    
    # Combine the coefficient and p-value text, separated onto two lines, for plotting
    stat.sum = paste(coeff, pval, sep="\n")
    
    p = ggplot(data=emp.eff.df, aes(x=InventoryID, y=Predicted)) +
      geom_line(data=cross.pred.dat, aes(x=InventoryID, y=Empirical, group=Iteration), color="grey", alpha=0.8) + 
      geom_line() +
      geom_ribbon(aes(ymin=lower, ymax=upper), color = "darkblue", fill="dodgerblue", alpha=0.2, lty=2) +
      geom_text(data=emp.eff.df, aes(x = max(InventoryID)/6, y = (5*max(Predicted))/6, label=stat.sum), size=5, family="serif") +
      ggtitle(title) +
      theme_bw() +
      theme(plot.title = element_text(hjust=0.5)) + 
      guides(color = "none") +
      xlab("Inventory")
    
  }
  
  else{
    stop("The phon.typ argument should match one of the following: 'true.vs.pred', 'heaps.pred', or 'heaps.emp'!")
  }
  
  return(p)
  
}

# Function to apply cross-prediction based on population
# df: (dataframe-like object) PHOIBLE
# n: (numeric) number of samples & models to generate
# status: (character) status(es) to select for training sample (e.g., "robust")
# pop.cut: (numeric) number indicating which quantile(s) to use to divide testing from training
#                 ex. 2 = 1st quantile, 3 = 2nd quantile, ...
# down.sample: (boolean) if true, all population categories will be downsampled to match the least frequent

pop.cross.pred = function(df, n, status=NULL, pop.cut=NULL, down.sample=F){
  
  # First, remove any NAs
  df = df %>% 
    filter(!is.na(simp.status)) 
  
  # Downsizing sample of each status to the smallest category
  if(down.sample){
    max.size = df %>%
      select(simp.status, InventoryID) %>%
      distinct() %>%
      group_by(simp.status) %>%
      summarize(count = n()) %>%
      pull(count) %>%
      min()
    
    inv.names = df %>%
      select(InventoryID, simp.status) %>%
      distinct() %>%
      group_by(simp.status) %>%
      slice_sample(n=max.size) 
    
    df = df %>%
      filter(InventoryID %in% inv.names$InventoryID)
  }
  
  pops = df %>%
    select(InventoryID, Pop.Size.Clean) %>%
    distinct() %>%
    select(Pop.Size.Clean) %>%
    pull()
  
  if(!is.null(pop.cut)){
    cut.point = quantile(pops)[pop.cut]
  }
  
  df.list = lapply(1:n, function(x){
    
    set.seed(x)
    
    # Randomize PHOIBLE
    rand.invs = unique(df$InventoryID) %>% sample()
    
    if(!is.null(status)){
      train.df = df %>%
        filter(simp.status %in% status) %>%
        arrange(match(InventoryID, rand.invs)) %>%
        mutate(Iteration = n)
      
      test.df = df %>% 
        filter(!simp.status %in% status) %>% 
        arrange(match(InventoryID, rand.invs)) %>%
        mutate(Iteration = n)
    }
    
    else if (!is.null(pop.cut)){
      
      # Take first X% of inventories as training data
      train.df = df %>% 
        filter(Pop.Size.Clean >= cut.point) %>% 
        arrange(match(InventoryID, rand.invs)) %>%
        mutate(Iteration = n)
      
      # Take first 1-X% of inventories as test data
      test.df = df %>% 
        filter(Pop.Size.Clean < cut.point) %>% 
        arrange(match(InventoryID, rand.invs)) %>%
        mutate(Iteration = n)
    }
    
    # Compute the cumulative frequencies for both
    train.df = get_cumulative_segment_counts_df(train.df) %>%
      mutate(InventoryID = seq(1, n()))
    
    test.df = get_cumulative_segment_counts_df(test.df) %>% 
      mutate(InventoryID = seq(1, n()))
    
    # Train a model on the training data
    train.mod = lm(log(results) ~ log(InventoryID), train.df)
    
    # Get model predictions for test data
    preds = exp(predict(train.mod, newdata = test.df))
    
    out.df = data.frame(Iteration = x,
                        Empirical = test.df$results,
                        Predicted = preds,
                        InventoryID = seq(1:length(preds)))
    
    return(out.df)})
  
  return(do.call(rbind, df.list) %>% mutate(Iteration=as.factor(Iteration)))
}


# Category walk function
# df: (dataframe-like object) PHOIBLE database with population info added
cat.walk = function(df, n){
  df.list = lapply(1:n, function(x){
    rand.order = df %>%
      filter(!is.na(simp.status)) %>%
      select(simp.status, ISO6393) %>%
      group_by(simp.status) %>%
      distinct() %>%
      sample_n(size=n()) %>%
      arrange(desc(simp.status)) %>%
      pull(ISO6393)
    
    df %>%
      arrange(match(ISO6393, rand.order)) %>%
      mutate(Iteration = x) %>%
      get_cumulative_segment_counts_df() %>%
      rename(Empirical = results)
  })
  
  output = do.call(rbind, df.list)
  
  preds = lmer(log(Empirical) ~ log(InventoryID) + (1|Iteration), data=output) %>% predict() %>% exp()
  
  output$Predicted = preds
  
  output = output[,c(1, 3, 4, 2)]
  
  return(output)
}

# Cross predicting within population class
cross.pred.cat = function(df, n, prop){
  df.list = lapply(1:n, function(x){
    
    df = df %>%
      mutate(Iteration = x)
    
    rand.order = df %>%
      filter(!is.na(simp.status)) %>%
      select(simp.status, ISO6393) %>%
      group_by(simp.status) %>%
      distinct() %>%
      sample_n(size=n()) %>%
      arrange(desc(simp.status))
    
    train.invs = rand.order %>%
      group_by(simp.status) %>%
      sample_n(size=floor(prop*n())) %>%
      arrange(desc(simp.status))
    
    train.df= df %>%
      filter(ISO6393 %in% train.invs$ISO6393) %>%
      get_cumulative_segment_counts_df()
    
    test.df = df %>% 
      filter(!ISO6393 %in% train.invs$ISO6393) %>%
      get_cumulative_segment_counts_df()
    
    # Train a model on the training data
    train.mod = lm(log(results) ~ log(InventoryID), train.df)
    
    # Get model predictions for test data
    preds = exp(predict(train.mod, newdata = test.df))
    
    out.df = data.frame(Iteration = as.factor(x),
                        Empirical = test.df$results,
                        Predicted = preds,
                        InventoryID = seq(1:length(preds)))
    
    return(out.df)})
  
  return(do.call(rbind, df.list))     
}


# Functoin to predict across classes, walking categrory by category, with diminishing population per category
cross.pred.cat.dim = function(df, n){
  
  df.list = lapply(1:n, function(x){
    
    set.seed(x)
    
    df = df %>%
      mutate(Iteration = x)
    
    rand.order = df %>%
      filter(!is.na(simp.status)) %>%
      select(simp.status, ISO6393) %>%
      distinct() %>%
      group_by(simp.status) %>%
      sample_n(size=n()) %>%
      arrange(desc(simp.status))
    
    rand.df = df %>%
      arrange(match(ISO6393, rand.order$ISO6393))
    
    sizes = as.data.frame(table(rand.order$simp.status)) %>%
      rename(simp.status = Var1) %>%
      arrange(desc(simp.status)) %>%
      mutate(sample.size = floor(Freq / seq(1:n()))) 
    
    train.invs = rand.order %>% 
      left_join(sizes, by = "simp.status") %>%
      group_by(simp.status) %>%
      mutate(samp = sample(n())) %>%
      filter(samp <= sample.size) %>%
      ungroup()
    
    train.df = rand.df %>%
      filter(ISO6393 %in% train.invs$ISO6393) %>%
      get_cumulative_segment_counts_df()
    
    test.df = rand.df %>%
      filter(! ISO6393 %in% train.invs$ISO6393) %>%
      get_cumulative_segment_counts_df()
    
    
    # Train a model on the training data
    train.mod = lm(log(results) ~ log(InventoryID), train.df)
    
    # Get model predictions for test data
    preds = exp(predict(train.mod, newdata = test.df))
    
    out.df = data.frame(Iteration = as.factor(x),
                        Empirical = test.df$results,
                        Predicted = preds,
                        InventoryID = seq(1:length(preds)))
    
    return(out.df)})
  
  output = do.call(rbind, df.list)
  
  return(output)
  
}
