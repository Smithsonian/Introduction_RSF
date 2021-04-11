# Adding anthropogenic risk covariate to each dataframe by using extract_covariates function
wb15<-lapply(wb15, function(x) extract_covariates(x,anth_risk))

# Adding woody distance variable to each dataframe by using extract_covariates function
wb15<-lapply(wb15, function(x) extract_covariates(x,woody_dist))

# Adding animal ID column to each dataframe in a list using the ID we subseted above
wb15<- mapply(cbind, wb15, "Animal_ID"=Mara.An, SIMPLIFY=F)

# Binding all dataframe in a list into a single dataframe for all animals
wb15<-do.call(rbind,wb15)

# *********** Fitting the data with generalised linear mixed model ************
# ******************************************************************************

require(lme4)
mixed<-glmer(case_ ~ scale(anth_risk) + scale(woody_dist) + (1|Animal_ID), 
             family = binomial(), data =wb15)

# Getting summary of the model
summary(mixed)

# ******************** Response curves ******************************************
# ******************************************************************************

# Plotting response curve of the generalised linear mixed model with binomial error
# distribution error for making inference at the general population level.
# You can see that the curve has changed in comparison to when we fitted using a 
# Single individual

# Anthropogenic risk
effect_plot(mixed,pred=anth_risk,interval = TRUE,
            x.label = "Anthropogenic risk",
            y.label = "Probability of selection")+
  theme(panel.grid.major = element_blank())

# Woody density
effect_plot(mixed,pred=woody_dist,interval = TRUE,
            x.label = "Woody density",
            y.label = "Probability of selection")+
  theme(panel.grid.major = element_blank()) 

# ******** Generating predictions across a range of raster values **************
# ******************************************************************************
# First we need to stack the raster layers together using the stack function
stackras<-stack(anth_risk,woody_dist) 

# Then we use the predict function to generate predictions across a range of value
Predmix <- predict(stackras, mixed,type="response", re.form =NA)

# Plotting predictions

plot(predmix)

#or 
plot(predmix,col=rev(topo.colors(5)))