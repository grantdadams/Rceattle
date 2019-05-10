

# We can plot both runs as well:
mod_list <- list(ss_run, ms_run)
mod_names <- c("SS", "MS")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names, incl_proj = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE, incl_proj = TRUE)


plot_index(ms_run)
plot_catch(ms_run)

plot_srv_comp(ms_run)
plot_fsh_comp(ms_run)

