# Jedynak et al., 2020 analysis main results
# Paulina Jedynak
# 23/09/20

# Load packages, functions and plan ----
source("R/packages.R") # analysis packages
source("R/variables_lists.R") # manual lists of compounds, variables, etc
source("R/plan.R") # analysis plan
source("R/tables.R") # table preparations for final outputs
source("R/figures.R") # analysis figures

# # List and visualize outdated targets ----
drake::outdated(plan)
drake::vis_drake_graph(plan,
                       targets_only = TRUE,
                       hover = TRUE,
                       main = "Jedynak et al., 2020")

# Run plan ----
drake::make(plan)


