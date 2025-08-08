# Exploratory Analysis to examine feasibility of GWAS

exp_dat <- read.csv("data/experiment/eelgrass_mesocosm_start_end_complete.csv")[,-1]

View(exp_dat)

# how would this work? Compare all pops in a given treatment for relevant traits
# difficult to disentangle the effects of salinity and temperature
# any way to do a GWAS looking at salinity tolerance?
# could do it by using control (W coast) control (E coast)
# these are at the same temp but one much lower salinity

# what information would we gain through a GWAS?
#   genetic basis behind low salinity tolerance?
#   if we find genes or genotypes better suited to low salinity, 
#   could be interesting to examine how these genes are distributed across clones 
#   + whether they are more frequently found in E coast clones 
#   (indirectly looking at why we see higher clonality on E coast)

# another option for a project
# development of local adaptation package!!

# another option would be to estimate quantitative genetic values, but I don't know if this is feasible

# could also be really cool to compare SDMs to GOs for these baltic eelgrass!

# additionally would be cool to sample + sequence some eelgrass here in the US, 
# run similar stress experiment but potentially just with salinity