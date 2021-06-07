require("purrr")
require("tidyr")
require("tibble")
require("TwoSampleMR")

get_mv_exposures <- function(exposure_list, full_gwas_list, clump_exposures=FALSE) {
  
  # Here we are using a re-written source code of `mv_extract_exposures` function from 2SMR package.
  # It was neccesary to do it, as it only works with MRBase database input at the moment (but we have external exposures)
  
  
  # Get effects of each instrument from each exposure
  # all tophit SNPs in both exposures, (clumped optionally)
  exposures <- exposure_list %>% purrr::reduce(bind_rows) 

  if (clump_exposures) {
    # ***optional*** : clump exposures 
    temp <- exposures
    temp$id.exposure <- 1
    temp <- clump_data(temp)
    #temp <- clump_data_local(temp, local_path)
    exposures <- filter(exposures, SNP %in% temp$SNP)
  }
  
  
  # merge exposures (in 'outcomes' ~ full gwas format)
  # extract all instruments from full GWAS of exposures
  for (i in 1:length(full_gwas_list)){
    full_gwas_list[[i]] <- full_gwas_list[[i]] %>% filter(SNP %in% exposures$SNP)
  }
  d1 <- full_gwas_list %>%
    purrr::reduce(bind_rows) %>% 
    distinct()

  # get ids
  id_exposure <- unique(d1$id.outcome) 
  
  # convert first trait to exposure format  -- exp1 is exposure
  tmp_exposure <- d1 %>% filter(id.outcome == id_exposure[1]) %>% convert_outcome_to_exposure()
  # keep other traits as outcome -- exp2+ are outcomes
  tmp_outcome <- d1 %>% filter(id.outcome != id_exposure[1])
  
  # Harmonise against the first trait
  d <- harmonise_data(exposure_dat = tmp_exposure, 
                      outcome_dat = tmp_outcome, action=2)
  
  # Only keep SNPs that are present in all
  snps_not_in_all <- d %>% dplyr::count(SNP)  %>% 
                    filter(n < length(exposure_list)-1) %>%
                    pull(SNP)
  d <- filter(d, !SNP %in% snps_not_in_all)

  
  # Subset and concat data
  
  # for exp1 get exposure cols
  dh1x <- d %>% filter(id.outcome == id.outcome[1]) %>% 
    select(SNP, contains("exposure"))
  # for exp2 get outcome cols
  dh2x <-d %>%  select(SNP, contains("outcome"))
  # rename outcome to exposure in these
  names(dh2x) <- gsub("outcome", "exposure", names(dh2x) )
  # join together (drop not needed cols)
  exposure_dat <- bind_rows(dh1x, dh2x) %>%  
    select(-c("samplesize.exposure" ,"mr_keep.exposure", "pval_origin.exposure")) %>% 
    distinct()
  
  return(exposure_dat)
}



## NB: this is a shortcut from exposure_dat in 2SMR to exposures_joined (from manual steps)
# `exposure_dat` is from  `exposure_dat <- rbind.fill(dh1x, dh2x) ...` above, 
# it is used to generate equivant of `exposures_joined` from manual steps, but here
# it is called `exposures_joined_auto`

#  function to convert 2SMR format into MVMR format
make_mvmr_input <- function(exposure_dat, exposure_list, outcome.id.mrbase="", outcome.data=""){
  # specify one or the other outcome agrument (provide MRbase ID or supply outcome data)
  
  exposures_joined_auto <- exposure_dat %>%
    select(SNP, beta.exposure, se.exposure, exposure) %>% 
    # conver to wider format
    pivot_wider(names_from = exposure, values_from = c(beta.exposure, se.exposure)) 

  
  # Outcomes
  
  # extract SNPs for both exposures from outcome dataset 
  if (outcome.id.mrbase != "") {
    # if mrbase.id is given
    outcome_dat <- extract_outcome_data(snps = exposures_joined_auto$SNP,
                                        outcomes = outcome.id.mrbase)
  } else if (outcome.data != ""){
    # if outcome df is given
    outcome_dat <- outcome.data %>% filter(SNP %in% exposures_joined_auto$SNP)
  }
  
  # harmonize datasets 
  exposures <- exposure_list %>% purrr::reduce(bind_rows) 
  outcome_harmonised <- harmonise_data(exposures, outcome_dat)
  
  
  # Create variables for MV
  # remove factors structure in SNPs, and sort by SNP (YGs and XGs must have SNPs in the same order)
  YG <- outcome_harmonised %>% 
    select("SNP", "beta.outcome", "se.outcome") %>%
    distinct() %>%
    mutate(SNP = as.character(SNP)) %>%
    arrange((SNP)) 
  XGs <- exposures_joined_auto %>% 
    filter(SNP %in% outcome_harmonised$SNP) %>%
    mutate(SNP = as.character(SNP)) %>%
    arrange((SNP)) 
  
  # some checks
  stopifnot(dim(XGs)[1]==dim(YG)[1])
  unique(YG$SNP %in% XGs$SNP)
  unique(XGs$SNP %in% YG$SNP)
  all.equal(YG$SNP, XGs$SNP)
  
  return(list(YG = YG,
              XGs = XGs))
  
}


tidy_mvmr_output <- function(mvmr_res) {
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    rename(b=Estimate) %>% 
    rename(se="Std. Error") %>% 
    rename(pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}






