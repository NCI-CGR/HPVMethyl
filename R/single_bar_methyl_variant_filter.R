single_bar_methyl_variant_filter <- function(variants, filteringTablePath, posConversionTable, pn_filters, manifest, control_defs,  control_freq_defs){
  
  require(fuzzyjoin)
  
  manifest <- manifest %>% filter(!is.na(Owner_Sample_ID))

  manifest %>%
  transform(BC1 = as.character(BC1)) %>% 
    rename(barcode = BC1) -> manifest
  
  
  filteringTable = read_tsv(filteringTablePath) %>%
    map_if(is.factor, as.character) %>%
    as_tibble() %>%
    rename(CHROM = Chr, POS = Base_num, REF = Base_ID, ALT = vcf_variant)
  
  #hotspot vars...
  
  GA_variants = variants %>% 
    filter(HS) %>%
    filter(!(ALT %in% c("C", "T"))) %>%
    transform(barcode = as.character(barcode)) %>%
    glimpse()

  manifest %>%
    inner_join(GA_variants) %>%
    select(-filename) %>%
    filter(!(is.na(Owner_Sample_ID))) %>%
    write_csv("lineage_variants_results.csv")
  
  pos_conversion = read_tsv(posConversionTable) %>%
    map_if(is.factor, as.character) %>%
    as_tibble()
  
  filtered_variants = variants %>% 
    filter(ALT %in% c("C", "T")) %>%  
    inner_join(filteringTable) %>% 
    transform(barcode = as.character(barcode)) %>%
    mutate(AF = as.double(AF)) %>%
    mutate(qc_reason = "") %>%
    mutate(qc_reason = ifelse(DP >= min_DP, qc_reason,
                              "min_DP")) %>% 
 #   filter(!(qc_reason == "min_DP")) %>%
    mutate(methyl_freq = case_when(REF == "T" & ALT == "C" ~ AF, REF == "C" & ALT == "T" ~ FRO/FDP)) %>%
    mutate(qc_reason = ifelse(SRF >= min_coverage_pos, qc_reason,
                              paste0(qc_reason, ";", "min_coverage_pos"))) %>%
    mutate(qc_reason = ifelse(SRR >= min_coverage_neg, qc_reason,
                              paste0(qc_reason, ";", "min_coverage_neg"))) %>%
    mutate(qc_reason = ifelse(SAF >= min_allele_coverage_pos, qc_reason,
                              paste0(qc_reason, ";", "min_allele_coverage_pos"))) %>%
    mutate(qc_reason = ifelse(SAR >= min_allele_coverage_neg, qc_reason,
                              paste0(qc_reason, ";", "min_allele_coverage_neg"))) %>%
    mutate(qc_reason = ifelse(QUAL >= min_qual, qc_reason,
                              paste0(qc_reason, ";", "min_qual"))) %>%
    mutate(qc_reason = ifelse(STB <= max_alt_strand_bias, qc_reason,
                              paste0(qc_reason, ";", "max_alt_strand_bias"))) %>%
    mutate(qc_reason = ifelse(methyl_freq >= min_freq, qc_reason,
                              paste0(qc_reason, ";", "min_freq"))) %>%
    mutate(qc_reason = ifelse(methyl_freq <= max_freq, qc_reason,
                              paste0(qc_reason, ";", "max_freq"))) %>%
    mutate(qc_reason = ifelse(FILTER == "PASS", qc_reason,
                              paste0(qc_reason, ";", FILTER))) %>%
    mutate(qc_reason = ifelse(qc_reason == "","Pass",paste0("Fail",";",qc_reason))) %>%
    mutate(status = ifelse(qc_reason == "Pass", "Pass", "Fail")) %>%
    rename(pos_amplicon = POS, chr_amplicon = CHROM) %>% 
    glimpse() %>% 
    inner_join(pos_conversion, by = c("chr_amplicon","pos_amplicon")) %>% 
    mutate(methyl_freq = ifelse(status == "Fail","NA",methyl_freq)) %>%
    select(chr_amplicon, pos_amplicon, pos_from_conversion_table = pos, DP, methyl_freq, QUAL, status, qc_reason, everything()) 
  
  pre_return_table = manifest %>% 
    left_join(filtered_variants %>% select(chr_amplicon,pos_amplicon,chr,pos_from_conversion_table,REF,ALT,DP,methyl_freq,status,qc_reason,QUAL,FILTER,CpG_Variant_Info,HS,TYPE,HRUN,everything())) %>% 
    filter(!(is.na(Owner_Sample_ID))) 
  
  ### focus on the depth (nothing to do with mehtyl_freq)
  coverage_matrix = pre_return_table %>% 
    group_by(Owner_Sample_ID, barcode, chr_amplicon) %>% 
    summarize(depth = max(DP)) %>% 
    ungroup() %>% 
   # mutate(chr_amplicon = str_replace_all(chr_amplicon,"-","_")) %>% 
    mutate(MASIC_reads = ifelse(str_detect(chr_amplicon,"MASIC*"),depth,0 )) %>% 
    mutate(HPV_reads = ifelse(str_detect(chr_amplicon,"HPV*"),depth,0 )) %>%
    mutate(Human_reads = ifelse(!str_detect(chr_amplicon,"MASIC") & !str_detect(chr_amplicon,"HPV*"),depth,0 )) %>% 
    group_by(barcode, Owner_Sample_ID) %>%
    mutate(total_MASIC_reads = sum(MASIC_reads),total_HPV_reads = sum(HPV_reads),total_human_reads = sum(Human_reads)) %>%
    select(-MASIC_reads,-HPV_reads,-Human_reads) %>%
    spread(chr_amplicon, depth) %>% 
    glimpse() 
  

  ### Get the column order from control_defs and use it to reorder coverage_majtrix
  control_defs %>%
    select(-control_code,-control_type) %>%
    colnames() -> con_sort_order
  
  coverage_matrix = coverage_matrix %>%
    select(Owner_Sample_ID, barcode, total_HPV_reads,total_human_reads, total_MASIC_reads, con_sort_order, everything())
  
  ### simple join manifest table to have the output for coverage_matrix_results.csv
  manifest %>%
    inner_join(coverage_matrix) %>%
    filter(!(is.na(Owner_Sample_ID))) %>%
    write_csv("coverage_matrix_results.csv") 
  
  
# Prepare the pn_filter
# chrom Min_reads_per_type
  pn_filter<-read.csv(pn_filters)
  pn_filter %>%
    rename(chrom = contig) -> pn_filters
 
# coverage_matrix + pn_filters => detailed_pn_matirx
# pos status is assigned for each amplicon.
 coverage_matrix %>% 
    gather(chrom,depth,-Owner_Sample_ID,-barcode,-total_HPV_reads,-total_MASIC_reads,-total_human_reads) %>%
    inner_join(pn_filters) %>%
    mutate(status = ifelse(depth >= Min_reads_per_type,"pos","neg")) %>%
    select(-Min_reads_per_type,-depth) %>% 
    spread(chrom,status) -> detailed_pn_matrix
 
#Calculate num_type_pos separately (the code was not well-written) 
# to add Num_Types_Pos 
    
    num_type_list <- detailed_pn_matrix %>%
    tidyr::gather("type_id", "type_status", starts_with("HPV")) %>%  #Only counting HPV contigs and not MASIC contigs here
      group_by(Owner_Sample_ID,barcode) %>% 
      mutate(count_pos = if_else(type_status == "pos", 1, 0)) %>%
      separate(type_id, into = c("grouped_type","extra"),sep = "_", remove = F) %>% 
      group_by(barcode,grouped_type) %>%
      mutate(total_pos_per_type = sum(count_pos)) %>% 
      ungroup() %>%
      mutate(new_status = ifelse(total_pos_per_type >= 2, "pos","neg")) %>%
      mutate(new_count = ifelse(new_status == "pos",1,0)) %>% 
      select(-type_id,-type_status, -extra, -count_pos) %>%
      unique() %>%
      group_by(Owner_Sample_ID,barcode) %>%
      mutate(Num_Types_Pos = sum(new_count)) %>% 
      ungroup() %>% 
      select(-grouped_type,-new_status,-new_count,-total_pos_per_type) %>%
      select(Owner_Sample_ID, barcode,Num_Types_Pos,everything()) %>%
      unique() %>% 
      inner_join(detailed_pn_matrix) %>% 
      transform(barcode = as.character(barcode)) 
    
    ### to add Num_Types_Pos to detailed_pn_matrix_results
    manifest %>%
      inner_join(num_type_list) %>% 
      filter(!(is.na(Owner_Sample_ID))) %>%
      write.csv("detailed_pn_matrix_results.csv", row.names = F)
       
    
    

  # Simple pn matrix
    
    simple_pn_matrix_long = detailed_pn_matrix %>%
      filter(!is.na(Owner_Sample_ID)) %>%
      gather("CHROM", "status", starts_with("HPV"), factor_key = TRUE) %>%
      separate(CHROM, sep = "_", into = c("type"), remove = FALSE, extra = "drop") %>%
      glimpse() %>% 
      mutate(status_as_integer = ifelse(status == "pos", 1, 0)) %>%
      group_by(barcode, type) %>%
      mutate(sum_status = sum(status_as_integer)) %>%
      mutate(simple_status = ifelse(sum_status >= 2, "pos", "neg")) %>%
      ungroup() %>%
      glimpse() %>%
      select(-status_as_integer, -CHROM, -sum_status, -status) %>%
      distinct() %>%
      glimpse() %>%
      group_by(barcode, type) %>%
      distinct() 
    
  simple_pn_matrix = simple_pn_matrix_long %>%
    spread(type, simple_status)
    
  ### return table is changed to pre_return_table and tvr_table (for the final return)
  ### update qc_reason and status too
    
    tvr_table = pre_return_table %>%
      inner_join(num_type_list %>% select(Owner_Sample_ID, barcode, Num_Types_Pos) %>% transform(Num_Types_Pos = as.integer(Num_Types_Pos)),by = c("Owner_Sample_ID","barcode")) %>% 
      left_join( simple_pn_matrix_long %>% ungroup %>% select(barcode, Gene=type, simple_status) )  %>%
      mutate(qc_reason = ifelse(!is.na(simple_status) & simple_status == 'neg' , ifelse(methyl_freq %in% c("0", "NA"), ifelse(qc_reason=="Pass", "Fail", qc_reason), "Fail;false-positive"),qc_reason)) %>%
      mutate(methyl_freq = ifelse(!is.na(simple_status) & simple_status == 'neg', "NA",methyl_freq)) %>%
      mutate(status = ifelse(!is.na(simple_status) & simple_status == 'neg', "Fail",status))

    write_csv(tvr_table,"target_variants_internal_results.csv")

  # to have a count table for masic1/2: Owner_Sample_ID     barcode num_MASIC_pos num_MASIC2_pos
  simple_pn_matrix_long %>%
      spread(type, simple_status) %>% 
      gather("MAS_type","value", starts_with("MASIC")) %>%
      separate(MAS_type, c("MAS","a")) %>%
    #  mutate(HPV_pos = ifelse(simple_status == "pos",1,0)) %>%
      mutate(masic1_total = ifelse(MAS == "MASIC1" & value == "pos",1,0)) %>% 
      mutate(masic2_total = ifelse(MAS == "MASIC2" & value == "pos",1,0)) %>%
      group_by(barcode) %>%
      mutate(num_MASIC_pos = sum(masic1_total)) %>%
      mutate(num_MASIC2_pos = sum(masic2_total)) %>%
      ungroup() %>%
      select(Owner_Sample_ID,barcode,num_MASIC_pos,num_MASIC2_pos) %>%
      distinct() -> count_table

  # to have HPV_count_table: Owner_Sample_ID,barcode,num_HPV_pos  
  simple_pn_matrix_long %>%
      mutate(HPV_pos = ifelse(simple_status == "pos",1,0)) %>%
      group_by(Owner_Sample_ID,barcode) %>%
      mutate(num_HPV_pos = sum(HPV_pos)) %>%
      select(Owner_Sample_ID,barcode,num_HPV_pos) %>%
      distinct() -> HPV_count_table
      
  ### to write simple_pn_matrix_results.csv 
  manifest %>%
    inner_join(count_table) %>%
    inner_join(HPV_count_table) %>% 
    inner_join(simple_pn_matrix) %>%
    filter(!(is.na(Owner_Sample_ID))) %>% 
     write.csv("simple_pn_matrix_results.csv", row.names = F)
    
  
  #freq_matrix
  # calculate mean methyl_freq of qualified sites in each amplicon
  # saved in a wide from
  freq_matrix = tvr_table %>% 
    transform(methyl_freq = as.numeric(methyl_freq)) %>%
    group_by(Owner_Sample_ID, barcode, chr_amplicon) %>%
    summarize(mean_freq = mean(methyl_freq[status == "Pass"], na.rm=T)) %>% 
    ungroup() %>% 
    group_by(barcode, Owner_Sample_ID) %>% 
    spread(chr_amplicon, mean_freq) %>% 
    select(Owner_Sample_ID,barcode,con_sort_order,everything()) %>%
    glimpse() 
  
  manifest %>%
    inner_join(freq_matrix ) %>% 
    filter(!(is.na(Owner_Sample_ID))) %>%
    write_csv("freq_matrix_results.csv") 
  
  ### Convert control_defs from wide form to the long form
  # with the 4 columns: 
# $ control_code <chr> "HEK", "human", "hg19_all-HR", "NTC", "HEK", "human", ...
# $ control_type <chr> "pos", "pos", "pos", "neg", "pos", "pos", "pos", "neg"...
# $ chrom        <chr> "MASIC1_A", "MASIC1_A", "MASIC1_A", "MASIC1_A", "MASIC...
# $ min_coverage <dbl> 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20...
  control_defs = control_defs %>%
    glimpse() %>%
    tidyr::gather("chrom", "min_coverage", -control_code,-control_type) %>%
    glimpse()
  
#   control_results[,1:6]
# # A tibble: 4 x 6
# # Groups:   Owner_Sample_ID, barcode [4]
#   Owner_Sample_ID     barcode num_pass_targets num_fail_targets ASCL1 EPB41L3
#   <chr>               <chr>              <dbl>            <dbl> <chr> <chr>  
# 1 0perc-hg19_all-HR   259                   55                2 Pass  Pass   
# 2 100perc-hg19_all-HR 229                   57                0 Pass  Pass   
# 3 50perc-hg19_all-HR  235                   54                3 Pass  Pass   
# 4 HEK                 287                   57                0 Pass  Pass   
  control_results = coverage_matrix %>% 
    tidyr::gather("chrom", "depth", -Owner_Sample_ID, -barcode, -starts_with("MASIC"), -total_MASIC_reads, -total_HPV_reads, -total_human_reads) %>% 
    mutate(depth = as.integer(depth)) %>%
    fuzzyjoin::fuzzy_join(control_defs, mode = "inner", by = c("Owner_Sample_ID" = "control_code"), match_fun = function(x, y) str_detect(x, fixed(y, ignore_case = TRUE))) %>%
    filter(chrom.x == chrom.y) %>% 
    mutate(control_results = ifelse(control_type == "pos" & depth >= min_coverage, "Pass","fail")) %>% 
    mutate(control_results = ifelse(control_type == "neg" & depth <= min_coverage, "Pass",control_results)) %>% 
   # mutate(control_result = ifelse(depth >= min_coverage, "pass", "fail")) %>%
    glimpse() %>%
    select(Owner_Sample_ID, barcode, chrom = chrom.x, control_results,starts_with("MASIC")) %>%
    arrange(Owner_Sample_ID, chrom) %>%
    mutate(pass_targets = ifelse(control_results == "Pass",1,0)) %>%
    mutate(fail_targets = ifelse(control_results == "fail",1,0)) %>%
    group_by(Owner_Sample_ID,barcode) %>% 
    mutate(num_pass_targets = sum(pass_targets), num_fail_targets = sum(fail_targets)) %>% 
    select(-pass_targets,-fail_targets,-starts_with("MASIC")) %>%
    spread(chrom, control_results) 
  
  control_results_final<- coverage_matrix %>%
    gather("MAS_chr","mas_depth",starts_with("MASIC"),-Owner_Sample_ID, -barcode) %>% 
    # inner_join(control_defs, by = c("MAS_chr" = "chrom","Owner_Sample_ID" = "control_code" )) %>%
    fuzzyjoin::fuzzy_join(control_defs, mode = "inner", by = c("Owner_Sample_ID" = "control_code"), match_fun = function(x, y) str_detect(x, fixed(y, ignore_case = TRUE))) %>%
    filter(MAS_chr == chrom) %>% 
    mutate(mas_results = ifelse(control_type == "pos" & mas_depth >= min_coverage, "Pass","fail")) %>% 
    mutate(mas_results = ifelse(control_type == "neg" & mas_depth < min_coverage, "Pass",mas_results)) %>% 
    mutate(mas_results = ifelse(Owner_Sample_ID == "NTC" & mas_depth >= min_coverage, "Pass",mas_results)) %>%    #Adding a filter for NTC 
    mutate(mas_results = ifelse(Owner_Sample_ID == "NTC" & mas_depth < min_coverage, "fail",mas_results))  %>% 
    mutate(mas_fail_targets = ifelse(mas_results == "fail",1,0)) %>% 
    group_by(Owner_Sample_ID,barcode) %>% 
    mutate(num_fail_MASIC = sum(mas_fail_targets)) %>%
    select(Owner_Sample_ID,barcode,num_fail_MASIC,MAS_chr,mas_results) %>%
    spread(MAS_chr,mas_results) %>% 
    # select(Owner_Sample_ID,barcode,num_fail_MASIC) %>%
    full_join(control_results, by = c("Owner_Sample_ID","barcode")) 
    
    control_results_final = control_results_final %>%
      select(Owner_Sample_ID,barcode,num_pass_targets,num_fail_targets,num_fail_MASIC, con_sort_order, everything())
    
  manifest %>%
    inner_join(control_results_final) %>%
    filter(!(is.na(Owner_Sample_ID))) %>% 
    write_csv("control_results.csv")
  
### control samples can be identified using Owner_Sample_ID in control_results_final 
write_csv(tvr_table %>% filter(! Owner_Sample_ID %in% control_results_final$Owner_Sample_ID),"target_variants_external_results.csv")


  #Control frequency QC
  
  
  control_freq_defs %>%
    filter(!is.na(control_code)) %>%
    tidyr::gather(chrom, value, -control_code,-Control_range) %>%
    transform(control_code = as.character(control_code)) %>%
    spread(Control_range, value) -> control_freq_defs
  
  control_freq_qc = freq_matrix %>%
    tidyr::gather("chrom", "freq", -Owner_Sample_ID, -barcode) %>%
    inner_join(control_freq_defs, by = c("Owner_Sample_ID" = "control_code", "chrom")) %>%
    mutate(status = ifelse(freq <= max & freq >= min, "Pass","Fail")) %>% 
    mutate(status = ifelse(freq == "neg", "NA",status)) %>% 
    mutate(pass_freq = ifelse(status == "Pass",1,0)) %>% 
    mutate(fail_freq = ifelse(status == "Fail",1,0)) %>%
    mutate(NA_freq = ifelse(is.na(status),1,0)) %>%
    group_by(Owner_Sample_ID,barcode) %>% 
    mutate(num_pass_freq = sum(pass_freq,na.rm = T),num_fail_freq = sum(fail_freq, na.rm = T),num_NA_freq = sum(NA_freq, na.rm = T)) %>% 
    select(Owner_Sample_ID, barcode, chrom,status, num_pass_freq,num_fail_freq,num_NA_freq) %>%
    spread(chrom,status) 
  
  control_freq_qc = control_freq_qc[,str_sort(colnames(control_freq_qc),numeric = T)] %>%
    select(Owner_Sample_ID,barcode,num_pass_freq,num_fail_freq,num_NA_freq,everything())
  
  manifest %>%
    inner_join(control_freq_qc) %>%
    filter(!(is.na(Owner_Sample_ID))) %>%
    write.csv("Control_frequency_results.csv", row.names = F)
  
  #non-hotspot vars...
  non_hotspot_vars = variants %>%
    filter(!HS)
  
  manifest %>%
    inner_join(non_hotspot_vars) %>%
    filter(!(is.na(Owner_Sample_ID))) %>%
    write_csv("non_target_variants_results.csv")
  
  ### rename
  prefix <- manifest$Assay_Batch_Code[1]
  cmd <- sprintf("rename 's/^/%s_/' *_results.csv", prefix)
  system(cmd)
  return(tvr_table)
  
  
}
