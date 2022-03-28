#+
vcf_to_dataframe <- function(vcf_files){
    require(tidyverse)
    require(fs)
    require(vcfR)

    temp = read.vcfR(vcf_files$vcf_out %>% unique())

    temp = vcfR2tidy(temp)

    temp = temp$fix %>%
    as_tibble() %>%
    mutate(filename = vcf_files$vcf_out %>% unique()) %>%
    select(filename, everything())

    variant_table_snv = temp %>%
      filter(!(str_detect(ALT, ","))) %>%
      glimpse()

    #This step is added to prevent selecting columns with NA which will cause problems in seperate rows step
    
    temp %>%
      filter(str_detect(ALT, ",")) -> get_names
    
    get_names[,complete.cases(t(get_names))] %>% colnames() -> id_all
    id_select<-c("ALT", "AO", "SAF", "SAR", "FAO", "AF", "FSAF", "FSAR", "TYPE", "LEN", "HRUN", "MLLD", 
                 "FWDB", "REVB", "REFB", "VARB", "STB", "STBP", "RBI", 
                 "FR", "SSSB", "SSEN", "SSEP", "PB", "PBP", "FDVR") 
    id_intersect<-intersect(id_all,id_select)
    
    
    variant_table_mult_split = temp %>%
      filter(str_detect(ALT, ",")) %>%
      separate_rows(id_intersect, sep = ",") %>%
      glimpse()

    variant_table_return = bind_rows(variant_table_snv, variant_table_mult_split)

    return(variant_table_return)


}
