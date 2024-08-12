# list of nts output by timo
ntlist = c('A','G','T','C','-')

# list of amino acids
aminoacids = c('G','A','L','M','F','W','K',
                'Q','E','S','P','V','I','C',
                'Y','H','R','N','D','T','*')

nsps = c('nsp1','nsp2',
        'nsp3','nsp4','nsp5','nsp6',
        'nsp7', 'nsp8','nsp9','nsp10',
        'nsp11','nsp12a','nsp12b',
        'nsp13', 'nsp14', 'nsp15','nsp16')

# Removing extra data: 
remove_patient = 13410 # Mirella had me remove - not a pro-longed infection

# a list of samples that were mislabeled when sequencing. The last 4 are from a patient (13410) that did not have a prolonged infection
removeit = c('22_CoV_888a_VTM2', '22_CoV_888_VTM2',
            "21_CoV_3403a_VTM2", "21-CoV-3403a-VTM2", 
            "21_CoV_3420_VTM2", "21_CoV_3420a_VTM2") # last four ^ are from patient that did not have prolonged infect. 13410

# These samples were collected from 2646 during their first infection, before remdesivir treatment. removing
remove_2646 = c('2646_2020-12-31', '2646_2021-01-04','2646_2021-01-07','2646_2021-01-11')

# Setting up plotting features: 
PlotTheme1 = theme_bw() +
            theme(axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              strip.background = element_rect(colour="black", fill="white"),
              text = element_text(size = 12))

# coding region colors: 
col_list = c('black','#bdcee0','#dac8b9',
             '#84bf3b','#66bace',
             '#4372d6','#949cf3','#6233e3',
             '#ba5fe5', '#b03766', '#db7032',
             '#eece50', '#cacaca','black')

gene_ord_list = c('5\'UTR','ORF1a', 'ORF1b',
                  'nsp1','nsp2','nsp3', 'nsp4', 'nsp5',
                  'nsp6','nsp7','nsp8', 'nsp9', 'nsp10',
                  'nsp11', 'nsp12a','nsp12b', 'nsp13','nsp14',
                  'nsp15', 'nsp16','S','ORF3a', 'E',
                  'M','ORF6', 'ORF7a','ORF7b','ORF8',
                  'N','ORF10', '3\'UTR','INTERGENIC')


gene_ord_list2 = c('5\'UTR',
                  'nsp1','nsp2', 'nsp3','nsp4','nsp5',
                  'nsp6', 'nsp7','nsp8','nsp9', 'nsp10',
                  'nsp11', 'nsp12', 'nsp12a','nsp12b','nsp13',
                  'nsp14', 'nsp15','nsp16', 'S','ORF3a',
                  'E', 'M', 'ORF6', 'ORF7a','ORF7b',
                  'ORF8','N','ORF10','3\'UTR', 'INTERGENIC')

aacolors = as.vector(cols25(21))
names(aacolors) = aminoacids
aa_colScale_fill <- scale_fill_manual(name = "aa",values = aacolors)
aa_colScale <- scale_colour_manual(name = "aa",values = aacolors)

short_gene_list = c('5\'UTR','ORF1a','ORF1b','S','ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10','3\'UTR')
names(col_list) = short_gene_list
gene_colScale_fill = scale_fill_manual(name = "gene id",values = col_list)
gene_colScale = scale_colour_manual(name = "gene id",values = col_list)

# filtering point shape: 
shape_list = c(19, 17)
names(shape_list) = c('pass', 'filter')
filter_shapes = scale_shape_manual(name = 'flag', values = shape_list)

# Prepping replicate information since naming isn't similar across seq runs
grab_replicates = function(df, percentage_type, coltype){
    pull_cols = c("Patient.number", "collection_date", "Order", "sample_id", c(coltype))
    return(df %>% filter(type == percentage_type) %>% select(all_of(pull_cols)) %>% rename("name" = coltype))
}

pull_treatments = c('dexamethasone','evusheld','methylprednisolone','paxlovid','prednisone','regen_cov',
                   'remdesivir','sotrovimab','steroid',
                   'rituximab','imbrutinib','tocilizumab',
                   'antithymocyte_globulin','bebtelovimab')


# columns to select for variant data: 
select_cols = c('name','segment','ntpos',
                'major','majorfreq','minor','minorfreq','binocheck',
                'totalcount','aapos','majoraa','majorcodon',
                'minoraa','minorcodon','refnt','gene_id',
                'refcodon','refaa','nonsyn')




# to set up one 'variant' dataframe made up of both high and low-freq variants, pull these cols
pull_major_columns2 = c("Patient.number", "collection_date", "sample_id", 
                        "segment", "totalcount", 
                        "gene_id", "ntpos", "major", "majorfreq",  "refnt", 
                        "aapos", "majoraa", "majorcodon", "refaa")     



pull_minor_columns2 = c("Patient.number", "collection_date", "sample_id", 
                        "segment", "totalcount", 
                        "gene_id", "ntpos", "minor", "minorfreq", "refnt", 
                        "aapos", "minoraa", "minorcodon", "refaa",'major')


# Setting up the correct reference information: 
cat_dfs = function(files, type){
    return_df = data.frame()
    
    for (f in files){   
            return_df = rbind(return_df, read.csv(f)) 
        }
    
    if (type == "aa"){
       return_df = return_df %>%
                    mutate(gene_id = gsub("-gene consensus", "", gene_id),
                           gene_id = gsub("variant ", "", gene_id))
        
    }
    return(return_df)
}


plot_isnv = function(pt, pt_med, ptdf, treat_colScale_fill, first_date){
    p6 =  
    ggplot() +
     geom_rect(data = pt_med, 
                  aes(xmin=start_first, 
                      xmax = end_first,
                     ymin = 0, 
                     ymax = max(ptdf$n + 2), 
                     fill = treatment_type), alpha = 0.2) +

         geom_line(data = ptdf, aes(x = to_diagnosis, y=n, group = Patient.number), linewidth = 1) + 
        geom_point(data = ptdf, size = 3, aes(x = to_diagnosis, y=n)) +
        labs(y="number of iSNV (2-98%)", x= "days to first diagnosis") + 
    PlotTheme1 + 
    treat_colScale_fill + 
    expand_limits(y = 0) +
    facet_grid(Patient.number + clade~.)

return(p6)
}


PrepDF = function(DF, ntlist = ntlist, aminoacids = aminoacids, coverage_cut = coverage_cut, min_varfreq, min_cov = varcov){
    DF = DF %>% filter(!sample_id %in% remove_2646) # remove samples from this patient where they did not recieve rdv treatment

    DF$totalcount.x[is.na(DF$totalcount.x)] = 0 

    DF$totalcount.y[is.na(DF$totalcount.y)] = 0
    
    DF$varfreq.x[is.na(DF$varfreq.x)] = 0 

    DF$varfreq.y[is.na(DF$varfreq.y)] = 0

    
    # setting up the info for each variant by comparing the two replicate sets of data: 
    DF = DF %>% 
        filter(varfreq.x > min_varfreq & varfreq.y > min_varfreq) %>% 
        unique() %>%
        rowwise() %>%
        mutate(varnt = ifelse(varnt %in% ntlist, varnt, 'N'), # if the majors are the same state that, if not then make 'n'

           varcodon = ifelse(varaa %in% aminoacids, varcodon, NA),

           varaa = ifelse(varaa %in% aminoacids, varaa, NA),
          
           varfreq = (varfreq.x + varfreq.y)/2,
           
           totalcount = (totalcount.x + totalcount.y)/2,
           
           vartype = ifelse(varfreq >= 0.50, 'major','minor'), 
           
           covflag = ifelse(vartype == 'major' & totalcount.x < coverage_cut |
                            vartype == 'major' & totalcount.y < coverage_cut | 
                            vartype == 'major' & totalcount < coverage_cut |
                            
                            vartype == 'minor' & totalcount.x < min_cov |
                            vartype == 'minor' & totalcount.y < min_cov | 
                            vartype == 'minor' & totalcount < min_cov, 'filter', 'pass'),
                       
           freqflag = ifelse(varfreq.x > min_varfreq & varfreq.y > min_varfreq, 'pass', 'filter')) %>% # generate a coverage flag
    ungroup()

    DF$varnt[is.na(DF$varnt)] = 'N' # change 'NA's to "N"

    DF = DF %>%
            rowwise() %>%
            mutate(filt_flag = ifelse(varnt %in% ntlist & 
                                      covflag == 'pass' &
                                      freqflag  == 'pass', 'pass','filter')#,
                  ) %>% # generate filt flag using nt list (which includes dels) and totalcount check
        ungroup() 
    
    return(DF)
}

FilterDF = function(DF, METADF, nsp_list ){
    DF_FILT = DF %>% 
                filter(filt_flag == "pass" & varnt != "N") %>% # the major nt needs to be agct
                unique() 
    
    # we are using nsp info rather than orf1a/b. But we need to merge with the reference seq data separating out that info to merge
    orf_info = DF_FILT %>% 
                filter(!gene_id %in% nsp_list) %>% 
                select(ntpos, aapos, gene_id) %>% 
                unique() %>%
                group_by(ntpos) %>%
                rename("merge_gene" = "gene_id", 
                       "merge_aapos" = "aapos") 

    DF_FILT = merge(DF_FILT, orf_info, by = c('ntpos'), all.x = TRUE) %>%
                filter(!gene_id %in% c('ORF1a', 'ORF1b')) %>%
                unique()

    DF_FILT = merge(DF_FILT, merge_meta, by = c('sample_id', 'Patient.number', 'collection_date')) 
    
    return(DF_FILT)
}


CombineVar = function(DF, ntdf, aadf, aminoacids){

    total_var = DF %>% filter(filt_flag == 'pass') %>% unique()


    total_var = merge(total_var, ntdf, by.x = c('ntpos','clade'), 
          by.y = c('ntpos','ref_clade'), all.x=TRUE)


    total_var = merge(total_var, aadf, 
                      by.x = c('merge_aapos','merge_gene','clade'), 
                      by.y=c('aapos','gene_id','ref_clade'), 
                      all.x = TRUE) %>%
        filter(varnt != nt & 
               varnt != '-' & 
               majnt.x != '-' &
              majnt.y != '-') %>% # can't merge bc of the 50/50 situation 
        mutate(aapos = as.numeric(as.character(aapos)), 
              varinfo = glue("{gene_id}: {aapos}")) %>%
        filter(varinfo != 'nsp8: 125') %>% 
        select(-varinfo) %>%
        ungroup() %>% 
        rowwise() %>%
        mutate(nonsyn = ifelse(varaa %in% aminoacids & varaa != aa & varaa == "*", 'stop',
                                ifelse(varaa %in% aminoacids & varaa != aa & varaa != "*", 'nonsyn',
                                ifelse(varaa %in% aminoacids & varaa == aa, 'syn', 'non-coding')))) %>%

        unique() 

    return(total_var)

}


plot_frequency = function(DF, pt, ptdf, pt_med, first_date, aa_colScale, g_list){

    temp = DF %>% 
                filter(gene_id %in% g_list & 
                            Patient.number == pt & 
                            ntpos %in% ptdf$ntpos 
                           ) %>% 
                select(Patient.number, collection_date,gene_id,varfreq, aapos, varnt, varaa, ntpos, filt_flag) %>% 
                unique() 


    temp$gene_id = factor(temp$gene_id, levels = g_list)
    
    
    temp = temp %>%
                mutate(aapos = as.numeric(as.character(aapos)), 
                       aapos = ifelse(gene_id == 'nsp12b', aapos + 9, aapos),
                      to_diagnosis = as.numeric(collection_date - first_date)) %>%
            unique() %>%
            arrange(Patient.number, ntpos, varaa, to_diagnosis) %>%
            drop_na(varaa, varfreq) 

        
        p8 =  
            ggplot() + 
                geom_vline(data = pt_med, aes(xintercept = end_first), color = 'gray', linetype = 2) + 
                geom_vline(data = pt_med, aes(xintercept = start_first), color = 'gray', linetype = 2) + 
                geom_line(data = temp, aes(x=collection_date-first_date, y = varfreq, color = glue("{varaa}"), group = varnt),
                          size = 1) + 
                geom_point(data = temp, aes(shape= filt_flag, x=collection_date-first_date, y = varfreq, color = glue("{varaa}"), group = varnt),
                           size = 3) + 
            PlotTheme1 +
            ylim(0, 1) +
            aa_colScale + 
            filter_shapes +
            scale_x_continuous(breaks = seq(0, max(ptdf$to_diagnosis), 10)) +    
            facet_grid(gene_id + aapos~Patient.number + clade) + 
            labs(y="relative frequency of SNV", x="days to first diagnosis")

        return(p8)
}

GeneratePatientDF = function(df, pt, minfreq, maxfreq){
    tempdf = df %>% 
        filter(Patient.number == pt & varfreq >= minfreq & varfreq <= maxfreq) %>% # filter for patient info
            group_by(sample_id, Patient.number, to_diagnosis) %>% # group by patient info and collection date
            add_tally() %>% # add tally to determine number of var / sample or date
            ungroup() %>% # ungroup
        group_by(ntpos, varnt, gene_id) %>% #group by nucleotide position
        add_tally(name = 'number_of_samples') %>% # determine number of patient samples with variant
        ungroup() %>% # ungroup
        mutate(freq_cutoff = ifelse(varfreq > 0.98, '>0.98', '<=0.98')) %>% # add a flag that indicates if 
        group_by(ntpos, varnt, gene_id, freq_cutoff) %>% # group by the flag
        add_tally(name = 'cutoff_count') %>% # count the number that are fixed or not
        ungroup() %>% 
        # add a filter flag that would remove variants that are greater than 98% (ie fixed and maybe not intrahost)
        mutate(filter_it = ifelse(freq_cutoff == '>0.98' & number_of_samples > 2 &  cutoff_count == number_of_samples, 'yes','no'))

    return(tempdf)  # return df
}

source_data = function(DF, pt, ptdf, pt_med, first_date, aa_colScale, g_list){

    temp = DF %>% 
                filter(gene_id %in% g_list & 
                            Patient.number == pt & 
                            ntpos %in% ptdf$ntpos 
                           ) %>% 
                select(Patient.number, collection_date,gene_id,varfreq, aapos, varnt, varaa, ntpos, filt_flag) %>% 
                unique() 


    temp$gene_id = factor(temp$gene_id, levels = g_list)
    
    
    temp = temp %>%
                mutate(aapos = as.numeric(as.character(aapos)), 
                       aapos = ifelse(gene_id == 'nsp12b', aapos + 9, aapos),
                      #gene_id = ifelse(gene_id == 'nsp12b', 'nsp12', gene_id),
                      to_diagnosis = as.numeric(collection_date - first_date)) %>%
            unique() %>%
            arrange(Patient.number, ntpos, varaa, to_diagnosis) %>%
            drop_na(varaa, varfreq) 

    print(levels(factor(temp$gene_id)))
    
    return(temp)
}