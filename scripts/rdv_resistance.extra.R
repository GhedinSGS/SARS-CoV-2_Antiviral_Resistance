# Filtering vectors: 
ntlist = c('A','G','T','C','-')

aminoacids = c('G','A','L','M','F','W','K',
                'Q','E','S','P','V','I','C',
                'Y','H','R','N','D','T','*')


nsps = c('nsp1',
                  'nsp2',
                  'nsp3',
                  'nsp4',
                  'nsp5',
                  'nsp6',
                  'nsp7',
                  'nsp8',
                  'nsp9',
                  'nsp10',
                  'nsp11',
                  'nsp12a',
                  'nsp12b',
                  'nsp13',
                  'nsp14',
                  'nsp15',
                  'nsp16')

# Removing extra data: 
remove_patient = 13410 # Mirella had me remove - not a pro-longed infection

# a list of samples that were mislabeled when sequencing. The last 4 are from a patient (13410) that did not have a prolonged infection
removeit = c('22_CoV_888a_VTM2', '22_CoV_888_VTM2',
            "21_CoV_3403a_VTM2", "21-CoV-3403a-VTM2", "21_CoV_3420_VTM2", "21_CoV_3420a_VTM2") # last four ^ are from patient that did not have prolonged infect. 13410

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
col_list = c('black',
             '#bdcee0',
             '#dac8b9',
             '#84bf3b',
             '#66bace',
             '#4372d6',
             '#949cf3',
             '#6233e3',
             '#ba5fe5',
             '#b03766',
             '#db7032',
             '#eece50',
             '#cacaca',
             'black')

gene_ord_list = c('5\'UTR',
                  'ORF1a',
                  'ORF1b',
                  'nsp1',
                  'nsp2',
                  'nsp3',
                  'nsp4',
                  'nsp5',
                  'nsp6',
                  'nsp7',
                  'nsp8',
                  'nsp9',
                  
                  'nsp10',
                  'nsp11',
                  'nsp12a',
                  'nsp12b',
                  'nsp13',
                  'nsp14',
                  'nsp15',
                  'nsp16',

                  'S',
                  'ORF3a',
                  'E',
                  'M',
                  'ORF6',
                  'ORF7a',
                  'ORF7b',
                  'ORF8',
                  'N',
                  'ORF10',
                  '3\'UTR',
                  'INTERGENIC')


gene_ord_list2 = c('5\'UTR',
                  'nsp1',
                  'nsp2',
                  'nsp3',
                  'nsp4',
                  'nsp5',
                  'nsp6',
                  'nsp7',
                  'nsp8',
                  'nsp9',
                  
                  'nsp10',
                  'nsp11',
                  'nsp12',
                  'nsp12a',
                  'nsp12b',
                  'nsp13',
                  'nsp14',
                  'nsp15',
                  'nsp16',

                  'S',
                  'ORF3a',
                  'E',
                  'M',
                  'ORF6',
                  'ORF7a',
                  'ORF7b',
                  'ORF8',
                  'N',
                  'ORF10',
                  '3\'UTR',
                  'INTERGENIC')

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

# treatment linetype - note used
treat_list = c(2,3)
names(treat_list) = c('paxlovid', 'remdesivir')
treat_lines_scale = scale_linetype_manual(name = 'antiviral', values = treat_list)

# Prepping replicate information since naming isn't similar across seq runs
grab_replicates = function(df, percentage_type, coltype){
    pull_cols = c("Patient.number", "collection_date", "Order", "sample_id", c(coltype))
    
    return(df %>% filter(type == percentage_type) %>% select(all_of(pull_cols)) %>% rename("name" = coltype))
}

# columns to select for variant data: 
select_cols = c('name','segment','ntpos',
                'major','majorfreq','minor','minorfreq','binocheck',
                'totalcount','aapos','majoraa','majorcodon',
                'minoraa','minorcodon','refnt','gene_id',
                'refcodon','refaa','nonsyn')

# columns to merge by replicate data: 
mergeit= c('segment','ntpos','aapos','refnt','gene_id',
            'refcodon','refaa','Patient.number', 'collection_date', 'sample_id')


# to set up one 'variant' dataframe made up of both high and low-freq variants, pull these cols
pull_major_columns = c("Patient.number", "collection_date", "sample_id", 
                        "segment", "totalcount", 
                        "gene_id", "ntpos", "major", "majorfreq",  "refnt", 
                        "aapos", "majoraa", "majorcodon", "refaa", "refcodon",
                        "covflag", "filt_flag")


pull_minor_columns = c("Patient.number", "collection_date", "sample_id", 
                        "segment", "totalcount", 
                        "gene_id", "ntpos", "minor", "minorfreq", "refnt", 
                        "aapos", "minoraa", "minorcodon", "refaa", "refcodon",
                        "covflag", "filt_flag", "major", "majoraa")


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

PrepConDF = function(CONDF, NYGC_DF, ntlist = ntlist, aminoacids = aminoacids, coverage_cut = coverage_cut){
    CONDF = CONDF %>% filter(!sample_id %in% remove_2646) # remove samples from this patient where they did not recieve rdv treatment

    CONDF$totalcount.x[is.na(CONDF$totalcount.x)] = 0 

    CONDF$totalcount.y[is.na(CONDF$totalcount.y)] = 0

    # If the sample was seq'd at NYGC, we do not have rep data, fill in the '2nd' rep with the 1st rep info (reduces errors downstream)
    CONDF = CONDF %>% 
        mutate(major.y = ifelse(name.x %in% NYGC_DF$name, major.x, major.y), 
           majorfreq.y = ifelse(name.x %in% NYGC_DF$name, majorfreq.x, majorfreq.y), 
           totalcount.y = ifelse(name.x %in% NYGC_DF$name, totalcount.x, totalcount.y), 
           majoraa.y = ifelse(name.x %in% NYGC_DF$name, majoraa.x, majoraa.y), 
           majorcodon.y = ifelse(name.x %in% NYGC_DF$name, majorcodon.x, majorcodon.y)
          ) 

    # setting up the info for each variant by comparing the two replicate sets of data: 
    CONDF = CONDF %>% 
        rowwise() %>%
        mutate(major = ifelse(major.x == major.y & 
                              major.x %in% ntlist, major.x, 'N'), # if the majors are the same state that, if not then make 'n'

           majorcodon = ifelse(majorcodon.x == majorcodon.y & 
                               majoraa.x %in% aminoacids, majorcodon.x, NA),

           majoraa = ifelse(majoraa.x == majoraa.y & majoraa.x %in% aminoacids, majoraa.x, NA),
          
           majorfreq = ifelse(major.x == major.y & major.x %in% ntlist, 
                            (majorfreq.x + majorfreq.y)/2, NA),
           
           totalcount = (totalcount.x + totalcount.y)/2,
           
           covflag = ifelse(totalcount.x < coverage_cut | totalcount.y < coverage_cut, 'filter', 'pass')
                           ) %>% # generate a coverage flag
    ungroup()

    CONDF$major[is.na(CONDF$major)] = 'N' # change 'NA's to "N"

    CONDF = CONDF %>%
            rowwise() %>%
            mutate(filt_flag = ifelse(major %in% ntlist & covflag == 'pass', 'pass','filter')#,
                  ) %>% # generate filt flag using nt list (which includes dels) and totalcount check
        ungroup() 
    return(CONDF)
}

FilterCon = function(CONDF, METADF, major_cols, nsp_list ){
    CON_FILT = CONDF %>% 
                filter(filt_flag == "pass" & major != "N") %>% # the major nt needs to be agct
                select(all_of(major_cols)) %>% # pull out the columns we care about
                unique() %>% 
                rename("varnt"="major",
                       "varfreq" = "majorfreq",
                       "varaa" = "majoraa",
                       "varcodon" = "majorcodon")

    # keeping the 'varnt' (which in the confile is the major) as 'major' column in case we need to compare later
    CON_FILT$majnt = CON_FILT$varnt
    CON_FILT$majaa = CON_FILT$varaa 
    CON_FILT$majorfreq = CON_FILT$varfreq 

    # we are using nsp info rather than orf1a/b. But we need to merge with the reference seq data separating out that info to merge
    orf_info = CON_FILT %>% 
                filter(!gene_id %in% nsp_list) %>% 
                select(ntpos, aapos, gene_id) %>% 
                unique() %>%
                group_by(ntpos) %>%
                rename("merge_gene" = "gene_id", 
                       "merge_aapos" = "aapos") 

    CON_FILT = merge(CON_FILT, orf_info, by = c('ntpos'), all.x = TRUE) %>%
                filter(!gene_id %in% c('ORF1a', 'ORF1b')) %>%
                unique()


    CON_FILT = merge(CON_FILT, METADF, by = c('sample_id','collection_date','Patient.number'))
    
    return(CON_FILT)
}


AdjustVarDF = function(VARDF, NYGC_DF, ntlist = ntlist, aminoacids = aminoacids, varcov = 200, varfreq = 0.02){
    VARDF = VARDF %>%
            filter(!sample_id %in% remove_2646)

    VARDF$minorfreq.x = as.numeric(as.character(VARDF$minorfreq.x))
    VARDF$minorfreq.y = as.numeric(as.character(VARDF$minorfreq.y))
    VARDF$minorfreq.x[is.na(VARDF$minorfreq.x)] = 0 # if they don't have any cov - add 0
    VARDF$minorfreq.y[is.na(VARDF$minorfreq.y)] = 0

    VARDF = VARDF %>% 
        mutate(major.y = ifelse(name.x %in% NYGC_DF$name, major.x, major.y), 
               majorfreq.y = ifelse(name.x %in% NYGC_DF$name, majorfreq.x, majorfreq.y), 
               totalcount.y = ifelse(name.x %in% NYGC_DF$name, totalcount.x, totalcount.y), 
               majoraa.y = ifelse(name.x %in% NYGC_DF$name, majoraa.x, majoraa.y), 
               majorcodon.y = ifelse(name.x %in% NYGC_DF$name, majorcodon.x, majorcodon.y),
               minor.y = ifelse(name.x %in% NYGC_DF$name, minor.x, minor.y), 
               minorfreq.y = ifelse(name.x %in% NYGC_DF$name, minorfreq.x, minorfreq.y), 
               minoraa.y = ifelse(name.x %in% NYGC_DF$name, minoraa.x, minoraa.y), 
               minorcodon.y = ifelse(name.x %in% NYGC_DF$name, minorcodon.x, minorcodon.y)
              )

    VARDF = VARDF %>% 
        rowwise() %>%
        mutate(major = ifelse(major.x == major.y & major.x %in% ntlist, major.x, 'N'), # if the majors are the same state that, if not then make 'n'

               majorcodon = ifelse(majorcodon.x == majorcodon.y & 
                                   majoraa.x %in% aminoacids, majorcodon.x, NA),

               majoraa = ifelse(majoraa.x == majoraa.y & majoraa.x %in% aminoacids, majoraa.x, NA),

               majorfreq = ifelse(major.x == major.y & major.x %in% ntlist, 
                                (majorfreq.x + majorfreq.y)/2, NA),

               totalcount = (totalcount.x + totalcount.y)/2,
               minor = ifelse(minor.x == minor.y & minor.x %in% ntlist, minor.x, 'N'), # if the majors are the same state that, if not then make 'n'
               minorcodon = ifelse(minorcodon.x == minorcodon.y & 
                                   minoraa.x %in% aminoacids, minorcodon.x, NA),

               minoraa = ifelse(minoraa.x == minoraa.y & minoraa.x %in% aminoacids, minoraa.x, NA),

               minorfreq = ifelse(minor.x == minor.y & minor.x %in% ntlist & minorfreq.x >= varfreq & minorfreq.y >= varfreq, 
                                (minorfreq.x + minorfreq.y)/2, NA),

               covflag = ifelse(totalcount.x < varcov | totalcount.y < varcov, 'filter', 'pass'),

               filt_flag = ifelse(minor %in% ntlist &
                                              major %in% ntlist & 
                                              minorfreq.x >= varfreq & 
                                              minorfreq.y >= varfreq & 
                                              covflag == 'pass', 'pass','filter')
              ) %>% # generate a coverage flag
        ungroup()

    VARDF$minor[is.na(VARDF$minor)] = 'N' # change 'NA's to "N"
    
    return(VARDF)
}


FilterVar = function(VARDF, METADF, minor_cols, nsp_list){
    VAR_FILT = VARDF %>%
                    filter(filt_flag == "pass" & 
                           major != "N" & 
                           minor != "N" 
                          ) %>%
                select(c(all_of(minor_cols), 'majorfreq')) %>%
                unique()

    VAR_FILT = VAR_FILT %>% 
        filter(filt_flag == "pass") %>% 
        rename("varnt"="minor",
              "varfreq" = "minorfreq",
              "varaa" = "minoraa",
              "varcodon" = "minorcodon",
              "majnt" = 'major',
              "majaa" = "majoraa",
              "majorfreq" = "majorfreq") %>%  
        drop_na(varfreq)
    
    
    orf_info = VAR_FILT %>% 
                filter(!gene_id %in% nsp_list) %>% 
                select(ntpos, aapos, gene_id) %>% unique() %>%
                group_by(ntpos) %>%
                rename("merge_gene" = "gene_id", 
                       "merge_aapos" = "aapos")  

    VAR_FILT = merge(VAR_FILT, orf_info, by = c('ntpos'), all.x = TRUE) %>%
                filter(!gene_id %in% c('ORF1a', 'ORF1b')) %>%
                unique() 



    VAR_FILT = merge(VAR_FILT, METADF, by = c('sample_id','collection_date','Patient.number'))

    return(VAR_FILT)
}

CombineVar = function(CONDF, VARDF, ntdf, aadf, aminoacids){

    total_var = rbind(CONDF, VARDF) %>% filter(filt_flag == 'pass') %>% unique()


    total_var = merge(total_var, ntdf, by.x = c('ntpos','clade'), 
          by.y = c('ntpos','ref_clade'), all.x=TRUE)


    total_var = merge(total_var, aadf, 
                      by.x = c('merge_aapos','merge_gene','clade'), 
                      by.y=c('aapos','gene_id','ref_clade'), 
                      all.x = TRUE) %>%
        filter(varnt != nt & 
               varnt != '-' & 
               majnt != '-' ) %>%
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


plot_frequency = function(CON, VAR, pt, ptdf, pt_med, first_date, aa_colScale, g_list){

        temp = rbind(
                CON %>% 
                     filter(gene_id %in% g_list & 
                            Patient.number == pt & 
                            ntpos %in% ptdf$ntpos 
                           ) %>% 

                     select(collection_date,gene_id,majorfreq, aapos, major, majoraa, ntpos, filt_flag) %>% unique() %>%
                     rename("varfreq" = "majorfreq", 
                               "varnt" = "major", 
                               "varaa" = "majoraa"), 

                     VAR %>% 
                     filter(gene_id %in% g_list & 
                            Patient.number == pt & 
                            ntpos %in% ptdf$ntpos 
                                        ) %>% 
                     select(collection_date, gene_id,
                            minorfreq, aapos, minor,minoraa, ntpos, filt_flag) %>% unique() %>%
                     rename("varfreq" = "minorfreq", 
                               "varnt" = "minor", 
                               "varaa" = "minoraa"), 

                    VAR %>% 
                     filter(gene_id %in% g_list & 
                            Patient.number == pt & 
                            ntpos %in% ptdf$ntpos 
                                        ) %>% 
                    mutate(filt_flag = ifelse(major %in% ntlist & totalcount > coverage_cut, 'pass','filter')) %>%
                     select(collection_date, gene_id,
                            majorfreq, aapos, major,majoraa, ntpos, filt_flag) %>% unique() %>%
                         rename("varfreq" = "majorfreq", 
                               "varnt" = "major", 
                               "varaa" = "majoraa")
                    ) %>%
        unique()


    temp$gene_id = factor(temp$gene_id, levels = g_list)
    
    
    temp = temp %>%
                mutate(aapos = as.numeric(as.character(aapos)), 
                       aapos = ifelse(gene_id == 'nsp12b', aapos + 9, aapos)) %>%
            unique() %>%
            drop_na(varaa, varfreq)

    print(levels(factor(temp$gene_id)))

        
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