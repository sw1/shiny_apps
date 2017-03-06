library(plotly)
library(shiny)
library(tidyr)
library(dplyr)
library(stringr)
library(lazyeval)
library(viridis)
library(rstanarm)


source('/data/sw1/Dropbox/stm_microbiome/code_active/stm_functions.R',local=TRUE)
source('/data/sw1/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R',local=TRUE)


pretty_taxa_names <- function(x){
  taxonomy <- c('Species','Genus','Family','Order','Class','Phylum','Kingdom')
  
  taxon <- gsub('[a-z]__','',x[,taxonomy[1]])
  taxon <- ifelse(taxon != '',paste(gsub('[a-z]__','',x[,taxonomy[2]]),taxon),'')
  bad_names <- taxon == ''
  
  for (tax in taxonomy[-1]){
    clean <- gsub('[a-z]__','',x[bad_names,tax])
    taxon[bad_names] <- ifelse(clean != '',paste(clean,stringr::str_to_lower(tax)),'')
    bad_names <- taxon == ''
  }
  
  taxon
}


shinyServer(function(input, output, session) {
  
  fitn <- 2
  K <- 25
  
  DAT <- readRDS(sprintf('/data/sw1/Dropbox/stm_microbiome/qiime_active/gevers/models/%s/stm_k_%s.rds',K,K))
  fit <- DAT$fits[[fitn]]
  META <- DAT$meta
  rownames(META) <- META$SampleID
  
  set.seed(453)
  eff <- stm::estimateEffect(as.formula(sprintf('1:%s ~ DIAGNOSIS',K)), 
                             fit, META, uncertainty='None')
  eff_plot <- stm::plot.estimateEffect(eff, 'DIAGNOSIS', model=fit, 
                                       topics=1:K, method='difference',cov.value1=1,cov.value2=0)
  
  topic_ord <- data.frame(topic=paste0('T',1:K)[order(unlist(eff_plot$means),decreasing=TRUE)],
                          topic_estimate=unlist(eff_plot$means)[order(unlist(eff_plot$means),decreasing=TRUE)],
                          topic_rank=1:K)
  topic_ord$topic <- factor(topic_ord$topic,levels=rev(topic_ord$topic),ordered=TRUE)
  
  
  p_est <- data.frame(topic=topic_ord$topic,est=topic_ord$topic_estimate) %>%
    mutate(sig=ifelse(est>0,'1','0')) %>%
    ggplot(aes(topic,est,color=sig)) +
    geom_hline(yintercept=0,linetype=3) +
    geom_point(size=4) +
    theme_void() + 
    labs(x='',y='Estimate') +
    theme(legend.position='none') +
    ggrepel::geom_label_repel(data=. %>% 
                                filter(est >= max(est) | est <= (min(est))) %>% 
                                mutate(group=ifelse(est < 0, 'CD-', 'CD+')),
                              aes(topic,est,fill=sig,label=group),
                              color='white',fontface='bold',size=10) +
    scale_color_brewer(type='qual',palette=6,direction=-1) +
    scale_fill_brewer(type='qual',palette=6,direction=-1) 
  
  
  BETA <- t(exp(fit$beta$logbeta[[1]]))
  rownames(BETA) <- fit$vocab
  colnames(BETA) <- paste0('T',1:K)
  
  BETA_RANK <- apply(dplyr::desc(BETA),2,dense_rank)
  BETA[BETA < 1e-5] <- 1e-5
  BETA <- log(BETA)
  BETA <- BETA[which(apply(BETA_RANK,1,min) <= 25),]
  
  dd_row <- as.dendrogram(hclust(vegan::vegdist(BETA,method='euclidean'),method='ward.D2'))
  otu_ord <- data.frame(otu=rownames(BETA)[order.dendrogram(dd_row)],
                        otu_rank=1:nrow(BETA))
  
  suppressWarnings(
    DF <- as.data.frame(BETA) %>% 
      mutate(otu=rownames(BETA)) %>%
      gather(topic,probability,-otu) %>%
      left_join(otu_ord,by='otu') %>%
      left_join(topic_ord,by='topic') %>%
      mutate(long=DAT$counts$ids[as.character(otu),'long']) %>%
      left_join(data.frame(taxon=pretty_taxa_names(DAT$taxa),
                           long=rownames(DAT$taxa)),by='long') %>%
      mutate(taxon=paste(taxon,gsub('^otu(.*)$','\\1',otu))) %>%
      arrange(dplyr::desc(topic_rank),otu_rank) %>%
      mutate(otu=factor(otu,levels=unique(otu),ordered=TRUE),
             taxon=factor(taxon,levels=unique(taxon),ordered=TRUE),
             topic=factor(topic,levels=unique(topic),ordered=TRUE))
  )
  
  p_tax <- DF %>%
    ggplot(aes(topic,taxon)) + 
    geom_raster(aes(fill=probability)) +
    scale_fill_viridis(name='Probability Rank') +
    labs(x='',y='',fill='Probability Rank') +
    theme(legend.position='none',
          axis.title.y=element_text(face='bold',size=20),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_text(angle=-90,hjust=0,vjust=.5))
  
  
  stan <- readRDS(sprintf('/data/sw1/Dropbox/stm_microbiome/qiime_active/gevers/models/%s/stan_lev3_k_%s_fit%s.rds',K,K,fitn-1))
  stan1 <- stan$stan
  biom_file <- biom::read_biom(sprintf('/data/sw1/Dropbox/stm_microbiome/qiime_active/gevers/models/%s/ko_k_%s_fit%s.biom',K,K,fitn-1))
  ko <- as.matrix(biom::biom_data(biom_file))
  ko <- ko[rowSums(ko)>0,]
  kegg_meta_2 <- kegg_metadata(biom_file,2)[rownames(ko)]
  kegg_meta <- sapply(kegg_meta_2,function(x) x$description)
  kegg_meta <- data.frame(ko=names(kegg_meta),description=kegg_meta)
  
  
  stan_sum <- stan1$stan_summary
  stan_sig_mean <- stan_sum[grepl('^b.* topic\\:pw',rownames(stan_sum)),c('mean')]
  stan_sig <- stan_sum[grepl('^b.* topic\\:pw',rownames(stan_sum)),c('10%','90%')]
  target_pws <- unique(gsub('_',' ',gsub('^b.*topic\\:pw\\:([0-9]+)\\:(.*)\\]$','\\2',rownames(stan_sig))))
  target_pws <- target_pws[!grepl('^b\\[\\(',target_pws)]
  stan_sig <- rownames(stan_sig)[rowSums(sign(stan_sig)) != 0]
  stan_sig_mean <- stan_sig_mean[stan_sig]
  stan_sig <- gsub('^b.*topic\\:pw\\:([0-9]+)\\:(.*)\\]$','\\1 \\2',stan_sig)
  stan_sig <- gsub('_',' ',stan_sig)
  stan_sig <- data.frame(do.call('rbind',str_split(stan_sig,' ',2)),
                         weight=stan_sig_mean)
  colnames(stan_sig) <- c('topic','pw','weight')
  stan_sig$topic <- paste0('T',stan_sig$topic)
  
  int_coef <- ranef(stan1)$`topic:pw`
  int_coef_names <- str_split(rownames(int_coef),':')
  int_mat <- matrix(0,length(target_pws),K,dimnames=list(target_pws,1:K))
  for (i in seq_len(nrow(int_coef))){
    k <- int_coef_names[[i]][1]
    pw <- int_coef_names[[i]][2]
    int_mat[pw,k] <- int_coef[i,1]
  }
  
  pw_row <- as.dendrogram(hclust(vegan::vegdist(int_mat,method='euclidean'),method='ward.D2'))
  pw_ord <- data.frame(pw=rownames(int_mat)[order.dendrogram(pw_row)],
                       pw_rank=1:nrow(int_mat))
  
  suppressWarnings(
    df <- data.frame(int_mat,pw=rownames(int_mat)) %>%
      gather(topic,weight,-pw) %>%
      mutate(topic=gsub('X','T',topic)) %>%
      left_join(topic_ord,by='topic') %>%
      left_join(pw_ord,by='pw') %>%
      arrange(topic_estimate,pw_rank) %>%
      mutate(topic=factor(topic,levels=unique(topic)),
             pw=factor(pw,levels=unique(pw)))
  )
  
  suppressWarnings(
    stan_sig2 <- stan_sig %>%
      left_join(df,by=c('topic','pw')) %>%
      rename(weight=weight.x)
  )
  
  df <- df %>% 
    filter(pw %in% stan_sig2$pw) %>%
    mutate(pw=factor(pw,levels=unique(pw)))
  
  p_fxn <- df %>%
    ggplot(aes(topic,pw,fill=weight)) +
    geom_raster() +
    geom_point(data=stan_sig2 %>% filter(weight.y>0),
               aes(topic,pw),color='red',shape=3) +
    geom_point(data=stan_sig2 %>% filter(weight.y<0),
               aes(topic,pw),color='cyan',shape=3) +
    scale_fill_viridis(name='Coefficient') +
    labs(x='',y='',fill='Coefficient') +
    theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=.5),
          legend.position='left') +
    geom_vline(xintercept=3.5,color='white',size=2) +
    geom_vline(xintercept=22.5,color='white',size=2) +
    theme(legend.position='none') 
  
  
  
  output$est <- renderPlot({
    p_est
  })
  
  output$tax <- renderPlotly({
    ggplotly(p_tax)
  })
  
  output$fxn <- renderPlotly({
    ggplotly(p_fxn,source='hm_fxn')
  })
  
  output$tbl <- renderTable({
    
    s <- event_data('plotly_click',source='hm_fxn')
    
    if (length(s)){
      k <- levels(df$topic)[s[['x']]]
      pw <- levels(df$pw)[s[['y']]]
      
      pw_desc <- pw
      pw_topics <- k
      
      suppressWarnings(
        stan1$data %>%
          mutate(topic=paste0('T',topic)) %>%
          filter(pw %in% pw_desc,
                 topic %in% pw_topics,
                 count > 0) %>%
          left_join(kegg_meta,by='ko') %>%
          dplyr::select(topic,count,description) %>%
          filter(!is.na(description)) %>%
          group_by(topic,description) %>%
          summarize(count=sum(count)) %>%
          spread(topic,count,fill=0) %>%
          arrange_(.dots=interp(~desc(x),x=as.name(pw_topics))) %>%
          rename(Description=description)
      )
    }else{
      NULL
    }
    
  },digits=0,width='900',align='l',rownames=FALSE,striped=TRUE,hover=TRUE,bordered=FALSE,spacing='s')
  
})
