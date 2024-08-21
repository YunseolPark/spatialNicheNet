plot_prioritization = function(original, updated, add_col=NULL){
  library(reshape2)
  library(ggplot2)
  library(ggnewscale)
  
  ## Steps:
  # given two NicheNet (NN) outputs (function: "prioritization_table")
  # in this case: prior_table , prior_table_update
  prior_table = original
  prior_table_update = updated
  
  
  # 1. generate two df with two cols: ligand, rank
  NN_ligand_auc_ranked_df = as.data.frame(prior_table %>%
                                          mutate(rank=c(1:nrow(.)))) %>%
    select(sender, ligand, receptor, rank)
  ENN_ligand_auc_ranked_df = as.data.frame(prior_table_update %>%
                                           mutate(rank=c(1:nrow(.)))) %>%
    select(sender, ligand, receptor, rank)
  
  # 2. merge over ligands and arrange using NN as "ranking ref" 
  merged_rank=merge(NN_ligand_auc_ranked_df, 
                    ENN_ligand_auc_ranked_df,
                    by=c("ligand", "receptor", "sender"),all=T) %>%
    arrange(rank.x)
  
  # 3. Create a dataframe with all necessary variables for plotting in ggplo2
  ranking_NN= merged_rank$rank.x
  ranking_ENN=merged_rank$rank.y
  receptor=merged_rank$receptor
  names=merged_rank$ligand
  sender = merged_rank$sender
  dat = data.frame(team=names, sender=sender, receptor=receptor,
                    rankA=ranking_NN,rankB=ranking_ENN)
  dat$color = ifelse(is.na(dat$rankB), FALSE, TRUE)
  dat_top = dat[(dat$rankA %in% c(1:30) | dat$rankB %in% c(1:30)),]
  
  #turn to df into a melted version
  dat_m = melt(dat_top,id.var=c("team", "receptor", "sender", 'color')) %>%
    arrange(desc(value))
  dat_m = dat_m[(dat_m$variable=='rankA')|
                  (dat_m$variable=='rankB'&!is.na(dat_m$value)) &
                dat_m$value<=30,]
  
  dat_m$sender %>% unique
  #colourCount = length(unique(dat_m$sender))
  #getPalette = colorRampPalette(brewer.pal(9, 'Dark2'))
  #library(randomcoloR)
  #palette = randomColor(15, luminosity='bright')
  #palette = c("#1eff3c", "#45d379", "#4aa7ba", "#53ba23", "#3d71d3", "#0ea010", "#10baac", "#2b2ed1", "#714bb7",
  #  "#a515f2", "#01b70a", "#f95f40", "#473ec9", "#6f00ef", "#00f2e6")
  #palette = distinctColorPalette(15, runTsne = TRUE)
  # For delauany
  #palette = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')
  # For squidpy
  #palette = c('#1b9e77','#984ea3','#e41a1c','#377eb8','#d95f02','#e6ab02','#7570b3','#e7298a','#66a61e')
  #palette = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','black','#a65628','#f781bf','#7570b3','#1b9e77')
  dat_m %>% mutate(sender=case_when(sender=='Central Vein Endothelial cells'~'CVECs',sender=='Portain Vein Endothelial cells'~'PVECs'))
  if (is.null(add_col)){
    palette = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c',  '#008080', '#9a6324', '#800000', '#808000', '#000075', '#808080', '#000000')
  }else{
    palette = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c',  '#008080', '#9a6324', '#800000', '#808000', '#000075', '#808080', '#000000','#e6ab02')
  }
    fig = ggplot(dat_m, aes(x=variable, y=value, 
                    group=interaction(team, receptor, sender))) +
        geom_line(arrow = arrow(length = unit(0.15, "cm"))) +
        geom_text(data=dat_m[dat_m$variable=="rankA",],
                  aes(label=paste(team,'-',receptor), color=sender),
                  hjust=  1.1, ,key_glyph='point')+
        scale_color_manual(values=palette)+
        geom_text(data=dat_m[dat_m$variable=="rankB"&dat_m$value<=30&
                               !is.na(dat_m$value),],
                  aes(label=paste(team,'-',receptor), color=sender),
                  hjust= -0.1,show.legend = F)+
        #scale_color_brewer(palette='Dark2')+
        #scale_color_manual(values=getPalette(colourCount))+
        geom_vline(xintercept = c(1,2)) +
        theme_classic() +
        theme(
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())+ylim(30,1)
  return(fig)
}
