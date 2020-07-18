source("functions.R")
require(ggtree)
dirw = file.path(dird, '01_tfbs')
diri = file.path(dirw, 'raw')
mdic = v3_to_v4(opt='dict_all')
ndic = read_symbol(opt = 'dict')

#{{{ compile motifs
ctag = 'cisbp'
#{{{ process cis-BP motifs
org = 'Zmays'
org = 'Athaliana'
diri1 = file.path(diri, 'cisbp', org, 'TF_Information.txt')
ti = read_tsv(diri1) %>%
    select(gid=DBID, name=TF_Name, motif=Motif_ID, status=TF_Status,
           method=Motif_Type, src_ref=MSource_Identifier, src_type=MSource_Type,
           src_id = DBID_1) %>%
    mutate(src_type = ifelse(src_type==method, src_type, str_c(src_type,method,sep="_"))) %>%
    select(-method) %>%
    filter(motif != '.') %>%
    mutate(fm = sprintf("%s/cisbp/%s/pwms_all_motifs/%s.txt", diri, org, motif)) %>%
    mutate(size=map_dbl(fm, file.size)) %>%
    arrange(size) %>% filter(size > 12) %>% select(-size) %>%
    #mutate(pwm = map(fm, read_cisbp)) %>%
    select(-fm) %>%
    mutate(org = !!org) %>% select(org, everything())
ti %>% print(width=Inf)
ti %>% count(src_type) %>% print(n=50)
ti %>% count(gid, src_id) %>% filter(n>1)
ti %>% count(status) %>% print(n=50)

#{{{ maize gene ID conversion
if (org == 'Athaliana') {
    tia = ti
} else if (org == 'Zmays') {
    tim = ti %>%
        mutate(gid = ifelse(gid %in% names(mdic), mdic[gid], gid)) %>%
        mutate(name = ifelse(name %in% names(mdic), mdic[name], name)) %>%
        mutate(name = ifelse(name %in% names(ndic), ndic[name], name))
}
#}}}

toc = rbind(tia, tim) %>%
    mutate(src_type = ifelse(src_type=='Dap-seq', 'DAP', src_type)) %>%
    mutate(src_type = ifelse(src_type %in% c('JASPAR_ChIP-chip','JASPAR_ChIP-seq'),'ChIP', src_type)) %>%
    mutate(src_type = ifelse(src_type=='JASPAR_PBM, CSA and or DIP-chip', 'JASPAR', src_type)) %>%
    mutate(src_type = ifelse(src_type=='JASPAR_SELEX', 'SELEX', src_type)) %>%
    mutate(name=sprintf("%s|%s[%s,%s]", str_sub(org,1,2), name, src_type, status)) %>%
    mutate(ctag = !!ctag) %>% select(ctag,motif,org,gid,everything())
toc %>% arrange(motif, org, gid) %>% print(n=20, width=Inf)
toc %>% count(src_type)
#}}}

ctag = 'planttfdb'
#{{{ planttfdb
org = 'Athaliana'
org = 'Zmays'
diri1 = file.path(diri, 'plant', org, 'TF_Information.txt')
fi = sprintf("%s/%s/%s/%s_TF_binding_motifs_information.txt", diri, ctag, org, str_sub(org,0,3))
ti = read_tsv(fi) %>%
    select(gid=Gene_id, name=Family, motif=Matrix_id, src_ref=Datasource,
           src_type=Method, src_id=Datasource_ID) %>%
    mutate(src_id = str_replace(src_id, 'transfer from ', '')) %>%
    mutate(src_id = str_replace(src_id, '\\(Arabidopsis thaliana\\)', '')) %>%
    mutate(ctag = !!ctag, org = !!org) %>%
    select(ctag, org, everything())
ti %>% count(src_type)

#{{{ maize gene ID conversion
if (org == 'Athaliana') {
    tia = ti
} else if (org == 'Zmays') {
    ti2 = ti %>% rename(ogid=gid) %>% inner_join(tm, by='ogid')
    ti2 %>% count(type)
    ti3 = ti2 %>% filter(type %in% c('1-to-1','extra')) %>%
        select(-ogid, -type) %>% select(org, gid, everything())
    tim = ti3
}
#}}}

top = rbind(tia, tim) %>% mutate(ctag = !!ctag) %>% select(ctag,everything())
top %>% arrange(motif, org, gid) %>% print(n=20, width=Inf)
#}}}

#gmtf = toc %>% bind_rows(top)
gmtf = toc
mtf = gmtf %>%
    group_by(ctag, motif) %>%
    summarise(name=str_c(name, collapse=' ')) %>% ungroup()

res = list(gmtf = gmtf, mtf = mtf)
fo = file.path(dirw, "01.rds")
saveRDS(res, fo)
#}}}

# check README.txt

#{{{ convert meme to universalmotif
fi = file.path(dirw, "01.rds")
ti = readRDS(fi)

fi = file.path(dirw, "05.motifs.meme")
tm = read_meme(fi)
x = tibble(i = 1:length(tm)) %>%
    mutate(pwm = map(i, mf <- function(i, tm) tm[[i]], tm)) %>%
    mutate(mid = map_chr(pwm, slot, name = 'name')) %>%
    select(i, mid, pwm) %>%
    full_join(ti$mtf, by=c('mid'='motif'))

to = ti
to$gmtf = ti$gmtf %>% rename(mid=motif)
to$mtf = x %>% select(ctag, mid, gene=name, pwm)
fo = file.path(dirw, '05.motifs.rds')
saveRDS(to, fo)

mids = to$mtf %>% pull(mid)
fo = file.path(dirw, '05.motifs.txt')
write(mids, file=fo)
#}}}

#{{{ collapse motifs
fi = file.path(dirw, '05.motifs.rds')
r = readRDS(fi)
fi = file.path(dirw, '05.motifs.meme')
mtfs = read_meme(fi)

cmp = compare_motifs(mtfs, method="PCC", min.mean.ic=0, score.strat="a.mean")
dst = as.dist(1 - cmp)
hcu = hclust(dst, method='ward.D')
tree = ape::as.phylo(hcu)

#{{{ distance heatmap
mids_sorted = hcu$labels[hcu$order]
tp = as_tibble(cmp) %>% mutate(mid = rownames(cmp)) %>%
    gather(mid2, score, -mid) %>%
    mutate(mid = factor(mid, levels=mids_sorted)) %>%
    mutate(mid2 = factor(mid2, levels=mids_sorted))

p = ggplot(tp) +
    geom_tile(aes(x=mid, y=mid2, fill=score)) +
    scale_x_discrete(position='bottom', expand=c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradientn(name='PCC', na.value='grey50', colors=cols100) +
    #scale_fill_viridis(name=leg) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           margin = c(.5,.5,.5,.5))
fo = file.path(dirw, '06.motif.heatmap.pdf')
ggsave(p, file=fo, width=8, height=8)
#}}}

#{{{ ggtree
tpt = ti %>% select(taxa=mtf, lgd=name)
p1 = ggtree(tree, ladderize=F) %<+%
    tpt +
    geom_tiplab(aes(label=lgd), size=1) +
    scale_x_continuous(expand=expansion(mult=c(0,.2))) +
    scale_y_continuous(expand=expansion(mult=c(.001,.001))) +
    scale_color_aaas() +
    theme_tree2() +
    theme(plot.margin = margin(.3,.5,.3,.5, 'lines')) +
    guides(color=F)
#
fo = file.path(dirw, '06.motifs.tree.pdf')
p1 %>% ggexport(filename = fo, width = 8, height = 60)
#}}}

x = cutree(hcu, h=.5)
tx = tibble(mid = names(x), grp = as.double(x)) %>%
    inner_join(r$mtf, by='mid') %>%
    mutate(icscore = map_dbl(pwm, 'icscore')) %>%
    arrange(grp, desc(icscore)) %>%
    select(mid, gene, icscore, pwm, grp) %>% print(n=30)

txs = tx %>% arrange(grp, desc(icscore)) %>%
    group_by(grp) %>%
    summarise(mids = str_c(mid, collapse=' '),
        grpname = str_c(unique(unlist(str_split(gene, ' '))), collapse=' ')) %>%
    ungroup()
tx2 = tx %>% inner_join(txs, by='grp')

fo = file.path(dirw, '09.motif.grp.rds')
res = list(gmtf = r$gmtf, mtf = tx2, grp = txs)
saveRDS(res, fo)

#{{{ distance heatmap
mids_sorted = hcu$labels[hcu$order]
tp = as_tibble(cmp) %>% mutate(mid = rownames(cmp)) %>%
    gather(mid2, score, -mid) %>%
    mutate(mid = factor(mid, levels=mids_sorted)) %>%
    mutate(mid2 = factor(mid2, levels=mids_sorted))

p = ggplot(tp) +
    geom_tile(aes(x=mid, y=mid2, fill=score)) +
    scale_x_discrete(position='bottom', expand=c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradientn(name='PCC', na.value='grey50', colors=cols100) +
    #scale_fill_viridis(name=leg) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           margin = c(.5,.5,.5,.5))
fo = file.path(dirw, '06.motif.heatmap.pdf')
ggsave(p, file=fo, width=8, height=8)
#}}}

#{{{ ggtree
x = tx2 %>% group_by(grp) %>% slice(1) %>% ungroup()
mtfs_nr = x$pwm
cmp = compare_motifs(mtfs_nr, method="PCC", min.mean.ic=0, score.strat="a.mean")
dst = as.dist(1 - cmp)
hcu = hclust(dst, method='ward.D')
tree = ape::as.phylo(hcu)

tpt = x %>% select(taxa=mid, grpname) %>% mutate(lgd=str_sub(grpname,1,50))
p1 = ggtree(tree, ladderize=F) %<+%
    tpt +
    geom_tiplab(aes(label=lgd), size=2) +
    scale_x_continuous(expand=expansion(mult=c(0,.2))) +
    scale_y_continuous(expand=expansion(mult=c(.01,.01))) +
    scale_color_aaas() +
    theme_tree2() +
    theme(plot.margin = margin(.3,.5,.3,.5, 'lines')) +
    guides(color=F)
#
fo = file.path(dirw, '09.motif.grp.tree.pdf')
p1 %>% ggexport(filename = fo, width = 8, height = 20)
#}}}


#}}}


