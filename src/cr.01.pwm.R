source("functions.R")
require(ggtree)
dirw = glue('{dird}/01_tfbs')
diri = glue('{dirw}/raw')
mdic = v3_to_v4(opt='dict_all')
ndic = read_symbol(opt = 'dict')
rename_mtf <- function(mtf, name) { mtf['name'] = name; mtf }

#{{{ cisbp
ctag = 'cisbp'
#{{{ read cis-BP motifs
org = 'Athaliana'
org = 'Zmays'
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
#
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
    select(motif,org,gid,everything())
toc %>% arrange(motif, org, gid) %>% print(n=20, width=Inf)
toc %>% count(src_type)
#}}}

fi = glue("{diri}/{ctag}/05.motifs.meme")
mtfs = read_meme(fi)
tm0 = tibble(i = 1:length(tm)) %>%
    mutate(pwm = map(i, mf <- function(i, mtfs) tm[[i]], tm)) %>%
    mutate(mid = map_chr(pwm, slot, name='name')) %>%
    select(i, mid, pwm)

gmtf = toc %>% rename(mid=motif)
tm = gmtf %>% group_by(mid) %>%
    summarise(name=str_c(name, collapse=' ')) %>% ungroup() %>%
    full_join(tm0, by=c('mid')) %>%
    select(mid, note=name, pwm)

r = list(gmtf = gmtf, mtf = tm)
fo = glue("{dirw}/01.{ctag}.rds")
saveRDS(r, fo)

#mids = to$mtf %>% pull(mid)
#fo = file.path(dirw, '05.motifs.txt')
#write(mids, file=fo)
#}}}

#{{{ planttfdb #-> don't use this
ctag = 'planttfdb'
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

#{{{ maize chipseq + dapseq
#{{{ find alias for 104 chipseq TFs
fi = glue('{dirw}/maize_tf_names.xlsx')
ti = read_xlsx(fi, sheet='Sheet2', col_names=c('tf'))
ta = read_symbol(opt='tibble')
to = ti %>% left_join(ta, by=c('tf'='gid'))
#write_tsv(to, glue('{dirw}/tmp.tsv'))
#}}}

fi = glue('{dirw}/maize_tf_names.xlsx')
ti1 = read_xlsx(fi, sheet='Sheet1', col_names=c('tf')) %>%
    mutate(note = str_to_lower(tf)) %>%
    mutate(ctag='ricci2019')
ti2 = read_xlsx(fi, sheet='Sheet2', col_names=c('tf','note')) %>%
    mutate(ctag='tu2020')
ti = rbind(ti1,ti2) %>% select(ctag, note, tf)

get_first_mtf <- function(mtfs) if(length(mtfs) > 1) mtfs[[1]] else mtfs
to = ti %>%
    mutate(fi=glue("{diri}/maize_tfs/02_transfac/{tf}.transfac")) %>%
    mutate(x = map(fi, read_transfac)) %>%
    mutate(pwm = map(x, get_first_mtf)) %>%
    select(ctag, note, pwm)

fo = glue("{dirw}/01.maize.rds")
saveRDS(to, fo)
#}}}

#{{{ integrate different sources of TFBS PWMs
fi = glue("{dirw}/01.cisbp.rds")
r = readRDS(fi)
tm1 = r$mtf %>% mutate(name = note) %>%
    mutate(name=str_replace_all(name, ".*\\|", "")) %>%
    mutate(name=str_replace_all(name, "\\[.*", "")) %>%
    #mutate(fname=str_replace_all(fname, "-", "")) %>%
    mutate(ctag='cisbp') %>% select(ctag, mid, name, pwm)

fi = glue("{dirw}/01.maize.rds")
r2 = readRDS(fi)
tm2 = r2 %>% mutate(ctag='maize') %>% select(ctag,name=note,pwm)

rename_mtf <- function(mtf, name) { mtf['name'] = name; mtf }
tm = tm1 %>% bind_rows(tm2) %>%
    mutate(mid = str_c('m', str_pad(1:n(), pad='0', width=4))) %>%
    mutate(pwm = map2(pwm, mid, rename_mtf))

fo = glue("{dirw}/05.motifs.rds")
saveRDS(tm, fo)
#}}}



#{{{ # distance heatmap
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

#{{{ collapse motifs (two-pass clustering)
#{{{ trimming
trim_motifs2 <- function(mtf)
    tryCatch(trim_motifs(mtf), error = function(c) NULL)
fi = glue('{dirw}/05.motifs.rds')
tm = readRDS(fi) %>%
    mutate(pwm = map(pwm, trim_motifs2)) %>%
    mutate(j = map_dbl(pwm, length)) %>% filter(j>0) %>% select(-j) %>%
    mutate(conseq=map_chr(pwm, 'consensus')) %>%
    mutate(nsites=map_int(conseq, nchar)) %>%
    filter(nsites >= 5)
#}}}

mtfs = tm$pwm
#{{{ # 1-pass clustering
cmp0 = compare_motifs(mtfs, method="PCC", min.mean.ic=.0,
                      min.overlap=5, score.strat="a.mean")
dst0 = as.dist(1 - cmp0)
hcu0 = hclust(dst0, method='average')

#x = cutreeHybrid(hcu, as.matrix(dst), deepSplit=2, minClusterSize=1, pamStage=T, maxPamDist=.15)
#tx0 = tibble(mid = hcu$labels, grp = as.double(x$labels))
x = cutree(hcu0, h=.05)
tx0 = tibble(mid = hcu0$labels, grp = as.integer(x))
tx0$grp[tx0$grp==0] = max(tx0$grp) + 1:sum(tx0$grp==0)
tx = tx0 %>%
    inner_join(tm, by='mid') %>%
    mutate(icscore = map_dbl(pwm, 'icscore')) %>%
    arrange(grp, desc(icscore)) %>%
    select(mid, name, icscore, grp, conseq, pwm) %>% print(n=30)
#
txs = tx %>% arrange(grp, desc(icscore)) %>%
    group_by(grp) %>%
    summarise(mids = str_c(mid, collapse=' '),
        grpname = str_c(unique(name), collapse=' ')) %>%
    ungroup()
tx = tx %>% inner_join(txs, by='grp')
tx %>% count(grp) %>% arrange(desc(grp))
#tx %>% filter(grp==4) %>% print(n=30)
#}}}
# tx txs

#{{{ # 2-pass clustering
tx2 = tx %>% arrange(grp, desc(icscore)) %>%
    group_by(grp) %>%
    summarise(mid=mid[1], pwm=pwm[1]) %>% ungroup()
#
cmp1 = compare_motifs(tx2$pwm, method="PCC", min.mean.ic=.0,
                      min.overlap=6, score.strat="a.mean")
dst = as.dist(1 - cmp1)
hcu = hclust(dst, method='average')

x = cutree(hcu, h=.05)
tx3 = tibble(mid = hcu$labels, grp2 = as.integer(x)) %>%
    inner_join(tx2, by='mid') %>% select(grp, grp2)
#
tr = tx %>% inner_join(tx3, by='grp') %>%
    select(-grp) %>% rename(grp=grp2) %>%
    select(-mids, -grpname) %>%
    arrange(grp, desc(icscore))
#}}}
# tr

rename_mtf <- function(mtf, name) { mtf['name'] = name; mtf }
tmc = tm %>% left_join(tr %>% select(mid,icscore,grp), by='mid') %>%
    arrange(grp, desc(icscore)) %>% rename(name0=name) %>%
    group_by(grp) %>%
    summarise(mid = mid[1], name = name0[1],
              conseq = conseq[1], pwm=pwm[1],
              n_mid = n(),
              n_cisbp=sum(ctag=='cisbp'), n_maize=sum(ctag=='maize'),
              mids = str_c(mid, collapse=' '),
              name_l = str_c(unique(name0), collapse=' ')) %>%
    ungroup()
tmc %>% filter(n_maize>0) %>% select(name,name_l,n_mid,n_cisbp,n_maize) %>% print(n=60)
#tmc %>% count(grp) %>% arrange(desc(grp))

#{{{ # similarity among motif clusters
mtfs = trs$pwm
cmp = compare_motifs(mtfs, method="PCC", min.mean.ic=0,
                     min.overlap = 6, score.strat="a.mean")
dst = as.dist(1 - cmp)
hcu = hclust(dst, method='average')

#{{{ distance heatmap
mids_sorted = hcu$labels[hcu$order]
tp = as_tibble(cmp) %>% mutate(mid = rownames(cmp)) %>%
    gather(mid2, score, -mid) %>%
    mutate(mid = factor(mid, levels=mids_sorted)) %>%
    mutate(mid2 = factor(mid2, levels=mids_sorted))
#
p = ggplot(tp) +
    geom_tile(aes(x=mid, y=mid2, fill=score)) +
    scale_x_discrete(position='bottom', expand=c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradientn(name='PCC', na.value='grey50', colors=cols100) +
    #scale_fill_viridis(name=leg) +
    otheme(legend.pos='none', legend.dir='h', legend.title=T,
           margin = c(.5,.5,.5,.5))
fo = file.path(dirw, '10.heatmap.pdf')
ggsave(p, file=fo, width=8, height=8)
#}}}

#{{{ ggtree for all motifs
tree = ape::as.phylo(hcu0)
tpt = tr %>% select(taxa=mid, grp,gene) %>%
    mutate(lgd=str_c(grp, str_sub(gene,1,20), sep=' '))
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
fo = file.path(dirw, '10.tree.all.pdf')
p1 %>% ggexport(filename = fo, width = 8, height = 60)
#}}}

#{{{ ggtree for motif centroids
tree = ape::as.phylo(hcu)
tpt = trs %>% select(taxa=mid, grp,grpname) %>%
    mutate(lgd=str_c(grp, str_sub(grpname,1,20), sep=' '))
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
fo = file.path(dirw, '10.tree.pdf')
p1 %>% ggexport(filename = fo, width = 8, height = 20)
#}}}
#}}}

fo = glue('{dirw}/10.fam.rds')
saveRDS(tmc, fo)
fo = glue('{dirw}/10.fam.meme')
write_meme(tmc$pwm, fo, overwrite=T)
fo = glue('{dirw}/10.fam.tsv')
fo = file.path(dirw, '10.fam.tsv')
write_tsv(tmc %>% select(-pwm), fo)
#}}}
