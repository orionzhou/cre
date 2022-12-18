source("functions.R")
#require(ggtree)
dirw = glue('{dird}/21_nam_pro')

panels = read_gt_panels()
panels
nam25 = get_gt_panel('nam25', panels)
maize31 = get_gt_panel('maize31', panels)
maize32_ph207 = get_gt_panel('maize32_ph207', panels) 

diri = '~/projects/dnaseq/data/Zmays_B73v5/'
fm = glue('{diri}/10_geno_lists/vt10.tsv')
tm = read_tsv(fm)
fi = glue("{diri}/01.meta.tsv")
th = read_tsv(fi)
ths = th %>% mutate(reads = mapped_uniq/1e6) %>% select(yid, gt, reads)

tm1 = tm %>% filter(gt %in% c(nam25, 'B73')) %>%
    inner_join(ths, by=c('yid','gt')) %>%
    arrange(gt, desc(reads)) %>%
    group_by(gt) %>% slice(1) %>% ungroup()

fo = glue("{dirw}/01.gts.txt")
to = tm1 %>% mutate(gt0 = glue("{yid}#{gt}")) %>%
    select(gt0)
write_tsv(to, fo, col_names=F)
