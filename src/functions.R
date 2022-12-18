require(devtools)
load_all('~/git/rmaize')
require(ggpubr)
require(ggforce)
require(knitr)
require(kableExtra)
options(knitr.table.format = "latex")
#require(pROC)
#require(PRROC)
require(universalmotif)
dirg = '~/data/genome'
dirp = '~/projects/cre'
dird = file.path(dirp, 'data')
gcfg = read_genome_conf()
#
#tsyn = read_syn(gcfg)
#symb = read_symbol(opt='') %>% mutate(symbol=ifelse(gid=='Zm00001d026147','R1',symbol))
cols100 = colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(100)
cols100v = viridis_pal(direction=-1,option='magma')(100)


