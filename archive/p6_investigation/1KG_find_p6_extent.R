# how can we look at the haplotypes in the 1KG data
# and find the extent of the p6 haplotype
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggpubr)
# i need random access to the file, reading certain rows by a key so that I can iterate over them
# that would be ideal, but I guess I can just start after the p6 and read rows iteratively

# sample_names <- c(read.table('~/data/1KG_data/individual_chromosomes/keep_p6_hom_hets.txt')[,1])
sample_names_p6_any <- c('HG00107','HG00110','HG00112','HG00113','HG00115','HG00118','HG00119','HG00127','HG00132','HG00136','HG00137','HG00139','HG00140','HG00141','HG00142','HG00143','HG00145','HG00151','HG00154','HG00155','HG00171','HG00174',
                         'HG00177','HG00178','HG00179','HG00185','HG00186','HG00187','HG00190','HG00232','HG00234','HG00235','HG00238','HG00239','HG00245','HG00246','HG00250','HG00251','HG00255','HG00256','HG00257','HG00258','HG00259','HG00261',
                         'HG00262','HG00265','HG00269','HG00271','HG00274','HG00275','HG00276','HG00281','HG00282','HG00284','HG00290','HG00308','HG00309','HG00315','HG00320','HG00321','HG00326','HG00327','HG00329','HG00331','HG00335','HG00337',
                         'HG00341','HG00343','HG00345','HG00350','HG00356','HG00358','HG00360','HG00364','HG00366','HG00367','HG00368','HG00369','HG00371','HG00376','HG00380','HG00381','HG00382','HG00383','HG00406','HG00407','HG00410','HG00422',
                         'HG00428','HG00436','HG00442','HG00446','HG00451','HG00452','HG00457','HG00464','HG00479','HG00524','HG00525','HG00534','HG00536','HG00551','HG00553','HG00556','HG00560','HG00565','HG00581','HG00584','HG00589','HG00592',
                         'HG00598','HG00599','HG00607','HG00610','HG00611','HG00614','HG00620','HG00623','HG00625','HG00626','HG00629','HG00637','HG00638','HG00651','HG00663','HG00672','HG00674','HG00675','HG00684','HG00689','HG00693','HG00701',
                         'HG00708','HG00729','HG00736','HG00743','HG00759','HG00851','HG00881','HG00956','HG00982','HG01031','HG01047','HG01051','HG01052','HG01054','HG01055','HG01060','HG01072','HG01085','HG01088','HG01089','HG01092','HG01095',
                         'HG01097','HG01098','HG01101','HG01102','HG01104','HG01105','HG01107','HG01110','HG01112','HG01119','HG01130','HG01131','HG01136','HG01137','HG01139','HG01148','HG01149','HG01161','HG01170','HG01173','HG01176','HG01177',
                         'HG01188','HG01191','HG01242','HG01247','HG01256','HG01257','HG01259','HG01260','HG01269','HG01271','HG01272','HG01277','HG01280','HG01305','HG01311','HG01325','HG01326','HG01334','HG01344','HG01348','HG01351','HG01365',
                         'HG01369','HG01375','HG01383','HG01389','HG01390','HG01393','HG01402','HG01412','HG01414','HG01432','HG01435','HG01440','HG01441','HG01443','HG01455','HG01456','HG01462','HG01464','HG01465','HG01468','HG01474','HG01479',
                         'HG01492','HG01494','HG01495','HG01497','HG01501','HG01504','HG01509','HG01510','HG01513','HG01518','HG01519','HG01522','HG01524','HG01525','HG01528','HG01530','HG01531','HG01551','HG01556','HG01565','HG01566','HG01571',
                         'HG01572','HG01577','HG01578','HG01583','HG01586','HG01595','HG01596','HG01599','HG01600','HG01605','HG01612','HG01613','HG01619','HG01620','HG01624','HG01625','HG01630','HG01632','HG01670','HG01673','HG01684','HG01685',
                         'HG01694','HG01704','HG01705','HG01756','HG01761','HG01766','HG01777','HG01783','HG01784','HG01790','HG01794','HG01797','HG01798','HG01802','HG01804','HG01808','HG01815','HG01817','HG01841','HG01844','HG01845','HG01846',
                         'HG01847','HG01848','HG01849','HG01850','HG01852','HG01853','HG01858','HG01861','HG01865','HG01866','HG01868','HG01869','HG01870','HG01872','HG01873','HG01874','HG01878','HG01880','HG01894','HG01896','HG01917','HG01918',
                         'HG01920','HG01921','HG01924','HG01926','HG01927','HG01932','HG01933','HG01936','HG01938','HG01939','HG01941','HG01944','HG01945','HG01947','HG01950','HG01953','HG01961','HG01967','HG01968','HG01971','HG01974','HG01979',
                         'HG01980','HG01985','HG01988','HG01990','HG02003','HG02006','HG02008','HG02009','HG02014','HG02016','HG02017','HG02019','HG02029','HG02035','HG02048','HG02051','HG02057','HG02058','HG02060','HG02064','HG02069','HG02070',
                         'HG02079','HG02084','HG02086','HG02089','HG02090','HG02102','HG02105','HG02113','HG02116','HG02121','HG02122','HG02130','HG02131','HG02134','HG02136','HG02137','HG02141','HG02146','HG02147','HG02150','HG02153','HG02154',
                         'HG02155','HG02164','HG02166','HG02178','HG02180','HG02184','HG02186','HG02187','HG02190','HG02219','HG02220','HG02232','HG02236','HG02250','HG02252','HG02253','HG02259','HG02260','HG02265','HG02266','HG02271','HG02272',
                         'HG02274','HG02275','HG02277','HG02278','HG02282','HG02285','HG02286','HG02292','HG02299','HG02301','HG02304','HG02308','HG02312','HG02318','HG02330','HG02334','HG02345','HG02348','HG02351','HG02373','HG02374','HG02379',
                         'HG02382','HG02383','HG02384','HG02386','HG02390','HG02394','HG02396','HG02398','HG02408','HG02409','HG02410','HG02427','HG02470','HG02471','HG02476','HG02489','HG02494','HG02521','HG02522','HG02537','HG02573','HG02574',
                         'HG02580','HG02595','HG02597','HG02600','HG02601','HG02604','HG02610','HG02628','HG02643','HG02649','HG02655','HG02660','HG02661','HG02666','HG02679','HG02682','HG02684','HG02685','HG02687','HG02690','HG02691','HG02694',
                         'HG02697','HG02699','HG02716','HG02722','HG02725','HG02727','HG02734','HG02736','HG02737','HG02763','HG02769','HG02771','HG02774','HG02778','HG02780','HG02784','HG02786','HG02789','HG02793','HG02811','HG02814','HG02817',
                         'HG02860','HG02861','HG02885','HG02923','HG02941','HG02946','HG03006','HG03018','HG03021','HG03022','HG03027','HG03028','HG03045','HG03049','HG03054','HG03066','HG03074','HG03085','HG03086','HG03103','HG03114','HG03136',
                         'HG03139','HG03157','HG03160','HG03175','HG03228','HG03229','HG03237','HG03238','HG03258','HG03301','HG03313','HG03351','HG03378','HG03394','HG03446','HG03455','HG03457','HG03472','HG03478','HG03484','HG03488','HG03490',
                         'HG03499','HG03520','HG03547','HG03558','HG03585','HG03593','HG03594','HG03595','HG03600','HG03604','HG03611','HG03615','HG03616','HG03629','HG03631','HG03634','HG03642','HG03643','HG03644','HG03646','HG03649','HG03652',
                         'HG03660','HG03668','HG03673','HG03679','HG03686','HG03687','HG03691','HG03696','HG03702','HG03705','HG03708','HG03709','HG03713','HG03714','HG03716','HG03718','HG03731','HG03733','HG03736','HG03738','HG03740','HG03743',
                         'HG03746','HG03755','HG03756','HG03760','HG03762','HG03767','HG03770','HG03771','HG03772','HG03774','HG03775','HG03777','HG03779','HG03785','HG03787','HG03789','HG03790','HG03792','HG03793','HG03802','HG03805','HG03809',
                         'HG03823','HG03829','HG03836','HG03854','HG03862','HG03868','HG03871','HG03872','HG03874','HG03875','HG03884','HG03890','HG03897','HG03899','HG03905','HG03907','HG03914','HG03917','HG03922','HG03926','HG03928','HG03940',
                         'HG03947','HG03949','HG03951','HG03960','HG03963','HG03967','HG03968','HG03971','HG03973','HG03991','HG04002','HG04014','HG04015','HG04018','HG04019','HG04020','HG04023','HG04035','HG04054','HG04059','HG04060','HG04061',
                         'HG04063','HG04098','HG04100','HG04106','HG04131','HG04134','HG04140','HG04141','HG04144','HG04152','HG04153','HG04155','HG04156','HG04162','HG04164','HG04173','HG04183','HG04185','HG04186','HG04189','HG04194','HG04195',
                         'HG04202','HG04209','HG04211','HG04219','HG04222','HG04227','HG04229','HG04235','NA06986','NA06989','NA07000','NA07037','NA07048','NA07051','NA07056','NA07347','NA11830','NA11840','NA11881','NA11892','NA11893','NA11894',
                         'NA11918','NA11919','NA11930','NA11933','NA11992','NA11994','NA12003','NA12043','NA12045','NA12058','NA12144','NA12154','NA12275','NA12283','NA12341','NA12348','NA12383','NA12414','NA12489','NA12546','NA12718','NA12748',
                         'NA12762','NA12775','NA12776','NA12777','NA12812','NA12813','NA12814','NA12815','NA12827','NA12828','NA12830','NA12842','NA12872','NA12873','NA12874','NA18499','NA18502','NA18504','NA18511','NA18519','NA18530','NA18532',
                         'NA18534','NA18538','NA18548','NA18553','NA18559','NA18567','NA18570','NA18577','NA18591','NA18595','NA18597','NA18599','NA18605','NA18608','NA18610','NA18611','NA18612','NA18615','NA18616','NA18618','NA18624','NA18628',
                         'NA18629','NA18630','NA18631','NA18632','NA18633','NA18634','NA18637','NA18640','NA18645','NA18646','NA18648','NA18740','NA18747','NA18865','NA18870','NA18876','NA18879','NA18923','NA18940','NA18943','NA18947','NA18949',
                         'NA18952','NA18954','NA18956','NA18959','NA18961','NA18962','NA18964','NA18966','NA18967','NA18971','NA18972','NA18973','NA18974','NA18976','NA18977','NA18981','NA18984','NA18985','NA18986','NA18991','NA18999','NA19000',
                         'NA19006','NA19007','NA19009','NA19010','NA19012','NA19025','NA19027','NA19030','NA19031','NA19037','NA19054','NA19058','NA19060','NA19068','NA19072','NA19074','NA19075','NA19079','NA19081','NA19083','NA19084','NA19085',
                         'NA19088','NA19091','NA19096','NA19107','NA19113','NA19121','NA19129','NA19141','NA19190','NA19201','NA19204','NA19206','NA19214','NA19223','NA19236','NA19238','NA19308','NA19315','NA19350','NA19360','NA19372','NA19391',
                         'NA19428','NA19445','NA19454','NA19455','NA19463','NA19467','NA19471','NA19648','NA19649','NA19651','NA19652','NA19654','NA19657','NA19663','NA19664','NA19669','NA19676','NA19679','NA19700','NA19701','NA19703','NA19717',
                         'NA19719','NA19720','NA19723','NA19726','NA19731','NA19732','NA19734','NA19740','NA19741','NA19746','NA19747','NA19750','NA19752','NA19755','NA19756','NA19758','NA19761','NA19764','NA19771','NA19773','NA19774','NA19777',
                         'NA19780','NA19783','NA19785','NA19792','NA19835','NA19900','NA19901','NA19922','NA20127','NA20281','NA20282','NA20299','NA20314','NA20339','NA20346','NA20348','NA20412','NA20502','NA20503','NA20504','NA20505','NA20507',
                         'NA20509','NA20518','NA20519','NA20522','NA20524','NA20525','NA20528','NA20529','NA20531','NA20533','NA20534','NA20535','NA20541','NA20542','NA20543','NA20581','NA20582','NA20585','NA20586','NA20753','NA20754','NA20757',
                         'NA20759','NA20760','NA20761','NA20762','NA20766','NA20768','NA20771','NA20772','NA20774','NA20775','NA20786','NA20787','NA20790','NA20792','NA20796','NA20799','NA20804','NA20808','NA20810','NA20812','NA20815','NA20818',
                         'NA20821','NA20826','NA20827','NA20832','NA20845','NA20854','NA20856','NA20859','NA20864','NA20866','NA20867','NA20868','NA20869','NA20872','NA20882','NA20886','NA20894','NA20897','NA20904','NA20905','NA21086','NA21090',
                         'NA21091','NA21094','NA21095','NA21098','NA21100','NA21101','NA21106','NA21108','NA21110','NA21111','NA21115','NA21117','NA21118','NA21119','NA21122','NA21126','NA21127','NA21129','NA21130','NA21133','NA21135','NA21144')
sample_names_p6_hom <- c('HG00107','HG00140','HG00154','HG00186','HG00232','HG00250','HG00255','HG00262','HG00281','HG00321','HG00326','HG00335','HG00341','HG00366','HG00367','HG00381','HG00382','HG00407','HG00436','HG00525','HG00536','HG00592','HG00607','HG00623','HG00708','HG01136','HG01148','HG01149','HG01277','HG01280','HG01311','HG01344','HG01383','HG01393','HG01464','HG01474','HG01479','HG01497','HG01501','HG01509','HG01530','HG01565','HG01572','HG01586','HG01599','HG01670','HG01704','HG01766','HG01790','HG01844','HG01845','HG01848','HG01850','HG01861','HG01868','HG01872','HG01896','HG01918','HG01921','HG01938','HG01941','HG01947','HG01961','HG01968','HG01974','HG01988','HG02003','HG02019','HG02070','HG02113','HG02131','HG02134','HG02146','HG02178','HG02219','HG02260','HG02265','HG02277','HG02278','HG02286','HG02292','HG02304','HG02308','HG02345','HG02379','HG02394','HG02410','HG02522','HG02601','HG02628','HG02660','HG02684','HG02694','HG02727','HG02778','HG02861','HG03228','HG03229','HG03237','HG03446','HG03585','HG03615','HG03631','HG03642','HG03731','HG03738','HG03755','HG03767','HG03771','HG03772','HG03777','HG03809','HG03862','HG03899','HG03967','HG04020','HG04059','NA06989','NA07056','NA11892','NA11919','NA11930','NA12043','NA12762','NA12776','NA18608','NA18616','NA18632','NA18962','NA18966','NA18977','NA18984','NA18985','NA19072','NA19085','NA19308','NA19372','NA19467','NA19471','NA19654','NA19676','NA19679','NA19752','NA19756','NA19758','NA19773','NA19777','NA20534','NA20542','NA20543','NA20757','NA20768','NA20799','NA20872','NA20894','NA20897','NA20904','NA20905','NA21086')
sample_names_p6_het <- setdiff(sample_names_p6_any, sample_names_p6_hom)
p6_variant_bed <- read.table('~/local_data/p6/p6_targets_onebp.bed', sep='\t', quote='', header=F)
GGCC_variant_bed <- read.table('~/local_data/p6/GGCC_targets_onebp.bed', sep='\t', quote='', header=F)

jak2_exon_bed <- read.table('~/local_data/p6/jak2_exons.bed', sep='\t',quote='')
colnames(jak2_exon_bed) <- c('chr', 'start', 'stop')

# get the starting position of the p6
found_p6 <- F
start_pos <- 0
read_nrow <- 10000

vcf_f <- gzfile('~/local_data/1KG/individual_chromosomes/1KG_allvar_p6_hom_hets.vcf.gz',open = 'rt')
while(! found_p6){
    print(start_pos)
    vcf <- read.table(textConnection(readLines(vcf_f, n = read_nrow)),
                      header = F, sep='\t', quote='')
    print(vcf[1,2])
    which_p6 <- which(vcf$V2 == p6_variant_bed[1,2])
    if(length(which_p6)>0){
        found_p6 <- T
        load_offset <- start_pos + which_p6 -1
        print(load_offset)
        close(vcf_f)
    } else{
        start_pos <- start_pos + read_nrow
    }
}

# result of the above 
start_pos <- 2e+05

# read the next n rows as a starting place
vcf_f <- gzfile('~/local_data/1KG/individual_chromosomes/1KG_allvar_p6_hom_hets.vcf.gz',open = 'rt')
vcf_post <- read.table(vcf_f,
                header = F, sep='\t', quote='',
                skip = load_offset, nrows = 20000)
close(vcf_f)
vcf_f <- gzfile('~/local_data/1KG/individual_chromosomes/1KG_allvar_p6_hom_hets.vcf.gz',open = 'rt')
vcf_pre <- read.table(vcf_f,
                      header = F, sep='\t', quote='',
                      skip = load_offset-20000, nrows = 20000)
close(vcf_f)
vcf <- rbind(vcf_pre, vcf_post)

vcf_meta <- vcf[,1:9]
vcf_data <- vcf[,10:ncol(vcf)]
colnames(vcf_data) <- sample_names
colnames(vcf_meta) <- c('chr', 'pos', 'id', 'ref','alt', 'gq', 'filter', 'm1','m2')
rownames(vcf_data) <- paste(vcf_meta$pos, vcf_meta$ref, vcf_meta$alt, sep='_')
rownames(vcf_meta) <- paste(vcf_meta$pos, vcf_meta$ref, vcf_meta$alt, sep='_')
vcf_data_hom <- vcf_data[, sample_names_p6_hom]
rownames(vcf_data_hom) <- paste(vcf_meta$pos, vcf_meta$ref, vcf_meta$alt, sep='_')
rm(vcf)

# eliminate rows that are all 0|0
vcf_data <- vcf_data[apply(vcf_data, 1, function(x) !all(x=='0|0')),]
vcf_data_hom <- vcf_data_hom[apply(vcf_data_hom, 1, function(x) !all(x=='0|0')),]

# get integers for each phased allele
a1_from_g <- function(x){
    if(x=='0|0' | x=='0|1'){
        return(0)
    } else {
        return(1)
    }
}
a2_from_g <- function(x){
    if(x=='0|0' | x=='1|0'){
        return(0)
    } else {
        return(1)
    }
}
# a1 <- apply(vcf_data, 1:2, function(x) as.numeric(strsplit(x, '|')[[1]][1]))
a1 <- apply(vcf_data, 1:2, a1_from_g)
a2 <- apply(vcf_data, 1:2, a2_from_g)
a1_hom <- apply(vcf_data_hom, 1:2, a1_from_g)
a2_hom <- apply(vcf_data_hom, 1:2, a2_from_g)

# For the hets, I need to find which chr their p6 is on, and look along that chr
# easier to start with the homs
# remove rows all zero
target_variant <- '5050706_C_T'
target_variant_GGCC <- '5068755_T_G'
min_var <- 2
a1_hom2 <- a1_hom[rowSums(a1_hom) >= min_var, ]
a1_hom_dist <- as.matrix(dist((a1_hom2), method = 'manh'))
a2_hom2 <- a2_hom[rowSums(a2_hom) >= min_var, ]
a2_hom_dist <- as.matrix(dist((a2_hom2), method = 'manh'))

plot_df_1 <- data.frame(id=colnames(a1_hom_dist), 
                      n=1:ncol(a1_hom_dist),
                      pos=vcf_meta[colnames(a1_hom_dist), "pos"],
                      dist=a1_hom_dist[,target_variant],
                      dist_rel=1-a1_hom_dist[,target_variant]/ncol(a1_hom))
plot_df_2 <- data.frame(id=colnames(a2_hom_dist), 
                      n=1:ncol(a2_hom_dist),
                      pos=vcf_meta[colnames(a2_hom_dist), "pos"],
                      dist=a2_hom_dist[,target_variant], 
                      dist_rel=1-a2_hom_dist[,target_variant]/ncol(a2_hom))
plot_df_1 <- plot_df_1[!duplicated(plot_df_1$pos), ]
plot_df_2 <- plot_df_2[!duplicated(plot_df_2$pos), ]

p1 <- ggplot(plot_df_1, aes(x=pos, y=dist_rel)) + 
    geom_point() + 
    labs(title='p6 homozygotes: allele 1', x='chr9 position', y='Fraction of samples identical') + 
    geom_vline(xintercept = p6_variant_bed$V2, col='firebrick')
p2 <- ggplot(plot_df_2, aes(x=pos, y=dist_rel)) + 
    geom_point() +
    labs(title='p6 homozygotes: allele 2', x='chr9 position', y='Fraction of samples identical') +
    geom_vline(xintercept = p6_variant_bed$V2, col='firebrick')

plot_grid(p1,p2, ncol=1)

common_points <- intersect(plot_df_1$pos, plot_df_2$pos)
plot_df_combined <- plot_df_1[plot_df_1$pos %in% common_points, ]
plot_df_combined$dist <- rowMeans(cbind(plot_df_1[plot_df_1$pos %in% common_points, "dist"],
                              plot_df_2[plot_df_2$pos %in% common_points, "dist"]))
plot_df_combined$dist_rel <- 1- plot_df_combined$dist/ncol(a2_hom)

p3 <- ggplot(plot_df_combined, aes(x=pos, y=dist_rel)) + 
    geom_point() + 
    labs(title='p6 homozygotes: haplotype extension', x='chr9 position', y='Fraction of samples identical') + 
    geom_vline(xintercept = p6_variant_bed$V2, col='firebrick') +
    geom_vline(xintercept = GGCC_variant_bed$V2, col='orange') +
    geom_vline(xintercept = c(4995000, 5265000), col='blue') +
    geom_segment(data=jak2_exon_bed, aes(y=-0.1, yend=-0.1, x=start, xend=stop+1000), lwd=2, col='darkgreen') +
    xlim(c(4980000, 5300000))
p3


# # # # # # # # # # # # # # #
# Find this for the p6 hets #
# # # # # # # # # # # # # # #

# identify for each subject, which allele the p6 is on
# just look at the 5 snps for now
p6_variant_strings <- rownames(vcf_meta)[vcf_meta$pos %in% p6_variant_bed$V2]
p6_variant_match <- c(1,1,1,0,1)
names(p6_variant_match) <- p6_variant_strings

a1_het <- a1[, sample_names_p6_het]
a2_het <- a2[, sample_names_p6_het]
a1_p6_inds <- which(apply(a1_het[p6_variant_strings, ],2, function(x) identical(x, p6_variant_match)))
a2_p6_inds <- which(apply(a2_het[p6_variant_strings, ],2, function(x) identical(x, p6_variant_match)))

# combine these 
combined_het <- cbind(a1_het[, a1_p6_inds], a2_het[, a2_p6_inds])

# do the same search
# remove rows all zero
target_variant <- '5050706_C_T'
min_var <- 2
combined_het2 <- combined_het[rowSums(combined_het) >= min_var, ]
combined_het_dist <- as.matrix(dist((combined_het2), method = 'manh'))

plot_df_het <- data.frame(id=colnames(combined_het_dist), 
                        n=1:ncol(combined_het_dist),
                        pos=vcf_meta[colnames(combined_het_dist), "pos"],
                        dist=combined_het_dist[,target_variant],
                        dist_rel=1-combined_het_dist[,target_variant]/ncol(combined_het))
plot_df_het <- plot_df_het[!duplicated(plot_df_het$pos), ]

p4 <- ggplot(plot_df_het, aes(x=pos, y=dist_rel)) + 
    geom_point() + 
    labs(title='p6 heterozygotes: haplotype extension', x='chr9 position', y='Fraction of samples identical') + 
    geom_vline(xintercept = p6_variant_bed$V2, col='firebrick') +
    geom_vline(xintercept = GGCC_variant_bed$V2, col='orange') +
    geom_vline(xintercept = c(4995000, 5265000), col='blue') + 
    geom_segment(data=jak2_exon_bed, aes(y=-0.1, yend=-0.1, x=start, xend=stop+1000), lwd=2, col='darkgreen') +
    xlim(c(4980000, 5300000))
p4

# save a combined plot
pdf('~/local_data/p6/p6_extent_1KG.pdf', height=6, width = 10)
plot_grid(p3, p4, nrow=2)
dev.off()


##################################################################################################
# Do the same but for the GGCC  ##################################################################
##################################################################################################
target_variant <- '5068755_T_G'
min_var <- 2
a1_hom2 <- a1_hom[rowSums(a1_hom) >= min_var, ]
a1_hom_dist <- as.matrix(dist((a1_hom2), method = 'manh'))
a2_hom2 <- a2_hom[rowSums(a2_hom) >= min_var, ]
a2_hom_dist <- as.matrix(dist((a2_hom2), method = 'manh'))

plot_df_1 <- data.frame(id=colnames(a1_hom_dist), 
                        n=1:ncol(a1_hom_dist),
                        pos=vcf_meta[colnames(a1_hom_dist), "pos"],
                        dist=a1_hom_dist[,target_variant],
                        dist_rel=1-a1_hom_dist[,target_variant]/ncol(a1_hom))
plot_df_2 <- data.frame(id=colnames(a2_hom_dist), 
                        n=1:ncol(a2_hom_dist),
                        pos=vcf_meta[colnames(a2_hom_dist), "pos"],
                        dist=a2_hom_dist[,target_variant], 
                        dist_rel=1-a2_hom_dist[,target_variant]/ncol(a2_hom))
plot_df_1 <- plot_df_1[!duplicated(plot_df_1$pos), ]
plot_df_2 <- plot_df_2[!duplicated(plot_df_2$pos), ]

p1 <- ggplot(plot_df_1, aes(x=pos, y=dist_rel)) + 
    geom_point() + 
    labs(title='GGCC homozygotes: allele 1', x='chr9 position', y='Fraction of samples identical') + 
    geom_vline(xintercept = GGCC_variant_bed$V2, col='firebrick')
p2 <- ggplot(plot_df_2, aes(x=pos, y=dist_rel)) + 
    geom_point() +
    labs(title='GGCC homozygotes: allele 2', x='chr9 position', y='Fraction of samples identical') +
    geom_vline(xintercept = GGCC_variant_bed$V2, col='firebrick')

plot_grid(p1,p2, ncol=1)

common_points <- intersect(plot_df_1$pos, plot_df_2$pos)
plot_df_combined <- plot_df_1[plot_df_1$pos %in% common_points, ]
plot_df_combined$dist <- rowMeans(cbind(plot_df_1[plot_df_1$pos %in% common_points, "dist"],
                                        plot_df_2[plot_df_2$pos %in% common_points, "dist"]))
plot_df_combined$dist_rel <- 1- plot_df_combined$dist/ncol(a2_hom)

p3 <- ggplot(plot_df_combined, aes(x=pos, y=dist_rel)) + 
    geom_point() + 
    labs(title='GGCC homozygotes: haplotype extension', x='chr9 position', y='Fraction of samples identical') + 
    geom_vline(xintercept = p6_variant_bed$V2, col='firebrick') +
    geom_vline(xintercept = GGCC_variant_bed$V2, col='orange') +
    geom_vline(xintercept = c(4995000, 5265000), col='blue') + 
    geom_segment(data=jak2_exon_bed, aes(y=-0.1, yend=-0.1, x=start, xend=stop+1000), lwd=2, col='darkgreen') +
    xlim(c(4980000, 5300000))
p3


# # # # # # # # # # # # # # #
# Find this for the GGCC hets #
# # # # # # # # # # # # # # #

# identify for each subject, which allele the GGCC is on
# just look at the 5 snps for now
GGCC_variant_strings <- rownames(vcf_meta)[vcf_meta$pos %in% GGCC_variant_bed$V2]
GGCC_variant_match <- c(1,1,1,1)
names(GGCC_variant_match) <- GGCC_variant_strings

a1_het <- a1[, sample_names_p6_het]
a2_het <- a2[, sample_names_p6_het]
a1_GGCC_inds <- which(apply(a1_het[GGCC_variant_strings, ],2, function(x) identical(x, GGCC_variant_match)))
a2_GGCC_inds <- which(apply(a2_het[GGCC_variant_strings, ],2, function(x) identical(x, GGCC_variant_match)))

# combine these 
combined_het <- cbind(a1_het[, a1_GGCC_inds], a2_het[, a2_GGCC_inds])

# do the same search
# remove rows all zero
min_var <- 2
combined_het2 <- combined_het[rowSums(combined_het) >= min_var, ]
combined_het_dist <- as.matrix(dist((combined_het2), method = 'manh'))

plot_df_het <- data.frame(id=colnames(combined_het_dist), 
                          n=1:ncol(combined_het_dist),
                          pos=vcf_meta[colnames(combined_het_dist), "pos"],
                          dist=combined_het_dist[,target_variant],
                          dist_rel=1-combined_het_dist[,target_variant]/ncol(combined_het))
plot_df_het <- plot_df_het[!duplicated(plot_df_het$pos), ]

p4 <- ggplot(plot_df_het, aes(x=pos, y=dist_rel)) + 
    geom_point() + 
    labs(title='GGCC heterozygotes: haplotype extension', x='chr9 position', y='Fraction of samples identical') + 
    geom_vline(xintercept = p6_variant_bed$V2, col='firebrick') +
    geom_vline(xintercept = GGCC_variant_bed$V2, col='orange') +
    geom_vline(xintercept = c(4995000, 5265000), col='blue') + 
    geom_segment(data=jak2_exon_bed, aes(y=-0.1, yend=-0.1, x=start, xend=stop+1000), lwd=2, col='darkgreen') +
    xlim(c(4980000, 5300000))
p4

# save a combined plot
pdf('~/local_data/p6/GGCC_extent_1KG.pdf', height=6, width = 10)
plot_grid(p3, p4, nrow=2)
dev.off()
