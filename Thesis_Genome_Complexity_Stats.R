setwd('/Users/alex/Thesis_Writing/Chapter5')
dex <- read.table('Genome_Complexity_Calculation_R.txt',sep='\t',header=T,row.names=1)

t.test(dex$Genome.Percentage[dex$Status == 'Biallelic'],dex$Genome.Percentage[dex$Status == 'Normal'])
t.test(dex$Genome.Percentage[dex$Status == 'Monoallelic'],dex$Genome.Percentage[dex$Status == 'Normal'])

t.test(dex$CN.Segments[dex$Status == 'Biallelic'],dex$CN.Segments[dex$Status == 'Normal'])
t.test(dex$CN.Segments[dex$Status == 'Monoallelic'],dex$CN.Segments[dex$Status == 'Normal'])

t.test(dex$Genome.Percentage[dex$TP53_Any == 'Yes'],dex$Genome.Percentage[dex$TP53_Any == 'No'])

t.test(dex$CN.Segments[dex$TP53_Any == 'Yes'],dex$CN.Segments[dex$TP53_Any == 'No'])

mean(dex$Genome.Percentage[dex$TP53_Any == 'Yes'])
mean(dex$Genome.Percentage[dex$TP53_Any == 'No'])
mean(dex$Genome.Percentage[dex$Status == 'Biallelic'])
mean(dex$Genome.Percentage[dex$Status == 'Monoallelic'])

mean(dex$CN.Segments[dex$TP53_Any == 'Yes'])
mean(dex$CN.Segments[dex$TP53_Any == 'No'])
mean(dex$CN.Segments[dex$Status == 'Biallelic'])
mean(dex$CN.Segments[dex$Status == 'Monoallelic'])

plex <- read.table('Ch5_Complexity_Summary.txt',sep='\t',header=T,row.names=1)

t.test(plex$N_Complex_abns[plex$Status == 'Biallelic'],plex$N_Complex_abns[plex$Status == 'No'])
t.test(plex$N_Complex_abns[plex$Status == 'Monoallelic'],plex$N_Complex_abns[plex$Status == 'No'])

t.test(plex$N_Complex_abns[plex$TP53_Any == 'Yes'], plex$N_Complex_abns[plex$TP53_Any == 'No'])

mean(plex$N_Complex_abns[plex$TP53_Any == 'Yes'])
mean(plex$N_Complex_abns[plex$TP53_Any == 'No'])
mean(plex$N_Complex_abns[plex$Status == 'Biallelic'])
mean(plex$N_Complex_abns[plex$Status == 'Monoallelic'])
