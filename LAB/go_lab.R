source('BEER.R')

#Load Demo Data (Oligodendroglioma, GSE70630)
#Download: https://sourceforge.net/projects/beergithub/files/

D1=readRDS('MGH36.RDS')
D2=readRDS('MGH53.RDS')
D3=readRDS('MGH54.RDS')
D4=readRDS('MGH60.RDS')
D5=readRDS('MGH93.RDS')
D6=readRDS('MGH97.RDS')

BATCH=c(rep('D1',ncol(D1)),
        rep('D2',ncol(D2)),
        rep('D3',ncol(D3)),
        rep('D4',ncol(D4)),
        rep('D5',ncol(D5)),
        rep('D6',ncol(D6)) )
        
D12=.simple_combine(D1,D2)$combine
D34=.simple_combine(D3,D4)$combine
D56=.simple_combine(D5,D6)$combine
D1234=.simple_combine(D12,D34)$combine
D123456=.simple_combine(D1234,D56)$combine

DATA=D123456   

rm(D1);rm(D2);rm(D3);rm(D4);rm(D5);rm(D6)
rm(D12);rm(D34);rm(D56);rm(D1234);rm(D123456)





