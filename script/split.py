
def buildItemDict(fitem):
    idic = dict()
    with open(fitem,'r') as fobj:
        fobj.readline()         # skip the first line
        csvobj = csv.reader(fobj)
        for line in csvobj:
            key = line[0]
            val1 = sline[1].replace(' ','+')
            val1 = val1.replace('/','+')
            val1 = val1.replace(',','+')
            val = key+'_'+val1+'_'+'_'.join(line[2:])
            idic[key]=val
    return idic
def splitTrain(trainf, idic, odir, field):
    fdic = set()
    with open(trainf,'r') as fobj:
        hline = fobj.readline()
        csvobj = csv.reader(fobj)
        for line in csvobj:
            if not itembn in fdic:
                fdic.add(itembn)
                fout = open(odir+'/'+idic[itembn],'w')
                fout.write(hline)
                fout.close()
            fout = open(odir+'/'+idic[itembn],'a')
            fout.write(line)
            fout.close()
            
