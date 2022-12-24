import pysam
import re
import numpy as np

samfile =  pysam.AlignmentFile("/media/wentao/d3294dee-a7f0-442b-a55c-1399de1c7729/PB000352/test.sam", "rb")
for read in samfile.fetch(until_eof = True):
    try:
        mmtag = read.get_tag("MZ")
        if(re.search('C\+m,\d+', mmtag)):
        #mltag = mltag.rstrip(';')
        #mltag = mltag.split(',')[1:]
            mmtag = mmtag.rstrip(';')
            mmtag = mmtag.split(',')[1:]
        #print(mltag)
            mmtag = np.array([int(x)+1 for x in mmtag])
            mmtagpos = mmtag.cumsum()-1
            mltag = [(int(i)+1)/256  for i in read.get_tag("ML")]
        #print(mmtagpos-1)
        #print(mmtagpos)
        seq = read.get_forward_sequence()
        cindex = re.finditer(pattern='C', string=seq)
        indices = np.array([index.start()+1 for index in cindex])
        relpos = indices[mmtagpos]
    #print(len(relpos))
        for i in range(0,len(relpos)):
            print(read.query_name,relpos[i], mltag[i], sep=",")
    except:
        print('not find')
    #rint(dir(read))

    
    