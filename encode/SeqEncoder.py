'''
Created on 22/giu/09

@author: andrea
'''
import re,copy

class seqfreq:
    '''Takes as input a protein seqrecord object, and computes
    the aminoacidic frequencies and parameter scales values for both
    whole sequence and part of it

    call examples:
    freq_obj=seqfreq(seqrec_obj)
    freq_obj.calculate_all()#calculate standard set of calculations
    print freq_obj.results # see the results

    #calfreq
    freq_obj.calfreq(fin=0,substring=1,Dfreq={'A':0, 'R':0 , 'N':0), ...,'W':0 },regex=0)
    freq_obj.results['AAfreq_0_1']
       {'A':0.15, 'R':0.03 , 'N':0.002), ...,'W':0 }
    freq_obj.calfreq(substring=2,Dfreq={'AA':0, 'AR':0 , 'AN':0), ...,'WW':0 })
    freq_obj.results['AAfreq_0_2']
       {'AA':0.015, 'AR':0.003 , 'AN':0.0002), ...,'WW':0 }

    #calscale
    freq_obj.calscale(20,'KD')#calculate the kite-dolittle scale for the first 20 residues
    freq_obj.results['KDscale_20']
       0.236780
    freq_obj.calscale(0,scale_name='My', scale={'A': 0.3, 'C': 1, ..., 'W': 0.})
    freq_obj.results['Myscale_0']
       0.236780

     #calscale_mobwin
     freq_obj.calscale_mobwin(fin=0,win_len=9,scale_name='KD',select_top=True)
         #calculate the highest score value for the KD scale on a window of lenght 19 of the
         #whole sequence.
     freq_obj.results['KDscale_mobwin_0_9']
         0.7832443
     freq_obj.calscale_mobwin(fin=0,win_len=9,scale_name='KD',select_top=False)
         #report the lower value
         -0.634203
     freq_obj.calscale_mobwin(fin=0,win_len=9,scale_name='KD',tresh_mode=True,thresh=0.8)
         #reports all the position in which the mobile windows scale was above the treshold
         {32:0.82431212,67:0.8677323,167:0.977781231}



    '''

    def __init__(self,seqrec):
        self.seq=str(seqrec.seq)
        self.results={}
        self.scales={}
        self.scales['KD']={'A':  1.800  ,'R': -4.500  ,'N': -3.500  ,'D': -3.500  ,'C':  2.500  ,'Q': -3.500  ,'E': -3.500  ,'G': -0.400  ,'H': -3.200  ,'I':  4.500  ,'L':  3.800  ,'K': -3.900  ,'M':  1.900  ,'F':  2.800  ,'P': -1.600  ,'S': -0.800  ,'T': -0.700  ,'W': -0.900  ,'Y': -1.300  ,'V':  4.200, 'X':0.0, 'Z': -3.500, 'B':-3.500, 'U':2.500}
        self.scales['pol']={'A':  8.100  ,'R': 10.500  ,'N': 11.600   ,'D': 13.000  ,'C':  5.500  ,'Q': 10.500  ,'E': 12.300  ,'G': 9.000  ,'H': 10.400  ,'I':  5.200  ,'L':  4.900  ,'K': 11.300  ,'M':  5.700  ,'F':  5.200  ,'P': 8.000  ,'S': 9.200  ,'T': 8.600  ,'W': 5.400  ,'Y': 6.200  ,'V':  5.900, 'X':8.325, 'Z': 11.400, 'B':12.300, 'U':5.500}
        return

    def calcall(self):
        '''calculate everything. To be implemented'''
        self.calfreq(0, 1, Dinput={}, regexp=0)
        self.calfreq(20, 1, Dinput={}, regexp=0)
        self.calfreq(40, 1, Dinput={}, regexp=0)
        self.calfreq(60, 1, Dinput={}, regexp=0)
        self.calfreq(-20, 1, Dinput={}, regexp=0)
        self.calfreq(-40, 1, Dinput={}, regexp=0)
        self.calfreq(-60, 1, Dinput={}, regexp=0)
        self.calfreq(0, 2, Dinput={}, regexp=0)
        self.calfreq(20, 2, Dinput={}, regexp=0)
        self.calfreq(50, 2, Dinput={}, regexp=0)
        self.calfreq(-20, 2, Dinput={}, regexp=0)
        self.calfreq(-50, 2, Dinput={}, regexp=0)
        self.calscale(0,'KD')
        self.calscale(20,'KD')
        self.calscale(40,'KD')
        self.calscale(-20,'KD')
        self.calscale(-40,'KD')
        self.calscale(0,'pol')
        self.calscale(20,'pol')
        self.calscale(40,'pol')
        self.calscale(-20,'pol')
        self.calscale(-40,'pol')
        self.calscale_mobwin(0, 9, 'KD',)
        self.calscale_mobwin(0, 9, 'pol',)
        self.calscale_mobwin(0, 9, 'KD',tresh_mode=True,tresh=0.5)
        self.calscale_mobwin(0, 9, 'KD',select_top=False,tresh_mode=True,tresh=-0.5)
        self.calscale_mobwin(0, 9, 'pol',tresh_mode=True,tresh=0.5)
        self.calscale_mobwin(0, 9, 'pol',select_top=False,tresh_mode=True,tresh=-0.5)




        return

    def _seqslice(self,position):
        ''' return a sliced sequence string given the cutting position
        slice sintax:
              0 = whole sequence
             20 = 20 aminoterminal residues
            -20 = 20 carboxyterminal residues'''
        seq=copy.deepcopy(self.seq)
        if  position> 0:
            seq=seq[:position]
        elif position < 0:
            seq=seq[position:]
        else:
            pass
        return seq



    def calfreq(self,fin=0,substring=1,Dinput={},regexp=0,completetag=''):
        '''Compute the frequency of a list of substring in a sequence

        fin= Sequence slice(int):
              0 = whole sequence
             20 = 20 aminoterminal residues
            -20 = 20 carboxyterminal residues
        substring=lenght of the substring to search in the sequence
        Dinput= dictinary containing each searched substring as key and a zero int as value
        regexp= force the use of RE module to search substrings

        Computed results are stored in seqfreq.results['AAfreq_'+str(fin)+'_'+str(substring)] as a dictionary'''

        if Dinput=={}:#Develop mode: please provide Dinput outside of the funtion for speed increase
            #BUG:id a dictionary is created without declaring the empy dic in the function
            #input it remains for furter iterations. ALWAYS provide at least an empy dic.
            Laa=['K','R','H',   'D','E','N','Q','F','Y','W',   'A','G','S','T','P','V','L','I','M','C']
            Lpep=['']
            for s in range(substring):
                Ltmp=Lpep[:]
                Lpep=[]
                for a in Laa:
                    for b in Ltmp:
                        Lpep.append(a+b)
            for key in Lpep:
                Dinput[key]=0

        sequenza=self._seqslice(fin)
        if not completetag:
            res_string='AAfreq_'+str(fin)+'_'+str(substring)
        else:
            res_string=completetag
        self.results[res_string]=copy.deepcopy(Dinput)


        incremento=1.0/len(sequenza)
        if (not regexp) and (substring>3):
            for i in xrange(len(sequenza)+1-substring):
                try:
                    self.results[res_string][sequenza[i:i+substring]]+=incremento
                except KeyError:
                    continue
        else:#Always use RE for substring>3 residues
            for key in Dinput.keys():
                regola=re.compile(key)
                self.results[res_string][key]+=(incremento*len(regola.findall(sequenza)))
        return


    def calscale(self,fin=0,scale_name='KD',scale={},completetag=''):
        '''Compute Protein scales values. default is Kite-Dolittle scale for whole seq
        scale values are normalized from -1 to 1 relating on the maximum and minimum values
        of the given scale applied to the analyzed sequence.
        Please note: scales need to be simmetric to fill the -1<-->1 range

        Scales:
        scale_name='KD'    -->Kite-Dolittle, Idrophobicity
        scale_name='Pol'    -->Grantham, polarity
        scale_name='myscale', scale={'A': 0.3, 'C': 1, ..., 'W': 0.} --> custom scale
        '''

        if scale_name=='KD':
            scale=self.scales['KD']
        elif scale_name=='pol':
            scale=self.scales['pol']
        res_string=scale_name+'scale_'+str(fin)


        sequence=self._seqslice(fin)

        maxvalue=float(len(sequence))*max(scale.values())
        minvalue=float(len(sequence))*min(scale.values())
        if abs(minvalue) >maxvalue:
            normalizer=abs(minvalue)
        else:
            normalizer=maxvalue

        scale_value=0.
        for residue in sequence:
            scale_value+=scale[residue]

        scale_value=scale_value/normalizer

        if completetag:
            res_string=completetag

        self.results[res_string]=scale_value
        return


    def calscale_mobwin(self,fin=0,win_len=9,scale_name='KD',scale={},select_top=True, tresh_mode=False,tresh=-1,completetag=''):
        '''compute the maximum value of a windows oftresh_mode=True,tresh=0tresh_mode=True,tresh=0 lenght win_len in the sequence
        for scales doc see calscale
        select_top=True returns the maximum value
        select_top=False returns the minimum value
        win_len=9 calculates for 19 residue long mobile window
        tresh_mode=True,tresh=0.56 ---> the number of  mobile windows above a given trheshold
        '''

        if scale_name=='KD':
            scale=self.scales['KD']
        elif scale_name=='pol':
            scale=self.scales['pol']
        res_string=scale_name+'scale_mobwin_'+str(fin)+'_'+str(win_len)
        if tresh_mode:
            res_string+='_tresh'
        if not select_top:
            res_string+='_seldown'

        sequence=self._seqslice(fin)

        maxvalue=float(win_len*2+1)*max(scale.values())
        minvalue=float(win_len*2+1)*min(scale.values())
        if abs(minvalue) >maxvalue:
            normalizer=abs(minvalue)
        else:
            normalizer=maxvalue

        Dres={}
        k=0
        mobwin=sequence[k:k+2*win_len+1]
        tot=0.0
        for A in mobwin:
            tot+=scale[A]
        Dres[k+win_len]=tot/normalizer

        while k <(len(sequence)-(2*win_len)-1):
            k+=1
            tot=tot-scale[sequence[k-1]]+scale[sequence[k+(2*win_len)]]
            Dres[k+win_len]=tot/normalizer

        if tresh_mode:
            if (select_top and tresh<=-1) or ((not select_top) and tresh>=1):#no threshold set, return all
                if completetag:
                    res_string=completetag
                self.results[res_string]=len(Dres)
            else:#calculate fraction of residues with surrounding enviroment beyond threshold
                Dtresh={}
                for k,v in Dres.items():
                    if select_top:
                        if v>=tresh:
                            Dtresh[k]=v
                    else:
                        if v<=tresh:
                            Dtresh[k]=v
                if completetag:
                    res_string=completetag
                self.results[res_string]=len(Dtresh)/float(len(sequence))

        else:
            if select_top:
                bestvalue=max(Dres.values())
            else:
                bestvalue=min(Dres.values())
            if completetag:
                res_string=completetag
            self.results[res_string]=bestvalue


        return
