'''
Created on 22/giu/09

@author: andrea
'''

import copy
import numpy

class profreq:
    '''Takes as input a protein spssm object, and computes
    the aminoacidic frequencies and parameter scales values for both
    whole sequence and part of it

    call examples:
    freq_obj=profreq(pssm_obj)


    #calfreq
    freq_obj.profreq(fin=0,)
    freq_obj.results['PROfreq_0_1']
       {'A':0.15, 'R':0.03 , 'N':0.002), ...,'W':0 }
    freq_obj.profreq(fin=0,completetag='test')
    freq_obj.results['test']
       {'A':0.15, 'R':0.03 , 'N':0.002), ...,'W':0 }




    '''

    def __init__(self,pssm,Laa_ord=['V','L','I','M','F','W','Y','G','A','P','S','T','C','H','R','K','Q','E','N','D'],):


        self.Laa_ord=Laa_ord
        self.pssm=pssm
        self.ord_pssm()
        self.results={}
        self.scales={}
        self.scales['KD']={'A':  1.800  ,'R': -4.500  ,'N': -3.500  ,'D': -3.500  ,'C':  2.500  ,'Q': -3.500  ,'E': -3.500  ,'G': -0.400  ,'H': -3.200  ,'I':  4.500  ,'L':  3.800  ,'K': -3.900  ,'M':  1.900  ,'F':  2.800  ,'P': -1.600  ,'S': -0.800  ,'T': -0.700  ,'W': -0.900  ,'Y': -1.300  ,'V':  4.200, 'X':0.0, 'Z': -3.500, 'B':-3.500, 'U':2.500}
        self.scales['pol']={'A':  8.100  ,'R': 10.500  ,'N': 11.600   ,'D': 13.000  ,'C':  5.500  ,'Q': 10.500  ,'E': 12.300  ,'G': 9.000  ,'H': 10.400  ,'I':  5.200  ,'L':  4.900  ,'K': 11.300  ,'M':  5.700  ,'F':  5.200  ,'P': 8.000  ,'S': 9.200  ,'T': 8.600  ,'W': 5.400  ,'Y': 6.200  ,'V':  5.900, 'X':8.325, 'Z': 11.400, 'B':12.300, 'U':5.500}
        return

    def ord_pssm(self):
        '''Reorder the PSSM with a given aminoacidic order, returns the pssm as an indented list
        DEFAULT ORDER:
        ['V','L','I','M','F','W','Y','G','A','P','S','T','C','H','R','K','Q','E','N','D']'''
        self.Lpssm=[]
        for pos in self.pssm:
            L=[]
            for aa in self.Laa_ord:
                if aa in pos:
                    L.append(pos[aa])
                else:
                    L.append(0.)
            self.Lpssm.append(L)
        self.Lpssm=numpy.array(self.Lpssm)

        return


    def _profslice(self,position,threshold=0):
        ''' return a sliced sequence string given the cutting position
        slice sintax:
              0 = whole sequence
             20 = 20 aminoterminal residues
            -20 = 20 carboxyterminal residues
        threshold, reports the position in the slice only in the trequency of top conserved residue for the position is
        above the given threshold'''

        Lpssm=copy.deepcopy(self.Lpssm)
        if  position> 0:
            if threshold:
                Lpss_conserved=[]
                for position_profile in Lpssm[:position]:
                    if max(position_profile)>=threshold:
                        Lpss_conserved.append(position_profile)
                Lpssm=numpy.array(Lpss_conserved)
            else:
                Lpssm=Lpssm[:position]
        elif position < 0:
            if threshold:
                Lpss_conserved=[]
                for position_profile in Lpssm[position:]:
                    if max(position_profile)>=threshold:
                        Lpss_conserved.append(position_profile)
                Lpssm=numpy.array(Lpss_conserved)
            else:
                Lpssm=Lpssm[position:]
        else:
            if threshold:
                Lpss_conserved=[]
                for position_profile in Lpssm:
                    if max(position_profile)>=threshold:
                        Lpss_conserved.append(position_profile)
                Lpssm=numpy.array(Lpss_conserved)
            else:
                pass
        return Lpssm



    def calfreq(self,fin=0,substring=1,completetag='',threshold=0):
        '''Compute the frequency of a list of substring in a sequence

        fin= Sequence slice(int):
              0 = whole sequence
             20 = 20 aminoterminal residues
            -20 = 20 carboxyterminal residues
        substring=lenght of the substring to search in the sequence - Just 1 for now  to use numpy
        Dinput= dictinary containing each searched substring as key and a zero int as value

        Computed results are stored in seqfreq.results['PROfreq_'+str(fin)+'_'+str(substring)] as a dictionary'''


        profilo=self._profslice(fin,threshold)
        if not completetag:
            res_string='PROfreq_'+str(fin)+'_'+str(substring)
        else:
            res_string=completetag

        #use numpy to calculate the mean frequency
        if len(profilo):
            self.results[res_string]=dict(zip(self.Laa_ord,profilo.mean(0)))
        else:
            self.results[res_string]=dict(zip(self.Laa_ord,[0.]*len(self.Laa_ord)))

        return
