'''
Created on 22/giu/09

@author: andrea
'''
import sys
from . import SeqEncoder, PSSMEncoder
from .PSSM import PSSM


class Protein_encoder:
    '''Takes a SeqRecod Object  and encode it in the requestested schemas '''
    def __init__(self,seqrec, coding_schema,align_db, biopyPSSM):

        self.seqrec=seqrec
        self.CODE=coding_schema
        self.align_db=align_db
        self.encoded_vectors={}
        self.encoded_protein={}
        self._gen_group()
        try:
            self.encode(biopyPSSM)
        except:
            raise
            sys.exit()
        print(self.encoded_protein)
        return


    def _gen_group(self):
        '''precalculate aminoacidic groups to calculate frequency from sequences and profile '''

        D_Lista_gruppi={
        20:[['K'],['R'],['H'],   ['D'],['E'],   ['N'],['Q'],['S'],['T'],['Y'],   ['A'],['V'],['L'],['I'],['P'],['F'],['M'],['W'],['G'],['C']],
        15:[['K','R'],['H'],   ['D'],['E'],   ['N'],['Q'],['S'],['T'],['F','Y'],   ['A'],['V','L','I','M'],['P'],['W'],['G'],['C']],
        13:[['K','R'],['H'],   ['D','E'],   ['N','Q'],['S'],['T'],['F','Y'],   ['A'],['V','L','I','M'],['P'],['W'],['G'],['C']],#aggiunta
        12:[['K','R','H'],   ['D','E'],   ['N','Q'],['S'],['T'],['F','Y'],   ['A'],['V','L','I','M'],['P'],['W'],['G'],['C']],#aggiunta
        10:[['K','R'],['H'],   ['D','E','N','Q'],['S','T'],['F','Y','W'],   ['A'],['V','L','I','M'],['P'],['G'],['C']],
        8:[['K','R'],['H'],   ['D','E','N','Q'],['S','T'],['F','Y','W'],   ['A','G'],['V','L','I','M','C'],['P']],
        4:[['K','R','H',   'D','E','N','Q'],['F','Y','W'],   ['A','G','S','T','P'],['V','L','I','M','C']],
        2:[['K','R','H',   'D','E','N','Q'],['F','Y','W',   'A','G','S','T','P','V','L','I','M','C']]}


        self._precalculated_grups={}#contiene per ogni indice di gruppo i corrispettivi dizionari calcolati solo se richiesti almeno una volta  es: {(20,1):{'A':0,'C':0., ... , Y:0.},...}}
        #==============copy and pasted form old script, rewrite!====#
        for params in self.CODE:
            if params['type']=='AAfreq':
                gruppo=params['aa grouping']
                if gruppo not in D_Lista_gruppi:
                    print ("Errore: gruppo non disponibile")
                    sys.exit
                Lgruppo=D_Lista_gruppi[gruppo]

                # calcolo sottostringhe con gruppi#
                lungh_substring=params['substring']
                if (params['substring'],params['lenght']) not in self._precalculated_grups:#se non e' gia' calcolato
                    Lposizioni=[]#non andare oltre i 5 con codifica a due gruppi  con 512Mb di ram
                    for k in xrange(lungh_substring):#e' la lunghezza delle sottostringhe da cercare
                        La_temp=Lposizioni[:]
                        Lposizioni=[]
                        for a in range(len(Lgruppo)):
                            if len(La_temp)==0:
                                Lposizioni.append([a])
                            else:
                                for b in La_temp:
                                    L=[a]
                                    L.extend(b)
                                    Lposizioni.append(L)
                    if lungh_substring <=3:
                        #calcolo Lpep come tutti i possibili peptidi
                        Lpep=[]
                        for gruppo in Lposizioni:
                            L=['']
                            for i in gruppo:
                                Ltmp=L[:]
                                L=[]
                                for elemento in Ltmp:
                                    for aa in Lgruppo[i]:
                                        L.append(elemento+aa)
                            Lpep.append(L)
                        Dposizioni_pep={}#contiene coppie substringa:indice
                        for k in xrange(len(Lpep)):
                            for stringa in Lpep[k]:
                                Dposizioni_pep[stringa]=0.

                    else:
                        #calcolo Lpep come RE
                        Lpep=[]
                        for gruppo in Lposizioni:
                            search_tot=''
                            for i in gruppo:
                                search='['
                                for aa in Lgruppo[i]:
                                    search+=aa
                                search+=']'
                                search_tot+=search
                            Lpep.append(search_tot)

                        Dposizioni_pep={}#contiene coppie substringa:zero, sara' modificato durante le conte
                        for k in xrange(len(Lpep)):
                                Dposizioni_pep[Lpep[k]]=0.


                    self._precalculated_grups[(params['substring'],params['aa grouping'])]=Dposizioni_pep




    def encode(self, biopyPSSM):

        freq_obj=SeqEncoder.seqfreq(self.seqrec)#create sequence aminoacidic frequency calculation object

        #PSSM_obj=PSSM(self.seqrec)#create PSSM
        #try:
        #    import multiprocessing
        #    CPU = multiprocessing.cpu_count()
        #except ImportError:
        #    CPU=2
        #speed up hack
        #PSSM_obj.build_pssm_biopy(self.align_db,nprocessors=CPU,evalue=1e-5,nalign=250)

        #prof_obj=PSSMEncoder.profreq(PSSM_obj.pssm)#create profile aminoacidic frequency calculation object
        prof_obj=PSSMEncoder.profreq(biopyPSSM)

        for Dcode in self.CODE:
            code=tuple(Dcode.items())
            vector=[]
            if Dcode['type']=='AAfreq':
                freq_obj.calfreq(fin=Dcode['lenght'],substring=Dcode['substring'],Dinput=self._precalculated_grups[(Dcode['substring'],Dcode['aa grouping'])],completetag=code)
                vector=freq_obj.results[code]
            elif Dcode['type']=='Scale':
                freq_obj.calscale(fin=Dcode['lenght'],scale_name=Dcode['scale_name'],completetag=code)
                vector=freq_obj.results[code]
            elif Dcode['type']=='ScaleMobWin':
                freq_obj.calscale_mobwin(fin=Dcode['lenght'],win_len=Dcode['window_lenght'],scale_name=Dcode['scale_name'],select_top=Dcode['select_top'], tresh_mode=Dcode['tresh_mode'],tresh=Dcode['tresh'],completetag=code)
                vector=freq_obj.results[code]
            elif Dcode['type']=='ProteinLenght':
                vector=len(self.seqrec.seq)/float(Dcode['normalize_value'])
            elif Dcode['type']=='Profile':
                if 'conservation_threshold' in Dcode:
                    threshold=Dcode['conservation_threshold']
                else:
                    threshold=0
                prof_obj.calfreq(fin=Dcode['lenght'],completetag=code,threshold=threshold)
                vector=prof_obj.results[code]

            self.encoded_protein[code]=vector#saves results persistently to local database
