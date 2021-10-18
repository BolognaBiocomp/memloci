'''
Created on 22/giu/09

@author: andrea
'''

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from subprocess import Popen,PIPE
import tempfile,os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
from Bio.Align import AlignInfo
import warnings

class PSSM:
    ''' builds a PSSM from a seqrecord object'''

    def __init__(self,seqrec):
        self.seqrec=seqrec

    def build_pssm_biopy(self,dbfile,include_gap=True,tmpfold='/tmp/', blastpgp_ex='blastpgp', normalize=1, evalue=1e-3, filter_SEG=False,nprocessors=1,nalign=500,back_compatible=0,weight=0,topnumber=5,):
        ''' build_pssm_biopy(dbfile,
        include_gap=True,
        tmpfold='./',
        blastpgp_ex='/usr/bin/blastpgp',
        normalize=1,
        evalue=1e-10,
        filter_SEG=False,
        nprocessors=1,
        nalign=500,
        back_compatible=0
        weight=0#weighted pssm: to be done
        topnumber=5 number of tophits to report
        )'''
        self.top_hits=[]

        tmpfasta=tempfile.NamedTemporaryFile(suffix='.fasta')
        open(tmpfasta.name,'w').write(self.seqrec.format('fasta'))
        cline = NcbiblastpCommandline(query=tmpfasta.name,
                                      db=dbfile,
                                      evalue=evalue,
                                      remote=False,
                                      outfmt=5,
                                      num_threads=nprocessors,
                                      #num_threads=1,
                                      )
        #print(cline)
        trytimes=0
        while trytimes<5:
            blastexec = Popen(str(cline), shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)#not using close_fds=Tru error: blastexec.stdout.read() ends before subprocess ends
            output=blastexec.stdout.read().decode('utf-8')

            blastexec.wait()
            while True:
                retcode=blastexec.poll()
                if retcode==0:
                    trytimes=5
                    break
                elif retcode==None:
                    pass
                else:
                    trytimes+=1
                    break
                    raise warnings.warn('BLAST Error, try=%i, return code= %i'%(trytimes,retcode))


        try:
            blast_results=NCBIXML.parse(StringIO(output))
        except:
            raise Exception('BLAST Error, stderr= %s'%blastexec.stderr.read())
        align=AlignIO.read(StringIO(self.seqrec.format('fasta')),'fasta')#initialize alignment
        for pair in blast_results.__next__().alignments:
            hps=pair.hsps[0]#consider just the best local alignment to include in the pssm
            string_4_pssm=['x']*(hps.query_start-1)
            for i in range(len(hps.query)):
                if hps.query[i]!='-':
                    string_4_pssm.append(hps.sbjct[i])
            string_4_pssm.extend(['x']*(len(self.seqrec.seq)-len(string_4_pssm)))
            string_4_pssm=''.join(string_4_pssm)
            align.append(SeqRecord(Seq(string_4_pssm), id=pair.hit_def))#new in bioython 1.56
            #align.add_sequence(pair.hit_def,string_4_pssm)
            if (len(self.top_hits)<topnumber) and ((hps.identities/float(hps.align_length))>=0.5) and ((hps.identities/float(hps.align_length)))>=0.5:
                self.top_hits.append([pair.hit_def.split('|')[1],'%.1f'%(hps.identities/float(hps.align_length)*100),hps.expect])
        summmary=AlignInfo.SummaryInfo(align)
        char_ign=['x']
        if not include_gap:
            char_ign.append('-')
        self.pssm=summmary.pos_specific_score_matrix(axis_seq=self.seqrec.seq,chars_to_ignore=char_ign)
        if normalize:
            for i in range(len(self.pssm.pssm)):
                tot=sum(self.pssm.pssm[i][1].values())
                for k,v in self.pssm.pssm[i][1].items():
                    if v:
                        if back_compatible:
                            self.pssm.pssm[i][1][k]=float(int(v/tot*100))/100
                        else:
                            self.pssm.pssm[i][1][k]=v/tot
                    else:
                        self.pssm.pssm[i][1][k]=v
        #print(self.pssm)
        return


    def ord_pssm(self,Laa_ord=['V','L','I','M','F','W','Y','G','A','P','S','T','C','H','R','K','Q','E','N','D']):
        '''Reorder the PSSM with a given aminoacidic order, returns the pssm as an indented list
        DEFAULT ORDER:
        ['V','L','I','M','F','W','Y','G','A','P','S','T','C','H','R','K','Q','E','N','D']'''
        self.Lpssm=[]
        for pos in self.pssm:
            L=[]
            for aa in Laa_ord:
                if pos.has_key(aa):
                    L.append(pos[aa])
                else:
                    L.append(0.)
            self.Lpssm.append(L)
        return
