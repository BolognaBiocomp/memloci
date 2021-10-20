#!/usr/bin/env python

'''
Created on 22/giu/09

@author: andrea
'''

import os
MEMLOCI_HOME=os.environ["MEMLOCI_HOME"]
import sys,types,shelve,os,datetime
import pickle as cPickle
import argparse, json
sys.path.append(MEMLOCI_HOME)
from Bio import SeqIO
from encode.ProteinEncoder import Protein_encoder
from SVM.SVMlike import *

from memlocilib import config
from memlocilib import utils
from memlocilib import workenv
from memlocilib import blast
from memlocilib import cpparser

class Parameters:
    def __init__(self):
        self.ALLSchemas=[   {'type':'ProteinLenght','normalize_value':2000,},
                            {'type':'Scale','lenght':40,'scale_name':'KD'},
                            #{'type':'Scale','lenght':40,'scale_name':'pol'},
                            {'type':'Scale','lenght':60,'scale_name':'KD'},
                            {'type':'Scale','lenght':-60,'scale_name':'KD'},
                            #{'type':'Scale','lenght':60,'scale_name':'pol'},
                            #{'type':'Scale','lenght':-60,'scale_name':'pol'},
                            {'type':'ScaleMobWin','lenght':0,'window_lenght':6,'scale_name':'KD','select_top':True,'tresh_mode':True,'tresh': 0.33},
                            {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':0},
                            {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':0,'conservation_threshold':0.75},
                            {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':40,'conservation_threshold':0.3},
                            {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':-40,'conservation_threshold':0.3},
                            {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':60,'conservation_threshold':0.3},
                            {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':60},
                            {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':-60,'conservation_threshold':0.3},
                            {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':40},
                            ]

        self.EncodingSchemas={'CM':[{'type':'Scale','lenght':40,'scale_name':'KD'},
                                    #{'type':'Scale','lenght':40,'scale_name':'pol'},
                                    {'type':'ProteinLenght','normalize_value':2000,},
                                    {'type':'ScaleMobWin','lenght':0,'window_lenght':6,'scale_name':'KD','select_top':True,'tresh_mode':True,'tresh': 0.33},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':0},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':0,'conservation_threshold':0.75},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':40,'conservation_threshold':0.3},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':-40,'conservation_threshold':0.3},
                                    ],
                             'ORG':[{'type':'ProteinLenght','normalize_value':2000,},
                                    {'type':'Scale','lenght':60,'scale_name':'KD'},
                                    {'type':'Scale','lenght':-60,'scale_name':'KD'},
                                    #{'type':'Scale','lenght':60,'scale_name':'pol'},
                                    #{'type':'Scale','lenght':-60,'scale_name':'pol'},
                                    {'type':'ScaleMobWin','lenght':0,'window_lenght':6,'scale_name':'KD','select_top':True,'tresh_mode':True,'tresh': 0.33},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':0},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':0,'conservation_threshold':0.75},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':60,'conservation_threshold':0.3},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':60},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':-60,'conservation_threshold':0.3},
                                    ],
                            'ENDO':[{'type':'ProteinLenght','normalize_value':2000,},
                                    {'type':'Scale','lenght':40,'scale_name':'KD'},
                                    #{'type':'Scale','lenght':40,'scale_name':'pol'},
                                    {'type':'ScaleMobWin','lenght':0,'window_lenght':6,'scale_name':'KD','select_top':True,'tresh_mode':True,'tresh': 0.33},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':0},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':0,'conservation_threshold':0.75},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':40},
                                    {'type':'Profile','profiletag':'PSSM_uniprot_1e-5','lenght':40,'conservation_threshold':0.3},
                                    ],
                             }
        self.SVM_model={'CM':'-t 2 -g 3 -c 6',
                    'ORG':'-t 2 -g 3 -c 6',
                    'ENDO':'-t 2 -g 3 -c 6',}


def Train(fasta_file,db_profile,IDs):
    '''Train the predictor with the new dataset '''

    from SVM.SVMtt import SVM_Modelizer,Post_SVM_learning
    params=Parameters()
    #Persist=shelve.open('training_data.shv',writeback=True)
    Persist={}
    Persist['EncodedSequences']={}


    inizio=datetime.datetime.now()
    def draw_advance(done,tot,wide):
        '''draws advance bar'''
        advance=float(done)/tot
        #percentage=fpformat.fix(advance*100,2)+'%'
        percentage = str(round(100*advance,2))+'%'
        position= int(advance*wide)
        before=position
        after=wide-position-1
        bar='['+'='*before+'>'+'-'*after+']'
        result=bar+' '+percentage
        if done==tot+1:
            result+=' DONE!'
        return result

    def Test(D_protein_vector,label):
        '''Test provided feature vector'''
        modpath='./SVM/'
        Dprediction={}
        svm=getSVMLight(modpath+label+'_MOD')
        for ID in D_protein_vector:
            if not Dprediction.has_key(ID):
                Dprediction[ID]={}
            Dprediction[ID]=svm.predict(D_protein_vector[ID])
        print (Dprediction[Dprediction.keys()[0]],len(Dprediction))
        return Dprediction


    '''Encode sequences '''
    c=0.
    print ("Encoding in progress:")
    totprot=len(list(SeqIO.parse(open(fasta_file), 'fasta')))
    for seqrec in SeqIO.parse(open(fasta_file), 'fasta'):
        try:
            if '|' in seqrec.id:
                seqrec.id=seqrec.id.split('|')[1]
            if '_' in seqrec.id:
                seqrec.id=seqrec.id.split('_')[1]
            c+=1
            try:
                ora=datetime.datetime.now()
                sys.stdout.write('\r'+draw_advance(c,totprot,50)+' ELAPSED TIME: '+str(ora-inizio).split('.')[0]+' REMAINING TIME: '+str(((ora-inizio)*(int((totprot/c)*1000))/1000)-(ora-inizio)).split('.')[0])
                sys.stdout.flush()
            except ZeroDivisionError:
                pass
            last_encoding = Protein_encoder(seqrec,params.ALLSchemas,db_profile).encoded_protein
            Persist['EncodedSequences'][seqrec.id]=last_encoding
        except:
            Persist['EncodedSequences'][seqrec.id]=last_encoding

    cPickle.dump(Persist, open('training_data.pkl','w',))
    #Persist.sync()
    '''build SVM models'''
    if not os.path.exists('./SVM/'):
        os.mkdir('./SVM/')
    D_Tested_proteins={}
    D_performance_gradients={}
    D_balancing_threshold={}
    for label in IDs:
        #create feature input  vectors
        Dcompletevectors={}
        Dlabel_input={}
        LID_pos=[]
        LID_neg=[]
        for class_label in IDs:
            for ID in IDs[class_label]:
                if ID in Persist['EncodedSequences'].keys():
                    if class_label == label:
                        Dlabel_input[ID]=1
                        LID_pos.append(ID)
                    else:
                        Dlabel_input[ID]=-1
                        LID_neg.append(ID)
                    protein_vector=[]
                    for CODE in params.EncodingSchemas[label]:
                        code=tuple(CODE.items())
                        value=Persist['EncodedSequences'][ID][code]
                        if type(value)== types.DictType:
                            protein_vector.extend(value.values())
                        elif (type(value)== types.TupleType) or (type(value)== types.ListType):
                            protein_vector.extend(value)
                        elif type(value)== types.FloatType:
                            protein_vector.append(value)
                        elif type(value)== types.IntType:
                            protein_vector.append(float(value))
                    Dcompletevectors[ID]=protein_vector
        SVM_Modelizer(Dlabel_input,Dcompletevectors,learning_parameters=params.SVM_model[label], path='./SVM/',modelname=label+'_MOD',)
        D_Tested_proteins[label]=Test(Dcompletevectors,label)#test training dataset

        D_balancing_threshold[label],D_performance_gradients[label]=Post_SVM_learning().balance(LID_pos,LID_neg,D_Tested_proteins[label])
    cPickle.dump(D_balancing_threshold, open('./SVM/thresholds.pk','w'),)
    cPickle.dump(D_performance_gradients, open('./SVM/gradients.pk','w'),)
    return



def Test(seqrec, biopyPSSM, align_db):
    '''Test provided dataset'''
    params=Parameters()
    EncodedSequence=Protein_encoder(seqrec,params.ALLSchemas,align_db, biopyPSSM).encoded_protein
    modpath=os.path.join(MEMLOCI_HOME, 'models')+'/'
    Dmodels={'CM':modpath+'CM_MOD','ORG':modpath+'ORG_MOD','ENDO':modpath+'ENDO_MOD',}
    Dprediction={}
    Dprediction_score={}
    DGradients=cPickle.load(open(modpath+'gradients.pk','rb'), encoding="latin-1")
    #print(DGradients.keys())
    for label in params.EncodingSchemas:
        protein_vector=[]
        for CODE in params.EncodingSchemas[label]:
            code=tuple(CODE.items())
            value=EncodedSequence[code]
            if isinstance(value, dict):
                protein_vector.extend(value.values())
            elif isinstance(value, tuple) or isinstance(value, list):
                protein_vector.extend(value)
            elif isinstance(value, float):
                protein_vector.append(value)
            elif isinstance(value, int):
                protein_vector.append(float(value))

        #predict
        svm=getSVMLight(Dmodels[label])
        Dprediction[label]=svm.predict(protein_vector)
        performance=1.
        points=list(DGradients[label].keys())
        points.sort()
        points.reverse()
        for point in points:
            if Dprediction[label] >=point:
                #performance=fpformat.fix((DGradients[label][point]['Cov'][0]-0.5)*200,0)+'%'
                performance = str(round((DGradients[label][point]['Cov'][0]-0.5)*200,0))+'%'
                break
        #performance='0%'
        Dprediction_score[label]=performance

    Lprediction=[]
    for k,v in Dprediction.items():
        Lprediction.append((v,k))
    Lprediction.sort()
    Lprediction.reverse()
    #print "\t".join(map(str, [seqrec.id,Lprediction[0][1],Dprediction_score,Dprediction]))
    return seqrec.id,Lprediction[0][1],Dprediction_score,Dprediction

def main():
    DESC = "MemLoci: Prediction of protein membrane localization"
    parser = argparse.ArgumentParser(description = DESC, prog = "memloci.py")
    parser.add_argument("-f", "--fasta", help = "The input FASTA file name", dest = "fasta", required = True)
    parser.add_argument("-d", "--dbfile",
                              help = "The PSIBLAST DB file",
                              dest = "dbfile", required= True)
    parser.add_argument("-o", "--outf", help = "The output file", dest = "outf", required = True)
    parser.add_argument("-m", "--outfmt", help = "The output format: json or gff3 (default)", choices=['json', 'gff3'], required = False, default = "gff3")
    parser.add_argument("-c", "--cache-dir", help="Cache dir for alignemnts", dest="cache_dir", required=False, default=None)
    parser.add_argument("-j", "--psiblast-iter", help="Number of PSIBLAST iterations (default 3)", dest="pbniter", required=False, default=3, type=int)
    parser.add_argument("-n", "--psiblast-nalign", help="PSIBLAST num_alignments parameter (default 5000)", dest="pbnalign", required=False, default=5000, type=int)
    parser.add_argument("-e", "--psiblast-evalue", help="PSIBLAST evalue parameter (default 0.001)", dest="pbeval", required=False, default=0.001, type=float)
    ns = parser.parse_args()

    workEnv = workenv.TemporaryEnv()
    data_cache = utils.get_data_cache(ns.cache_dir)
    i=0
    protein_jsons = []
    ofs=open(ns.outf, 'w')
    for record in SeqIO.parse(ns.fasta, 'fasta'):
        seqid = record.id
        seq = str(record.seq)
        prefix="seq%d"%i
        fastaSeq  = workEnv.createFile(prefix+".", ".fasta")
        fsofs=open(fastaSeq,'w')
        print(">%s" % seqid, file=fsofs)
        print(seq, file=fsofs)
        fsofs.close()
        pssmFile = blast.runPsiBlast(prefix, ns.dbfile, fastaSeq, workEnv, data_cache=data_cache,
                                     num_alignments=ns.pbnalign, num_iterations=ns.pbniter, evalue=ns.pbeval)
        profile_matrix = cpparser.BlastCheckPointProfile(pssmFile)
        seqrec = SeqIO.read(open(fastaSeq),'fasta')
        biopyPSSM = utils.get_biopy_pssm(str(seqrec), profile_matrix)
        memloci_pred = Test(seqrec, biopyPSSM, ns.dbfile)
        if ns.outfmt == "gff3":
            loc = memloci_pred[1]
            score = float(memloci_pred[2][loc][:-1])/100.0
            utils.write_gff_output(seqid, seq, ofs, loc, score)
        else:
            i_json = {'accession': seqid, 'comments': [], "dbReferences": []}
            i_json['sequence'] = {
                                "length": len(seq),
                                "sequence": seq
                             }
            acc_json = utils.get_json_output(i_json, memloci_pred)
            protein_jsons.append(acc_json)
        i = i + 1
    if ns.outfmt == "json":
        json.dump(protein_jsons, ofs, indent=5)
    ofs.close()
    workEnv.destroy()

    """
    ifs = open(ns.i_json)
    input_json = json.load(ifs)
    ifs.close()
    protein_jsons = []
    i=0
    for i_json in input_json:
        seqid = i_json['accession']
        #seq = i_json['sequence']['sequence']
        seq, cleavage = utils.cut_peptide(i_json)

        prefix="seq%d"%i
        fastaSeq  = workEnv.createFile(prefix+".", ".fasta")
        fsofs=open(fastaSeq,'w')
        #SeqIO.write([fasta], fsofs, 'fasta')
        print(">%s" % seqid, file=fsofs)
        print(seq, file=fsofs)
        fsofs.close()
        seqrec = SeqIO.read(open(fastaSeq),'fasta')
        memloci_pred = Test(seqrec, ns.dbfile)
        #print(memloci_pred)
        acc_json = utils.get_json_output(i_json, memloci_pred)
        protein_jsons.append(acc_json)
        i = i + 1
    ofs=open(ns.outf, 'w')
    json.dump(protein_jsons, ofs, indent=5)
    ofs.close()
    workEnv.destroy()
    """

if __name__ == '__main__':
    """
    if len(sys.argv)==1:
        USAGE='''
        TRAIN: memloci.py -train   Traindataset.fasta  IDfile:CM  IDfile:ORG  IDfile:ENDO   database_to_build_profile
        TEST:  memloci.py file.fasta database_to_build_profile
        '''
        print USAGE
    else:
        if sys.argv[1]=='-train':#train
            #try:
                DID={}
                for filename in sys.argv[3:-1]:
                    DID[filename.split(':')[1]]=open(filename.split(':')[0]).read().split()
                Train(fasta_file=sys.argv[2],db_profile=sys.argv[-1], IDs=DID)
                print '''Memloci was trained correctly with new dataset'''
            #except Exception, error:
            #    print '''Training procedure went wrong... ''',error
        else:# test
            for seqrec in SeqIO.parse(open(sys.argv[1]), 'fasta' ):
                Test(seqrec, sys.argv[-1])
    """
    main()
#                if seqrec.id.split('_')[1] in '''P32790 ['Cell_memb', 'ER']
#Q9H4E7 ['Cell_memb', 'ER']
#P41318 ['Cell_memb', 'ER']
#Q00583 ['ER', 'Mito']
#Q60766 ['Cell_memb', 'ER']
#Q15642 ['Cell_memb', 'ER']
#Q96RU3 ['Cell_memb', 'ER']
#P23377 ['Cell_memb', 'ER']
#Q92797 ['Cell_memb', 'ER']
#Q9NY26 ['Cell_memb', 'ER']
#P20398 ['Cell_memb', 'ER']
#P49597 ['Cell_memb', 'ER']
#P25618 ['Cell_memb', 'ER']
#Q969V5 ['ER', 'Mito']
#Q924S8 ['Cell_memb', 'ER']
#Q9EQD0 ['Cell_memb', 'ER']
#P08011 ['ER', 'Mito']
#Q99JE4 ['Cell_memb', 'ER']
#O00750 ['Cell_memb', 'ER']
#Q93052 ['Cell_memb', 'ER']
#O95996 ['Cell_memb', 'ER']
#Q96IF1 ['Cell_memb', 'ER']
#Q3ZCQ8 ['ER', 'Mito']
#''':
#
#
#                    Test(seqrec, sys.argv[-1])
