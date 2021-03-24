'''
Created on 23/giu/09

@author: andrea
'''
import os,sys
from subprocess import Popen,PIPE

class SVM_Modelizer:
    
    def __init__(self,data_label,data_vector,learning_parameters='-t 2 -g 3 -c 100',modelname='svm_model',path='./',write_input=True,feature_filter=False,learn=True):
        '''
        data_label: dictionary containing  for each protein as key a name and the corrisponding class label (1,0,-1,...etc)
                    must contain only the protein to be included in the model
        data_vector: dictionary containing  for each protein as key a name and the corrisponding feature vector as value
                    can contain also vetors non included in the model
        path: the path to write the input file and the model file
        modelname: the name of the file containing the generated model
        learning_parameters: svm parameter to be passed to svm_light
        '''
        self.labels=data_label
        self.vectors=data_vector
        self.path=path
        if not os.path.exists(path):
            os.mkdir(path)
        self.modname=modelname
        self.svm_param=learning_parameters
        self.input_file=self.path+'svm_in_learn_'+self.modname
        
        if feature_filter:
            self.ff_flag=True
            self.feature_filter=feature_filter
        else:
            self.ff_flag=False
        
        if write_input:
            self.write_input()
        if learn:
            self.learn()
        
    def _encode_vector(self,ID):
        encoded_vector=[str(self.labels[ID])]
        if self.ff_flag:
            if not self.vectors.has_key(ID):
                print 'WARNING:',ID,'''not present in encoded vector for modeling, total encoded vectors:''',len(self.vectors)
                sys.exit()
            else:
                for k,v in zip(xrange(len(self.vectors[ID])),self.vectors[ID]):
                    if k in self.feature_filter:
                        if v !=0.:#lower disk space usage
                            encoded_vector.append(str(k+1)+':'+str(v))
        else:
            if not self.vectors.has_key(ID):
                print 'WARNING:',ID,'''not present in encoded vector for modeling, total encoded vectors:''',len(self.vectors)
                sys.exit()
            else:            
                for k,v in zip(xrange(len(self.vectors[ID])),self.vectors[ID]):
                    if v !=0.:#lower disk space usage
                        encoded_vector.append(str(k+1)+':'+str(v))
        encoded_vector.append('#'+ID)
        encoded_vector=' '.join(encoded_vector)
        return encoded_vector
        
    
    def write_input(self):
        
        out=open(self.input_file,'w')
        for ID in self.labels.keys():
            out.write(self._encode_vector(ID)+'\n')
        out.close()
        
    def learn(self):
        cmd=' svm_learn '+self.svm_param+' '+self.input_file+' '+self.path+self.modname
        trainer = Popen(str(cmd), shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)#not using close_fds=Tru error:.read() ends before subprocess ends
        output=trainer.stdout.read()
        trainer.wait()
    


class SVM_Tester:
    '''takes in input a data vector dict and returns a dict with predicted value
    it performs at default:
         self.write_input    write the input vector to disk
         self.test            launch svm_classify to predict given a model
         self.retrieve        read predictions and returns them as a dict
         
    self.feature_filter is a list of allowed features to be written in the svm input
        alternative faster implementation with numpy, but little to gain, do it just for production
        
        self.Dresults contains the retrieved results'''
         
    
    def __init__(self,data_vector,modelname='svm_model',path='./',feature_filter=False,test_label='0',input_file='',output_file='',write_input=True,test=True,retrieve=True,tag=''):
        self.vectors=data_vector
        self.path=path
        self.modname=modelname
        self.ID_ordered_list=self.vectors.keys()
        self.label=test_label
        if feature_filter:
            self.ff_flag=True
            self.feature_filter=feature_filter
        else:
            self.ff_flag=False
        if input_file:
            self.input_file=input_file
        else:
            self.input_file=self.path+'svm_in_test_'+tag+self.modname
        if output_file:
            self.output_file=output_file
        else:
            self.output_file=self.path+'svm_in_test_'+tag+self.modname+'_predicted'
        if write_input:
            self.write_input()
        if test:
            self.test()
        if retrieve:
            self.retrieve()
            return 
        
    def _encode_vector(self,ID):
        encoded_vector=[self.label]
        if self.ff_flag:
            for k,v in zip(xrange(len(self.vectors[ID])),self.vectors[ID]):
                if k in self.feature_filter:
                    if v !=0.:#lower disk space usage
                        encoded_vector.append(str(k+1)+':'+str(v))
        else:
            for k,v in zip(xrange(len(self.vectors[ID])),self.vectors[ID]):
                if v !=0.:#lower disk space usage
                    encoded_vector.append(str(k+1)+':'+str(v))
        encoded_vector.append('#'+ID)
        encoded_vector=' '.join(encoded_vector)
        return encoded_vector

    def write_input(self):
        out=open(self.input_file,'w')
        for ID in self.ID_ordered_list:
            out.write(self._encode_vector(ID)+'\n')
        out.close()
        
    def test(self):
        cmd=' svm_classify '+self.input_file+' '+self.path+self.modname+' '+self.output_file
        r,w=os.popen2(cmd)
        w.readlines()
        r.close()
        w.close()
        
    def retrieve(self):
        '''reads and returns the predicted values'''
        
        self.Dresults={}
        out=open(self.output_file).readlines()
        for i in xrange(len(out)):
            self.Dresults[self.ID_ordered_list[i]]=float(out[i])

class Post_SVM_learning:
    
    def __init__(self):
        pass
    
    def balance(self,LID_pos,LID_neg,Dprediction):
        '''bacello-like '''
        import operator
        
        Lpositive=[]
        for ID in LID_pos:
            Lpositive.append(Dprediction[ID])

        Lnegative=[]
        for ID in LID_neg:
            Lnegative.append(Dprediction[ID])
        
        def _calculate_matrix(Lpos,Lneg,thresh):
            M=[[0.,0.],[0.,0.]]
            
            Lpos_mod=map(operator.add,Lpos,[thresh]*len(Lpos))
            M[0][0]=sum(map(operator.ge,Lpos_mod,[0]*len(Lpos_mod)))
            M[0][1]=len(Lpos_mod)-M[0][0]
            Lneg_mod=map(operator.add,Lneg,[thresh]*len(Lneg))
            M[1][0]=sum(map(operator.ge,Lneg_mod,[0]*len(Lneg_mod)))
            M[1][1]=len(Lneg_mod)-M[1][0]
            
            #if (M[0][0]!=0 and  M[1][0]!=0) or ((M[0][1]!=0 and  M[1][1]!=0)):
            try:
                PP=Calculate_Prediction_Performance(M)
                PP._calc_all()
                
                return {'Cov':PP.Cov,
                        'Acc':PP.Acc,
                        'nAcc':PP.nAcc,
                        'MCC':PP.MCC,
                        'IC':PP.IC,
                        'Spec':PP.spec,
                        'FPrate':PP.LFPrate,
                        'Q':PP.Q,
                        'nQ':PP.nQ,
                        'GC':PP.GC,
                        'GC_norm':PP.GC_norm,
                        'Matrix':M,
                        'valid':True,
                        'threshold':thresh}                        
            except:
                return {'Cov':[0.,0.],
                    'Acc':[0.,0.],
                    'nAcc':[0.,0.],
                    'MCC':[0.,0.],
                    'IC':0.,
                    'Spec':[0.,0.],
                    'FPrate':[0.,0.],
                    'Q':0.,
                    'nQ':0.,
                    'GC':0.,
                    'GC_norm':0.,
                    'Matrix':[[0,0],[0,0]],
                    'valid':False,
                    'threshold':thresh}
            
        #search range

        Lint=map(operator.div,range(-500,505,5),[100.0]*len(range(-500,505,5)))
        Lperf=map(_calculate_matrix,[Lpositive]*len(Lint),[Lnegative]*len(Lint),Lint)
        #choose balancing parameter
        Lperf_treshold=[]
        for perf in Lperf:
            if perf['valid']:
                Lperf_treshold.append((perf['GC_norm'],perf['threshold']))
        threshold=max(Lperf_treshold)[1]
        self.performance_scale={}
        for perfD in Lperf:
            if perfD['valid']:
                self.performance_scale[perfD['threshold']]=perfD
        self.balancing_thresh=threshold
        
        return self.balancing_thresh,self.performance_scale


class Calculate_Prediction_Performance: 
    '''prende ininput una matrice di confusione NxN'''
    import warnings
    warnings.simplefilter("ignore","all")


    def __init__(self,M):
        ''' inizilizza la patrice calcolandone le cov,Acc,nAcc piu' GC e IC'''
        self.matrice=M
        
        
        
    def _calc_all(self):
        self.calCovAcc()
        self.calGC()
        self.calIC()
        self.calMCC()
        self.calGC_norm()
        self.calSpec()
        
        

    def calCovAcc(self):#calcola Coverage, Accuratezza e Accuratezza Normalizzata
        '''Prende in input la matrice di confuzione e calcola Coverage, Accuratezza e Accuratezza Normalizzata
        risultati:
        self.Cov=Lcov------>lista contenente le Coverage per colonna
        self.Acc=Lacc------>lista contenente le Accuracy per colonna
        self.nAcc=Laccper-->lista contenente le Normalized Accuracy per colonna
        self.MCC=LMMC------>lista contenente le MCC per colonna
        self.spec=spec----->lista contenente le tuple di specificity e FP rate '''

        K=len(self.matrice)
        Lcov=[]
        for i in range(K):
            corr=self.matrice[i][i]
            tot=0
            for j in range(K):
                tot=tot+self.matrice[i][j]
            try:
                cov=float(corr)/tot
            except ZeroDivisionError:
                cov=0.0
            Lcov.append(cov)
            
        Lacc=[]
        for j in range(K):
            corr=self.matrice[j][j]
            tot=0
            for i in range(K):
                tot=tot+self.matrice[i][j]
            try:
                acc=float(corr)/tot
            except ZeroDivisionError:
                acc=0.0
            Lacc.append(acc)
            
        Laccper=[]
        Ltotc=[]
        for i in range(K):
            totclasse=0
            for j in range(K):
                totclasse=totclasse+self.matrice[i][j]
            Ltotc.append(totclasse)
        for j in range(K):
            corr=self.matrice[j][j]
            tot=0
            for i in range(K):
                try:
                    tot=tot+(float(self.matrice[i][j])/Ltotc[i])
                except ZeroDivisionError:
                    None
            try:
                acc=((float(corr)/Ltotc[j])/tot)
            except ZeroDivisionError:
                acc=0.0
            Laccper.append(acc)
        
        self.Cov=Lcov
        self.Acc=Lacc
        self.nAcc=Laccper
        
        
        tot=0.0
        for k in range(len(self.matrice)):
            for s in range(len(self.matrice[k])):
                tot=tot+self.matrice[k][s]
        
        c=0.0
        for k in range(len(self.matrice)):
            c=c+self.matrice[k][k]
        
        Qper=0.0
        for cc in Lcov:
            Qper=Qper+cc
        Qper=Qper/len(Lcov)
        
        self.elementi=tot
        self.corretti=c
        self.Q=c/tot
        self.nQ=Qper
        
        
        return 
        
    def calMCC(self): #CONTROLLARE A VOLTE VIENE > 1
        '''prende in input una matrice e calcola l'MCC per ogni classe considerandola 1wsAll'''
        import math
        K=len(self.matrice)
        LMCC=[]
        for k in range(K):
            #generazione matrice 2x2:
            M=[0,0,0,0]#TP,FN,FP,TN
            for i in range(K):
                for j in range(K):
                    if i==j==k:#TP
                        M[0]=self.matrice[i][j]
                    elif i==k:#FN
                        M[1]+=self.matrice[i][j]
                    elif j==k:#FP
                        M[2]+=self.matrice[i][j]
                    else:#TN
                        M[3]+=self.matrice[i][j]
            try:
                TP=M[0]
                TN=M[3]
                FP=M[2]
                FN=M[1]
                MCC=(TP*TN-FP*FN)/math.sqrt((TP+FP)*(TN+FN)*(TP+FN)*(TN+FP))
            except ZeroDivisionError:
                MCC=0.0
            LMCC.append(MCC)
        self.MCC=LMCC

    
    def calSpec(self):
        '''prende in input una matrice e calcola FP rate e specificity per ogni classe considerandola 1wsAll'''
        K=len(self.matrice)
        Lspec=[]
        LFPrate=[]
        for k in range(K):
            #generazione matrice 2x2:
            M=[0,0,0,0]#TP,FN,FP,TN
            for i in range(K):
                for j in range(K):
                    if i==j==k:#TP
                        M[0]=self.matrice[i][j]
                    elif i==k:#FN
                        M[1]+=self.matrice[i][j]
                    elif j==k:#FP
                        M[2]+=self.matrice[i][j]
                    else:#TN
                        M[3]+=self.matrice[i][j]
            try:
                TN=M[3]
                FP=M[2]
                FPrate=(float(FP)/(FP+TN))
                
            except ZeroDivisionError:
                FPrate=0.0
            spec=1.0-FPrate
                
            #Lspec.append((spec,FPrate))
            Lspec.append(spec)
            LFPrate.append(FPrate)
        self.spec=Lspec
        self.LFPrate=LFPrate
        
    def calGC(self):
        '''prende in input una matrice e ne calcola la Generalized correlation:
        risultati:
        self.GC=generalized Correlation'''
        
        
        
        import math
        K=len(self.matrice)
        #calcolo N
        N=0.0
        for i in range(K):
            for j in range(K):
                N=N+self.matrice[i][j]
                
        #calcolo e_ij
        e=[]#matrice che conterra' gli e_ij
        for i in range(K):
            eriga=[]
            for j in range(K):
                #x_i
                x=0.0
                for ix in range(K):
                    x=x+self.matrice[i][ix]
                #y_i
                y=0.0
                for iy in range(K):
                    y=y+self.matrice[iy][j]
                try:
                    epos=(x*y)/N
                except ZeroDivisionError:
                    epos=0.0
                eriga.append(epos)
            e.append(eriga)    
        
        #calcolo CG^2
        
        somNUM=0.0
        DEN=(N*(K-1))
        for i in range(K):
            for j in range(K):
                if e[i][j] != 0:
                    try:
                        NUM=((self.matrice[i][j]- e[i][j])**2)/e[i][j]
                    except ZeroDivisionError:
                        NUM=0.0
                    somNUM=somNUM+NUM
        try:
            GC2=somNUM/DEN
        except ZeroDivisionError:
            GC2=0.0
        GC=math.sqrt(GC2)
        
        self.GC=GC
        
        return


    def calIC(self):
        '''prende in input una matrice e ne calcola la il coefficiente di mutua informazione:
        risultati:
        self.IC=Mutual Information Coefficient'''
        
        import math
        K=len(self.matrice)
        #calcolo N
        N=0.0
        for i in range(K):
            for j in range(K):
                N=N+self.matrice[i][j]
        I=0.0
        H=0.0
        for i in range(K):
            #x_i
            x=0.0
            for ix in range(K):
                x=x+self.matrice[i][ix]
            
            A=0.0
            try:
                A=(0.0-((x/N)*math.log(x/N)))
            except :
                A=0.0
            
            B=0.0
            for j in range(K):
                
                #y_i
                y=0.0
                for iy in range(K):
                    y=y+self.matrice[iy][j]
                    



                try:
                    Bj=(self.matrice[i][j]/N)*math.log(self.matrice[i][j]/y)
                except :
                    Bj=0.0
                
                B=B+Bj
            Id=A+B
            I=I+Id
        
            Hd=A
            H=H+Hd
        
        try:
            IC=I/H
        except ZeroDivisionError:
            IC=0.0
        self.IC=IC
            
        return 


    def redmat(self,Lind):
    
        '''prende in input una matrice di confusione e una lista contenente i raggruppamenti
        Lind=[[0,1],[2,3,],[4]] #primo livello: nuovo elemento della matrice formato dai gruppi del secondo livello che saranno accorpati
        Da in output una matrice di confusione con numero di cruppi ridotto alla 
        lunghezza di Lind accorpando correttamente i TP e i FP
        NOTA: i FP che occorrono tra le gli elementi raggruppati in un gruppo (es. 0 e 1 o 2 e 3) vengono aggiunti ai TP di quel gruppo'''
        ''' una volta ridotta modifica l'oggetto facendolo diventare la matrice ridotta e ne ricalcola le cov,Acc,nAcc piu' GC e IC
        risultati:
            self=Mred[:]
            self.calCovAcc()
            self.GC()
            self.IC()'''
        
        LTP=[]
        LFP=[]
        #calcolo TP e FP
        for k in range(len(Lind)):#gruppo di indici
            lTP=[]
            TP=0 #somma diagonali
            for t in range(len(Lind[k])):
                for s in range(len(Lind[k])):
                    i=Lind[k][t]
                    j=Lind[k][s]
                    TP=TP+self.matrice[i][j]
            lTP.append(TP)
            
            
            lFP=[]
            for s in range(len(Lind)):
                if s != k: #se non e' un TP
                    Lidfp=Lind[s]
                    FP=0 #somma colonne
                    for t in range(len(Lind[k])):
                        i=Lind[k][t] #indice gruppo
                        for v in range(len(Lidfp)):
                            j=Lidfp[v]
                            FP=FP+self.matrice[j][i]
                    lFP.append(FP)
            LTP.append(lTP)
            LFP.append(lFP)
        #ricostruzione matrice ridotta
        Mred=[]
        for k in range(len(Lind)):
            lriga=[]
            for t in range(len(Lind)):
                if k==t:
                    lriga.append(LTP[k][0])
                else:
                    if k>t:
                        lriga.append(LFP[t][k-1])
                    else:
                        lriga.append(LFP[t][k])
    
            Mred.append(lriga)
        
        self.matrice_rid=Mred[:]
        self.calCovAcc()
        self.calGC()
        self.calIC()
        
        
        return 


    def calGC_norm(self,):
        '''prende in input una matrice e ne calcola la Generalized correlation dopo averla normalizzata:
        risultati:
        self.GC=generalized Correlation'''
        
        import math,copy


        mat_norm=copy.deepcopy(self.matrice)
        
            

        #normalizza matrice                    
        for i in range(len(mat_norm)):
            riga1=mat_norm[i]
            tot=0.0
            for val in riga1:
                tot=tot+val
            for k in range(len(riga1)):
                try:
                    mat_norm[i][k]=float(riga1[k])/tot
                except ZeroDivisionError:
                    mat_norm[i][k]=0.0
        #calcolo
        K=len(mat_norm)
        #calcolo N
        N=0.0
        for i in range(K):
            for j in range(K):
                N=N+mat_norm[i][j]
                
        #calcolo e_ij
        e=[]#matrice che conterra' gli e_ij
        for i in range(K):
            eriga=[]
            for j in range(K):
                #x_i
                x=0.0
                for ix in range(K):
                    x=x+mat_norm[i][ix]
                #y_i
                y=0.0
                for iy in range(K):
                    y=y+mat_norm[iy][j]
                try:
                    epos=(x*y)/N
                except ZeroDivisionError:
                    epos=0.0
                eriga.append(epos)
            e.append(eriga)    
        
        #calcolo CG^2
        
        somNUM=0.0
        DEN=(N*(K-1))
        for i in range(K):
            for j in range(K):
                if e[i][j] != 0:
                    try:
                        NUM=((mat_norm[i][j]- e[i][j])**2)/e[i][j]
                    except ZeroDivisionError:
                        NUM=0.0
                    somNUM=somNUM+NUM
        try:
            GC2=somNUM/DEN
        except ZeroDivisionError:
            GC2=0.0
        GC=math.sqrt(GC2)
        
        self.GC_norm=GC
