from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Corrfunc.theory.DD import DD
from Corrfunc.io import read_catalog
from Corrfunc.utils import convert_3d_counts_to_cf



pth='/home/ignacio/Documents/HAB/Catalogs/dfG_halos.csv'
pth2='/home/ignacio/Documents/HAB/Catalogs/'
cat=pd.read_csv(pth)
cat.sort_values(by=['Group_M_Crit200'],ascending=True,inplace=True)
bins=[]
ls=[]
for i in cat.index:
    ls.append(np.log10(cat['Group_M_Crit200'][i]))
cat['logmass']=ls
for i in range(23): #22 bins
    b=min(cat['logmass'])+ (max(cat['logmass'])-min(cat['logmass'])) * (i) /22  
    b=round(b,13)
    bins.append(b)
    print(i)
cat['mass_bins'] = pd.cut(cat['logmass'], bins)
#cat.to_csv(pth2+'dfG_halos.csv',index=False)

###Lo de arriba es un código con el que generamos bin de masa. Archivos en mass_bins.


pth='/home/ignacio/Documents/HAB/Catalogs/Illustris_data/bins_spin/'
num=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
bns=[]
for i in num:
    data=pd.read_csv(pth+'bins_'+str(i)+'.csv')
    bns.append([data,'bins_'+str(i)]) #contiene 22 elementos

quants25=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
quants75=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
for i in bns:
    if i[1]=='bins_1':
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[0].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[0].append([i[0]['Spin'][ind],i[1],'<q25'])
    elif i[1]=='bins_2':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[1].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[1].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_3':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[2].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[2].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_4':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[3].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[3].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_5':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[4].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[4].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_6':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[5].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[5].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_7':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[6].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[6].append([i[0]['Spin'][ind],i[1],'<q25']) 
    elif i[1]=='bins_8':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[7].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[7].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_9':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[8].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[8].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_10':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[9].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[9].append([i[0]['Spin'][ind],i[1],'<q25']) 
    elif i[1]=='bins_11':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[10].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[10].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_12':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[11].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[11].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins13':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[12].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[12].append([i[0]['Spin'][ind],i[1],'<q25']) 
    elif i[1]=='bins_14':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[13].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[13].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins15':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[14].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[14].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_16':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[15].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[15].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_17':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[16].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[16].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_18':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[17].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[17].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_19':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[18].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[18].append([i[0]['Spin'][ind],i[1],'<q25']) 
    elif i[1]=='bins_20':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[19].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[19].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_21':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[20].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[20].append([i[0]['Spin'][ind],i[1],'<q25'])  
    elif i[1]=='bins_22':   
        q25,q75=np.percentile(i[0]['Spin'],[25,75])
        for ind in i[0].index:
            if i[0]['Spin'][ind]>=q75:
                quants75[21].append([i[0]['Spin'][ind],i[1],'>q75'])
            elif i[0]['Spin'][ind]<=q25:
                quants25[21].append([i[0]['Spin'][ind],i[1],'<q25']) 



sp75=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]] 
for bin in range(len(quants75)):
    for num_bin in quants75[bin]:
        numero_bin=num_bin[1]
        if numero_bin=='bins_1':
            sp75[0].append(num_bin[0])
        elif numero_bin=='bins_2':
            sp75[1].append(num_bin[0])
        elif numero_bin=='bins_3':
            sp75[2].append(num_bin[0])
        elif numero_bin=='bins_4':
            sp75[3].append(num_bin[0])
        elif numero_bin=='bins_5':
            sp75[4].append(num_bin[0])
        elif numero_bin=='bins_6':
            sp75[5].append(num_bin[0])
        elif numero_bin=='bins_7':
            sp75[6].append(num_bin[0])
        elif numero_bin=='bins_8':
            sp75[7].append(num_bin[0])
        elif numero_bin=='bins_9':
            sp75[8].append(num_bin[0])
        elif numero_bin=='bins_10':
            sp75[9].append(num_bin[0])
        elif numero_bin=='bins_11':
            sp75[10].append(num_bin[0])
        elif numero_bin=='bins_12':
            sp75[11].append(num_bin[0])
        elif numero_bin=='bins_13':
            sp75[12].append(num_bin[0])
        elif numero_bin=='bins_14':
            sp75[13].append(num_bin[0])
        elif numero_bin=='bins_15':
            sp75[14].append(num_bin[0])
        elif numero_bin=='bins_16':
            sp75[15].append(num_bin[0])
        elif numero_bin=='bins_17':
            sp75[16].append(num_bin[0])
        elif numero_bin=='bins_18':
            sp75[17].append(num_bin[0])
        elif numero_bin=='bins_19':
            sp75[18].append(num_bin[0])
        elif numero_bin=='bins_20':
            sp75[19].append(num_bin[0])
        elif numero_bin=='bins_21':
            sp75[20].append(num_bin[0])
        elif numero_bin=='bins_22':
            sp75[21].append(num_bin[0])

min75=[]
for i in range(len(sp75)): #Arroja el minimo del quantil 75
    print(min(sp75[i]))
    min75.append(min(sp75[i]))
    
sp25=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]] 
for bin in range(len(quants25)):
    for num_bin in quants25[bin]:
        numero_bin=num_bin[1]
        if numero_bin=='bins_1':
            sp25[0].append(num_bin[0])
        elif numero_bin=='bins_2':
            sp25[1].append(num_bin[0])
        elif numero_bin=='bins_3':
            sp25[2].append(num_bin[0])
        elif numero_bin=='bins_4':
            sp25[3].append(num_bin[0])
        elif numero_bin=='bins_5':
            sp25[4].append(num_bin[0])
        elif numero_bin=='bins_6':
            sp25[5].append(num_bin[0])
        elif numero_bin=='bins_7':
            sp25[6].append(num_bin[0])
        elif numero_bin=='bins_8':
            sp25[7].append(num_bin[0])
        elif numero_bin=='bins_9':
            sp25[8].append(num_bin[0])
        elif numero_bin=='bins_10':
            sp25[9].append(num_bin[0])
        elif numero_bin=='bins_11':
            sp25[10].append(num_bin[0])
        elif numero_bin=='bins_12':
            sp25[11].append(num_bin[0])
        elif numero_bin=='bins_13':
            sp25[12].append(num_bin[0])
        elif numero_bin=='bins_14':
            sp25[13].append(num_bin[0])
        elif numero_bin=='bins_15':
            sp25[14].append(num_bin[0])
        elif numero_bin=='bins_16':
            sp25[15].append(num_bin[0])
        elif numero_bin=='bins_17':
            sp25[16].append(num_bin[0])
        elif numero_bin=='bins_18':
            sp25[17].append(num_bin[0])
        elif numero_bin=='bins_19':
            sp25[18].append(num_bin[0])
        elif numero_bin=='bins_20':
            sp25[19].append(num_bin[0])
        elif numero_bin=='bins_21':
            sp25[20].append(num_bin[0])
        elif numero_bin=='bins_22':
            sp25[21].append(num_bin[0])
max25=[]
for i in range(len(sp25)): #Arroja el maximo del quantil 25
    print(max(sp25[i]))
    max25.append(max(sp25[i]))


pth='/home/ignacio/Documents/HAB/Catalogs/Illustris_data/bins_spin/'
num=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
upper75=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]] 
for i in num:
    data=pd.read_csv(pth+'bins_'+str(i)+'.csv')
    data=pd.DataFrame(data)
    if i==1:
        df_mask=data['Spin']>=min75[0]
        df_split=data[df_mask]
        upper75[0].append(df_split) 
    elif i==2:
        df_mask=data['Spin']>=min75[1]
        df_split=data[df_mask]
        upper75[1].append(df_split) 
    elif i==3:
        df_mask=data['Spin']>=min75[2]
        df_split=data[df_mask]
        upper75[2].append(df_split)
    elif i==4:
        df_mask=data['Spin']>=min75[3]
        df_split=data[df_mask]
        upper75[3].append(df_split)
    elif i==5:
        df_mask=data['Spin']>=min75[4]
        df_split=data[df_mask]
        upper75[4].append(df_split)
    elif i==6:
        df_mask=data['Spin']>=min75[5]
        df_split=data[df_mask]
        upper75[5].append(df_split)
    elif i==7:
        df_mask=data['Spin']>=min75[6]
        df_split=data[df_mask]
        upper75[6].append(df_split)
    elif i==8:
        df_mask=data['Spin']>=min75[7]
        df_split=data[df_mask]
        upper75[7].append(df_split) 
    elif i==9:
        df_mask=data['Spin']>=min75[8]
        df_split=data[df_mask]
        upper75[8].append(df_split)
    elif i==10:
        df_mask=data['Spin']>=min75[9]
        df_split=data[df_mask]
        upper75[9].append(df_split)
    elif i==11:
        df_mask=data['Spin']>=min75[10]
        df_split=data[df_mask]
        upper75[10].append(df_split)
    elif i==12:
        df_mask=data['Spin']>=min75[11]
        df_split=data[df_mask]
        upper75[11].append(df_split)
    elif i==13:
        df_mask=data['Spin']>=min75[12]
        df_split=data[df_mask]
        upper75[12].append(df_split)
    elif i==14:
        df_mask=data['Spin']>=min75[13]
        df_split=data[df_mask]
        upper75[13].append(df_split) 
    elif i==15:
        df_mask=data['Spin']>=min75[14]
        df_split=data[df_mask]
        upper75[14].append(df_split)
    elif i==16:
        df_mask=data['Spin']>=min75[15]
        df_split=data[df_mask]
        upper75[15].append(df_split)
    elif i==17:
        df_mask=data['Spin']>=min75[16]
        df_split=data[df_mask]
        upper75[16].append(df_split)
    elif i==18:
        df_mask=data['Spin']>=min75[17]
        df_split=data[df_mask]
        upper75[17].append(df_split)
    elif i==19:
        df_mask=data['Spin']>=min75[18]
        df_split=data[df_mask]
        upper75[18].append(df_split)
    elif i==20:
        df_mask=data['Spin']>=min75[19]
        df_split=data[df_mask]
        upper75[19].append(df_split)
    elif i==21:
        df_mask=data['Spin']>=min75[20]
        df_split=data[df_mask]
        upper75[20].append(df_split)
    elif i==22:
        df_mask=data['Spin']>=min75[21]
        df_split=data[df_mask]
        upper75[21].append(df_split)
        
pth='/home/ignacio/Documents/HAB/Catalogs/Illustris_data/bins_spin/'
num=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
lower25=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]] 
for i in num:
    data=pd.read_csv(pth+'bins_'+str(i)+'.csv')
    data=pd.DataFrame(data)
    if i==1:
        df_mask=data['Spin']<=max25[0]
        df_split=data[df_mask]
        lower25[0].append(df_split) 
    elif i==2:
        df_mask=data['Spin']>=max25[1]
        df_split=data[df_mask]
        lower25[1].append(df_split) 
    elif i==3:
        df_mask=data['Spin']>=max25[2]
        df_split=data[df_mask]
        lower25[2].append(df_split)
    elif i==4:
        df_mask=data['Spin']>=max25[3]
        df_split=data[df_mask]
        lower25[3].append(df_split)
    elif i==5:
        df_mask=data['Spin']>=max25[4]
        df_split=data[df_mask]
        lower25[4].append(df_split)
    elif i==6:
        df_mask=data['Spin']>=max25[5]
        df_split=data[df_mask]
        lower25[5].append(df_split)
    elif i==7:
        df_mask=data['Spin']>=max25[6]
        df_split=data[df_mask]
        lower25[6].append(df_split)
    elif i==8:
        df_mask=data['Spin']>=max25[7]
        df_split=data[df_mask]
        lower25[7].append(df_split) 
    elif i==9:
        df_mask=data['Spin']>=max25[8]
        df_split=data[df_mask]
        lower25[8].append(df_split)
    elif i==10:
        df_mask=data['Spin']>=max25[9]
        df_split=data[df_mask]
        lower25[9].append(df_split)
    elif i==11:
        df_mask=data['Spin']>=max25[10]
        df_split=data[df_mask]
        lower25[10].append(df_split)
    elif i==12:
        df_mask=data['Spin']>=max25[11]
        df_split=data[df_mask]
        lower25[11].append(df_split)
    elif i==13:
        df_mask=data['Spin']>=max25[12]
        df_split=data[df_mask]
        lower25[12].append(df_split)
    elif i==14:
        df_mask=data['Spin']>=max25[13]
        df_split=data[df_mask]
        lower25[13].append(df_split) 
    elif i==15:
        df_mask=data['Spin']>=max25[14]
        df_split=data[df_mask]
        lower25[14].append(df_split)
    elif i==16:
        df_mask=data['Spin']>=max25[15]
        df_split=data[df_mask]
        lower25[15].append(df_split)
    elif i==17:
        df_mask=data['Spin']>=max25[16]
        df_split=data[df_mask]
        lower25[16].append(df_split)
    elif i==18:
        df_mask=data['Spin']>=max25[17]
        df_split=data[df_mask]
        lower25[17].append(df_split)
    elif i==19:
        df_mask=data['Spin']>=max25[18]
        df_split=data[df_mask]
        lower25[18].append(df_split)
    elif i==20:
        df_mask=data['Spin']>=max25[19]
        df_split=data[df_mask]
        lower25[19].append(df_split)
    elif i==21:
        df_mask=data['Spin']>=max25[20]
        df_split=data[df_mask]
        lower25[20].append(df_split)
    elif i==22:
        df_mask=data['Spin']>=max25[21]
        df_split=data[df_mask]
        lower25[21].append(df_split)
data_bins=[]
for i in range(22):
    data_bins.append([lower25[i],upper75[i],bns[i][0]])
data_bins[1]
