import numpy as np
import matplotlib.pyplot as plt

a=np.array([3,2])
b=np.array([18,42,24])
c=np.array([(2,1),(2,3),(3,1)])
a1=np.array([1,3])
b1=np.array([10,5,4])
c1=np.array([(1,2),(1,0),(0,1)])

class LPSsolver():
    
    def solve(self,a,b,c):
        
        la=len(a)   #no. of variables
        lb=len(b)   #no. of equations
        x=[0 for i in range(la)]        #final answer
        graph=[]
        b=b[:,np.newaxis]
        arr=np.identity(lb)
        arr=np.concatenate((c,arr,b),axis=1)
        C=np.concatenate([a,np.array([0 for i in range(lb)])])
        B=np.array([0 for i in range(lb)])
        bx=list(range(la,la+lb))
        minpos=[0 for i in range(lb)]
        maxpos=[0 for i in range(la+lb)]
        z=np.array([[0.0 for i in range(la+lb+1)]])
        
        for i in range(la+lb+1):
            if i<la: z[0][i]=-a[i]
        
        arr=np.concatenate((arr,z),axis=0)
        z1,k=0,0
        ann1=[]
        
        while 1:
            i,j=0,0
            for i in range(la+lb):
                maxpos[i]=C[i]-sum(arr[:-1,i]*B)
        
            maxmaxpos=max(maxpos)
           
            
            if maxmaxpos<=0: break
            imax=maxpos.index(maxmaxpos)
            for i in range(lb):
                if arr[i,imax]==0:
                    minpos[i]=1e40
                else:
                    minpos[i]=arr[i,-1]/arr[i,imax]
                    
                    
            
            qminpos=filter(lambda x:x>0,minpos)
            if qminpos==[]: break
            minminpos=min(qminpos)
            imin=minpos.index(minminpos)
            B[imin]=C[imax]
            bx[imin]=imax
            arr[imin]=arr[imin]/arr[imin,imax]
            
            for i in range(0,imin)+range(imin+1,lb+1):
                arr[i]=arr[i]-arr[i,imax]*arr[imin]/arr[imin,imax]
            
            graph.append([0]*la)
            for i in range(la):
                if i in bx:
                    graph[k][j]=arr[bx.index(i),-1]
                else:
                    graph[k][j]=0
                j+=1
            k+=1
                    
            
            ann1.append(arr[-1,-1])        
            if z1<arr[-1,-1]:
                z1=arr[-1,-1]
            else: break
            
            
         
          
        for i in range(lb):
            if bx[i]<la:
                x[bx[i]]=arr[i,-1]
        x.append(sum([x[i]*a[i] for i in range(la)]))
        
        
        graph1=np.array(graph)
        
        extral=[]
        if 0 in graph1[:,1]:
            for i in range(lb):
                if c[i,0]!=0:
                    extral.append(b[i,0]/c[i,1])
            extra=min(extral)
            graph.append([0,extra])
            ann1.append(a[1]*extra)
        elif 0 in graph1[:,0]:
            for i in range(lb):
                if c[i,0]!=0:
                    extral.append(b[i,0]/c[i,0])
            extra=min(extral)
            graph.append([extra,0])
            ann1.append(a[0]*extra)
                
        graph=np.array(graph)
        
        xmax=max(graph[:,0])
        ymax=max(graph[:,1])
        
        #Defining Subplots
        fig=plt.figure()
        fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)
        constraints=fig.add_subplot(221)
        sol=fig.add_subplot(223)
        final=fig.add_subplot(122)
        
        #Defining line names
        linename=[0]*3
        for i in range(3):
            tempx=str(c[i,0])
            tempy=str(c[i,1])
            if tempx=='1': tempx=''
            if tempy=='1': tempy=''
            if tempx=='0':
                linename[i]="%sy-%d=0"%(tempy,b[i])
            elif tempy=='0':
                linename[i]="%sx-%d=0"%(tempx,b[i])
            else:
                linename[i]="%sy+%sx-%d=0"%(tempy,tempx,b[i])
        
        
        
        
        #Plotting constraints
        t=np.arange(0,xmax+5,0.1)
        ycon=[]
        xcon=[]
        yconmax=[]
        linecon=[0]*3
        for i in range(3):
            if c[i][1]==0:
                ycon.append(np.arange(0,int(ymax)+10,0.1))
                xcon.append([float(b[i]/c[i][0]) for j in ycon[i]])
            elif c[i][0]==0:
                xcon.append(t)
                ycon.append([float(b[i]/c[i][1]) for j in t])
            else:
                xcon.append(t)
                ycon.append((b[i]-c[i][0]*t)/c[i][1])
                yconmax.append(max(ycon[i]))
        yconmax=max(yconmax)
        
        st=['red','blue','green']
        for i in range(3):
            print(xcon[i])
            print(ycon[i])
            linecon[i],=constraints.plot(xcon[i],ycon[i])
            constraints.fill_between(xcon[i],ycon[i],facecolor=st[i], alpha=0.5)
        constraints.axis([0,xmax+4,0,yconmax])
        constraints.annotate('Constraint Equations',xy=(xmax-3.7,yconmax-2), xytext=(xmax-3.7,yconmax-3))
        constraints.legend()

            
        #Plotting sol
        sol.plot(graph[:,0],graph[:,1],'ro',markersize=8,label="Z Achieved At Each Iteration")
        sol.plot(graph[:,0],graph[:,1],'k-',marker='o',markerfacecolor='black', markersize=8)
        sol.plot(graph[-2,0],graph[-2,1],'ro',markersize=10)
        sol.fill_between(graph[:,0],graph[:,1],)
        sol.axis([0,xmax+2.5,0,ymax+3])
        sol.annotate('Z at each iteration',xy=(xmax-3.8,ymax+1), xytext=(xmax-3.7,ymax+1))
        for i in range(0,len(graph)):
            ann="Z = %d"%(ann1[i])
            
            sol.annotate(ann, xy=(graph[i,0], graph[i,1]), xytext=(graph[i,0]+.2, graph[i,1]+.3))
        
        #Defining Final
        znormal,=final.plot(graph[:,0],graph[:,1],'ko',markersize=8,)
        final.plot(graph[:,0],graph[:,1],'k-',marker='o',markerfacecolor='black', markersize=8)
        zmax,=final.plot(graph[-2,0],graph[-2,1],'ro',markersize=10)
        final.fill_between(graph[:,0],graph[:,1],)
        final.axis([0,xmax+1,0,yconmax+3])
        
        for i in range(3):
            final.plot(xcon[i],ycon[i])
            final.fill_between(xcon[i],ycon[i],facecolor=st[i], alpha=0.2)
        
    
        
        #Defining legends
        fig.legend(linecon+[znormal,zmax],linename+['Z at each Iteration','Maximum Z'],numpoints=1,loc='upper center',prop={'size':12},fancybox=True,shadow=True,ncol=5,title='SIMPLEX METHOD').get_frame().set_alpha(0.5)
        
        plt.axis([0,xmax+5,0,ymax+5])

        
            
        
        
        
        
        ann="Maximum Z = %d"%(x[-1])
        plt.annotate(ann, xy=(graph[-2,0]+.1, graph[-2,1]+.1), xytext=(graph[-2,0]+2, graph[-2,1]+2),arrowprops=dict(facecolor='black', shrink=0.05),)
        for i in range(0,len(graph)-2)+range(len(graph)-1,len(graph)):
            ann="Z = %d"%(ann1[i])
            
            plt.annotate(ann, xy=(graph[i,0], graph[i,1]), xytext=(graph[i,0]+.2, graph[i,1]+.3))
        
        plt.show()
        x=map(int,x)
        return x
    
lps=LPSsolver()
solution=lps.solve(a,b,c)
print(solution)
            
            
                
        
